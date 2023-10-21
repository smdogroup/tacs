"""
==============================================================================
TACS Nonlinear Continuation Solver
==============================================================================
This continuation solver uses a predictor-corrector scheme to increment the load scale in nonlinear problems.
Each iteration of the continuation solver consists of three steps:

#. Predictor computation: If enabled, the solver will use the solutions from previous continuation steps to extrapolate the equilibrium path to the current load scale, which should provide a good initial guess for the inner solver.
#. Corrector computation: The continuation solver calls an inner solver to solve the nonlinear problem at the current load scale.
#. Load scale update: The continuation solver increments the load scale, the size of the step taken is adapted each iteration to try and achieve a target number of inner solver iterations.
"""

# ==============================================================================
# Standard Python modules
# ==============================================================================
from typing import Optional, Callable, Any, Dict

# ==============================================================================
# External Python modules
# ==============================================================================
import numpy as np
import mpi4py

# ==============================================================================
# Extension modules
# ==============================================================================
import tacs.TACS
from tacs.solvers import BaseSolver


def lagrangeInterp(xKnown, yKnown, xQuery, yQuery):
    """Interpolate an array using lagrange polynomials

    Parameters
    ----------
    xKnown : iterable of length n
        scalar x values of known points
    yKnown : iterable of n np.ndarrays
        arrays at known points
    xQuery : float
        x value to interpolate at
    yQuery : np.ndarray
        array to store interpolated result in
    """
    numPoints = len(xKnown)
    yQuery[:] = 0.0
    yTemp = np.zeros_like(yQuery)
    for jj in range(numPoints):
        yTemp[:] = 1.0
        for mm in range(numPoints):
            if mm != jj:
                yTemp[:] = yTemp[:] * (xQuery - xKnown[mm]) / (xKnown[jj] - xKnown[mm])
        yQuery[:] = yQuery[:] + yKnown[jj] * yTemp


class ContinuationSolver(BaseSolver):
    defaultOptions = {
        "MaxLambda": [
            float,
            1.0,
            "Final continuation parameter value to aim for.",
        ],
        "AbsTol": [
            float,
            1e-8,
            "Convergence criteria for the nonlinear residual norm.",
        ],
        "RelTol": [
            float,
            1e-8,
            "Relative convergence criteria for the nonlinear residual norm, norm is measured relative to that of the external load vector.",
        ],
        "CoarseAbsTol": [
            float,
            1e-4,
            "Residual norm criteria for intermediate continuation steps, making this larger may speed up the nonlinear solver by allowing it to only partially converge intermediate steps.",
        ],
        "CoarseRelTol": [
            float,
            1e-4,
            "Relative residual norm criteria for intermediate load increments.",
        ],
        "TargetIter": [
            int,
            8,
            "Target number of Newton iterations for each continuation increment.",
        ],
        "MaxIter": [int, 30, "Maximum number of continuation steps."],
        "InitialStep": [float, 0.2, "Initial continuation step size."],
        "MinStep": [float, 1e-4, "Minimum continuation step size."],
        "MaxStep": [float, np.inf, "Maximum continuation step size."],
        "MinStepFactor": [
            float,
            0.5,
            "The minimum factor by which the continuation step size can decrease in a single step.",
        ],
        "MaxStepFactor": [
            float,
            2.0,
            "The maximum factor by which the continuation step size can increase in a single step.",
        ],
        "RetractionFactor": [
            float,
            0.5,
            "The factor by which the continuation step size is reduced when the Newton solver fails to converge.",
        ],
        # Predictor step options
        "UsePredictor": [
            bool,
            False,
            "Flag for using predictor step in continuation.",
        ],
        "NumPredictorStates": [
            int,
            2,
            "Number of previous equilibrium states to use in computing the predictor step.",
        ],
        # "predictorUseDerivative": [
        #     bool,
        #     False,
        #     "Whether to use the equilibrium path slope in the computation of the predictor step. This requires a linear solve and thus greatly increases the cost of the predictor step computation.",
        # ],
    }

    def __init__(
        self,
        jacFunc: Callable,
        pcUpdateFunc: Callable,
        linearSolver: tacs.TACS.KSM,
        setLambdaFunc: Callable,
        getLambdaFunc: Callable,
        innerSolver: BaseSolver,
        options: Optional[dict] = None,
        comm: Optional[mpi4py.MPI.Comm] = None,
    ) -> None:
        """Create a continuation solver instance

        Parameters
        ----------
        jacFunc : function
            Function to update the residual Jacobian at the current state, with signature `jacFunc() -> None`
        pcUpdateFunc : function
            Function to update the residual Jacobian preconditioner at the current state, with signature `pcUpdateFunc() -> None`
        linearSolver : tacs.TACS.KSM
            TACS linear solver object to use for the Newton solve, the linear solver owns the matrix and preconditioner
        setLambdaFunc : function
            Function to set the continuation parameter, with signature `setLambdaFunc(lambda:float) -> None`
        setLambdaFunc : function
            Function to get the current continuation parameter, with signature `getLambdaFunc() -> float`
        innerSolver : TACS Nonlinear Solver
            A Solver object to use for the corrector solve in each increment (e.g a NewtonSolver object)
        options : dict, optional
            Dictionary holding solver-specific option parameters (case-insensitive)., by default None
        comm : mpi4py.MPI.Intracomm, optional
            The comm object on which to create the pyTACS object., by default mpi4py.MPI.COMM_WORLD
        """

        self.jacFunc = jacFunc
        self.pcUpdateFunc = pcUpdateFunc
        self.linearSolver = linearSolver
        self.setLambdaFunc = setLambdaFunc
        self.getLambdaFunc = getLambdaFunc
        self.innerSolver = innerSolver

        self.equilibriumPathStates = []
        self.equilibriumPathLoadScales = []

        BaseSolver.__init__(
            self,
            assembler=self.innerSolver.assembler,
            setStateFunc=self.innerSolver.setStateFunc,
            resFunc=self.innerSolver.resFunc,
            stateVec=self.innerSolver.stateVec,
            resVec=self.innerSolver.resVec,
            options=options,
            comm=comm,
        )

        # Create additional vectors
        self.fInt = self.assembler.createVec()
        self.fExt = self.assembler.createVec()
        self.du_e = self.assembler.createVec()
        self.du_i = self.assembler.createVec()
        self.predictorStep = self.assembler.createVec()
        self.incStartState = self.assembler.createVec()

    @property
    def resFunc(self) -> Callable:
        return self.innerSolver.resFunc

    @resFunc.setter
    def resFunc(self, resFunc: Callable) -> None:
        self.innerSolver.resFunc = resFunc

    def getHistoryVariables(self) -> Dict[str, Dict]:
        """Get the variables to be stored in the solver history

        This method allows for implementation of any logic that dictates any changes in the stored variables depending on the current options.

        Returns
        -------
        Dict[str, Dict]
            Dictionary of solver variables, keys are the variable names, value is another dictionary with keys "type" and "print", where "type" is the data type of the variable and "print" is a boolean indicating whether or not to print the variable to the screen
        """
        variables = {}
        variables["Increment"] = {"type": int, "print": True}
        variables["Lambda"] = {"type": float, "print": True}
        variables["SubIter"] = {"type": int, "print": True}

        # Add the variables from the inner solver
        variables.update(self.innerSolver.getHistoryVariables())

        return variables

    def _setupPredictorVectors(self) -> None:
        """Setup the structures containing the data for computing the predictor steps"""
        self.equilibriumPathStates = []
        self.equilibriumPathLoadScales = []
        if self.getOption("UsePredictor"):
            for _ in range(self.getOption("NumPredictorStates")):
                self.equilibriumPathStates.append(self.assembler.createVec())
                self.equilibriumPathLoadScales.append(None)

    def setOption(self, name: str, value: Any) -> None:
        BaseSolver.setOption(self, name, value)

        # Update the predictor computation data structures if the relevant options are changed
        if name.lower() in [
            "usepredictor",
            "numpredictorstates",
        ]:
            self._setupPredictorVectors()

    def setConvergenceTolerance(
        self, absTol: Optional[float] = None, relTol: Optional[float] = None
    ) -> None:
        """Set the convergence tolerance of the solver

        Parameters
        ----------
        absTol : float, optional
            Absolute tolerance, not changed if no value is provided
        relTol : float, optional
            Relative tolerance, not changed if no value is provided
        """
        if absTol is not None:
            self.setOption("AbsTol", absTol)
        if relTol is not None:
            self.setOption("RelTol", relTol)

        return

    def initializeSolve(self, u0: Optional[tacs.TACS.Vec] = None) -> None:
        """Perform any initialization required before the solve"""
        BaseSolver.initializeSolve(self, u0)
        if self.getOption("UsePredictor"):
            for ii in range(self.getOption("NumPredictorStates")):
                self.equilibriumPathLoadScales[ii] = None
                self.equilibriumPathStates[ii].zeroEntries()

    def solve(
        self, u0: Optional[tacs.TACS.Vec] = None, result: Optional[tacs.TACS.Vec] = None
    ) -> None:
        MAX_LAMBDA = self.getOption("MaxLambda")
        TARGET_ITERS = self.getOption("TargetIter")
        INIT_STEP = self.getOption("InitialStep")
        MIN_STEP = self.getOption("MinStep")
        MAX_STEP = self.getOption("MaxStep")
        MAX_INCREMENTS = self.getOption("MaxIter")
        MIN_STEP_FACTOR = self.getOption("MinStepFactor")
        MAX_STEP_FACTOR = self.getOption("MaxStepFactor")
        STEP_RETRACT_FACTOR = self.getOption("RetractionFactor")

        ABS_TOL = self.getOption("AbsTol")
        REL_TOL = self.getOption("RelTol")
        COARSE_ABS_TOL = self.getOption("CoarseAbsTol")
        COARSE_REL_TOL = self.getOption("CoarseRelTol")

        USE_PREDICTOR = self.getOption("UsePredictor")
        # PREDICTOR_USE_DERIVATIVE = self.getOption("predictorUseDerivative")

        self.initializeSolve()

        if u0 is not None:
            self.stateVec.copyValues(u0)

        # Compute the internal and external forcing vectors at the current point
        self.computeForceVectors()
        if np.real(self.fExt.norm()) == 0:
            self.du_e.copyValues(self.stateVec)
            self.stateVec.zeroEntries()
            self.setStateFunc(self.stateVec)
            self.resFunc(self.resVec)
            self.setRefNorm(np.real(self.resVec.norm()))
        else:
            self.setRefNorm(np.real(self.fExt.norm()))

        # ==============================================================================
        # Compute the initial load scale
        # ==============================================================================
        self.setLambdaFunc(INIT_STEP)
        loadStepDirection = 1

        # If we're restarting from a previous solution we should compute the optimum load scale
        # to restart from. This is done by computing the load scale that minimizes the work
        # done by the resulting Newton step:
        # optLoadScale = (Fe^T dUi + Fi^T dUe) / (-2 Fe^T dUe)
        # Where: Fe = external force, Fi = internal force, dUi = inv(K) * Fi, dUe = inv(K) * Fe
        isRestartIncrement = False
        if np.real(self.stateVec.norm()) > 0 and np.real(self.fExt.norm()) > 0.0:
            self.jacFunc()
            self.pcUpdateFunc()
            self.linearSolver.solve(self.fExt, self.du_e)
            self.linearSolver.solve(self.fInt, self.du_i)
            FeUe = np.real(self.fExt.dot(self.du_e))
            FeUi = np.real(self.fExt.dot(self.du_i))
            FiUe = np.real(self.fInt.dot(self.du_e))
            optLoadScale = (FeUi + FiUe) / (-2 * FeUe)

            if optLoadScale > 2 * MAX_LAMBDA or optLoadScale < 0.0:
                # If the optimum load scale is more than double the max load scale we're aiming for, or if it's
                # negative then the loading/structure has changed so much that we'll be closer to the final
                # solution if we just reset the displacements to zero and start the solver from there
                self.stateVec.zeroEntries()
                self.setStateFunc(self.stateVec)
                optLoadScale = INIT_STEP
            elif np.abs(optLoadScale - MAX_LAMBDA) < 1e-2:
                # If the optimum load scale is close to the max load scale then we'll just use the max load scale
                optLoadScale = MAX_LAMBDA
                isRestartIncrement = True
            else:
                # Otherwise choose the maximum of the ideal load scale and the default initial load scale
                optLoadScale = max(optLoadScale, INIT_STEP)
                isRestartIncrement = True
                # If the optimum load scale is greater than the max we want to get to then we need to reverse the
                # direction of load incrementation
                if optLoadScale > MAX_LAMBDA:
                    loadStepDirection = -1

            self.setLambdaFunc(optLoadScale)
        elif np.real(self.fExt.norm()) == 0.0:
            # If the external force is zero then the load scale doesn't mean anything and we can just set it to 1
            self.setLambdaFunc(1.0)

        # If starting from zero, we can assume that u=0, lambda=0 is an equilibrium state
        if USE_PREDICTOR and self.stateVec.norm() == 0:
            self.equilibriumPathLoadScales[-1] = 0.0
            self.equilibriumPathStates[-1].copyValues(self.stateVec)

        stepSize = INIT_STEP
        currentLambda = self.getLambdaFunc()

        for increment in range(MAX_INCREMENTS):
            # Save displacement at start of this increment, this is what
            # we'll reset to if the increment diverges
            self.incStartState.copyValues(self.stateVec)

            # --- Compute predictor step ---
            # TODO: Adapt this to enable use of equilibrium path slope
            # We only compute a predictor step if we have at least 3 equilibrium states so that we can do a nonlinear extrapolation
            numValidPredictorStates = sum(
                [1 for x in self.equilibriumPathLoadScales if x is not None]
            )
            if numValidPredictorStates > 2:
                stateArrays = [
                    x.getArray()
                    for x in self.equilibriumPathStates[-numValidPredictorStates:]
                ]
                lagrangeInterp(
                    self.equilibriumPathLoadScales[-numValidPredictorStates:],
                    stateArrays,
                    currentLambda,
                    self.stateVec.getArray(),
                )
                self.setStateFunc(self.stateVec)

            # --- Call inner solver for corrector step ---
            if currentLambda == MAX_LAMBDA:
                rtol = REL_TOL
                atol = ABS_TOL
            else:
                rtol = COARSE_REL_TOL
                atol = COARSE_ABS_TOL
            self.innerSolver.setConvergenceTolerance(absTol=atol, relTol=rtol)

            # Before calling the inner solver we need to create a callback function so that we can store data in this solver's history file at every iteration of the inner solver
            def continuationcallBack(solver, u, res, monitorVars):
                monitorVars["SubIter"] = solver.iterationCount
                monitorVars["Increment"] = increment
                monitorVars["Lambda"] = currentLambda

                if self.rank == 0:
                    self.history.write(monitorVars)

                if self.userCallback is not None:
                    self.userCallback(solver, u, res, monitorVars)

            self.innerSolver.setCallback(continuationcallBack)

            self.innerSolver.setRefNorm(self.refNorm * currentLambda)
            self.innerSolver.solve()
            success = self.innerSolver.hasConverged
            numIters = self.innerSolver.iterationCount
            self._iterationCount += numIters

            # --- Check convergence ---
            isLastIncrement = increment == MAX_INCREMENTS - 1
            if not success:
                # If this increment failed then we should probably erase any saved states used for the predictor steps because they're clearly not working well
                if USE_PREDICTOR:
                    for ii in range(len(self.equilibriumPathLoadScales)):
                        self.equilibriumPathLoadScales[ii] = None
                        self.equilibriumPathStates[ii].zeroEntries()
                # If this was the first increment restarting from a previous solution then we don't have a safe state to reset to, so we just have to do a full reset
                if isRestartIncrement:
                    self.stateVec.zeroEntries()
                    self.incStartState.zeroEntries()
                    self.setStateFunc(self.stateVec)
                    stepSize = INIT_STEP
                    currentLambda = 0.0
                    loadStepDirection = 1
                # If the inner solve failed then we'll reduce the step size and try again, unless we've hit the increment or step size limits
                elif not isLastIncrement and stepSize > MIN_STEP:
                    self.setStateFunc(self.incStartState)
                    currentLambda -= stepSize * loadStepDirection
                    self.setLambdaFunc(currentLambda)
                    stepSize *= STEP_RETRACT_FACTOR
                else:
                    self._fatalFailure = self.innerSolver.fatalFailure
                    break
            else:
                # If inner solver converged and we're at the max load scale then we're done
                if currentLambda == MAX_LAMBDA:
                    self._hasConverged = True
                    break
                else:
                    if numIters != 0:
                        stepChangeFactor = np.sqrt(TARGET_ITERS / numIters)
                        stepSize *= np.clip(
                            stepChangeFactor, MIN_STEP_FACTOR, MAX_STEP_FACTOR
                        )
                    if USE_PREDICTOR:
                        stateToOverwrite = self.equilibriumPathStates.pop(0)
                        stateToOverwrite.copyValues(self.stateVec)
                        self.equilibriumPathStates.append(stateToOverwrite)

                        self.equilibriumPathLoadScales.pop(0)
                        self.equilibriumPathLoadScales.append(currentLambda)

            maxStep = min(np.abs(MAX_LAMBDA - currentLambda), MAX_STEP)
            stepSize = np.clip(stepSize, MIN_STEP, maxStep)
            currentLambda += loadStepDirection * stepSize
            self.setLambdaFunc(currentLambda)

            isRestartIncrement = False

        if result is not None:
            result.copyValues(self.stateVec)

        return

    def computeForceVectors(self) -> None:
        """Compute the current forcing vector

        The continuation solver is based on the assumption that the residual takes the following form:

        r(u, lambda) = F_int(u) + lambda * F_ext(u, lambda)

        This function computes fInt and fExt using two residual evaluations at lambda=0 and lambda=1:
        f_int = r(u, 0)
        f_ext = r(u, 1) - f_int
        """
        currentLambda = self.getLambdaFunc()
        self.setLambdaFunc(0.0)
        self.resFunc(self.fInt)
        self.setLambdaFunc(1.0)
        self.resFunc(self.fExt)
        self.fExt.axpy(-1.0, self.fInt)
        self.setLambdaFunc(currentLambda)
        return
