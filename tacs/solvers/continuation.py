"""
==============================================================================
TACS Nonlinear Continuation Solver
==============================================================================
@Author : Alasdair Christison Gray
@Description : A predictor-corrector force continuation solver for nonlinear TACS problems
"""

# ==============================================================================
# Standard Python modules
# ==============================================================================
from typing import Optional, Callable

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
        "continuationMaxLambda": [
            float,
            1.0,
            "Final continuation parameter value to aim for.",
        ],
        "continuationAbsTol": [
            float,
            1e-8,
            "Convergence criteria for the nonlinear residual norm.",
        ],
        "continuationRelTol": [
            float,
            1e-8,
            "Relative convergence criteria for the nonlinear residual norm, norm is measured relative to that of the external load vector.",
        ],
        "continuationCoarseAbsTol": [
            float,
            1e-4,
            "Residual norm criteria for intermediate continuation steps, making this larger may speed up the nonlinear solver by allowing it to only partially converge intermediate steps.",
        ],
        "continuationCoarseRelTol": [
            float,
            1e-4,
            "Relative residual norm criteria for intermediate load increments.",
        ],
        "continuationTargetIter": [
            int,
            8,
            "Target number of Newton iterations for each continuation increment.",
        ],
        "continuationMaxIter": [int, 30, "Maximum number of continuation steps."],
        "continuationInitialStep": [float, 0.2, "Initial continuation step size."],
        "continuationMinStep": [float, 1e-4, "Minimum continuation step size."],
        "continuationMaxStep": [float, np.inf, "Maximum continuation step size."],
        "continuationMinStepFactor": [
            float,
            0.5,
            "The minimum factor by which the continuation step size can decrease in a single step.",
        ],
        "continuationMaxStepFactor": [
            float,
            2.0,
            "The maximum factor by which the continuation step size can increase in a single step.",
        ],
        "continuationRetractionFactor": [
            float,
            0.5,
            "The factor by which the continuation step size is reduced when the Newton solver fails to converge.",
        ],
        # Predictor step options
        "continuationUsePredictor": [
            bool,
            False,
            "Flag for using predictor step in continuation.",
        ],
        "continuationNumPredictorStates": [
            int,
            2,
            "Number of previous equilibrium states to use in computing the predictor step.",
        ],
        "predictorUseDerivative": [
            bool,
            False,
            "Whether to use the equilibrium path slope in the computation of the predictor step. This requires a linear solve and thus greatly increases the cost of the predictor step computation.",
        ],
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
            Function to update the residual Jacobian at the current state, with signature jacFunc()
        pcUpdateFunc : function
            Function to update the residual Jacobian preconditioner at the current state, with signature pcUpdateFunc()
        linearSolver : tacs.TACS.KSM
            TACS linear solver object to use for the Newton solve, the linear solver owns the matrix and preconditioner
        setLambdaFunc : function
            Function to set the continuation parameter, with signature setLambdaFunc(lambda:float) -> None
        setLambdaFunc : function
            Function to get the current continuation parameter, with signature getLambdaFunc() -> float
        innerSolver : TACS Nonlinear Solver
            A Solver object to use for the corrector solve in each increment (e.g a NewtonSolver object)
        options : dict, optional
            Dictionary holding solver-specific option parameters (case-insensitive)., by default None
        comm : mpi4py.MPI.Intracomm, optional
            The comm object on which to create the pyTACS object., by default MPI.COMM_WORLD
        """

        self.jacFunc = jacFunc
        self.pcUpdateFunc = pcUpdateFunc
        self.linearSolver = linearSolver
        self.setLambdaFunc = setLambdaFunc
        self.getLambdaFunc = getLambdaFunc
        self.innerSolver = innerSolver
        self.defaultOptions.update(innerSolver.defaultOptions)

        self.equilibriumPathStates = None
        self.equilibriumPathLoadScales = None

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

    def _setupPredictorVectors(self) -> None:
        """Setup the structures containing the data for computing the predictor steps"""
        self.equilibriumPathStates = []
        self.equilibriumPathLoadScales = []
        if self.getOption("continuationUsePredictor"):
            for _ in range(self.getOption("continuationNumPredictorStates")):
                self.equilibriumPathStates.append(self.assembler.createVec())
                self.equilibriumPathLoadScales.append(None)

    def setOption(self, name, value) -> None:
        BaseSolver.setOption(self, name, value)

        # Pass option to inner solver if it is a inner solver option
        if name.lower() in [opt.lower() for opt in self.innerSolver.defaultOptions]:
            self.innerSolver.setOption(name, value)

        # Update the predictor computation data structures if the relevant options are changed
        if name.lower() in [
            "continuationusepredictor",
            "continuationnumpredictorstates",
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
            self.setOption("continuationAbsTol", absTol)
        if relTol is not None:
            self.setOption("continuationRelTol", relTol)

        return

    def initializeSolve(self) -> None:
        """Perform any initialization required before the solve"""
        BaseSolver.initializeSolve(self)
        if self.getOption("continuationUsePredictor"):
            for ii in range(self.getOption("continuationNumPredictorStates")):
                self.equilibriumPathLoadScales[ii] = None
                self.equilibriumPathStates[ii].zeroEntries()

    def solve(
        self, u0: Optional[tacs.TACS.Vec] = None, result: Optional[tacs.TACS.Vec] = None
    ) -> None:
        MAX_LAMBDA = self.getOption("continuationMaxLambda")
        TARGET_ITERS = self.getOption("continuationTargetIter")
        INIT_STEP = self.getOption("continuationInitialStep")
        MIN_STEP = self.getOption("continuationMinStep")
        MAX_STEP = self.getOption("continuationMaxStep")
        MAX_INCREMENTS = self.getOption("continuationMaxIter")
        MIN_STEP_FACTOR = self.getOption("continuationMinStepFactor")
        MAX_STEP_FACTOR = self.getOption("continuationMaxStepFactor")
        STEP_RETRACT_FACTOR = self.getOption("continuationRetractionFactor")

        ABS_TOL = self.getOption("continuationAbsTol")
        REL_TOL = self.getOption("continuationRelTol")
        COARSE_ABS_TOL = self.getOption("continuationCoarseAbsTol")
        COARSE_REL_TOL = self.getOption("continuationCoarseRelTol")

        USE_PREDICTOR = self.getOption("continuationUsePredictor")
        # PREDICTOR_USE_DERIVATIVE = self.getOption("predictorUseDerivative")

        self.initializeSolve()

        if u0 is not None:
            self.stateVec.copyValues(u0)

        # Compute the internal and external forcing vectors at the current point
        self.computeForceVectors()
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
        if np.real(self.stateVec.norm()) > 0:
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

            # Before calling the inner solver we need to create a wrapper around the user's callback function so that we can add some continuation information to the monitor variables
            if self.callback is not None:

                def continuationcallBack(solver, u, res, monitorVars):
                    monitorVars["SubIter"] = solver.iterationCount
                    monitorVars["Increment"] = increment
                    monitorVars["Lambda"] = currentLambda
                    self.callback(solver, u, res, monitorVars)

                self.innerSolver.setCallback(continuationcallBack)

            self.innerSolver.setRefNorm(self.refNorm * currentLambda)
            self.innerSolver.solve()
            success = self.innerSolver.hasConverged
            numIters = self.innerSolver.iterationCount
            self._iterationCount += numIters

            # --- Check convergence ---
            isLastIncrement = increment == MAX_INCREMENTS - 1
            if not success:
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
