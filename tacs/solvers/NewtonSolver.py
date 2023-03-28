"""
==============================================================================
TACS Nonlinear Newton Solver
==============================================================================
@Author : Alasdair Christison Gray
@Description : A Newton solver for nonlinear TACS problems
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


class NewtonSolver(BaseSolver):
    defaultOptions = {
        "newtonSolverMaxIter": [int, 40, "Maximum number of Newton iterations."],
        "newtonSolverForceFirstIter": [
            bool,
            False,
            "Force the solver to perform the first Newton iteration, even if the convergence criteria are satisfied at the initial point.",
        ],
        "newtonSolverAbsTol": [
            float,
            1e-8,
            "Convergence criteria for the nonlinear residual norm.",
        ],
        "newtonSolverRelTol": [
            float,
            1e-8,
            "Relative convergence criteria for the nonlinear residual norm, norm is measured relative to that of the external load vector.",
        ],
        "newtonSolverDivergenceTol": [
            float,
            1e10,
            "Residual norm at which the nonlinear solver is jugded to have diverged",
        ],
        "newtonSolverAbsLinTol": [float, 1e-12, "Linear solver residual tolerance."],
        "newtonSolverRelLinTol": [
            float,
            1e-12,
            "Linear solver relative residual tolerance.",
        ],
        "newtonSolverMaxLinIters": [
            int,
            0,
            "If the linear solver takes more than this number of iterations to converge, the preconditioner is updated.",
        ],
        "newtonSolverUseEW": [
            bool,
            False,
            "Flag for enabling use of variable linear solver convergence using the Eisenstat-Walker method.",
        ],
        "newtonSolverEWMaxTol": [
            float,
            0.01,
            "Eisenstat-Walker max allowable linear solver tolerance.",
        ],
        "newtonSolverEWGamma": [float, 1.0, "Eisenstat-Walker gamma parameter."],
        "newtonSolverEWAlpha": [
            float,
            0.5 * (1.0 + np.sqrt(5)),
            "Eisenstat-Walker alpha parameter.",
        ],
        # Line search options
        "useLineSearch": [
            bool,
            True,
            "Flag for using line search in the nonlinear solver.",
        ],
        "lineSearchMonitor": [
            bool,
            False,
            "Flag for printing out line search information.",
        ],
        "skipFirstNLineSearch": [
            int,
            0,
            "Skip the first N line searches. Setting this to 1 can improve the convergence speed of Newton solver, but also decreases robustness",
        ],
        "lineSearchMaxIter": [int, 25, "Maximum number of linesearch iterations."],
        "lineSearchExpectedDecrease": [
            float,
            1e-4,
            "Minimum fraction of the expected decrease in the energy gradient during the linesearch. Should be between 0 and 1. Higher values should improve robustness at the expense of solution time.",
        ],
        "lineSearchMaxStep": [
            float,
            2.0,
            "Maximum step size for the linesearch, as a fraction of the Newton step",
        ],
        "lineSearchMinStep": [
            float,
            1e-2,
            "Minimum step size for the linesearch, as a fraction of the Newton step",
        ],
        "lineSearchMaxStepChange": [
            float,
            0.5,
            "Maximum change in the step size from one linesearch iteration to the next, can be useful in cases where secant method bounces between upper and lower step bounds.",
        ],
        "lineSearchFallbackStepLimit": [
            float,
            0.9,
            "Often, the value of the merit function at the Newton step (alpha = 1.0), is orders of magnitude greater than at the start point. In these situations, the linesearch then tries to evaluate a point with a very small step size, which usually meets the expected decrease criteria but results in very slow progress of the Newton solver. To combat this, this value limits how far the linesearch can backtrack on the first iteration after evaluating alpha = 1. This has the effect of encouraging the linesearch to find larger steps that meet the expected decrease criterion, which results in faster convergence of the Newton solver.",
        ],
    }

    def __init__(
        self,
        assembler: tacs.TACS.Assembler,
        setStateFunc: Callable,
        resFunc: Callable,
        jacFunc: Callable,
        pcUpdateFunc: Callable,
        linearSolver: tacs.TACS.KSM,
        stateVec: Optional[tacs.TACS.Vec] = None,
        resVec: Optional[tacs.TACS.Vec] = None,
        options: Optional[dict] = None,
        comm: Optional[mpi4py.MPI.Comm] = None,
    ) -> None:
        """Create a Newton solver instance

        Parameters
        ----------
        assembler : tacs.TACS.Assembler
            TACS assembler object related to the problem being solved, required in order for the solver to create it's own vectors
        setStateFunc : function
            Function to set the state vector, with signature setStateFunc(stateVec: tacs.TACS.Vec) -> None
        resFunc : function
            Function to evaluate the residual at the current state, with signature resFunc(resVec: tacs.TACS.Vec) -> None
        jacFunc : function
            Function to update the residual Jacobian at the current state, with signature jacFunc() -> None
        pcUpdateFunc : function
            Function to update the residual Jacobian preconditioner at the current state, with signature pcUpdateFunc() -> None
        linearSolver : tacs.TACS.KSM
            TACS linear solver object to use for the Newton solve, the linear solver owns the matrix and preconditioner
        stateVec : tacs.TACS.Vec, optional
            Vector to store the state in, by default the solver will create it's own but these can be passed to save additional allocations
        resVec : tacs.TACS.Vec, optional
            Vector to store the residual in, by default the solver will create it's own but these can be passed to save additional allocations
        options : dict, optional
            Dictionary holding solver-specific option parameters (case-insensitive)., by default None
        comm : mpi4py.MPI.Intracomm, optional
            The comm object on which to create the pyTACS object., by default MPI.COMM_WORLD
        """
        BaseSolver.__init__(
            self,
            assembler=assembler,
            setStateFunc=setStateFunc,
            resFunc=resFunc,
            stateVec=stateVec,
            resVec=resVec,
            options=options,
            comm=comm,
        )
        self.jacFunc = jacFunc
        self.pcUpdateFunc = pcUpdateFunc
        self.linearSolver = linearSolver

        # Create additional vectors
        self.update = self.assembler.createVec()

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
            self.setOption("newtonSolverAbsTol", absTol)
        if relTol is not None:
            self.setOption("newtonSolverRelTol", relTol)
        return

    def solve(
        self, u0: Optional[tacs.TACS.Vec] = None, result: Optional[tacs.TACS.Vec] = None
    ) -> None:
        """Solve the residual equations r(u) = 0

        Parameters
        ----------
        u0 : TACS vector, optional
            Initial guess, by default uses the current state
        result : TACS vector, optional
            Vector in which to store the solution, by default None.
            The problem's state is updated with the solution whether or not this is provided.
        """
        USE_LINESEARCH = self.getOption("useLineSearch")
        LINESEARCH_SKIP_ITERS = self.getOption("skipFirstNLineSearch")
        MAX_ITERS = self.getOption("newtonSolverMaxIter")
        MAX_RES = self.getOption("newtonSolverDivergenceTol")
        MAX_LIN_ITERS = self.getOption("newtonSolverMaxLinIters")
        FORCE_FIRST_ITER = self.getOption("newtonSolverForceFirstIter")

        # Linear solver convergence options
        USE_EW = self.getOption("newtonSolverUseEW")
        LIN_SOLVE_TOL_MAX = self.getOption("newtonSolverEWMaxTol")
        LIN_SOLVE_TOL_MIN = self.getOption("newtonSolverRelLinTol")
        EW_ALPHA = self.getOption("newtonSolverEWAlpha")
        EW_GAMMA = self.getOption("newtonSolverEWGamma")
        linCovergenceRel = LIN_SOLVE_TOL_MAX if USE_EW else LIN_SOLVE_TOL_MIN

        ABS_TOL = self.getOption("newtonSolverAbsTol")
        REL_TOL = self.getOption("newtonSolverRelTol")

        self.initializeSolve()

        flags = ""

        if u0 is not None:
            self.stateVec.copyValues(u0)

        for iteration in range(MAX_ITERS):
            self._iterationCount = iteration
            prevLinCovergenceRel = linCovergenceRel

            # Compute residual and norms
            self.resFunc(self.resVec)
            if iteration > 0:
                prevResNorm = resNorm
            resNorm = self.resVec.norm()
            uNorm = self.stateVec.norm()

            # Test convergence
            if not FORCE_FIRST_ITER or iteration > 0:
                self._hasConverged = (
                    np.real(resNorm) / np.real(self.refNorm) < REL_TOL
                    or np.real(resNorm) < ABS_TOL
                )
                self._fatalFailure = np.real(resNorm) >= MAX_RES or np.isnan(resNorm)
                if self._hasConverged:
                    flags += "C"
                elif self.fatalFailure:
                    flags += "D"

            # Write data to history
            monitorVars = {
                "SubIter": iteration,
                "Res norm": resNorm,
                "Rel res norm": resNorm / self.refNorm,
                "U norm": uNorm,
                "Flags": flags,
            }
            if iteration > 0:
                monitorVars["Lin iters"] = linearSolveIterations
                monitorVars["Lin res"] = np.abs(linearSolveResNorm)
                monitorVars["LS step"] = alpha
                monitorVars["LS iters"] = lineSearchIters
                if USE_EW:
                    monitorVars["EW Tol"] = prevLinCovergenceRel
            if self.callback is not None:
                self.callback(self, self.stateVec, self.resVec, monitorVars)

            flags = ""

            # exit if converged/diverged
            if self.hasConverged or self.fatalFailure:
                break

            # Update Jacobian
            self.jacFunc()

            # Update preconditioner, or skip if last linear solve converged in few enough iterations
            if iteration > 0 and linearSolveIterations <= MAX_LIN_ITERS:
                pass
            else:
                flags += "P"
                self.pcUpdateFunc()

            # Update linear solver convergence tolerance using Einstat-Walker method b)
            if USE_EW:
                if iteration > 0:
                    zeta = EW_GAMMA * np.real(resNorm / prevResNorm) ** EW_ALPHA
                    threshold = EW_GAMMA * prevLinCovergenceRel**EW_ALPHA
                    if threshold <= 0.1:
                        linCovergenceRel = zeta
                    else:
                        linCovergenceRel = max(zeta, threshold)
                linCovergenceRel = np.clip(
                    linCovergenceRel, LIN_SOLVE_TOL_MIN, LIN_SOLVE_TOL_MAX
                )
                self.linearSolver.setTolerances(
                    float(linCovergenceRel), self.getOption("newtonSolverAbsLinTol")
                )

            # Compute Newton step
            linSolveConverged = self.linearSolver.solve(self.resVec, self.update)
            linSolveConverged = linSolveConverged == 1

            # If we didn't converge because the preconditioner isn't up to date, update preconditioner and try again
            if not linSolveConverged and "P" not in flags:
                self.pcUpdateFunc()
                flags += "P"
                self.linearSolver.solve(self.resVec, self.update)

            self.update.scale(-1.0)

            # Check data from linear solve
            linearSolveIterations = self.linearSolver.getIterCount()
            linearSolveResNorm = self.linearSolver.getResidualNorm()

            if USE_LINESEARCH and iteration >= LINESEARCH_SKIP_ITERS:
                # Do linesearch
                alpha, lineSearchIters = self.energyLineSearch(
                    self.stateVec, self.update
                )
            else:
                alpha = 1.0
                lineSearchIters = 1
            self.stateVec.axpy(alpha, self.update)
            self.setStateFunc(self.stateVec)

        if result is not None:
            result.copyValues(self.stateVec)

    def energyLineSearch(self, u, stepDir, slope=None):
        MAX_LINESEARCH_ITERS = self.getOption("lineSearchMaxIter")
        LINESEARCH_MU = self.getOption("lineSearchExpectedDecrease")
        LINESEARCH_ALPHA_MIN = self.getOption("lineSearchMinStep")
        LINESEARCH_ALPHA_MAX = self.getOption("lineSearchMaxStep")
        LINESEARCH_MAX_STEP_CHANGE = self.getOption("lineSearchMaxStepChange")
        PRINT_LINESEARCH_ITERS = self.getOption("lineSearchMonitor")
        if slope is None:
            slope = 1.0

        # Compute residual and merit function at u0
        self.setStateFunc(u)
        self.resFunc(self.resVec)
        f0 = np.real(self.resVec.dot(stepDir))
        fOld = f0
        alphaOld = 0.0
        uNorm = u.norm()
        if self.rank == 0 and PRINT_LINESEARCH_ITERS:
            print(
                f"Line search iter  0: alpha = {0: 11e},   f0 = {(f0): 11e}, uNorm = {uNorm: 11e}"
            )

        # 3. Set $\alpha = 1$
        alpha = 1.0
        alphaNew = alpha
        for iteration in range(MAX_LINESEARCH_ITERS):
            # 4. Increment state, $u = u + \alpha \Delta u$
            u.axpy(alpha, stepDir)
            self.setStateFunc(u)

            # 5. Compute residual, $r = r(u)$
            self.resFunc(self.resVec)

            # 6. Compute merit function,  $f(\alpha)=f(u, r, \Delta u)$
            fNew = np.real(self.resVec.dot(stepDir))

            # 7. if $abs(f(\alpha)) \leq \mu f_0 + \alpha f'_0$:
            #     exit
            uNorm = u.norm()
            if self.rank == 0 and PRINT_LINESEARCH_ITERS:
                print(
                    f"Line search iter {(iteration+1):2d}: alpha = {alpha: 11e}, f/f0 = {(fNew/f0): 11e}, uNorm = {uNorm: 11e}"
                )
            u.axpy(-alpha, stepDir)
            fReduction = np.abs(fNew / f0)
            if fReduction <= 1 - LINESEARCH_MU * min(alpha, 1.0) * slope:
                break
            else:
                # 8. Update $\alpha$ (based on search method)
                if iteration == 0:
                    alphaMin = 0.9
                else:
                    alphaMin = LINESEARCH_ALPHA_MIN
                if fNew == fOld:
                    alphaNew = alpha + LINESEARCH_ALPHA_MIN
                else:
                    alphaNew = np.clip(
                        alpha - fNew * (alpha - alphaOld) / (fNew - fOld),
                        alphaMin,
                        LINESEARCH_ALPHA_MAX,
                    )
                if iteration > 0 and abs(alphaNew - alpha) > LINESEARCH_MAX_STEP_CHANGE:
                    alphaNew = (
                        alpha + np.sign(alphaNew - alpha) * LINESEARCH_MAX_STEP_CHANGE
                    )
                alphaOld = alpha
                alpha = alphaNew
                fOld = fNew
            # 9. return to step 4
        return alpha, iteration + 1
