"""
==============================================================================
TACS: Base Solver Class
==============================================================================
@Description : Base class for pyTACS nonlinear solvers
"""

# ==============================================================================
# Standard Python modules
# ==============================================================================
import abc
from typing import Optional, Callable

# ==============================================================================
# External Python modules
# ==============================================================================
import mpi4py

# ==============================================================================
# Extension modules
# ==============================================================================
import tacs.TACS
from tacs.utilities import BaseUI


class BaseSolver(BaseUI):
    def __init__(
        self,
        assembler: tacs.TACS.Assembler,
        setStateFunc: Callable,
        resFunc: Callable,
        stateVec: Optional[tacs.TACS.Vec] = None,
        resVec: Optional[tacs.TACS.Vec] = None,
        options: Optional[dict] = None,
        comm: Optional[mpi4py.MPI.Comm] = None,
    ) -> None:
        """Create a solver instance

        Parameters
        ----------
        assembler : tacs.TACS.Assembler
            TACS assembler object related to the problem being solved, required in order for the solver to create it's own vectors
        setStateFunc : function
            Function to set the state vector, with signature setStateFunc(stateVec: tacs.TACS.Vec) -> None
        resFunc : function
            Function to evaluate the residual at the current state, with signature resFunc(resVec: tacs.TACS.Vec) -> None
        stateVec : tacs.TACS.Vec, optional
            Vector to store the state in, by default the solver will create it's own but these can be passed to save additional allocations
        resVec : tacs.TACS.Vec, optional
            Vector to store the residual in, by default the solver will create it's own but these can be passed to save additional allocations
        options : dict, optional
            Dictionary holding solver-specific option parameters (case-insensitive)., by default None
        comm : mpi4py.MPI.Intracomm, optional
            The comm object on which to create the pyTACS object., by default MPI.COMM_WORLD
        """
        self.assembler = assembler
        self.setStateFunc = setStateFunc
        self.resFunc = resFunc
        self.stateVec = stateVec if stateVec is not None else self.assembler.createVec()
        self.resVec = resVec if resVec is not None else self.assembler.createVec()
        self.refNorm = 1.0

        self._hasConverged = False
        self._fatalFailure = False
        self._iterationCount = 0
        self.callback = None
        BaseUI.__init__(self, options, comm)

    @property
    def hasConverged(self) -> bool:
        """Whether the solver has converged, set as a property rather than an attribute so that it is read-only"""
        return self._hasConverged

    @property
    def fatalFailure(self) -> bool:
        """Whether the solver has failed, set as a property rather than an attribute so that it is read-only

        Note that a fatalFailure is not the same as not converging, this flag is meant to reflect that there has been a fatal failure in the solver which requires a full reset
        """
        return self._fatalFailure

    @property
    def iterationCount(self) -> int:
        """Number of iterations performed, set as a property rather than an attribute so that it is read-only"""
        return self._iterationCount

    @abc.abstractmethod
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

    @abc.abstractmethod
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

    def initializeSolve(self) -> None:
        """Perform any initialization required before the solve

        In the base solver class, this simply involves resetting the iteration counter and convergence flags
        """
        self._iterationCount = 0
        self._hasConverged = False
        self._fatalFailure = False

    def reset(self) -> None:
        """Reset the solver

        Currently this just zeros out the state vector, but more functionality may be added in future
        """
        self.stateVec.zeroEntries()
        self.setStateFunc(self.stateVec)

    def setRefNorm(self, norm: float) -> None:
        """Set the reference norm used to compute relative convergence measures

        Parameters
        ----------
        norm : float
            Reference norm
        """
        self.refNorm = norm

    def setCallback(self, callback: Callable) -> None:
        """Set a callback function to be called at each iteration

        Parameters
        ----------
        callback : callable, optional
            Callback function, should have the following signature:
            callback(solver: tacs.solver, u: tacs.TACS.Vec, res: tacs.TACS.Vec, monitorVars: dict) -> None
            Where:
                - ``solver`` is the solver object
                - ``u`` is the current state vector
                - ``res`` is the current residual vector
                - ``monitorVars`` is a dictionary of variables to monitor, which can be specified through the ``"nonlinearSolverMonitorVars"`` option
        """
        self.callback = callback
