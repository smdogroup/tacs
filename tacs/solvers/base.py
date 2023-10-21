"""
==============================================================================
TACS: Base Solver Class
==============================================================================
This is the base class from which all other TACS solvers are derived.
It is an abstract class, so cannot be used directly. Instead, it defines
the methods that all solvers must implement.
"""

# ==============================================================================
# Standard Python modules
# ==============================================================================
import abc
from typing import Optional, Callable, Dict

# ==============================================================================
# External Python modules
# ==============================================================================
import mpi4py

# ==============================================================================
# Extension modules
# ==============================================================================
import tacs.TACS
from tacs.utilities import BaseUI, SolverHistory


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
            Function to evaluate the residual at the current state, with signature `resFunc(resVec: tacs.TACS.Vec) -> None`
        stateVec : tacs.TACS.Vec, optional
            Vector to store the state in, by default the solver will create it's own but these can be passed to save additional allocations
        resVec : tacs.TACS.Vec, optional
            Vector to store the residual in, by default the solver will create it's own but these can be passed to save additional allocations
        options : dict, optional
            Dictionary holding solver-specific option parameters (case-insensitive)., by default None
        comm : mpi4py.MPI.Intracomm, optional
            The comm object on which to create the pyTACS object., by default mpi4py.MPI.COMM_WORLD
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
        BaseUI.__init__(self, options, comm)

        # No solver history for now
        self.history = None
        self.userCallback = None

    def _createSolverHistory(self):
        """Setup the solver history object based on the current options

        The solver history is only created on the root processor.
        """
        histVars = self.getHistoryVariables()
        if self.comm.rank == 0:
            self.history = SolverHistory()
            for varName, variable in histVars.items():
                self.history.addVariable(
                    varName, variable["type"], printVar=variable["print"]
                )

    @abc.abstractmethod
    def getHistoryVariables(self) -> Dict[str, Dict]:
        """Get the variables to be stored in the solver history

        This method allows for implementation of any logic that dictates any changes in the stored variables depending on the current options.

        Returns
        -------
        Dict[str, Dict]
            Dictionary of solver variables, keys are the variable names, value is another dictionary with keys "type" and "print", where "type" is the data type of the variable and "print" is a boolean indicating whether or not to print the variable to the screen
        """

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

    def initializeSolve(self, u0: Optional[tacs.TACS.Vec] = None) -> None:
        """Perform any initialization required before the solve

        Parameters
        ----------
        u0 : TACS vector, optional
            Initial guess, by default uses the current state
        """
        self._iterationCount = 0
        self._hasConverged = False
        self._fatalFailure = False

        # Initialize the state vector and set the boundary conditions
        if u0 is not None:
            self.stateVec.copyValues(u0)
        self.assembler.setBCs(self.stateVec)
        self.setStateFunc(self.stateVec)

        # Reset the solver history and store the solver options as metadata
        if self.rank == 0:
            if self.history is None:
                self._createSolverHistory()
            else:
                self.history.reset()
            if "options" not in self.history.getMetadata():
                self.history.addMetadata("options", self.options)

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
        self.userCallback = callback
