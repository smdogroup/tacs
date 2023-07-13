"""
The main purpose of this class is to represent all relevant
information for a static analysis. This will include
information defining the loading condition as well as various
other pieces of information.

.. note:: This class should be created using the
    :meth:`pyTACS.createStaticProblem <tacs.pytacs.pyTACS.createStaticProblem>` method.
"""

# =============================================================================
# Imports
# =============================================================================
import copy
import os
import time
from collections import OrderedDict

import numpy as np
import pyNastran.bdf as pn

import tacs.TACS
import tacs.elements
import tacs.solvers
from tacs.problems.base import TACSProblem
from tacs.utilities import SolverHistory


class StaticProblem(TACSProblem):
    # Default options for class
    defaultOptions = {
        "outputDir": [str, "./", "Output directory for F5 file writer."],
        # Solution Options
        "KSMSolver": [
            str,
            "GMRES",
            "Krylov subspace method to use for linear solver. Currently only supports 'GMRES'",
        ],
        "orderingType": [
            int,
            tacs.TACS.ND_ORDER,
            "Ordering type to use for matrix partitioning.\n"
            "\t Acceptable values are:\n"
            f"\t\t tacs.TACS.NATURAL_ORDER = {tacs.TACS.NATURAL_ORDER}\n"
            f"\t\t tacs.TACS.RCM_ORDER = {tacs.TACS.RCM_ORDER}\n"
            f"\t\t tacs.TACS.ND_ORDER = {tacs.TACS.ND_ORDER}\n"
            f"\t\t tacs.TACS.TACS_AMD_ORDER = {tacs.TACS.TACS_AMD_ORDER}\n"
            f"\t\t tacs.TACS.MULTICOLOR_ORDER = {tacs.TACS.MULTICOLOR_ORDER}",
        ],
        "PCFillLevel": [int, 1000, "Preconditioner fill level."],
        "PCFillRatio": [float, 20.0, "Preconditioner fill ratio."],
        "subSpaceSize": [int, 10, "Subspace size for Krylov solver."],
        "nRestarts": [int, 15, "Max number of restarts for Krylov solver."],
        "flexible": [
            bool,
            True,
            "Flag for whether the preconditioner is flexible.",
        ],
        "L2Convergence": [
            float,
            1e-12,
            "Absolute convergence tolerance for linear solver based on l2 norm of residual.",
        ],
        "L2ConvergenceRel": [
            float,
            1e-12,
            "Relative convergence tolerance for linear solver based on l2 norm of residual.",
        ],
        "RBEStiffnessScaleFactor": [
            float,
            1e3,
            "Constraint matrix scaling factor used in RBE Lagrange multiplier stiffness matrix.",
        ],
        "RBEArtificialStiffness": [
            float,
            1e-3,
            "Artificial constant added to diagonals of RBE Lagrange multiplier stiffness matrix \n"
            "\t to stabilize preconditioner.",
        ],
        "useMonitor": [
            bool,
            False,
            "Flag for whether to attach a debug monitor to the linear solver.",
        ],
        "monitorFrequency": [
            int,
            10,
            "Print frequency for sub iterations of linear solver.",
        ],
        "nonlinearSolverMonitorVars": [
            list,
            [
                "linSolverIters",
                "linSolverRes",
                "lineSearchStep",
                "lineSearchIters",
            ],
            "List of variables to include in nonlinear solver monitor output. Choose from 'linSolverIters', 'linSolverRes', 'loadScale', 'lineSearchStep', 'EWTol', and 'lineSearchIters'.",
        ],
        # Output Options
        "writeSolution": [
            bool,
            True,
            "Flag for suppressing all f5 file writing.",
        ],
        "writeNLIterSolutions": [
            bool,
            False,
            "Flag to save a solution file at every nonlinear solver iterations.",
        ],
        "numberSolutions": [
            bool,
            True,
            "Flag for attaching solution counter index to f5 files.",
        ],
        "printTiming": [
            bool,
            False,
            "Flag for printing out timing information for class procedures.",
        ],
    }
    defaultOptions.update(tacs.solvers.NewtonSolver.defaultOptions)
    defaultOptions.update(tacs.solvers.ContinuationSolver.defaultOptions)

    def __init__(
        self,
        name,
        assembler,
        comm,
        outputViewer=None,
        meshLoader=None,
        isNonlinear=False,
        options=None,
    ):
        """
        NOTE: This class should not be initialized directly by the user.
        Use pyTACS.createStaticProblem instead.

        Parameters
        ----------
        name : str
            Name of this tacs problem

        assembler : tacs.TACS.Assembler
            Cython object responsible for creating and setting tacs objects used to solve problem

        comm : mpi4py.MPI.Intracomm
            The comm object on which to create the pyTACS object.

        outputViewer : tacs.TACS.TACSToFH5
            Cython object used to write out f5 files that can be converted and used for postprocessing.

        meshLoader : tacs.pymeshloader.pyMeshLoader
            pyMeshLoader object used to create the assembler.

        options : dict
            Dictionary holding problem-specific option parameters (case-insensitive).

        """

        # Problem name
        self.name = name

        # Set linear solver to None, until we set it up later
        self.KSM = None
        self.history = None
        self.newtonSolver = None
        self.nonlinearSolver = None

        # Default setup for common problem class objects, sets up comm and options
        TACSProblem.__init__(
            self, assembler, comm, options, outputViewer, meshLoader, isNonlinear
        )

        # Create problem-specific variables
        self._createVariables()

        # Setup solver and solver history objects for nonlinear problems
        if self.isNonlinear:
            self._createSolverHistory()

            # Create Newton solver, the inner solver for the continuation solver
            solverOptions = {}
            solverOptionNames = [
                name.lower() for name in tacs.solvers.NewtonSolver.defaultOptions
            ]
            for key in self.options:
                if key.lower() in solverOptionNames:
                    solverOptions[key] = self.getOption(key)

            self.newtonSolver = tacs.solvers.NewtonSolver(
                assembler=self.assembler,
                setStateFunc=self.setVariables,
                resFunc=self.getResidual,
                jacFunc=self.updateJacobian,
                pcUpdateFunc=self.updatePreconditioner,
                linearSolver=self.KSM,
                stateVec=self.u,
                resVec=self.res,
                options=solverOptions,
                comm=self.comm,
            )

            # And now create the continuation solver
            solverOptions = {}
            solverOptionNames = [
                name.lower() for name in tacs.solvers.ContinuationSolver.defaultOptions
            ]
            for key in self.options:
                if key.lower() in solverOptionNames:
                    solverOptions[key] = self.getOption(key)

            def getLoadScale():
                return self.loadScale

            self.nonlinearSolver = tacs.solvers.ContinuationSolver(
                jacFunc=self.updateJacobian,
                pcUpdateFunc=self.updatePreconditioner,
                linearSolver=self.KSM,
                setLambdaFunc=self.setLoadScale,
                getLambdaFunc=getLoadScale,
                innerSolver=self.newtonSolver,
                options=solverOptions,
                comm=self.comm,
            )
            self.nonlinearSolver.setCallback(self._nonlinearCallback)

    def _createSolverHistory(self):
        """Setup the solver history object based on the current options

        The solver history is only created on the root processor.
        """
        monitorVars = [s.lower() for s in self.getOption("nonlinearSolverMonitorVars")]
        numType = float if self.dtype == np.float64 else complex
        if self.comm.rank == 0:
            history = SolverHistory()

            # Define the variables to be stored in the history
            # Continuation increment
            history.addVariable("Increment", int, printVar=True)
            # Load scale
            history.addVariable("Lambda", float, printVar=True)
            # Newton solve iteration number
            history.addVariable("SubIter", int, printVar=True)
            # Einstat walker linear solver tolerance
            if self.getOption("newtonSolverUseEW"):
                history.addVariable("EW Tol", float, printVar="ewtol" in monitorVars)
            # Number of linear solver iterations
            history.addVariable(
                "Lin iters", int, printVar="linsolveriters" in monitorVars
            )
            # Linear solver residual norm
            history.addVariable(
                "Lin res", numType, printVar="linsolverres" in monitorVars
            )
            # Residual norm (absolute and relative)
            history.addVariable("Res norm", numType, printVar=True)
            history.addVariable("Rel res norm", numType, printVar=True)
            # state norm
            history.addVariable("U norm", numType, printVar=True)
            # Line search step size
            history.addVariable(
                "LS step",
                float,
                printVar="linesearchstep" in monitorVars
                and self.getOption("useLineSearch"),
            )
            # Num line search iterations
            history.addVariable(
                "LS iters",
                int,
                printVar="linesearchiters" in monitorVars
                and self.getOption("useLineSearch"),
            )
            # Flags
            history.addVariable("Flags", str, printVar=True)

            self.history = history

    def _createVariables(self):
        """Internal to create the variable required by TACS"""

        opt = self.getOption

        # Generic residual vector
        self.res = self.assembler.createVec()
        self.rhs = self.assembler.createVec()

        # Dictionaries to hold adjoint/sens vectors for each evalFunc
        self.adjointList = OrderedDict()
        self.dIduList = OrderedDict()
        self.dvSensList = OrderedDict()
        self.xptSensList = OrderedDict()

        # Temporary vector for adjoint solve
        self.phi = self.assembler.createVec()
        self.adjRHS = self.assembler.createVec()

        # Load vector
        self.F = self.assembler.createVec()
        self.F_array = self.F.getArray()

        # State variable vector
        self.u = self.assembler.createVec()
        self.u_array = self.u.getArray()

        # Vectors used to decompose residual into external and internal forces
        self.externalForce = self.assembler.createVec()
        self.internalForce = self.assembler.createVec()

        if self.isNonlinear:
            self.u_inc_start = self.assembler.createVec()

        # Auxiliary element object for applying tractions/pressure
        self.auxElems = tacs.TACS.AuxElements()
        # Counter for number of calls to `solve` method
        self.callCounter = -1

        # Norms
        self.initNorm = 0.0
        self.startNorm = 0.0
        self.finalNorm = 0.0

        # Load scaling factor
        self._loadScale = 1.0

        # Tangent Stiffness --- process the ordering option here:
        ordering = opt("orderingType")

        # True stiffness matrix
        self.K = self.assembler.createSchurMat(ordering)
        # Artificial stiffness for RBE numerical stabilization to stabilize PC
        self.rbeArtificialStiffness = self.assembler.createSchurMat(ordering)

        # Additional Vecs for updates
        self.update = self.assembler.createVec()

        # Setup PCScMat and KSM solver
        self.alpha = 1.0
        self.beta = 0.0
        self.gamma = 0.0

        # Computes stiffness matrix w/o art. terms
        # Set artificial stiffness factors in rbe class to zero
        tacs.elements.RBE2.setScalingParameters(opt("RBEStiffnessScaleFactor"), 0.0)
        tacs.elements.RBE3.setScalingParameters(opt("RBEStiffnessScaleFactor"), 0.0)
        self.assembler.assembleJacobian(
            self.alpha,
            self.beta,
            self.gamma,
            self.res,
            self.K,
            loadScale=self.loadScale,
        )

        # Now isolate art. terms
        # Recompute stiffness with artificial terms included
        tacs.elements.RBE2.setScalingParameters(
            opt("RBEStiffnessScaleFactor"), opt("RBEArtificialStiffness")
        )
        tacs.elements.RBE3.setScalingParameters(
            opt("RBEStiffnessScaleFactor"), opt("RBEArtificialStiffness")
        )
        self.assembler.assembleJacobian(
            self.alpha, self.beta, self.gamma, None, self.rbeArtificialStiffness
        )
        # Subtract full stiffness w/o artificial terms from full stiffness w/ terms
        # to isolate  artificial stiffness terms
        self.rbeArtificialStiffness.axpy(-1.0, self.K)

        reorderSchur = 1
        self.PC = tacs.TACS.Pc(
            self.K,
            lev_fill=opt("PCFillLevel"),
            ratio_fill=opt("PCFillRatio"),
            reorder=reorderSchur,
        )

        # Operator, fill level, fill ratio, msub, rtol, ataol
        if opt("KSMSolver").upper() == "GMRES":
            self.KSM = tacs.TACS.KSM(
                self.K,
                self.PC,
                opt("subSpaceSize"),
                opt("nRestarts"),
                opt("flexible"),
            )
        # TODO: Fix this
        # elif opt('KSMSolver').upper() == 'GCROT':
        #    self.KSM = tacs.TACS.GCROT(
        #        self.K, self.PC, opt('subSpaceSize'), opt('subSpaceSize'),
        #        opt('nRestarts'), opt('flexible'))
        else:
            raise self._TACSError(
                "Unknown KSMSolver option. Valid options are " "'GMRES' or 'GCROT'"
            )

        self.KSM.setTolerances(
            self.getOption("L2ConvergenceRel"), self.getOption("L2Convergence")
        )

        if opt("useMonitor"):
            self.KSM.setMonitor(
                self.comm,
                _descript=opt("KSMSolver").upper(),
                freq=opt("monitorFrequency"),
            )

        # Linear solver factor flag
        self._jacobianUpdateRequired = True
        self._preconditionerUpdateRequired = True

        # Give new vectors and linear solver to nonlinear solver
        for solver in [self.newtonSolver, self.nonlinearSolver]:
            if solver is not None:
                solver.linearSolver = self.KSM
                solver.stateVec = self.u
                solver.resVec = self.res

    def setOption(self, name, value):
        """
        Set a solver option value. The name is not case sensitive.

        Parameters
        ----------
        name : str
            Name of option to modify

        value : depends on option
            New option value to set
        """
        # Default setOption for common problem class objects
        TACSProblem.setOption(self, name, value)

        createVariables = True

        if self.KSM is not None:
            # Update tolerances
            if "l2convergence" in name.lower():
                createVariables = False
                self.KSM.setTolerances(
                    self.getOption("L2ConvergenceRel"),
                    self.getOption("L2Convergence"),
                )
                if self.nonlinearSolver is not None:
                    self.nonlinearSolver.setOption(
                        "newtonSolverRelLinTol", self.getOption("L2ConvergenceRel")
                    )

            # No need to reset solver for output options
            if name.lower() in [
                "writesolution",
                "writenlitersolutions",
                "printtiming",
                "numbersolutions",
                "outputdir",
            ]:
                createVariables = False

            if self.isNonlinear:
                # Pass option to nonlinear solver if it is a nonlinear solver option
                if name.lower() in [
                    opt.lower() for opt in self.nonlinearSolver.defaultOptions
                ]:
                    createVariables = False
                    if self.nonlinearSolver is not None:
                        self.nonlinearSolver.setOption(name, value)

                # We need to create a new solver history object if the monitor variables have updated
                if (
                    name.lower() in ["nonlinearsolvermonitorvars", "newtonsolveruseew"]
                    and self.history is not None
                ):
                    createVariables = False
                    self._createSolverHistory()

            # Reset solver for all other option changes
            if createVariables:
                self._createVariables()

    @property
    def loadScale(self):
        """This is a scaling factor applied to all forcing terms

        Forcing terms includes both the user supplied force vector and the forcing terms coming from aux elements in
        the TACS assembler (e.g inertial, centrifugal forces)

        Returns
        -------
        float or complex
            The current load scale
        """
        return self._loadScale

    @loadScale.setter
    def loadScale(self, value):
        """Set the scaling applied to external loads

        This function exists so that calling `problem.loadScale = value` has the same effect as calling::

            `problem.setLoadScale(value)`

        This is important in case we want to update other things when the load scale is changed in future

        Parameters
        ----------
        value : float or complex
            Value to set the load scale to
        """
        self.setLoadScale(value)

    def setLoadScale(self, value):
        """Set the scaling applied to external loads

        Parameters
        ----------
        value : float or complex
            Value to set the load scale to
        """
        if value != self._loadScale:
            self._jacobianUpdateRequired = True
            self._loadScale = value

    def addFunction(self, funcName, funcHandle, compIDs=None, **kwargs):
        """
        Generic method to add a function for TACS. It is intended to
        be reasonably generic since the user supplies the actual
        function handle to use. See the :py:mod:`~tacs.functions` module
        for supported TACS eval functions.

        Parameters
        ----------
        funcName : str
            The user-supplied name for the function. This will
            typically be a string that is meaningful to the user

        funcHandle : tacs.TACS.Function
            The function handle to use for creation. This must come
            from the functions module in tacs.

        compIDs: list
            List of compIDs to select.

        **kwargs:
            Any keyword arguments to be passed to the TACS function during setup.
        """
        success = TACSProblem.addFunction(self, funcName, funcHandle, compIDs, **kwargs)
        if success:
            # Create additional tacs BVecs to hold adjoint and sens info
            self.adjointList[funcName] = self.assembler.createVec()
            self.dIduList[funcName] = self.assembler.createVec()
            self.dvSensList[funcName] = self.assembler.createDesignVec()
            self.xptSensList[funcName] = self.assembler.createNodeVec()
        return success

    def setDesignVars(self, x):
        """
        Update the design variables used by tacs.

        Parameters
        ----------
        x : numpy.ndarray or dict or tacs.TACS.Vec
            The variables (typically from the optimizer) to set. It
            looks for variable in the ``self.varName`` attribute.

        """
        TACSProblem.setDesignVars(self, x)
        self._jacobianUpdateRequired = True

    def setNodes(self, coords):
        """
        Set the mesh coordinates of the structure.

        Parameters
        ----------
        coords : numpy.ndarray
            Structural coordinate in array of size (N * 3) where N is
            the number of structural nodes on this processor.
        """
        TACSProblem.setNodes(self, coords)
        self._jacobianUpdateRequired = True

    ####### Load adding methods ########

    def addLoadToComponents(self, compIDs, F, averageLoad=False):
        """
        This method is used to add a *FIXED TOTAL LOAD* on one or more
        components, defined by COMPIDs. The purpose of this routine is to add loads that
        remain fixed throughout an optimization. An example would be an engine load.
        This routine determines all the unique nodes in the FE model that are part of the
        requested components, then takes the total 'force' by F and divides by the
        number of nodes. This average load is then applied to the nodes.

        Parameters
        ----------

        compIDs : list[int] or int
            The components with added loads. Use pyTACS selectCompIDs method
            to determine this.

        F : numpy.ndarray 1d or 2d length (varsPerNodes) or (numCompIDs, varsPerNodes)
            Vector(s) of 'force' to apply to each component.  If only one force vector is provided,
            force will be copied uniformly across all components.

        averageLoad : bool
            Flag to determine whether load should be split evenly across all components (True)
            or copied and applied individually to each component (False). Defaults to False.

        Notes
        -----

        The units of the entries of the 'force' vector F are not
        necessarily physical forces and their interpretation depends
        on the physics problem being solved and the dofs included
        in the model.

        A couple of examples of force vector components for common problems are listed below:

            In Heat Conduction with varsPerNode = 1
                F = [Qdot] # heat rate
            In Elasticity with varsPerNode = 3,
                F = [fx, fy, fz] # forces
            In Elasticity with varsPerNode = 6,
                F = [fx, fy, fz, mx, my, mz] # forces + moments
            In Thermoelasticity with varsPerNode = 4,
                F = [fx, fy, fz, Qdot] # forces + heat rate
            In Thermoelasticity with varsPerNode = 7,
                F = [fx, fy, fz, mx, my, mz, Qdot] # forces + moments + heat rate
        """
        self._addLoadToComponents(self.F, compIDs, F, averageLoad)

    def addLoadToNodes(self, nodeIDs, F, nastranOrdering=False):
        """
        This method is used to add a fixed point load of F to the
        selected node IDs.

        Parameters
        ----------

        nodeIDs : list[int]
            The nodes IDs with added loads.

        F : Numpy 1d or 2d array length (varsPerNodes) or (numNodeIDs, varsPerNodes)
            Array of force vectors, one for each node. If only one force vector is provided,
            force will be copied uniformly across all nodes.

        nastranOrdering : bool
            Flag signaling whether nodeIDs are in TACS (default)
            or NASTRAN (grid IDs in bdf file) ordering

        Notes
        -----

        The units of the entries of the 'force' vector F are not
        necessarily physical forces and their interpretation depends
        on the physics problem being solved and the dofs included
        in the model.

        A couple of examples of force vector components for common problems are listed below:

            In Heat Conduction with varsPerNode = 1
                F = [Qdot] # heat rate
            In Elasticity with varsPerNode = 3,
                F = [fx, fy, fz] # forces
            In Elasticity with varsPerNode = 6,
                F = [fx, fy, fz, mx, my, mz] # forces + moments
            In Thermoelasticity with varsPerNode = 4,
                F = [fx, fy, fz, Qdot] # forces + heat rate
            In Thermoelasticity with varsPerNode = 7,
                F = [fx, fy, fz, mx, my, mz, Qdot] # forces + moments + heat rate
        """

        self._addLoadToNodes(self.F, nodeIDs, F, nastranOrdering)

    def addLoadToRHS(self, Fapplied):
        """
        This method is used to add a *FIXED TOTAL LOAD* directly to the
        right hand side vector given the equation below:

            K*u = f

        Where:
            - K : Stiffness matrix for problem
            - u : State variables for problem
            - f : Right-hand side vector to add loads to

        Parameters
        ----------

        Fapplied : numpy.ndarray or tacs.TACS.Vec
            Distributed array containing loads to applied to RHS of the problem.

        """
        self._addLoadToRHS(self.F, Fapplied)

    def addTractionToComponents(self, compIDs, tractions, faceIndex=0):
        """
        This method is used to add a *FIXED TOTAL TRACTION* on one or more
        components, defined by COMPIDs. The purpose of this routine is
        to add loads that remain fixed throughout an optimization.

        Parameters
        ----------

        compIDs : list[int] or int
            The components with added loads. Use pyTACS selectCompIDs method
            to determine this.

        tractions : Numpy array length 1 or compIDs
            Array of traction vectors for each component

        faceIndex : int
            Indicates which face (side) of element to apply traction to.
            Note: not required for certain elements (i.e. shells)
        """
        self._addTractionToComponents(self.auxElems, compIDs, tractions, faceIndex)

    def addTractionToElements(
        self, elemIDs, tractions, faceIndex=0, nastranOrdering=False
    ):
        """
        This method is used to add a fixed traction to the
        selected element IDs. Tractions can be specified on an
        element by element basis (if tractions is a 2d array) or
        set to a uniform value (if tractions is a 1d array)

        Parameters
        ----------

        elemIDs : list[int]
            The global element ID numbers for which to apply the traction.

        tractions : Numpy 1d or 2d array length varsPerNodes or (elemIDs, varsPerNodes)
            Array of traction vectors for each element

        faceIndex : int
            Indicates which face (side) of element to apply traction to.
            Note: not required for certain elements (i.e. shells)

        nastranOrdering : bool
            Flag signaling whether elemIDs are in TACS (default)
            or NASTRAN ordering
        """

        self._addTractionToElements(
            self.auxElems, elemIDs, tractions, faceIndex, nastranOrdering
        )

    def addPressureToComponents(self, compIDs, pressures, faceIndex=0):
        """
        This method is used to add a *FIXED TOTAL PRESSURE* on one or more
        components, defined by COMPIds. The purpose of this routine is
        to add loads that remain fixed throughout an optimization. An example
        would be a fuel load.

        Parameters
        ----------

        compIDs : list[int] or int
            The components with added loads. Use pyTACS selectCompIDs method
            to determine this.

        pressures : Numpy array length 1 or compIDs
            Array of pressure values for each component

        faceIndex : int
            Indicates which face (side) of element to apply pressure to.
            Note: not required for certain elements (i.e. shells)
        """
        self._addPressureToComponents(self.auxElems, compIDs, pressures, faceIndex)

    def addPressureToElements(
        self, elemIDs, pressures, faceIndex=0, nastranOrdering=False
    ):
        """
        This method is used to add a fixed presure to the
        selected element IDs. Pressures can be specified on an
        element by element basis (if pressures is an array) or
        set to a uniform value (if pressures is a scalar)

        Parameters
        ----------

        elemIDs : list[int]
            The global element ID numbers for which to apply the pressure.

        pressures : Numpy array length 1 or elemIDs
            Array of pressure values for each element

        faceIndex : int
            Indicates which face (side) of element to apply pressure to.
            Note: not required for certain elements (i.e. shells)

        nastranOrdering : bool
            Flag signaling whether elemIDs are in TACS (default)
            or NASTRAN ordering
        """

        self._addPressureToElements(
            self.auxElems, elemIDs, pressures, faceIndex, nastranOrdering
        )

    def addInertialLoad(self, inertiaVector):
        """
        This method is used to add a fixed inertial load due to
        a uniform acceleration over the entire model.
        This is most commonly used to model gravity loads on a model.

        Parameters
        ----------
        inertiaVector : numpy.ndarray
            Acceleration vector used to define inertial load.
        """
        self._addInertialLoad(self.auxElems, inertiaVector)

    def addCentrifugalLoad(self, omegaVector, rotCenter, firstOrder=False):
        """
        This method is used to add a fixed centrifugal load due to a
        uniform rotational velocity over the entire model.
        This is most commonly used to model rotors, rolling aircraft, etc.

        Parameters
        ----------

        omegaVector : numpy.ndarray
            Rotational velocity vector (rad/s) used to define centrifugal load.

        rotCenter : numpy.ndarray
            Location of center of rotation used to define centrifugal load.

        firstOrder : bool, optional
            Whether to use first order approximation for centrifugal load,
            which computes the force in the displaced position. By default False
        """
        self._addCentrifugalLoad(self.auxElems, omegaVector, rotCenter, firstOrder)

    def addLoadFromBDF(self, loadID, scale=1.0):
        """
        This method is used to add a fixed load set defined in the BDF file to the problem.
        Currently, only supports LOAD, FORCE, MOMENT, GRAV, RFORCE, PLOAD2, and PLOAD4.

        Parameters
        ----------

        loadID : int
            Load identification number of load set in BDF file user wishes to add to problem.

        scale : float
            Factor to scale the BDF loads by before adding to problem.
        """
        self._addLoadFromBDF(self.F, self.auxElems, loadID, scale)

    ####### Static solver methods ########

    def _updateAssemblerVars(self):
        """
        Make sure that the assembler is using
        the input variables associated with this problem
        """

        self.assembler.setDesignVars(self.x)
        self.assembler.setNodes(self.Xpts)
        self.assembler.setAuxElements(self.auxElems)
        # Set state variables
        self.assembler.setVariables(self.u)
        # Zero any time derivative terms
        self.assembler.zeroDotVariables()
        self.assembler.zeroDDotVariables()
        # Set artificial stiffness factors in rbe class to zero
        tacs.elements.RBE2.setScalingParameters(
            self.getOption("RBEStiffnessScaleFactor"), 0.0
        )
        tacs.elements.RBE3.setScalingParameters(
            self.getOption("RBEStiffnessScaleFactor"), 0.0
        )

    def _initializeSolve(self):
        """
        Initialize the solution of the structural system for the
        loadCase. The stiffness matrix is assembled and factored.
        """

        self.updateJacobian()
        self.updatePreconditioner()

        if self.isNonlinear:
            # Reset the solver history
            if self.rank == 0:
                self.history.reset(clearMetadata=True)
                self.history.addMetadata("Options", self.options)
                self.history.addMetadata("Name", self.name)

    def solve(self, Fext=None):
        """
        Solution of the static problem for current load set. The
        forces must already be set.

        Parameters
        ----------
        Optional Arguments:

        Fext : numpy.ndarray or tacs.TACS.Vec
            Distributed array containing additional loads (ex. aerodynamic forces for aerostructural coupling)
            to applied to RHS of the static problem.

        """
        startTime = time.time()

        self.callCounter += 1
        setupProblemTime = time.time()

        # Set problem vars to assembler
        self._updateAssemblerVars()

        # Check if we need to initialize
        self._initializeSolve()

        initSolveTime = time.time()

        if self.isNonlinear:
            hasConverged = self.solveNonlinear(Fext)
        else:
            # Get current residual
            self.getResidual(self.res, Fext=Fext)

            # Get rhs vector
            self.K.mult(self.u, self.rhs)
            self.rhs.axpy(-1.0, self.res)

            # Set initnorm as the norm of rhs
            self.initNorm = np.real(self.rhs.norm())

            # Starting Norm for this computation
            self.startNorm = np.real(self.res.norm())

            initNormTime = time.time()

            # Solve Linear System for the update
            hasConverged = self.solveLinear(self.res, self.update)

            self.update.scale(-1.0)

            solveTime = time.time()

            # Update State Variables
            self.assembler.getVariables(self.u)
            self.u.axpy(1.0, self.update)
            self.assembler.setVariables(self.u)

            stateUpdateTime = time.time()

            # Get updated residual
            self.getResidual(self.res, Fext)
            self.finalNorm = np.real(self.res.norm())

            finalNormTime = time.time()

            # If timing was was requested print it, if the solution is nonlinear
            # print this information automatically if prinititerations was requested.
            if self.getOption("printTiming"):
                self._pp("+--------------------------------------------------+")
                self._pp("|")
                self._pp("| TACS Solve Times:")
                self._pp("|")
                self._pp(
                    "| %-30s: %10.3f sec"
                    % ("TACS Setup Time", setupProblemTime - startTime)
                )
                self._pp(
                    "| %-30s: %10.3f sec"
                    % ("TACS Solve Init Time", initSolveTime - setupProblemTime)
                )
                self._pp(
                    "| %-30s: %10.3f sec"
                    % ("TACS Init Norm Time", initNormTime - initSolveTime)
                )
                self._pp(
                    "| %-30s: %10.3f sec"
                    % ("TACS Solve Time", solveTime - initNormTime)
                )
                self._pp(
                    "| %-30s: %10.3f sec"
                    % ("TACS State Update Time", stateUpdateTime - solveTime)
                )
                self._pp(
                    "| %-30s: %10.3f sec"
                    % ("TACS Final Norm Time", finalNormTime - stateUpdateTime)
                )
                self._pp("|")
                self._pp(
                    "| %-30s: %10.3f sec"
                    % ("TACS Total Solution Time", finalNormTime - startTime)
                )
                self._pp("+--------------------------------------------------+")

        return hasConverged

    def solveNonlinear(self, Fext=None):
        # Compute the internal and external force components of the residual at the current point
        self.getForces(
            externalForceVec=self.externalForce,
            internalForceVec=self.internalForce,
            Fext=Fext,
        )
        self.initNorm = np.real(self.externalForce.norm())

        # We need to update the residual function handle used by the nonlinear solver based on the current external force vector
        def resFunc(res):
            self.getResidual(res, Fext=Fext)

        self.nonlinearSolver.resFunc = resFunc
        self.newtonSolver.resFunc = resFunc
        self.nonlinearSolver.setRefNorm(self.initNorm)
        self.nonlinearSolver.solve()

        # Since the state has changed, we need to flag that the jacobian and preconditioner should be before the next primal or adjoint solve
        self._preconditionerUpdateRequired = True
        self._jacobianUpdateRequired = True

        # We should also reset the linear solver convergence so that any loosening of the linear convergence criteria are undone before a possible adjoint solve
        self.KSM.setTolerances(
            self.getOption("L2ConvergenceRel"),
            self.getOption("L2Convergence"),
        )

        # Finally return a bool indicating whether the solve was successful
        return self.nonlinearSolver.hasConverged

    def _nonlinearCallback(self, solver, u, res, monitorVars):
        """Callback function to be called by the nonlinear solver at each iteration

        Parameters
        ----------
        solver : pyTACS solver object
            The solver
        u : tacs.TACS.Vec
            Current state vector
        res : tacs.TACS.Vec
            Current residual vector
        monitorVars : dict
            Dictionary of variables to monitor, the values the solver should include can be
            specified through the ``"nonlinearSolverMonitorVars"`` option.
        """
        iteration = 0
        if self.rank == 0:
            # Figure out the iteration number
            iteration = self.history.getIter()
            if iteration % 50 == 0:
                self.history.printHeader()
            self.history.write(monitorVars)
            self.history.printData()

        self.comm.bcast(iteration, root=0)

        if self.getOption("writeNLIterSolutions"):
            # import pdb

            # pdb.set_trace()
            self.writeSolution(
                baseName=f"{self.name}-{self.callCounter:03d}-NLIter", number=iteration
            )

    def solveLinear(self, res, sol):
        """Solve the linear system J * sol = res using the current Jacobian matrix

        Parameters
        ----------
        res : tacs.TACS.Vec
            Right hand side of the linear system
        sol : tacs.TACS.Vec
            Vector to store the solution in

        Returns
        -------
        bool
            Whether the linear solve converged
        """
        success = self.KSM.solve(res, sol)
        success = success == 1

        if not success:
            self._TACSWarning(
                "Linear solver failed to converge. "
                "This is likely a sign that the problem is ill-conditioned. "
                "Check that the model is properly restrained."
            )
        return success

    def updateJacobian(self, res=None):
        """Update the Jacobian (a.k.a stiffness) matrix

        The Jacobian will only actually be updated if the
        ``_jacobianUpdateRequired`` flag is set to True.

        Parameters
        ----------
        res : tacs.TACS.Vec, optional
            If provided, the residual is also computed and stored in this vector
        """
        if self._jacobianUpdateRequired:
            # Assemble residual and stiffness matrix (w/o artificial terms)
            self.assembler.assembleJacobian(
                self.alpha,
                self.beta,
                self.gamma,
                res,
                self.K,
                loadScale=self._loadScale,
            )
            self._preconditionerUpdateRequired = True

    def updatePreconditioner(self):
        """Update the Jacobian (a.k.a stiffness) matrix preconditioner

        By default, the static problem uses a full LU factorization of the
        Jacobian matrix as a preconditioner, which means that this step is
        typically the most expensive operation in the solution process.

        The preconditioner will only actually be updated if the
        ``_preconditionerUpdateRequired`` flag is set to True. This occurs
        whenever the Jacobian is updated.
        """
        if self._preconditionerUpdateRequired:
            # Stiffness matrix must include artificial terms before pc factor
            # to prevent factorization issues w/ zero-diagonals
            self.K.axpy(1.0, self.rbeArtificialStiffness)
            self.PC.factor()
            # Remove artificial stiffness terms to get true stiffness mat
            self.K.axpy(-1.0, self.rbeArtificialStiffness)
            self._preconditionerUpdateRequired = False

    ####### Function eval/sensitivity methods ########

    def evalFunctions(self, funcs, evalFuncs=None, ignoreMissing=False):
        """
        This is the main routine for returning useful information from
        pytacs. The functions corresponding to the strings in
        evalFuncs are evaluated and updated into the provided
        dictionary.

        Parameters
        ----------
        funcs : dict
            Dictionary into which the functions are saved.
        evalFuncs : iterable object containing strings.
            If not none, use these functions to evaluate.
        ignoreMissing : bool
            Flag to supress checking for a valid function. Please use
            this option with caution.

        Examples
        --------
        >>> funcs = {}
        >>> staticProblem.solve()
        >>> staticProblem.evalFunctions(funcs, ['mass'])
        >>> funcs
        >>> # Result will look like (if StaticProblem has name of 'c1'):
        >>> # {'cl_mass':12354.10}
        """
        startTime = time.time()

        # Set problem vars to assembler
        self._updateAssemblerVars()

        if evalFuncs is None:
            evalFuncs = sorted(list(self.functionList))
        else:
            evalFuncs = sorted(list(evalFuncs))

        if not ignoreMissing:
            for f in evalFuncs:
                if f not in self.functionList:
                    raise self._TACSError(
                        f"Supplied function '{f}' has not been added "
                        "using addFunction()."
                    )

        setupProblemTime = time.time()

        # Fast parallel function evaluation of structural funcs:
        handles = [self.functionList[f] for f in evalFuncs if f in self.functionList]
        funcVals = self.assembler.evalFunctions(handles)

        functionEvalTime = time.time()

        # Assign function values to appropriate dictionary
        i = 0
        for f in evalFuncs:
            if f in self.functionList:
                key = self.name + "_%s" % f
                funcs[key] = funcVals[i]
                i += 1

        dictAssignTime = time.time()

        if self.getOption("printTiming"):
            self._pp("+--------------------------------------------------+")
            self._pp("|")
            self._pp("| TACS Function Times:")
            self._pp("|")
            self._pp(
                "| %-30s: %10.3f sec"
                % ("TACS Function Setup Time", setupProblemTime - startTime)
            )
            self._pp(
                "| %-30s: %10.3f sec"
                % (
                    "TACS Function Eval Time",
                    functionEvalTime - setupProblemTime,
                )
            )
            self._pp(
                "| %-30s: %10.3f sec"
                % ("TACS Dict Time", dictAssignTime - functionEvalTime)
            )
            self._pp("|")
            self._pp(
                "| %-30s: %10.3f sec"
                % ("TACS Function Time", dictAssignTime - startTime)
            )
            self._pp("+--------------------------------------------------+")

    def evalFunctionsSens(self, funcsSens, evalFuncs=None):
        """
        This is the main routine for returning useful (sensitivity)
        information from problem. The derivatives of the functions
        corresponding to the strings in evalFuncs are evaluated and
        updated into the provided dictionary. The derivitives with
        respect to all design variables and node locations are computed.

        Parameters
        ----------
        funcsSens : dict
            Dictionary into which the derivatives are saved.
        evalFuncs : iterable object containing strings
            The functions the user wants returned

        Examples
        --------
        >>> funcsSens = {}
        >>> staticProblem.evalFunctionsSens(funcsSens, ['mass'])
        >>> funcsSens
        >>> # Result will look like (if StaticProblem has name of 'c1'):
        >>> # {'c1_mass':{'struct':[1.234, ..., 7.89], 'Xpts':[3.14, ..., 1.59]}}
        """

        startTime = time.time()

        # Set problem vars to assembler
        self._updateAssemblerVars()

        if evalFuncs is None:
            evalFuncs = sorted(list(self.functionList))
        else:
            evalFuncs = sorted(list(evalFuncs))
        # Check that the functions are all ok.
        # and prepare tacs vecs for adjoint procedure
        dvSenses = []
        xptSenses = []
        dIdus = []
        adjoints = []
        for f in evalFuncs:
            if f not in self.functionList:
                raise self._TACSError(
                    "Supplied function has not been added " "using addFunction()"
                )
            else:
                # Populate the lists with the tacs bvecs
                # we'll need for each adjoint/sens calculation
                dvSens = self.dvSensList[f]
                dvSens.zeroEntries()
                dvSenses.append(dvSens)

                xptSens = self.xptSensList[f]
                xptSens.zeroEntries()
                xptSenses.append(xptSens)

                dIdu = self.dIduList[f]
                dIdu.zeroEntries()
                dIdus.append(dIdu)

                adjoint = self.adjointList[f]
                adjoint.zeroEntries()
                adjoints.append(adjoint)

        setupProblemTime = time.time()

        adjointStartTime = {}
        adjointEndTime = {}

        # Next we will solve all the adjoints
        # Set adjoint rhs
        self.addSVSens(evalFuncs, dIdus)
        adjointRHSTime = time.time()
        for i, f in enumerate(evalFuncs):
            adjointStartTime[f] = time.time()
            self.solveAdjoint(dIdus[i], adjoints[i])
            adjointEndTime[f] = time.time()

        adjointFinishedTime = time.time()
        # Evaluate all the adoint res prooduct at the same time for
        # efficiency:
        self.addDVSens(evalFuncs, dvSenses)
        self.addAdjointResProducts(adjoints, dvSenses)
        self.addXptSens(evalFuncs, xptSenses)
        self.addAdjointResXptSensProducts(adjoints, xptSenses)

        # Recast sensititivities into dict for user
        for i, f in enumerate(evalFuncs):
            key = self.name + "_%s" % f
            # Return sensitivities as array in sens dict
            funcsSens[key] = {
                self.varName: dvSenses[i].getArray().copy(),
                self.coordName: xptSenses[i].getArray().copy(),
            }

        totalSensitivityTime = time.time()

        if self.getOption("printTiming"):
            self._pp("+--------------------------------------------------+")
            self._pp("|")
            self._pp("| TACS Adjoint Times:")
            self._pp("|")
            self._pp(
                "| %-30s: %10.3f sec"
                % ("TACS Sens Setup Problem Time", setupProblemTime - startTime)
            )
            self._pp(
                "| %-30s: %10.3f sec"
                % ("TACS Adjoint RHS Time", adjointRHSTime - setupProblemTime)
            )
            for f in evalFuncs:
                self._pp(
                    "| %-30s: %10.3f sec"
                    % (
                        "TACS Adjoint Solve Time - %s" % (f),
                        adjointEndTime[f] - adjointStartTime[f],
                    )
                )
            self._pp(
                "| %-30s: %10.3f sec"
                % (
                    "Total Sensitivity Time",
                    totalSensitivityTime - adjointFinishedTime,
                )
            )
            self._pp("|")
            self._pp(
                "| %-30s: %10.3f sec"
                % (
                    "Complete Sensitivity Time",
                    totalSensitivityTime - startTime,
                )
            )
            self._pp("+--------------------------------------------------+")

    def addSVSens(self, evalFuncs, svSensList):
        """
        Add the state variable partial sensitivity to the ADjoint RHS for given evalFuncs

        Parameters
        ----------
        evalFuncs : list[str]
            The functions the user wants returned

        svSensList : list[tacs.TACS.Vec] or list[numpy.ndarray]
            List of sensitivity vectors to add partial sensitivity to
        """
        # Set problem vars to assembler
        self._updateAssemblerVars()

        # Get list of TACS function handles from evalFuncs
        funcHandles = [
            self.functionList[f] for f in evalFuncs if f in self.functionList
        ]

        # Create a tacs BVec copy for the operation if the output is a numpy array
        if isinstance(svSensList[0], np.ndarray):
            svSensBVecList = [
                self._arrayToVec(svSensArray) for svSensArray in svSensList
            ]
        # Otherwise the input is already a BVec and we can do the operation in place
        else:
            svSensBVecList = svSensList

        self.assembler.addSVSens(
            funcHandles, svSensBVecList, self.alpha, self.beta, self.gamma
        )

        # Update from the BVec values, if the input was a numpy array
        if isinstance(svSensList[0], np.ndarray):
            for svSensArray, svSensBVec in zip(svSensList, svSensBVecList):
                svSensArray[:] = svSensBVec.getArray()

    def addDVSens(self, evalFuncs, dvSensList, scale=1.0):
        """
        Add partial sensitivity contribution due to design vars for evalFuncs

        Parameters
        ----------
        evalFuncs : list[str]
            The functions the user wants returned

        dvSensList : list[tacs.TACS.Vec] or list[numpy.ndarray]
            List of sensitivity vectors to add partial sensitivity to

        scale : float
            Scalar to multiply partial sensitivity by. Defaults to 1.0
        """
        # Set problem vars to assembler
        self._updateAssemblerVars()

        # Get list of TACS function handles from evalFuncs
        funcHandles = [
            self.functionList[f] for f in evalFuncs if f in self.functionList
        ]

        # Create a tacs BVec copy for the operation if the output is a numpy array
        if isinstance(dvSensList[0], np.ndarray):
            dvSensBVecList = [
                self._arrayToDesignVec(dvSensArray) for dvSensArray in dvSensList
            ]
        # Otherwise the input is already a BVec and we can do the operation in place
        else:
            dvSensBVecList = dvSensList

        self.assembler.addDVSens(funcHandles, dvSensBVecList, scale)

        # Finalize sensitivity arrays across all procs
        for dvSensBVec in dvSensBVecList:
            dvSensBVec.beginSetValues()
            dvSensBVec.endSetValues()

        # Update the BVec values, if the input was a numpy array
        if isinstance(dvSensList[0], np.ndarray):
            for dvSensArray, dvSensBVec in zip(dvSensList, dvSensBVecList):
                # Copy values to numpy array
                dvSensArray[:] = dvSensBVec.getArray()

    def addAdjointResProducts(self, adjointlist, dvSensList, scale=-1.0):
        """
        Add the adjoint product contribution to the design variable sensitivity arrays

        Parameters
        ----------
        adjointlist : list[tacs.TACS.Vec] or list[numpy.ndarray]
            List of adjoint vectors for residual sensitivity product

        dvSensList : list[tacs.TACS.Vec] or list[numpy.ndarray]
            List of sensitivity vectors to add product to

        scale : float
            Scalar to multiply product by. Defaults to -1.0
        """
        # Set problem vars to assembler
        self._updateAssemblerVars()

        # Create a tacs BVec copy for the operation if the output is a numpy array
        if isinstance(adjointlist[0], np.ndarray):
            adjointBVeclist = [
                self._arrayToVec(adjointArray) for adjointArray in adjointlist
            ]
        # Otherwise the input is already a BVec and we can do the operation in place
        else:
            adjointBVeclist = adjointlist

        # Make sure BC terms are zeroed out in adjoint
        for adjoint in adjointBVeclist:
            self.assembler.applyBCs(adjoint)

        # Create a tacs BVec copy for the operation if the output is a numpy array
        if isinstance(dvSensList[0], np.ndarray):
            dvSensBVecList = [
                self._arrayToDesignVec(dvSensArray) for dvSensArray in dvSensList
            ]
        # Otherwise the input is already a BVec and we can do the operation in place
        else:
            dvSensBVecList = dvSensList

        self.assembler.addAdjointResProducts(adjointBVeclist, dvSensBVecList, scale)

        # Finalize sensitivity arrays across all procs
        for dvSensBVec in dvSensBVecList:
            dvSensBVec.beginSetValues()
            dvSensBVec.endSetValues()

        # Update the BVec values, if the input was a numpy array
        if isinstance(dvSensList[0], np.ndarray):
            for dvSensArray, dvSensBVec in zip(dvSensList, dvSensBVecList):
                # Copy values to numpy array
                dvSensArray[:] = dvSensBVec.getArray()

    def addXptSens(self, evalFuncs, xptSensList, scale=1.0):
        """
        Add partial sensitivity contribution due to nodal coordinates for evalFuncs

        Parameters
        ----------
        evalFuncs : list[str]
            The functions the user wants returned

        xptSensList : list[tacs.TACS.Vec] or list[numpy.ndarray]
            List of sensitivity vectors to add partial sensitivity to

        scale : float
            Scalar to multiply partial sensitivity by. Defaults to 1.0
        """
        # Set problem vars to assembler
        self._updateAssemblerVars()

        # Get list of TACS function handles from evalFuncs
        funcHandles = [
            self.functionList[f] for f in evalFuncs if f in self.functionList
        ]

        # Create a tacs BVec copy for the operation if the output is a numpy array
        if isinstance(xptSensList[0], np.ndarray):
            xptSensBVecList = [
                self._arrayToNodeVec(xptSensArray) for xptSensArray in xptSensList
            ]
        # Otherwise the input is already a BVec and we can do the operation in place
        else:
            xptSensBVecList = xptSensList

        self.assembler.addXptSens(funcHandles, xptSensBVecList, scale)

        # Finalize sensitivity arrays across all procs
        for xptSensBVec in xptSensBVecList:
            xptSensBVec.beginSetValues()
            xptSensBVec.endSetValues()

        # Update from the BVec values, if the input was a numpy array
        if isinstance(xptSensList[0], np.ndarray):
            for xptSensArray, xptSensBVec in zip(xptSensList, xptSensBVecList):
                # Copy values to numpy array
                xptSensArray[:] = xptSensBVec.getArray()

    def addAdjointResXptSensProducts(self, adjointlist, xptSensList, scale=-1.0):
        """
        Add the adjoint product contribution to the nodal coordinates sensitivity arrays

        Parameters
        ----------
        adjointlist : list[tacs.TACS.Vec] or list[numpy.ndarray]
            List of adjoint vectors for residual sensitivity product

        xptSensList : list[tacs.TACS.Vec] or list[numpy.ndarray]
            List of sensitivity vectors to add product to

        scale : float
            Scalar to multiply product by. Defaults to -1.0
        """
        # Set problem vars to assembler
        self._updateAssemblerVars()

        # Create a tacs BVec copy for the operation if the output is a numpy array
        if isinstance(adjointlist[0], np.ndarray):
            adjointBVeclist = [
                self._arrayToVec(adjointArray) for adjointArray in adjointlist
            ]
        # Otherwise the input is already a BVec and we can do the operation in place
        else:
            adjointBVeclist = adjointlist

        # Make sure BC terms are zeroed out in adjoint
        for adjoint in adjointBVeclist:
            self.assembler.applyBCs(adjoint)

        # Create a tacs BVec copy for the operation if the output is a numpy array
        if isinstance(xptSensList[0], np.ndarray):
            xptSensBVecList = [
                self._arrayToNodeVec(xptSensArray) for xptSensArray in xptSensList
            ]
        # Otherwise the input is already a BVec and we can do the operation in place
        else:
            xptSensBVecList = xptSensList

        self.assembler.addAdjointResXptSensProducts(
            adjointBVeclist, xptSensBVecList, scale
        )

        # Finalize sensitivity arrays across all procs
        for xptSensBVec in xptSensBVecList:
            xptSensBVec.beginSetValues()
            xptSensBVec.endSetValues()

        if isinstance(xptSensList[0], np.ndarray):
            for xptSensArray, xptSensBVec in zip(xptSensList, xptSensBVecList):
                # Copy values to numpy array
                xptSensArray[:] = xptSensBVec.getArray()

    def getResidual(self, res, Fext=None):
        """
        This routine is used to evaluate directly the structural
        residual. Only typically used with aerostructural analysis.

        Parameters
        ----------
        res : tacs.TACS.Vec or numpy.ndarray
            If res is not None, place the residuals into this array.

        Fext : tacs.TACS.Vec or numpy.ndarray, optional
            Distributed array containing additional loads (ex. aerodynamic forces for aerostructural coupling)
            to applied to RHS of the static problem.

        """
        # Make sure assembler variables are up-to-date
        self._updateAssemblerVars()

        # Determine if the user vector is a BVec or numpy array
        if isinstance(res, tacs.TACS.Vec):
            resArray = None
        else:  # Input is a numpy array
            resArray = res
            res = self.res

        # Sum the forces from the loads not handled by TACS
        self.rhs.copyValues(self.F)  # Fixed loads

        # Add external loads, if specified
        if Fext is not None:
            if isinstance(Fext, tacs.TACS.Vec):
                self.rhs.axpy(1.0, Fext)
            elif isinstance(Fext, np.ndarray):
                rhsArray = self.rhs.getArray()
                rhsArray[:] = rhsArray[:] + Fext[:]

        # Zero out forces on DOF that are subject to BCs
        self.assembler.applyBCs(self.rhs)

        # Assemble the TACS residual and subtract the externally handled loads
        self.assembler.assembleRes(res, self._loadScale)
        res.axpy(-self._loadScale, self.rhs)

        # If requested, copy the residual to the output array
        if resArray is not None:
            resArray[:] = res.getArray()

    def getForces(self, externalForceVec, internalForceVec, Fext=None):
        """Compute the internal and external forces acting on the structure

        The computations here are based on the residual equation:
            r(u, loadScale) = -(Fint(u) + loadScale * Fext(u))
        Thus, the internal forces are given by:
            Fint(u) = -r(u, 0)
        And the external forces are given by:
            Fext(u) = -r(u, 1) - Fi(u)

        Parameters
        ----------
        externalForceVec : TACS BVec or numpy array
            Vector/array to store external forces in
        internalForceVec : TACS BVec or numpy array
            Vector/array to store internal forces in
        Fext : TACS BVec or numpy array, optional
            Distributed array containing additional loads (ex. aerodynamic forces for aerostructural coupling)
            to applied to RHS of the static problem.
        """
        loadScale = self._loadScale
        self.setLoadScale(0.0)
        self.getResidual(internalForceVec, Fext)
        self.setLoadScale(1.0)
        self.getResidual(externalForceVec, Fext)
        self.setLoadScale(loadScale)

        # Compute internal forces
        if isinstance(internalForceVec, tacs.TACS.Vec):
            internalForceVec.scale(-1.0)
        elif isinstance(internalForceVec, np.ndarray):
            internalForceVec[:] = -internalForceVec[:]

        # Compute external forces
        if isinstance(externalForceVec, tacs.TACS.Vec):
            externalForceVec.scale(-1.0)
        elif isinstance(externalForceVec, np.ndarray):
            externalForceVec[:] = -externalForceVec[:]

        if isinstance(internalForceVec, tacs.TACS.Vec):
            if isinstance(externalForceVec, tacs.TACS.Vec):
                externalForceVec.axpy(-1.0, internalForceVec)
            elif isinstance(externalForceVec, np.ndarray):
                externalForceVec[:] = externalForceVec[:] - internalForceVec.getArray()
        elif isinstance(internalForceVec, np.ndarray):
            if isinstance(externalForceVec, np.ndarray):
                externalForceVec[:] = externalForceVec[:] - internalForceVec[:]
            elif isinstance(externalForceVec, tacs.TACS.Vec):
                externalForceVec.axpy(-1.0, self._arrayToVec(internalForceVec))

    def getJacobian(self):
        """Get the problem's Jacobian in sciPy sparse matrix format

        Returns
        -------
        K : (scipy.sparse.bsr_matrix, scipy.sparse.bsr_matrix) or (scipy.sparse.bsr_matrix, scipy.sparse.bsr_matrix, scipy.sparse.bsr_matrix, scipy.sparse.bsr_matrix)
            A tuple of 2 scipy.sparse.bsr_matrices (A, B) if Jacobian is a TACSParallelMat, or 4
            scipy.sparse.bsr_matrices (A, B, C, D) if Jacobian is a TACSSchurMat
        """
        # Make sure stiffness mat is up-to-date
        self._updateAssemblerVars()
        self._initializeSolve()
        # Return copy of scipy mat
        return copy.deepcopy(self.K.getMat())

    def addTransposeJacVecProduct(self, phi, prod, scale=1.0):
        """
        Adds product of transpose Jacobian and input vector into output vector as shown below:
        prod += scale * J^T . phi

        Parameters
        ----------
        phi : tacs.TACS.Vec or numpy.ndarray
            Input vector to product with the transpose Jacobian.

        prod : tacs.TACS.Vec or numpy.ndarray
            Output vector to add Jacobian product to.

        scale : float
            Scalar used to scale Jacobian product by.
        """
        # Create a tacs bvec copy of the adjoint vector
        if isinstance(phi, tacs.TACS.Vec):
            self.phi.copyValues(phi)
        elif isinstance(phi, np.ndarray):
            self.phi.getArray()[:] = phi

        # Tacs doesn't actually transpose the matrix here so keep track of
        # RHS entries that TACS zeros out for BCs.
        bcTerms = self.update
        bcTerms.copyValues(self.phi)
        self.assembler.applyBCs(self.phi)
        bcTerms.axpy(-1.0, self.phi)

        # Set problem vars to assembler
        self._updateAssemblerVars()

        self.K.mult(self.phi, self.res)
        # Add bc terms back in
        self.res.axpy(1.0, bcTerms)

        # Output residual
        if isinstance(prod, tacs.TACS.Vec):
            prod.axpy(scale, self.res)
        else:
            prod[:] = prod + scale * self.res.getArray()

    def zeroVariables(self):
        """
        Zero all the tacs solution b-vecs
        """
        self.res.zeroEntries()
        self.u.zeroEntries()
        self.assembler.setVariables(self.u)
        self.update.zeroEntries()

    def zeroLoads(self):
        """
        Zero all applied loads
        """
        self.F.zeroEntries()
        self.auxElems = tacs.TACS.AuxElements()

    def solveAdjoint(self, rhs, phi):
        """
        Solve the structural adjoint.

        Parameters
        ----------
        rhs : tacs.TACS.Vec or numpy.ndarray
            right hand side vector for adjoint solve
        phi : tacs.TACS.Vec or numpy.ndarray
            BVec or numpy array into which the adjoint is saved
        """

        # Set problem vars to assembler
        self._updateAssemblerVars()

        # Check if we need to initialize
        self._initializeSolve()

        # Create a copy of the adjoint/rhs guess
        if isinstance(phi, tacs.TACS.Vec):
            self.phi.copyValues(phi)
        elif isinstance(phi, np.ndarray):
            self.phi.getArray()[:] = phi

        if isinstance(rhs, tacs.TACS.Vec):
            self.adjRHS.copyValues(rhs)
        elif isinstance(rhs, np.ndarray):
            self.adjRHS.getArray()[:] = rhs

        # Tacs doesn't actually transpose the matrix here so keep track of
        # RHS entries that TACS zeros out for BCs.
        bcTerms = self.update
        bcTerms.copyValues(self.adjRHS)
        self.assembler.applyBCs(self.adjRHS)
        bcTerms.axpy(-1.0, self.adjRHS)

        # Solve Linear System
        self.solveLinear(self.adjRHS, self.phi)
        self.assembler.applyBCs(self.phi)
        # Add bc terms back in
        self.phi.axpy(1.0, bcTerms)

        # Copy output values back to user vectors
        if isinstance(phi, tacs.TACS.Vec):
            phi.copyValues(self.phi)
        elif isinstance(phi, np.ndarray):
            phi[:] = self.phi.getArray()

    def getVariables(self, states=None):
        """
        Return the current state values for the
        problem

        Parameters
        ----------
        states : tacs.TACS.Vec or numpy.ndarray
            Vector to place current state variables into (optional)

        Returns
        ----------
        states : numpy.ndarray
            current state vector
        """

        if isinstance(states, tacs.TACS.Vec):
            states.copyValues(self.u)
        elif isinstance(states, np.ndarray):
            states[:] = self.u_array[:]

        return self.u_array.copy()

    def setVariables(self, states):
        """
        Set the structural states for current load case.

        Parameters
        ----------
        states : numpy.ndarray
            Values to set. Must be the size of getNumVariables()
        """
        # Copy array values
        if isinstance(states, tacs.TACS.Vec):
            self.u.copyValues(states)
        elif isinstance(states, np.ndarray):
            self.u_array[:] = states[:]
        # Apply boundary conditions
        self.assembler.applyBCs(self.u)
        # Set states to assembler
        self.assembler.setVariables(self.u)

        # If this is a nonlinear problem then changing the state will change the jacobian
        if self.isNonlinear:
            self._jacobianUpdateRequired = True

    def getOutputFileName(self, outputDir=None, baseName=None, number=None):
        """Figure out a base path/name for output files

        Parameters
        ----------
        outputDir : str, optional
            Directory to write file to, by default uses the 'outputDir' option
        baseName : str, optional
            Case name, by default uses self.name
        number : int, optional
            A number to append to the filename, by default uses the call-count of the problem

        Returns
        -------
        str
            Full path to output file (excluding any extension)
        """
        # Check input
        if outputDir is None:
            outputDir = self.getOption("outputDir")

        if baseName is None:
            baseName = self.name

        # If we are numbering solution, it saving the sequence of
        # calls, add the call number
        if number is not None:
            # We need number based on the provided number:
            baseName = baseName + "_%3.3d" % number
        else:
            # if number is none, i.e. standalone, but we need to
            # number solutions, use internal counter
            if self.getOption("numberSolutions"):
                baseName = baseName + "_%3.3d" % self.callCounter
        return os.path.join(outputDir, baseName)

    def writeSolution(self, outputDir=None, baseName=None, number=None):
        """
        This is a generic shell function that writes the output
        file(s).  The intent is that the user or calling program can
        call this function and pyTACS writes all the files that the
        user has defined. It is recommended that this function is used
        along with the associated logical flags in the options to
        determine the desired writing procedure

        Parameters
        ----------
        outputDir : str or None
            Use the supplied output directory
        baseName : str or None
            Use this supplied string for the base filename. Typically
            only used from an external solver.
        number : int or None
            Use the user supplied number to index solution. Again, only
            typically used from an external solver
        """
        # Make sure assembler variables are up to date
        self._updateAssemblerVars()

        # Figure out the output file base name
        baseName = self.getOutputFileName(outputDir, baseName, number)

        # Unless the writeSolution option is off write actual file:
        if self.getOption("writeSolution"):
            self.outputViewer.writeToFile(baseName + ".f5")

    def writeSolutionHistory(self, outputDir=None, baseName=None, number=None):
        """Write the nonlinear solver history to a file

        Parameters
        ----------
        outputDir : str or None
            Use the supplied output directory
        baseName : str or None
            Use this supplied string for the base filename. Typically
            only used from an external solver.
        number : int or None
            Use the user supplied number to index solution. Again, only
            typically used from an external solver
        """
        # Figure out the output file base name
        baseName = self.getOutputFileName(outputDir, baseName, number)
        if self.history is not None:
            self.history.save(baseName)

    def writeLoadToBDF(self, bdfFile, loadCaseID):
        """
        Write loads from problem to NASTRAN BDF file.
        NOTE: To get correct loads, `solve` method should be called before this method.

        Parameters
        ----------
        bdfFile: str or pyNastran.bdf.bdf.BDF or None
            Name of file to write BDF file to. Only required on root proc,
            can be None otherwise.
        loadCaseID: int
            NASTARAN loadcase ID to assign loads to in BDF.
        """

        # Grab RHS vector from previous solve
        F = self.rhs
        F_array = np.real(F.getArray())

        # Get local force info for each processor
        multNodes = self.meshLoader.getLocalMultiplierNodeIDs()
        globalToLocalNodeIDDict = self.meshLoader.getGlobalToLocalNodeIDDict()

        # Gather local info to root processor
        allMultNodes = self.comm.gather(multNodes, root=0)
        allGlobalToLocalNodeIDDict = self.comm.gather(globalToLocalNodeIDDict, root=0)
        allF = self.comm.gather(F_array, root=0)

        vpn = self.getVarsPerNode()

        # Assemble new BDF file on root
        if self.comm.rank == 0:
            if isinstance(bdfFile, str):
                newBDFInfo = pn.bdf.BDF(debug=False)
            elif isinstance(bdfFile, pn.bdf.BDF):
                newBDFInfo = bdfFile

            # Save subcase info to bdf
            if newBDFInfo.case_control_deck is not None:
                newBDFInfo.case_control_deck.create_new_subcase(loadCaseID)
                newBDFInfo.case_control_deck.add_parameter_to_local_subcase(
                    loadCaseID, f"SUBTITLE={self.name}"
                )
                newBDFInfo.case_control_deck.add_parameter_to_local_subcase(
                    loadCaseID, f"ANALYSIS=STATICS"
                )
                newBDFInfo.case_control_deck.add_parameter_to_local_subcase(
                    loadCaseID, f"LOAD={loadCaseID}"
                )

            # Tolerance for writing out point loads
            zero_tol = 1e-6
            # Write out force values
            nastranNodeIDs = list(self.bdfInfo.node_ids)
            # Loop through each proc and pull out nodal forces
            for proc_i in range(self.comm.size):
                Fxyz = allF[proc_i].reshape(-1, vpn)
                for tacsGNodeID in allGlobalToLocalNodeIDDict[proc_i]:
                    # Get local node ID
                    tacsLNodeID = allGlobalToLocalNodeIDDict[proc_i][tacsGNodeID]
                    # Get Global nastran ID
                    nastranGNodeID = nastranNodeIDs[tacsGNodeID]
                    # Add force to bdf file (if its not a multiplier node)
                    if tacsLNodeID not in allMultNodes[proc_i]:
                        # Check if force is above tolerance before adding to bdf
                        if (
                            vpn >= 3
                            and np.linalg.norm(Fxyz[tacsLNodeID][:3]) > zero_tol
                        ):
                            f = np.zeros(3)
                            for i in range(3):
                                if abs(Fxyz[tacsLNodeID][i]) > zero_tol:
                                    f[i] = Fxyz[tacsLNodeID][i]
                            newBDFInfo.add_force(loadCaseID, nastranGNodeID, 1.0, f)
                        if (
                            vpn >= 6
                            and np.linalg.norm(Fxyz[tacsLNodeID][3:6]) > zero_tol
                        ):
                            m = np.zeros(3)
                            for i in range(3):
                                if abs(Fxyz[tacsLNodeID][i + 3]) > zero_tol:
                                    m[i] = Fxyz[tacsLNodeID][i + 3]
                            newBDFInfo.add_moment(loadCaseID, nastranGNodeID, 1.0, m)

            # If bdf file was provided as a file name save it directly
            if isinstance(bdfFile, str):
                newBDFInfo.write_bdf(
                    bdfFile, size=16, is_double=True, write_header=False
                )

        # All procs should wait for root
        self.comm.barrier()
