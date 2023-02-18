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

import tacs.TACS
import tacs.elements
from .base import TACSProblem
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
        # Output Options
        "writeSolution": [
            bool,
            True,
            "Flag for suppressing all f5 file writing.",
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
        # Nonlinear continuation options
        "continuationTargetIter": [
            int,
            8,
            "Target number of Newton iterations for each continuation increment.",
        ],
        "continuationMaxIter": [int, 100, "Maximum number of continuation steps."],
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
        "usePredictor": [bool, True, "Flag for using predictor step in continuation."],
        "predictorNumStates": [
            int,
            2,
            "Number of previous equilibrium states to use in computing the predictor step.",
        ],
        "predictorUseDerivative": [
            bool,
            False,
            "Whether to use the equilibrium path slope in the computation of the predictor step. This requires a linear solve and thus greatly increases the cost of the predictor step computation.",
        ],
        # Newton solver options
        "newtonSolverMonitorVars": [
            list,
            [
                "linSolverIters",
                "linSolverRes",
                "loadScale",
                "lineSearchStep",
                "lineSearchIters",
            ],
            "List of variables to include in nonlinear solver monitor output. Choose from 'linSolverIters', 'linSolverRes', 'loadScale', 'lineSearchStep' and 'lineSearchIters'.",
        ],
        "newtonSolverMaxIter": [int, 40, "Maximum number of Newton iterations."],
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
        "newtonSolverCoarseAbsTol": [
            float,
            1e-4,
            "Residual norm criteria for intermediate continuation steps, making this larger may speed up the nonlinear solver by allowing it to only partially converge intermediate steps.",
        ],
        "newtonSolverCoarseRelTol": [
            float,
            1e-4,
            "Relative residual norm criteria for intermediate load increments.",
        ],
        "newtonSolverDivergenceTol": [
            float,
            1e10,
            "Residual norm at which the nonlinear solver is jugded to have diverged",
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
        name,
        assembler,
        comm,
        outputViewer=None,
        meshLoader=None,
        isNonlinear=False,
        options={},
    ):
        """
        NOTE: This class should not be initialized directly by the user.
        Use pyTACS.createStaticProblem instead.

        Parameters
        ----------
        name : str
            Name of this tacs problem

        assembler : TACS.Assembler
            Cython object responsible for creating and setting tacs objects used to solve problem

        comm : mpi4py.MPI.Intracomm
            The comm object on which to create the pyTACS object.

        outputViewer : TACS.TACSToFH5
            Cython object used to write out f5 files that can be converted and used for postprocessing.

        meshLoader : pymeshloader.pyMeshLoader
            pyMeshLoader object used to create the assembler.

        options : dict
            Dictionary holding problem-specific option parameters (case-insensitive).

        """

        # Problem name
        self.name = name

        # Default setup for common problem class objects
        TACSProblem.__init__(
            self, assembler, comm, outputViewer, meshLoader, isNonlinear=isNonlinear
        )

        # Process the default options which are added to self.options
        # under the 'defaults' key. Make sure the key are lower case
        def_keys = self.defaultOptions.keys()
        self.options["defaults"] = {}
        for key in def_keys:
            self.options["defaults"][key.lower()] = self.defaultOptions[key]
            self.options[key.lower()] = self.defaultOptions[key]

        # Set user-defined options
        for key in options:
            TACSProblem.setOption(self, key, options[key])

        # Setup solver history object for nonlinear problems
        self.history = None
        if self.isNonlinear:
            self._createSolverHistory()

        # Create problem-specific variables
        self._createVariables()

    def _createSolverHistory(self):
        """Setup the solver history object based on the current options

        The solver history is only created on the root processor.
        """
        monitorVars = [s.lower() for s in self.getOption("newtonSolverMonitorVars")]
        numType = float if self.dtype == np.float64 else complex
        if self.comm.rank == 0:
            history = SolverHistory()

            # Define the variables to be stored in the history
            # Continuation increment number
            history.addVariable("Increment", int, printVar=True)
            # Load scale
            history.addVariable(
                "Load scale", float, printVar="loadscale" in monitorVars
            )
            # Newton solve iteration number
            history.addVariable("SubIter", int, printVar=True)
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

        # Vectors used to compute extrapolate the equilibrium path during nonlinear solutions
        self.equilibriumPathStates = []
        self.equilibriumPathLoadScales = []
        if self.isNonlinear and opt("usePredictor"):
            for _ in range(opt("predictorNumStates")):
                self.equilibriumPathStates.append(self.assembler.createVec())
                self.equilibriumPathLoadScales.append(None)

        if self.isNonlinear:
            self.u_inc_start = self.assembler.createVec()

        # Auxiliary element object for applying tractions/pressure
        self.auxElems = tacs.TACS.AuxElements()
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
        self._stiffnessUpdateRequired = True
        self._factorOnNext = True

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

        # Update tolerances
        if "l2convergence" in name.lower():
            self.KSM.setTolerances(
                self.getOption("L2ConvergenceRel"),
                self.getOption("L2Convergence"),
            )
        # No need to reset solver for output options
        elif (
            name.lower()
            in [
                "writesolution",
                "printtiming",
                "numbersolutions",
                "outputdir",
                "usePredictor",
                "predictorUseDerivative",
            ]
            or "linesearch" in name.lower()
            or "newtonsolver" in name.lower()
        ):
            pass
        # Reset solver for all other option changes
        else:
            self._createVariables()

        # We need to create a new solver history object if the monitor variables have updated
        if name.lower() == "newtonsolvermonitorvars":
            self._createSolverHistory()

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
            self._stiffnessUpdateRequired = True
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

        funcHandle : TACS.Function
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
        x : numpy.ndarray
            The variables (typically from the optimizer) to set. It
            looks for variable in the ``self.varName`` attribute.

        """
        TACSProblem.setDesignVars(self, x)
        self._stiffnessUpdateRequired = True

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
        self._stiffnessUpdateRequired = True

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
            Vector(s) of 'force' to apply to each components.  If only one force vector is provided,
            force will be copied uniformly across all components.

        averageLoad : bool
            Flag to determine whether load should be split evenly across all components (True)
            or copied and applied individually to each component (False). Defaults to False.

        Notes
        ----------

        The units of the entries of the 'force' vector F are not
        necessarily physical forces and their interpretation depends
        on the physics problem being solved and the dofs included
        in the model.

        A couple of examples of force vector components for common problem are listed below:

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
        ----------

        The units of the entries of the 'force' vector F are not
        necessarily physical forces and their interpretation depends
        on the physics problem being solved and the dofs included
        in the model.

        A couple of examples of force vector components for common problem are listed below:

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

        Fapplied : numpy.ndarray or TACS.Vec
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
            Array of traction vectors for each components

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
            Array of pressure values for each components

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

    def solve(self, Fext=None):
        """
        Solution of the static problem for current load set. The
        forces must already be set.

        Parameters
        ----------
        Optional Arguments:

        Fext : numpy.ndarray or TACS.Vec
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
            self.solveNonlinear(Fext)
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
            self.solveJacLinear(self.res, self.update)

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

        return

    def solveNonlinear(self, Fext=None, maxLoadScale=1.0):
        TARGET_ITERS = self.getOption("continuationTargetIter")
        INIT_STEP = self.getOption("continuationInitialStep")
        MIN_STEP = self.getOption("continuationMinStep")
        MAX_STEP = self.getOption("continuationMaxStep")
        MAX_INCREMENTS = self.getOption("continuationMaxIter")
        MIN_STEP_FACTOR = self.getOption("continuationMinStepFactor")
        MAX_STEP_FACTOR = self.getOption("continuationMaxStepFactor")
        STEP_RETRACT_FACTOR = self.getOption("continuationRetractionFactor")

        USE_PREDICTOR = self.getOption("usePredictor")
        NUM_PREDICTOR_STATES = self.getOption("predictorNumStates")
        PREDICTOR_USE_DERIVATIVE = self.getOption("predictorUseDerivative")

        # Compute the internal and external force components of the residual at the current point
        self.getForces(
            externalForceVec=self.externalForce,
            internalForceVec=self.internalForce,
            Fext=Fext,
        )
        self.initNorm = np.real(self.externalForce.norm())

        # ==============================================================================
        # Compute the initial load scale
        # ==============================================================================
        self.setLoadScale(min(maxLoadScale, INIT_STEP))
        loadStepDirection = 1

        # If we're restarting from a previous solution we should compute the optimum load scale
        # to restart from. This is done by computing the load scale that minimizes the work
        # done by the resulting Newton step:
        # optLoadScale = (Fe^T dUi + Fi^T dUe) / (-2 Fe^T dUe)
        # Where: Fe = external force, Fi = internal force, dUi = inv(K) * Fi, dUe = inv(K) * Fe
        if np.real(self.u.norm()) > 0:
            if self._stiffnessUpdateRequired:
                self.updateJacobian()
            if self._factorOnNext:
                self.updatePreconditioner()
            du_i = self.u
            du_e = self.update
            self.solveJacLinear(self.externalForce, du_e)
            self.solveJacLinear(self.internalForce, du_i)
            FeUe = self.externalForce.dot(du_e)
            FeUi = self.externalForce.dot(du_i)
            FiUe = self.internalForce.dot(du_e)
            optLoadScale = (FeUi + FiUe) / (-2 * FeUe)

            if optLoadScale > 2 * maxLoadScale or optLoadScale < 0.0:
                # If the optimum load scale is more than double the max load scale we're aiming for, or if it's
                # negative then the loading/structure has changed so much that we'll be closer to the final
                # solution if we just reset the displacements to zero and start the solver from there
                self.zeroVariables()
                optLoadScale = self.loadScale
            elif np.abs(optLoadScale - self.loadScale) < 1e-2:
                # If the optimum load scale is close to the max load scale then we'll just use the max load scale
                optLoadScale = maxLoadScale
            else:
                # Otherwise choose the maximum of the ideal load scale and the default initial load scale
                optLoadScale = max(optLoadScale, self.loadScale)
                # If the optimum load scale is greater than the max we want to get to then we need to reverse the
                # direction of load incrementation
                if optLoadScale > maxLoadScale:
                    loadStepDirection = -1

            self.setLoadScale(optLoadScale)

        stepSize = INIT_STEP

        # Reset the solver history
        if self.rank == 0:
            self.history.reset(clearMetadata=True)
            self.history.addMetadata("Options", self.options)
            self.history.addMetadata("Name", self.name)

        for increment in range(MAX_INCREMENTS):

            # Save displacement at start of this increment, this is what
            # we'll reset to if the increment diverges
            self.u_inc_start.copyValues(self.u)

            # --- Compute predictor step ---
            # TODO: Add predictor computation here

            success, numIters = self.newtonSolve(Fext=Fext)

            # --- Check convergence ---
            if not success:
                # If the Newton solve failed then we'll reduce the step size and try again
                self.setVariables(self.u_inc_start)
                self.setLoadScale(self.loadScale - stepSize * loadStepDirection)
                stepSize *= STEP_RETRACT_FACTOR
            else:
                if self.loadScale == maxLoadScale:
                    break
                else:
                    stepChangeFactor = np.sqrt(TARGET_ITERS / numIters)
                    stepSize *= np.clip(
                        stepChangeFactor, MIN_STEP_FACTOR, MAX_STEP_FACTOR
                    )
                    if USE_PREDICTOR:
                        stateToOverwrite = self.equilibriumPathStates.pop(0)
                        stateToOverwrite.copyValues(self.u)
                        self.equilibriumPathStates.append(stateToOverwrite)

                        self.equilibriumPathLoadScales.pop(0)
                        self.equilibriumPathLoadScales.append(self.loadScale)

            maxStep = min(np.abs(maxLoadScale - self.loadScale), MAX_STEP)
            stepSize = np.clip(stepSize, MIN_STEP, maxStep)
            self.setLoadScale(self.loadScale + loadStepDirection * stepSize)

        # ==============================================================================
        # End of nonlinear solution
        # ==============================================================================

    def newtonSolve(self, Fext=None):
        USE_LINESEARCH = self.getOption("useLineSearch")
        LINESEARCH_SKIP_ITERS = self.getOption("skipFirstNLineSearch")
        MAX_ITERS = self.getOption("newtonSolverMaxIter")


        for iteration in range(MAX_ITERS):
            # TODO: Write output file here based on option
            # self.writeSolution(baseName=f"{self.name}-NLIter", number=iteration)

            # Compute residual
            self.getResidual(self.res, Fext=Fext)
            resNorm = self.res.norm()
            uNorm = self.u.norm()

            # Write data to history
            histData = {
                "Increment":0,
                "SubIter": iteration,
                "Load scale": self.loadScale,
                "Res norm": resNorm,
                "Rel res norm":resNorm/self.initNorm,
                "U norm":uNorm,
            }
            if iteration > 0:
                histData["Lin iters"] = linearSolveIterations
                histData["Lin res"] = linearSolveResNorm
                histData["LS step"] = alpha
                histData["LS iters"] = lineSearchIters
            if self.rank == 0:
                self.history.write(histData)
                if iteration % 50 == 0:
                    self.history.printHeader()
                self.history.printData()


            # Test convergence (exit if converged/diverged)
            hasConverged = self.checkConvergence(resNorm)
            hasDiverged = self.checkDivergence(resNorm)
            if hasConverged or hasDiverged:
                break

            # Update Jacobian
            self.updateJacobian()
            self.updatePreconditioner()

            # Compute Newton step
            self.solveJacLinear(self.res, self.update)
            self.update.scale(-1.0)

            # Check data from linear solve
            linearSolveIterations = self.KSM.getIterCount()
            linearSolveResNorm = self.KSM.getResidualNorm()

            if USE_LINESEARCH and iteration >= LINESEARCH_SKIP_ITERS:
                # Do linesearch
                alpha, lineSearchIters = self.energyLineSearch(self.u, self.update, Fext=Fext)
            else:
                alpha = 1.0
                lineSearchIters = 1
            self.u.axpy(alpha, self.update)
            self.assembler.setVariables(self.u)
            self._stiffnessUpdateRequired = True

        return hasConverged, iteration

    def energyLineSearch(self, u, stepDir, Fext=None):
        MAX_LINESEARCH_ITERS = self.getOption("lineSearchMaxIter")
        LINESEARCH_MU = self.getOption("lineSearchExpectedDecrease")
        LINESEARCH_ALPHA_MIN = self.getOption("lineSearchMinStep")
        LINESEARCH_ALPHA_MAX = self.getOption("lineSearchMaxStep")
        LINESEARCH_MAX_STEP_CHANGE = self.getOption("lineSearchMaxStepChange")
        PRINT_LINESEARCH_ITERS = self.getOption("lineSearchMonitor")

        # Compute residual and merit function at u0
        self.assembler.setVariables(u)
        self.getResidual(self.res, Fext=Fext)
        f0 = np.real(self.res.dot(stepDir))
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
            self.assembler.setVariables(u)

            # 5. Compute residual, $r = r(u)$
            self.getResidual(self.res, Fext=Fext)

            # 6. Compute merit function,  $f(\alpha)=f(u, r, \Delta u)$
            fNew = np.real(self.res.dot(stepDir))

            # 7. if $abs(f(\alpha)) \leq \mu f_0 + \alpha f'_0$:
            #     1. exit
            uNorm = u.norm()
            if self.rank == 0 and PRINT_LINESEARCH_ITERS:
                print(
                    f"Line search iter {(iteration+1):2d}: alpha = {alpha: 11e}, f/f0 = {(fNew/f0): 11e}, uNorm = {uNorm: 11e}"
                )
            u.axpy(-alpha, stepDir)
            fReduction = np.abs(fNew / f0)
            if fReduction <= 1 - LINESEARCH_MU * min(alpha, 1.0):
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
        return alpha, iteration

    def solveJacLinear(self, res, sol):
        success = self.KSM.solve(res, sol)

        if not success:
            self._TACSWarning(
                "Linear solver failed to converge. "
                "This is likely a sign that the problem is ill-conditioned. "
                "Check that the model is properly restrained."
            )
        return success

    def updateJacobian(self, res=None):
        if self._stiffnessUpdateRequired:
            # Assemble residual and stiffness matrix (w/o artificial terms)
            self.assembler.assembleJacobian(
                self.alpha,
                self.beta,
                self.gamma,
                res,
                self.K,
                loadScale=self._loadScale,
            )
            self._factorOnNext = True

    def updatePreconditioner(self):
        if self._factorOnNext:
            # Stiffness matrix must include artificial terms before pc factor
            # to prevent factorization issues w/ zero-diagonals
            self.K.axpy(1.0, self.rbeArtificialStiffness)
            self.PC.factor()
            # Remove artificial stiffness terms to get true stiffness mat
            self.K.axpy(-1.0, self.rbeArtificialStiffness)
            self._factorOnNext = False

    def checkConvergence(self, resNorm):
        """Check whether the residual is sufficiently converged

        Returns
        -------
        _type_
            _description_
        """
        return resNorm / self.initNorm < 1e-8

    def checkDivergence(self, resNorm):
        """Check whether the residual has diverged

        Returns
        -------
        _type_
            _description_
        """
        return resNorm > 1e10 or np.isnan(resNorm)

    ####### Function eval/sensitivity methods ########

    def evalFunctions(self, funcs, evalFuncs=None, ignoreMissing=False):
        """
        This is the main routine for returning useful information from
        pytacs. The functions corresponding to the strings in
        EVAL_FUNCS are evaluated and updated into the provided
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
        corresponding to the strings in EVAL_FUNCS are evaluated and
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

        if self.getOption("printTiming") and self.rank == 0:
            self._pp("+--------------------------------------------------+")
            self._pp("|")
            self._pp("| TACS Adjoint Times:")
            print("|")
            print(
                "| %-30s: %10.3f sec"
                % ("TACS Sens Setup Problem Time", setupProblemTime - startTime)
            )
            print(
                "| %-30s: %10.3f sec"
                % ("TACS Adjoint RHS Time", adjointRHSTime - setupProblemTime)
            )
            for f in evalFuncs:
                print(
                    "| %-30s: %10.3f sec"
                    % (
                        "TACS Adjoint Solve Time - %s" % (f),
                        adjointEndTime[f] - adjointStartTime[f],
                    )
                )
            print(
                "| %-30s: %10.3f sec"
                % (
                    "Total Sensitivity Time",
                    totalSensitivityTime - adjointFinishedTime,
                )
            )
            print("|")
            print(
                "| %-30s: %10.3f sec"
                % (
                    "Complete Sensitivity Time",
                    totalSensitivityTime - startTime,
                )
            )
            print("+--------------------------------------------------+")

    def addSVSens(self, evalFuncs, svSensList):
        """
        Add the state variable partial sensitivity to the ADjoint RHS for given evalFuncs

        Parameters
        ----------
        evalFuncs : list[str]
            The functions the user wants returned

        svSensList : list[TACS.Vec] or list[numpy.ndarray]
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

        dvSensList : list[BVec] or list[numpy.ndarray]
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
        adjointlist : list[BVec] or list[numpy.ndarray]
            List of adjoint vectors for residual sensitivity product

        dvSensList : list[BVec] or list[numpy.ndarray]
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

        xptSensList : list[BVec] or list[numpy.ndarray]
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
        adjointlist : list[BVec] or list[numpy.ndarray]
            List of adjoint vectors for residual sensitivity product

        xptSensList : list[BVec] or list[numpy.ndarray]
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
        res : TACS BVec or numpy array
            If res is not None, place the residuals into this array.

        Fext : TACS BVec or numpy array, optional
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
        tuple of 2 or 4 scipy.sparse.bsr_matrices
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
        phi : TACS BVec or numpy array
            Input vector to product with the transpose Jacobian.

        prod : TACS BVec or numpy array
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
        rhs : TACS BVec or numpy array
            right hand side vector for adjoint solve
        phi : TACS BVec or numpy array
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
        self.KSM.solve(self.adjRHS, self.phi)
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
        states : TACS.Vec or numpy.ndarray
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

    def getOutputFileName(self, outputDir=None, baseName=None, number=None):
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
        return baseName

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
            base = os.path.join(outputDir, baseName) + ".f5"
            self.outputViewer.writeToFile(base)

    def writeSolverHistory(self, outputDir=None, baseName=None, number=None):
        # Figure out the output file base name
        baseName = self.getOutputFileName(outputDir, baseName, number)
        if self.history is not None:
            self.history.save(baseName)
