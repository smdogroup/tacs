"""
The main purpose of this class is to represent all relevant
information for a linearized buckling analysis.

.. note:: This class should be created using the
    :meth:`pyTACS.createBucklingProblem <tacs.pytacs.pyTACS.createBucklingProblem>` method.
"""

# =============================================================================
# Imports
# =============================================================================
import os
import time

import numpy as np

import tacs.TACS
from .base import TACSProblem


class BucklingProblem(TACSProblem):
    # Default Option List
    defaultOptions = {
        "outputDir": [str, "./", "Output directory for F5 file writer."],
        # Solution Options
        "L2Convergence": [
            float,
            1e-12,
            "Absolute convergence tolerance for Eigenvalue solver based on l2 norm of residual.",
        ],
        "L2ConvergenceRel": [
            float,
            1e-12,
            "Relative convergence tolerance for Eigenvalue solver based on l2 norm of residual.",
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
        "subSpaceSize": [
            int,
            10,
            "Subspace size for Krylov solver used by Eigenvalue solver.",
        ],
        "nRestarts": [
            int,
            15,
            "Max number of resets for Krylov solver used by Eigenvalue solver.",
        ],
        # Output Options
        "writeSolution": [bool, True, "Flag for suppressing all f5 file writing."],
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
        "printLevel": [
            int,
            0,
            "Print level for Eigenvalue solver.\n"
            "\t Accepts:\n"
            "\t\t   0 : No printing.\n"
            "\t\t   1 : Print major iterations.\n"
            "\t\t > 1 : Print major + minor iterations.",
        ],
    }

    def __init__(
        self,
        name,
        sigma,
        numEigs,
        assembler,
        comm,
        outputViewer=None,
        meshLoader=None,
        options=None,
    ):
        """
        NOTE: This class should not be initialized directly by the user.
        Use pyTACS.createBucklingProblem instead.

        Parameters
        ----------
        name : str
            Name of this tacs problem

        sigma : float
            Guess for the lowest eigenvalue. This corresponds to the lowest buckling load factor.

        numEigs : int
            Number of eigenvalues to solve for

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

        # Default setup for common problem class objects, sets up comm and options
        TACSProblem.__init__(self, assembler, comm, options, outputViewer, meshLoader)

        # Set time eigenvalue parameters
        self.sigma = sigma
        self.numEigs = numEigs

        # String name used in evalFunctions
        self.valName = "eigsb"
        self._initializeFunctionList()

        # Create problem-specific variables
        self._createVariables()

    def _createVariables(self):
        """
        Internal to create the objects required by TACS Integrator
        """

        self.callCounter = -1

        # Buckling load state
        self.u0 = self.assembler.createVec()
        # Load vector
        self.F = self.assembler.createVec()
        self.F_array = self.F.getArray()
        # RHS vector
        self.rhs = self.assembler.createVec()
        # Auxiliary element object for applying tractions/pressure
        self.auxElems = tacs.TACS.AuxElements()

        self.aux = self.assembler.createSchurMat()
        self.G = self.assembler.createSchurMat()
        self.K = self.assembler.createSchurMat()

        self.pc = tacs.TACS.Pc(self.K)

        # Set artificial stiffness factors in rbe class
        c1 = self.getOption("RBEStiffnessScaleFactor")
        c2 = self.getOption("RBEArtificialStiffness")
        tacs.elements.RBE2.setScalingParameters(c1, c2)
        tacs.elements.RBE3.setScalingParameters(c1, c2)

        # Assemble and factor the stiffness/Jacobian matrix. Factor the
        # Jacobian and solve the linear system for the displacements
        self.assembler.assembleMatType(tacs.TACS.STIFFNESS_MATRIX, self.K)
        self.assembler.assembleMatType(tacs.TACS.GEOMETRIC_STIFFNESS_MATRIX, self.G)

        subspace = self.getOption("subSpaceSize")
        restarts = self.getOption("nRestarts")
        atol = self.getOption("L2Convergence")
        rtol = self.getOption("L2ConvergenceRel")
        self.gmres = tacs.TACS.KSM(self.aux, self.pc, subspace, restarts)
        self.gmres.setTolerances(rtol, atol)

        # Create the buckling analysis object
        self.buckleSolver = tacs.TACS.BucklingAnalysis(
            self.assembler,
            self.sigma,
            self.G,
            self.K,
            self.gmres,
            num_eigs=self.numEigs,
            eig_tol=rtol,
        )

    def _initializeFunctionList(self):
        """
        Create FunctionList dict which maps eigenvalue strings
        to mode indices used in evalFunctions method.
        """
        self.functionList = {}
        for mode_i in range(self.numEigs):
            self.functionList[f"{self.valName}.{mode_i}"] = mode_i

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

        # No need to reset solver for output options
        if name.lower() in [
            "writesolution",
            "printtiming",
            "numbersolutions",
            "outputdir",
        ]:
            pass
        # Reset solver for all other option changes
        else:
            self._createVariables()

    def setValName(self, valName):
        """
        Set a name for the eigenvalues. Only needs
        to be changed if more than 1 pytacs object is used in an
        optimization

        Parameters
        ----------
        valName : str
            Name of the eigenvalue output used in evalFunctions().
        """
        self.valName = valName
        self._initializeFunctionList()

    def getNumEigs(self):
        """
        Get the number of eigenvalues requested from solver for this problem.

        Returns
        ----------
        numEigs : int
            Number of eigenvalues.
        """
        return self.numEigs

    def addFunction(self, funcName, funcHandle, compIDs=None, **kwargs):
        """
        NOT SUPPORTED FOR THIS PROBLEM
        """
        self._TACSWarning("addFunction method not supported for this class.")

    def evalFunctions(self, funcs, evalFuncs=None, ignoreMissing=False):
        """
        Evaluate eigenvalues for problem. The functions corresponding to
        the integers in evalFuncs are evaluated and updated into
        the provided dictionary.

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
        >>> bucklingProblem.solve()
        >>> bucklingProblem.evalFunctions(funcs, 'eigsm.0')
        >>> funcs
        >>> # Result will look like (if bucklingProblem has name of 'c1'):
        >>> # {'c1_eigsm.0':12354.10}
        """
        # Make sure assembler variables are up to date
        self._updateAssemblerVars()

        # Check if user specified which eigvals to output
        # Otherwise, output them all
        if evalFuncs is None:
            evalFuncs = self.functionList
        else:
            userFuncs = sorted(list(evalFuncs))
            evalFuncs = {}
            for func in userFuncs:
                if func in self.functionList:
                    evalFuncs[func] = self.functionList[func]

        if not ignoreMissing:
            for f in evalFuncs:
                if f not in self.functionList:
                    raise self._TACSError(
                        f"Supplied function '{f}' has not been added "
                        "using addFunction()."
                    )

        # Loop through each requested eigenvalue
        for funcName in evalFuncs:
            mode_i = evalFuncs[funcName]
            key = f"{self.name}_{funcName}"
            funcs[key], _ = self.getVariables(mode_i)

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
        >>> bucklingProblem.evalFunctionsSens(funcsSens, 'eigsm.0')
        >>> funcsSens
        >>> # Result will look like (if bucklingProblem has name of 'c1'):
        >>> # {'c1_eigsm.0':{'struct':[1.234, ..., 7.89], 'Xpts':[3.14, ..., 1.59]}}
        """
        # Make sure assembler variables are up to date
        self._updateAssemblerVars()

        # Check if user specified which eigvals to output
        # Otherwise, output them all
        if evalFuncs is None:
            evalFuncs = self.functionList
        else:
            userFuncs = sorted(list(evalFuncs))
            evalFuncs = {}
            for func in userFuncs:
                if func in self.functionList:
                    evalFuncs[func] = self.functionList[func]

        dvSens = self.assembler.createDesignVec()
        xptSens = self.assembler.createNodeVec()

        indices = [evalFuncs[funcName] for funcName in evalFuncs]
        dvSensList = [self.assembler.createDesignVec() for funcName in evalFuncs]
        xptSensList = [self.assembler.createNodeVec() for funcName in evalFuncs]
        svSensList = [self.assembler.createVec() for funcName in evalFuncs]
        adjointList = [self.assembler.createVec() for funcName in evalFuncs]

        self.addDVSens(indices, dvSensList, scale=1.0)
        self.addXptSens(indices, xptSensList, scale=1.0)
        self.evalSVSens(indices, svSensList)

        self.aux.copyValues(self.K)
        self.pc.factor()

        # Loop through each requested eigenvalue
        for i, funcName in enumerate(evalFuncs):
            rhs = svSensList[i]
            adjoint = adjointList[i]
            self.gmres.solve(rhs, adjoint)

            # Evaluate adjoint contribution to nodal sens
            xptSens = xptSensList[i]
            self.assembler.addAdjointResXptSensProducts([adjoint], [xptSens], -1.0)
            xptSens.beginSetValues()
            xptSens.endSetValues()
            # Evaluate adjoint contribution to dv sens
            dvSens = dvSensList[i]
            self.assembler.addAdjointResProducts([adjoint], [dvSens], -1.0)
            dvSens.beginSetValues()
            dvSens.endSetValues()

            key = f"{self.name}_{funcName}"
            funcsSens[key] = {}
            funcsSens[key][self.varName] = dvSens.getArray().copy()
            funcsSens[key][self.coordName] = xptSens.getArray().copy()

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

    def zeroLoads(self):
        """
        Zero all applied loads
        """
        self.F.zeroEntries()
        self.auxElems = tacs.TACS.AuxElements()

    ####### Buckling solver methods ########

    def _updateAssemblerVars(self):
        """
        Make sure that the assembler is using
        the input variables associated with this problem
        """

        self.assembler.setDesignVars(self.x)
        self.assembler.setNodes(self.Xpts)
        # Set state variables
        self.assembler.setVariables(self.u0)
        # Zero any time derivative terms
        self.assembler.zeroDotVariables()
        self.assembler.zeroDDotVariables()
        # Make sure previous auxiliary loads are removed
        self.assembler.setAuxElements(self.auxElems)
        # Set artificial stiffness factors in rbe class
        c1 = self.getOption("RBEStiffnessScaleFactor")
        c2 = self.getOption("RBEArtificialStiffness")
        tacs.elements.RBE2.setScalingParameters(c1, c2)
        tacs.elements.RBE3.setScalingParameters(c1, c2)

    def solve(self, Fext=None, u0=None):
        """
        Solve the eigenvalue problem. The
        forces must already be set.

        Parameters
        ----------
        Fext : tacs.TACS.Vec or numpy.ndarray, optional
            Distributed array containing additional loads (ex. aerodynamic forces for aerostructural coupling)
            to applied to RHS of the static problem.

        u0 : tacs.TACS.Vec or numpy.ndarray, optional
            Distributed array containing additional loads (ex. aerodynamic forces for aerostructural coupling)
            to applied to RHS of the static problem. Alternate to Fext.
        """
        startTime = time.time()

        self.callCounter += 1

        setupProblemTime = time.time()

        # Set problem vars to assembler
        self._updateAssemblerVars()

        # states were not prescribed, pass in forces
        if u0 is None:
            # Sum the forces from the loads not handled by TACS
            self.rhs.copyValues(self.F)  # Fixed loads

            # Add external loads, if specified
            if Fext is not None:
                if isinstance(Fext, tacs.TACS.Vec):
                    self.rhs.axpy(1.0, Fext)
                elif isinstance(Fext, np.ndarray):
                    rhsArray = self.rhs.getArray()
                    rhsArray[:] = rhsArray[:] + Fext[:]

            force = self.rhs
            path = None

        # states already prescribed, we don't need forces
        else:
            # Convert to bvec
            if isinstance(u0, np.ndarray):
                path = self._arrayToVec(u0)
            else:
                path = u0
            force = None

        initSolveTime = time.time()

        # Solve the buckling analysis problem
        self.buckleSolver.solve(
            force=force, path=path, print_flag=self.getOption("printLevel")
        )

        # Save state vars
        self.assembler.getVariables(self.u0)

        solveTime = time.time()

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
                "| %-30s: %10.3f sec" % ("TACS Solve Time", solveTime - initSolveTime)
            )
            self._pp("|")
            self._pp(
                "| %-30s: %10.3f sec"
                % ("TACS Total Solution Time", solveTime - startTime)
            )
            self._pp("+--------------------------------------------------+")

        return

    def getVariables(self, index, states=None):
        """
        Return the current state values for one mode of the current problem

        Parameters
        ----------
        index : int
            Mode index to return solution for.

        states : tacs.TACS.Vec or numpy.ndarray or None
            Place eigenvector for mode into this array (optional).

        Returns
        --------
        eigVal: float
            Eigenvalue for mode corresponds to buckling load factor

        states : numpy.ndarray
            Eigenvector for mode
        """
        eigVal, err = self.buckleSolver.extractEigenvalue(index)
        eigVector = self.assembler.createVec()
        self.buckleSolver.extractEigenvector(index, eigVector)
        # Inplace assignment if vectors were provided
        if isinstance(states, tacs.TACS.Vec):
            states.copyValues(eigVector)
        elif isinstance(states, np.ndarray):
            states[:] = eigVector.getArray()
        return eigVal, eigVector.getArray()

    def addXptSens(self, indices, xptSensList, scale=1.0):
        """
        Add partial sensitivity contribution due to nodal coordinates for eigenvalue

        Parameters
        ----------
        indices : list[int]
            Mode indices to return sensitivity for.

        xptSensList : list[tacs.TACS.Vec] or list[numpy.ndarray]
            List of sensitivity vectors to add partial sensitivity to

        scale : float
            Scalar to multiply partial sensitivity by. Defaults to 1.0
        """
        # Set problem vars to assembler
        self._updateAssemblerVars()

        for index, xptSens in zip(indices, xptSensList):
            # Create a tacs BVec copy for the operation if the output is a numpy array
            if isinstance(xptSens, np.ndarray):
                xptSensBVec = self._arrayToNodeVec(xptSens)
            # Otherwise the input is already a BVec and we can do the operation in place
            else:
                xptSensBVec = xptSens

            self.buckleSolver.addEigenXptSens(scale, index, xptSensBVec)

            # Finalize sensitivity arrays across all procs
            xptSensBVec.beginSetValues()
            xptSensBVec.endSetValues()

            # Update from the BVec values, if the input was a numpy array
            if isinstance(xptSens, np.ndarray):
                # Copy values to numpy array
                xptSens[:] = xptSensBVec.getArray()

    def addDVSens(self, indices, dvSensList, scale=1.0):
        """
        Add partial sensitivity contribution due to design variables for eigenvalue

        Parameters
        ----------
        indices : list[int]
            Mode indices to return sensitivity for.

        dvSensList : list[tacs.TACS.Vec] or list[numpy.ndarray]
            List of sensitivity vectors to add partial sensitivity to

        scale : float
            Scalar to multiply partial sensitivity by. Defaults to 1.0
        """
        # Set problem vars to assembler
        self._updateAssemblerVars()

        for index, dvSens in zip(indices, dvSensList):
            # Create a tacs BVec copy for the operation if the output is a numpy array
            if isinstance(dvSens, np.ndarray):
                dvSensBVec = self._arrayToDesignVec(dvSens)
            # Otherwise the input is already a BVec and we can do the operation in place
            else:
                dvSensBVec = dvSens

            self.buckleSolver.addEigenDVSens(scale, index, dvSensBVec)

            # Finalize sensitivity arrays across all procs
            dvSensBVec.beginSetValues()
            dvSensBVec.endSetValues()

            # Update from the BVec values, if the input was a numpy array
            if isinstance(dvSens, np.ndarray):
                # Copy values to numpy array
                dvSens[:] = dvSensBVec.getArray()

    def evalSVSens(self, indices, svSensList):
        """
        Add partial sensitivity contribution due to state variables for eigenvalue

        Parameters
        ----------
        indices : list[int]
            Mode indices to return sensitivity for.

        svSensList : list[tacs.TACS.Vec] or list[numpy.ndarray]
            List of sensitivity vectors to add partial sensitivity to
        """
        # Set problem vars to assembler
        self._updateAssemblerVars()

        for index, svSens in zip(indices, svSensList):
            # Create a tacs BVec copy for the operation if the output is a numpy array
            if isinstance(svSens, np.ndarray):
                svSensBVec = self._arrayToVec(svSens)
            # Otherwise the input is already a BVec and we can do the operation in place
            else:
                svSensBVec = svSens

            self.buckleSolver.evalEigenSVSens(index, svSensBVec)

            # Finalize sensitivity arrays across all procs
            svSensBVec.beginSetValues()
            svSensBVec.endSetValues()

            # Update from the BVec values, if the input was a numpy array
            if isinstance(svSens, np.ndarray):
                # Copy values to numpy array
                svSens[:] = svSensBVec.getArray()

    def writeSolution(self, outputDir=None, baseName=None, number=None, indices=None):
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
        indices : int or list[int] or None
            Mode index or indices to get state variables for.
            If None, returns a solution for all modes.
            Defaults to None.
        """
        # Make sure assembler variables are up-to-date
        self._updateAssemblerVars()

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

        # Unless the writeSolution option is off write actual file:
        if self.getOption("writeSolution"):
            # If indices is None, output all modes
            if indices is None:
                indices = np.arange(self.numEigs)

            # Write out each specified mode
            indices = np.atleast_1d(indices)
            vec = self.assembler.createVec()
            for index in indices:
                # Extract eigenvector
                eig, _ = self.getVariables(index, states=vec)
                # Set eigen mode in assembler
                self.assembler.setVariables(vec)
                # Write out mode shape as f5 file
                modeName = baseName + "_%3.3d" % index
                fileName = os.path.join(outputDir, modeName) + ".f5"
                self.outputViewer.writeToFile(fileName)
