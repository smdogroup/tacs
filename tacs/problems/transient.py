"""
The main purpose of this class is to represent all relevant
information for a transient analysis. This will include
information defining the loading condition as well as various
other pieces of information.

.. note:: This class should be created using the
    :meth:`pyTACS.createTransientProblem <tacs.pytacs.pyTACS.createTransientProblem>` method.
"""

# =============================================================================
# Imports
# =============================================================================
import os
import numpy as np
import time
from .base import TACSProblem
import tacs.TACS


class TransientProblem(TACSProblem):

    # Default Option List
    defaultOptions = {
        'outputdir': [str, './', 'Output directory for F5 file writer.'],

        # Solution Options
        'timeIntegrator': [str, 'BDF', "Time integration scheme to use. Currently supports 'BDF' and 'DIRK'."],
        'integrationOrder': [int, 2, "Integration order for time marching scheme."],
        'L2Convergence': [float, 1e-12, 'Absolute convergence tolerance for integrator based on l2 norm of residual.'],
        'L2ConvergenceRel': [float, 1e-12,
                             'Relative convergence tolerance for integrator based on l2 norm of residual.'],
        'jacAssemblyFreq': [int, 1, 'How frequently to reassemble Jacobian during time integration process.'],

        # Output Options
        'writeSolution': [bool, True, 'Flag for suppressing all f5 file writing.'],
        'numberSolutions': [bool, True, 'Flag for attaching solution counter index to f5 files.'],
        'printTiming': [bool, False, 'Flag for printing out timing information for class procedures.'],
        'printLevel': [int, 0, 'Print level for integration solver.\n'
                               '\t Accepts:\n'
                               '\t\t   0 : No printing.\n'
                               '\t\t   1 : Print major iterations.\n'
                               '\t\t > 1 : Print major + minor iterations.'],

    }

    def __init__(self, name, tInit, tFinal, numSteps,
                 assembler, comm, outputViewer=None, meshLoader=None,
                 options={}):
        """
        NOTE: This class should not be initialized directly by the user.
        Use pyTACS.createTransientProblem instead.

        Parameters
        ----------
        name : str
            Name of this tacs problem

        tInit : float
            Starting time for transient problem integration

        tFinal : float
            Ending time for transient problem integration

        numSteps : int
            Number of time steps for transient problem integration

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

        # Defualt setup for common problem class objects
        super().__init__(assembler, comm, outputViewer, meshLoader)

        # Set time interval parmeters
        self.tInit = tInit
        self.tFinal = tFinal
        self.numSteps = numSteps

        # Process the default options which are added to self.options
        # under the 'defaults' key. Make sure the key are lower case
        def_keys = self.defaultOptions.keys()
        self.options['defaults'] = {}
        for key in def_keys:
            self.options['defaults'][key.lower()] = self.defaultOptions[key]
            self.options[key.lower()] = self.defaultOptions[key]

        # Set user-defined options
        for key in options:
            super().setOption(key, options[key])

        # Create problem-specific variables
        self._createVariables()

    def _createVariables(self):
        """
        Internal to create the objects required by TACS Integrator
        """

        self.callCounter = -1

        # Create a force vector for each time step
        self.F = [self.assembler.createVec() for i in range(self.numSteps + 1)]
        # Auxillary element object for applying tractions/pressure
        self.auxElems = [tacs.TACS.AuxElements() for i in range(self.numSteps + 1)]
        # Initialize the initial conditions tacs vectors
        self.vars0 = self.assembler.createVec()
        self.dvars0 = self.assembler.createVec()
        self.ddvars0 = self.assembler.createVec()

        # Create the BDF integrator solver
        order = self.getOption('integrationOrder')
        solverType = self.getOption('timeIntegrator')
        # Chose solver type
        if solverType.upper() == 'BDF':
            self.integrator = tacs.TACS.BDFIntegrator(self.assembler, self.tInit, self.tFinal,
                                                      float(self.numSteps), order)
        elif solverType.upper() == 'DIRK':
            self.integrator = tacs.TACS.DIRKIntegrator(self.assembler, self.tInit, self.tFinal,
                                                       float(self.numSteps), order)

        printLevel = self.getOption('printLevel')
        self.integrator.setPrintLevel(printLevel)
        # Set solver tolerances
        atol = self.getOption('L2Convergence')
        self.integrator.setAbsTol(atol)
        rtol = self.getOption('L2ConvergenceRel')
        self.integrator.setRelTol(rtol)
        # Jacobian assembly frequency
        jacFreq = self.getOption('jacAssemblyFreq')
        self.integrator.setJacAssemblyFreq(jacFreq)

        # Set output viewer for integrator
        self.integrator.setFH5(self.outputViewer)
        outputDir = self.getOption('outputDir')
        self.integrator.setOutputPrefix(outputDir)

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
        # Defualt setOption for common problem class objects
        super().setOption(name, value)

        # Update tolerances
        if 'l2convergence' in name.lower():
            # Set solver tolerances
            atol = self.getOption('L2Convergence')
            self.integrator.setAbsTol(atol)
            rtol = self.getOption('L2ConvergenceRel')
            self.integrator.setRelTol(rtol)
        elif name.lower() == 'printlevel':
            printLevel = self.getOption('printLevel')
            self.integrator.setPrintLevel(printLevel)
        elif name.lower() == 'jacassemblyfreq':
            # Jacobian assembly frequency
            jacFreq = self.getOption('jacAssemblyFreq')
            self.integrator.setJacAssemblyFreq(jacFreq)
        # No need to reset solver for output options
        elif name.lower() in ['writesolution', 'printtiming',
                              'numbersolutions', 'outputdir']:
            pass
        # Reset solver for all other option changes
        else:
            self._createVariables()

    def getNumTimeSteps(self):
        """
        Get the number of timesteps used in time integration for this problem.

        Returns
        ----------
        numSteps : int
            Number of time steps.
        """
        # Should this really be numSteps + 1 ?
        return self.numSteps

    def getTimeSteps(self):
        """
        Get the discrete time step slices used in time integration.

        Returns
        ----------
        numSteps : int
            Number of time steps.
        """
        timeSteps = np.linspace(self.tInit, self.tFinal, self.numSteps + 1)
        return timeSteps

    ####### Load adding methods ########

    def addLoadToComponents(self, timeStep, compIDs, F, averageLoad=False):
        """"
        This method is used to add a *FIXED TOTAL LOAD* on one or more
        components, defined by COMPIDs, at a specifc time instance.
        The purpose of this routine is to add loads that remain fixed throughout
        an optimization. An example would be an engine load. This routine determines
        all the unqiue nodes in the FE model that are part of the the requested components,
        then takes the total 'force' by F and divides by the number of nodes.
        This average load is then applied to the nodes.

        Parameters
        ----------

        timeStep : int
            Time step index to apply load to.

        compIDs : list[int] or int
            The components with added loads. Use pyTACS.selectCompIDs method
            to determine this.

        F : Numpy 1d or 2d array length (varsPerNodes) or (numNodeIDs, varsPerNodes)
            Vector(s) of 'force' to apply to each components.  If only one force vector is provided,
            force will be copied uniformly across all components.

        averageLoad : bool
            Flag to determine whether load should be split evenly across all components (True)
            or copied and applied individually to each component (False). Defaults to False.

        Notes
        -----

        The units of the entries of the 'force' vector F are not
        necesarily physical forces and their interpretation depends
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
        self._addLoadToComponents(self.F[timeStep], compIDs, F, averageLoad)

    def addLoadToNodes(self, timeStep, nodeIDs, F, nastranOrdering=False):
        """
        This method is used to add a fixed point load of F to the
        selected node IDs at a specified time instance.

        Parameters
        ----------

        timeStep : int
            Time step index to apply load to.

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
        necesarily physical forces and their interpretation depends
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

        self._addLoadToNodes(self.F[timeStep], nodeIDs, F, nastranOrdering)

    def addLoadToRHS(self, timeStep, Fapplied):
        """"
        This method is used to add a *FIXED TOTAL LOAD* directly to the
        right hand side vector given the equation below:

            M*udotdot + K*u = f

        Where:
            - K : Stiffness matrix for problem
            - u : State variables for problem
            - M : Mass matrix for problem
            - udotdot : Second time derivitive of state variables for problem
            - f : Right-hand side vector to add loads to

        Parameters
        ----------

        timeStep : int
            Time step index to apply load to.

        Fapplied : numpy.ndarray or TACS.Vec
            Distributed array containing loads to applied to RHS of the problem.

        """
        self._addLoadToRHS(self.F[timeStep], Fapplied)

    def addTractionToComponents(self, timeStep, compIDs, tractions,
                                faceIndex=0):
        """
        This method is used to add a *FIXED TOTAL TRACTION* on one or more
        components, defined by COMPIDs, at specified time instance. The purpose of
        this routine is to add loads that remain fixed throughout an optimization.

        Parameters
        ----------

        timeStep : int
            Time step index to apply load to.

        compIDs : list[int] or int
            The components with added loads. Use pyTACS.selectCompIDs method
            to determine this.

        tractions : numpy.ndarray length 1 or compIDs
            Array of traction vectors for each components

        faceIndex : int
            Indicates which face (side) of element to apply traction to.
            Note: not required for certain elements (i.e. shells)
        """
        self._addTractionToComponents(self.auxElems[timeStep], compIDs, tractions, faceIndex)

    def addTractionToElements(self, timeStep, elemIDs, tractions,
                              faceIndex=0, nastranOrdering=False):
        """
        This method is used to add a fixed traction to the
        selected element IDs at specified time instance.
        Tractions can be specified on an element by element basis
        (if tractions is a 2d array) or set to a uniform value (if tractions is a 1d array)

        Parameters
        ----------

        timeStep : int
            Time step index to apply load to.

        elemIDs : list[int]
            The global element ID numbers for which to apply the traction.

        tractions : numpy.ndarray 1d or 2d length varsPerNodes or (elemIDs, varsPerNodes)
            Array of traction vectors for each element

        faceIndex : int
            Indicates which face (side) of element to apply traction to.
            Note: not required for certain elements (i.e. shells)

        nastranOrdering : bool
            Flag signaling whether elemIDs are in TACS (default)
            or NASTRAN ordering
        """

        self._addTractionToElements(self.auxElems[timeStep], elemIDs, tractions, faceIndex, nastranOrdering)

    def addPressureToComponents(self, timeStep, compIDs, pressures,
                                faceIndex=0):
        """
        This method is used to add a *FIXED TOTAL PRESSURE* on one or more
        components, defined by COMPIDs, at specified time instance. The purpose of this routine is
        to add loads that remain fixed throughout an optimization. An example
        would be a fuel load.

        Parameters
        ----------

        timeStep : int
            Time step index to apply load to.

        compIDs : list[int] or int
            The components with added loads. Use pyTACS.selectCompIDs method
            to determine this.

        pressures : Numpy array length 1 or compIDs
            Array of pressure values for each components

        faceIndex : int
            Indicates which face (side) of element to apply pressure to.
            Note: not required for certain elements (i.e. shells)
        """
        self._addPressureToComponents(self.auxElems[timeStep], compIDs, pressures, faceIndex)

    def addPressureToElements(self, timeStep, elemIDs, pressures,
                              faceIndex=0, nastranOrdering=False):
        """
        This method is used to add a fixed presure to the
        selected element IDs at specified time instance.
        Pressures can be specified on an element by element
        basis (if pressures is an array) or set to a uniform value (if pressures is a scalar)

        Parameters
        ----------

        timeStep : int
            Time step index to apply load to.

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

        self._addPressureToElements(self.auxElems[timeStep], elemIDs, pressures,
                                    faceIndex, nastranOrdering)

    def addInertialLoad(self, timeStep, inertiaVector):
        """
        This method is used to add a fixed inertial load  at a specified time step
        due to a uniform acceleration over the entire model.
        This is most commonly used to model gravity loads on a model.

        Parameters
        ----------

        timeStep : int
            Time step index to apply load to.

        inertiaVector : numpy.ndarray
            Acceleration vector used to define inertial load.
        """
        self._addInertialLoad(self.auxElems[timeStep], inertiaVector)

    ####### Transient solver methods ########

    def setInitConditions(self, vars=None, dvars=None, ddvars=None):
        """
        Set the initial conditions associated with this problem

        Parameters
        ----------
        vars : float or numpy.ndarray or TACS.Vec
            Initial conditions of the state variables
        dvars : float or numpy.ndarray or TACS.Vec
            Initial conditions of the first time-derivative of the state variables
        ddvars : float or numpy.ndarray or TACS.Vec
            Initial conditions of the second time-derivative of the state variables
        """

        if vars is not None:
            if isinstance(vars, np.ndarray):
                vars0Array = self.vars0.getArray()
                vars0Array[:] = vars[:]
            elif isinstance(vars, tacs.TACS.Vec):
                self.vars0.copyValues(vars)
            else:  # assume type=float
                vars0Array = self.vars0.getArray()
                vars0Array[:] = vars

        if dvars is not None:
            if isinstance(dvars, np.ndarray):
                dvars0Array = self.dvars0.getArray()
                dvars0Array[:] = dvars[:]
            elif isinstance(dvars, tacs.TACS.Vec):
                self.dvars0.copyValues(dvars)
            else:  # assume type=float
                dvars0Array = self.dvars0.getArray()
                dvars0Array[:] = dvars

        if ddvars is not None:
            if isinstance(ddvars, np.ndarray):
                ddvars0Array = self.ddvars0.getArray()
                ddvars0Array[:] = ddvars[:]
            elif isinstance(ddvars, tacs.TACS.Vec):
                self.ddvars0.copyValues(ddvars)
            else:  # assume type=float
                ddvars0Array = self.ddvars0.getArray()
                ddvars0Array[:] = ddvars

    def _updateAssemblerVars(self):
        """
        Make sure that the assembler is using
        the input variables associated with this problem
        """

        self.assembler.setDesignVars(self.x)
        self.assembler.setNodes(self.Xpts)
        self.assembler.setInitConditions(vec=self.vars0, dvec=self.dvars0, ddvec=self.ddvars0)

    def solve(self):
        """
        Solve the time integrated transient problem. The
        forces must already be set.
        """
        startTime = time.time()

        self.callCounter += 1

        setupProblemTime = time.time()

        # Set problem vars to assembler
        self._updateAssemblerVars()

        initSolveTime = time.time()

        # Loop over every time step and solve transient problem
        for i in range(self.numSteps + 1):
            # Set the auxiliary elements for this time step (tractions/pressures)
            self.assembler.setAuxElements(self.auxElems[i])
            self.integrator.iterate(i, forces=self.F[i])

        solveTime = time.time()

        # If timing was was requested print it, if the solution is nonlinear
        # print this information automatically if printTiming was requested.
        if self.getOption('printTiming'):
            self._pp('+--------------------------------------------------+')
            self._pp('|')
            self._pp('| TACS Solve Times:')
            self._pp('|')
            self._pp('| %-30s: %10.3f sec' % ('TACS Setup Time', setupProblemTime - startTime))
            self._pp('| %-30s: %10.3f sec' % ('TACS Solve Init Time', initSolveTime - setupProblemTime))
            self._pp('| %-30s: %10.3f sec' % ('TACS Solve Time', solveTime - initSolveTime))
            self._pp('|')
            self._pp('| %-30s: %10.3f sec' % ('TACS Total Solution Time', solveTime - startTime))
            self._pp('+--------------------------------------------------+')

        return

    ####### Function eval/sensitivity methods ########

    def evalFunctions(self, funcs, evalFuncs=None,
                      ignoreMissing=False):
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
        >>> transientProblem.solve()
        >>> transientProblem.evalFunctions(funcs, ['mass'])
        >>> funcs
        >>> # Result will look like (if TransientProblem has name of 'c1'):
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
                    raise self._TACSError(f"Supplied function '{f}' has not been added "
                                          "using addFunction().")

        setupProblemTime = time.time()

        # Fast parallel function evaluation of structural funcs:
        handles = [self.functionList[f] for f in evalFuncs if
                   f in self.functionList]
        # Set functions for integrator
        self.integrator.setFunctions(handles)
        # Evaluate functions
        funcVals = self.integrator.evalFunctions(handles)

        functionEvalTime = time.time()

        # Assign function values to appropriate dictionary
        i = 0
        for f in evalFuncs:
            if f in self.functionList:
                key = self.name + '_%s' % f
                funcs[key] = funcVals[i]
                i += 1

        dictAssignTime = time.time()

        if self.getOption('printTiming'):
            self._pp('+--------------------------------------------------+')
            self._pp('|')
            self._pp('| TACS Function Times:')
            self._pp('|')
            self._pp('| %-30s: %10.3f sec' % ('TACS Function Setup Time', setupProblemTime - startTime))
            self._pp('| %-30s: %10.3f sec' % ('TACS Function Eval Time', functionEvalTime - setupProblemTime))
            self._pp('| %-30s: %10.3f sec' % ('TACS Dict Time', dictAssignTime - functionEvalTime))
            self._pp('|')
            self._pp('| %-30s: %10.3f sec' % ('TACS Function Time', dictAssignTime - startTime))
            self._pp('+--------------------------------------------------+')

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
        >>> transientProblem.evalFunctionsSens(funcsSens, ['mass'])
        >>> funcs
        >>> # Result will look like (if TransientProblem has name of 'c1'):
        >>> # {'c1_mass':{'struct':[1.234, ..., 7.89], 'Xpts':[3.14, ..., 1.59]}}
        """

        startTime = time.time()

        # Set problem vars to assembler
        self._updateAssemblerVars()

        if evalFuncs is None:
            evalFuncs = sorted(list(self.functionList))
        else:
            evalFuncs = sorted(list(evalFuncs))

        for f in evalFuncs:
            if f not in self.functionList:
                raise self._TACSError("Supplied function has not been added "
                                      "using addFunction()")

        # Fast parallel function evaluation of structural funcs:
        handles = [self.functionList[f] for f in evalFuncs if
                   f in self.functionList]
        # Set functions for integrator
        self.integrator.setFunctions(handles)

        # integrate adjoints backwards in time from step = numSteps
        # to step = 0
        for i in range(self.numSteps, -1, -1):
            self.assembler.setAuxElements(self.auxElems[i])
            self.integrator.initAdjoint(i)
            self.integrator.iterateAdjoint(i)
            self.integrator.postAdjoint(i)

        adjointFinishedTime = time.time()

        # Recast sensitivities into dict for user
        for i, f in enumerate(evalFuncs):
            key = self.name + '_%s' % f
            # Finalize sensitivity arrays across all procs
            dvSens = self.integrator.getGradient(i)
            dvSens.beginSetValues()
            dvSens.endSetValues()
            xptSens = self.integrator.getXptGradient(i)
            xptSens.beginSetValues()
            xptSens.endSetValues()
            # Return sensitivities as array in sens dict
            funcsSens[key] = {self.varName: dvSens.getArray().copy(),
                              self.coordName: xptSens.getArray().copy()}

        totalSensitivityTime = time.time()

        if self.getOption('printTiming'):
            self._pp('+--------------------------------------------------+')
            self._pp('|')
            self._pp('| TACS Adjoint Times:')
            print('|')
            print('| %-30s: %10.3f sec' % ('Adjoint solve time', adjointFinishedTime - startTime))
            print('|')
            print('| %-30s: %10.3f sec' % ('Complete Sensitivity Time', totalSensitivityTime - startTime))
            print('+--------------------------------------------------+')

    ####### Post processing methods ########

    def getVariables(self, timeStep, states=None, dstates=None, ddstates=None):
        """
        Return the current state values for the current problem

        Parameters
        ----------
        timeStep : int
            Time step index to get state variables for.

        states : TACS.Vec or numpy.ndarray or None
            If states is not None, place the state variables into this array (optional).

        dstates : TACS.Vec or numpy.ndarray or None
            If dstates is not None, place the time derivitive of the state variables into this array (optional).

        ddstates : TACS.Vec or numpy.ndarray or None
            If ddstates is not None, place the second time derivitive of the state variables into this array (optional).

        Returns
        --------
        time: float
            The time at specified step

        states : TACS.Vec or numpy.ndarray
            The state variables.

        dstates : TACS.Vec or numpy.ndarray or None
            The time derivitive of the state variables.

        ddstates : TACS.Vec or numpy.ndarray or None
            The second time derivitive of the state variables.

        """

        # Get BVecs for time instance
        time, q, qdot, qddot = self.integrator.getStates(timeStep)

        # Convert to arrays
        qArray = q.getArray()
        qdotArray = qdot.getArray()
        qddotArray = qddot.getArray()

        # Inplace assignment if vectors were provided
        if isinstance(states, tacs.TACS.Vec):
            states.copyValues(q)
        elif isinstance(states, np.ndarray):
            states[:] = qArray

        if isinstance(dstates, tacs.TACS.Vec):
            dstates.copyValues(qdot)
        elif isinstance(dstates, np.ndarray):
            dstates[:] = qdotArray

        if isinstance(ddstates, tacs.TACS.Vec):
            ddstates.copyValues(qddot)
        elif isinstance(ddstates, np.ndarray):
            ddstates[:] = qddotArray

        # Return arrays
        return time, qArray, qdotArray, qddotArray

    def writeSolution(self, outputDir=None, baseName=None, number=None, timeSteps=None):
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
            Use the user spplied number to index solution. Again, only
            typically used from an external solver
        timeSteps : int or list[int] or None
            Time step index or indices to get state variables for.
            If None, returns a solution for all time steps.
            Defaults to None.
        """
        # Make sure assembler variables are up to date
        self._updateAssemblerVars()

        # Check input
        if outputDir is None:
            outputDir = self.getOption('outputDir')

        if baseName is None:
            baseName = self.name

        # If we are numbering solution, it saving the sequence of
        # calls, add the call number
        if number is not None:
            # We need number based on the provided number:
            baseName = baseName + '_%3.3d' % number
        else:
            # if number is none, i.e. standalone, but we need to
            # number solutions, use internal counter
            if self.getOption('numberSolutions'):
                baseName = baseName + '_%3.3d' % self.callCounter

        # Unless the writeSolution option is off write actual file:
        if self.getOption('writeSolution'):

            # If timeSteps is None, output all modes
            if timeSteps is None:
                timeSteps = np.arange(self.numSteps + 1)

            # Write out each specified timestep
            timeSteps = np.atleast_1d(timeSteps)
            vec = self.assembler.createVec()
            for timeStep in timeSteps:
                # Extract eigenvector
                self.getVariables(timeStep, states=vec)
                # Set eigen mode in assembler
                self.assembler.setVariables(vec)
                # Write out mode shape as f5 file
                modeName = baseName + '_%3.3d' % timeStep
                fileName = os.path.join(outputDir, modeName) + '.f5'
                self.outputViewer.writeToFile(fileName)
