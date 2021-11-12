"""
pyStatic_problem
"""

# =============================================================================
# Imports
# =============================================================================
import warnings
import os
import numpy as np
from collections import OrderedDict
import time
from .base import BaseProblem
import tacs.TACS, tacs.constitutive, tacs.elements, tacs.functions, tacs.problems.static

class TransientProblem(BaseProblem):
    """
    The main purpose of this class is to represent all relevant
    information for a static analysis. This will include
    information defining the loading condition as well as various
    other pieces of information.

    Parameters
    ----------
    name : str
        Name of this tacs problem

    Examples
    --------
    >>> sp = TransientProblem('lc0')
    """

    def __init__(self, name, tInit, tFinal, numSteps,
                 assembler, comm, outputViewer=None, meshLoader=None,
                 options={}):
        # python object name
        self.objectName = 'TransientProblem'

        # Problem name
        self.name = name

        # Defualt setup for common problem class objects
        super().__init__(assembler, comm, outputViewer, meshLoader)

        # Set time interval parmeters
        self.tInit = tInit
        self.tFinal = tFinal
        self.numSteps = numSteps

        # Default Option List
        defOpts = {
            'outputdir': [str, './'],

            # Solution Options
            'timeIntegrator': [str, 'BDF'],
            'integrationOrder': [int, 2],
            'L2Convergence': [float, 1e-12],
            'L2ConvergenceRel': [float, 1e-12],
            'jacAssemblyFreq': [int, 1],
            'useMonitor': [bool, False],
            'monitorFrequency': [int, 10],

            # Output Options
            'outputFrequency': [int, 0],
            'writeSolution': [bool, True],
            'numberSolutions': [bool, True],
            'printTiming': [bool, False],
            'printIterations': [bool, True],

        }

        # Process the default options which are added to self.options
        # under the 'defaults' key. Make sure the key are lower case
        self.options = {}
        def_keys = defOpts.keys()
        self.options['defaults'] = {}
        for key in def_keys:
            self.options['defaults'][key.lower()] = defOpts[key]
            self.options[key.lower()] = defOpts[key]

        # Set user-defined options
        for key in options:
            self.setOption(key, options[key])

        # Create problem-specific variables
        self._createVariables()

        # Set output viewer for integrator
        self.integrator.setFH5(self.outputViewer)
        outputFreq = self.getOption('outputFrequency')
        self.integrator.setOutputFrequency(outputFreq)
        outputDir = self.getOption('outputDir')
        self.integrator.setOutputPrefix(outputDir)

    def _createVariables(self):
        """Internal to create the objects required by TACS Integrator"""

        self.callCounter = -1

        # Create a force vector for each time step
        self.F = [self.assembler.createVec() for i in range(self.numSteps + 1)]
        # Auxillary element object for applying tractions/pressure
        self.auxElems = [tacs.TACS.AuxElements() for i in range(self.numSteps + 1)]

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

        if self.getOption('printIterations'):
            self.integrator.setPrintLevel(1)

        # Set solver tolerances
        atol = self.getOption('L2Convergence')
        self.integrator.setAbsTol(atol)
        rtol = self.getOption('L2ConvergenceRel')
        self.integrator.setRelTol(rtol)
        # Jacobian assembly frequency
        jacFreq = self.getOption('jacAssemblyFreq')
        self.integrator.setJacAssemblyFreq(jacFreq)

    def getNumTimeSteps(self):
        return self.numSteps

    def getTimeSteps(self):
        timeSteps = np.linspace(self.tInit, self.tFinal, self.numSteps + 1)
        return timeSteps

####### Static load methods ########

    def addLoadToComponents(self, timeStep, compIDs, F, averageLoad=False):
        """"
        The function is used to add a *FIXED TOTAL LOAD* on one or more
        components, defined by COMPIDs. The purpose of this routine is
        to add loads that remain fixed throughout an optimization. An example
        would be an engine load. This routine determines all the unqiue nodes
        in the FE model that are part of the the requested components, then
        takes the total 'force' by F and divides by the number of nodes.
        This average load is then applied to the nodes.

        NOTE: The units of the entries of the 'force' vector F are not
        necesarily physical forces and their interpretation depends
        on the physics problem being solved and the dofs included
        in the model.

        A couple of examples of force vector components for common problem are listed below:

        In Elasticity with varsPerNode = 3,
        F = [fx, fy, fz] # forces
        In Elasticity with varsPerNode = 6,
        F = [fx, fy, fz, mx, my, mz] # forces + moments
        In Thermoelasticity with varsPerNode = 4,
        F = [fx, fy, fz, Q] # forces + heat
        In Thermoelasticity with varsPerNode = 7,
        F = [fx, fy, fz, mx, my, mz, Q] # forces + moments + heat

        Parameters
        ----------

        compIDs : The components with added loads. Use selectCompIDs()
            to determine this.
        F : Numpy array length varsPerNode
            Vector of 'force' components
        """
        self._addLoadToComponents(self.F[timeStep], compIDs, F, averageLoad)

    def addLoadToNodes(self, timeStep, nodeIDs, F, nastranOrdering=False):
        """
        The function is used to add a fixed point load of F to the
        selected node IDs. This is similar to the addLoadToPoints method,
        except we select the load points based on node ID rather than
        physical location.

        NOTE: This should be the prefered method (over addLoadToPoints) for adding forces to
        specific nodes for the following reasons:
            1. This method is more efficient, as it does not require a
            closest point search to locate the node.
            2. In the case where the mesh features coincident nodes
            it is impossible to uniquely specify which node gets the load
            through x,y,z location, however the points can be specified uniquely by node ID.

        A couple of examples of force vector components for common problem are listed below:

        In Elasticity with varsPerNode = 3,
        F = [fx, fy, fz] # forces
        In Elasticity with varsPerNode = 6,
        F = [fx, fy, fz, mx, my, mz] # forces + moments
        In Thermoelasticity with varsPerNode = 4,
        F = [fx, fy, fz, Q] # forces + heat
        In Thermoelasticity with varsPerNode = 7,
        F = [fx, fy, fz, mx, my, mz, Q] # forces + moments + heat

        Parameters
        ----------

        nodeIDs : list[int]
            The nodes with added loads.
        F : Numpy 1d or 2d array length (varsPerNodes) or (numNodeIDs, varsPerNodes)
            Array of force vectors, one for each node. If only one force vector is provided,
            force will be copied uniformly across all nodes.
        nastranOrdering : bool
            Flag signaling whether nodeIDs are in TACS (default)
            or NASTRAN (grid IDs in bdf file) ordering
        """

        self._addLoadToNodes(self.F[timeStep], nodeIDs, F, nastranOrdering)

    def addTractionToComponents(self, timeStep, compIDs, tractions,
                                faceIndex=0):
        """
        The function is used to add a *FIXED TOTAL TRACTION* on one or more
        components, defined by COMPIDs. The purpose of this routine is
        to add loads that remain fixed throughout an optimization.

        Parameters
        ----------

        compIDs : The components with added loads. Use selectCompIDs()
            to determine this.
        tractions : Numpy array length 1 or compIDs
            Array of traction vectors for each components
        faceIndex : int
            Indicates which face (side) of element to apply traction to.
            Note: not required for certain elements (i.e. shells)
        """
        self._addTractionToComponents(self.auxElems[timeStep], compIDs, tractions, faceIndex)

    def addTractionToElements(self, timeStep, elemIDs, tractions,
                              faceIndex=0, nastranOrdering=False):
        """
        The function is used to add a fixed traction to the
        selected element IDs. Tractions can be specified on an
        element by element basis (if tractions is a 2d array) or
        set to a uniform value (if tractions is a 1d array)

        Parameters
        ----------

        elemIDs : List
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

        self._addTractionToElements(self.auxElems[timeStep], elemIDs, tractions, faceIndex, nastranOrdering)

    def addPressureToComponents(self, timeStep, compIDs, pressures,
                                faceIndex=0):
        """
        The function is used to add a *FIXED TOTAL PRESSURE* on one or more
        components, defined by COMPIds. The purpose of this routine is
        to add loads that remain fixed throughout an optimization. An example
        would be a fuel load.

        Parameters
        ----------

        compIDs : The components with added loads. Use selectCompIDs()
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
        The function is used to add a fixed presure to the
        selected element IDs. Pressures can be specified on an
        element by element basis (if pressures is an array) or
        set to a uniform value (if pressures is a scalar)

        Parameters
        ----------

        elemIDs : List
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

    ####### Static solver methods ########

    def _setProblemVars(self):
        """
        Make sure that the assembler is using
        the input variables associated with this problem
        """

        self.assembler.setDesignVars(self.x)
        self.assembler.setNodes(self.Xpts)

    def solve(self):
        """
        Solution of the structural system for loadCase. The
        forces must already be set.

        Parameters
        ----------
        Optional Arguments:

        Fext : ndarray
            Distributed array containing additional loads (i.e. aerodynamic)
            to applied to RHS of static problem.

        """
        startTime = time.time()

        self.callCounter += 1

        setupProblemTime = time.time()

        # Set problem vars to assembler
        self._setProblemVars()

        initSolveTime = time.time()

        # Loop over every time step and solve transient problem
        for i in range(self.numSteps + 1):
            # Set the auxilliary elements for this time step (tractions/pressures)
            self.assembler.setAuxElements(self.auxElems[i])
            self.integrator.iterate(i, forces=self.F[i])

        solveTime = time.time()

        # If timing was was requested print it, if the solution is nonlinear
        # print this information automatically if prinititerations was requested.
        if self.getOption('printTiming') or self.getOption('printIterations'):
            self.pp('+--------------------------------------------------+')
            self.pp('|')
            self.pp('| TACS Solve Times:')
            self.pp('|')
            self.pp('| %-30s: %10.3f sec' % ('TACS Setup Time', setupProblemTime - startTime))
            self.pp('| %-30s: %10.3f sec' % ('TACS Solve Init Time', initSolveTime - setupProblemTime))
            self.pp('| %-30s: %10.3f sec' % ('TACS Solve Time', solveTime - initSolveTime))
            self.pp('|')
            self.pp('| %-30s: %10.3f sec' % ('TACS Total Solution Time', solveTime - startTime))
            self.pp('+--------------------------------------------------+')

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
        >>> problem.solve()
        >>> problem.evalFunctions(funcs, ['mass'])
        >>> funcs
        >>> # Result will look like (if structProblem, sp, has name of 'c1'):
        >>> # {'cl_mass':12354.10}
        """
        startTime = time.time()

        # Set problem vars to assembler
        self._setProblemVars()

        if evalFuncs is None:
            evalFuncs = sorted(list(self.functionList))
        else:
            evalFuncs = sorted(list(evalFuncs))

        if not ignoreMissing:
            for f in evalFuncs:
                if not f in self.functionList:
                    raise Error("Supplied function '%s' has not been added "
                                "using addFunction()." % f)

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
            self.pp('+--------------------------------------------------+')
            self.pp('|')
            self.pp('| TACS Function Times:')
            self.pp('|')
            self.pp('| %-30s: %10.3f sec' % ('TACS Function Setup Time', setupProblemTime - startTime))
            self.pp('| %-30s: %10.3f sec' % ('TACS Function Eval Time', functionEvalTime - setupProblemTime))
            self.pp('| %-30s: %10.3f sec' % ('TACS Dict Time', dictAssignTime - functionEvalTime))
            self.pp('|')
            self.pp('| %-30s: %10.3f sec' % ('TACS Function Time', dictAssignTime - startTime))
            self.pp('+--------------------------------------------------+')

    def evalFunctionsSens(self, funcsSens, evalFuncs=None):
        """
        This is the main routine for returning useful (sensitivity)
        information from pytacs. The derivatives of the functions
        corresponding to the strings in EVAL_FUNCS are evaluated and
        updated into the provided dictionary.

        Parameters
        ----------
        structProblem : pyStructProblem class
            Structural problem to get the solution for
        funcsSens : dict
            Dictionary into which the derivatives are saved.
        evalFuncs : iterable object containing strings
            The functions the user wants returned

        Examples
        --------
        >>> funcsSens = {}
        >>> FEAsolver.evalFunctionsSens(sp, funcsSens, ['mass'])
        >>> funcs
        >>> # Result will look like (if structProblem, sp, has name of 'c1'):
        >>> # {'c1_mass':{'struct':[1.234, ..., 7.89]}
        """

        startTime = time.time()

        # Set problem vars to assembler
        self._setProblemVars()

        if evalFuncs is None:
            evalFuncs = sorted(list(self.functionList))
        else:
            evalFuncs = sorted(list(evalFuncs))

        for f in evalFuncs:
            if f not in self.functionList:
                raise self.TACSError("Supplied function has not been added "
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

        # Recast sensititivities into dict for user
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
            self.pp('+--------------------------------------------------+')
            self.pp('|')
            self.pp('| TACS Adjoint Times:')
            print('|')
            print('| %-30s: %10.3f sec' % ('Adjoint solve time', adjointFinishedTime - startTime))
            print('|')
            print('| %-30s: %10.3f sec' % ('Complete Sensitivity Time', totalSensitivityTime - startTime))
            print('+--------------------------------------------------+')

    ####### Post processing methods ########

    def getVariables(self, states=None):
        """Return the current state values for the current
        structProblem"""
        self.setStructProblem(structProblem)

        if states is None:
            states = self.u.getArray().copy()
        else:
            states[:] = self.u.getArray()

        return states

    def setVariables(self, states):
        """ Set the structural states for current load case. Typically
        only used for aerostructural analysis

        Parameters
        ----------
        states : array
            Values to set. Must be the size of getNumVariables()
        """
        self.setStructProblem(structProblem)
        self.u.setValues(states)
        self.assembler.setVariables(self.u)

    def writeSolution(self, timeSteps=None):
        """This is a generic shell function that writes the output
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
        """
        # Unless the writeSolution option is off write actual file:
        if self.getOption('writeSolution'):
            # Check input
            outputDir = self.getOption('outputDir')
            baseName = self.name
            if self.getOption('numberSolutions'):
                baseName = baseName + '_%3.3d' % self.callCounter

            base = os.path.join(outputDir, baseName) + '.f5'
            # Write specific time step out
            if timeSteps is not None:
                timeSteps = np.atleast_1d(timeSteps)
                for timeStep in timeSteps:
                    self.integrator.writeStepToF5(timeStep)
            # Write all time steps out
            else:
                self.integrator.writeSolutionToF5()