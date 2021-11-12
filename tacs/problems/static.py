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

class StaticProblem(BaseProblem):
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
    >>> sp = StaticProblem('lc0')
    """

    def __init__(self, name, assembler, comm, outputViewer=None, meshLoader=None, options={}):
        # python object name
        self.objectName = 'StaticProblem'

        # Problem name
        self.name = name

        # Defualt setup for common problem class objects
        super().__init__(assembler, comm, outputViewer, meshLoader)

        # Default Option List
        defOpts = {
            'outputdir': [str, './'],

            # Solution Options
            'KSMSolver': [str, 'GMRES'],
            'orderingType': [str, 'ND'],
            'PCFillLevel': [int, 1000],
            'PCFillRatio': [float, 20.0],
            'subSpaceSize': [int, 10],
            'nRestarts': [int, 15],
            'flexible': [int, 1],
            'L2Convergence': [float, 1e-12],
            'L2ConvergenceRel': [float, 1e-12],
            'useMonitor': [bool, False],
            'monitorFrequency': [int, 10],
            'resNormUB': [float, 1e20],

            # Output Options
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

        # Linear solver factor flag
        self._factorOnNext = True

        # Create problem-specific variables
        self._createVariables()

    def _createVariables(self):
        """Internal to create the variable required by TACS"""

        # Generic residual vector
        self.res = self.assembler.createVec()

        # Dictionaries to hold adjoint/sens vectors for each evalFunc
        self.adjointList = OrderedDict()
        self.dIduList = OrderedDict()
        self.dvSensList = OrderedDict()
        self.xptSensList = OrderedDict()

        # Load vector
        self.F = self.assembler.createVec()
        self.F_array = self.F.getArray()
        # State variable vector
        self.u = self.assembler.createVec()
        self.u_array = self.u.getArray()
        # Auxillary element object for applying tractions/pressure
        self.auxElems  = tacs.TACS.AuxElements()
        self.callCounter = -1

        # Norms
        self.initNorm = 0.0
        self.startNorm = 0.0
        self.finalNorm = 0.0

        opt = self.getOption

        # Tangent Stiffness --- process the ordering option here:
        tmp = opt('orderingType').lower()
        if tmp == 'natural':
            ordering = tacs.TACS.NATURAL_ORDER
        elif tmp == 'nd':
            ordering = tacs.TACS.ND_ORDER
        elif tmp == 'rcm':
            ordering = tacs.TACS.RCM_ORDER
        elif tmp == 'tacs_amd':
            ordering = tacs.TACS.TACS_AMD_ORDER
        elif tmp == 'multicolor':
            ordering = tacs.TACS.MULTICOLOR_ORDER
        else:
            raise self.TACSError("Unrecognized 'orderingType' option value: "
                        "'%s'. Valid values are: 'natural', 'nd', 'rcm', "
                        "'tacs_amd', or 'multicolor'." % tmp)

        self.K = self.assembler.createSchurMat(ordering)

        # Additional Vecs for updates
        self.update = self.assembler.createVec()

        # Setup PCScMat and KSM solver
        self.alpha = 1.0
        self.beta = 0.0
        self.gamma = 0.0
        self.assembler.assembleJacobian(self.alpha, self.beta, self.gamma, self.res, self.K)

        reorderSchur = 1
        self.PC = tacs.TACS.Pc(self.K, lev_fill=opt('PCFillLevel'),
                               ratio_fill=opt('PCFillRatio'), reorder=reorderSchur)

        # Operator, fill level, fill ratio, msub, rtol, ataol
        if opt('KSMSolver').upper() == 'GMRES':
            self.KSM = tacs.TACS.KSM(
                self.K, self.PC, opt('subSpaceSize'), opt('nRestarts'),
                opt('flexible'))
        # TODO: Fix this
        # elif opt('KSMSolver').upper() == 'GCROT':
        #    self.KSM = tacs.TACS.GCROT(
        #        self.K, self.PC, opt('subSpaceSize'), opt('subSpaceSize'),
        #        opt('nRestarts'), opt('flexible'))
        else:
            raise self.TACSError("Unknown KSMSolver option. Valid options are "
                        "'GMRES' or 'GCROT'")

        self.KSM.setTolerances(self.getOption('L2ConvergenceRel'),
                               self.getOption('L2Convergence'))

        if opt('useMonitor'):
            self.KSM.setMonitor(tacs.TACS.KSMPrintStdout(
                opt('KSMSolver'), self.comm.rank, opt('monitorFrequency')))

    def addFunction(self, funcName, funcHandle, compIDs=None, **kwargs):
        """
            Generic function to add a function for TACS. It is intended to
            be reasonably generic since the user supplies the actual
            function handle to use. The following functions can be used:
            KSFailure, KSBuckling, MaxBuckling, AverageKSFailure,
            MaxFailure, AverageMaxFailure, AverageKSBuckling,
            StructuralMass, Compliance, AggregateDisplacement.

            Parameters
            ----------
            funcName : str
                The user-supplied name for the function. This will
                typically be a string that is meanful to the user
            funcHandle : tacs.functions
                The fucntion handle to use for creation. This must come
                from the functions module in tacs.
            compIDs: list
                List of compIDs to select. Alternative to selectCompIDs
                arguments.
        """
        success = super().addFunction(funcName, funcHandle, compIDs, **kwargs)
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
        x : ndarray
            The variables (typically from the optimizer) to set. It
            looks for variable in the ``self.varName`` attribute.

        """
        super().setDesignVars(x)
        self._factorOnNext = True

    def setCoordinates(self, coords):
        """
        Set the mesh coordinates of the structure.

        Returns
        -------
        coords : array
            Structural coordinate in array of size (N, 3) where N is
            the number of structural nodes on this processor.
        """
        super().setCoordinates(coords)
        self._factorOnNext = True

####### Static load methods ########

    def addLoadToComponents(self, compIDs, F, averageLoad=False):
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
        self._addLoadToComponents(self.F, compIDs, F, averageLoad)

    def addLoadToNodes(self, nodeIDs, F, nastranOrdering=False):
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

        self._addLoadToNodes(self.F, nodeIDs, F, nastranOrdering)

    def addTractionToComponents(self, compIDs, tractions,
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
        self._addTractionToComponents(self.auxElems, compIDs, tractions, faceIndex)

    def addTractionToElements(self, elemIDs, tractions,
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

        self._addTractionToElements(self.auxElems, elemIDs, tractions, faceIndex, nastranOrdering)

    def addPressureToComponents(self, compIDs, pressures,
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
        self._addPressureToComponents(self.auxElems, compIDs, pressures, faceIndex)

    def addPressureToElements(self, elemIDs, pressures,
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

        self._addPressureToElements(self.auxElems, elemIDs, pressures,
                                    faceIndex, nastranOrdering)

    ####### Static solver methods ########

    def _setProblemVars(self):
        """
        Make sure that the assembler is using
        the input variables associated with this problem
        """

        self.assembler.setDesignVars(self.x)
        self.assembler.setNodes(self.Xpts)
        self.assembler.setAuxElements(self.auxElems)
        # Set state variables
        self.assembler.setVariables(self.u)
        # Zero any time derivitive terms
        self.assembler.zeroDotVariables()
        self.assembler.zeroDDotVariables()

    def _initializeSolve(self):
        """
        Initialze the solution of the structural system for the
        loadCase. The stiffness matrix is assembled and factored.
        """

        if self._factorOnNext:
            self.assembler.assembleJacobian(self.alpha, self.beta, self.gamma, self.res, self.K)
            self.PC.factor()
            self._factorOnNext = False

    def solve(self, Fext=None):
        """
        Solution of the structural system for loadCase. The
        forces must already be set.

        Parameters
        ----------
        Optional Arguments:

        Fext : numpy array
            Distributed array containing additional loads (ex. aerodynamic forces for aerostructural coupling)
            to applied to RHS of the static problem.

        """
        startTime = time.time()

        self.callCounter += 1

        setupProblemTime = time.time()

        # Set problem vars to assembler
        self._setProblemVars()

        # Check if we need to initialize
        self._initializeSolve()

        initSolveTime = time.time()

        # Compute the RHS
        self.assembler.assembleRes(self.res)
        # Zero out bc terms in F
        self.assembler.applyBCs(self.F)
        # Add the -F
        self.res.axpy(-1.0, self.F)
        # Add external forces, if necessary
        if Fext is not None:
            resArray = self.res.getArray()
            resArray[:] -= Fext[:]

        # Set initnorm as the norm of F
        self.initNorm = np.real(self.F.norm())

        # Starting Norm for this compuation
        self.startNorm = np.real(self.res.norm())

        initNormTime = time.time()

        # Solve Linear System for the update
        self.KSM.solve(self.res, self.update)

        self.update.scale(-1.0)

        solveTime = time.time()

        # Update State Variables
        self.assembler.getVariables(self.u)
        self.u.axpy(1.0, self.update)
        self.assembler.setVariables(self.u)

        stateUpdateTime = time.time()

        # Compute final FEA Norm
        self.assembler.assembleRes(self.res)
        self.res.axpy(-1.0, self.F)  # Add the -F
        self.finalNorm = np.real(self.res.norm())

        finalNormTime = time.time()

        # If timing was was requested print it, if the solution is nonlinear
        # print this information automatically if prinititerations was requested.
        if self.getOption('printTiming') or self.getOption('printIterations'):
            self.pp('+--------------------------------------------------+')
            self.pp('|')
            self.pp('| TACS Solve Times:')
            self.pp('|')
            self.pp('| %-30s: %10.3f sec' % ('TACS Setup Time', setupProblemTime - startTime))
            self.pp('| %-30s: %10.3f sec' % ('TACS Solve Init Time', initSolveTime - setupProblemTime))
            self.pp('| %-30s: %10.3f sec' % ('TACS Init Norm Time', initNormTime - initSolveTime))
            self.pp('| %-30s: %10.3f sec' % ('TACS Solve Time', solveTime - initNormTime))
            self.pp('| %-30s: %10.3f sec' % ('TACS State Update Time', stateUpdateTime - solveTime))
            self.pp('| %-30s: %10.3f sec' % ('TACS Final Norm Time', finalNormTime - stateUpdateTime))
            self.pp('|')
            self.pp('| %-30s: %10.3f sec' % ('TACS Total Solution Time', finalNormTime - startTime))
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
        >>> FEAsolver(sp)
        >>> FEAsolver.evalFunctions(funcs, ['mass'])
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
        funcVals = self.assembler.evalFunctions(handles)

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
        # Check that the functions are all ok.
        # and prepare tacs vecs for adjoint procedure
        dvSenses = []
        xptSenses = []
        dIdus = []
        adjoints = []
        for f in evalFuncs:
            if f not in self.functionList:
                raise Error("Supplied function has not beed added "
                            "using addFunction()")
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
            key = self.name + '_%s' % f
            # Finalize sensitivity arrays across all procs
            dvSenses[i].beginSetValues()
            dvSenses[i].endSetValues()
            xptSenses[i].beginSetValues()
            xptSenses[i].endSetValues()
            # Return sensitivities as array in sens dict
            funcsSens[key] = {self.varName: dvSenses[i].getArray().copy(),
                              self.coordName: xptSenses[i].getArray().copy()}

        totalSensitivityTime = time.time()

        if self.getOption('printTiming'):
            self.pp('+--------------------------------------------------+')
            self.pp('|')
            self.pp('| TACS Adjoint Times:')
            print('|')
            print('| %-30s: %10.3f sec' % ('TACS Sens Setup Problem Time', setupProblemTime - startTime))
            print('| %-30s: %10.3f sec' % (
                'TACS Adjoint RHS Time', adjointRHSTime - setupProblemTime))
            for f in evalFuncs:
                print('| %-30s: %10.3f sec' % (
                'TACS Adjoint Solve Time - %s' % (f), adjointEndTime[f] - adjointStartTime[f]))
            print('| %-30s: %10.3f sec' % ('Total Sensitivity Time', totalSensitivityTime - adjointFinishedTime))
            print('|')
            print('| %-30s: %10.3f sec' % ('Complete Sensitivity Time', totalSensitivityTime - startTime))
            print('+--------------------------------------------------+')

    def addSVSens(self, evalFuncs, dIduList):
        """ Add the state variable sensitivity to the ADjoint RHS for given evalFuncs"""
        funcHandles = [self.functionList[f] for f in evalFuncs if
                       f in self.functionList]
        self.assembler.addSVSens(funcHandles, dIduList, self.alpha, self.beta, self.gamma)

    def addDVSens(self, evalFuncs, dvSensList, scale=1.0):
        """ Add pratial sensitivity contribution due to design vars for evalFuncs"""
        funcHandles = [self.functionList[f] for f in evalFuncs if
                       f in self.functionList]
        self.assembler.addDVSens(funcHandles, dvSensList, scale)

    def addAdjointResProducts(self, adjointlist, dvSensList, scale=-1.0):
        """ Add the adjoint product contribution to the design variable sensitivity arrays"""
        self.assembler.addAdjointResProducts(adjointlist, dvSensList, scale)

    def addXptSens(self, evalFuncs, xptSensList, scale=1.0):
        """ Add pratial sensitivity contribution due to nodal coordinates for evalFuncs"""
        funcHandles = [self.functionList[f] for f in evalFuncs if
                       f in self.functionList]
        self.assembler.addXptSens(funcHandles, xptSensList, scale)

    def addAdjointResXptSensProducts(self, adjointlist, xptSensList, scale=-1.0):
        """ Add the adjoint product contribution to the nodal coordinates sensitivity arrays"""
        self.assembler.addAdjointResXptSensProducts(adjointlist, xptSensList, scale)

    def getResidual(self, res=None, Fext=None):
        """
        This routine is used to evaluate directly the structural
        residual. Only typically used with aerostructural analysis.

        Parameters
        ----------
        res : numpy array
            If res is not None, place the residuals into this array.

        Fext : numpy array
            Distributed array containing additional loads (ex. aerodynamic forces for aerostructural coupling)
            to applied to RHS of the static problem.

        Returns
        -------
        res : array
            The same array if res was provided, (otherwise a new
            array) with evaluated residuals
        """
        # Make sure assembler variables are up to date
        self._setProblemVars()
        # Assemble residual
        self.assembler.assembleRes(self.res)
        # Add the -F
        self.res.axpy(1.0, self.F)

        # Add external loads, if specified
        if Fext is not None:
            resArray = self.res.getArray()
            resArray[:] -= Fext[:]

        # Output residual
        if res is None:
            res = self.res.getArray().copy()
        else:
            res[:] = self.res.getArray()

        return res

    def zeroVectors(self):
        """Zero all the tacs solution b-vecs"""
        self.res.zeroEntries()
        self.u.zeroEntries()
        self.assembler.setVariables(self.u)
        self.update.zeroEntries()

    def solveAdjoint(self, rhs, phi, damp=1.0):
        """
        Solve the structural adjoint.

        Parameters
        ----------
        rhs : TACS BVec
            right hand side vector for adjoint solve
        phi : TACS BVec
            BVec into which the adjoint is saved
        damp : float
            A damping variable for adjoint update. Typically only used
            in multidisciplinary analysis
        """

        # First compute the residual
        self.K.mult(phi, self.res)
        self.res.axpy(-1.0, rhs)  # Add the -RHS

        # Starting Norm for this compuation
        self.startNorm = np.real(self.res.norm())

        # Solve Linear System
        zeroGuess = 0
        self.update.zeroEntries()
        self.KSM.solve(self.res, self.update, zeroGuess)

        # Update the adjoint vector with the (damped) update
        phi.axpy(-damp, self.update)

        # Compute actual final FEA Norm
        self.K.mult(phi, self.res)
        self.res.axpy(-1.0, rhs)  # Add the -RHS
        self.finalNorm = np.real(self.res.norm())

    def getVariables(self, states=None):
        """Return the current state values for the
        problem"""

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
        self.u_array[:] = states[:]
        self.assembler.setVariables(self.u)

    def writeSolution(self, outputDir=None, baseName=None, number=None):
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
        # Make sure assembler variables are up to date
        self._setProblemVars()

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
            base = os.path.join(outputDir, baseName) + '.f5'
            self.outputViewer.writeToFile(base)