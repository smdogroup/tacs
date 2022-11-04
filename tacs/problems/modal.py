"""
The main purpose of this class is to represent all relevant
information for a modal analysis.

.. note:: This class should be created using the
    :meth:`pyTACS.createModalProblem <tacs.pytacs.pyTACS.createModalProblem>` method.
"""

# =============================================================================
# Imports
# =============================================================================
import os
import numpy as np
import time
from .base import TACSProblem
import tacs.TACS


class ModalProblem(TACSProblem):

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
        options={},
    ):
        """
        NOTE: This class should not be initialized directly by the user.
        Use pyTACS.createModalProblem instead.

        Parameters
        ----------
        name : str
            Name of this tacs problem

        sigma : float
            Guess for the lowest eigenvalue. This corresponds to the lowest frequency squared. (rad^2/s^2)

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

        # Default setup for common problem class objects
        TACSProblem.__init__(self, assembler, comm, outputViewer, meshLoader)

        # Set time eigenvalue parameters
        self.sigma = sigma
        self.numEigs = numEigs

        # String name used in evalFunctions
        self.valName = "eigsm"
        self._initializeFunctionList()

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

        # Create problem-specific variables
        self._createVariables()

    def _createVariables(self):
        """
        Internal to create the objects required by TACS Integrator
        """

        self.callCounter = -1

        # Solve the eigenvalue problem
        self.M = self.assembler.createSchurMat()
        self.K = self.assembler.createSchurMat()

        self.pc = tacs.TACS.Pc(self.K)

        # Assemble and factor the stiffness/Jacobian matrix. Factor the
        # Jacobian and solve the linear system for the displacements
        alpha = 1.0
        beta = 0.0
        gamma = 0.0
        self.assembler.assembleJacobian(alpha, beta, gamma, None, self.K)
        self.pc.factor()  # LU factorization of stiffness matrix

        subspace = self.getOption("subSpaceSize")
        restarts = self.getOption("nRestarts")
        self.gmres = tacs.TACS.KSM(self.K, self.pc, subspace, restarts)

        eigTol = self.getOption("L2Convergence")

        # Create the frequency analysis object
        self.freqSolver = tacs.TACS.FrequencyAnalysis(
            self.assembler,
            self.sigma,
            self.M,
            self.K,
            self.gmres,
            num_eigs=self.numEigs,
            eig_tol=eigTol,
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
        the integers in EVAL_FUNCS are evaluated and updated into
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
        >>> modalProblem.solve()
        >>> modalProblem.evalFunctions(funcs, 'eigsm.0')
        >>> funcs
        >>> # Result will look like (if ModalProblem has name of 'c1'):
        >>> # {'c1_eigsm.0':12354.10}
        """
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
        >>> modalProblem.evalFunctionsSens(funcsSens, 'eigsm.0')
        >>> funcsSens
        >>> # Result will look like (if ModalProblem has name of 'c1'):
        >>> # {'c1_eigsm.0':{'struct':[1.234, ..., 7.89], 'Xpts':[3.14, ..., 1.59]}}
        """
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

        # Loop through each requested eigenvalue
        for funcName in evalFuncs:
            mode_i = evalFuncs[funcName]
            key = f"{self.name}_{funcName}"
            funcsSens[key] = {}
            # Evaluate dv sens
            self.freqSolver.evalEigenDVSens(mode_i, dvSens)
            funcsSens[key][self.varName] = dvSens.getArray().copy()
            # Evaluate nodal sens
            self.freqSolver.evalEigenXptSens(mode_i, xptSens)
            funcsSens[key][self.coordName] = xptSens.getArray().copy()

    ####### Modal solver methods ########

    def _updateAssemblerVars(self):
        """
        Make sure that the assembler is using
        the input variables associated with this problem
        """

        self.assembler.setDesignVars(self.x)
        self.assembler.setNodes(self.Xpts)

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

        # Solve the frequency analysis problem
        self.freqSolver.solve(
            print_flag=self.getOption("printLevel"),
            print_level=self.getOption("printLevel"),
        )

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

        states : TACS.Vec or numpy.ndarray or None
            Place eigenvector for mode into this array (optional).

        Returns
        --------
        eigVal: float
            Eigenvalue for mode corresponds to square of eigenfrequency (rad^2/s^2)

        states : numpy.ndarray
            Eigenvector for mode
        """
        eigVal, err = self.freqSolver.extractEigenvalue(index)
        eigVector = self.assembler.createVec()
        self.freqSolver.extractEigenvector(index, eigVector)
        # Inplace assignment if vectors were provided
        if isinstance(states, tacs.TACS.Vec):
            states.copyValues(eigVector)
        elif isinstance(states, np.ndarray):
            states[:] = eigVector.getArray()
        return eigVal, eigVector.getArray()

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
            Use the user spplied number to index solution. Again, only
            typically used from an external solver
        indices : int or list[int] or None
            Mode index or indices to get state variables for.
            If None, returns a solution for all modes.
            Defaults to None.
        """
        # Make sure assembler variables are up to date
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
