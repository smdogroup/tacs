"""
The main purpose of this class is to place linear constraints on design variables within the same component.
The constraints are of the form:

    c = a_0 * dv_0 + a_1 * dv_1 + ... + a_n * dv_n

Where which design variables to include (dv_0, dv_1, etc.) and the corresponding weights (a_0, a_1, etc.) are
defined by the user.

As an example, this constraint can be used to enforce ply fraction constraints for composite design optimization:

    pf_0 + pf_45 + pf_m45 + pf_90 = 1.0

Or to enforce that stiffener and panel thicknesses of a blade-stiffened panel do not differ by too much:

    -delta_t < st - pt < delta_t

For applying constraints on design variables across components, see :class:`AdjacencyConstraint <tacs.constraints.AdjacencyConstraint>`.

.. note:: This class should be created using the
    :meth:`pyTACS.createDVConstraint <tacs.pytacs.pyTACS.createDVConstraint>` method.
"""

# =============================================================================
# Imports
# =============================================================================
import scipy as sp

from tacs.constraints.base import TACSConstraint, SparseLinearConstraint


class DVConstraint(TACSConstraint):
    def __init__(
        self,
        name,
        assembler,
        comm,
        outputViewer=None,
        meshLoader=None,
        options=None,
    ):
        """
        NOTE: This class should not be initialized directly by the user.
        Use pyTACS.createDVConstraint instead.

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

        # Default setup for common constraint class objects, sets up comm and options
        TACSConstraint.__init__(
            self, assembler, comm, options, outputViewer, meshLoader
        )

    def addConstraint(
        self, conName, compIDs=None, lower=-1e20, upper=1e20, dvIndices=0, dvWeights=1.0
    ):
        """
        Generic method to adding a new constraint set for TACS.

        Parameters
        ----------
        conName : str
            The user-supplied name for the constraint set. This will
            typically be a string that is meaningful to the user

        compIDs: list[int] or None
            List of compIDs to select. One constraint will be added for each component.
            If None, all compIDs will be selected. Defaults to None.

        lower: float or complex
            lower bound for constraint. Defaults to -1e20.

        upper: float or complex
            upper bound for constraint. Defaults to 1e20.

        dvIndices : int or array-like[int]
            Index numbers of element DVs to be used in constraint.
            Defaults to 0.

        dvWeights : float or complex or array-like[float] or array-like[complex]
            Linear scaling factors for each DV used in constraint definition.
            If list, should match length of dvIndices. Defaults to 1's.

        """
        if compIDs is not None:
            # Make sure CompIDs is flat and get element numbers on each proc corresponding to specified compIDs
            compIDs = self._flatten(compIDs)
        else:
            nComps = self.meshLoader.getNumComponents()
            compIDs = list(range(nComps))

        if hasattr(dvIndices, "__iter__"):
            dvIndices = list(dvIndices)
        elif isinstance(dvIndices, int):
            dvIndices = [dvIndices]

        if hasattr(dvWeights, "__iter__"):
            dvWeights = list(dvWeights)
        elif isinstance(dvWeights, float) or isinstance(dvWeights, complex):
            dvWeights = [dvWeights]

        constrObj = self._createConstraint(dvIndices, dvWeights, compIDs, lower, upper)
        if constrObj.nCon > 0:
            self.constraintList[conName] = constrObj
            success = True
        else:
            self._TACSWarning(f"No valid `compIDs` provided. Skipping {conName}.")
            success = False

        return success

    def _createConstraint(self, dvIndices, dvWeights, compIDs, lbound, ubound):
        """
        Create a new constraint object for TACS.

        Parameters
        ----------
        dvIndices : list[int]
            Index numbers of element DVs to be used in constraint.
            Defaults to 0.

        dvWeights : list[float or complex]
            Linear scaling factors for each DV used in constraint definition.
            If list, should match length of dvIndices. Defaults to 1's.

        compIDs: list[int]
            List of compIDs to select.

        lbound: float or complex
            lower bound for constraint. Defaults to 0.0.

        ubound: float or complex
            upper bound for constraint. Defaults to 1e20.

        Returns
        -------
        constraint : tacs.constraints.base.SparseLinearConstraint or None
            Constraint object if successful, None otherwise.

        """
        # Assemble constraint info
        conCount = 0
        rows = []
        cols = []
        vals = []
        for comp in compIDs:
            # Get the TACS element object associated with this compID
            elemObj = self.meshLoader.getElementObject(comp, 0)
            elemIndex = 0
            # Get the dvs owned by this element
            globalDvNums = elemObj.getDesignVarNums(elemIndex)
            # Check if each specified dv num is owned by this proc
            for dvIndex, dvWeight in zip(dvIndices, dvWeights):
                if globalDvNums[dvIndex] in self.globalToLocalDVNums:
                    globalDVNum = globalDvNums[dvIndex]
                    localDVNum = self.globalToLocalDVNums[globalDVNum]
                    rows.append(conCount)
                    cols.append(localDVNum)
                    vals.append(dvWeight)
            conCount += 1

        nLocalDVs = self.getNumDesignVars()

        return SparseLinearConstraint(
            self.comm, rows, cols, vals, conCount, nLocalDVs, lbound, ubound
        )

    def evalConstraints(self, funcs, evalCons=None, ignoreMissing=False):
        """
        Evaluate values for constraints. The constraints corresponding to the strings in
        evalCons are evaluated and updated into the provided
        dictionary.

        Parameters
        ----------
        funcs : dict
            Dictionary into which the constraints are saved.
        evalCons : iterable object containing strings.
            If not none, use these constraints to evaluate.
        ignoreMissing : bool
            Flag to supress checking for a valid constraint. Please use
            this option with caution.

        Examples
        --------
        >>> funcs = {}
        >>> dvConstraint.evalConstraints(funcs, 'LE_SPAR')
        >>> funcs
        >>> # Result will look like (if DVConstraint has name of 'c1'):
        >>> # {'c1_LE_SPAR': array([12354.10])}
        """
        # Check if user specified which constraints to output
        # Otherwise, output them all
        evalCons = self._processEvalCons(evalCons, ignoreMissing)

        # Loop through each requested constraint set
        for conName in evalCons:
            key = f"{self.name}_{conName}"
            funcs[key] = self.constraintList[conName].evalCon(self.x.getArray())

    def evalConstraintsSens(self, funcsSens, evalCons=None):
        """
        This is the main routine for returning useful (sensitivity)
        information from constraint. The derivatives of the constraints
        corresponding to the strings in evalCons are evaluated and
        updated into the provided dictionary. The derivitives with
        respect to all design variables and node locations are computed.

        Parameters
        ----------
        funcsSens : dict
            Dictionary into which the derivatives are saved.
        evalCons : iterable object containing strings
            The constraints the user wants returned

        Examples
        --------
        >>> funcsSens = {}
        >>> dvConstraint.evalConstraintsSens(funcsSens, 'LE_SPAR')
        >>> funcsSens
        >>> # Result will look like (if DVConstraint has name of 'c1'):
        >>> # {'c1_LE_SPAR':{'struct':<50x242 sparse matrix of type '<class 'numpy.float64'>' with 100 stored elements in Compressed Sparse Row format>}}
        """
        # Check if user specified which constraints to output
        # Otherwise, output them all
        evalCons = self._processEvalCons(evalCons)

        # Get number of nodes coords on this proc
        nCoords = self.getNumCoordinates()

        # Loop through each requested constraint set
        for conName in evalCons:
            key = f"{self.name}_{conName}"
            # Get sparse Jacobian for dv sensitivity
            funcsSens[key] = {}
            funcsSens[key][self.varName] = self.constraintList[conName].evalConSens(
                self.x.getArray()
            )

            # Nodal sensitivities are always zero for this constraint,
            # Add an empty sparse matrix
            nCon = self.constraintList[conName].nCon
            funcsSens[key][self.coordName] = sp.sparse.csr_matrix(
                (nCon, nCoords), dtype=self.dtype
            )
