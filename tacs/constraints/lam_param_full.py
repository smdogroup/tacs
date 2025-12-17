"""
Lamination parameter constraints for composite design variables.

This module provides a TACS constraint implementation that enforces
non-linear relationships between lamination parameters for composite
components. For each selected component, three scalar constraints are
evaluated to enforce in-plane and bending feasibility limits.

Lamination parameter ordering
-----------------------------
Internally, lamination parameters are stored per-component in a fixed
6-entry ordering (indexed 1..6 for user-facing DV indices):

    1 -> V1
    2 -> V3
    3 -> W1
    4 -> W2
    5 -> W3
    6 -> W4

Constraint definitions (per component)
-------------------------------------
- Constraint 1 (in-plane relation):
    f1 = 2 * V1^2 - V3

- Constraint 2 (bending magnitude):
    f2 = W1^2 + W2^2

- Constraint 3 (bending feasibility):
    There are two cases depending on which W-parameters are present
    for a given component:

    * If W1 and W3 exist but W2 and W4 do not:
        f3 = 2 * W1^2 - W3

    * Otherwise (general case, missing W values treated as zero):
        f3 = 2*W1^2*(1 - W3) + 2*W2^2*(1 + W3) + W3^2 + W4^2 - 4*W1*W2*W4

Usage
-----
Create the constraint via :meth:`pyTACS.createLamParamFullConstraint`.

.. note:: This class should be created using the
    :meth:`pyTACS.createLamParamFullConstraint <tacs.pytacs.pyTACS.createLamParamFullConstraint>` method.
"""

# =============================================================================
# Imports
# =============================================================================
import scipy as sp
import numpy as np

from tacs.constraints.base import TACSConstraint


class LamParamFullConstraint(TACSConstraint):
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
        Use pyTACS.createLamParamFullConstraint instead.

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
        self, conName, compIDs=None, lower=-1e20, upper=1.0, dvIndices=0, dvWeights=1.0
    ):
        """
        Generic method to adding a new constraint set for TACS.

        Parameters
        ----------
        conName : str
            The user-supplied name for the constraint set. This will
            typically be a string that is meaningful to the user

        compIDs: list[int] or None
            List of compIDs to select. Three constraints will be added for each component.
            If None, all compIDs will be selected. Defaults to None.

        lower: float or complex
            Lower bound for constraints. Defaults to -1e20 for all lamination parameter constraints.

        upper: float or complex
            Upper bound for constraints. Defaults to 1.0 for all lamination parameter constraints.

        dvIndices : int or array-like[int]
            Index numbers of lamination parameter DVs to be used in constraint.
            Indices permitted correspond to the internal object DV numbers
            of the lamination parameters, that is [1,2,3,4,5,6] = [V1, V3, W1, W2, W3, W4]
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
            Index numbers of lamination parameter DVs to be used in constraint.
            Indices permitted correspond to the internal object DV numbers
            of the lamination parameters, that is [1,2,3,4,5,6] = [V1, V3, W1, W2, W3, W4]
            Defaults to 0.

        dvWeights : list[float or complex]
            Linear scaling factors for each DV used in constraint definition.
            If list, should match length of dvIndices. Defaults to 1's.

        compIDs: list[int]
            List of compIDs to select.

        lbound: float or complex
            Lower bound for constraints. Defaults to -1e20.

        ubound: float or complex
            Upper bound for constraints. Defaults to 1.

        Returns
        -------
        constraint : tacs.constraints.base.SparseLinearConstraint or None
            Constraint object if successful, None otherwise.

        """

        # We map lamination parameters into the fixed 6-entry ordering
        # indices: [1..6] -> [V1, V3, W1, W2, W3, W4]
        TOTAL_LP = 6
        nComps = len(compIDs)

        # Build weight lookup for provided dvIndices
        # dvIndices are element-local indices (e.g. 1..6). Map to weights
        weightByDv = {}
        for dv, w in zip(dvIndices, dvWeights):
            weightByDv[int(dv)] = w

        # Rows are indexed as (compIndex * TOTAL_LP + lpIndex)
        rows = []
        cols = []
        vals = []
        for compIndex, comp in enumerate(compIDs):
            # Get the TACS element object associated with this compID
            elemObj = self.meshLoader.getElementObject(comp, 0)
            elemIndex = 0
            # Get the dvs owned by this element
            globalDvNums = elemObj.getDesignVarNums(elemIndex)

            # For each of the 6 lamination parameter slots, check if the
            # corresponding element DV index (1..6) was requested and owned
            for lpIndex in range(TOTAL_LP):
                dvIndex = lpIndex + 1
                if dvIndex in weightByDv:
                    # Ensure the element actually has this DV slot
                    if (
                        dvIndex < len(globalDvNums)
                        and globalDvNums[dvIndex] in self.globalToLocalDVNums
                    ):
                        globalDVNum = globalDvNums[dvIndex]
                        localDVNum = self.globalToLocalDVNums[globalDVNum]
                        rows.append(compIndex * TOTAL_LP + lpIndex)
                        cols.append(localDVNum)
                        vals.append(weightByDv[dvIndex])

        nLocalDVs = self.getNumDesignVars()

        return SparseLamParamFullConstraint(
            self.comm,
            rows,
            cols,
            vals,
            nComps,
            TOTAL_LP,
            nLocalDVs,
            lbound,
            ubound,
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
            Flag to suppress checking for a valid constraint. Please use
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
        updated into the provided dictionary. The derivatives with
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


class SparseLamParamFullConstraint(object):
    """
    Sparse constraint object for lamination-parameter-based constraints.

    This object holds a sparse mapping matrix `A_lp` from local design
    variables to per-component lamination parameters, evaluates the
    three nonlinear constraints per component, and returns analytic
    sparse sensitivities.

    Attributes
    ----------
    A_lp : scipy.sparse.csr_matrix
        Sparse mapping of shape (nComps * nLP, nLocalDVs) mapping local
        DVs -> lamination parameters.
    nLP : int
        Number of lamination-parameter slots per component (6).
    nComps : int
        Number of components selected for the constraint.
    nCon : int
        Number of scalar constraints (3 * nComps).
    present : ndarray of bool
        Boolean mask of shape (nComps, nLP) indicating which LP rows
        are present (non-zero) for each component.
    """

    dtype = TACSConstraint.dtype

    def __init__(self, comm, rows, cols, vals, nComps, nLP, ncols, lb=-1e20, ub=1e20):
        """
        Initialize the sparse lamination-parameter constraint object.

        Parameters
        ----------
        comm : MPI.Intracomm
            MPI communicator used for allreduce operations.
        rows, cols, vals : sequence
            Row indices, column indices and values for constructing the
            CSR mapping matrix `A_lp` of shape (nComps*nLP, ncols).
        nComps : int
            Number of components selected.
        nLP : int
            Number of lamination-parameter slots per component.
        ncols : int
            Number of local design variables (columns of `A_lp`).
        lb, ub : float or array-like
            Lower and upper bounds for each scalar constraint. Can be a
            scalar (applies to all constraints), an array of length
            `nComps` (applies to each component), or an array of length
            `3*nComps` (one bound per scalar constraint).
        """
        # Sparse mapping from design vars to lamination parameters
        self.A_lp = sp.sparse.csr_matrix(
            (vals, (rows, cols)), shape=(nComps * nLP, ncols), dtype=self.dtype
        )
        # Number of lamination parameters per component
        self.nLP = nLP
        # Number of components
        self.nComps = nComps
        # Number of nonlinear constraints: 3 per component
        self.nCon = 3 * nComps
        # MPI comm
        self.comm = comm

        # Determine which lamination-parameter rows are present (nonzero)
        # row_nnz is length nComps * nLP
        rowNnz = np.diff(self.A_lp.indptr)
        # Get a boolean array of shape (nComps, nLP) indicating presence
        self.present = rowNnz.reshape(self.nComps, self.nLP) > 0

        # Save bound information (apply to each constraint)
        if isinstance(lb, np.ndarray) and len(lb) == self.nCon:
            self.lb = lb.astype(self.dtype)
        elif isinstance(lb, float) or isinstance(lb, complex):
            self.lb = np.array([lb] * self.nCon, dtype=self.dtype)
        else:
            # If an array of length nComps provided, expand to 3 per comp
            lbArr = np.asarray(lb)
            if lbArr.ndim == 1 and len(lbArr) == nComps:
                self.lb = np.repeat(lbArr.astype(self.dtype), 3)
            else:
                self.lb = np.array([lb] * self.nCon, dtype=self.dtype)

        if isinstance(ub, np.ndarray) and len(ub) == self.nCon:
            self.ub = ub.astype(self.dtype)
        elif isinstance(ub, float) or isinstance(ub, complex):
            self.ub = np.array([ub] * self.nCon, dtype=self.dtype)
        else:
            ubArr = np.asarray(ub)
            if ubArr.ndim == 1 and len(ubArr) == nComps:
                self.ub = np.repeat(ubArr.astype(self.dtype), 3)
            else:
                self.ub = np.array([ub] * self.nCon, dtype=self.dtype)

    def evalCon(self, x):
        """
        Evaluate the nonlinear constraints for the provided local design
        variable vector `x`.

        Parameters
        ----------
        x : ndarray
            Local design variable vector.

        Returns
        -------
        ndarray
            An array of length `3 * nComps` containing the three scalar
            constraint values per component.
        """
        # Compute lamination parameters globally. 
        # Compute the local contribution and then sum across MPI ranks.
        lpGlobal = self.comm.allreduce(self.A_lp.dot(x))
        # lpGlobal is length (nComps * nLP) so reshape to (nComps, nLP) for easier access
        lpMat = lpGlobal.reshape(self.nComps, self.nLP)

        conList = []
        for compIndex in range(self.nComps):
            lp = lpMat[compIndex]

            # Map indices (values default to 0.0 if not present)
            V1 = lp[0]
            V3 = lp[1]
            W1 = lp[2]
            W2 = lp[3]
            W3 = lp[4]
            W4 = lp[5]

            hasV1 = bool(self.present[compIndex, 0])
            hasV3 = bool(self.present[compIndex, 1])
            hasW1 = bool(self.present[compIndex, 2])
            hasW2 = bool(self.present[compIndex, 3])
            hasW3 = bool(self.present[compIndex, 4])
            hasW4 = bool(self.present[compIndex, 5])

            # Constraint 1: in-plane relation (only if V1 or V3 present)
            if hasV1 or hasV3:
                con1 = 2.0 * V1 * V1 - V3
            else:
                con1 = 0.0

            # Constraint 2: bending magnitude (requires W1 and W2)
            if hasW1 and hasW2:
                con2 = W1 * W1 + W2 * W2
            else:
                con2 = 0.0

            # Constraint 3: bending feasibility
            if hasW1 and hasW3 and not (hasW2 or hasW4):
                # Case where W1 and W3 exist but W2 and W4 do not
                con3 = 2.0 * W1 * W1 - W3
            elif hasW1 or hasW2 or hasW3 or hasW4:
                # General case (missing values are zero)
                con3 = (
                    2.0 * W1 * W1 * (1.0 - W3)
                    + 2.0 * W2 * W2 * (1.0 + W3)
                    + W3 * W3
                    + W4 * W4
                    - 4.0 * W1 * W2 * W4
                )
            else:
                con3 = 0.0

            conList.extend([con1, con2, con3])

        return np.asarray(conList, dtype=self.dtype)

    def evalConSens(self, x=None):
        """
        Compute analytic sparse sensitivities of the constraints w.r.t
        the local design variables `x`.

        Parameters
        ----------
        x : ndarray
            Local design variable vector.

        Returns
        -------
        scipy.sparse.csr_matrix
            Sparse Jacobian of shape `(nCon, nLocalDVs)`.
        """
        # Compute lamination parameters globally
        lpGlobal = self.comm.allreduce(self.A_lp.dot(x))
        lpMat = lpGlobal.reshape(self.nComps, self.nLP)

        # Build Jacobian of constraints w.r.t lamination parameters (sparse)
        # J_lp has shape (nCon, nComps*nLP)
        data = []
        rowIdx = []
        colIdx = []

        for compIndex in range(self.nComps):
            lp = lpMat[compIndex]
            baseRow = compIndex * 3
            baseLp = compIndex * self.nLP

            # Values
            V1 = lp[0]
            V3 = lp[1]
            W1 = lp[2]
            W2 = lp[3]
            W3 = lp[4]
            W4 = lp[5]

            hasV1 = bool(self.present[compIndex, 0])
            hasV3 = bool(self.present[compIndex, 1])
            hasW1 = bool(self.present[compIndex, 2])
            hasW2 = bool(self.present[compIndex, 3])
            hasW3 = bool(self.present[compIndex, 4])
            hasW4 = bool(self.present[compIndex, 5])

            # d(con1)/d(lp): con1 = 2 V1^2 - V3 --> dV1 = 4 V1, dV3 = -1
            if hasV1:
                rowIdx.append(baseRow)
                colIdx.append(baseLp + 0)
                data.append(4.0 * V1)
            if hasV3:
                rowIdx.append(baseRow)
                colIdx.append(baseLp + 1)
                data.append(-1.0)

            # d(con2)/d(lp): con2 = W1^2 + W2^2 --> dW1 = 2 W1, dW2 = 2 W2
            if hasW1:
                rowIdx.append(baseRow + 1)
                colIdx.append(baseLp + 2)
                data.append(2.0 * W1)
            if hasW2:
                rowIdx.append(baseRow + 1)
                colIdx.append(baseLp + 3)
                data.append(2.0 * W2)

            # d(con3)/d(lp): handle special and general expression
            if hasW1 and hasW3 and not (hasW2 or hasW4):
                # special case: con3 = 2 W1^2 - W3
                rowIdx.append(baseRow + 2)
                colIdx.append(baseLp + 2)
                data.append(4.0 * W1)
                rowIdx.append(baseRow + 2)
                colIdx.append(baseLp + 4)
                data.append(-1.0)
            elif hasW1 or hasW2 or hasW3 or hasW4:
                # general derivatives (only add entries for present LPs)
                if hasW1:
                    rowIdx.append(baseRow + 2)
                    colIdx.append(baseLp + 2)
                    data.append(4.0 * W1 * (1.0 - W3) - 4.0 * W2 * W4)
                if hasW2:
                    rowIdx.append(baseRow + 2)
                    colIdx.append(baseLp + 3)
                    data.append(4.0 * W2 * (1.0 + W3) - 4.0 * W1 * W4)
                if hasW3:
                    rowIdx.append(baseRow + 2)
                    colIdx.append(baseLp + 4)
                    data.append(2.0 * (W2 * W2 - W1 * W1) + 2.0 * W3)
                if hasW4:
                    rowIdx.append(baseRow + 2)
                    colIdx.append(baseLp + 5)
                    data.append(2.0 * W4 - 4.0 * W1 * W2)

        # Construct sparse Jacobian matrix
        J_lp = sp.sparse.csr_matrix(
            (data, (rowIdx, colIdx)),
            shape=(self.nCon, self.nComps * self.nLP),
            dtype=self.dtype,
        )

        # Chain rule: J_x = J_lp * A_lp
        J_x = J_lp.dot(self.A_lp)

        return J_x

    def getBounds(self):
        """
        Return lower and upper bounds for the assembled constraints.

        Returns
        -------
        (ndarray, ndarray)
            Tuple of (lb, ub) arrays of length `nCon`.
        """
        return self.lb.copy(), self.ub.copy()
