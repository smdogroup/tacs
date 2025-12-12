"""
The purpose of this module is to define a class to handle design variable constraints related to laminate parameters in composite structures.

Constraints considered are placed on membrane and bending lamination parameters.

TODO: Add definitions



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
        weight_by_dv = {}
        for dv, w in zip(dvIndices, dvWeights):
            weight_by_dv[int(dv)] = w

        # Rows are indexed as (comp_index*TOTAL_LP + lp_index)
        rows = []
        cols = []
        vals = []
        for comp_i, comp in enumerate(compIDs):
            # Get the TACS element object associated with this compID
            elemObj = self.meshLoader.getElementObject(comp, 0)
            elemIndex = 0
            # Get the dvs owned by this element
            globalDvNums = elemObj.getDesignVarNums(elemIndex)

            # For each of the 6 lamination parameter slots, check if the
            # corresponding element DV index (1..6) was requested and owned
            for lp_pos in range(TOTAL_LP):
                dv_index = lp_pos + 1
                if dv_index in weight_by_dv:
                    # Ensure the element actually has this DV slot
                    if dv_index < len(globalDvNums) and globalDvNums[dv_index] in self.globalToLocalDVNums:
                        globalDVNum = globalDvNums[dv_index]
                        localDVNum = self.globalToLocalDVNums[globalDVNum]
                        rows.append(comp_i * TOTAL_LP + lp_pos)
                        cols.append(localDVNum)
                        vals.append(weight_by_dv[dv_index])

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

class SparseLamParamFullConstraint(object):
    dtype = TACSConstraint.dtype

    def __init__(self, comm, rows, cols, vals, nComps, nLP, ncols, lb=-1e20, ub=1e20):
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
        row_nnz = np.diff(self.A_lp.indptr)
        try:
            self.present = (row_nnz.reshape(self.nComps, self.nLP) > 0)
        except Exception:
            # If reshape fails, fall back to all-false presence
            self.present = np.zeros((self.nComps, self.nLP), dtype=bool)

        # Save bound information (apply to each constraint)
        if isinstance(lb, np.ndarray) and len(lb) == self.nCon:
            self.lb = lb.astype(self.dtype)
        elif isinstance(lb, float) or isinstance(lb, complex):
            self.lb = np.array([lb] * self.nCon, dtype=self.dtype)
        else:
            # If an array of length nComps provided, expand to 3 per comp
            lb_arr = np.asarray(lb)
            if lb_arr.ndim == 1 and len(lb_arr) == nComps:
                self.lb = np.repeat(lb_arr.astype(self.dtype), 3)
            else:
                self.lb = np.array([lb] * self.nCon, dtype=self.dtype)

        if isinstance(ub, np.ndarray) and len(ub) == self.nCon:
            self.ub = ub.astype(self.dtype)
        elif isinstance(ub, float) or isinstance(ub, complex):
            self.ub = np.array([ub] * self.nCon, dtype=self.dtype)
        else:
            ub_arr = np.asarray(ub)
            if ub_arr.ndim == 1 and len(ub_arr) == nComps:
                self.ub = np.repeat(ub_arr.astype(self.dtype), 3)
            else:
                self.ub = np.array([ub] * self.nCon, dtype=self.dtype)

    def evalCon(self, x):
        # Compute lamination parameters globally
        lp_global = self.comm.allreduce(self.A_lp.dot(x))
        # lp_global is length (nComps * nLP)
        lp_mat = lp_global.reshape(self.nComps, self.nLP)

        print(lp_mat)

        con_list = []
        for i in range(self.nComps):
            lp = lp_mat[i]

            # Map indices (values default to 0.0 if not present)
            V1 = lp[0]
            V3 = lp[1]
            W1 = lp[2]
            W2 = lp[3]
            W3 = lp[4]
            W4 = lp[5]

            hasV1 = bool(self.present[i, 0])
            hasV3 = bool(self.present[i, 1])
            hasW1 = bool(self.present[i, 2])
            hasW2 = bool(self.present[i, 3])
            hasW3 = bool(self.present[i, 4])
            hasW4 = bool(self.present[i, 5])

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

            con_list.extend([con1, con2, con3])

        print("final con_list")
        print(con_list)
        return np.asarray(con_list, dtype=self.dtype)

    def evalConSens(self, x=None):
        # Compute lamination parameters globally
        lp_global = self.comm.allreduce(self.A_lp.dot(x))
        lp_mat = lp_global.reshape(self.nComps, self.nLP)

        # Build Jacobian of constraints w.r.t lamination parameters (sparse)
        # J_lp has shape (nCon, nComps*nLP)
        data = []
        row_idx = []
        col_idx = []

        for comp_i in range(self.nComps):
            lp = lp_mat[comp_i]
            base_row = comp_i * 3
            base_lp = comp_i * self.nLP
            # Values
            V1 = lp[0]
            V3 = lp[1]
            W1 = lp[2]
            W2 = lp[3]
            W3 = lp[4]
            W4 = lp[5]

            hasV1 = bool(self.present[comp_i, 0])
            hasV3 = bool(self.present[comp_i, 1])
            hasW1 = bool(self.present[comp_i, 2])
            hasW2 = bool(self.present[comp_i, 3])
            hasW3 = bool(self.present[comp_i, 4])
            hasW4 = bool(self.present[comp_i, 5])

            # d(con1)/d(lp): con1 = 2 V1^2 - V3 --> dV1 = 4 V1, dV3 = -1
            if hasV1:
                row_idx.append(base_row)
                col_idx.append(base_lp + 0)
                data.append(4.0 * V1)
            if hasV3:
                row_idx.append(base_row)
                col_idx.append(base_lp + 1)
                data.append(-1.0)

            # d(con2)/d(lp): con2 = W1^2 + W2^2 --> dW1 = 2 W1, dW2 = 2 W2
            if hasW1:
                row_idx.append(base_row + 1)
                col_idx.append(base_lp + 2)
                data.append(2.0 * W1)
            if hasW2:
                row_idx.append(base_row + 1)
                col_idx.append(base_lp + 3)
                data.append(2.0 * W2)

            # d(con3)/d(lp): handle special and general expression
            if hasW1 and hasW3 and not (hasW2 or hasW4):
                # special case: con3 = 2 W1^2 - W3
                row_idx.append(base_row + 2)
                col_idx.append(base_lp + 2)
                data.append(4.0 * W1)
                row_idx.append(base_row + 2)
                col_idx.append(base_lp + 4)
                data.append(-1.0)
            elif hasW1 or hasW2 or hasW3 or hasW4:
                # general derivatives (only add entries for present LPs)
                if hasW1:
                    row_idx.append(base_row + 2)
                    col_idx.append(base_lp + 2)
                    data.append(4.0 * W1 * (1.0 - W3) - 4.0 * W2 * W4)
                if hasW2:
                    row_idx.append(base_row + 2)
                    col_idx.append(base_lp + 3)
                    data.append(4.0 * W2 * (1.0 + W3) - 4.0 * W1 * W4)
                if hasW3:
                    row_idx.append(base_row + 2)
                    col_idx.append(base_lp + 4)
                    data.append(2.0 * (W2 * W2 - W1 * W1) + 2.0 * W3)
                if hasW4:
                    row_idx.append(base_row + 2)
                    col_idx.append(base_lp + 5)
                    data.append(2.0 * W4 - 4.0 * W1 * W2)

        J_lp = sp.sparse.csr_matrix((data, (row_idx, col_idx)), shape=(self.nCon, self.nComps * self.nLP), dtype=self.dtype)

        # Chain rule: J_x = J_lp * A_lp
        J_x = J_lp.dot(self.A_lp)

        return J_x

    def getBounds(self):
        return self.lb.copy(), self.ub.copy()
