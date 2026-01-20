"""
This class implements a constraint which enforces the stiffener
length design variable values passed to elements using the IsoRectangularBeamStiffness
constitutive model to be consistent with the true length of the stiffener they are
a part of.

.. note:: This class should be created using the
    :meth:`pyTACS.createStiffenerLengthConstraint <tacs.pytacs.pyTACS.createStiffenerLengthConstraint>` method.
"""

# ==============================================================================
# Standard Python modules
# ==============================================================================
from collections import defaultdict, deque

# ==============================================================================
# External Python modules
# ==============================================================================
import numpy as np
import scipy as sp

# ==============================================================================
# Extension modules
# ==============================================================================
from tacs.constraints.base import TACSConstraint
from tacs import constitutive


class StiffenerLengthConstraint(TACSConstraint):
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
        Use pyTACS.createStiffenerLengthConstraint instead.

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

        return

    def addConstraint(self, conName, compIDs=None, lower=None, upper=None, dvIndex=0):
        """
        Generic method to adding a new constraint set for TACS.

        Parameters
        ----------
        conName : str
            The user-supplied name for the constraint set. This will
            typically be a string that is meaningful to the user

        compIDs: list[int] or None
            List of compIDs to apply constraints to. If None, all compIDs will be used. Defaults to None.

        lower: float or complex
            lower bound for constraint. Not used.

        upper: float or complex
            upper bound for constraint. Not used.

        dvIndex : int
            Index number of the panel length DV's. Defaults to 0.

        """
        if compIDs is not None:
            # Make sure CompIDs is flat and get element numbers on each proc corresponding to specified compIDs
            compIDs = self._flatten(compIDs)
        else:
            nComps = self.meshLoader.getNumComponents()
            compIDs = list(range(nComps))

        if lower is None:
            lower = 0.0
        if upper is None:
            upper = 0.0

        constrObj = self._createConstraint(dvIndex, compIDs, lower, upper)
        if constrObj.nCon > 0:
            self.constraintList[conName] = constrObj
            success = True
        else:
            self._TACSWarning(f"No valid `compIDs` provided. Skipping {conName}.")
            success = False

        return success

    def _createConstraint(self, dvIndex, compIDs, lbound, ubound):
        """
        Create a new constraint object for TACS.

        Parameters
        ----------
        dvIndex : int
            Index number of element DV to be used in constraint.
            Defaults to 0.

        compIDs: list[int]
            List of compIDs to select.

        lbound: float or complex
            lower bound for constraint. Defaults to 0.0.

        ubound: float or complex
            upper bound for constraint. Defaults to 1e20.

        Returns
        -------
        constraint : SparseLengthConstraint or None
            Constraint object if successful, None otherwise.

        """
        # Assemble constraint info for each component
        lengthDVs = np.zeros(len(compIDs), dtype=int)
        lengthDVsOwnerProc = np.zeros(len(compIDs), dtype=int)
        allEndNodeLocalIDs = np.zeros([len(compIDs), 2], dtype=int)
        allEndNodeOwnerProc = np.zeros([len(compIDs), 2], dtype=int)

        for conCount, comp in enumerate(compIDs):
            # Get the TACS element object associated with this compID
            elemObj = self.meshLoader.getElementObject(comp, 0)
            conObj = elemObj.getConstitutive()
            if not isinstance(conObj, constitutive.IsoRectangleBeamConstitutive):
                raise ValueError("Only Beams!")
            elemIndex = 0
            # Get the design variables owned by this element
            globalDvNums = elemObj.getDesignVarNums(elemIndex)
            globalDVNum = globalDvNums[dvIndex]
            # Map global DV number to local DV number if owned by this processor
            if globalDVNum in self.globalToLocalDVNums:
                localDVNum = self.globalToLocalDVNums[globalDVNum]
                lengthDVs[conCount] = localDVNum
                lengthDVsOwnerProc[conCount] = self.comm.rank

            # Get connectivity and find end nodes of the stiffener chain
            compConn = self.meshLoader.getConnectivityForComp(
                comp, nastranOrdering=False
            )
            endNodeGlobalIDs = self._findChainEnds(compConn)
            endNodeLocalIDs = self.meshLoader.getLocalNodeIDsFromGlobal(
                endNodeGlobalIDs, nastranOrdering=False
            )
            # Store end node information for constraint evaluation
            for end_i, nodeID in enumerate(endNodeLocalIDs):
                if nodeID >= 0:  # Valid local node ID
                    allEndNodeLocalIDs[conCount, end_i] = nodeID
                    allEndNodeOwnerProc[conCount, end_i] = self.comm.rank

        # Gather information from all processors
        lengthDVs = self.comm.allreduce(lengthDVs)
        allEndNodeLocalIDs = self.comm.allreduce(allEndNodeLocalIDs)
        allEndNodeOwnerProc = self.comm.allreduce(allEndNodeOwnerProc)

        nLocalDVs = self.getNumDesignVars()

        return SparseLengthConstraint(
            self.comm,
            lengthDVs,
            lengthDVsOwnerProc,
            allEndNodeLocalIDs,
            allEndNodeOwnerProc,
            nLocalDVs,
            lbound,
            ubound,
        )

    def _findChainEnds(self, connectivity):
        """
        Takes a list of elements (each element is a list of node IDs) and verifies
        that they form exactly one continuous topological body with exactly two ends.
        Finally, it returns the two end nodes of the chain.

        Parameters
        ----------
        connectivity : list[list[int]]
            List of element nodal connectivities

        Returns
        -------
        ends : tuple[int, int]
            (endNodeID1, endNodeID2) - the two end node IDs of the chain

        Raises
        ------
        ValueError
            If the connectivity does not form exactly one continuous body or
            if the body does not have exactly two ends
        """

        # Build adjacency map from connectivity list
        adj = defaultdict(set)
        for element in connectivity:
            # Connect consecutive nodes in each element
            for i in range(len(element) - 1):
                n1, n2 = element[i], element[i + 1]
                adj[n1].add(n2)
                adj[n2].add(n1)

        if not adj:
            raise ValueError("No nodes found in connectivity list.")

        # --- 1. Check there is exactly one continuous body ---
        visited = set()
        startNode = next(iter(adj))  # pick any node

        queue = deque([startNode])
        while queue:
            node = queue.popleft()
            if node in visited:
                continue
            visited.add(node)
            queue.extend(adj[node] - visited)

        if len(visited) != len(adj):
            raise ValueError("Connectivity list describes more than one body.")

        # --- 2. Find degree (number of connected neighbors) of each node ---
        degrees = {node: len(neighbors) for node, neighbors in adj.items()}

        # Ends should have degree == 1 (only one neighbor)
        ends = [node for node, deg in degrees.items() if deg == 1]

        if len(ends) != 2:
            raise ValueError(
                f"Body does not have exactly two ends (found {len(ends)} ends)."
            )

        return tuple(ends)

    def evalConstraints(self, funcs, evalCons=None, ignoreMissing=False):
        """
        Evaluate values for constraints. The constraints corresponding to the strings in
        evalCons are evaluated and updated into the provided dictionary.

        The same constraint arrays are returned on every proc

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
        >>> plConstraint.evalConstraints(funcs, 'LE_SPAR')
        >>> funcs
        >>> # Result will look like (if PanelLengthConstraint has name of 'c1'):
        >>> # {'c1_LE_SPAR': array([1.325, 2.1983645, 3.1415926, ...])}
        """

        # Check if user specified which constraints to output
        # Otherwise, output them all
        evalCons = self._processEvalCons(evalCons, ignoreMissing)

        # Loop through each requested constraint set
        for conName in evalCons:
            key = f"{self.name}_{conName}"
            funcs[key] = self.constraintList[conName].evalCon(
                self.x.getArray(), self.Xpts.getArray()
            )

    def evalConstraintsSens(self, funcsSens, evalCons=None):
        """This is the main routine for returning useful (sensitivity)
        information from constraint. The derivatives of the constraints
        corresponding to the strings in evalCons are evaluated and
        updated into the provided dictionary. The derivitives with
        respect to all design variables and node locations are computed.

        The sensitivities returned on each proc are a sparse m x n matrix
        where m is the number of constraints and n is the number of design
        variables or 3x the number of nodes on this proc. The matrix contains
        the sensitivities of all constraints w.r.t only the design
        variables/node coordinates on this proc.

        Parameters
        ----------
        funcsSens : dict
            Dictionary into which the derivatives are saved.
        evalCons : iterable object containing strings
            The constraints the user wants returned

        Examples
        --------
        >>> funcsSens = {}
        >>> adjConstraint.evalConstraintsSens(funcsSens, 'LE_SPAR')
        >>> funcsSens
        >>> # Result will look like (if AdjacencyConstraint has name of 'c1'):
        >>> # {'c1_LE_SPAR':{'struct':<50x242 sparse matrix of type '<class 'numpy.float64'>' with 100 stored elements in Compressed Sparse Row format>}}
        """

        # Check if user specified which constraints to output
        # Otherwise, output them all
        evalCons = self._processEvalCons(evalCons)

        # Loop through each requested constraint set
        for conName in evalCons:
            key = f"{self.name}_{conName}"
            xSens, XptSens = self.constraintList[conName].evalConSens(
                self.x.getArray(), self.Xpts.getArray()
            )
            # Get sparse Jacobian for dv sensitivity
            funcsSens[key] = {self.varName: xSens, self.coordName: XptSens}


class SparseLengthConstraint(object):
    """
    A class for handling sparse length constraints in parallel.

    This class manages the evaluation of constraints that enforce the consistency
    between design variable length values and the actual geometric length of
    stiffener elements. It handles parallel computation across multiple processors
    where different processors may own different design variables and nodes.

    The constraint ensures that the design variable representing the length of
    a stiffener matches the actual geometric distance between the end nodes
    of the stiffener chain.
    """

    dtype = TACSConstraint.dtype

    def __init__(
        self,
        comm,
        lengthLocalDVNums,
        lengthDVsOwnerProc,
        allEndNodeLocalIDs,
        allEndNodeOwnerProc,
        nLocalDVs,
        lbound,
        ubound,
    ):
        """
        Initialize the SparseLengthConstraint object.

        Parameters
        ----------
        comm : mpi4py.MPI.Intracomm
            MPI communicator for parallel operations
        lengthLocalDVNums : array_like
            Local design variable numbers for length DVs
        lengthDVsOwnerProc : array_like
            Processor ranks that own each length DV
        allEndNodeLocalIDs : array_like
            Local node IDs for end nodes of each stiffener (nCon x 2)
        allEndNodeOwnerProc : array_like
            Processor ranks that own each end node (nCon x 2)
        nLocalDVs : int
            Number of local design variables on this processor
        lbound : float or array_like
            Lower bounds for constraints
        ubound : float or array_like
            Upper bounds for constraints
        """
        # Number of constraints
        self.nCon = len(lengthLocalDVNums)
        self.lengthLocalDVNums = lengthLocalDVNums
        self.lengthDVsOwnerProc = lengthDVsOwnerProc
        self.comm = comm

        self.allEndNodeLocalIDs = allEndNodeLocalIDs
        self.allEndNodeOwnerProc = allEndNodeOwnerProc
        self.nLocalDVs = nLocalDVs

        # Save bound information
        if isinstance(lbound, np.ndarray) and len(lbound) == self.nCon:
            self.lbound = lbound.astype(self.dtype)
        elif isinstance(lbound, float) or isinstance(lbound, complex):
            self.lbound = np.array([lbound] * self.nCon, dtype=self.dtype)

        if isinstance(ubound, np.ndarray) and len(ubound) == self.nCon:
            self.ubound = ubound.astype(self.dtype)
        elif isinstance(ubound, float) or isinstance(ubound, complex):
            self.ubound = np.array([ubound] * self.nCon, dtype=self.dtype)

    def evalCon(self, x, Xpts):
        """
        Evaluate the constraint values.

        Parameters
        ----------
        x : array_like
            Design variable values
        Xpts : array_like
            Node coordinates

        Returns
        -------
        array_like
            Constraint values (difference between DV lengths and exact lengths)
        """
        Lexact = self._computeExactLength(Xpts)
        Ldv = self._getDVLengths(x)
        Ldiff = Ldv - Lexact
        return Ldiff

    def _computeExactLength(self, Xpts):
        """
        Compute the exact geometric length of each stiffener.

        Parameters
        ----------
        Xpts : array_like
            Node coordinates

        Returns
        -------
        array_like
            Exact lengths of each stiffener
        """
        stiffenerEndLocations = np.zeros([self.nCon, 2, 3], dtype=self.dtype)
        for con_i in range(self.nCon):
            for end_j in range(2):
                if self.allEndNodeOwnerProc[con_i, end_j] == self.comm.rank:
                    localNodeID = self.allEndNodeLocalIDs[con_i, end_j]
                    stiffenerEndLocations[con_i, end_j, :] = Xpts[
                        3 * localNodeID : 3 * localNodeID + 3
                    ]
        stiffenerEndLocations = self.comm.allreduce(stiffenerEndLocations)
        dX = stiffenerEndLocations[:, 1, :] - stiffenerEndLocations[:, 0, :]
        Lexact = np.sqrt(np.sum(dX * dX, axis=1))
        return Lexact

    def _getDVLengths(self, x):
        """
        Get the design variable length values for each stiffener.

        Parameters
        ----------
        x : array_like
            Design variable values

        Returns
        -------
        array_like
            Length values from design variables
        """
        Ldv = np.zeros(self.nCon, dtype=self.dtype)
        for con_i in range(self.nCon):
            if self.lengthDVsOwnerProc[con_i] == self.comm.rank:
                Ldv[con_i] = x[self.lengthLocalDVNums[con_i]]
        Ldv = self.comm.allreduce(Ldv)
        return Ldv

    def evalConSens(self, x, Xpts):
        """
        Evaluate the constraint sensitivities.

        Parameters
        ----------
        x : array_like
            Design variable values
        Xpts : array_like
            Node coordinates

        Returns
        -------
        tuple
            (xSens, XptSens) - sensitivities w.r.t. design variables and coordinates
        """
        # Initialize sparse matrix components for Jacobian
        dvJacVals = []
        dvJacRows = []
        dvJacCols = []
        coordJacVals = []
        coordJacRows = []
        coordJacCols = []
        Lexact = self._computeExactLength(Xpts)

        for con_i in range(self.nCon):
            # Sensitivity w.r.t. design variables (derivative = 1.0)
            if self.lengthDVsOwnerProc[con_i] == self.comm.rank:
                dvJacVals.append(1.0)
                dvJacRows.append(con_i)
                dvJacCols.append(self.lengthLocalDVNums[con_i])

            # Sensitivity w.r.t. node coordinates
            for end_j in range(2):
                if self.allEndNodeOwnerProc[con_i, end_j] == self.comm.rank:
                    localNodeID = self.allEndNodeLocalIDs[con_i, end_j]
                    # Derivative of length w.r.t. node coordinates
                    val = -Xpts[3 * localNodeID : 3 * localNodeID + 3] / Lexact[con_i]
                    if end_j == 0:
                        val *= -1
                    coordJacVals.extend(val)
                    coordJacRows.extend([con_i, con_i, con_i])
                    coordJacCols.extend(
                        np.arange(3 * localNodeID, 3 * localNodeID + 3, dtype=int)
                    )
        xSens = sp.sparse.csr_matrix(
            (dvJacVals, (dvJacRows, dvJacCols)),
            (self.nCon, self.nLocalDVs),
            dtype=self.dtype,
        )
        XptSens = sp.sparse.csr_matrix(
            (coordJacVals, (coordJacRows, coordJacCols)),
            (self.nCon, Xpts.size),
            dtype=self.dtype,
        )
        return xSens, XptSens

    def getBounds(self):
        """
        Get the constraint bounds.

        Returns
        -------
        tuple
            (lbound, ubound) - lower and upper bounds for constraints
        """
        return self.lbound.copy(), self.ubound.copy()
