"""
The main purpose of this class is to constrain design variables step sizes across adjacent components.
The formulation is a linear constraint that takes the following form:

    c = dv_i - dv_j

Where dv_i and dv_j are two design variables in adjacent components.

A common example of this is ensuring that thicknesses vary to abruptly across panel boundaries:

    -delta_t < t_i - t_j < delta_t

.. note:: This class should be created using the
    :meth:`pyTACS.createAdjacencyConstraint <tacs.pytacs.pyTACS.createAdjacencyConstraint>` method.
"""

# =============================================================================
# Imports
# =============================================================================
import scipy as sp

from tacs.constraints.base import TACSConstraint, SparseLinearConstraint


class AdjacencyConstraint(TACSConstraint):
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
        Use pyTACS.createAdjacencyConstraint instead.

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

        # Create a list of all adjacent components on root proc
        self._initializeAdjacencyList()

    def _initializeAdjacencyList(self):
        """
        Create a list of all components with common edges.
        """

        if self.comm.rank == 0:
            edgeToFace = {}
            nComp = self.meshLoader.getNumComponents()
            for compID in range(nComp):
                compConn = self.meshLoader.getConnectivityForComp(
                    compID, nastranOrdering=False
                )
                for elemConn in compConn:
                    nnodes = len(elemConn)
                    if nnodes >= 2:
                        for j in range(nnodes):
                            nodeID1 = elemConn[j]
                            nodeID2 = elemConn[(j + 1) % nnodes]

                            if nodeID1 < nodeID2:
                                key = (nodeID1, nodeID2)
                            else:
                                key = (nodeID2, nodeID1)

                            if key not in edgeToFace:
                                edgeToFace[key] = [compID]
                            elif compID not in edgeToFace[key]:
                                edgeToFace[key].append(compID)

            # Now we loop back over each element and each edge. By
            # using the edgeToFace dictionary, we can now determine
            # which components IDs (jComp) are connected to the
            # current component ID (iComp).
            self.adjacentComps = []

            for edgeKey in edgeToFace:
                if len(edgeToFace[edgeKey]) >= 2:
                    for i, iComp in enumerate(edgeToFace[edgeKey][:-1]):
                        for jComp in edgeToFace[edgeKey][i + 1 :]:
                            if iComp < jComp:
                                dvKey = (iComp, jComp)
                            else:
                                dvKey = (jComp, iComp)
                            if dvKey not in self.adjacentComps:
                                self.adjacentComps.append(dvKey)

        else:
            self.adjacentComps = None

        # Wait for root
        self.comm.barrier()

    def addConstraint(self, conName, compIDs=None, lower=-1e20, upper=1e20, dvIndex=0):
        """
        Generic method to adding a new constraint set for TACS.

        Parameters
        ----------
        conName : str
            The user-supplied name for the constraint set. This will
            typically be a string that is meaningful to the user

        compIDs: list[int] or None
            List of compIDs to select. If None, all compIDs will be selected. Defaults to None.

        lower: float or complex
            lower bound for constraint. Defaults to -1e20.

        upper: float or complex
            upper bound for constraint. Defaults to 1e20.

        dvIndex : int
            Index number of element DV to be used in constraint. Defaults to 0.

        """
        if compIDs is not None:
            # Make sure CompIDs is flat and get element numbers on each proc corresponding to specified compIDs
            compIDs = self._flatten(compIDs)
        else:
            nComps = self.meshLoader.getNumComponents()
            compIDs = range(nComps)

        constrObj = self._createConstraint(dvIndex, compIDs, lower, upper)
        if constrObj.nCon > 0:
            self.constraintList[conName] = constrObj
            success = True
        else:
            self._TACSWarning(
                f"No adjacent components found in `compIDs`. Skipping {conName}."
            )
            success = False

        return success

    def _createConstraint(self, dvIndex, compIDs, lbound, ubound):
        size = self.comm.size
        rank = self.comm.rank
        # Gather the dv mapping from each proc
        globalToLocalDVNumsOnProc = self.comm.gather(self.globalToLocalDVNums, root=0)
        # Assemble constraint info on root proc
        if rank == 0:
            # Create a list of lists that will hold the sparse data info on each proc
            rowsOnProc = [[] for _ in range(size)]
            colsOnProc = [[] for _ in range(size)]
            valsOnProc = [[] for _ in range(size)]
            conCount = 0
            # Loop through all adjacent component pairs
            for compPair in self.adjacentComps:
                # Check if they are in the user provided compIDs
                if compPair[0] in compIDs and compPair[1] in compIDs:
                    # We found a new constraint
                    for i, comp in enumerate(compPair):
                        # Get the TACS element object associated with this compID
                        elemObj = self.meshLoader.getElementObject(comp, 0)
                        elemIndex = 0
                        # Get the dvs owned by this element
                        globalDvNums = elemObj.getDesignVarNums(elemIndex)
                        # Check if specified dv num is owned by each proc
                        for proc_i in range(size):
                            globalToLocalDVNums = globalToLocalDVNumsOnProc[proc_i]
                            if globalDvNums[dvIndex] in globalToLocalDVNums:
                                globalDVNum = globalDvNums[dvIndex]
                                localDVNum = globalToLocalDVNums[globalDVNum]
                                rowsOnProc[proc_i].append(conCount)
                                colsOnProc[proc_i].append(localDVNum)
                                if i == 0:
                                    valsOnProc[proc_i].append(1.0)
                                else:
                                    valsOnProc[proc_i].append(-1.0)
                                break
                    conCount += 1

        else:
            rowsOnProc = None
            colsOnProc = None
            valsOnProc = None
            conCount = 0

        # Scatter local sparse indices/values to remaining procs
        rows = self.comm.scatter(rowsOnProc, root=0)
        cols = self.comm.scatter(colsOnProc, root=0)
        vals = self.comm.scatter(valsOnProc, root=0)

        # Get local sparse matrix dimensions
        conCount = self.comm.bcast(conCount, root=0)
        nLocalDVs = self.getNumDesignVars()

        # Create linear constraint object
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
        >>> adjConstraint.evalConstraints(funcs, 'LE_SPAR')
        >>> funcs
        >>> # Result will look like (if AdjacencyConstraint has name of 'c1'):
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
        >>> adjConstraint.evalConstraintsSens(funcsSens, 'LE_SPAR')
        >>> funcsSens
        >>> # Result will look like (if AdjacencyConstraint has name of 'c1'):
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
