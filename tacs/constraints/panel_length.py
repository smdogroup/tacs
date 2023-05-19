"""
==============================================================================
Panel Length Constraint
==============================================================================
@File    :   panel_length.py
@Date    :   2023/04/23
@Author  :   Alasdair Christison Gray
@Description : This class implements a constraint which enforces the panel
length design variable values passed to elements using the BladeStiffenedShell
constitutive model to be consistent with the true length of the panel they are
a part of.
"""

# ==============================================================================
# Standard Python modules
# ==============================================================================

# ==============================================================================
# External Python modules
# ==============================================================================
from numba import njit
import numpy as np
import scipy as sp

# ==============================================================================
# Extension modules
# ==============================================================================
from tacs.constraints.base import TACSConstraint


@njit(cache=True, fastmath=True)
def computePanelLength(points, direction):
    """Given the sorted points around the perimeter of a panel, compute the length of the panel in a given direction

    Note: This function is approximate, it works best when the length direction is close to parallel with the panel

    Parameters
    ----------
    points : n x 3 array
        Coordinates of the perimeter points of the panel, in sorted order, so that points[i] - points[i-1] is a vector
        along the perimeter
    direction : length 3 array
        Direction in which to compute the panel length

    Returns
    -------
    float or complex
        The panel length in the given direction
    """
    numPoints = points.shape[0]
    length = 0.0
    for pointInd in range(numPoints):
        for edgeInd in range(numPoints):
            startInd = edgeInd
            endInd = (edgeInd + 1) % numPoints
            # We only need to check edges that are not adjacent to the current point
            if not (startInd == pointInd or endInd == pointInd):
                # Find the intersection of a line through the point in the given
                # direction with a line along the current edge
                edge = points[endInd] - points[startInd]
                mat = np.stack((direction, -edge), axis=1)
                rhs = points[startInd] - points[pointInd]
                sol, res, _, _ = np.linalg.lstsq(mat, rhs)
                alpha = sol[0]
                beta = sol[1]
                # Only compute the length if the intersection occurs within the true bounds of the edge
                if np.real(beta) - 1e-12 <= 1 and np.real(beta) + 1e-12 >= 0:
                    # intersectionPoint = points[startInd] + beta * edge
                    # newLength = np.sqrt(np.sum((intersectionPoint - points[pointInd]) ** 2))
                    newLength = np.sqrt(np.sum((alpha * direction) ** 2))
                    if (
                        np.real(newLength) > np.real(length)
                        and np.real(newLength) > res
                    ):
                        length = newLength
    return length


@njit(cache=True, fastmath=True)
def computePanelLengthSens(points, direction):
    """Given the sorted points around the perimeter of a panel, compute the
    sensitivitiy of the length of the panel in a given direction with
    respect to the coordinates of the points

    Parameters
    ----------
    points : n x 3 array
        Coordinates of the perimeter points of the panel, in sorted order, so that points[i] - points[i-1] is a vector
        along the perimeter
    direction : length 3 array
        Direction in which to compute the panel length

    Returns
    -------
    n x 3 array
        Panel length sensitivities
    """
    sens = np.zeros_like(points)
    pointsPert = np.zeros_like(points, dtype=np.complex128)
    pointsPert[:] = points[:]
    for ii in range(points.shape[0]):
        for jj in range(points.shape[1]):
            pointsPert[ii, jj] += 1e-200j
            sens[ii, jj] = np.imag(computePanelLength(pointsPert, direction)) * 1e200
            pointsPert[ii, jj] -= 1e-200j
    return sens


def simplifyPoly(nodeIDs, nodes, angleTol=18.0):
    """
    Take a (closed) chain of nodes and remove any nodes that turn by less than angleTol degrees.
    This simplifies the polygon by leaving only "sharp" corners
    """
    cont = True
    while cont:
        newNodes = []
        newNodeIDs = []
        for i in range(len(nodes)):
            im1 = i - 1
            ip1 = i + 1
            if i == 0:
                im1 = len(nodes) - 1
            if i == len(nodes) - 1:
                ip1 = 0

            v1 = nodes[ip1] - nodes[i]
            v2 = nodes[im1] - nodes[i]
            # Angle between vectors
            arg = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
            arg = max(-1, arg)
            arg = min(1, arg)
            theta = np.arccos(arg)

            if theta < np.pi * (1 - angleTol / 180):
                newNodes.append(nodes[i])
                newNodeIDs.append(nodeIDs[i])

        if len(newNodes) == len(nodes):
            cont = False

        nodes = np.array(newNodes)
        nodeIDs = np.array(newNodeIDs)

    return list(nodeIDs), nodes


class PanelLengthConstraint(TACSConstraint):
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
        Use pyTACS.createPanelLengthConstraint instead.

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

        # Create a map from the global DV index to the proc index that owns it and the local index on that proc
        self.globalToLocalDVNumsOnProc = self.comm.gather(
            self.globalToLocalDVNums, root=0
        )
        self.DVMap = {}
        if self.rank == 0:
            for procInd in range(self.comm.size):
                for globalInd, localInd in self.globalToLocalDVNumsOnProc[
                    procInd
                ].items():
                    self.DVMap[globalInd] = {"proc": procInd, "localInd": localInd}

        # Now create the same thing for the nodes
        nodeDict = self.meshLoader.getGlobalToLocalNodeIDDict()
        nodeDicts = self.comm.gather(nodeDict, root=0)
        self.nodeMap = {}
        if self.rank == 0:
            for procInd in range(self.comm.size):
                for globalInd, localInd in nodeDicts[procInd].items():
                    self.nodeMap[globalInd] = {"proc": procInd, "localInd": localInd}

        # Store the number of DVs and nodes on each proc
        self.numLocalDVs = self.getNumDesignVars()

        self.computePanelLength = computePanelLength
        self.computePanelLengthSens = computePanelLengthSens

        # Store flags for whether or not we need to recompute the constraints and derivatives
        self.constraintsUpToDate = {}
        self.funcs = {}
        self.constraintsSensUpToDate = {}
        self.funcsSens = {}

    def setDesignVars(self, x):
        """Wrapper around the base class setDesignVars method that also sets
        the flag to indicate that the constraints need to be recomputed
        """
        TACSConstraint.setDesignVars(self, x)
        for key in self.constraintsUpToDate:
            self.constraintsUpToDate[key] = False

    def setNodes(self, Xpts):
        """Wrapper around the base class setNodes method that also sets
        the flag to indicate that the constraints and their derivatives
        need to be recomputed
        """
        TACSConstraint.setNodes(self, Xpts)
        for key in self.constraintsUpToDate:
            self.constraintsUpToDate[key] = False
            self.constraintsSensUpToDate[key] = False

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
            compIDs = range(nComps)

        dvGlobalInds = []
        dvProcs = []
        dvLocalInds = []
        boundaryNodeGlobalInds = []
        boundaryNodeLocalInds = []
        boundaryNodeLocalProcs = []
        refAxes = []
        dvJacRows = []
        dvJacCols = []
        dvJacVals = []
        coordJacRows = []
        coordJacCols = []
        for ii in range(self.comm.size):
            dvJacRows.append([])
            dvJacCols.append([])
            dvJacVals.append([])
            coordJacRows.append([])
            coordJacCols.append([])

        # Get the boundary node IDs for each component
        boundaryNodeIDs, boundaryNodeCoords = self._getComponentBoundaryNodes(compIDs)

        if self.rank == 0:
            constraintInd = 0
            for compID in compIDs:
                # Get the TACS element object associated with this compID to
                # get the ref axis
                elemObj = self.meshLoader.getElementObject(compID, 0)
                transObj = elemObj.getTransform()
                try:
                    refAxis = transObj.getRefAxis()
                except AttributeError as e:
                    raise AttributeError(
                        f"The elements in component {self.meshLoader.compDescripts[compID]} do not have a reference axis. Please define one by using the 'ShellRefAxisTransform' class with your elements"
                    ) from e

                # For a more accurate length calculation, roject the ref axis
                # onto the "average" plane of the baseline panel geometry by
                # using an SVD to compute a normal vector
                centroid = np.mean(boundaryNodeCoords[compID], axis=0, keepdims=True)
                centredPoints = boundaryNodeCoords[compID] - centroid
                _, _, VT = np.linalg.svd(centredPoints, full_matrices=False)
                panelNormal = VT[-1]
                refAxis -= np.dot(refAxis, panelNormal) * panelNormal
                refAxis /= np.linalg.norm(refAxis)
                refAxes.append(refAxis)

                # Now figure out where the DV for this component lives
                globalDvNums = elemObj.getDesignVarNums(0)
                dvGlobalInds.append(globalDvNums[dvIndex])
                dvProcs.append(self.DVMap[globalDvNums[dvIndex]]["proc"])
                dvLocalInds.append(self.DVMap[globalDvNums[dvIndex]]["localInd"])

                # Do the same for the boundary nodes, this is a little more
                # complicated because each node may be on a different proc
                GlobalInds = []
                LocalInds = []
                LocalProcs = []
                for nodeID in boundaryNodeIDs[compID]:
                    GlobalInds.append(nodeID)
                    LocalInds.append(self.nodeMap[nodeID]["localInd"])
                    LocalProcs.append(self.nodeMap[nodeID]["proc"])
                boundaryNodeGlobalInds.append(GlobalInds)
                boundaryNodeLocalInds.append(LocalInds)
                boundaryNodeLocalProcs.append(LocalProcs)

                # Figure out the jacobian sparsity for each proc
                # The DV jacobian on the proc that owns this component's DV will have a -1 in the row corresponding to this constraint and the column corresponding to the local DV index
                dvJacRows[dvProcs[-1]].append(constraintInd)
                dvJacCols[dvProcs[-1]].append(dvLocalInds[-1])
                dvJacVals[dvProcs[-1]].append(-1.0)

                for ii in range(len(boundaryNodeIDs[compID])):
                    # the coordinate jacobian on the proc that owns this node will have 3 entries in the row corresponding to this constraint and the columns corresponding to the local node index on the proc
                    proc = boundaryNodeLocalProcs[-1][ii]
                    localNodeInd = boundaryNodeLocalInds[-1][ii]
                    coordJacRows[proc] += [constraintInd] * 3
                    coordJacCols[proc] += [
                        3 * localNodeInd,
                        3 * localNodeInd + 1,
                        3 * localNodeInd + 2,
                    ]

                constraintInd += 1

            # Add the constraint to the constraint list, we need to store:
            # - The size of this constraint
            # - The global index of each DV
            # - The proc number that is in charge of each DV
            # - The local index of each DV on the proc that is in charge of it
            # - The global boundary node IDs for each component
            # - The proc that owns each boundary node for each component
            # - The local index of each boundary node on the proc that owns it
            # - The reference axis direction for each component
            # - The DV jacobian for each proc (because it's constant)
            # - The coordinate jacobian sparsity pattern for each proc (because it's also constant)
            self.constraintList[conName] = {
                "nCon": len(compIDs),
                "compIDs": compIDs,
                "dvGlobalInds": dvGlobalInds,
                "dvProcs": dvProcs,
                "dvLocalInds": dvLocalInds,
                "boundaryNodeGlobalInds": boundaryNodeGlobalInds,
                "boundaryNodeLocalInds": boundaryNodeLocalInds,
                "boundaryNodeLocalProcs": boundaryNodeLocalProcs,
                "refAxes": refAxes,
            }
        else:
            self.constraintList[conName] = {"nCon": len(compIDs)}
            dvJacRows = None
            dvJacCols = None
            dvJacVals = None
            coordJacRows = None
            coordJacCols = None

        # These constraints are linear w.r.t the DVs so we can just precompute
        # the DV jacobian for each proc
        dvJacRows = self.comm.scatter(dvJacRows, root=0)
        dvJacCols = self.comm.scatter(dvJacCols, root=0)
        dvJacVals = self.comm.scatter(dvJacVals, root=0)
        coordJacRows = self.comm.scatter(coordJacRows, root=0)
        coordJacCols = self.comm.scatter(coordJacCols, root=0)
        self.constraintList[conName]["dvJac"] = sp.sparse.csr_matrix(
            (dvJacVals, (dvJacRows, dvJacCols)),
            shape=(len(compIDs), self.numLocalDVs),
        )

        # The constraints aren't linear w.r.t the coordinates, but the
        # sparsity pattern is constant so we can store that
        self.constraintList[conName]["coordJacRows"] = coordJacRows
        self.constraintList[conName]["coordJacCols"] = coordJacCols

        self.constraintsUpToDate[conName] = False
        self.constraintsSensUpToDate[conName] = False
        success = True

        return success

    def getConstraintBounds(self, bounds, evalCons=None):
        """
        Get bounds for constraints. The constraints corresponding to the strings in
        `evalCons` are evaluated and updated into the provided
        dictionary.

        The panel length constraints are equality constraints so both the upper and lower bounds are zero

        Parameters
        ----------
        bounds : dict
            Dictionary into which the constraint bounds are saved.
            Bounds will be saved as a tuple: (lower, upper)
        evalCons : iterable object containing strings.
            If not none, use these constraints to evaluate.

        Examples
        --------
        >>> conBounds = {}
        >>> tacsConstraint.getConstraintBounds(conBounds, 'LE_SPAR')
        >>> conBounds
        >>> # Result will look like (if TACSConstraint has name of 'c1'):
        >>> # {'c1_LE_SPAR': (array([-1e20]), array([1e20]))}
        """
        # Check if user specified which constraints to output
        # Otherwise, output them all
        evalCons = self._processEvalCons(evalCons)

        # Loop through each requested constraint set
        for conName in evalCons:
            key = f"{self.name}_{conName}"
            nCon = self.constraintList[conName]["nCon"]
            bounds[key] = (np.zeros(nCon), np.zeros(nCon))

    def getConstraintSizes(self, sizes, evalCons=None):
        """
        Get number for constraint equations in each set.
        The constraints corresponding to the strings in `evalCons`
        are evaluated and updated into the provided dictionary.

        Parameters
        ----------
        sizes : dict
            Dictionary into which the constraint sizes are saved.
        evalCons : iterable object containing strings.
            If not none, use these constraints to evaluate.

        Examples
        --------
        >>> conSizes = {}
        >>> tacsConstraint.getConstraintSizes(conSizes, 'LE_SPAR')
        >>> funconSizescs
        >>> # Result will look like (if TACSConstraint has name of 'c1'):
        >>> # {'c1_LE_SPAR': 10}
        """
        # Check if user specified which constraints to output
        # Otherwise, output them all
        evalCons = self._processEvalCons(evalCons)

        # Loop through each requested constraint set
        for conName in evalCons:
            key = f"{self.name}_{conName}"
            sizes[key] = self.constraintList[conName]["nCon"]

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
        >>> adjConstraint.evalConstraints(funcs, 'LE_SPAR')
        >>> funcs
        >>> # Result will look like (if PanelLengthConstraint has name of 'c1'):
        >>> # {'c1_LE_SPAR': array([1.325, 2.1983645, 3.1415926, ...])}
        """

        # Check if user specified which constraints to output
        # Otherwise, output them all
        evalCons = self._processEvalCons(evalCons, ignoreMissing)

        # We will compute everything on the root proc and then broadcast the
        # results, this is definitely not the most efficient way to do this
        # so we may want to re-do it in future, but this works for now

        DVs = None
        nodes = None

        # Loop through each requested constraint set
        for conName in evalCons:
            key = f"{self.name}_{conName}"
            nCon = self.constraintList[conName]["nCon"]
            if self.constraintsUpToDate[conName]:
                funcs[key] = np.copy(self.funcs[key])
            else:
                constraintValues = np.zeros(nCon, dtype=self.dtype)
                if nodes is None:
                    # Get all of the DVs and nodes on the root proc
                    DVs = self.comm.gather(self.x.getArray(), root=0)
                    nodes = self.comm.gather(self.Xpts.getArray(), root=0)
                    if self.rank == 0:
                        for ii in range(self.comm.size):
                            nodes[ii] = nodes[ii].reshape(-1, 3)
                if self.rank == 0:
                    for ii in range(nCon):
                        numPoints = len(
                            self.constraintList[conName]["boundaryNodeGlobalInds"][ii]
                        )
                        boundaryPoints = np.zeros((numPoints, 3), dtype=self.dtype)
                        for jj in range(numPoints):
                            nodeProc = self.constraintList[conName][
                                "boundaryNodeLocalProcs"
                            ][ii][jj]
                            localInd = self.constraintList[conName][
                                "boundaryNodeLocalInds"
                            ][ii][jj]
                            boundaryPoints[jj] = nodes[nodeProc][localInd]
                        refAxis = self.constraintList[conName]["refAxes"][ii]
                        DVProc = self.constraintList[conName]["dvProcs"][ii]
                        DVLocalInd = self.constraintList[conName]["dvLocalInds"][ii]
                        constraintValues[ii] = (
                            self.computePanelLength(boundaryPoints, refAxis)
                            - DVs[DVProc][DVLocalInd]
                        )
                funcs[key] = self.comm.bcast(constraintValues, root=0)
                self.funcs[key] = np.copy(funcs[key])
                self.constraintsUpToDate[conName] = True

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

        # We will compute everything on the root proc and then broadcast the results,
        # this is definitely not the most efficient way to do this so we may want to
        # re-do it in future, but this works for now

        nodes = None

        # Get number of nodes coords on this proc
        nCoords = self.getNumCoordinates()

        # Loop through each requested constraint set
        for conName in evalCons:
            key = f"{self.name}_{conName}"
            nCon = self.constraintList[conName]["nCon"]
            funcsSens[key] = {}
            funcsSens[key][self.varName] = self.constraintList[conName]["dvJac"]
            if self.constraintsSensUpToDate[conName]:
                funcsSens[key][self.coordName] = self.funcsSens[key].copy()
            else:
                # Get all of the nodes on the root proc
                if nodes is None:
                    nodes = self.comm.gather(self.Xpts.getArray(), root=0)
                    if self.rank == 0:
                        for ii in range(self.comm.size):
                            nodes[ii] = nodes[ii].reshape(-1, 3)
                funcsSens[key][self.coordName] = sp.sparse.csr_matrix(
                    (nCon, nCoords), dtype=self.dtype
                )
                if self.rank == 0:
                    coordJacVals = []
                    for ii in range(self.comm.size):
                        coordJacVals.append([])
                    for ii in range(nCon):
                        numPoints = len(
                            self.constraintList[conName]["boundaryNodeGlobalInds"][ii]
                        )
                        boundaryPoints = np.zeros((numPoints, 3), dtype=self.dtype)

                        for jj in range(numPoints):
                            nodeProc = self.constraintList[conName][
                                "boundaryNodeLocalProcs"
                            ][ii][jj]
                            localInd = self.constraintList[conName][
                                "boundaryNodeLocalInds"
                            ][ii][jj]
                            boundaryPoints[jj] = nodes[nodeProc][localInd]
                        refAxis = self.constraintList[conName]["refAxes"][ii]
                        LSens = self.computePanelLengthSens(boundaryPoints, refAxis)
                        for jj in range(numPoints):
                            nodeProc = self.constraintList[conName][
                                "boundaryNodeLocalProcs"
                            ][ii][jj]
                            coordJacVals[nodeProc] += LSens[jj].tolist()
                else:
                    coordJacVals = None
                coordJacVals = self.comm.scatter(coordJacVals, root=0)
                coordJacRows = self.constraintList[conName]["coordJacRows"]
                coordJacCols = self.constraintList[conName]["coordJacCols"]
                self.funcsSens[key] = sp.sparse.csr_matrix(
                    (coordJacVals, (coordJacRows, coordJacCols)),
                    (nCon, nCoords),
                    dtype=self.dtype,
                )
                funcsSens[key][self.coordName] = self.funcsSens[key].copy()
                self.constraintsSensUpToDate[conName] = True

    def _getComponentBoundaryNodes(self, compIDs):
        """For a given list of components, find the nodes on the boundaries of
        each of the components.

        The results are broadcast to all procs

        Parameters
        ----------
        compIDs: list[int]
            List of compIDs to find boundaries of.

        Returns
        --------
        dict[int, list[int]]
            Dictionary where dict[compID] = sorted list of nodeIDs on the
            boundary of the component
        dict[int, np.array]
            Dictionary where dict[compID] = array of node coordinates on the
            boundary of the component
        """

        boundaryNodeIDs = {}
        boundaryNodeCoords = {}

        if self.rank == 0:
            for compID in compIDs:
                allEdges = set()
                dupEdges = set()
                compConn = self.meshLoader.getConnectivityForComp(
                    compID, nastranOrdering=False
                )
                # Go over all the elements in the component and add their edges to the sets of all and possibly duplicate edges
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

                            # Skip degenerate edges
                            if key[0] != key[1]:
                                # Either add to allEdges or dupEdges depending on whether we've seen this edge before
                                if key not in allEdges:
                                    allEdges.add(key)
                                else:
                                    dupEdges.add(key)
                # Now get a list of all the edges that aren't duplicated, these are the boundary edges
                boundaryEdges = list(allEdges - dupEdges)

                # Create a nodeToElem Pointer using a dictionary:
                nodeToElem = {}
                for iEdge in range(len(boundaryEdges)):
                    edge = boundaryEdges[iEdge]
                    for ii in range(2):
                        if edge[ii] in nodeToElem:
                            nodeToElem[edge[ii]].append(iEdge)
                        else:
                            nodeToElem[edge[ii]] = [iEdge]

                # Now check that each nodeToElem has a length of
                # 2. This means we have a chance it is a
                # closed curve
                for key in nodeToElem:
                    if len(nodeToElem[key]) != 2:
                        raise ValueError(
                            "The topology of the geometry associated with "
                            "a constitutive object is not manifold "
                            "(There is a node with three or more edges "
                            "attached. This constitutive object cannot "
                            "use a panel-type constitutive object. "
                            f"CompIDs are: {repr(compIDs)}"
                        )

                # Now we will "order" the edges if possible. This
                # will also allow us to detect multiple loops
                # which isn't allowed, or a non-manifold local
                # geometry - ie. in this context a node connected
                # to three edges. This is also not allowed.

                nodeChain = [boundaryEdges[0][0], boundaryEdges[0][1]]
                cont = True
                curElem = 0
                while cont:
                    # We arbitrarily pick the first 'element'
                    # (edge) containing the first two nodes of our
                    # chain. Next step is to find the next element
                    # and node in the chain:
                    nextElems = nodeToElem[nodeChain[-1]]
                    # Get the 'nextElem' that isn't the current
                    # one
                    if nextElems[0] == curElem:
                        nextElem = nextElems[1]
                    else:
                        nextElem = nextElems[0]

                    # Now nextElem is the next in the chain. Get
                    # the nodes for this elem:
                    nextNodes = boundaryEdges[nextElem]

                    # Append the node that isn't the last one
                    # (that is already in the chain)
                    if nextNodes[0] == nodeChain[-1]:
                        nodeChain.append(nextNodes[1])
                    else:
                        nodeChain.append(nextNodes[0])

                    # Exit condition:
                    if nodeChain[-1] == nodeChain[0]:
                        # We've made it all the way around!
                        cont = False

                    # Set current element
                    curElem = nextElem

                # Now check that we've *actually* used all of our
                # nodes. Since we've determined it is manifold,
                # this must mean we have multiple loops which
                # *also* isn't allowed.
                if len(nodeChain) - 1 != len(boundaryEdges):
                    raise ValueError(
                        "Detected more than one closed loop for "
                        "constitutive object. This is not allowed. "
                        "This constitutive object cannot use a "
                        "panel-type constitutive object. "
                        f"CompIDs are: {repr(compIDs)}"
                    )
                nodeChain = nodeChain[:-1]
                nodeChainCoords = self.meshLoader.getBDFNodes(
                    nodeChain, nastranOrdering=False
                )
                nodeIDs, coords = simplifyPoly(nodeChain, nodeChainCoords)
                boundaryNodeIDs[compID] = nodeIDs
                boundaryNodeCoords[compID] = coords

        return self.comm.bcast(boundaryNodeIDs, root=0), self.comm.bcast(
            boundaryNodeCoords, root=0
        )
