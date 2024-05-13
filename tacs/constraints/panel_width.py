"""
Author:
    - Sean P Engelstad, Alasdair Christison Gray

This class implements a constraint which enforces the panel
width design variable values passed to elements using the GPBladeStiffenedShell
constitutive model to be consistent with the true width of the panel they are
a part of.

.. note:: This class should be created using the
    :meth:`pyTACS.createPanelLengthConstraint <tacs.pytacs.pyTACS.createPanelWidthConstraint>` method.
"""

# ==============================================================================
# Standard Python modules
# ==============================================================================

# ==============================================================================
# External Python modules
# ==============================================================================
import numpy as np
import scipy as sp

# ==============================================================================
# Extension modules
# ==============================================================================
from .panel_length import *


class PanelWidthConstraint(PanelLengthConstraint):
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
        Use pyTACS.createPanelWidthConstraint instead.

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

        super(PanelWidthConstraint, self).__init__(
            name, assembler, comm, outputViewer, meshLoader, options
        )

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
        widthAxes = []
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

                # now this refAxis is the panel length axis, compute the width axis which is normal to the length + panel normal axis
                widthAxis = np.cross(refAxis, panelNormal)
                widthAxis /= np.linalg.norm(widthAxis)
                widthAxes.append(widthAxes)

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
                # The DV jacobian on the proc that owns this component's DV
                # will have a -1 in the row corresponding to this constraint
                # and the column corresponding to the local DV index
                dvJacRows[dvProcs[-1]].append(constraintInd)
                dvJacCols[dvProcs[-1]].append(dvLocalInds[-1])
                dvJacVals[dvProcs[-1]].append(-1.0)

                for ii in range(len(boundaryNodeIDs[compID])):
                    # the coordinate jacobian on the proc that owns this node
                    # will have 3 entries in the row corresponding to this
                    # constraint and the columns corresponding to the local
                    # node index on the proc
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
                "lengthAxes": refAxes,
                "widthAxes": widthAxes,
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
        >>> # Result will look like (if PanelWidthConstraint has name of 'c1'):
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
                        widthAxis = self.constraintList[conName]["widthAxes"][ii]
                        DVProc = self.constraintList[conName]["dvProcs"][ii]
                        DVLocalInd = self.constraintList[conName]["dvLocalInds"][ii]
                        constraintValues[ii] = (
                            self.computePanelLength(boundaryPoints, widthAxis)
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
                        widthAxis = self.constraintList[conName]["widthAxes"][ii]
                        LSens = self.computePanelLengthSens(boundaryPoints, widthAxis)
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
