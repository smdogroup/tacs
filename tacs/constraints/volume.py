"""
The main purpose of this class is to constrain the size of a closed volume.
Only shell and solid elements are supported for this constraint.
For shell elements, the enclosed volume MUST be manifold and water-tight (no missing/internal faces).
The formulation is a nonlinear constraint based on the nodal coordinates.

A common example of this is ensuring enough volume in the wingbox for fuel:

    vol_wing >= vol_fuel
"""

# =============================================================================
# Imports
# =============================================================================
import os
import copy

import numpy as np
import scipy as sp
import pyNastran.bdf as pn
import pyNastran.converters.tecplot.tecplot as tp

from tacs.constraints.base import TACSConstraint
from tacs.functions import EnclosedVolume


class VolumeConstraint(TACSConstraint):
    # Default options for class
    defaultOptions = {
        "volCheckTol": [
            float,
            1e-12,
            "Relative tolerance for surface closure check for shell elements.",
        ],
        "outputDir": [str, "./", "Output directory for F5 file writer."],
        "numberSolutions": [
            bool,
            True,
            "Flag for attaching solution counter index to f5 files.",
        ],
    }

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
        Use pyTACS.createVolumeConstraint instead.

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

        # Set call counter
        self.callCounter = -1

    def _initializeAdjacencyList(self):
        """
        Create a list of all components with common edges.
        """

        if self.comm.rank == 0:
            # First, create a dictionary of common edges shared by components
            edgeToFace = {}
            for elemID in self.bdfInfo.elements:
                elemInfo = self.bdfInfo.elements[elemID]
                elemConn = elemInfo.nodes
                compID = self.meshLoader.nastranToTACSCompIDDict[elemInfo.pid]
                nnodes = len(elemConn)
                if nnodes >= 2:
                    for j in range(nnodes):
                        nodeID1 = elemConn[j]
                        nodeID2 = elemConn[(j + 1) % nnodes]

                        key = (nodeID1, nodeID2)

                        if key not in edgeToFace:
                            edgeToFace[key] = [compID]
                        elif compID not in edgeToFace[key]:
                            edgeToFace[key].append(compID)

            # Now we loop back over each element and each edge. By
            # using the edgeToFace dictionary, we can now determine
            # which components IDs (jComp) are connected to the
            # current component ID (iComp).
            nComp = self.meshLoader.getNumComponents()
            self.adjacentOrientationMatch = {compID: {} for compID in range(nComp)}

            # Loop through each shared edge and determine if
            # components sharing the edge have matching orientations or not
            for edgeKey in edgeToFace:
                if len(edgeToFace[edgeKey]) >= 2:
                    # If this component shares an edge with another component
                    # that uses the same edge ordering, the normal orientations
                    # will be inconsistent
                    for i, iComp in enumerate(edgeToFace[edgeKey][:-1]):
                        for jComp in edgeToFace[edgeKey][i + 1 :]:
                            self.adjacentOrientationMatch[iComp][
                                jComp
                            ] = self.adjacentOrientationMatch[jComp][iComp] = False
                # If this component shares an edge with another component
                # but the edge ordering is reversed, the normal orientations
                # will be consistent
                flippedEdgeKey = (edgeKey[1], edgeKey[0])
                if flippedEdgeKey in edgeToFace:
                    for iComp in edgeToFace[edgeKey]:
                        for jComp in edgeToFace[flippedEdgeKey]:
                            self.adjacentOrientationMatch[iComp][
                                jComp
                            ] = self.adjacentOrientationMatch[jComp][iComp] = True

        else:
            self.adjacentOrientationMatch = None

        # Wait for root
        self.comm.barrier()

    def addConstraint(self, conName, compIDs=None, lower=0, upper=1e20):
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
            lower bound for constraint. Defaults to 0.0.

        upper: float or complex
            upper bound for constraint. Defaults to 1e20.

        """
        if compIDs is not None:
            # Make sure CompIDs is flat and get element numbers on each proc corresponding to specified compIDs
            compIDs = self._flatten(compIDs)
        else:
            nComps = self.meshLoader.getNumComponents()
            compIDs = list(range(nComps))

        constrObj = self._createConstraint(compIDs, lower, upper)
        if constrObj.nCon is not None:
            self.constraintList[conName] = constrObj
            success = True
        else:
            self._TACSWarning(
                f"No closed volume found in `compIDs`. Skipping {conName}."
            )
            success = False

        return success

    def _createConstraint(self, compIDs, lbound, ubound):
        """
        Create a new constraint object for TACS.

        Parameters
        ----------
        compIDs: list[int]
            List of compIDs to select.

        lbound: float or complex
            lower bound for constraint. Defaults to 0.0.

        ubound: float or complex
            upper bound for constraint. Defaults to 1e20.

        Returns
        -------
        constraint : ParallelVolumeConstraint or None
            Constraint object if successful, None otherwise.

        """
        # Check if elements in supplied compIDs are all shell elements or all solid elements
        elemIDs = self.meshLoader.getGlobalElementIDsForComps(
            compIDs, nastranOrdering=True
        )
        allShells = False
        allSolids = False
        for elemID in elemIDs:
            elemInfo = self.bdfInfo.elements[elemID]
            if isinstance(elemInfo, pn.cards.elements.shell.ShellElement):
                if allShells is False:
                    allShells = True

            elif isinstance(elemInfo, pn.cards.elements.solid.SolidElement):
                if allSolids is False:
                    allSolids = True

        if not (allShells or allSolids):
            self._TACSWarning("No shell or solid elements found in provided compIDs.")
            return None
        elif allShells and allSolids:
            self._TACSWarning(
                "Both shell and solid elements were found in compIDs. "
                "Only all-shell or all-solid element meshes are supported."
            )
            return None

        # Assemble constraint info on root proc
        if self.comm.rank == 0:
            properNormalCompIDs = []
            flippedNormalCompIDs = []
            # For shell element meshes we have to take an additional step.
            # All components must have consistent normal directions (normals pointing out)
            # for the volume calculation to work correctly. In general, this won't be the case.
            # So we'll need to loop through each component and figure out if their normals need
            # to be flipped. We'll also need to make sure the provided volume is closed.
            if allShells:
                # First, make sure that every component is connected to at least one other component
                for compID in compIDs:
                    adjIDs = self.adjacentOrientationMatch[compID].keys()
                    atLeastOneConnection = False
                    for adjID in adjIDs:
                        if adjID in compIDs:
                            atLeastOneConnection = True
                            break
                    if len(compIDs) > 1 and atLeastOneConnection is False:
                        self._TACSWarning(
                            "Provided compID's do not form a closed volume."
                        )
                        return None

                # Assume that the first component is proper (we'll test this assumption later)
                currCompID = compIDs[0]
                properNormalCompIDs.append(currCompID)
                # Sort the remaining components based on their connectivity and adjacentOrientationMatch dict
                self._sortIDs(
                    currCompID, compIDs, properNormalCompIDs, flippedNormalCompIDs
                )

                # Now we check to make sure the volume is closed
                # This can be checked by integrating the normal over the surface of the volume.
                # If the volume is closed, this integral should come out to the zero vector
                avgNormal = np.zeros(3)
                volSurfArea = 0.0

                # Make sure bdf is cross-referenced
                if self.bdfInfo.is_xrefed is False:
                    self.bdfInfo.cross_reference()
                    self.bdfInfo.is_xrefed = True

                elemIDs = self.meshLoader.getGlobalElementIDsForComps(
                    compIDs, nastranOrdering=True
                )
                for elemID in elemIDs:
                    elemInfo = self.bdfInfo.elements[elemID]
                    compID = self.meshLoader.nastranToTACSCompIDDict[elemInfo.pid]
                    if compID in properNormalCompIDs:
                        factor = 1.0
                    else:  # compID in flippedNormalCompIDs
                        factor = -1.0
                    n = elemInfo.Normal()
                    A = elemInfo.Area()
                    volSurfArea += A
                    avgNormal[:] += factor * n * A

                avgNormal /= volSurfArea
                if np.linalg.norm(avgNormal) > self.getOption("volCheckTol"):
                    self._TACSWarning("Specified volume may not be closed manifold.")

            # We can treat all components as being proper for solids,
            # since their volumes don't depend on their normals
            else:  # allSolids is True
                properNormalCompIDs = compIDs

        else:  # self.comm.rank != 0
            properNormalCompIDs = None
            flippedNormalCompIDs = None

        # Broadcast CompID lists from root
        properNormalCompIDs = self.comm.bcast(properNormalCompIDs, root=0)
        flippedNormalCompIDs = self.comm.bcast(flippedNormalCompIDs, root=0)

        return ParallelVolumeConstraint(
            self.assembler,
            self.meshLoader,
            properNormalCompIDs,
            flippedNormalCompIDs,
            self.Xpts,
            lb=lbound,
            ub=ubound,
        )

    def _sortIDs(self, currID, compIDs, orient1, orient2):
        """
        Sort all compIDs into one of two lists (orient1 or orient2), depending on
        the orientation of their neighbor components. This operation is performed recursively.
        """
        for adjID in self.adjacentOrientationMatch[currID]:
            sameOrientation = self.adjacentOrientationMatch[currID][adjID]
            if adjID not in compIDs:
                continue
            elif sameOrientation:
                if currID in orient1 and adjID not in orient1:
                    orient1.append(adjID)
                elif currID in orient2 and adjID not in orient2:
                    orient2.append(adjID)
                else:
                    continue
            else:  # flipped orientations
                if currID in orient1 and adjID not in orient2:
                    orient2.append(adjID)
                elif currID in orient2 and adjID not in orient1:
                    orient1.append(adjID)
                else:
                    continue
            self._sortIDs(adjID, compIDs, orient1, orient2)

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
        >>> volConstraint.evalConstraints(funcs, 'LE_SPAR')
        >>> funcs
        >>> # Result will look like (if VolumeConstraint has name of 'c1'):
        >>> # {'c1_wing': array([12354.10])}
        """
        # Check if user specified which constraints to output
        # Otherwise, output them all
        evalCons = self._processEvalCons(evalCons, ignoreMissing)

        # Loop through each requested constraint set
        for conName in evalCons:
            key = f"{self.name}_{conName}"
            funcs[key] = self.constraintList[conName].evalCon(self.Xpts)

        # Update call counter
        self.callCounter += 1

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
        >>> volConstraint.evalConstraintsSens(funcsSens, 'LE_SPAR')
        >>> funcsSens
        >>> # Result will look like (if VolumeConstraint has name of 'c1'):
        >>> # {'c1_wing':{'struct':<50x242 sparse matrix of type '<class 'numpy.float64'>' with 100 stored elements in Compressed Sparse Row format>}}
        """
        # Check if user specified which constraints to output
        # Otherwise, output them all
        evalCons = self._processEvalCons(evalCons)

        # Get number of dvs/coords on this proc
        nDVs = self.getNumDesignVars()
        ncoords = self.getNumCoordinates()

        # Loop through each requested constraint set
        for conName in evalCons:
            key = f"{self.name}_{conName}"
            # DV sensitivities are always zero for this constraint,
            funcsSens[key] = {}
            funcsSens[key][self.varName] = sp.sparse.csr_matrix(
                (1, nDVs), dtype=self.dtype
            )

            # Get nodal sensitivity
            xptSens = self.constraintList[conName].evalConSens(self.Xpts)
            xpt_array = xptSens.getArray()
            # Reshape from 1d array to a 1 by N column vector
            funcsSens[key][self.coordName] = xpt_array.reshape(1, ncoords)

    def writeVisualization(self, outputDir=None, baseName=None, number=None):
        """
        This function can be used to write a tecplot file for
        the purposes of visualization.

        Parameters
        ----------
        outputDir : str or None
            Use the supplied output directory
        baseName : str or None
            Use this supplied string for the base filename. Typically
            only used from an external solver.
        number : int or None
            Use the user supplied number to index output. Again, only
            typically used from an external solver
        """

        # Check input
        if outputDir is None:
            outputDir = self.getOption("outputDir")

        if baseName is None:
            baseName = self.name

        # If we are numbering output, it saving the sequence of
        # calls, add the call number
        if number is not None:
            # We need number based on the provided number:
            baseName = baseName + "_%3.3d" % number
        else:
            # if number is none, i.e. standalone, but we need to
            # number solutions, use internal counter
            if self.getOption("numberSolutions"):
                baseName = baseName + "_%3.3d" % self.callCounter

        base = os.path.join(outputDir, baseName) + ".dat"

        # Get current TACS nodes bvec
        localXpts = self.Xpts.getArray()
        localXpts = localXpts.reshape(-1, 3)

        # Setup tecplot writer on root
        if self.comm.rank == 0:
            tecInfo = tp.Tecplot(debug=None)

        # Loop through each constraint group and write to tecplot zone
        for constrName, constrObj in self.constraintList.items():
            # Identify all constrained compIDs
            posCompIDs = constrObj.posCompIDs
            negCompIDs = constrObj.negCompIDs
            allCompIDs = posCompIDs + negCompIDs
            # Get all global nodeIDs for this constraint
            allNodeIDsGlobal = self.meshLoader.getGlobalNodeIDsForComps(
                allCompIDs, nastranOrdering=True
            )
            # Find the local node ID corresponding to each global ID on each proc
            allNodeIDsLocal = self.meshLoader.getLocalNodeIDsFromGlobal(
                allNodeIDsGlobal, nastranOrdering=True
            )
            # Pull the updated xyz values of each node from each proc
            nodeXYZDict = {}
            for gID, lID in zip(allNodeIDsGlobal, allNodeIDsLocal):
                # Node is owned by this proc, add values to dict
                if lID >= 0:
                    nodeXYZDict[gID] = localXpts[lID]
            # Gather node xyz values to root
            nodeDicts = self.comm.gather(nodeXYZDict, root=0)

            # Write zone data for this constraint
            if self.comm.rank == 0:
                # Combine all node dicts received from each proc into a single dict
                allNodeXYZDict = {}
                for nodeXYZDict in nodeDicts:
                    allNodeXYZDict.update(nodeXYZDict)

                # Get all global elemIDs for this constraint
                allElemIDsGlobal = self.meshLoader.getGlobalElementIDsForComps(
                    allCompIDs, nastranOrdering=True
                )

                # Create new zone
                zoneInfo = Zone(tecInfo.log)

                # Setup header info
                zoneInfo.title = constrName
                zoneInfo.headers_dict["VARIABLES"] = ["X", "Y", "Z"]

                # Initialize connectivity/node locations
                nodeXYZ = []
                quadConn = []
                triConn = []
                hexConn = []
                tetConn = []

                # Loop through and set up connectivity/node locations
                currNodeID = 0
                oldToNewNodeID = {}
                for elemID in allElemIDsGlobal:
                    elemInfo = self.bdfInfo.elements[elemID]
                    compID = self.meshLoader.nastranToTACSCompIDDict[elemInfo.pid]
                    # Reverse the connectivity if normal is facing inward
                    if compID in negCompIDs:
                        elemInfo = copy.deepcopy(elemInfo)
                        elemInfo.flip_normal()
                    # We'll need to renumber the nodes to be sequential
                    oldConn = elemInfo.nodes
                    newConn = []
                    for oldNodeID in oldConn:
                        if oldNodeID not in oldToNewNodeID:
                            oldToNewNodeID[oldNodeID] = currNodeID
                            nodeXYZ.append(allNodeXYZDict[oldNodeID])
                            currNodeID += 1
                        newNodeID = oldToNewNodeID[oldNodeID]
                        newConn.append(newNodeID)

                    # Append new connectivity to the appropriate list
                    if "CQUAD" in elemInfo.type:
                        quadConn.append(newConn)
                    elif "CTRI" in elemInfo.type:
                        triConn.append(newConn)
                    elif "CHEXA" in elemInfo.type:
                        hexConn.append(newConn)
                    elif "CTETRA" in elemInfo.type:
                        tetConn.append(newConn)

                # Set all node/connectivity info to zone
                zoneInfo.xyz = np.array(nodeXYZ, dtype=float)
                zoneInfo.quad_elements = np.array(quadConn, dtype=np.intc)
                zoneInfo.tri_elements = np.array(triConn, dtype=np.intc)
                zoneInfo.hexa_elements = np.array(hexConn, dtype=np.intc)
                zoneInfo.tet_elements = np.array(tetConn, dtype=np.intc)

                # Add constraint zone to tecplot file
                tecInfo.zones.append(zoneInfo)

            # Wait for root
            self.comm.barrier()

        # Write out tecplot file
        if self.comm.rank == 0:
            tecInfo.write_tecplot(base)

        # Wait for root
        self.comm.barrier()


# Simple class for handling sparse volume constraints in parallel
class ParallelVolumeConstraint(object):
    dtype = TACSConstraint.dtype

    def __init__(
        self, assembler, meshLoader, posCompIDs, negCompIDs, Xpts0, lb=0.0, ub=1e20
    ):
        # Number of constraints
        self.nCon = 1
        self.assembler = assembler

        # Create two TACS Volume function domains,
        # One for the properly oriented components...
        if len(posCompIDs) > 0:
            posVolFunc = EnclosedVolume(self.assembler)
            elemIDs = meshLoader.getLocalElementIDsForComps(posCompIDs)
            posVolFunc.setDomain(elemIDs)
        else:
            posVolFunc = None

        # and one for the flipped components.
        # The volumes computed by these components will be negative
        if len(negCompIDs) > 0:
            negVolFunc = EnclosedVolume(self.assembler)
            elemIDs = meshLoader.getLocalElementIDsForComps(negCompIDs)
            negVolFunc.setDomain(elemIDs)
        else:
            negVolFunc = None

        # Compute the volume and check that we sorted the correct orientations
        self.posCompIDs = posCompIDs
        self.negCompIDs = negCompIDs
        self.posVolFunc = posVolFunc
        self.negVolFunc = negVolFunc
        vol0 = self.evalCon(Xpts0)

        # If the volume is negative, we need to flip the orientations
        if vol0 < 0.0:
            self.posVolFunc, self.negVolFunc = self.negVolFunc, self.posVolFunc
            self.posCompIDs, self.negCompIDs = self.negCompIDs, self.posCompIDs

        # Save bound information
        if isinstance(lb, np.ndarray) and len(lb) == self.nCon:
            self.lb = lb.astype(self.dtype)
        elif isinstance(lb, float) or isinstance(lb, complex):
            self.lb = np.array([lb] * self.nCon, dtype=self.dtype)

        if isinstance(ub, np.ndarray) and len(ub) == self.nCon:
            self.ub = ub.astype(self.dtype)
        elif isinstance(ub, float) or isinstance(ub, complex):
            self.ub = np.array([ub] * self.nCon, dtype=self.dtype)

    def evalCon(self, x):
        self.assembler.setNodes(x)
        posVol, negVol = self.assembler.evalFunctions(
            [self.posVolFunc, self.negVolFunc]
        )
        totVol = posVol - negVol
        return totVol

    def evalConSens(self, x):
        self.assembler.setNodes(x)
        totSens = self.assembler.createNodeVec()
        # Compute positive and negative volumes contributions to sens
        self.assembler.addXptSens([self.posVolFunc], [totSens], 1.0)
        self.assembler.addXptSens([self.negVolFunc], [totSens], -1.0)
        # Distribute non-local vals across procs
        totSens.beginSetValues()
        totSens.endSetValues()
        return totSens

    def getBounds(self):
        return self.lb.copy(), self.ub.copy()


# Slight modification of pyNASTRAN's tecplot Zone class
# By default this class doesn't write out the zone title in the tecplot file for some reason
# We make a new class here with an overriden method that includes this in the zone header
class Zone(tp.Zone):
    def write_unstructured_zone(
        self,
        tecplot_file,
        ivars,
        is_points,
        nnodes,
        nelements,
        zone_type,
        log,
        is_tris,
        is_quads,
        is_tets,
        is_hexas,
        adjust_nids=True,
    ):
        msg = "ZONE "
        self.log.info("is_points = %s" % is_points)
        datapacking = "POINT" if is_points else "BLOCK"
        # Make sure to include title
        msg += f' T="{self.title}", n={nnodes:d}, e={nelements:d}, ZONETYPE={zone_type}, DATAPACKING={datapacking}\n'
        tecplot_file.write(msg)

        self._write_xyz_results(tecplot_file, is_points, ivars)
        self._write_elements(
            tecplot_file,
            nnodes,
            is_tris,
            is_quads,
            is_tets,
            is_hexas,
            adjust_nids=adjust_nids,
        )
