#!/usr/bin/python
"""
pymeshloader - A python-based mesh loader for setting up TACS

This python interface is designed to provide a easier interface to the
c-layer of TACS. This module uses the pyNastran library for its bdf parsing
functionality.

"""
# =============================================================================
# Imports
# =============================================================================
import numpy as np
import tacs.TACS, tacs.elements
from pyNastran.bdf.bdf import read_bdf

class pyMeshLoader(object):

    def __init__(self, comm, dtype, printDebug=False):
        # MPI communicator
        self.comm = comm
        # TACS scalara data type (float or complex)
        self.dtype = dtype
        # Debug printing flag
        self.printDebug = printDebug
        self.bdfInfo = None

    def scanBdfFile(self, fileName):
        """
        Scan nastran bdf file using pyNastran's bdf parser.
        We also set up arrays that will be require later to build tacs.
        """

        # Only print debug info on root, if requested
        if self.comm.rank == 0:
            debugPrint = self.printDebug
        else:
            debugPrint = False

        # Read in bdf file as pynsatran object
        # By default we avoid cross-referencing unless we actually need it,
        # since its expensive for large models
        self.bdfInfo = read_bdf(fileName, validate=False, xref=False, debug=debugPrint)
        # Set flag letting us know model is not xrefed yet
        self.bdfInfo.is_xrefed = False

        # If any property cards are msiing in the bdf, we have to add dummy cards
        # so pynastran doesn't through errors when cross-referencing
        # Loop through all elements and add dummy property, as necessary
        self.bdfInfo.missing_properties = False
        for element_id in self.bdfInfo.elements:
            element = self.bdfInfo.elements[element_id]
            if element.pid not in self.bdfInfo.property_ids:
                # If no material properties were found,
                # add dummy properties and materials
                matID = 1
                E = 70.0
                G = 35.0
                nu = 0.3
                self.bdfInfo.add_mat1(matID, E, G, nu)
                self.bdfInfo.add_pbar(element.pid, matID)
                # Warn the user that the property card is missing
                # and should not be read in using pytacs elemCallBackFromBDF method
                self.bdfInfo.missing_properties = True
                if self.printDebug:
                    TACSWarning('Element ID %d references undefined property ID %d in bdf file. '
                                'A user-defined elemCalBack function will need to be provided.'% (element_id, element.pid),
                                self.comm)

        # Create dictionaries for mapping between tacs and nastran id numbering
        self._updateNastranToTACSDicts()

        # Try to get the node x,y,z locations from bdf file
        try:
            self.bdfXpts = self.bdfInfo.get_xyz_in_coord()
        # If this fails, the file may reference multiple coordinate systems
        # and will have to be cross-refernced to work
        except:
            self.bdfInfo.cross_reference()
            self.bdfInfo.is_xrefed = True
            self.bdfXpts = self.bdfInfo.get_xyz_in_coord()

        # element card contained within each property group (may contain multiple per group)
        # Each entry will eventually have its own tacs element object assigned to it
        # (ex. self.elemDescripts = [['CQUAD4', 'CTRIA3'], ['CQUAD9', 'CQUAD4'], ['CBAR']])
        self.elemDescripts = []
        # Pointer index used to map elemDescript entries into flat list
        # (ex. self.elementObjectNumByComp = [[0, 1], [2, 3], [4]])
        self.elemObjectNumByComp = []
        # User-defined name for property/component group,
        # this may be read in through ICEMCFD/FEMAP format comments in BDF
        self.compDescripts = []

        # Populate list entries with default values
        for pID in self.bdfInfo.property_ids:
            self.elemDescripts.append([])
            self.elemObjectNumByComp.append([])
            # Check if there is a Femap label for this component
            propComment = self.bdfInfo.properties[pID].comment
            if '$ Femap Property' in propComment:
                # Pick off last word from comment, this is the name
                propName = propComment.split()[-1]
                self.compDescripts.append(propName)
            #  Default component name
            else:
                self.compDescripts.append('untitled')

        # Element connectivity information
        self.elemConnectivity = [None] * self.bdfInfo.nelements
        self.elemConnectivityPointer = [None] * (self.bdfInfo.nelements + 1)
        self.elemConnectivityPointer[0] = 0
        elementObjectCounter = 0
        # List specifying which tacs element object each element in bdf should point to
        self.elemObjectNumByElem = [None] * (self.bdfInfo.nelements)

        # Loop through every element and record information needed for tacs
        for tacsElementID, nastranElementID in enumerate(self.bdfInfo.element_ids):
            element = self.bdfInfo.elements[nastranElementID]
            elementType = element.type.upper()
            propertyID = element.pid
            componentID = self.idMap(propertyID, self.nastranToTACSCompIDDict)

            # This element type has not been added to the list for the component group yet, so we append it
            if elementType not in self.elemDescripts[componentID]:
                self.elemDescripts[componentID].append(elementType)
                self.elemObjectNumByComp[componentID].append(elementObjectCounter)
                elementObjectCounter += 1

            # Find the index number corresponding to the element object number for this component
            componentTypeIndex = self.elemDescripts[componentID].index(elementType)
            self.elemObjectNumByElem[tacsElementID] = self.elemObjectNumByComp[componentID][componentTypeIndex]

            # We've identified a ICEM property label
            if 'Shell element data for family' in element.comment:
                componentName = element.comment.split()[-1]
                self.compDescripts[componentID] = componentName

            conn = element.nodes.copy()

            # TACS has a different node ordering than Nastran for certain elements,
            # we now perform the reordering (if necessary)
            if elementType in ['CQUAD4', 'CQUADR']:
                conn = [conn[0], conn[1], conn[3], conn[2]]
            elif elementType in ['CQUAD9', 'CQUAD']:
                conn = [conn[0], conn[4], conn[1], conn[7], conn[8], conn[5], conn[3], conn[6], conn[2]]
            elif elementType in ['CHEXA8', 'CHEXA']:
                conn = [conn[0], conn[1], conn[3], conn[2], conn[4], conn[5], conn[7], conn[6]]

            # Map node ids in connectivity from Nastan numbering to TACS numbering
            self.elemConnectivity[tacsElementID] = self.idMap(conn, self.nastranToTACSNodeIDDict)
            self.elemConnectivityPointer[tacsElementID + 1] = self.elemConnectivityPointer[
                                                                  tacsElementID] + len(element.nodes)

        # Allocate list for user-specified tacs element objects
        self.elemObjects = [None] * elementObjectCounter

    def _updateNastranToTACSDicts(self):
        '''
        Create dictionaries responsible for mapping over
        global node, global element, and property/component ID numbers
        from NASTRAN (as read in from BDF) to TACS (contiguous, 0-based).
        The keys of each dictionary are the NASTRAN ID and the entry the TACS ID,
        such that:
            tacsNodeID = self.nastranToTACSNodeIDDict[nastranNodeID]
            tacsComponentID = self.nastranToTACSCompIDDict[nastranPropertyID]
            tacsElementID = self.nastranToTACSNodeIDDict[nastranElementID]
        The dictionaries contain all elements/nodes found in the BDF,
        not just those *owned* by this processor
        '''
        # Create Node ID map
        nastranIDs = self.bdfInfo.node_ids
        tacsIDs = range(self.bdfInfo.nnodes)
        nodeTuple = zip(nastranIDs, tacsIDs)
        self.nastranToTACSNodeIDDict = dict(nodeTuple)

        # Create Property/Component ID map
        nastranIDs = self.bdfInfo.property_ids
        tacsIDs = range(self.bdfInfo.nproperties)
        propTuple = zip(nastranIDs, tacsIDs)
        self.nastranToTACSCompIDDict = dict(propTuple)

        # Create Element ID map
        nastranIDs = self.bdfInfo.element_ids
        tacsIDs = range(self.bdfInfo.nelements)
        elemTuple = zip(nastranIDs, tacsIDs)
        self.nastranToTACSElemIDDict = dict(elemTuple)

    def getBDFInfo(self):
        """
        Return pynastran bdf object.
        """
        return self.bdfInfo

    def getNumComponents(self):
        '''
        Return number of nodes found in bdf.
        '''
        return self.bdfInfo.nproperties

    def getNumComponents(self):
        '''
        Return number of nodes found in bdf.
        '''
        return self.bdfInfo.nproperties

    def getNumBDFNodes(self):
        '''
        Return number of nodes found in bdf.
        '''
        return self.bdfInfo.nnodes

    def getNumOwnedNodes(self):
        '''
        Return number of nodes owned by this processor.
        '''
        return self.assembler.getNumOwnedNodes()

    def getNumBDFElements(self):
        '''
        Return number of elements found in bdf.
        '''
        return self.bdfInfo.nelements

    def getBDFNodes(self, nodeIDs, nastranOrdering=False):
        '''
        Return x,y,z location of specified node in bdf file.
        '''
        # Convert to tacs numbering, if necessary
        if nastranOrdering:
            nodeIDs = self.idMap(nodeIDs, self.nastranToTACSNodeIDDict)
        return self.bdfXpts[nodeIDs]

    def getElementComponents(self):
        elements = self.bdfInfo.elements
        propertyIDList = [elements[eID].pid for eID in self.bdfInfo.element_ids]
        compIDList = self.idMap(propertyIDList, self.nastranToTACSCompIDDict)
        return compIDList

    def getConnectivityForComp(self, componentID, nastranOrdering=False):
        # Find all of the element IDs belonging to this property group
        propertyID = list(self.bdfInfo.property_ids)[componentID]
        elementIDs = self.bdfInfo.get_element_ids_dict_with_pids()[propertyID]
        compConn = []
        for elementID in elementIDs:
            # We've now got the connectivity for this element, but it is in nastrans node numbering
            nastranConn = self.bdfInfo.elements[elementID].nodes
            if nastranOrdering:
                compConn.append(nastranConn)
            else:
                # Convert Nastran node numbering back to tacs numbering
                tacsConn = [self.bdfInfo.nid_map[nID] for nID in nastranConn]
                # Append element connectivity to list for component
                compConn.append(tacsConn)
        return compConn

    def getElementObjectNumsForComp(self, componentID):
        """
        Return tacs element object number for each element type
        belonging to this component.
        """
        return self.elemObjectNumByComp[componentID][:]

    def getElementDescripts(self):
        """
        Get nested list containing all element types owned by each component group
        example: [['CQUAD4', 'CTRIA3], ['CQUAD4'], ['CQUAD4', CQUAD9']]
        """
        return self.elemDescripts

    def getComponentDescripts(self):
        return self.compDescripts

    def getLocalNodeIDsFromGlobal(self, globalIDs, nastranOrdering=False):
        """
        given a list of node IDs in global (non-partitioned) ordering
        returns the local (partitioned) node IDs on each processor.
        If a requested node is not included on this processor,
        an entry of -1 will be returned.
        """
        # Convert to tacs node numbering, if necessary
        if nastranOrdering:
            globalIDs = self.idMap(globalIDs, self.nastranToTACSNodeIDDict)

        # Ensure input is list-like
        globalIDs = np.atleast_1d(globalIDs)

        # Get the node id offset for this processor
        OwnerRange = self.assembler.getOwnerRange()
        nodeOffset = OwnerRange[self.comm.rank]

        # Get the local ID numbers for this proc
        tacsLocalIDs = []
        for gID in globalIDs:
            lIDs = self.creator.getAssemblerNodeNums(self.assembler,
                                                     np.array([gID], dtype=np.intc))
            # Node was not found on this proc, return -1
            if len(lIDs) == 0:
                tacsLocalIDs.append(-1)
            # Node was found on this proc, shift by nodeOffset to get local index for node
            else:
                tacsLocalIDs.append(lIDs[0] - nodeOffset)

        return tacsLocalIDs

    def getLocalElementIDsFromGlobal(self, globalIDs, nastranOrdering=False):
        """
        given a list of element IDs in global (non-partitioned) ordering
        returns the local (partitioned) element IDs on each processor.
        If a requested element is not included on this processor,
        an entry of -1 will be returned.
        """
        # Convert to tacs node numbering, if necessary
        if nastranOrdering:
            globalIDs = self.idMap(globalIDs, self.nastranToTACSElemIDDict)

        # Ensure input is list-like
        globalIDs = np.atleast_1d(globalIDs)

        # Get the node id offset for this processor
        OwnerRange = self.assembler.getOwnerRange()
        nodeOffset = OwnerRange[self.comm.rank]

        # Get the local ID numbers for this proc
        tacsLocalIDs = []
        for gID in globalIDs:
            # element was found on this proc, get local ID num
            if gID in self.globalToLocalNodeIDDict:
                lID = self.globalToLocalNodeIDDict[gID]
            # element was not found on this proc, return -1
            else:
                lID = -1
            tacsLocalIDs.append(lID)

        return tacsLocalIDs

    def getGlobalElementIDsForComps(self, componentIDs, nastranOrdering=False):
        """
        Returns a list of element IDs belonging to specified components
        """
        # Make sure list is flat
        componentIDs = self._flatten(componentIDs)
        # Convert tacs component IDs to nastran property IDs
        propertyIDs = [0] * len(componentIDs)
        for i, componentID in enumerate(componentIDs):
            propertyIDs[i] = list(self.bdfInfo.property_ids)[componentID]
        # Get dictionary whose values are the element ids we are looking for
        elementIDDict = self.bdfInfo.get_element_ids_dict_with_pids(propertyIDs)
        # Convert to list
        elementIDs = list(elementIDDict.values())
        # Make sure list is flat
        elementIDs = self._flatten(elementIDs)
        # Convert to tacs element numbering, if necessary
        if not nastranOrdering:
            elementIDs = self.idMap(elementIDs, self.nastranToTACSElemIDDict)
        return elementIDs

    def getLocalElementIDsForComps(self, componentIDs):
        """
        Get the local element numbers on each proc used by tacs
        corresponding to the component groups in componentIDs.
        """
        if self.creator is None:
            raise Error("TACS assembler has not been created. "
                        "Assembler must created first by running 'createTACS' method.")
        # Make sure list is flat
        componentIDs = self._flatten(componentIDs)
        # Get the element object IDs belonging to each of these components
        objIDs = []
        for componentID in componentIDs:
            tmp = self.getElementObjectNumsForComp(componentID)
            objIDs.extend(tmp)
        objIDs = np.array(objIDs, dtype=np.intc)
        # Get the local element IDs corresponding to the object IDs on this processor (if any)
        localElemIDs = self.creator.getElementIdNums(objIDs)
        return list(localElemIDs)

    def getGlobalToLocalNodeIDDict(self):
        """
        Creates a dictionary who's keys correspond to the global ID of each node (tacs ordering)
        owned by this processor and whose values correspond to the local node index for this proc.
        The dictionary can be used to convert from a global node ID to local using the assignment below*:
        localNodeID = globalToLocalNodeIDDict[globalNodeID]

        * assuming globalNodeID is owned on this processor
        """
        # Do sorting on root proc
        if self.comm.rank == 0:
            # List containing which poc every element belongs to
            elemPartition = self.creator.getElementPartition()
            # Create an empty list that we will use to sort what elements are on what procs
            allOwnedElementIDs = [[] for proc_i in range(self.comm.size)]
            for elemGlobalID in range(len(elemPartition)):
                ownerProc = elemPartition[elemGlobalID]
                allOwnedElementIDs[ownerProc].append(elemGlobalID)
        else:
            allOwnedElementIDs = None

        # Scatter the list from the root so each proc knows what element ID it owns
        ownedElementIDs = self.comm.scatter(allOwnedElementIDs, root=0)
        # Create dictionary that gives the corresponding local ID for each global ID owned by this proc
        globalToLocalNodeIDDict = {gID: lID for lID, gID in enumerate(ownedElementIDs)}

        return globalToLocalNodeIDDict

    def setElementObject(self, componentID, objectIndex, elemObject):
        pointer = self.elemObjectNumByComp[componentID][objectIndex]
        self.elemObjects[pointer] = elemObject

    def getElementObjectForElemID(self, elemID, nastranOrdering=False):
        # Convert to tacs numbering, if necessary
        if nastranOrdering:
            elemID = self.idMap(elemID, self.nastranToTACSElemIDDict)
        # Get the pointer for the tacs element object for this element
        elemObjNum = self.elemObjectNumByElem[elemID]
        elemObj = self.elemObjects[elemObjNum]
        return elemObj

    def createTACSAssembler(self, varsPerNode):
        """
        Setup TACSCreator object responsible for creating TACSAssembler
        """
        self.creator = tacs.TACS.Creator(self.comm, varsPerNode)

        # Append RBE elements to element list, these are not setup by the user
        for rbe in self.bdfInfo.rigid_elements.values():
            if rbe.type == 'RBE2':
                self._addTACSRBE2(rbe, varsPerNode)
            elif rbe.type == 'RBE3':
                self._addTACSRBE3(rbe, varsPerNode)
            else:
                raise NotImplementedError("Rigid element of type '{}' is not supported".format(rbe.type))

        # Setup element connectivity and boundary condition info on root processor
        if self.comm.rank == 0:
            # Set connectivity for all elements
            ptr = np.array(self.elemConnectivityPointer, dtype=np.intc)
            conn = np.array(self._flatten(self.elemConnectivity), dtype=np.intc)
            objectNums = np.array(self.elemObjectNumByElem, dtype=np.intc)
            self.creator.setGlobalConnectivity(self.bdfInfo.nnodes, ptr, conn, objectNums)

            # Set up the boundary conditions
            bcDict = {}
            for spc_id in self.bdfInfo.spcs:
                for spc in self.bdfInfo.spcs[spc_id]:
                    # Loop through every node specifed in this spc and record bc info
                    for j, nastranNode in enumerate(spc.nodes):
                        tacsNode = self.idMap(nastranNode, self.nastranToTACSNodeIDDict)
                        # If node hasn't been added to bc dict yet, add it
                        if tacsNode not in bcDict:
                            bcDict[tacsNode] = {}

                        # Loop through each dof and record bc info if it is included in this spc
                        for dof in range(varsPerNode):
                            # Add 1 to get nastran dof number
                            nastranDOF = dof + 1
                            if spc.type == 'SPC':
                                # each node may have its own dofs uniquely constrained
                                constrainedDOFs = spc.components[j]
                                # The boundary condition may be forced to a non-zero value
                                constrainedVal = spc.enforced[j]
                            else:  # SPC1?
                                # All nodes always have the same dofs constrained
                                constrainedDOFs = spc.components
                                # This boundary condition is always 0
                                constrainedVal = 0.0
                            # if nastran dof is in spc components string, add it to the bc dict
                            if self._isDOFInString(constrainedDOFs, nastranDOF):
                                bcDict[tacsNode][dof] = constrainedVal


            # Convert bc information from dict to list
            bcnodes = []
            bcdofs = []
            bcptr = [0]
            bcvals = []
            numbcs = 0
            for tacsNode in bcDict:
                bcnodes.append(tacsNode)
                # Store constrained dofs for this node
                dofs = bcDict[tacsNode].keys()
                bcdofs.extend(dofs)
                # Store enforced bc value
                vals = bcDict[tacsNode].values()
                bcvals.extend(vals)
                # Increment bc pointer with how many constraints weve added for this node
                numbcs += len(bcDict[tacsNode])
                bcptr.append(numbcs)

            # Recast lists as numpy arrays
            bcnodes = np.array(bcnodes, dtype=np.intc)
            bcdofs = np.array(bcdofs, dtype=np.intc)
            bcptr = np.array(bcptr, dtype=np.intc)
            bcvals = np.array(bcvals, dtype=self.dtype)
            # Set boundary conditions in tacs
            self.creator.setBoundaryConditions(bcnodes, bcptr, bcdofs, bcvals)

            # Set node locations
            Xpts = self.bdfInfo.get_xyz_in_coord().astype(tacs.TACS.dtype)
            self.creator.setNodes(Xpts.flatten())

        # Set the elements for each component
        self.creator.setElements(self.elemObjects)

        self.assembler = self.creator.createTACS()

        self.globalToLocalNodeIDDict = self.getGlobalToLocalNodeIDDict()

        return self.assembler

    def _isDOFInString(self, constrained_dofs, dof):
        """
        Find if dof number (nastran numbering) occurs in constraint string.

        Parameters
        ----------
        constrained_dofs : string
            String containing list of dofs (ex. '123456')
        dof : int or string
            nastran dof number to check for
        """
        # Convert to string, if necessary
        if isinstance(dof, int):
            dof = '%d' % (dof)
        # pyNastran only supports 0,1,2,3,4,5,6 as valid dof components
        # For this reason, we'll treat 0 as if its 7, since it's traditionally never used in nastran
        if dof == '7':
            dof = '0'
        location = constrained_dofs.find(dof)
        # if dof is found, return true
        if location > -1:
            return True
        else:
            return False

    def _addTACSRBE2(self, rbeInfo, varsPerNode):
        """
        Method to automatically set up RBE2 element from bdf file for user.
        User should *NOT* set these up in their elemCallBack function.
        """
        indepNode = rbeInfo.independent_nodes
        depNodes = []
        depConstrainedDOFs = []
        dummyNodes = []
        dofsAsString = rbeInfo.cm
        dofsAsList = self.isDOFInString(dofsAsString, varsPerNode)
        for node in rbeInfo.dependent_nodes:
            depNodes.append(node)
            depConstrainedDOFs.extend(dofsAsList)
            # add dummy nodes for all lagrange multiplier
            dummyNodeNum = list(self.bdfInfo.node_ids)[-1] + 1  # Next available node number
            # Add the dummy node coincident to the dependent node in x,y,z
            self.bdfInfo.add_grid(dummyNodeNum, self.bdfInfo.nodes[node].xyz)
            dummyNodes.append(dummyNodeNum)

        conn = indepNode + depNodes + dummyNodes
        nTotalNodes = len(conn)
        # Update Nastran to TACS ID mapping dicts, since we just added new nodes to model
        self._updateNastranToTACSDicts()
        # Append RBE information to the end of the element lists
        self.elemConnectivity.append(self.idMap(conn, self.nastranToTACSNodeIDDict))
        self.elemConnectivityPointer.append(self.elemConnectivityPointer[-1] + nTotalNodes)
        rbeObj = tacs.elements.RBE2(nTotalNodes, np.array(depConstrainedDOFs, dtype=np.intc))
        self.elemObjectNumByElem.append(self.elemObjectNumByElem[-1] + 1)
        self.elemObjects.append(rbeObj)
        return

    def _addTACSRBE3(self, rbeInfo, varsPerNode):
        """
        Method to automatically set up RBE3 element from bdf file for user.
        User should *NOT* set these up in their elemCallBack function.
        """
        depNode = rbeInfo.dependent_nodes
        depConstrainedDOFs = self.isDOFInString(rbeInfo.refc, varsPerNode)

        # add dummy node for lagrange multipliers
        dummyNodeNum = list(self.bdfInfo.node_ids)[-1] + 1  # Next available node number
        # Add the dummy node coincident to the dependent node in x,y,z
        self.bdfInfo.add_grid(dummyNodeNum, self.bdfInfo.nodes[depNode[0]].xyz)
        dummyNodes = [dummyNodeNum]
        # Update Nastran to TACS ID mapping dicts, since we just added new nodes to model
        self._updateNastranToTACSDicts()

        indepNodes = []
        indepWeights = []
        indepConstrainedDOFs = []
        for depNodeGroup in rbeInfo.wt_cg_groups:
            wt = depNodeGroup[0]
            dofsAsString = depNodeGroup[1]
            dofsAsList = self.isDOFInString(dofsAsString, varsPerNode)
            for node in depNodeGroup[2]:
                indepNodes.append(node)
                indepWeights.append(wt)
                indepConstrainedDOFs.extend(dofsAsList)

        conn = depNode + indepNodes + dummyNodes
        nTotalNodes = len(conn)
        # Append RBE information to the end of the element lists
        self.elemConnectivity.append(self.idMap(conn, self.nastranToTACSNodeIDDict))
        self.elemConnectivityPointer.append(self.elemConnectivityPointer[-1] + nTotalNodes)
        rbeObj = tacs.elements.RBE3(nTotalNodes,
                                    np.array(depConstrainedDOFs, dtype=np.intc),
                                    np.array(indepWeights),
                                    np.array(indepConstrainedDOFs, dtype=np.intc))
        self.elemObjects.append(rbeObj)
        self.elemObjectNumByElem.append(self.elemObjectNumByElem[-1] + 1)
        return

    def isDOFInString(self, dofString, numDOFs):
        """
        Converts a dof string to a boolean list.
        Examples:
            '123' -> [1, 1, 1, 0, 0, 0]
            '1346' -> [1, 0, 1, 1, 0, 1]
        """
        dofList = []
        for dof in range(numDOFs):
            dof = str(dof + 1)
            loc = dofString.find(dof)
            if loc > -1:
                dofList.append(1)
            else:
                dofList.append(0)
        return dofList

    def idMap(self, fromIDs, tacsIDDict):
        """
        Translate fromIDs numbering from nastran numbering to tacs numbering.
        If node ID doesn't exist in nastranIDList, return -1 for entry
        """
        # Input is a list return a list
        if hasattr(fromIDs, '__iter__'):
            toIDs = [None] * len(fromIDs)
            # Iterate through list and call function recursively one element at a time
            for i, id in enumerate(fromIDs):
                toIDs[i] = self.idMap(id, tacsIDDict)
            return toIDs
        # Input is a int, return an int
        else:
            if fromIDs in tacsIDDict:
                return tacsIDDict[fromIDs]
            else:
                return -1

    def _flatten(self, l, ltypes=(list, tuple)):
        """
        Flattens nested list into single 1d list
        """
        ltype = type(l)
        l = list(l)
        i = 0
        while i < len(l):
            while isinstance(l[i], ltypes):
                if not l[i]:
                    l.pop(i)
                    i -= 1
                    break
                else:
                    l[i:i + 1] = l[i]
            i += 1
        return ltype(l)

class Error(Exception):
    """
    Format the error message in a box to make it clear this
    was a expliclty raised exception.
    """

    def __init__(self, message):
        msg = '\n+' + '-' * 78 + '+' + '\n' + '| pyMeshLoader Error: '
        i = 21
        for word in message.split():
            if len(word) + i + 1 > 78:  # Finish line and start new one
                msg += ' ' * (78 - i) + '|\n| ' + word + ' '
                i = 1 + len(word) + 1
            else:
                msg += word + ' '
                i += len(word) + 1
        msg += ' ' * (78 - i) + '|\n' + '+' + '-' * 78 + '+' + '\n'
        print(msg)
        Exception.__init__(self)

class TACSWarning(object):
    """
    Format a warning message
    """

    def __init__(self, message, comm):
        if comm.rank == 0:
            msg = '\n+' + '-' * 78 + '+' + '\n' + '| pyMeshLoader Warning: '
            i = 23
            for word in message.split():
                if len(word) + i + 1 > 78:  # Finish line and start new one
                    msg += ' ' * (78 - i) + '|\n| ' + word + ' '
                    i = 1 + len(word) + 1
                else:
                    msg += word + ' '
                    i += len(word) + 1
            msg += ' ' * (78 - i) + '|\n' + '+' + '-' * 78 + '+' + '\n'
            print(msg)