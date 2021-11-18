"""
pyBase_problem
"""

# =============================================================================
# Imports
# =============================================================================
import warnings
import numpy as np
from mpi4py import MPI
from ..utilities import BaseUI
from collections import OrderedDict
import tacs.TACS, tacs.constitutive, tacs.elements, tacs.functions, tacs.problems.static

class BaseProblem(BaseUI):
    """
    Base class for TACS problem types. Contains methods common to all TACS problems.
    """
    def __init__(self, assembler, comm, outputViewer=None, meshLoader=None):

        # TACS assembler object
        self.assembler = assembler
        # TACS F5 output writer
        self.outputViewer = outputViewer
        # TACS pyMeshLoader object
        self.meshLoader = meshLoader
        # MPI communicator object
        self.comm = comm

        # Data type to use for TACS scalars (float or complex)
        self.dtype = tacs.TACS.dtype

        # Create Design variable vector
        self.x = self.assembler.createDesignVec()
        self.assembler.getDesignVars(self.x)
        self.varName = 'struct'
        # Create Nodal coordinate vector
        self.Xpts = self.assembler.createNodeVec()
        self.assembler.getNodes(self.Xpts)
        self.coordName = 'Xpts'
        # List of functions
        self.functionList = OrderedDict()

        return

    ####### Design variable methods ########

    def setVarName(self, varName):
        """
        Set a name for the structural variables in pyOpt. Only needs
        to be changed if more than 1 pytacs object is used in an
        optimization

        Parameters
        ----------
        varName : str
            Name of the structural variable used in addVarGroup().
            """
        self.varName = varName

    def getDesignVars(self):
        """
        get the design variables that were specified with
        addVariablesPyOpt.

        Returns
        ----------
        x : array
            The current design variable vector set in tacs.

        Notes
        -----
        This routine **can** also accept a list or vector of
        variables. This is used internally in pytacs, but is not
        recommended to used externally.
        """
        return self.x.getArray().copy()

    def setDesignVars(self, x):
        """
        Update the design variables used by tacs.

        Parameters
        ----------
        x : ndarray
            The variables (typically from the optimizer) to set. It
            looks for variable in the ``self.varName`` attribute.

        """
        # Check if the design variables are being handed in a dict
        if isinstance(x, dict):
            if self.varName in x:
                self.x.getArray()[:] = x[self.varName]
        # or array
        elif isinstance(x, np.ndarray):
            self.x.getArray()[:] = x
        # Or TACS BVec
        elif isinstance(x, tacs.TACS.Vec):
            self.x.copyValues(x)
        else:
            raise ValueError("setDesignVars must be called with either a numpy array, dict, or TACS Vec as input.")

        # Set the variables in tacs
        self.assembler.setDesignVars(self.x)

    def _arrayToDesignVec(self, dvArray):
        """
        Converts a distributed numpy array into a TACS design variable BVec.
        NOTE: dvArray must have correct size on each processor
        """
        xVec = self.assembler.createDesignVec()

        # Set values
        xVec.getArray()[:] = dvArray

        # Return as tacs bvec object
        return xVec

    def getNumDesignVars(self):
        """
        Return the number of design variables on this processor.
        """
        return self.x.getSize()

    # TODO: Change below to getNodes/setNodes for consistency
    def getCoordinates(self):
        """
        Return the mesh coordiantes of this problem.

        Returns
        -------
        coords : array
            Structural coordinate in array of size (N, 3) where N is
            the number of structural nodes on this processor.
        """
        return self.Xpts.getArray().copy()

    def setCoordinates(self, Xpts):
        """
        Set the mesh coordinates of this problem.

        Returns
        -------
        coords : array
            Structural coordinate in array of size (N, 3) where N is
            the number of structural nodes on this processor.
        """
        # Check if the design variables are being handed in a dict
        if isinstance(Xpts, dict):
            if self.coordName in Xpts:
                self.Xpts.getArray()[:] = Xpts[self.coordName]
        # or array
        elif isinstance(Xpts, np.ndarray):
            self.Xpts.getArray()[:] = Xpts
        # Or TACS BVec
        elif isinstance(Xpts, tacs.TACS.Vec):
            self.Xpts.copyValues(Xpts)
        else:
            raise ValueError("setCoordinates must be called with either a numpy array, dict, or TACS Vec as input.")
        self.assembler.setNodes(self.Xpts)

    def _arrayToNodeVec(self, xptsArray):
        """
        Converts a distributed numpy array into a TACS nodal BVec.
        NOTE: xptsArray must have correct size on each processor
        """
        Xptsvec = self.assembler.createNodeVec()

        # Set values
        Xptsvec.getArray()[:] = xptsArray

        # Return as tacs bvec object
        return Xptsvec

    def getNumCoordinates(self):
        """
        Return the number of mesh coordinates on this processor.
        """
        return self.Xpts.getSize()

    ####### Variable methods ########

    def getVarsPerNode(self):
        """
        Get the number of variables per node for the model.
        """
        return self.assembler.getVarsPerNode()

    def getNumOwnedNodes(self):
        """
        Get the number of nodes owned by this processor.
        """
        return self.assembler.getNumOwnedNodes()

    def _arrayToVec(self, varArray):
        """
        Converts a distributed numpy array into a TACS state variable BVec.
        NOTE: dvArray must have correct size on each processor
        """
        varVec = self.assembler.createVec()

        # Set values
        varVec.getArray()[:] = varArray

        # Return as tacs bvec object
        return varVec

    def getNumVariables(self):
        """Return the number of degrees of freedom (states) that are
        on this processor

        Returns
        -------
        nstate : int
            number of states.
        """
        vpn =  self.getVarsPerNode()
        nnodes =  self.getNumOwnedNodes()
        return vpn * nnodes

    ####### Eval function methods ########

    def addFunction(self, funcName, funcHandle, compIDs=None, **kwargs):
        """
        Generic function to add a function for TACS. It is intended to
        be reasonably generic since the user supplies the actual
        function handle to use. The following functions can be used:
        KSFailure, KSBuckling, MaxBuckling, AverageKSFailure,
        MaxFailure, AverageMaxFailure, AverageKSBuckling,
        StructuralMass, Compliance, AggregateDisplacement.

        Parameters
        ----------
        funcName : str
            The user-supplied name for the function. This will
            typically be a string that is meanful to the user
        funcHandle : tacs.functions
            The fucntion handle to use for creation. This must come
            from the functions module in tacs.
        compIDs: list
            List of compIDs to select. Alternative to selectCompIDs
            arguments.
        """

        # We try to setup the function, if it fails it may not be implimented:
        try:
            # pass assembler an function-specific kwargs straight to tacs function
            self.functionList[funcName] = funcHandle(self.assembler, **kwargs)
        except:
            self.TACSWarning("Function type %s is not currently supported "
                        "in pyTACS. Skipping function." % funcHandle)
            return False

        # First we will get the required domain, but only if compIDs
        # was specified. If not, just use the entire domain by default
        if compIDs is not None:
            # Make sure CompIDs is flat and get element numbers on each proc corresponding to specified compIDs
            compIDs = self._flatten(compIDs)
            elemIDs = self.meshLoader.getLocalElementIDsForComps(compIDs)
            # Finally set the domain information
            self.functionList[funcName].setDomain(elemIDs)

        return True

    def getFunctionKeys(self):
        """Return a list of the current function key names"""
        return list(self.functionList.keys())

####### Load adding methods ########

    def _addLoadToComponents(self, FVec, compIDs, F, averageLoad=False):
        """"
        The function is used to add a *FIXED TOTAL LOAD* on one or more
        components, defined by COMPIDs. The purpose of this routine is
        to add loads that remain fixed throughout an optimization. An example
        would be an engine load. This routine determines all the unqiue nodes
        in the FE model that are part of the the requested components, then
        takes the total 'force' by F and divides by the number of nodes.
        This average load is then applied to the nodes.

        NOTE: The units of the entries of the 'force' vector F are not
        necesarily physical forces and their interpretation depends
        on the physics problem being solved and the dofs included
        in the model.

        A couple of examples of force vector components for common problem are listed below:

        In Elasticity with varsPerNode = 3,
        F = [fx, fy, fz] # forces
        In Elasticity with varsPerNode = 6,
        F = [fx, fy, fz, mx, my, mz] # forces + moments
        In Thermoelasticity with varsPerNode = 4,
        F = [fx, fy, fz, Q] # forces + heat
        In Thermoelasticity with varsPerNode = 7,
        F = [fx, fy, fz, mx, my, mz, Q] # forces + moments + heat

        Parameters
        ----------

        compIDs : The components with added loads. Use selectCompIDs()
            to determine this.
        F : Numpy array length varsPerNode
            Vector of 'force' components
        """
        # Make sure CompIDs are flat
        compIDs = self._flatten([compIDs])

        # Apply a unique force vector to each component
        if not averageLoad:
            F = np.atleast_2d(F)

            # If the user only specified one force vector,
            # we assume the force should be the same for each component
            if F.shape[0] == 1:
                F = np.repeat(F, [len(compIDs)], axis=0)
            # If the dimensions still don't match, raise an error
            elif F.shape[0] != len(compIDs):
                raise Error("Number of forces must match number of compIDs,"
                            " {} forces were specified for {} compIDs".format(F.shape[0], len(compIDs)))

            # Call addLoadToComponents again, once for each compID
            for i, compID in enumerate(compIDs):
                self._addLoadToComponents(FVec, compID, F[i], averageLoad=True)

        # Average one force vector over all components
        else:
            F = np.atleast_1d(F)

            # First determine the actual physical nodal location in the
            # original BDF ordering of the nodes we want to add forces
            # to. Only the root rank need do this:
            uniqueNodes = None
            if self.comm.rank == 0:
                allNodes = []
                compIDs = set(compIDs)
                for cID in compIDs:
                    tmp = self.meshLoader.getConnectivityForComp(cID, nastranOrdering=True)
                    allNodes.extend(self._flatten(tmp))

                # Now just unique all the nodes:
                uniqueNodes = np.unique(allNodes)

            uniqueNodes = self.comm.bcast(uniqueNodes, root=0)

            # Now generate the final average force vector
            Favg = F / len(uniqueNodes)

            self._addLoadToNodes(FVec, uniqueNodes, Favg, nastranOrdering=True)

            # Write out a message of what we did:
            self._info("Added a fixed load of %s to %d components, "
                       "distributed over %d nodes." % (
                           repr(F), len(compIDs), len(uniqueNodes)),
                       maxLen=80, box=True)

    def _addLoadToNodes(self, FVec, nodeIDs, F, nastranOrdering=False):
        """
        The function is used to add a fixed point load of F to the
        selected node IDs. This is similar to the addLoadToPoints method,
        except we select the load points based on node ID rather than
        physical location.

        NOTE: This should be the prefered method (over addLoadToPoints) for adding forces to
        specific nodes for the following reasons:
            1. This method is more efficient, as it does not require a
            closest point search to locate the node.
            2. In the case where the mesh features coincident nodes
            it is impossible to uniquely specify which node gets the load
            through x,y,z location, however the points can be specified uniquely by node ID.

        A couple of examples of force vector components for common problem are listed below:

        In Elasticity with varsPerNode = 3,
        F = [fx, fy, fz] # forces
        In Elasticity with varsPerNode = 6,
        F = [fx, fy, fz, mx, my, mz] # forces + moments
        In Thermoelasticity with varsPerNode = 4,
        F = [fx, fy, fz, Q] # forces + heat
        In Thermoelasticity with varsPerNode = 7,
        F = [fx, fy, fz, mx, my, mz, Q] # forces + moments + heat

        Parameters
        ----------

        nodeIDs : list[int]
            The nodes with added loads.
        F : Numpy 1d or 2d array length (varsPerNodes) or (numNodeIDs, varsPerNodes)
            Array of force vectors, one for each node. If only one force vector is provided,
            force will be copied uniformly across all nodes.
        nastranOrdering : bool
            Flag signaling whether nodeIDs are in TACS (default)
            or NASTRAN (grid IDs in bdf file) ordering
        """

        # Make sure the inputs are the correct shape
        nodeIDs = np.atleast_1d(nodeIDs)
        F = np.atleast_2d(F)

        numNodes = len(nodeIDs)

        # If the user only specified one force vector,
        # we assume the force should be the same for each node
        if F.shape[0] == 1:
            F = np.repeat(F, [numNodes], axis=0)
        # If the dimensions still don't match, raise an error
        elif F.shape[0] != numNodes:
            raise self.TACSError("Number of forces must match number of nodes,"
                        " {} forces were specified for {} node IDs".format(F.shape[0], numNodes))

        vpn = self.assembler.getVarsPerNode()
        if len(F[0]) != vpn:
            raise self.TACSError("Length of force vector must match varsPerNode specified "
                        "for problem, which is {}, "
                        "but length of vector provided was {}".format(vpn, len(F[0])))

        # First find the cooresponding local node ID on each processor
        localNodeIDs = self.meshLoader.getLocalNodeIDsFromGlobal(nodeIDs, nastranOrdering)

        # Flag to make sure we find all user-specified nodes
        nodeFound = np.zeros(numNodes, dtype=int)

        F_array = FVec.getArray()
        nnodes = self.assembler.getNumOwnedNodes()
        vpn = self.assembler.getVarsPerNode()
        F_array = F_array.reshape(nnodes, vpn)

        # Loop through every node and if it's owned by this processor, add the load
        for i, nodeID in enumerate(localNodeIDs):
            # The node was found on this proc
            if nodeID >= 0:
                # Add contribution to global force array
                F_array[nodeID, :] += F[i]
                nodeFound[i] = 1

        # Reduce the node flag and make sure that every node was found on exactly 1 proc
        nodeFound = self.comm.allreduce(nodeFound, op=MPI.SUM)

        # Warn the user if any nodes weren't found
        if nastranOrdering:
            orderString = 'Nastran'
        else:
            orderString = 'TACS'

        for i in range(numNodes):
            if not nodeFound[i]:
                self.TACSWarning("Can't add load to node ID {} ({} ordering), node not found in model. "
                            "Double check BDF file.".format(nodeIDs[i], orderString))

    def _addTractionToComponents(self, auxElems, compIDs, tractions,
                                faceIndex=0):
        """
        The function is used to add a *FIXED TOTAL TRACTION* on one or more
        components, defined by COMPIDs. The purpose of this routine is
        to add loads that remain fixed throughout an optimization.

        Parameters
        ----------

        compIDs : The components with added loads. Use selectCompIDs()
            to determine this.
        tractions : Numpy array length 1 or compIDs
            Array of traction vectors for each components
        faceIndex : int
            Indicates which face (side) of element to apply traction to.
            Note: not required for certain elements (i.e. shells)
        """
        # Make sure compIDs is flat and unique
        compIDs = set(self._flatten(compIDs))
        tractions = np.atleast_1d(tractions)

        # Get global element IDs for the elements we're applying tractions to
        elemIDs = self.meshLoader.getGlobalElementIDsForComps(compIDs, nastranOrdering=False)
        # Add tractions element by element
        self._addTractionToElements(auxElems, elemIDs, tractions, faceIndex, nastranOrdering=False)

        # Write out a message of what we did:
        self._info("Added a fixed traction of %s to %d components, "
                   "distributed over %d elements." % (
                       repr(tractions), len(compIDs), len(elemIDs)),
                   maxLen=80, box=True)

    def _addTractionToElements(self, auxElems, elemIDs, tractions,
                              faceIndex=0, nastranOrdering=False):
        """
        The function is used to add a fixed traction to the
        selected element IDs. Tractions can be specified on an
        element by element basis (if tractions is a 2d array) or
        set to a uniform value (if tractions is a 1d array)

        Parameters
        ----------

        elemIDs : List
            The global element ID numbers for which to apply the traction.
        tractions : Numpy 1d or 2d array length varsPerNodes or (elemIDs, varsPerNodes)
            Array of traction vectors for each element
        faceIndex : int
            Indicates which face (side) of element to apply traction to.
            Note: not required for certain elements (i.e. shells)
        nastranOrdering : bool
            Flag signaling whether elemIDs are in TACS (default)
            or NASTRAN ordering
        """

        # Make sure the inputs are the correct shape
        elemIDs = np.atleast_1d(elemIDs)
        tractions = np.atleast_2d(tractions).astype(dtype=self.dtype)

        numElems = len(elemIDs)

        # If the user only specified one traction vector,
        # we assume the force should be the same for each element
        if tractions.shape[0] == 1:
            tractions = np.repeat(tractions, [numElems], axis=0)
        # If the dimensions still don't match, raise an error
        elif tractions.shape[0] != numElems:
            raise Error("Number of tractions must match number of elements,"
                        " {} tractions were specified for {} element IDs".format(tractions.shape[0], numElems))

        # First find the coresponding local element ID on each processor
        localElemIDs = self.meshLoader.getLocalElementIDsFromGlobal(elemIDs, nastranOrdering=nastranOrdering)

        # Flag to make sure we find all user-specified elements
        elemFound = np.zeros(numElems, dtype=int)

        # Loop through every element and if it's owned by this processor, add the traction
        for i, elemID in enumerate(localElemIDs):
            # The element was found on this proc
            if elemID >= 0:
                # Mark element as found
                elemFound[i] = 1
                # Get the pointer for the tacs element object for this element
                elemObj = self.meshLoader.getElementObjectForElemID(elemIDs[i], nastranOrdering=nastranOrdering)
                # Create appropriate traction object for this element type
                tracObj = elemObj.createElementTraction(faceIndex, tractions[i])
                # Traction not implemented for element
                if tracObj is None:
                    self.TACSWarning("TACS element of type {} does not hav a traction implimentation. "
                                "Skipping element in addTractionToElement procedure.".format(elemObj.getObjectName()))
                # Traction implemented
                else:
                    # Add new traction to auxiliary element object
                    auxElems.addElement(elemID, tracObj)

        # Reduce the element flag and make sure that every element was found on exactly 1 proc
        elemFound = self.comm.allreduce(elemFound, op=MPI.SUM)

        # Warn the user if any elements weren't found
        if nastranOrdering:
            orderString = 'Nastran'
        else:
            orderString = 'TACS'

        for i in range(numElems):
            if not elemFound[i]:
                self.TACSWarning("Can't add traction to element ID {} ({} ordering), element not found in model. "
                            "Double check BDF file.".format(elemIDs[i], orderString))

    def _addPressureToComponents(self, auxElems, compIDs, pressures,
                                faceIndex=0):
        """
        The function is used to add a *FIXED TOTAL PRESSURE* on one or more
        components, defined by COMPIds. The purpose of this routine is
        to add loads that remain fixed throughout an optimization. An example
        would be a fuel load.

        Parameters
        ----------

        compIDs : The components with added loads. Use selectCompIDs()
            to determine this.
        pressures : Numpy array length 1 or compIDs
            Array of pressure values for each components
        faceIndex : int
            Indicates which face (side) of element to apply pressure to.
            Note: not required for certain elements (i.e. shells)
        """
        # Make sure compIDs is flat and unique
        compIDs = set(self._flatten(compIDs))
        pressures = np.atleast_1d(pressures)

        # Get global element IDs for the elements we're applying pressure to
        elemIDs = self.meshLoader.getGlobalElementIDsForComps(compIDs, nastranOrdering=False)
        # Add pressure element by element
        self._addPressureToElements(auxElems, elemIDs, pressures, faceIndex, nastranOrdering=False)

        # Write out a message of what we did:
        self._info("Added a fixed pressure of %s to %d components, "
                   "distributed over %d elements." % (
                       repr(pressures), len(compIDs), len(elemIDs)),
                   maxLen=80, box=True)

    def _addPressureToElements(self, auxElems, elemIDs, pressures,
                              faceIndex=0, nastranOrdering=False):
        """
        The function is used to add a fixed presure to the
        selected element IDs. Pressures can be specified on an
        element by element basis (if pressures is an array) or
        set to a uniform value (if pressures is a scalar)

        Parameters
        ----------

        elemIDs : List
            The global element ID numbers for which to apply the pressure.
        pressures : Numpy array length 1 or elemIDs
            Array of pressure values for each element
        faceIndex : int
            Indicates which face (side) of element to apply pressure to.
            Note: not required for certain elements (i.e. shells)
        nastranOrdering : bool
            Flag signaling whether elemIDs are in TACS (default)
            or NASTRAN ordering
        """

        # Make sure the inputs are the correct shape
        elemIDs = np.atleast_1d(elemIDs)
        pressures = np.atleast_1d(pressures)

        numElems = len(elemIDs)

        # If the user only specified one pressure,
        # we assume the force should be the same for each element
        if pressures.shape[0] == 1:
            pressures = np.repeat(pressures, [numElems], axis=0)
        # If the dimensions still don't match, raise an error
        elif pressures.shape[0] != numElems:
            raise Error("Number of pressures must match number of elements,"
                        " {} pressures were specified for {} element IDs".format(pressures.shape[0], numElems))

        # First find the coresponding local element ID on each processor
        localElemIDs = self.meshLoader.getLocalElementIDsFromGlobal(elemIDs, nastranOrdering=nastranOrdering)

        # Flag to make sure we find all user-specified elements
        elemFound = np.zeros(numElems, dtype=int)

        # Loop through every element and if it's owned by this processor, add the pressure
        for i, elemID in enumerate(localElemIDs):
            # The element was found on this proc
            if elemID >= 0:
                elemFound[i] = 1
                # Get the pointer for the tacs element object for this element
                elemObj = self.meshLoader.getElementObjectForElemID(elemIDs[i], nastranOrdering=nastranOrdering)
                # Create appropriate pressure object for this element type
                pressObj = elemObj.createElementPressure(faceIndex, pressures[i])
                # Pressure not implemented for element
                if pressObj is None:
                    self.TACSWarning("TACS element of type {} does not hav a pressure implimentation. "
                                "Skipping element in addPressureToElement procedure.".format(elemObj.getObjectName()))
                # Pressure implemented
                else:
                    # Add new pressure to auxiliary element object
                    auxElems.addElement(elemID, pressObj)

        # Reduce the element flag and make sure that every element was found on exactly 1 proc
        elemFound = self.comm.allreduce(elemFound, op=MPI.SUM)

        # Warn the user if any elements weren't found
        if nastranOrdering:
            orderString = 'Nastran'
        else:
            orderString = 'TACS'

        for i in range(numElems):
            if not elemFound[i]:
                self.TACSWarning("Can't add pressure to element ID {} ({} ordering), element not found in model. "
                            "Double check BDF file.".format(elemIDs[i], orderString))