"""
pyBase_problem
"""

# =============================================================================
# Imports
# =============================================================================
import numpy as np
from mpi4py import MPI
from ..utilities import BaseUI
from collections import OrderedDict
import tacs.TACS


class TACSProblem(BaseUI):
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
        # pyNastran BDF object
        if self.meshLoader:
            self.bdfInfo = self.meshLoader.getBDFInfo()
        # MPI communicator object
        self.comm = comm

        # Create Design variable vector
        self.x = self.assembler.createDesignVec()
        self.assembler.getDesignVars(self.x)
        self.varName = "struct"
        # Create Nodal coordinate vector
        self.Xpts = self.assembler.createNodeVec()
        self.assembler.getNodes(self.Xpts)
        self.coordName = "Xpts"
        # List of functions
        self.functionList = OrderedDict()

        # Empty options dict, should be filled out by child class
        self.options = {}

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
        Get the current set of  design variables for this problem.

        Returns
        ----------
        x : array
            The current design variable vector set in tacs.

        """
        return self.x.getArray().copy()

    def setDesignVars(self, x):
        """
        Update the design variables used by tacs.

        Parameters
        ----------
        x : numpy.ndarray or dict or TACS.Vec
            The variables (typically from the optimizer) to set. It
            looks for variable in the ``self.varName`` attribute if in dict.

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
            raise ValueError(
                "setDesignVars must be called with either a numpy array, dict, or TACS Vec as input."
            )

        # Set the variables in tacs
        self.assembler.setDesignVars(self.x)

    def _arrayToDesignVec(self, dvArray):
        """
        Converts a distributed numpy array into a TACS design variable BVec.

        Parameters
        ----------
        dvArray : numpy.ndarray
                  Numpy array for which to convert to TACS designVec.

        Returns
        -------
        xVec : TACS.Vec
               Converted TACS designVec.

        Notes
        -----
        dvArray must have correct size on each processor.
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

    def getNodes(self):
        """
        Return the mesh coordiantes of this problem.

        Returns
        -------
        coords : array
            Structural coordinate in array of size (N * 3) where N is
            the number of structural nodes on this processor.
        """
        return self.Xpts.getArray().copy()

    def setNodes(self, Xpts):
        """
        Set the mesh coordinates of the structure.

        Parameters
        ----------
        coords : numpy.ndarray
            Structural coordinate in array of size (N * 3) where N is
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
            raise ValueError(
                "setNodes must be called with either a numpy array, dict, or TACS Vec as input."
            )
        self.assembler.setNodes(self.Xpts)

    def _arrayToNodeVec(self, xptsArray):
        """
        Converts a distributed numpy array into a TACS node BVec.

        Parameters
        ----------
        xptsArray : numpy.ndarray
                    Numpy array for which to convert to TACS nodeVec.

        Returns
        -------
        Xptsvec : TACS.Vec
                  Converted TACS nodeVec.

        Notes
        -----
        xptsArray must have correct size on each processor.
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

        Parameters
        ----------
        varArray : numpy.ndarray
                   Numpy array for which to convert to TACS Vec.

        Returns
        -------
        varVec : TACS.Vec
                 Converted TACS Vec.

        Notes
        -----
        varArray must have correct size on each processor.
        """
        varVec = self.assembler.createVec()

        # Set values
        varVec.getArray()[:] = varArray

        # Return as tacs bvec object
        return varVec

    def getNumVariables(self):
        """
        Return the number of degrees of freedom (states) that are
        on this processor

        Returns
        -------
        nstate : int
            number of states.
        """
        vpn = self.getVarsPerNode()
        nnodes = self.getNumOwnedNodes()
        return vpn * nnodes

    ####### Eval function methods ########

    def addFunction(self, funcName, funcHandle, compIDs=None, **kwargs):
        """
        Generic method to add a function for TACS. It is intended to
        be reasonably generic since the user supplies the actual
        function handle to use. See the :py:mod:`~tacs.functions` module
        for supported TACS eval functions.

        Parameters
        ----------
        funcName : str
            The user-supplied name for the function. This will
            typically be a string that is meaningful to the user

        funcHandle : TACS.Function
            The function handle to use for creation. This must come
            from the functions module in tacs.

        compIDs: list
            List of compIDs to select.

        **kwargs:
            Any keyword arguments to be passed to the TACS function during setup.
        """

        # We try to setup the function, if it fails it may not be implemented:
        try:
            # pass assembler an function-specific kwargs straight to tacs function
            self.functionList[funcName] = funcHandle(self.assembler, **kwargs)
        except:
            self._TACSWarning(
                f"Function type {funcHandle} is not currently supported. "
                "in pyTACS. Skipping function."
            )
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
        """
        Return a list of the current function key names
        """
        return list(self.functionList.keys())

    ####### Load adding methods ########

    def _addLoadToComponents(self, FVec, compIDs, F, averageLoad=False):
        """ "
        This is an internal helper function for doing the addLoadToComponents method for
        inherited TACSProblem classes. The function should NOT be called by the user should
        use the addLoadToComponents method for the respective problem class. The function is
        used to add a *FIXED TOTAL LOAD* on one or more components, defined by COMPIDs.
        The purpose of this routine is to add loads that remain fixed throughout an optimization.
        An example would be an engine load. This routine determines all the unqiue nodes in the
        FE model that are part of the requested components, then takes the total 'force' by F and
        divides by the number of nodes. This average load is then applied to the nodes.

        Parameters
        ----------

        FVec : TACS.Vec
            TACS BVec to add loads to.

        compIDs : list[int] or int
            The components with added loads. Use pyTACS.selectCompIDs method
            to determine this.

        F : numpy.ndarray 1d or 2d length (varsPerNodes) or (numCompIDs, varsPerNodes)
            Vector(s) of 'force' to apply to each components.  If only one force vector is provided,
            force will be copied uniformly across all components.

        averageLoad : bool
            Flag to determine whether load should be split evenly across all components (True)
            or copied and applied individually to each component (False). Defaults to False.

        Notes
        -----

        The units of the entries of the 'force' vector F are not
        necesarily physical forces and their interpretation depends
        on the physics problem being solved and the dofs included
        in the model.

        A couple of examples of force vector components for common problem are listed below:

            In Heat Conduction with varsPerNode = 1
                F = [Qdot] # heat rate
            In Elasticity with varsPerNode = 3,
                F = [fx, fy, fz] # forces
            In Elasticity with varsPerNode = 6,
                F = [fx, fy, fz, mx, my, mz] # forces + moments
            In Thermoelasticity with varsPerNode = 4,
                F = [fx, fy, fz, Qdot] # forces + heat rate
            In Thermoelasticity with varsPerNode = 7,
                F = [fx, fy, fz, mx, my, mz, Qdot] # forces + moments + heat rate
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
                raise self._TACSError(
                    "Number of forces must match number of compIDs,"
                    " {} forces were specified for {} compIDs".format(
                        F.shape[0], len(compIDs)
                    )
                )

            # Call addLoadToComponents again, once for each compID
            for i, compID in enumerate(compIDs):
                self._addLoadToComponents(FVec, compID, F[i], averageLoad=True)

        # Average one force vector over all components
        else:
            F = np.atleast_1d(F)

            # First determine the unique global node IDs corresponding to components:
            uniqueNodes = self.meshLoader.getGlobalNodeIDsForComps(
                compIDs, nastranOrdering=False
            )

            # Now generate the final average force vector
            Favg = F / len(uniqueNodes)

            self._addLoadToNodes(FVec, uniqueNodes, Favg, nastranOrdering=False)

            # Write out a message of what we did:
            self._info(
                "Added a fixed load of %s to %d components, "
                "distributed over %d nodes."
                % (repr(F), len(compIDs), len(uniqueNodes)),
                maxLen=80,
                box=True,
            )

    def _addLoadToNodes(self, FVec, nodeIDs, F, nastranOrdering=False):
        """
        This is an internal helper function for doing the addLoadToNodes method for
        inherited TACSProblem classes. The function should NOT be called by the user should
        use the addLoadToNodes method for the respective problem class. The function is
        used to add a fixed point load of F to the selected node IDs.

        Parameters
        ----------

        FVec : TACS.Vec
            TACS BVec to add loads to.

        nodeIDs : list[int]
            The nodes IDs with added loads.

        F : Numpy 1d or 2d array length (varsPerNodes) or (numNodeIDs, varsPerNodes)
            Array of force vectors, one for each node. If only one force vector is provided,
            force will be copied uniformly across all nodes.

        nastranOrdering : bool
            Flag signaling whether nodeIDs are in TACS (default)
            or NASTRAN (grid IDs in bdf file) ordering

        Notes
        ----------

        The units of the entries of the 'force' vector F are not
        necesarily physical forces and their interpretation depends
        on the physics problem being solved and the dofs included
        in the model.

        A couple of examples of force vector components for common problem are listed below:

            In Heat Conduction with varsPerNode = 1
                F = [Qdot] # heat rate
            In Elasticity with varsPerNode = 3,
                F = [fx, fy, fz] # forces
            In Elasticity with varsPerNode = 6,
                F = [fx, fy, fz, mx, my, mz] # forces + moments
            In Thermoelasticity with varsPerNode = 4,
                F = [fx, fy, fz, Qdot] # forces + heat rate
            In Thermoelasticity with varsPerNode = 7,
                F = [fx, fy, fz, mx, my, mz, Qdot] # forces + moments + heat rate
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
            raise self._TACSError(
                "Number of forces must match number of nodes,"
                " {} forces were specified for {} node IDs".format(F.shape[0], numNodes)
            )

        vpn = self.assembler.getVarsPerNode()
        if len(F[0]) != vpn:
            raise self._TACSError(
                "Length of force vector must match varsPerNode specified "
                "for problem, which is {}, "
                "but length of vector provided was {}".format(vpn, len(F[0]))
            )

        # First find the cooresponding local node ID on each processor
        localNodeIDs = self.meshLoader.getLocalNodeIDsFromGlobal(
            nodeIDs, nastranOrdering
        )

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
            orderString = "Nastran"
        else:
            orderString = "TACS"

        for i in range(numNodes):
            if not nodeFound[i]:
                self._TACSWarning(
                    "Can't add load to node ID {} ({} ordering), node not found in model. "
                    "Double check BDF file.".format(nodeIDs[i], orderString)
                )

    def _addLoadToRHS(self, Frhs, Fapplied):
        """ "
        This is an internal helper function for doing the addLoadToRHS method for
        inherited TACSProblem classes. The function should NOT be called by the user should
        use the addLoadToRHS method for the respective problem class.
        The function is used to add a *FIXED TOTAL LOAD* directly to the
        right hand side vector given the equation below:

            K*u = f

        Where:
            K : Stiffness matrix for problem
            u : State variables for problem
            f : Right-hand side vector to add loads to

        Parameters
        ----------

        Fapplied : numpy.ndarray or TACS.Vec
            Distributed array containing loads to applied to RHS of the problem.

        """
        if isinstance(Fapplied, tacs.TACS.Vec):
            Frhs.axpy(1.0, Fapplied)
        elif isinstance(Fapplied, np.ndarray):
            if len(Fapplied) != Frhs.getSize():
                raise self._TACSError(
                    "User-supplied distributed vector not correct length, "
                    "expected size of {} on processor {}, but got length {}.".format(
                        Frhs.getSize(), self.comm.rank, len(Fapplied)
                    )
                )
            rhsArray = Frhs.getArray()
            rhsArray[:] = rhsArray[:] + Fapplied[:]

    def _addTractionToComponents(self, auxElems, compIDs, tractions, faceIndex=0):
        """
        This is an internal helper function for doing the addTractionToComponents method for
        inherited TACSProblem classes. The function should NOT be called by the user should
        use the addTractionToComponents method for the respective problem class. The function is used
        to add a *FIXED TOTAL TRACTION* on one or more components, defined by COMPIDs. The purpose of
        this routine is to add loads that remain fixed throughout an optimization.

        Parameters
        ----------

         auxElems : TACS AuxElements object
            AuxElements object to add loads to.

        compIDs : list[int] or int
            The components with added loads. Use pyTACS.selectCompIDs method
            to determine this.

        tractions : TACS AuxElements object
            Array of traction vectors for each components

        faceIndex : int
            Indicates which face (side) of element to apply traction to.
            Note: not required for certain elements (i.e. shells)
        """
        # Make sure compIDs is flat and unique
        compIDs = set(self._flatten(compIDs))
        tractions = np.atleast_1d(tractions)

        # Get global element IDs for the elements we're applying tractions to
        elemIDs = self.meshLoader.getGlobalElementIDsForComps(
            compIDs, nastranOrdering=False
        )
        # Add tractions element by element
        self._addTractionToElements(
            auxElems, elemIDs, tractions, faceIndex, nastranOrdering=False
        )

        # Write out a message of what we did:
        self._info(
            "Added a fixed traction of %s to %d components, "
            "distributed over %d elements."
            % (repr(tractions), len(compIDs), len(elemIDs)),
            maxLen=80,
            box=True,
        )

    def _addTractionToElements(
        self, auxElems, elemIDs, tractions, faceIndex=0, nastranOrdering=False
    ):
        """
        This is an internal helper function for doing the addTractionToElements method for
        inherited TACSProblem classes. The function should NOT be called by the user should
        use the addTractionToElements method for the respective problem class. The function
        is used to add a fixed traction to the selected element IDs. Tractions can be specified on an
        element by element basis (if tractions is a 2d array) or set to a uniform value (if tractions is a 1d array)

        Parameters
        ----------

        auxElems : TACS AuxElements object
            AuxElements object to add loads to.

        elemIDs : list[int]
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
            raise self._TACSError(
                "Number of tractions must match number of elements,"
                " {} tractions were specified for {} element IDs".format(
                    tractions.shape[0], numElems
                )
            )

        # First find the coresponding local element ID on each processor
        localElemIDs = self.meshLoader.getLocalElementIDsFromGlobal(
            elemIDs, nastranOrdering=nastranOrdering
        )

        # Flag to make sure we find all user-specified elements
        elemFound = np.zeros(numElems, dtype=int)

        # Loop through every element and if it's owned by this processor, add the traction
        for i, elemID in enumerate(localElemIDs):
            # The element was found on this proc
            if elemID >= 0:
                # Mark element as found
                elemFound[i] = 1
                # Get the pointer for the tacs element object for this element
                elemObj = self.meshLoader.getElementObjectForElemID(
                    elemIDs[i], nastranOrdering=nastranOrdering
                )
                # Create appropriate traction object for this element type
                tracObj = elemObj.createElementTraction(faceIndex, tractions[i])
                # Traction not implemented for element
                if tracObj is None:
                    self._TACSWarning(
                        "TACS element of type {} does not hav a traction implimentation. "
                        "Skipping element in addTractionToElement procedure.".format(
                            elemObj.getObjectName()
                        )
                    )
                # Traction implemented
                else:
                    # Add new traction to auxiliary element object
                    auxElems.addElement(elemID, tracObj)

        # Reduce the element flag and make sure that every element was found on exactly 1 proc
        elemFound = self.comm.allreduce(elemFound, op=MPI.SUM)

        # Warn the user if any elements weren't found
        if nastranOrdering:
            orderString = "Nastran"
        else:
            orderString = "TACS"

        for i in range(numElems):
            if not elemFound[i]:
                self._TACSWarning(
                    "Can't add traction to element ID {} ({} ordering), element not found in model. "
                    "Double check BDF file.".format(elemIDs[i], orderString)
                )

    def _addPressureToComponents(self, auxElems, compIDs, pressures, faceIndex=0):
        """
        This is an internal helper function for doing the addPressureToComponents method for
        inherited TACSProblem classes. The function should NOT be called by the user should
        use the addPressureToComponents method for the respective problem class. The function
        is used to add a *FIXED TOTAL PRESSURE* on one or more components, defined by COMPIds.
        The purpose of this routine is to add loads that remain fixed throughout an optimization.
        An example would be a fuel load.

        Parameters
        ----------

        auxElems : TACS AuxElements object
            AuxElements object to add loads to.

        compIDs : list[int] or int
            The components with added loads. Use pyTACS.selectCompIDs method
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
        elemIDs = self.meshLoader.getGlobalElementIDsForComps(
            compIDs, nastranOrdering=False
        )
        # Add pressure element by element
        self._addPressureToElements(
            auxElems, elemIDs, pressures, faceIndex, nastranOrdering=False
        )

        # Write out a message of what we did:
        self._info(
            "Added a fixed pressure of %s to %d components, "
            "distributed over %d elements."
            % (repr(pressures), len(compIDs), len(elemIDs)),
            maxLen=80,
            box=True,
        )

    def _addPressureToElements(
        self, auxElems, elemIDs, pressures, faceIndex=0, nastranOrdering=False
    ):
        """
        This is an internal helper function for doing the addPressureToElements method for
        inherited TACSProblem classes. The function should NOT be called by the user should
        use the addPressureToElements method for the respective problem class. The function
        is used to add a fixed presure to the selected element IDs. Pressures can be specified on an
        element by element basis (if pressures is an array) or set to a uniform value (if pressures is a scalar)

        Parameters
        ----------

        auxElems : TACS AuxElements object
            AuxElements object to add loads to.

        elemIDs : list[int]
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
            raise self._TACSError(
                "Number of pressures must match number of elements,"
                " {} pressures were specified for {} element IDs".format(
                    pressures.shape[0], numElems
                )
            )

        # First find the coresponding local element ID on each processor
        localElemIDs = self.meshLoader.getLocalElementIDsFromGlobal(
            elemIDs, nastranOrdering=nastranOrdering
        )

        # Flag to make sure we find all user-specified elements
        elemFound = np.zeros(numElems, dtype=int)

        # Loop through every element and if it's owned by this processor, add the pressure
        for i, elemID in enumerate(localElemIDs):
            # The element was found on this proc
            if elemID >= 0:
                elemFound[i] = 1
                # Get the pointer for the tacs element object for this element
                elemObj = self.meshLoader.getElementObjectForElemID(
                    elemIDs[i], nastranOrdering=nastranOrdering
                )
                # Create appropriate pressure object for this element type
                pressObj = elemObj.createElementPressure(faceIndex, pressures[i])
                # Pressure not implemented for element
                if pressObj is None:
                    self._TACSWarning(
                        "TACS element of type {} does not hav a pressure implimentation. "
                        "Skipping element in addPressureToElement procedure.".format(
                            elemObj.getObjectName()
                        )
                    )
                # Pressure implemented
                else:
                    # Add new pressure to auxiliary element object
                    auxElems.addElement(elemID, pressObj)

        # Reduce the element flag and make sure that every element was found on exactly 1 proc
        elemFound = self.comm.allreduce(elemFound, op=MPI.SUM)

        # Warn the user if any elements weren't found
        if nastranOrdering:
            orderString = "Nastran"
        else:
            orderString = "TACS"

        for i in range(numElems):
            if not elemFound[i]:
                self._TACSWarning(
                    "Can't add pressure to element ID {} ({} ordering), element not found in model. "
                    "Double check BDF file.".format(elemIDs[i], orderString)
                )

    def _addInertialLoad(self, auxElems, inertiaVector):
        """
        This is an internal helper function for doing the addInertialLoad method for
        inherited TACSProblem classes. The function should NOT be called by the user should
        use the addInertialLoad method for the respective problem class. The function
        is used to add a fixed inertial load due to a uniform acceleration over the entire model.
        This is most commonly used to model gravity loads on a model.

        Parameters
        ----------

         auxElems : TACS AuxElements object
            AuxElements object to add loads to.

        inertiaVector : numpy.ndarray
            Acceleration vector used to define inertial load.
        """
        # Make sure vector is right type
        inertiaVector = np.atleast_1d(inertiaVector).astype(self.dtype)
        # Get elements on this processor
        localElements = self.assembler.getElements()
        # Loop through every element and apply inertial load
        for elemID, elemObj in enumerate(localElements):
            # Create appropriate inertial force object for this element type
            inertiaObj = elemObj.createElementInertialForce(inertiaVector)
            # Inertial force is implemented for element
            if inertiaObj is not None:
                # Add new inertial force to auxiliary element object
                auxElems.addElement(elemID, inertiaObj)

    def _addCentrifugalLoad(self, auxElems, omegaVector, rotCenter):
        """
        This is an internal helper function for doing the addCentrifugalLoad method for
        inherited TACSProblem classes. The function should NOT be called by the user should
        use the addCentrifugalLoad method for the respective problem class. The function
        is used to add a fixed centrifugal load due to a uniform rotational velocity over the entire model.
        This is most commonly used to model rotors, rolling aircraft, etc.

        Parameters
        ----------

        auxElems : TACS AuxElements object
            AuxElements object to add loads to.

        omegaVector : numpy.ndarray
            Rotational velocity vector (rad/s) used to define centrifugal load.

        rotCenter : numpy.ndarray
            Location of center of rotation used to define centrifugal load.
        """
        # Make sure vector is right type
        omegaVector = np.atleast_1d(omegaVector).astype(self.dtype)
        rotCenter = np.atleast_1d(rotCenter).astype(self.dtype)
        # Get elements on this processor
        localElements = self.assembler.getElements()
        # Loop through every element and apply centrifugal load
        for elemID, elemObj in enumerate(localElements):
            # Create appropriate centrifugal force object for this element type
            centrifugalObj = elemObj.createElementCentrifugalForce(
                omegaVector, rotCenter
            )
            # Centrifugal force is implemented for element
            if centrifugalObj is not None:
                # Add new centrifugal force to auxiliary element object
                auxElems.addElement(elemID, centrifugalObj)

    def _addLoadFromBDF(self, FVec, auxElems, loadID, setScale=1.0):
        """
        This is an internal helper function for doing the addLoadFromBDF method for
        inherited TACSProblem classes. The function should NOT be called by the user should
        use the addLoadFromBDF method for the respective problem class. This method is
        used to add a fixed load set defined in the BDF file to the problem.
        Currently, only supports LOAD, FORCE, MOMENT, GRAV, RFORCE, PLOAD2, and PLOAD4.

        Parameters
        ----------

        FVec : TACS.Vec
            TACS BVec to add loads to.

        auxElems : TACS AuxElements object
            AuxElements object to add loads to.

        loadID : int
            Load identification number of load set in BDF file user wishes to add to problem.

        setScale : float
            Factor to scale the BDF loads by before adding to problem.
        """
        vpn = self.assembler.getVarsPerNode()
        # Get loads and scalers for this load case ID
        loadSet, loadScale, _ = self.bdfInfo.get_reduced_loads(loadID)
        # Loop through every load in set and add it to problem
        for loadInfo, scale in zip(loadSet, loadScale):
            scale *= setScale
            # Add any point force or moment cards
            if loadInfo.type == "FORCE" or loadInfo.type == "MOMENT":
                nodeID = loadInfo.node_ref.nid

                loadArray = np.zeros(vpn)
                if loadInfo.type == "FORCE" and vpn >= 3:
                    F = scale * loadInfo.scaled_vector
                    loadArray[:3] += loadInfo.cid_ref.transform_vector_to_global(F)
                elif loadInfo.type == "MOMENT" and vpn >= 6:
                    M = scale * loadInfo.scaled_vector
                    loadArray[3:6] += loadInfo.cid_ref.transform_vector_to_global(M)
                self._addLoadToNodes(FVec, nodeID, loadArray, nastranOrdering=True)

            # Add any gravity loads
            elif loadInfo.type == "GRAV":
                inertiaVec = np.zeros(3, dtype=self.dtype)
                inertiaVec[:3] = scale * loadInfo.scale * loadInfo.N
                # Convert acceleration to global coordinate system
                inertiaVec = loadInfo.cid_ref.transform_vector_to_global(inertiaVec)
                self._addInertialLoad(auxElems, inertiaVec)

            # Add any centrifugal loads
            elif loadInfo.type == "RFORCE":
                omegaVec = np.zeros(3, dtype=self.dtype)
                if loadInfo.nid_ref:
                    rotCenter = loadInfo.nid_ref.get_position()
                else:
                    rotCenter = np.zeros(3, dtype=self.dtype)
                omegaVec[:3] = scale * loadInfo.scale * np.array(loadInfo.r123)
                # Convert omega from rev/s to rad/s
                omegaVec *= 2 * np.pi
                # Convert omega to global coordinate system
                omegaVec = loadInfo.cid_ref.transform_vector_to_global(omegaVec)
                self._addCentrifugalLoad(auxElems, omegaVec, rotCenter)

            # Add any pressure loads
            # Pressure load card specific to shell elements
            elif loadInfo.type == "PLOAD2":
                elemIDs = loadInfo.eids
                pressure = scale * loadInfo.pressure
                self._addPressureToElements(
                    auxElems, elemIDs, pressure, nastranOrdering=True
                )

            # Alternate more general pressure load type
            elif loadInfo.type == "PLOAD4":
                self._addPressureFromPLOAD4(auxElems, loadInfo, scale)

            else:
                self._TACSWarning(
                    "Unsupported load type "
                    f" '{loadInfo.type}' specified for load set number {loadInfo.sid}, skipping load"
                )

    def _addPressureFromPLOAD4(self, auxElems, loadInfo, scale=1.0):
        """
        Add pressure to tacs static/transient problem from pynastran PLOAD4 card.
        Should only be called by createTACSProbsFromBDF and not directly by user.
        """
        # Dictionary mapping nastran element face indices to TACS equivilent numbering
        nastranToTACSFaceIDDict = {
            "CTETRA4": {1: 1, 2: 3, 3: 2, 4: 0},
            "CTETRA": {2: 1, 4: 3, 3: 2, 1: 0},
            "CHEXA": {1: 4, 2: 2, 3: 0, 4: 3, 5: 0, 6: 5},
        }

        # We don't support pressure variation across elements, for now just average it
        pressure = scale * np.mean(loadInfo.pressures)
        for elemInfo in loadInfo.eids_ref:
            elemID = elemInfo.eid

            # Get the correct face index number based on element type
            if "CTETRA" in elemInfo.type:
                for faceIndex in elemInfo.faces:
                    if (
                        loadInfo.g1 in elemInfo.faces[faceIndex]
                        and loadInfo.g34 not in elemInfo.faces[faceIndex]
                    ):
                        # For some reason CTETRA4 is the only element that doesn't
                        # use ANSYS face numbering convention by default
                        if len(elemInfo.nodes) == 4:
                            faceIndex = nastranToTACSFaceIDDict["CTETRA4"][faceIndex]
                        else:
                            faceIndex = nastranToTACSFaceIDDict["CTETRA"][faceIndex]
                        # Positive pressure is inward for solid elements, flip pressure if necessary
                        # We don't flip it for face 0, because the normal for that face points inward by convention
                        # while the rest point outward
                        if faceIndex != 0:
                            pressure *= -1.0
                        break

            elif "CHEXA" in elemInfo.type:
                for faceIndex in elemInfo.faces:
                    if (
                        loadInfo.g1 in elemInfo.faces[faceIndex]
                        and loadInfo.g34 in elemInfo.faces[faceIndex]
                    ):
                        faceIndex = nastranToTACSFaceIDDict["CHEXA"][faceIndex]
                        # Pressure orientation is flipped for solid elements per Nastran convention
                        pressure *= -1.0
                        break

            elif "CQUAD" in elemInfo.type or "CTRIA" in elemInfo.type:
                # Face index doesn't matter for shells, just use 0
                faceIndex = 0

            else:
                raise self._TACSError(
                    "Unsupported element type "
                    f"'{elemInfo.type}' specified for PLOAD4 load set number {loadInfo.sid}."
                )

            # Figure out whether or not this is a traction based on if a vector is defined
            if np.linalg.norm(loadInfo.nvector) == 0.0:
                self._addPressureToElements(
                    auxElems, elemID, pressure, faceIndex, nastranOrdering=True
                )
            else:
                trac = pressure * loadInfo.nvector
                self._addTractionToElements(
                    auxElems, elemID, trac, faceIndex, nastranOrdering=True
                )
