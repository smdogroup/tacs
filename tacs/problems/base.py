"""
pyBase_problem
"""

# =============================================================================
# Imports
# =============================================================================
from collections import OrderedDict

import numpy as np
from mpi4py import MPI

import tacs.TACS
from ..system import TACSSystem


class TACSProblem(TACSSystem):
    """
    Base class for TACS problem types. Contains methods common to all TACS problems.
    """

    def __init__(
        self,
        assembler,
        comm=None,
        options=None,
        outputViewer=None,
        meshLoader=None,
        isNonlinear=False,
    ):
        """
        Parameters
        ----------
        assembler : tacs.TACS.Assembler
            Cython object responsible for creating and setting tacs objects used to solve problem

        comm : mpi4py.MPI.Intracomm
            The comm object on which to create the pyTACS object.

        options : dict
            Dictionary holding problem-specific option parameters (case-insensitive).

        outputViewer : tacs.TACS.TACSToFH5
            Cython object used to write out f5 files that can be converted and used for postprocessing.

        meshLoader : tacs.pymeshloader.pyMeshLoader
            pyMeshLoader object used to create the assembler.
        """

        # Set nonlinear flag
        self._isNonlinear = isNonlinear

        # Set attributes and options
        TACSSystem.__init__(self, assembler, comm, options, outputViewer, meshLoader)

        # List of functions
        self.functionList = OrderedDict()

        return

    @property
    def isNonlinear(self):
        """The public interface for the isNonlinear attribute. Implemented as a property so that it is read-only."""
        return self._isNonlinear

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

        Returns
        -------
        funcNames : list[str]
            List containing user-defined names for functions added so far.
        """
        return list(self.functionList.keys())

    ####### Load adding methods ########

    def _addLoadToComponents(self, FVec, compIDs, F, averageLoad=False):
        """
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

        FVec : tacs.TACS.Vec
            TACS BVec to add loads to.

        compIDs : list[int] or int
            The components with added loads. Use pyTACS.selectCompIDs method
            to determine this.

        F : numpy.ndarray 1d or 2d length (varsPerNodes) or (numCompIDs, varsPerNodes)
            Vector(s) of 'force' to apply to each component.  If only one force vector is provided,
            force will be copied uniformly across all components.

        averageLoad : bool
            Flag to determine whether load should be split evenly across all components (True)
            or copied and applied individually to each component (False). Defaults to False.

        Notes
        -----

        The units of the entries of the 'force' vector F are not
        necessarily physical forces and their interpretation depends
        on the physics problem being solved and the dofs included
        in the model.

        A couple of examples of force vector components for common problems are listed below:

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

        FVec : tacs.TACS.Vec
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
        -----

        The units of the entries of the 'force' vector F are not
        necessarily physical forces and their interpretation depends
        on the physics problem being solved and the dofs included
        in the model.

        A couple of examples of force vector components for common problems are listed below:

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

        # First find the corresponding local node ID on each processor
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
        """
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

        Fapplied : numpy.ndarray or tacs.TACS.Vec
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

         auxElems : tacs.TACS.AuxElements
            AuxElements object to add loads to.

        compIDs : list[int] or int
            The components with added loads. Use pyTACS.selectCompIDs method
            to determine this.

        tractions : numpy.ndarray
            Array of traction vectors for each component

        faceIndex : int
            Indicates which face (side) of the element to apply traction to.
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

        auxElems : tacs.TACS.AuxElements
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

        # First find the corresponding local element ID on each processor
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
                # Create an appropriate traction object for this element type
                tracObj = elemObj.createElementTraction(faceIndex, tractions[i])
                # Traction not implemented for this element
                if tracObj is None:
                    self._TACSWarning(
                        "TACS element of type {} does not hav a traction implementation. "
                        "Skipping element in addTractionToElement procedure.".format(
                            elemObj.getObjectName()
                        )
                    )
                # Traction implemented
                else:
                    # Add new traction to the auxiliary element object
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

        auxElems : tacs.TACS.AuxElements
            AuxElements object to add loads to.

        compIDs : list[int] or int
            The components with added loads. Use pyTACS.selectCompIDs method
            to determine this.

        pressures : Numpy array length 1 or compIDs
            Array of pressure values for each component

        faceIndex : int
            Indicates which face (side) of the element to apply pressure to.
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

        auxElems : tacs.TACS.AuxElements
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

        # First find the corresponding local element ID on each processor
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
                        "TACS element of type {} does not hav a pressure implementation. "
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

         auxElems : tacs.TACS.AuxElements
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

    def _addCentrifugalLoad(self, auxElems, omegaVector, rotCenter, firstOrder=False):
        """
        This is an internal helper function for doing the addCentrifugalLoad method for
        inherited TACSProblem classes. The function should NOT be called by the user should
        use the addCentrifugalLoad method for the respective problem class. The function
        is used to add a fixed centrifugal load due to a uniform rotational velocity over the entire model.
        This is most commonly used to model rotors, rolling aircraft, etc.

        Parameters
        ----------

        auxElems : tacs.TACS.AuxElements
            AuxElements object to add loads to.

        omegaVector : numpy.ndarray
            Rotational velocity vector (rad/s) used to define centrifugal load.

        rotCenter : numpy.ndarray
            Location of center of rotation used to define centrifugal load.

        firstOrder : bool, optional
            Whether to use first order approximation for centrifugal load,
            which computes the force in the displaced position. By default False
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
                omegaVector, rotCenter, firstOrder=firstOrder
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

        FVec : tacs.TACS.Vec
            TACS BVec to add loads to.

        auxElems : tacs.TACS.AuxElements
            AuxElements object to add loads to.

        loadID : int
            Load identification number of load set in BDF file user wishes to add to problem.

        setScale : float
            Factor to scale the BDF loads by before adding to problem.
        """
        # Make sure bdf is cross-referenced
        if self.bdfInfo.is_xrefed is False:
            self.bdfInfo.cross_reference()
            self.bdfInfo.is_xrefed = True

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
        # Dictionary mapping nastran element face indices to TACS equivalent numbering
        nastranToTACSFaceIDDict = {
            "CTETRA4": {1: 1, 2: 3, 3: 2, 4: 0},
            "CTETRA": {2: 1, 4: 3, 3: 2, 1: 0},
            "CHEXA": {1: 4, 2: 2, 3: 1, 4: 3, 5: 0, 6: 5},
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

    def writeSensFile(self, evalFuncs, tacsAim, proc: int = 0, root=0):
        """
        write an ESP/CAPS .sens file from the tacs aim
        Optional tacs_aim arg for TacsAim wrapper class object in root/tacs/caps2tacs/

        Parameters
        ----------
        evalFuncs : list[str]
            names of TACS functions to be evaluated
        tacsAim : tacs.caps2tacs.TacsAIM
            class which handles the sensitivity file writing for ESP/CAPS shape derivatives
        proc: int
            which processor (in case of parallel tacsAIM instances) to write the sens file to
        root: int
            which processor we get the data from (usually proc 0)
        """

        is_dummy_file = evalFuncs is None
        if is_dummy_file:
            evalFuncs = ["dummy-func"]

        # obtain the functions and sensitivities from TACS assembler
        tacs_funcs = {}
        if not (is_dummy_file):
            local_tacs_sens = {}
            self.evalFunctions(tacs_funcs, evalFuncs=evalFuncs)
            self.evalFunctionsSens(local_tacs_sens, evalFuncs=evalFuncs)

            # get nastran=>local tacs id map for this processor
            # this local_tacs_ids is len(num_nodes) globally with -1 for nodes
            # not on this processor and the local tacs_ids for nodes owned by this processor
            num_nodes = self.meshLoader.bdfInfo.nnodes
            bdfNodes = range(num_nodes)
            local_tacs_ids = self.meshLoader.getLocalNodeIDsFromGlobal(
                bdfNodes, nastranOrdering=False
            )

            # only need to combine xpts sens across processors, the struct sens is only root proc (usually rank 0)
            local_xpt_sens = {}
            for tacs_key in local_tacs_sens:
                local_xpt_sens[tacs_key] = local_tacs_sens[tacs_key]["Xpts"]

            # gather the tacs ids and tacs xpts sens globally
            mpi_tacs_ids = self.comm.gather(local_tacs_ids, root=root)
            mpi_xpt_sens = self.comm.gather(local_xpt_sens, root=root)

        num_funcs = len(evalFuncs)
        assert tacsAim is not None
        num_struct_dvs = len(tacsAim.thickness_variables)
        num_nodes = self.meshLoader.bdfInfo.nnodes

        # uses other proc and broadcast so needed before root-proc check
        sens_file_path = tacsAim.sens_file_path(proc)

        if self.comm.rank == root:
            # open the sens file nastran_CAPS.sens and write coordinate derivatives
            # and any other struct derivatives to it
            with open(sens_file_path, "w") as hdl:
                for func_name in evalFuncs:
                    hdl.write(f"{num_funcs} {num_struct_dvs}\n")

                    # for each function write the values and coordinate derivatives
                    for func_name in evalFuncs:
                        # get the tacs key
                        for tacs_key in tacs_funcs:
                            if func_name in tacs_key:
                                break

                        if not (is_dummy_file):
                            # write the func name, value and nnodes
                            hdl.write(f"{func_name}\n")
                            hdl.write(f"{tacs_funcs[tacs_key].real}\n")
                            hdl.write(f"{num_nodes}\n")

                            # write the coordinate derivatives for the given function
                            for nastran_id_m1 in range(num_nodes):
                                dfdxyz = None
                                # find the xpt-sens among those of each processor in the mpi_tacs_sens list
                                for ilocal, local_tacs_ids in enumerate(mpi_tacs_ids):
                                    if local_tacs_ids[nastran_id_m1] != -1:
                                        local_tacs_id = local_tacs_ids[nastran_id_m1]
                                        dfdxyz = mpi_xpt_sens[ilocal][tacs_key][
                                            3 * local_tacs_id : 3 * local_tacs_id + 3
                                        ]
                                        break

                                nastran_id = nastran_id_m1 + 1
                                hdl.write(
                                    f"{nastran_id} {dfdxyz[0].real} {dfdxyz[1].real} {dfdxyz[2].real}\n"
                                )

                            # write any struct derivatives if there are struct derivatives
                            if num_struct_dvs > 0:
                                struct_sens = local_tacs_sens[tacs_key]["struct"]
                                for idx, thick_var in enumerate(
                                    tacsAim.thickness_variables
                                ):
                                    # assumes these are sorted in tacs aim wrapper
                                    hdl.write(f"{thick_var.name}\n")
                                    hdl.write("1\n")
                                    hdl.write(f"{struct_sens[idx].real}\n")

                        else:  # is a dummy sens file for animating the structure shape
                            # write the func name, value and nnodes
                            hdl.write(f"{func_name}\n")
                            hdl.write(f"{0.0}\n")
                            hdl.write(f"{num_nodes}\n")

                            # write the coordinate derivatives for the given function
                            for bdf_ind in range(num_nodes):
                                nastran_node = bdf_ind + 1
                                hdl.write(f"{nastran_node} 0.0 0.0 0.0\n")

                            # write any struct derivatives if there are struct derivatives
                            if num_struct_dvs > 0:
                                for idx, thick_var in enumerate(
                                    tacsAim.thickness_variables
                                ):
                                    # assumes these are sorted in tacs aim wrapper
                                    hdl.write(f"{thick_var.name}\n")
                                    hdl.write("1\n")
                                    hdl.write("0.0\n")

            return
