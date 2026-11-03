"""
pySystem
"""

# =============================================================================
# Imports
# =============================================================================

import numpy as np

from tacs.utilities import BaseUI


class TACSSystem(BaseUI):
    """
    Base class for TACS problem/constraint types. Contains methods common to all TACS systems dealing with design variables.
    """

    def __init__(
        self, assembler, comm=None, options=None, outputViewer=None, meshLoader=None
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
        # TACS assembler object
        self.assembler = assembler
        # TACS F5 output writer
        self.outputViewer = outputViewer
        # TACS pyMeshLoader object
        self.meshLoader = meshLoader
        # pyNastran BDF object
        if self.meshLoader:
            self.bdfInfo = self.meshLoader.getBDFInfo()

        # Create Design variable vector
        self.x = self.assembler.createDesignVec()
        self.assembler.getDesignVars(self.x)
        self.varName = "struct"
        # Create Nodal coordinate vector
        self.Xpts = self.assembler.createNodeVec()
        self.assembler.getNodes(self.Xpts)
        self.coordName = "Xpts"

        # Setup comm and options
        BaseUI.__init__(self, options=options, comm=comm)

        return

    ####### Design variable methods ########

    def setVarName(self, varName):
        """
        Set a name for the design variables in pyOpt. Only needs
        to be changed if more than 1 pytacs object is used in an
        optimization

        Parameters
        ----------
        varName : str
            Name of the design variables used in setDesignVars() dict.
        """
        self.varName = varName

    def getVarName(self):
        """
        Get name for the design variables in pyOpt. Only needed
        if more than 1 pytacs object is used in an optimization

        Returns
        -------
        varName : str
            Name of the design variables used in setDesignVars() dict.
        """
        return self.varName

    def getDesignVars(self):
        """
        Get the current set of  design variables for this problem.

        Returns
        -------
        x : numpy.ndarray
            The current design variable vector set in tacs.

        """
        return self.x.getArray().copy()

    def setDesignVars(self, x):
        """
        Update the design variables used by tacs.

        Parameters
        ----------
        x : numpy.ndarray or dict or tacs.TACS.Vec
            The variables (typically from the optimizer) to set. It
            looks for variable in the ``self.varName`` attribute if in dict.

        """
        try:
            # Check if the design variables are being handed in a dict
            if isinstance(x, dict):
                if self.varName in x:
                    self.copyToTACSVec(x[self.varName], self.x)
            # or array or TACSBVec
            else:
                self.copyToTACSVec(x, self.x)
        except ValueError as err:
            raise ValueError(
                "setDesignVars must be called with either a numpy array, TACS Vec, or dict containing one of the two, as input."
            ) from err

        # Set the variables in tacs
        self.assembler.setDesignVars(self.x)

    def getDesignVarRange(self):
        """
        Get the lower/upper bounds for the design variables.

        Returns
        -------
        xlb : numpy.ndarray
            The design variable lower bound.
        xub : numpy.ndarray
            The design variable upper bound.

        """
        xlb = self.assembler.createDesignVec()
        xub = self.assembler.createDesignVec()
        self.assembler.getDesignVarRange(xlb, xub)
        return xlb.getArray(), xub.getArray()

    def getComponentDesignVars(self):
        """
        Get the design variable groups and their current values for every component.

        The returned dictionary is keyed by component description and each entry maps design
        variable group names to their current values. Group names and value types match the
        keyword arguments of the corresponding constitutive class constructor, so the values
        can be passed directly back to the constructor inside an elemCallBack function to
        recreate the current sizing state in another TACS execution.

        One exception: the ``lp`` group of ``LamParamFullShellConstitutive`` cannot be set
        through the constructor; restore it with ``setLaminationParameters`` after
        construction.

        All design variable groups are always returned, whether or not their entries are
        active design variables. Active entries reflect this problem's or constraint's
        current design variable values, i.e. the values last passed to its
        :meth:`setDesignVars` method; inactive entries hold the values the constitutive
        object was constructed with. Use :meth:`getComponentDesignVarNums` to determine
        which entries are active.

        This method is collective on this object's comm; the returned dictionary is
        identical on every processor.

        Returns
        -------
        compDVs : dict[str, dict[str, Union[float, np.ndarray]]]
            Dictionary of design variable group values for each component. Scalar groups
            are returned as scalars and array groups as 1D numpy arrays.
        """
        if self.meshLoader is None:
            raise self._TACSError(
                "getComponentDesignVars requires a meshLoader. Create this "
                "problem/constraint through the pyTACS factory methods (e.g. "
                "createStaticProblem)."
            )
        self._checkUniqueComponentDescripts()
        compDescripts = self.meshLoader.getComponentDescripts()
        nComp = self.meshLoader.getNumComponents()
        globalDVValues = self._gatherGlobalDesignVarValues()
        compDVs = {}
        for compID in range(nComp):
            descript = compDescripts[compID]
            con = self._getComponentConstitutive(compID)
            groupValues = {}
            if con is not None:
                self._checkDVGroupsImplemented(compID, con)
                groupValues = con.getDesignVarGroups()
                groupNums = con.getDesignVarGroupDVNums()
                # Overlay this object's current values onto the active entries
                for name, nums in groupNums.items():
                    if isinstance(nums, np.ndarray):
                        for ii in range(len(nums)):
                            if nums[ii] >= 0:
                                groupValues[name][ii] = globalDVValues[nums[ii]]
                    elif nums >= 0:
                        groupValues[name] = globalDVValues[nums].item()
            compDVs[descript] = groupValues
        return compDVs

    def getComponentDesignVarNums(self):
        """
        Get the global design variable numbers of every design variable group for every
        component.

        The returned dictionary has the same structure as the one returned by
        :meth:`getComponentDesignVars`, but contains global design variable numbers instead
        of values. Entries that are not active design variables have a design variable
        number of -1. Note that design variable numbers are specific to this TACS execution
        and should not be stored across executions.

        Design variable numbers are shared by all problems and constraints within one TACS
        execution, so this dictionary does not depend on which problem or constraint it is
        called on.

        Returns
        -------
        compDVNums : dict[str, dict[str, Union[int, np.ndarray]]]
            Dictionary of design variable group numbers for each component. Scalar groups
            are returned as integers and array groups as 1D numpy integer arrays.
        """
        if self.meshLoader is None:
            raise self._TACSError(
                "getComponentDesignVarNums requires a meshLoader. Create this "
                "problem/constraint through the pyTACS factory methods (e.g. "
                "createStaticProblem)."
            )
        self._checkUniqueComponentDescripts()
        compDescripts = self.meshLoader.getComponentDescripts()
        nComp = self.meshLoader.getNumComponents()
        compDVNums = {}
        for compID in range(nComp):
            descript = compDescripts[compID]
            con = self._getComponentConstitutive(compID)
            groupNums = {}
            if con is not None:
                self._checkDVGroupsImplemented(compID, con)
                groupNums = con.getDesignVarGroupDVNums()
            compDVNums[descript] = groupNums
        return compDVNums

    def _checkUniqueComponentDescripts(self):
        """Raise an error if any two components share a description, since the component
        design variable dictionaries are keyed by description.
        """
        compDescripts = self.meshLoader.getComponentDescripts()
        seen = set()
        duplicates = set()
        for descript in compDescripts:
            if descript in seen:
                duplicates.add(descript)
            seen.add(descript)
        if len(duplicates) > 0:
            raise self._TACSError(
                "Component descriptions must be unique to use getComponentDesignVars/"
                "getComponentDesignVarNums, but the following descriptions are shared by "
                f"multiple components: {sorted(duplicates)}. Rename the property groups in "
                "the BDF file so that every component has a unique description."
            )

    def _getComponentConstitutive(self, compID):
        """Get the single constitutive object shared by all element objects of a component.

        Returns None if the component has no constitutive object. Raises an error if the
        component's element objects hold more than one distinct constitutive object, since
        the component design variable dictionaries cannot represent that faithfully.
        """
        compDescripts = self.meshLoader.getComponentDescripts()
        conObj = None
        numObjs = len(self.meshLoader.getElementObjectNumsForComp(compID))
        for objIndex in range(numObjs):
            elemObj = self.meshLoader.getElementObject(compID, objIndex)
            elemCon = elemObj.getConstitutive()
            if elemCon is None:
                continue
            if conObj is None:
                conObj = elemCon
            elif elemCon is not conObj:
                raise self._TACSError(
                    f"Component {compID} ('{compDescripts[compID]}') uses more than "
                    "one distinct constitutive object across its element types. "
                    "getComponentDesignVars/getComponentDesignVarNums require all element "
                    "types in a component to share a single constitutive object."
                )
        return conObj

    def _gatherGlobalDesignVarValues(self):
        """Gather this problem's/constraint's current design variable values into a
        rank-identical numpy array indexed by global design variable number.

        The design vector is distributed in contiguous blocks of design variable numbers,
        so concatenating the local arrays in rank order yields the global array.
        """
        localValues = self.x.getArray().copy()
        return np.concatenate(self.comm.allgather(localValues))

    def _checkDVGroupsImplemented(self, compID, con):
        """Warn if a constitutive object's design variable groups do not cover all of its
        design variables, since the missing DVs will be absent from the output.
        """
        compDescripts = self.meshLoader.getComponentDescripts()
        numDVs = len(con.getDesignVarNums())
        if con.getNumDesignVarGroups() == 0 and numDVs > 0:
            self._TACSWarning(
                f"Constitutive class '{type(con).__name__}' used by component {compID} "
                f"('{compDescripts[compID]}') has design variables but does not "
                "implement the design variable group interface. Its design variables will "
                "be missing from the component design variable dictionaries."
            )
        else:
            numActiveGroupEntries = 0
            for nums in con.getDesignVarGroupDVNums().values():
                numActiveGroupEntries += int(np.count_nonzero(np.atleast_1d(nums) >= 0))
            if numActiveGroupEntries != numDVs:
                self._TACSWarning(
                    f"Constitutive class '{type(con).__name__}' used by component {compID} "
                    f"('{compDescripts[compID]}') has {numDVs} design variables but its "
                    f"design variable groups only cover {numActiveGroupEntries} of them. The "
                    "remaining design variables will be missing from the component design "
                    "variable dictionaries."
                )

    def _arrayToDesignVec(self, dvArray):
        """
        Converts a distributed numpy array into a TACS design variable BVec.

        Parameters
        ----------
        dvArray : numpy.ndarray
                  Numpy array for which to convert to TACS designVec.

        Returns
        -------
        xVec : tacs.TACS.Vec
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

        Returns
        -------
        ndvs : int
            Number of design variables on this processor.
        """
        return self.x.getSize()

    def setCoordName(self, coordName):
        """
        Set a name for the nodal coordinates in pyOpt.

        Parameters
        ----------
        coordName : str
            Name of the nodal coordinates used in setNodes() dict.
        """
        self.coordName = coordName

    def getCoordName(self):
        """
        Get name for the nodal coordinates in pyOpt.

        Returns
        -------
        coordName : str
            Name of the nodal coordinates used in setNodes() dict.
        """
        return self.coordName

    def getNodes(self):
        """
        Return the mesh coordinates of this problem.

        Returns
        -------
        coords : numpy.ndarray
            Structural coordinate in array of size (N * 3) where N is
            the number of structural nodes on this processor.
        """
        return self.Xpts.getArray().copy()

    def setNodes(self, Xpts):
        """
        Set the mesh coordinates of the structure.

        Parameters
        ----------
        Xpts : numpy.ndarray
            Structural coordinate in array of size (N * 3) where N is
            the number of structural nodes on this processor.
        """
        try:
            # Check if the coordinates are being handed in a dict
            if isinstance(Xpts, dict):
                if self.coordName in Xpts:
                    self.copyToTACSVec(Xpts[self.coordName], self.Xpts)
            # or array or TACSBVec
            else:
                self.copyToTACSVec(Xpts, self.Xpts)
        except ValueError as err:
            raise ValueError(
                "setNodes must be called with either a numpy array, TACS Vec, or dict containing one of the two, as input."
            ) from err
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
        Xptsvec : tacs.TACS.Vec
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

        Returns
        -------
        ncoords : int
            Number of mesh coordinates on this processor.
        """
        return self.Xpts.getSize()

    ####### Variable methods ########

    def getVarsPerNode(self):
        """
        Get the number of variables per node for the model.

        Returns
        -------
        vpn : int
            Number of variables per node.
        """
        return self.assembler.getVarsPerNode()

    def getNumOwnedNodes(self):
        """
        Get the number of nodes owned by this processor.

        Returns
        -------
        nnodes : int
            Number of nodes on this processor.
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
        varVec : tacs.TACS.Vec
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
