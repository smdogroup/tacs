#!/usr/bin/python
"""
pytacs - The Python wrapper for the TACS assembler

This python interface is designed to provide a easier interface to the
C++ layer of TACS. User-supplied hooks allow for nearly complete
customization of any or all parts in the problem setup. There are two
main parts of this module: The first deals with setting up the TACS
model including reading the mesh, setting elements and design variables.
The second part deals with creating problem instances that are responsible
for setting loads and functions, performing analysis, and gradient computations.

Developers:
    - Dr. G.K.W. Kenway (GKK)
    - Dr. T.R Brooks

History:
    - v. 1.0 pyTACS initial implementation
    - v. 3.0 updated TACS 3.0 pyTACS implementation
"""
# =============================================================================
# Imports
# =============================================================================
from __future__ import print_function

import copy
import numbers
import time
import warnings
from functools import wraps

import numpy as np
import pyNastran.bdf as pn

import tacs.TACS
import tacs.constitutive
import tacs.constraints
import tacs.elements
import tacs.functions
import tacs.problems
from tacs.pymeshloader import pyMeshLoader
from tacs.utilities import BaseUI

warnings.simplefilter("default")


# Define decorator functions for methods that must be called before initialize
def preinitialize_method(method):
    @wraps(method)
    def wrapped_method(self, *args, **kwargs):
        if self.assembler is not None:
            raise self._TACSError(
                f"`{method.__name__}` is a pre-initialize method. "
                "It may only be called before the 'initialize' method has been called."
            )
        else:
            return method(self, *args, **kwargs)

    return wrapped_method


# Define decorator functions for methods that must be called after initialize
def postinitialize_method(method):
    @wraps(method)
    def wrapped_method(self, *args, **kwargs):
        if self.assembler is None:
            raise self._TACSError(
                f"`{method.__name__}` is a post-initialize method. "
                "It may only be called after the 'initialize' method has been called."
            )
        else:
            return method(self, *args, **kwargs)

    return wrapped_method


class pyTACS(BaseUI):
    """
    The class for working with a TACS structure
    """

    # Default class options
    defaultOptions = {
        # Meshloader options
        "printDebug": [
            bool,
            False,
            "Flag for whether to print debug information while loading file.",
        ],
        # Output Options
        "outputElement": [
            int,
            None,
            "Specifies which element type should be written out in the f5 file.\n"
            "\t If None, the type will be inferred from the first element in the model.\n"
            "\t Acceptable values are:\n"
            f"\t\t tacs.TACS.ELEMENT_NONE = {tacs.TACS.ELEMENT_NONE}\n"
            f"\t\t tacs.TACS.SCALAR_2D_ELEMENT = {tacs.TACS.SCALAR_2D_ELEMENT}\n"
            f"\t\t tacs.TACS.SCALAR_3D_ELEMENT = {tacs.TACS.SCALAR_3D_ELEMENT}\n"
            f"\t\t tacs.TACS.BEAM_OR_SHELL_ELEMENT = {tacs.TACS.BEAM_OR_SHELL_ELEMENT}\n"
            f"\t\t tacs.TACS.PLANE_STRESS_ELEMENT = {tacs.TACS.PLANE_STRESS_ELEMENT}\n"
            f"\t\t tacs.TACS.SOLID_ELEMENT = {tacs.TACS.SOLID_ELEMENT}\n"
            f"\t\t tacs.TACS.RIGID_ELEMENT = {tacs.TACS.RIGID_ELEMENT}\n"
            f"\t\t tacs.TACS.MASS_ELEMENT = {tacs.TACS.MASS_ELEMENT}\n"
            f"\t\t tacs.TACS.SPRING_ELEMENT = {tacs.TACS.SPRING_ELEMENT}\n"
            f"\t\t tacs.TACS.PCM_ELEMENT = {tacs.TACS.PCM_ELEMENT}",
        ],
        "writeConnectivity": [
            bool,
            True,
            "Flag for whether to include element connectivity in f5 file.",
        ],
        "writeNodes": [bool, True, "Flag for whether to include nodes in f5 file."],
        "writeDisplacements": [
            bool,
            True,
            "Flag for whether to include nodal displacements in f5 file.",
        ],
        "writeStrains": [
            bool,
            True,
            "Flag for whether to include element strains in f5 file.",
        ],
        "writeStresses": [
            bool,
            True,
            "Flag for whether to include element stresses in f5 file.",
        ],
        "writeExtras": [
            bool,
            True,
            "Flag for whether to include element extra variables in f5 file.",
        ],
        "writeLoads": [
            bool,
            True,
            "Flag for whether to include external nodal loads in f5 file.",
        ],
        "writeCoordinateFrame": [
            bool,
            False,
            "Flag for whether to include element coordinate frames in f5 file.",
        ],
        "familySeparator": [
            str,
            "/",
            "Family separator character used for condensing groups in f5 file.",
        ],
        "printTiming": [
            bool,
            False,
            "Flag for printing out timing information for class procedures.",
        ],
        "linearityTol": [
            float,
            1e-14,
            "When created, pyTACS will check if the model is linear or nonlinear by checking whether (res(2*u) - res(0)) - 2 * (res(u) - res(0)) == 0 this tolerance controls how close to zero the residual must be to be considered linear.",
        ],
    }

    def __init__(self, bdf, comm=None, dvNum=0, scaleList=None, options=None):
        """

        Parameters
        ----------
        bdf : str or pyNastran.bdf.bdf.BDF
            The BDF file or a pyNastran BDF object to load.

        comm : mpi4py.MPI.Intracomm
            The comm object on which to create the pyTACS object.

        dvNum : int
            A user-supplied offset to the design variable
            numbering. This is typically used with tacs+tripan when
            geometric variables have already been added and assigned
            global tacs numberings.

        scaleList: list
            when dvNum is non-zero, the scaleList must be same size
            as the number of design variables already added. i.e.
            len(scaleList) = dvNum

        options : dict
            Dictionary holding model-specific option parameters (case-insensitive).
        """

        startTime = time.time()

        # Setup comm and options
        BaseUI.__init__(self, options=options, comm=comm)

        importTime = time.time()

        # Create and load mesh loader object.
        debugFlag = self.getOption("printDebug")
        self.meshLoader = pyMeshLoader(self.comm, debugFlag)
        self.meshLoader.scanBdfFile(bdf)
        # Save pynastran bdf object
        self.bdfInfo = self.meshLoader.getBDFInfo()
        self.bdfName = self.bdfInfo.bdf_filename

        meshLoadTime = time.time()

        # Retrieve the number of components. This is the maximum
        # number of unique constitutive objects possible in this model.
        self.nComp = self.meshLoader.getNumComponents()

        # Load all the component descriptions
        self.compDescripts = self.meshLoader.getComponentDescripts()
        self.elemDescripts = self.meshLoader.getElementDescripts()

        # Set the starting dvNum and scaleList
        self.dvNum = dvNum
        self.scaleList = scaleList
        if scaleList is None:
            self.scaleList = []

        DVPreprocTime = time.time()

        # List of DV groups
        self.globalDVs = {}
        self.massDVs = {}
        self.compIDBounds = {}
        self.addedCompIDs = set()

        # List of initial coordinates
        self.Xpts0 = None
        # List of initial designvars
        self.x0 = None
        # Design var upper/lower-bounds
        self.xub = None
        self.xlb = None

        # Variables per node for model
        self.varsPerNode = None

        # TACS assembler object
        self.assembler = None

        # Nonlinear flag
        self._isNonlinear = None

        initFinishTime = time.time()
        if self.getOption("printTiming"):
            self._pp("+--------------------------------------------------+")
            self._pp("|")
            self._pp("| TACS Init Times:")
            self._pp("|")
            self._pp(
                "| %-30s: %10.3f sec" % ("TACS Module Time", importTime - startTime)
            )
            self._pp(
                "| %-30s: %10.3f sec"
                % ("TACS Meshload Time", meshLoadTime - importTime)
            )
            self._pp(
                "| %-30s: %10.3f sec"
                % ("TACS DV Processing Time", DVPreprocTime - meshLoadTime)
            )
            self._pp(
                "| %-30s: %10.3f sec"
                % ("TACS Finalize Initialization Time", initFinishTime - DVPreprocTime)
            )
            self._pp("|")
            self._pp(
                "| %-30s: %10.3f sec"
                % ("TACS Total Initialization Time", initFinishTime - startTime)
            )
            self._pp("+--------------------------------------------------+")

    @property
    def isNonlinear(self):
        """The public interface for the isNonlinear attribute. Implemented as a property so that it is read-only."""
        return self._isNonlinear

    @preinitialize_method
    def addGlobalDV(self, descript, value, lower=None, upper=None, scale=1.0):
        """
        Add a global design variable that can affect multiple components.

        This function allows adding design variables that are not
        cleanly associated with a particular constitutive object. One
        example is the pitch of the stiffeners for blade-stiffened
        panels. It is often the same for many different constitutive
        objects. By calling this function, the internal dvNum counter
        is incremented, and the user doesn't have to worry about
        it.

        Parameters
        ----------
        descript : str
            A user-supplied string that can be used to retrieve the
            variable number and value elemCallBackFunction.
        value : float
            Initial value for variable.
        lower : float
            Lower bound. This may be None for unbounded
        upper : float
            Upper bound. This may be None for unbounded
        scale : float
            Scale factor for variable
        """
        self.globalDVs[descript] = {
            "num": self.dvNum,
            "value": value,
            "lowerBound": lower,
            "upperBound": upper,
            "isMassDV": False,
        }
        self.dvNum += 1
        self.scaleList.append(scale)

    def getGlobalDVs(self):
        """
        Return dict holding info about all current global DVs.

        Returns
        -------
        globalDVs : dict
            Dictionary holding global dv information.
        """
        return self.globalDVs.copy()

    def getGlobalDVKeys(self):
        """
        Get key names for all current global DVs.

        Returns
        -------
        globalDVKeys : list[str]
            List holding global dv names.
        """
        return list(self.globalDVs.keys())

    def getGlobalDVNums(self):
        """
        Get the dv nums corresponding to global DVs.

        Returns
        -------
        globalDVNums : list[int]
            List holding dv nums corresponding to global DVs.
        """
        return [self.globalDVs[descript]["num"] for descript in self.globalDVs]

    def getTotalNumGlobalDVs(self):
        """
        Get the total number of global DVs across all processors.

        Returns
        -------
        globalDVs : dict
            Dictionary holding global dv information.
        """
        return len(self.globalDVs)

    @preinitialize_method
    def assignMassDV(self, descript, eIDs, dvName="m"):
        """
        Assign a global DV to a point mass element.

        Parameters
        ----------
        descript : str
            Global DV key to assign mass design variable to. If the key is does not exist,
            it will automatically be created and added to global DVs.

        eIDs : int or list[int]
            Element IDs of concentrated mass to assign DV to (NASTRAN ordering)

        dvName : str
            Name of mass property to apply DV to.
            May be `m` for mass, `I11`, `I22`, `I12`, etc. for moment of inertia components.
            Defaults to `m` (mass).

        Notes
        -----
        Currently only CONM2 cards are supported.
        """
        # Make sure eID is an array
        eIDs = np.atleast_1d(eIDs)

        # Check if referenced element ID is a CONM2 element
        for eID in eIDs:
            is_mass_element = False
            if eID in self.bdfInfo.masses:
                if self.bdfInfo.masses[eID].type in ["CONM2"]:
                    is_mass_element = True

            if not is_mass_element:
                raise self._TACSError(
                    f"Element ID '{eID}' does not correspond to a `CONM2` element. "
                    "Only `CONM2` elements are supported for this method."
                )

        # Check if descript already exists in global dvs, if not add it
        if descript not in self.globalDVs:
            self.addGlobalDV(descript, None)

        dv_dict = self.globalDVs[descript]

        # Flag this global dv as being a mass dv
        dv_dict["isMassDV"] = True

        massDV = dv_dict["num"]
        value = dv_dict["value"]
        ub = dv_dict["upperBound"]
        lb = dv_dict["lowerBound"]

        for eID in eIDs:
            # If the element ID hasn't already been added to massDVs, add it
            if eID not in self.massDVs:
                self.massDVs[eID] = {}

            # Update the element entry with the dv num
            self.massDVs[eID][f"{dvName}Num"] = massDV

            # Update the element entry with the dv name
            if value is not None:
                self.massDVs[eID][dvName] = value
            # If value was defined from previous call, remove it
            elif dvName in self.massDVs[eID]:
                self.massDVs[eID].pop(dvName)

            # Update the element entry with the dv upper bound
            if ub is not None:
                self.massDVs[eID][f"{dvName}ub"] = ub
            # If upper bound was defined from previous call, remove it
            elif f"{dvName}ub" in self.massDVs[eID]:
                self.massDVs[eID].pop(f"{dvName}ub")

            # Update the element entry with the dv lower bound
            if lb is not None:
                self.massDVs[eID][f"{dvName}lb"] = lb
            # If lower bound was defined from previous call, remove it
            elif f"{dvName}lb" in self.massDVs[eID]:
                self.massDVs[eID].pop(f"{dvName}lb")

    def selectCompIDs(
        self,
        include=None,
        exclude=None,
        includeBounds=None,
        nGroup=1,
        includeOp="or",
        excludeOp="or",
        projectVector=None,
        **kwargs,
    ):
        """
        This is the most important function of the entire setup
        process.
        The basic idea is as follows: We have a list of nComp
        which are the component descriptions.
        What we need is a way of
        generating subgroups of these for the purposes of adding
        design variables, constitutive objects, KS domains, and mass
        domains.
        All of these operations boil down to selecting a
        subset of the compIDs.

        This function attempts to support as many ways as possible to
        select parts of the structure.
        Easy and efficient selection of
        parts is critical to the end user.

        Methods of selection:

        1. include, integer, string, list of integers and/or strings: The
        simplest and most direct way of selecting a component.
        The
        user supplies the index of the componentID, a name or partial
        name, or a list containing a combination of both.

        For example::

            # Select the 11th component
            selectCompIDs(include=10)

            # Select the first and fifth component
            selectCompIDs(include=[0, 4])

            # Select any component containing 'rib.00'
            selectCompIDs(include='rib.00')

            # Select any components containing 'rib.00' and 'rib.10'
            selectCompIDs(include=['rib.00', 'rib.10'])

            # Select any component containing 'rib.00', the 11th
            # component and any component containing 'spar'
            # (This is probably not advisable!)
            selectCompIDs(include=['rib.00', 10, 'spar'])

        2. Exclude, operates similarly to 'include'.
        The behaviour of exclude is identical to include above, except that
        component ID's that are found using 'exclude' are
        'subtracted' from those found using include.
        A special case is treated if 'include' is NOT given: if only an
        exclude list is given, this implies the selection of all
        compID's EXCEPT the those in exclude.

        For example::

            # This will return will [0, 1, 2, 3, 5, ..., nComp-1]
            selectCompIDs(exclude = 4)

            # This will return [0, 1, 4, 5, ..., nComp-1]
            selectCompIDs(exclude = [2, 3]) will return

            # This will return components that have 'ribs' in the
            # component ID, but not those that have 'le_ribs' in the
            # component id.
            selectCompIDs(include='ribs', exclude='le_ribs')

        3. includeBounds, list of components defining a region inside
        which 'include' components will be selected.
        This functionality uses a geometric approach to select the compIDs.
        All components within the project 2D convex hull are included.
        Therefore, it is essential to split up concave include regions
        into smaller convex regions.
        Use multiple calls to selectCompIDs to accumulate multiple regions.

        For example::

            # This will select upper skin components between the
            # leading and trailing edge spars and between ribs 1 and 4.
            selectCompIDs(include='U_SKIN', includeBound=
                ['LE_SPAR', 'TE_SPAR', 'RIB.01', 'RIB.04'])

        4. nGroup: The number of groups to divide the found components
        into.
        Generally this will be 1. However, in certain cases, it
        is convenient to create multiple groups in one pass.

        For example::

            # This will 'evenly' create 10 groups on all components
            # containing LE_SPAR.
            Note that once the components are
            # selected, they are sorted **alphabetically** and assigned
            # sequentially.
            selectCompIDs(include='LE_SPAR', nGroup=10)

        nGroup can also be negative.
        If it is negative, then a single
        design variable group is added to each of the found
        components.

        For example::

            # will select all components and assign a design variable
            # group to each one.
            selectCompIDs(nGroup=-1)

        includeOp, str: 'and' or 'or'.
        Selects the logical operation
        used for item in 'include' option.
        For example:

        selectCompIDs(include=['LE_SPAR', 'TE_SPAR'],
        includeOpt='or') will select the LE_SPAR and TE_SPAR
        components (default behaviour).

        selectCompIDs(include=['RIB', 'SEG.01'], includeOpt='and')
        will select any component with 'RIB' in the description AND
        'SEG.01' in the description.
        """

        # Defaults
        includeIDs = np.arange(self.nComp)
        excludeIDs = []
        includeBoundIDs = None

        if include is not None:
            includeIDs = self._getCompIDs(includeOp, include)

        if exclude is not None:
            excludeIDs = self._getCompIDs(excludeOp, exclude)

        iSet = set(includeIDs)
        eSet = set(excludeIDs)

        # First take the intersection of iSet and ibSet
        if includeBoundIDs is not None:
            tmp = iSet.intersection(set(includeBoundIDs))
        else:
            tmp = iSet

        # Next take the difference between tmp and eSet
        compIDs = tmp.difference(eSet)

        # Convert back to a list:
        compIDs = list(compIDs)

        # If we only want a single group, we're done, otherwise, we
        # have a bit more work to do...
        if nGroup > 1:
            # The user wants to have nGroups returned from compIDs.

            # First check that nGroup <= len(compIDs), print warning
            # and clip if not
            if nGroup > len(compIDs):
                self._TACSWarning(
                    f"nGroup={nGroup} is larger than the number of\
                selected components={len(compIDs)}. nGroup will be clipped to {nGroup}"
                )
                nGroup = len(compIDs)

            # Pluck out the component descriptions again and we will
            # sort them
            compDescript = []
            for i in range(len(compIDs)):
                compDescript.append(self.compDescripts[compIDs[i]])

            # define a general argsort
            def argsort(seq):
                return sorted(range(len(seq)), key=seq.__getitem__)

            # ind is the index that would result in a sorted list.
            ind = argsort(compDescript)

            # Now simply divide 'ind' into 'nGroups' as evenly as
            # possible, in the integer sense.
            def split_list(alist, wanted_parts=1):
                length = len(alist)
                return [
                    alist[i * length // wanted_parts : (i + 1) * length // wanted_parts]
                    for i in range(wanted_parts)
                ]

            ind = split_list(ind, nGroup)

            # Finally assemble the nested list of component IDs
            tmp = []
            for i in range(len(ind)):
                tmp.append([])
                for j in range(len(ind[i])):
                    tmp[-1].append(compIDs[ind[i][j]])
            compIDs = tmp
        elif nGroup < 0:
            # Negative number signifies 'add one dv to each component'
            tmp = []
            for comp in compIDs:
                tmp.append([comp])
            compIDs = tmp
        else:
            # Otherwise, just put the current list of compIDs in a
            # list of length 1.
            compIDs = [compIDs]

        return compIDs

    def getBDFInfo(self):
        """
        Return a pynastran bdf object.
        This object can be used interactively
        to parse information (nodes, elements, loads, etc.) included in the bdf file.

        Returns
        -------
        bdfInfo : pyNastran.bdf.bdf.BDF
            pyNastran bdf object.
        """
        return self.bdfInfo

    def getCompNames(self, compIDs=None):
        """
        Return a list of component descriptions for the given component
        IDs. compIDs should come from a call to selectCompIDs

        Parameters
        ----------
        compIDs : int or list[int] or None
            List of integers containing the compIDs numbers. If None, returns names for all components.
            Defaults to None.

        Returns
        -------
        compDescript : list[str]
            List of strings containing the names of the corresponding compIDs
        """
        # Return all component names
        if compIDs is None:
            return copy.deepcopy(self.compDescripts)
        # Convert to list
        elif isinstance(compIDs, (int, np.integer)):
            compIDs = [compIDs]
        # Make sure list is flat
        else:
            compIDs = self._flatten(compIDs)

        compDescripts = []
        for i in range(len(compIDs)):
            compDescripts.append(self.compDescripts[compIDs[i]])

        return compDescripts

    def getGlobalNodeIDsForComps(self, compIDs, nastranOrdering=False):
        """
        Return the global (non-partitioned) node IDs belonging to a given list of component IDs

        Parameters
        ----------
        compIDs : int or list[int] or None
            List of integers containing the compIDs numbers.
            If None, returns nodeIDs for all components.
            Defaults to None.

        nastranOrdering : bool
            Flag signaling whether nodeIDs are in TACS (default) or NASTRAN (grid IDs in bdf file) ordering
            Defaults to False.

        Returns
        -------
        nodeIDs : list[int]
            List of unique nodeIDs that belong to the given list of compIDs
        """
        # Return all component ids
        if compIDs is None:
            compIDs = list(range(self.nComp))

        return self.meshLoader.getGlobalNodeIDsForComps(compIDs, nastranOrdering)

    @postinitialize_method
    def getLocalNodeIDsForComps(self, compIDs):
        """
        Return the local (partitioned) node IDs belonging to a given list of component IDs

        Parameters
        ----------
         compIDs : int or list[int] or None
            List of integers containing the compIDs numbers.
            If None, returns nodeIDs for all components.
            Defaults to None.

        Returns
        -------
        nodeIDs : list[int]
            List of unique nodeIDs that belong to the given list of compIDs
        """
        # Return all component ids
        if compIDs is None:
            compIDs = list(range(self.nComp))

        return self.meshLoader.getLocalNodeIDsForComps(compIDs)

    def getLocalNodeIDsFromGlobal(self, globalIDs, nastranOrdering=False):
        """
        Given a list of node IDs in global (non-partitioned) ordering
        returns the local (partitioned) node IDs on each processor.
        If a requested node is not included on this processor,
        an entry of -1 will be returned.

        Parameters
        ----------
        globalIDs : int or list[int]
            List of global node IDs.

        nastranOrdering : bool
            Flag signaling whether globalIDs is in TACS (default) or NASTRAN (grid IDs in bdf file) ordering
            Defaults to False.

        Returns
        -------
        localIDs : list[int]
            List of local node IDs for each entry in globalIDs.
            If the node is not owned by this processor, its index is filled with a value of -1.
        """

        return self.meshLoader.getLocalNodeIDsFromGlobal(globalIDs, nastranOrdering)

    def initialize(self, elemCallBack=None):
        """
        This is the 'last' method to be called during the setup. The
        user should have already added all the design variables,
        domains, etc. Before this function is called. This function
        finalizes the problem initialization and cannot be changed at
        later time. If the user does not provide an elemCallBack function,
        we will use pyNastran to generate one automatically from element
        properties provided in the BDF file.

        Parameters
        ----------
        elemCallBack : collections.abc.Callable or None

           The calling sequence for elemCallBack **must** be as
           follows::

             def elemCallBack(dvNum, compID, compDescript, elemDescripts,
                             globalDVs, **kwargs):

           The dvNum is the current counter which must be used by the
           user when creating a constitutive object with design
           variables.

           compID is the ID number used by tacs to reference this property group.
           Use kwargs['propID'] to get the corresponding Nastran property ID that
           is read in from the BDF.

           compDescript is the component description label read in from optional
           formatted comments in BDF file

           elemDescripts are the name of the elements belonging to this group
           (e.g. CQUAD4, CTRIA3, CTETRA, etc). This value will be a list since
           one component may contain multiple compatible element types.
           Example: ['CQUAD4', CTRIA3']

           globalDVs is a dictionary containing information about any
           global DVs that have been added.

           elemCallBack must return a list containing as many TACS element
           objects as there are element types in elemDescripts (one for each).

        """

        if elemCallBack is None:
            elemCallBack = self._elemCallBackFromBDF()
        self._createOutputGroups()
        self._createElements(elemCallBack)

        self.assembler = self.meshLoader.createTACSAssembler(
            self.varsPerNode, self.massDVs
        )
        self._createOutputViewer()

        # Store original node locations read in from bdf file
        self.Xpts0 = self.assembler.createNodeVec()
        self.assembler.getNodes(self.Xpts0)

        # Store initial design variable values
        self.x0 = self.assembler.createDesignVec()
        self.assembler.getDesignVars(self.x0)

        # Store design variable upper/lower-bounds
        self.xub = self.assembler.createDesignVec()
        self.xlb = self.assembler.createDesignVec()
        self.assembler.getDesignVarRange(self.xlb, self.xub)

        self._isNonlinear = self._checkNonlinearity()

    @postinitialize_method
    def _checkNonlinearity(self) -> bool:
        """Check if the finite element model is nonlinear

        This check works by checking whether the residual is nonlinear w.r.t the states using 3 residual evaluations.

        Returns
        -------
        bool
            True if the problem is nonlinear, False otherwise.
        """
        res0 = self.assembler.createVec()
        res1 = self.assembler.createVec()
        res2 = self.assembler.createVec()
        state = self.assembler.createVec()

        # Evaluate r(0)
        state.zeroEntries()
        self.assembler.setVariables(state, state, state)
        self.assembler.assembleRes(res0)

        # Evaluate r(u) - r(0)
        state.initRand()
        state.setRand()
        self.setBCsInVec(state)
        self.assembler.setVariables(state, state, state)
        self.assembler.assembleRes(res1)
        res1.axpy(-1.0, res0)

        # Evaluate r(2u) -  r(0)
        state.scale(2.0)
        self.setBCsInVec(state)
        self.assembler.setVariables(state, state, state)
        self.assembler.assembleRes(res2)
        res2.axpy(-1.0, res0)

        # Reset the state variables
        state.zeroEntries()
        self.assembler.setVariables(state, state, state)

        # Check if (res2-res0) - 2 * (res1 - res0) is zero (or very close to it)
        resNorm = np.real(res1.norm())
        res2.axpy(-2.0, res1)
        if resNorm == 0.0 or (np.real(res2.norm()) / resNorm) <= self.getOption(
            "linearityTol"
        ):
            return False  # not nonlinear case
        else:
            return True  # nonlinear case

    def _elemCallBackFromBDF(self):
        """
        Automatically setup elemCallBack using information contained in BDF file.
        This function assumes all material properties are specified in the BDF.
        """

        # Check if any properties are in the BDF
        if self.bdfInfo.missing_properties:
            raise self._TACSError(
                f"BDF file '{self.bdfName}' has missing properties cards. "
                "Set 'printDebug' option to True for more information. "
                "User must define own elemCallBack function."
            )

        # Make sure cross-referencing is turned on in pynastran
        if self.bdfInfo.is_xrefed is False:
            self.bdfInfo.cross_reference()
            self.bdfInfo.is_xrefed = True

        # Create a dictionary to sort all elements by property number
        elemDict = {}
        for elementID in self.bdfInfo.elements:
            element = self.bdfInfo.elements[elementID]
            propertyID = element.pid
            if propertyID not in elemDict:
                elemDict[propertyID] = {}
                elemDict[propertyID]["elements"] = []
                elemDict[propertyID]["dvs"] = {}
            elemDict[propertyID]["elements"].append(element)

        # Create a dictionary to sort all design variables
        for dv in self.bdfInfo.dvprels:
            propertyID = self.bdfInfo.dvprels[dv].pid
            dvName = self.bdfInfo.dvprels[dv].pname_fid
            self.dvNum = max(self.dvNum, self.bdfInfo.dvprels[dv].dvids[0])
            elemDict[propertyID]["dvs"][dvName] = self.bdfInfo.dvprels[dv]
        # Create option for user to specify scale values in BDF
        self.scaleList = [1.0] * self.dvNum

        # Callback function to return appropriate tacs MaterialProperties object
        # For a pynastran mat card
        def matCallBack(matInfo):
            # Nastran isotropic material card
            if matInfo.type == "MAT1":
                mat = tacs.constitutive.MaterialProperties(
                    rho=matInfo.rho,
                    E=matInfo.e,
                    nu=matInfo.nu,
                    ys=matInfo.St,
                    alpha=matInfo.a,
                )
            # Nastran orthotropic material card
            elif matInfo.type == "MAT8":
                E1 = matInfo.e11
                E2 = matInfo.e22
                nu12 = matInfo.nu12
                G12 = matInfo.g12
                G13 = matInfo.g1z
                G23 = matInfo.g2z
                # If out-of-plane shear values are 0, Nastran defaults them to the in-plane
                if G13 == 0.0:
                    G13 = G12
                if G23 == 0.0:
                    G23 = G12
                rho = matInfo.rho
                Xt = matInfo.Xt
                Xc = matInfo.Xc
                Yt = matInfo.Yt
                Yc = matInfo.Yc
                S12 = matInfo.S

                if (
                    S12 == 0 or Xt == 0 or Xc == 0 or Yt == 0 or Yc == 0
                ):
                    self._TACSWarning(
                        f"MAT8 card {matInfo.mid} has a zero strength, check Xc, Xt, Yc, Yt, and S12."
                        "Otherwise Tsai-Wu Failure criterion is undefined or infinity."
                    )

                # TODO: add alpha
                mat = tacs.constitutive.MaterialProperties(
                    rho=rho,
                    E1=E1,
                    E2=E2,
                    nu12=nu12,
                    G12=G12,
                    G13=G13,
                    G23=G23,
                    Xt=Xt,
                    Xc=Xc,
                    Yt=Yt,
                    Yc=Yc,
                    S12=S12,
                )
            # Nastran 2D anisotropic material card
            elif matInfo.type == "MAT2":
                C11 = matInfo.G11
                C12 = matInfo.G12
                C22 = matInfo.G22
                C13 = matInfo.G13
                C23 = matInfo.G23
                C33 = matInfo.G33
                rho = matInfo.rho
                # See if this card features anisotropic coupling terms (which we don't support yet)
                if (
                    np.abs(C13) / (C11 + C22) >= 1e-8
                    or np.abs(C23) / (C11 + C22) >= 1e-8
                ):
                    self._TACSWarning(
                        f"MAT2 card {matInfo.mid} has anisotropic stiffness components that are not currently supported. "
                        "These terms will be dropped and the material treated as orthotropic. "
                        "Result accuracy may be affected."
                    )
                nu12 = C12 / C22
                nu21 = C12 / C11
                E1 = C11 * (1 - nu12 * nu21)
                E2 = C22 * (1 - nu12 * nu21)
                G12 = G13 = G23 = C33
                # TODO: add alpha
                mat = tacs.constitutive.MaterialProperties(
                    rho=rho, E1=E1, E2=E2, nu12=nu12, G12=G12, G13=G13, G23=G23
                )

            else:
                raise self._TACSError(
                    f"Unsupported material type '{matInfo.type}' for material number {matInfo.mid}."
                )

            return mat

        def elemCallBack(
            dvNum, compID, compDescript, elemDescripts, globalDVs, **kwargs
        ):
            # Initialize scale list for design variables we will add
            scaleList = []

            # Get the Nastran property ID
            propertyID = kwargs["propID"]
            propInfo = self.bdfInfo.properties[propertyID]
            elemInfo = elemDict[propertyID]["elements"][0]

            # First we define the material object
            mat = None
            # This property only references one material
            if hasattr(propInfo, "mid_ref"):
                matInfo = propInfo.mid_ref
                mat = matCallBack(matInfo)
            # This property references multiple materials (maybe a laminate)
            elif hasattr(propInfo, "mids_ref"):
                mat = []
                for matInfo in propInfo.mids_ref:
                    mat.append(matCallBack(matInfo))

            # Next we define the constitutive object
            if propInfo.type == "PSHELL":  # Nastran isotropic shell
                kcorr = propInfo.tst

                if "T" in elemDict[propertyID]["dvs"]:
                    thickness = elemDict[propertyID]["dvs"]["T"].dvids_ref[0].xinit
                    tNum = elemDict[propertyID]["dvs"]["T"].dvids[0] - 1
                    minThickness = elemDict[propertyID]["dvs"]["T"].dvids_ref[0].xlb
                    maxThickness = elemDict[propertyID]["dvs"]["T"].dvids_ref[0].xub
                    name = elemDict[propertyID]["dvs"]["T"].dvids_ref[0].label
                    self.scaleList[tNum - 1] = elemDict[propertyID]["dvs"]["T"].coeffs[
                        0
                    ]
                else:
                    thickness = propInfo.t
                    tNum = -1
                    minThickness = 0.0
                    maxThickness = 1e20

                con = tacs.constitutive.IsoShellConstitutive(
                    mat, t=thickness, tlb=minThickness, tub=maxThickness, tNum=tNum
                )

            elif propInfo.type == "PCOMP":  # Nastran composite shell
                numPlies = propInfo.nplies
                plyThicknesses = []
                plyAngles = []
                plyMats = []

                # if the laminate is symmetric, mirror the ply indices
                if propInfo.lam == "SYM":
                    plyIndices = list(range(numPlies // 2))
                    plyIndices.extend(plyIndices[::-1])
                else:
                    plyIndices = range(numPlies)

                # Loop through plies and setup each entry in layup
                for ply_i in plyIndices:
                    plyThicknesses.append(propInfo.thicknesses[ply_i])
                    plyMat = tacs.constitutive.OrthotropicPly(
                        plyThicknesses[ply_i], mat[ply_i]
                    )
                    plyMats.append(plyMat)
                    plyAngles.append(np.deg2rad(propInfo.thetas[ply_i]))

                # Convert thickness/angles to appropriate numpy array
                plyThicknesses = np.array(plyThicknesses, dtype=self.dtype)
                plyAngles = np.array(plyAngles, dtype=self.dtype)

                # Get the total laminate thickness
                lamThickness = propInfo.Thickness()
                # Get the offset distance from the ref plane to the midplane
                tOffset = -(propInfo.z0 / lamThickness + 0.5)

                if propInfo.lam is None or propInfo.lam in ["SYM", "MEM"]:
                    # Discrete laminate class (not for optimization)
                    con = tacs.constitutive.CompositeShellConstitutive(
                        plyMats, plyThicknesses, plyAngles, tOffset=tOffset
                    )

                elif propInfo.lam == "SMEAR":
                    plyFractions = plyThicknesses / lamThickness
                    con = tacs.constitutive.SmearedCompositeShellConstitutive(
                        plyMats, lamThickness, plyAngles, plyFractions, t_offset=tOffset
                    )

                # Need to add functionality to consider only membrane in TACS for type = MEM
                else:
                    raise self._TACSError(
                        f"Unrecognized LAM type '{propInfo.lam}' for PCOMP number {propertyID}."
                    )

            elif propInfo.type == "PSOLID":  # Nastran solid property
                if "T" in elemDict[propertyID]["dvs"]:
                    thickness = elemDict[propertyID]["dvs"]["T"].dvids_ref[0].xinit
                    tNum = elemDict[propertyID]["dvs"]["T"].dvids[0] - 1
                    minThickness = elemDict[propertyID]["dvs"]["T"].dvids_ref[0].xlb
                    maxThickness = elemDict[propertyID]["dvs"]["T"].dvids_ref[0].xub
                    name = elemDict[propertyID]["dvs"]["T"].dvids_ref[0].label
                    self.scaleList[tNum - 1] = elemDict[propertyID]["dvs"]["T"].coeffs[
                        0
                    ]
                else:
                    thickness = 1.0
                    tNum = -1
                    minThickness = 0.0
                    maxThickness = 10.0

                con = tacs.constitutive.SolidConstitutive(
                    mat, t=thickness, tlb=minThickness, tub=maxThickness, tNum=tNum
                )

            elif propInfo.type == "PBUSH":  # Nastran spring
                k = np.zeros(6)
                for j in range(len(k)):
                    if propInfo.Ki[j]:
                        k[j] = propInfo.Ki[j]
                con = tacs.constitutive.DOFSpringConstitutive(k=k)

            elif propInfo.type == "PBAR":  # Nastran bar
                area = propInfo.A
                I1 = propInfo.i1
                I2 = propInfo.i2
                # Nastran uses negative convention for POI's
                I12 = -propInfo.i12
                J = propInfo.j
                k1 = propInfo.k1
                k2 = propInfo.k2

                # pynastran defaults these values to 1e8,
                # which can lead to scaling issues in the stiffness matrix
                # We truncate this value to 1e3 to prevent this
                if k1 is None or k1 > 1e3:
                    k1 = 1e3

                if k2 is None or k2 > 1e3:
                    k2 = 1e3

                con = tacs.constitutive.BasicBeamConstitutive(
                    mat, A=area, Iy=I2, Iz=I1, Iyz=I12, J=J, ky=k1, kz=k2
                )

            elif propInfo.type == "PROD":  # Nastran rod
                area = propInfo.A
                J = propInfo.j
                k1 = 0.0
                k2 = 0.0

                con = tacs.constitutive.BasicBeamConstitutive(
                    mat, A=area, J=J, ky=k1, kz=k2
                )

            else:
                raise self._TACSError(
                    f"Unsupported property type '{propInfo.type}' for property number {propertyID}. "
                )

            # Set up transform object which may be required for certain elements
            transform = None
            if propInfo.type in ["PSHELL", "PCOMP"]:
                mcid = elemDict[propertyID]["elements"][0].theta_mcid_ref
                if mcid:
                    if mcid.type == "CORD2R":
                        refAxis = mcid.i
                        transform = tacs.elements.ShellRefAxisTransform(refAxis)
                    else:  # Don't support spherical/cylindrical yet
                        raise self._TACSError(
                            "Unsupported material coordinate system type "
                            f"'{mcid.type}' for property number {propertyID}."
                        )
            elif propInfo.type in ["PBAR"]:
                refAxis = elemDict[propertyID]["elements"][0].g0_vector
                transform = tacs.elements.BeamRefAxisTransform(refAxis)
            elif propInfo.type == "PROD":
                refAxis = np.array(
                    [1.0, -1.0, 1.0]
                )  # dummy ref_axis, not really needed for rods
                transform = tacs.elements.BeamRefAxisTransform(refAxis)
            elif propInfo.type == "PBUSH":
                if elemDict[propertyID]["elements"][0].cid_ref:
                    refAxis_i = elemDict[propertyID]["elements"][0].cid_ref.i
                    refAxis_j = elemDict[propertyID]["elements"][0].cid_ref.j
                    transform = tacs.elements.SpringRefFrameTransform(
                        refAxis_i, refAxis_j
                    )
                elif elemDict[propertyID]["elements"][0].x[0]:
                    refAxis = (
                        np.array(elemDict[propertyID]["elements"][0].x)
                        - elemDict[propertyID]["elements"][0]
                        .nodes_ref[0]
                        .get_position()
                    )
                    transform = tacs.elements.SpringRefAxisTransform(refAxis)
                elif elemDict[propertyID]["elements"][0].g0_ref:
                    refAxis = (
                        elemDict[propertyID]["elements"][0].g0_ref.get_position()
                        - elemDict[propertyID]["elements"][0]
                        .nodes_ref[0]
                        .get_position()
                    )
                    transform = tacs.elements.SpringRefAxisTransform(refAxis)

            # Finally set up the element objects belonging to this component
            elemList = []
            for descript in elemDescripts:
                if descript in ["CQUAD4", "CQUADR"]:
                    elem = tacs.elements.Quad4Shell(transform, con)
                elif descript in ["CQUAD9", "CQUAD"]:
                    elem = tacs.elements.Quad9Shell(transform, con)
                elif descript in ["CTRIA3", "CTRIAR"]:
                    elem = tacs.elements.Tri3Shell(transform, con)
                elif descript in ["CBAR", "CROD"]:
                    elem = tacs.elements.Beam2(transform, con)
                elif "CTETRA" in descript:
                    # May have variable number of nodes in card
                    nnodes = len(elemInfo.nodes)
                    if nnodes == 4:
                        basis = tacs.elements.LinearTetrahedralBasis()
                    elif nnodes == 10:
                        basis = tacs.elements.QuadraticTetrahedralBasis()
                    else:
                        raise self._TACSError(
                            f"TACS does not currently support CTETRA elements with {nnodes} nodes."
                        )
                    model = tacs.elements.LinearElasticity3D(con)
                    elem = tacs.elements.Element3D(model, basis)
                elif descript in ["CHEXA8", "CHEXA"]:
                    basis = tacs.elements.LinearHexaBasis()
                    model = tacs.elements.LinearElasticity3D(con)
                    elem = tacs.elements.Element3D(model, basis)
                elif descript == "CBUSH":
                    elem = tacs.elements.SpringElement(transform, con)
                else:
                    raise self._TACSError(
                        "Unsupported element type "
                        f"'{descript}' specified for property number {propertyID}."
                    )
                elemList.append(elem)

            return elemList, scaleList

        return elemCallBack

    @postinitialize_method
    def getOrigDesignVars(self):
        """
        get the original design variables that were specified with
        during assembler creation.

        Returns
        -------
        x : numpy.ndarray
            The original design variable vector set in tacs.

        """
        return self.x0.getArray().copy()

    @postinitialize_method
    def getDesignVarRange(self):
        """
        get the lower/upper bounds for the design variables.

        Returns
        -------
        xlb : numpy.ndarray
            The design variable lower bound.
        xub : numpy.ndarray
            The design variable upper bound.

        """
        return self.xlb.getArray().copy(), self.xub.getArray().copy()

    @postinitialize_method
    def createDesignVec(self, asBVec=False):
        """
        Create a new tacs distributed design vector.
        Values are initialized to zero.

        Parameters
        ----------
        asBVec : bool
            Flag that determines whether to return
            design vector as tacs :class:`~TACS.Vec` (True) or numpy array (False).
            Defaults to False.

        Returns
        -------
        x : numpy.ndarray or tacs.TACS.Vec
            Distributed design variable vector
        """
        xVec = self.assembler.createDesignVec()
        if asBVec:
            return xVec
        else:
            return xVec.getArray()

    @postinitialize_method
    def getNumDesignVars(self):
        """
        Return the number of design variables on this processor.

        Returns
        -------
        ndvs : int
            Number of design variables on this processor.
        """
        return self.x0.getSize()

    @postinitialize_method
    def getTotalNumDesignVars(self):
        """
        Return the number of design variables across all processors.

        Returns
        -------
        ndvs : int
            Total number of design variables across all processors.
        """
        return self.dvNum

    @postinitialize_method
    def getOrigNodes(self):
        """
        Return the original mesh coordinates read in from the meshLoader.

        Returns
        -------
        coords : numpy.ndarray
            Structural coordinate in array of size (N * 3) where N is
            the number of structural nodes on this processor.
        """
        return self.Xpts0.getArray().copy()

    @postinitialize_method
    def createNodeVec(self, asBVec=False):
        """
        Create a new tacs distributed node vector.
        Values are initialized to zero.

        Parameters
        ----------
        asBVec : bool
            Flag that determines whether to return
            node vector as tacs :class:`~TACS.Vec` (True) or numpy array (False).
            Defaults to False.

        Returns
        -------
        xpts : numpy.ndarray or tacs.TACS.Vec
            Distributed node coordinate vector
        """

        xptVec = self.assembler.createNodeVec()
        if asBVec:
            return xptVec
        else:
            return xptVec.getArray()

    @postinitialize_method
    def getNumOwnedNodes(self):
        """
        Get the number of nodes owned by this processor.

        Returns
        -------
        nNodes : int
            Number of nodes owned by this proc.
        """
        return self.assembler.getNumOwnedNodes()

    @postinitialize_method
    def getNumOwnedMultiplierNodes(self):
        """
        Get the number of lagrange multiplier nodes owned by this processor.

        Returns
        -------
        nMultNodes : int
            Number of multiplier nodes owned by this proc.
        """
        return len(self.meshLoader.getLocalMultiplierNodeIDs())

    @postinitialize_method
    def getLocalMultiplierNodeIDs(self):
        """
        Get the tacs indices of multiplier nodes used to hold lagrange multipliers on this processor.

        Returns
        -------
        nodeIDs : list[int]
            List of multiplier node ID's owned by this proc.
        """
        return self.meshLoader.getLocalMultiplierNodeIDs()

    @postinitialize_method
    def createVec(self, asBVec=False):
        """
        Create a new tacs distributed state variable vector.
        Values are initialized to zero.

        Parameters
        ----------
        asBVec : bool
            Flag that determines whether to return
            state vector as tacs :class:`~TACS.Vec` (True) or numpy array (False).
            Defaults to False.

        Returns
        -------
        vars : numpy.ndarray or tacs.TACS.Vec
            Distributed state variable vector
        """
        vars = self.assembler.createVec()
        if asBVec:
            return vars
        else:
            return vars.getArray()

    @postinitialize_method
    def getVarsPerNode(self):
        """
        Get the number of variables per node for the model.

        Returns
        -------
        vpn : int
            Number of variables per node.
        """
        return self.assembler.getVarsPerNode()

    @postinitialize_method
    def applyBCsToVec(self, vec):
        """
        Applies zeros to boundary condition DOFs in input vector.

        Parameters
        ----------
        vec : numpy.ndarray or tacs.TACS.Vec
            Vector to apply boundary conditions to.
        """
        # Check if input is a BVec or numpy array
        if isinstance(vec, tacs.TACS.Vec):
            self.assembler.applyBCs(vec)
        elif isinstance(vec, np.ndarray):
            array = vec
            # Create temporary BVec
            vec = self.assembler.createVec()
            # Copy array values to BVec
            vec.getArray()[:] = array
            # Apply BCs
            self.assembler.applyBCs(vec)
            # Copy values back to array
            array[:] = vec.getArray()

    @postinitialize_method
    def setBCsInVec(self, vec):
        """
        Sets dirichlet boundary condition values in the input vector.

        Parameters
        ----------
        vec : numpy.ndarray or tacs.TACS.Vec
            Vector to set boundary conditions in.
        """
        # Check if input is a BVec or numpy array
        if isinstance(vec, tacs.TACS.Vec):
            self.assembler.setBCs(vec)
        elif isinstance(vec, np.ndarray):
            array = vec
            # Create temporary BVec
            vec = self.assembler.createVec()
            # Copy array values to BVec
            vec.getArray()[:] = array
            # Apply BCs
            self.assembler.setBCs(vec)
            # Copy values back to array
            array[:] = vec.getArray()

    @postinitialize_method
    def createStaticProblem(self, name, options=None):
        """
        Create a new staticProblem for modeling a static load cases.
        This object can be used to set loads, evalFunctions as well as perform
        solutions and sensitivities related to static problems

        Parameters
        ----------
        name : str
            Name to assign problem.
        options : dict
            Problem-specific options to pass to StaticProblem instance (case-insensitive).
            Defaults to None.

        Returns
        -------
        problem : tacs.problems.StaticProblem
            StaticProblem object used for modeling and solving static cases.
        """
        problem = tacs.problems.static.StaticProblem(
            name,
            self.assembler,
            self.comm,
            self.outputViewer,
            self.meshLoader,
            self.isNonlinear,
            options,
        )
        # Set with original design vars and coordinates, in case they have changed
        problem.setDesignVars(self.x0)
        problem.setNodes(self.Xpts0)
        return problem

    @postinitialize_method
    def createTransientProblem(self, name, tInit, tFinal, numSteps, options=None):
        """
        Create a new TransientProblem for modeling a transient load cases.
        This object can be used to set loads, evalFunctions as well as perform
        solutions and sensitivities related to transient problems

        Parameters
        ----------
        name : str
            Name to assign problem.
        tInit : float
            Starting time for transient time integration
        tFinal : float
            Ending time for transient time integration
        numSteps : int
            Number of time steps for transient time integration
        options : dict
            Problem-specific options to pass to TransientProblem instance (case-insensitive).
            Defaults to None.

        Returns
        -------
        problem : tacs.problems.TransientProblem
            TransientProblem object used for modeling and solving transient cases.
        """
        problem = tacs.problems.transient.TransientProblem(
            name,
            tInit,
            tFinal,
            numSteps,
            self.assembler,
            self.comm,
            self.outputViewer,
            self.meshLoader,
            self.isNonlinear,
            options,
        )
        # Set with original design vars and coordinates, in case they have changed
        problem.setDesignVars(self.x0)
        problem.setNodes(self.Xpts0)
        return problem

    @postinitialize_method
    def createModalProblem(self, name, sigma, numEigs, options=None):
        """
        Create a new ModalProblem for performing modal analysis.
        This problem can be used to identify the natural frequencies and mode
        shapes of the model through eigenvalue analysis.

        Parameters
        ----------
        name : str
            Name to assign problem.
        sigma : float
            Guess for the lowest eigenvalue.
            This corresponds to the lowest expected frequency squared. (rad^2/s^2)
        numEigs : int
            Number of eigenvalues to solve for.
        options : dict
            Problem-specific options to pass to ModalProblem instance (case-insensitive).
            Defaults to None.

        Returns
        -------
        problem : tacs.problems.ModalProblem
            ModalProblem object used for performing modal eigenvalue analysis.
        """
        problem = tacs.problems.modal.ModalProblem(
            name,
            sigma,
            numEigs,
            self.assembler,
            self.comm,
            self.outputViewer,
            self.meshLoader,
            self.isNonlinear,
            options,
        )
        # Set with original design vars and coordinates, in case they have changed
        problem.setDesignVars(self.x0)
        problem.setNodes(self.Xpts0)
        return problem

    @postinitialize_method
    def createBucklingProblem(self, name, sigma, numEigs, options=None):
        """
        Create a new BucklingProblem for performing linearized buckling analysis.
        This problem can be used to identify the buckling load factors and mode
        shapes of the model through eigenvalue analysis.

        Parameters
        ----------
        name : str
            Name to assign problem.
        sigma : float
            Guess for the lowest eigenvalue.
            This corresponds to the lowest expected buckling load factor.
        numEigs : int
            Number of eigenvalues to solve for.
        options : dict
            Problem-specific options to pass to ModalProblem instance (case-insensitive).
            Defaults to None.

        Returns
        -------
        problem : tacs.problems.BucklingProblem
            BucklingProblem object used for performing buckling eigenvalue analysis.
        """
        problem = tacs.problems.buckling.BucklingProblem(
            name,
            sigma,
            numEigs,
            self.assembler,
            self.comm,
            self.outputViewer,
            self.meshLoader,
            self.isNonlinear,
            options,
        )
        # Set with original design vars and coordinates, in case they have changed
        problem.setDesignVars(self.x0)
        problem.setNodes(self.Xpts0)
        return problem

    @postinitialize_method
    def createTACSProbsFromBDF(self):
        """
        Automatically define tacs problem classes with loads using information contained in BDF file.
        This function assumes all loads are specified in the BDF and allows users to
        skip setting loads in Python.

        Returns
        -------
        structProblems : dict[int, tacs.problems.TACSProblem]
            Dictionary containing a predefined TACSProblem for every loadcase found in the BDF.
            The dictionary keys are the loadcase IDs from the BDF.

        Notes
        -----
        Currently only supports LOAD, FORCE, MOMENT, GRAV, RFORCE, PLOAD2, PLOAD4, TLOAD1, TLOAD2, and DLOAD cards.
        Currently only supports staticProblem (SOL 101), transientProblem (SOL 109), and modalProblems (SOL 103)
        """
        # Make sure cross-referencing is turned on in pynastran
        if self.bdfInfo.is_xrefed is False:
            self.bdfInfo.cross_reference()
            self.bdfInfo.is_xrefed = True

        structProblems = {}

        # If subcases have been added in Nastran, then subCase 0 should not be run
        if len(self.bdfInfo.subcases) > 1:
            skipCaseZero = True
        else:
            skipCaseZero = False

        # Loop through every load set and create a corresponding structural problem
        for subCase in self.bdfInfo.subcases.values():
            if skipCaseZero and subCase.id == 0:
                continue

            if "SUBTITLE" in subCase:
                name = subCase["SUBTITLE"][0]
            else:
                name = "load_set_%.3d" % (subCase.id)

            if self.bdfInfo.sol == 103:
                methodID = subCase["METHOD"][0]
                methodInfo = self.bdfInfo.methods[methodID]
                if methodInfo.v1 is not None:
                    sigma = (2 * np.pi * methodInfo.v1) ** 2
                elif methodInfo.v2 is not None:
                    sigma = (2 * np.pi * methodInfo.v2) ** 2
                else:
                    sigma = 1.0
                if methodInfo.nd is not None:
                    nEigs = methodInfo.nd
                else:
                    nEigs = 20
                problem = self.createModalProblem(name, sigma, nEigs)

            elif self.bdfInfo.sol == 109:
                # Get time step info
                if "TSTEP" in subCase:
                    tStepID = subCase["TSTEP"][0]
                    tStep = self.bdfInfo.tsteps[tStepID]
                    nSteps = tStep.N[0]
                    dt = tStep.DT[0]
                # If no time step info was included, we'll skip this case
                else:
                    self._TACSWarning(
                        f"No TSTEP entry found in control deck for subcase number {subCase.id}, "
                        "skipping case."
                    )
                    continue
                problem = self.createTransientProblem(
                    name, tInit=0.0, tFinal=dt * nSteps, numSteps=nSteps
                )

                # Find dynamic load specified for this subcase
                if "DLOAD" in subCase:
                    dloadsID = subCase["DLOAD"][0]
                    dloadSet, dloadScale = self.bdfInfo.get_reduced_dloads(dloadsID)
                    for dloadInfo, dscale in zip(dloadSet, dloadScale):
                        timeSteps = problem.getTimeSteps()
                        if dloadInfo.type in ["TLOAD1", "TLOAD2"]:
                            if dloadInfo.type == "TLOAD1":
                                loadScales = dloadInfo.get_load_at_time(
                                    timeSteps, dscale
                                )
                            elif dloadInfo.type == "TLOAD2":
                                loadScales = dloadInfo.get_load_at_time(
                                    timeSteps, dscale
                                )
                            if dloadInfo.Type != "LOAD":
                                self._TACSWarning(
                                    "Only 'LOAD' types are supported for "
                                    f"'{dloadInfo.type}' card, but '{dloadInfo.type}' {dloadInfo.sid}, "
                                    f"was specified as {dloadInfo.Type} type"
                                )
                            loadsID = dloadInfo.excite_id
                        else:
                            self._TACSWarning(
                                "Unsupported dload type "
                                f"'{dloadInfo.type}' specified for load set number {dloadInfo.sid},"
                                f" skipping load"
                            )
                            continue
                        # Loop through each time step and add loads to problem
                        for timeIndex, scale in enumerate(loadScales):
                            problem.addLoadFromBDF(timeIndex, loadsID, scale=scale)

            else:
                problem = self.createStaticProblem(name)

                # Find the static load specified for this test case
                if "LOAD" in subCase:
                    # Add loads to problem
                    loadsID = subCase["LOAD"][0]
                    problem.addLoadFromBDF(loadsID)

            # append to list of structural problems
            structProblems[subCase.id] = problem

        return structProblems

    @postinitialize_method
    def writeBDF(self, fileName, problems):
        """
        Write NASTRAN BDF file from problem class.
        Assumes all supplied Problems share the same nodal and design variable values.

        NOTE: Only supports writing loads from StaticProblem types.

        Parameters
        ----------
        fileName: str
            Name of file to write BDF file to.
        problems: tacs.problems.TACSProblem or list[tacs.problems.TACSProblem]
            List of pytacs Problem classes to write BDF file from.
        """
        # Make sure problems is in a list
        if hasattr(problems, "__iter__") == False:
            problems = [problems]
        elif isinstance(problems, dict):
            problems = list(problems.values())
        else:
            problems = list(problems)

        # Check that each problem was created by this pyTACS instance
        for problem in problems:
            if problem.assembler != self.assembler:
                raise self._TACSError(
                    f"This problem instance ({problem.name}) is not associated with this instance of pyTACS."
                )

        # Make sure design variables are up-to-date
        dv_bvec = self.createDesignVec(asBVec=True)
        dv_bvec.getArray()[:] = problems[0].getDesignVars()
        # Transfer all non-local dvs
        dv_bvec.beginDistributeValues()
        dv_bvec.endDistributeValues()

        # Get local node info for each processor
        multNodes = self.getLocalMultiplierNodeIDs()
        globalToLocalNodeIDDict = self.meshLoader.getGlobalToLocalNodeIDDict()
        Xpts_bvec = np.real(problems[0].getNodes())

        # Gather local info to root processor
        allMultNodes = self.comm.gather(multNodes, root=0)
        allGlobalToLocalNodeIDDict = self.comm.gather(globalToLocalNodeIDDict, root=0)
        allXpts = self.comm.gather(Xpts_bvec, root=0)

        # Assemble new BDF file for mesh on root
        if self.comm.rank == 0:
            newBDFInfo = pn.bdf.BDF(debug=False)

            # Write out updated node locations
            nastranNodeIDs = list(self.bdfInfo.node_ids)
            # Loop through each proc and pull out new node locations
            for proc_i in range(self.comm.size):
                xyz = allXpts[proc_i].reshape(-1, 3)
                for tacsGNodeID in allGlobalToLocalNodeIDDict[proc_i]:
                    # Get local node ID
                    tacsLNodeID = allGlobalToLocalNodeIDDict[proc_i][tacsGNodeID]
                    # Get Global nastran ID
                    nastranGNodeID = nastranNodeIDs[tacsGNodeID]
                    # Add node to bdf file (if its not a multiplier node)
                    if tacsLNodeID not in allMultNodes[proc_i]:
                        newBDFInfo.add_grid(nastranGNodeID, xyz[tacsLNodeID])

            # Copy over boundary conditions
            # Set all con IDs to one
            newBDFInfo.spcs[1] = []
            for spcID in self.bdfInfo.spcs:
                for spcCard in self.bdfInfo.spcs[spcID]:
                    newCard = copy.deepcopy(spcCard)
                    newCard.conid = 1
                    newBDFInfo.spcs[1].append(newCard)

            # Write updated properties and elements
            transObjs = {}
            matObjs = []
            conObjs = []
            for compID, propID in enumerate(self.bdfInfo.properties):
                # Get TACS element object
                elemObj = self.meshLoader.getElementObject(compID, 0)
                # get dv nums for element
                dvNums = elemObj.getDesignVarNums(0)
                # Update design variable values
                dvVals = dv_bvec.getValues(dvNums)
                elemObj.setDesignVars(0, dvVals)
                # Get TACS constitutive object for element (if applicable)
                conObj = elemObj.getConstitutive()
                if conObj is not None:
                    # Set the property ID number for the class to be used in the Nastran card
                    conObj.setNastranID(propID)
                    conObjs.append(conObj)
                    # Get TACS material properties object for constitutive (if applicable)
                    matObj = conObj.getMaterialProperties()
                    # May be a single object...
                    if isinstance(matObj, tacs.constitutive.MaterialProperties):
                        if matObj not in matObjs:
                            matObjs.append(matObj)
                    # or a list (plys for composite classes)
                    elif isinstance(matObj, list):
                        for mat_i in matObj:
                            if mat_i not in matObjs:
                                matObjs.append(mat_i)
                # Get TACS transform object for element (if applicable)
                transObj = elemObj.getTransform()
                if transObj is not None:
                    transObjs[compID] = transObj

            # Write material cards from TACS MaterialProperties class
            for i, matObj in enumerate(matObjs):
                matID = i + 1
                matObj.setNastranID(matID)
                newBDFInfo.materials[matID] = matObj.generateBDFCard()

            # Write property/element cards from TACSConstitutive/TACSElement classes
            curCoordID = 1
            for compID, conObj in enumerate(conObjs):
                propID = conObj.getNastranID()
                propCard = conObj.generateBDFCard()
                if propCard is not None:
                    # Copy property comment (may include component name info)
                    # Make sure to remove comment `$` from string
                    propCard.comment = self.bdfInfo.properties[propID].comment[1:]
                    # Add property card to BDF
                    newBDFInfo.properties[propID] = propCard
                elemIDs = self.meshLoader.getGlobalElementIDsForComps(
                    [compID], nastranOrdering=True
                )
                # Convert any transform objects to nastran COORD2R cards, if necessary
                transObj = transObjs.get(compID, None)
                if isinstance(
                    transObj, tacs.elements.ShellRefAxisTransform
                ) or isinstance(transObj, tacs.elements.SpringRefFrameTransform):
                    coordID = curCoordID
                    origin = np.zeros(3)
                    if isinstance(transObj, tacs.elements.SpringRefFrameTransform):
                        vec1, vec2 = transObj.getRefAxes()
                    else:
                        vec1 = transObj.getRefAxis()
                        vec2 = np.random.random(3)
                    # Add COORD2R card to BDF
                    pn.cards.coordinate_systems.define_coord_e123(
                        newBDFInfo,
                        "CORD2R",
                        coordID,
                        origin,
                        xaxis=np.real(vec1),
                        xzplane=np.real(vec2),
                        add=True,
                    )
                    curCoordID += 1
                # We just need the ref vector for these types
                elif isinstance(
                    transObj, tacs.elements.BeamRefAxisTransform
                ) or isinstance(transObj, tacs.elements.SpringRefAxisTransform):
                    vec = transObj.getRefAxis()
                    vec = np.real(vec)
                # Otherwise, there's no transform associated with this element, use default
                else:
                    coordID = None
                # Copy and update element cards
                for elemID in elemIDs:
                    # Create copy of card
                    newCard = copy.deepcopy(self.bdfInfo.elements[elemID])
                    # Copy element comment (may include component name info)
                    # Make sure to remove comment `$` from string
                    newCard.comment = self.bdfInfo.elements[elemID].comment[1:]
                    # Update element coordinate frame info, if necessary
                    if "CQUAD" in newCard.type or "CTRI" in newCard.type:
                        newCard.theta_mcid = coordID
                    elif "CBAR" in newCard.type:
                        newCard.x = vec
                        newCard.g0 = None
                    elif "CBEAM" in newCard.type:
                        newCard.x = vec
                        newCard.g0 = None
                        if propCard.type != "PBEAM":
                            # TACS wrote out a PBAR card that we must convert
                            newPropCard = (
                                pn.cards.properties.beam.PBEAM_init_from_empty()
                            )
                            newPropCard.A[0] = propCard.Area()
                            newPropCard.i1[0] = propCard.I11()
                            newPropCard.i2[0] = propCard.I22()
                            newPropCard.i12[0] = propCard.I12()
                            if hasattr(propCard, "J"):
                                newPropCard.j[0] = propCard.J()
                            else:
                                newPropCard.j[0] = propCard.j
                            newPropCard.comment = propCard.comment
                            propCard = newPropCard
                    elif "CROD" in newCard.type and propCard.type != "PROD":
                        # TACS wrote out a PBAR card that we must convert
                        if hasattr(propCard, "J"):
                            J = propCard.J()
                        else:
                            J = propCard.j
                        newPropCard = pn.cards.properties.rods.PROD(
                            propCard.pid, propCard.mid, propCard.Area(), J
                        )
                        newBDFInfo.properties[propID] = newPropCard
                        newPropCard.comment = propCard.comment
                        propCard = newPropCard
                    elif newCard.type == "CBUSH":
                        if isinstance(transObj, tacs.elements.SpringRefAxisTransform):
                            newCard.x = vec
                            newCard.g0 = None
                        else:
                            newCard.cid = coordID
                    # Add element card to bdf
                    newBDFInfo.elements[elemID] = newCard

            # Copy over masses elements
            for massCard in self.bdfInfo.masses.values():
                elemID = massCard.eid
                # We'll have to create a new CONM2 card in case the point mass is associated with tacs dvs
                if massCard.type == "CONM2":
                    nodeID = massCard.nid
                    elemObj = self.meshLoader.getElementObjectForElemID(
                        elemID, nastranOrdering=True
                    )
                    conObj = elemObj.getConstitutive()
                    M = conObj.evalMassMatrix()
                    mass = np.real(M[0])
                    I11 = np.real(M[15])
                    I22 = np.real(M[18])
                    I33 = np.real(M[20])
                    # Nastran uses negative convention for POI's
                    I12 = -np.real(M[16])
                    I13 = -np.real(M[17])
                    I23 = -np.real(M[19])
                    newBDFInfo.add_conm2(
                        elemID, nodeID, mass, I=[I11, I12, I22, I13, I23, I33]
                    )
                # CONM1's can't be updated by TACS, so we can just copy the original value
                else:
                    newBDFInfo.masses[elemID] = copy.deepcopy(massCard)
                # Copy over comments
                newBDFInfo.masses[elemID].comment = massCard.comment

            # Copy over rigid elements
            newBDFInfo.rigid_elements.update(self.bdfInfo.rigid_elements)

            # Add case control deck for loads
            caseConLines = [
                "TITLE = TACS Analysis Set",
                "ECHO = NONE",
                "DISPLACEMENT(PLOT) = ALL",
                "SPCFORCE(PLOT) = ALL",
                "OLOAD(PLOT) = ALL",
                "FORCE(PLOT,CORNER) = ALL",
                "STRESS(PLOT,CORNER) = ALL",
                "SPC = 1",
            ]
            newBDFInfo.case_control_deck = pn.case_control_deck.CaseControlDeck(
                caseConLines
            )
            # Set solution type to static (101)
            newBDFInfo.sol = 101

        else:
            newBDFInfo = None

        # All procs should wait for root
        self.comm.barrier()

        # Append forces from problem classes
        for i, problem in enumerate(problems):
            if isinstance(problem, tacs.problems.StaticProblem):
                loadCase = i + 1
                problem.writeLoadToBDF(newBDFInfo, loadCase)

        # Write out BDF file
        if self.comm.rank == 0:
            newBDFInfo.write_bdf(
                fileName, size=16, is_double=True, write_header=False, enddata=True
            )

        # All procs should wait for root
        self.comm.barrier()

    @postinitialize_method
    def createAdjacencyConstraint(self, name, options=None):
        """
        Create a new AdjacencyConstraint for calculating
        design variable differences across adjacent components.
        This constraint can be used to ensure that the design variables
        do not change too abruptly across components.
        The formulation is a linear constraint that takes the following form:

        c = dv_i - dv_j

        Where dv_i and dv_j are two design variables in adjacent components.

        Parameters
        ----------
        name : str
            Name to assign constraint.
        options : dict
            Class-specific options to pass to AdjacencyConstraint instance (case-insensitive).
            Defaults to None.

        Returns
        -------
        constraint : tacs.constraints.AdjacencyConstraint
            AdjacencyConstraint object used for calculating constraints.
        """
        constr = tacs.constraints.AdjacencyConstraint(
            name,
            self.assembler,
            self.comm,
            self.outputViewer,
            self.meshLoader,
            options,
        )
        # Set with original design vars and coordinates, in case they have changed
        constr.setDesignVars(self.x0)
        constr.setNodes(self.Xpts0)
        return constr

    @postinitialize_method
    def createDVConstraint(self, name, options=None):
        """
        Create a new DVConstraint for calculating linear constraints based
        on design variables within the same component.

        The constraints are of the form:

            c = a_0 * dv_0 + a_1 * dv_1 + ... + a_n * dv_n

        Where which design variables to include (dv_0, dv_1, etc.)
        and the corresponding weights (a_0, a_1, etc.) are defined by the user.

        Parameters
        ----------
        name : str
            Name to assign constraint.
        options : dict
            Class-specific options to pass to DVConstraint instance (case-insensitive).
            Defaults to None.

        Returns
        -------
        constraint : tacs.constraints.DVConstraint
            DVConstraint object used for calculating constraints.
        """
        constr = tacs.constraints.DVConstraint(
            name,
            self.assembler,
            self.comm,
            self.outputViewer,
            self.meshLoader,
            options,
        )
        # Set with original design vars and coordinates, in case they have changed
        constr.setDesignVars(self.x0)
        constr.setNodes(self.Xpts0)
        return constr

    @postinitialize_method
    def createPanelLengthConstraint(self, name, options=None):
        """Create a new PanelLengthConstraint for enforcing that the panel
        length DV values passed to components match the actual panel lengths.

        Parameters
        ----------
        name : str
            Name to assign constraint.
        options : dict
            Class-specific options to pass to DVConstraint instance (case-insensitive).

        Returns
        ----------
        constraint : tacs.constraints.PanelLengthConstraint
            PanelLengthConstraint object used for calculating constraints.
        """
        constr = tacs.constraints.PanelLengthConstraint(
            name,
            self.assembler,
            self.comm,
            self.outputViewer,
            self.meshLoader,
            options,
        )
        # Set with original design vars and coordinates, in case they have changed
        constr.setDesignVars(self.x0)
        constr.setNodes(self.Xpts0)
        return constr

    @postinitialize_method
    def createVolumeConstraint(self, name, options=None):
        """
        Create a new VolumeConstraint for constraining the size of a closed volume.
        Only shell and solid elements are supported for this constraint.
        For shell elements, the enclosed volume MUST be manifold and water-tight (no missing/internal faces).
        The formulation is a nonlinear constraint based on the nodal coordinates.

        A common example of this is ensuring enough volume in the wingbox for fuel:

            vol_wing >= vol_fuel

        Parameters
        ----------
        name : str
            Name to assign constraint.
        options : dict
            Class-specific options to pass to VolumeConstraint instance (case-insensitive).
            Defaults to None.

        Returns
        -------
        constraint : tacs.constraints.VolumeConstraint
            VolumeConstraint object used for calculating constraints.
        """
        constr = tacs.constraints.VolumeConstraint(
            name,
            self.assembler,
            self.comm,
            self.outputViewer,
            self.meshLoader,
            options,
        )
        # Set with original design vars and coordinates, in case they have changed
        constr.setDesignVars(self.x0)
        constr.setNodes(self.Xpts0)
        return constr

    def getNumComponents(self):
        """
        Return number of components (property) groups found in bdf.

        Returns
        -------
        nComp : int
            Number of components in model
        """
        return self.nComp

    def _createOutputGroups(self):
        """Automatically determine how to split out the output file
        for easier viewing"""

        self.fam = []
        for i in range(self.nComp):
            aux = self.compDescripts[i].split(self.getOption("familySeparator"))
            self.fam.append(aux[0])

        # Uniqify them and sort
        self.fam = sorted(np.unique(self.fam))

        self.compFam = np.zeros(self.nComp, dtype="intc")
        for i in range(self.nComp):
            aux = self.compDescripts[i].split(self.getOption("familySeparator"))
            self.compFam[i] = self.fam.index(aux[0])

    def _createOutputViewer(self):
        """
        Internal method to create the appropriate output viewer
        (TACSToFH5 object) for TACS.
        """

        # Depending on the user-supplied options generate the
        # write_flag
        write_flag = 0
        if self.getOption("writeConnectivity"):
            write_flag |= tacs.TACS.OUTPUT_CONNECTIVITY
        if self.getOption("writeNodes"):
            write_flag |= tacs.TACS.OUTPUT_NODES
        if self.getOption("writeDisplacements"):
            write_flag |= tacs.TACS.OUTPUT_DISPLACEMENTS
        if self.getOption("writeStrains"):
            write_flag |= tacs.TACS.OUTPUT_STRAINS
        if self.getOption("writeStresses"):
            write_flag |= tacs.TACS.OUTPUT_STRESSES
        if self.getOption("writeExtras"):
            write_flag |= tacs.TACS.OUTPUT_EXTRAS
        if self.getOption("writeLoads"):
            write_flag |= tacs.TACS.OUTPUT_LOADS
        if self.getOption("writeCoordinateFrame"):
            write_flag |= tacs.TACS.OUTPUT_COORDINATES

        # Create actual viewer
        if self.getOption("outputElement") is not None:
            elementType = self.getOption("outputElement")
        else:
            # Set the output type based on the first element in the model
            elem = self.meshLoader.getElementObjectForElemID(0, nastranOrdering=False)
            elementType = elem.getElementType()

        self.outputViewer = tacs.TACS.ToFH5(self.assembler, elementType, write_flag)

        # Set the names of each of the output families
        for i in range(len(self.fam)):
            self.outputViewer.setComponentName(i, self.fam[i])

    def _getCompIDs(self, op, *inList):
        """Internal method to return the component IDs mathing
        information in inList"""

        # First recursively flatten the inList in case it was nested:
        inList = self._flatten(inList)

        # Neste list container for compIDs
        compIDs = []

        # Look at each item in list (which is a list because of the *)
        for item in inList:
            compIDs.append([])
            if isinstance(item, int):
                # Integers are easy, just check if in bounds and add:
                if item >= 0 and item < self.nComp:
                    compIDs[-1].append(item)
                else:
                    self._TACSWarning(
                        f"Trying to add component ID of {item}, which\
                    is out of the range 0 <= compID < {self.nComp}"
                    )

            elif isinstance(item, str):
                # This is a little inefficient here; loop over
                # self.compDescripts and see if 'item' (a string) in
                # part of the description. if so add.
                item = item.upper()
                for i in range(self.nComp):
                    if item in self.compDescripts[i].upper():
                        compIDs[-1].append(i)
            else:
                self._TACSWarning(
                    f"Unidentifiable information given for 'include'\
                or 'exclude'. Valid data are integers 0 <= i < {self.nComp}, or \
                strings."
                )

        if op == "and":
            # First convert each entry to a set:
            for i in range(len(compIDs)):
                compIDs[i] = set(compIDs[i])

            # We want to go through and take only the intersection of
            # each of the sets we have found:
            tmp = copy.deepcopy(compIDs[0])

            for i in range(1, len(compIDs)):
                tmp = tmp.intersection(compIDs[i])
            compIDs = tmp

        # Finally, convert to a list
        compIDs = self._flatten(list(compIDs))

        return compIDs

    def _createElements(self, elemCallBack):
        """
        Create all the constitutive objects by calling the
        userSupplied or default callback function

        Parameters
        ----------
        elemCallBack : callable
            Element callback function provided by user or pyTACS
            to set up TACS element objects.
        """

        for i in range(self.nComp):
            # Get a list of compDescripts to help the user
            compDescript = self.compDescripts[i]
            numElements = len(self.elemDescripts[i])
            # TACS component ID
            compID = i
            # Nastran property ID
            propID = list(self.bdfInfo.property_ids)[i]

            # Call the user function
            result = elemCallBack(
                self.dvNum,
                compID,
                compDescript,
                self.elemDescripts[i],
                self.globalDVs,
                propID=propID,
            )

            # For maximum flexibility, multiple pieces of information
            # can be returned. At a minimum, the element objects
            # must be returned!

            # Note: If two objects are returned, the
            # first one is used as the element object list and the
            # second one is treated as a scale list for the added dvs.

            # Check that result is an element object instance or .
            numFoundElements = 0
            scaleList = None

            if isinstance(result, tuple):
                elemObjects = result[0]
                if hasattr(result[1], "__iter__"):
                    # Iterable item, the scale list:
                    scaleList = result[1]
                elif isinstance(result[1], numbers.Number):
                    scaleList = [result[1]]
                else:
                    print(result[1])
                    # Don't know what it is:
                    self._TACSWarning(
                        "Could not identify objects returned \
                    from elemCallBack. Valid return objects are: \
                    A list of TACS element objects (required, first), \
                    an iterable object \
                    (eg, list or array) containing the scaling parameters \
                    for the added design variables (optional, second). The \
                    string representation of the offending object is: \
                    '%s'"
                        % repr(result[1])
                    )

            else:
                elemObjects = result

            if isinstance(elemObjects, tacs.TACS.Element):
                # There was only one element, recast it as a list and continue
                elemObjects = [elemObjects]
                numFoundElements += 1
            elif isinstance(elemObjects, list):
                # Multiple elements were returned, count how many
                for object in elemObjects:
                    if isinstance(object, tacs.TACS.Element):
                        numFoundElements += 1
                    else:
                        self._TACSError(
                            f"Object of type {type(object)} returned in elemCallBack function "
                            f"is not a valid TACS element object. The \
                               string representation of the offending object is: \
                               '{repr(object)}'"
                        )

            if numFoundElements != numElements:
                raise self._TACSError(
                    f"Unexpected number of element objects \
                    returned from user-supplied elemCallBack function. \
                    {numElements} element types ({repr(self.elemDescripts[i])}) are contained in Component {i}, \
                    but {numFoundElements} element objects were returned by elemCallback."
                )

            # Now determine the number of design variables. This is
            # NOT as simple as just getting the number of design
            # variables; Not all variables added in the conObject are
            # 'new' variables, some of the variable number may have
            # been already used.
            newVars = []
            for elemObject in elemObjects:
                dvs = elemObject.getDesignVarNums(0)

                if len(dvs) > 0:
                    # We will also check if the user screwed up. That is
                    # make sure that for added variables, the are
                    # continuous starting at self.dvNum
                    for var in dvs:
                        if var >= self.dvNum:
                            newVars.append(var)

            # Remove repeated dv nums from list
            newVars = np.unique(newVars)
            newVars.sort()

            if len(newVars) > 0:
                # Now the length of newVars must the same as
                # newVars[-1]-newVars[0] + 1
                if not len(newVars) == newVars[-1] - newVars[0] + 1:
                    raise self._TACSError(
                        "Inconsistent design variables detected. "
                        "The added design variables are not continuous."
                        f" The added design variables are {repr(newVars)}."
                    )

            # Finally, increment the dv counter
            self.dvNum += len(newVars)

            if len(newVars) > 0:
                if scaleList is None:
                    self.scaleList.extend(np.ones(len(newVars)))
                else:
                    # Make sure that the scaleList is the correct length.
                    if len(scaleList) != len(newVars):
                        self._TACSWarning(
                            f"An incorrect number of scale variables \
                        were returned. There were {len(newVars)} variables added, but only \
                        {len(scaleList)} scale variables returned. The scale for these \
                        variables will be set to 1.0. The scale variables are {repr(scaleList)}."
                        )
                        self.scaleList.extend(np.ones(len(newVars)))
                    else:
                        self.scaleList.extend(scaleList)

            # Loop through every element type in this component,
            # there may be multiple (e.g CQUAD4 + CTRIA3)
            for j, elemObject in enumerate(elemObjects):
                # Set component-specific family id
                elemObject.setComponentNum(self.compFam[i])
                # Set each of the elements for this component
                self.meshLoader.setElementObject(i, j, elemObject)
                # set varsPerNode
                elemVarsPerNode = elemObject.getVarsPerNode()
                if self.varsPerNode is None:
                    self.varsPerNode = elemVarsPerNode
                elif self.varsPerNode != elemVarsPerNode:
                    raise self._TACSError(
                        "Model references elements with differing numbers of variables per node "
                        f"({self.varsPerNode} and {elemVarsPerNode}). "
                        "All elements must use same number of variables to be compatible."
                    )

        # If varsPerNode still hasn't been set (because there were no elements added in the callback)
        # Default to 6
        if self.varsPerNode is None:
            self.varsPerNode = 6
