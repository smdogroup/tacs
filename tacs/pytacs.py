#!/usr/bin/python
"""
pytacs - The Python wrapper for the TACS solver

This python interface is designed to provide a easier interface to the
c-layer of TACS. It combines all the functionality of the old pyTACS
and pyTACS_Mesh. User-supplied hooks allow for nearly complete
customization of any or all parts of the problem setup. There are two
main parts of this module: The first deals with setting up the TACS
problem including reading the mesh, setting design variables,
functions, constraints etc (Functionality in the former
pyTACS_Mesh). The second part deals with solution of the structural
analysis and gradient computations.

Copyright (c) 2013 by Dr. G.K.W. Kenway
All rights reserved. Not to be used for commercial purposes.

Developers:
-----------
- Dr. G.K.W. Kenway (GKK)

History
-------
    v. 1.0  - pyTACS initial implementation
"""
# =============================================================================
# Imports
# =============================================================================
from __future__ import print_function
import copy
import os
import numbers
import numpy
import time

import numpy as np
from mpi4py import MPI
import warnings
import tacs.TACS, tacs.constitutive, tacs.elements, tacs.functions, tacs.problems.static
from .utilities import BaseUI
from tacs.pymeshloader import pyMeshLoader

DEG2RAD = np.pi / 180.0

warnings.simplefilter('default')
try:
    from collections import OrderedDict
except ImportError:
    try:
        from ordereddict import OrderedDict
    except ImportError:
        print("Could not find any OrderedDict class. "
              "For python 2.6 and earlier, use:"
              "\n pip install ordereddict")

class pyTACS(BaseUI):

    def __init__(self, fileName, comm=None, dvNum=0,
                 scaleList=None, **kwargs):
        """
        The class for working with a TACS structure

        Parameters
        ----------
        fileName : str
            The filename of the BDF file to load.

        comm : MPI Intracomm
            The comm object on which to create the pyTACS object.

        dvNum : int
            An user supplied offset to the design variable
            numbering. This is typically used with tacs+tripan when
            geometric variables have already been added and assigned
            global tacs numberings.

        scaleList: list
            when dvNum is non zero, the scaleList must be same size
            as the number of design variables already added. i.e.
            len(scaleList) = dvNum
        """

        self.objectName = 'pyTACS'

        startTime = time.time()

        # Default Option List
        defOpts = {
            # Meshloader options
            'printDebug': [bool, False],

            # Output Options
            'outputElement': [int, None],
            'writeBDF': [bool, False],
            'writeConnectivity': [bool, True],
            'writeNodes': [bool, True],
            'writeDisplacements': [bool, True],
            'writeStrains': [bool, True],
            'writeStresses': [bool, True],
            'writeExtras': [bool, True],
            'writeCoordinateFrame': [bool, False],
            'familySeparator': [str, '/'],
            'printTiming': [bool, False],
            'printIterations': [bool, True],

        }

        # Data type (real or complex)
        self.dtype = tacs.TACS.dtype

        # Set the communicator and rank -- defaults to MPI_COMM_WORLD
        if comm is None:
            comm = MPI.COMM_WORLD
        self.comm = comm
        self.rank = comm.rank

        # Process the default options which are added to self.options
        # under the 'defaults' key. Make sure the key are lower case
        self.options = {}
        def_keys = defOpts.keys()
        self.options['defaults'] = {}
        for key in def_keys:
            self.options['defaults'][key.lower()] = defOpts[key]
            self.options[key.lower()] = defOpts[key]

        # Process the user-supplied options
        koptions = kwargs.pop('options', {})
        kopt_keys = koptions.keys()
        for key in kopt_keys:
            self.setOption(key, koptions[key])

        importTime = time.time()

        # Create and load mesh loader object.
        debugFlag = self.getOption('printDebug')
        self.meshLoader = pyMeshLoader(self.comm, self.dtype, debugFlag)
        self.meshLoader.scanBdfFile(fileName)
        self.bdfName = fileName
        # Save pynastran bdf object
        self.bdfInfo = self.meshLoader.getBDFInfo()

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
        self.compIDBounds = {}
        self.addedCompIDs = set()

        # List of initial coordinates
        self.Xpts0 = None
        # List of initial designvars
        self.x0 = None

        # Variables per node for model
        self.varsPerNode = None

        # TACS assembler object
        self.assembler = None

        initFinishTime = time.time()
        if self.getOption('printTiming'):
            self.pp('+--------------------------------------------------+')
            self.pp('|')
            self.pp('| TACS Init Times:')
            self.pp('|')
            self.pp('| %-30s: %10.3f sec' % ('TACS Module Time', importTime - startTime))
            self.pp('| %-30s: %10.3f sec' % ('TACS Meshload Time', meshLoadTime - importTime))
            self.pp('| %-30s: %10.3f sec' % ('TACS DV Processing Time', DVPreprocTime - meshLoadTime))
            self.pp('| %-30s: %10.3f sec' % ('TACS Finalize Initialization Time', initFinishTime - DVPreprocTime))
            self.pp('|')
            self.pp('| %-30s: %10.3f sec' % ('TACS Total Initialization Time', initFinishTime - startTime))
            self.pp('+--------------------------------------------------+')

    def addGlobalDV(self, descript, value,
                    lower=None, upper=None, scale=1.0):
        """
        This function allows adding design variables that are not
        cleanly associated with a particular constiutive object. One
        example is the pitch of the stiffeners for blade stiffened
        panels; It often is the same for many different constitutive
        objects. By calling this function, the internal dvNum counter
        is incremented and the user doesn\'t have to worry about
        it.

        Parameters
        ----------
        descript : str
            A user supplied string that can be used to retrieve the
            variable number and value elemCallBackFunction.

        value : float
            Initial value for variable.
        lower : float
            Lower bound. May be None for unbounded
        upper : float
            Upper bound. May be None for unbounded
        scale : float
            Scale factor for variable

        Returns
        -------
        None, but the information is provided to the user in the
        elemCallBack function
        """
        self.globalDVs[descript] = {'num': self.dvNum,
                                     'value': value,
                                     'lowerBound': lower,
                                     'upperBound': upper}
        self.dvNum += 1
        self.scaleList.append(scale)

    def selectCompIDs(self, include=None, exclude=None,
                      includeBounds=None, nGroup=1, includeOp='or',
                      excludeOp='or', projectVector=None, **kwargs):
        """
        This is the most important function of the entire setup
        process. The basic idea is as follow: We have a list of nComp
        which are the component descriptions. What we need is a way of
        generating subgroups of these for the purposes of adding
        design variables, constitutive objects, KS domains and mass
        domains. All of these operations boil down to selecting a
        subset of the compIDs.

        This function attemps to support as many ways as possible to
        select parts of the structure. Easy and efficient selection of
        parts is critical to the end user.

        Methods of selction:

        1. include, integer, string, list of integers and/or strings: The
           simpliest and most direct way of selecting a component. The
           user supplies the index of the componentID, a name or partial
           name, or a list of a combination of both.

           For exammple::

            # Select the 11th component
            selectCompIDs(include=10)

            # Select the first and fifth component
            selectCompIDs(include=[0, 4])

            # Select any component containing 'rib.00'
            selectCompIDs(include='rib.00')

            # Select any components containg 'rib.00' and 'rib.10'
            selectCompIDs(include=['rib.00', 'rib.10'])

            # Select any componet containing 'rib.00', the 11th
            # component and any component containing 'spar'
            # (This is probably not advisable!)
            selectCompIDs(include=['rib.00', 10, 'spar'])

        2. Exclude, operates similarally to 'include'. The behaviour
           of exclude is identical to include above, except that
           component ID's that are found using 'exclude' are
           'subtracted' from those found using include. A special
           case is treated if 'include' is NOT given: if only an
           exclude list is given, this implies the selection of all
           compID's EXCEPT the those in exclude.

           For example::

               # This will return will [0, 1, 2, 3, 5, ..., nComp-1]
               selectCompIDs(exclude = 4)

               # This will return [0, 1, 4, 5, ..., nComp-1]
               selectCompIDs(exclude = [2, 3]) will return

               # This will return components that have 'ribs' in the
               # componet ID, but not those that have 'le_ribs' in the
               # componet id.
               selectCompIDs(include='ribs', exclude='le_ribs')

        3. includeBounds, list of componets defining a region inside
           of which 'include' components will be selected. This
           functionality uses a geometric approach to select the compIDs.
           All components within the project 2D convex hull are included.
           Therefore it is essential to split up concave include regions
           into smaller convex regions. Use multiple calls to selectCompIDs to
           accumulate multiple regions.

           For example::

               # This will select upper skin components between the
               # leading and trailing edge spars and between ribs 1 and 4.
               selectCompIDs(include='U_SKIN', includeBound=
                   ['LE_SPAR', 'TE_SPAR', 'RIB.01', 'RIB.04'])

        4. nGroup: The number of groups to divde the found componets
           into. Generally this will be 1. However, in certain cases, it
           is convient to create multiple groups in one pass.

           For example::

             # This will 'evenly' create 10 groups on all components
             # containing LE_SPAR. Note that once the componets are
             # selected, they are sorted **alphetically** and assigned
             # sequentially.
             selectCompIDs(include='LE_SAPR', nGroup=10)

           nGroup can also be negative. If it is negative, then a single
           design variable group is added to each of the found
           components.

           For example::

             # will select all components and assign a design variable
             # group to each one.
             selectCompIDs(nGroup=-1)

        includeOp, str: 'and' or 'or'. Selects the logical operation
        used for item in 'include' option. For example:

        selectCompIDs(include=['LE_SPAR', 'TE_SPAR'],
        includeOpt='or') will select the LE_SPAR and TE_SPAR
        components (default behaviour).

        selectCompIDs(include=['RIB', 'SEG.01'], includeOpt='and')
        will select any componet with 'RIB' in the description AND
        'SEG.01' in the description.
        """

        # Defaults
        includeIDs = numpy.arange(self.nComp)
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
                self.TACSWarning('nGroup=%d is larger than the number of\
                selected components=%d. nGroup will be clipped to %d' %
                            (nGroup, len(compIDs), nGroup))
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
                return [alist[i * length // wanted_parts:
                              (i + 1) * length // wanted_parts]
                        for i in range(wanted_parts)]

            ind = split_list(ind, nGroup)

            # Finally assemble the nested list of componet IDs
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
        Return pynastran bdf object. This object can be used interactively
        to parse information (nodes, elements, loads etc) included in the bdf file.

        Returns
        -------
        bdfInfo : BDF object
            pyNastran bdf object.
        """
        return self.bdfInfo

    def getCompNames(self, compIDs):
        """
        Return a list of component descriptions for the given component
        IDs. compIDs should come from a call to selectCompIDs

        Parameters
        ----------
        compIDs : list
            List of integers of the compIDs numbers

        Returns
        -------
        compDescript : list
            List of strings of the names of the corresponding compIDs
        """
        compIDs = self._flatten(compIDs)
        compDescripts = []
        for i in range(len(compIDs)):
            compDescripts.append(self.compDescripts[compIDs[i]])

        return compDescripts

    def initialize(self, elemCallBack=None):
        """
        This is the 'last' function to be called during the setup. The
        user should have already added all the design variables,
        domains ect. before this function is call. This function
        finializes the problem initialization and cannot be changed at
        later time. If a elemCallBack function is not provided by the user,
        we will use pyNastran to generate one automatically from element
        properties provided in the BDF file.

        Parameters
        ----------
        elemCallBack : python function handle

           The calling sequence for elemCallBack **must** be as
           follows::

             def elemCallBack(dvNum, compID, compDescript, elemDescripts,
                             globalDVs, **kwargs):

           The dvNum is the current counter which must be used by the
           user when creating constitutive object with design
           variables.

           compID is the ID number used by tacs to reference this property group.
           Use kwargs['propID'] to get the corresponding Nastran property ID that
           is read in from the BDF.

           compDescript is the component descriptions read in from the BDF file

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

        self.assembler = self.meshLoader.createTACSAssembler(self.varsPerNode)

        self._createOutputViewer()

        # Initial design variable vec if necessary
        self.dvs0 = self.assembler.createDesignVec()
        self.assembler.getDesignVars(self.dvs0)

        # Store original node locations read in from bdf file
        self.Xpts0 = self.assembler.createNodeVec()
        self.assembler.getNodes(self.Xpts0)

        # Store initial design variable values
        self.x0 = self.assembler.createDesignVec()
        self.assembler.getDesignVars(self.x0)

    def _elemCallBackFromBDF(self):
        """
        Automatically setup elemCallBack using information contained in BDF file.
        This function assumes all material properties are specified in the BDF.
        """

        # Check if any properties are in the BDF
        if self.bdfInfo.missing_properties:
            raise self.TACSError("BDF file '%s' has missing properties cards. "
                        "Set 'debugPrint' option to True for more information."
                        "User must define own elemCallBack function." % (self.bdfName))

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
                elemDict[propertyID]['elements'] = []
                elemDict[propertyID]['dvs'] = {}
            elemDict[propertyID]['elements'].append(element)

        # Create a dictionary to sort all design variables
        for dv in self.bdfInfo.dvprels:
            propertyID = self.bdfInfo.dvprels[dv].pid
            dvName = self.bdfInfo.dvprels[dv].pname_fid
            self.dvNum = max(self.dvNum, self.bdfInfo.dvprels[dv].dvids[0])
            elemDict[propertyID]['dvs'][dvName] = self.bdfInfo.dvprels[dv]
        # Create option for user to specify scale values in BDF
        self.scaleList = [1.0] * self.dvNum

        # Callback function to return appropriate tacs MaterialProperties object
        # For a pynastran mat card
        def matCallBack(matInfo):
            # First we define the material property object
            if matInfo.type == 'MAT1':
                mat = tacs.constitutive.MaterialProperties(rho=matInfo.rho, E=matInfo.e,
                                                           nu=matInfo.nu, ys=matInfo.St,
                                                           alpha=matInfo.a)
            elif matInfo.type == 'MAT8':
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
                # TODO: add alpha
                mat = tacs.constitutive.MaterialProperties(rho=rho, E1=E1, E2=E2, nu12=nu12, G12=G12, G13=G13, G23=G23,
                                                           Xt=Xt, Xc=Xc, Yt=Yt, Yc=Yc, S12=S12)
            else:
                raise self.TACSError("Unsupported material type '%s' for material number %d. " % (matInfo.type, matInfo.mid))

            return mat

        def elemCallBack(dvNum, compID, compDescript, elemDescripts, globalDVs, **kwargs):
            # Initialize scale list for design variables we will add
            scaleList = []

            # Get the Nastran property ID
            propertyID = kwargs['propID']
            propInfo = self.bdfInfo.properties[propertyID]
            elemInfo = elemDict[propertyID]['elements'][0]

            # First we define the material object
            mat = None
            # This property only references one material
            if hasattr(propInfo, 'mid_ref'):
                matInfo = propInfo.mid_ref
                mat = matCallBack(matInfo)
            # This property references multiple materials (maybe a laminate)
            elif hasattr(propInfo, 'mids_ref'):
                mat = []
                for matInfo in propInfo.mids_ref:
                    mat.append(matCallBack(matInfo))

            # Next we define the constitutive object
            if propInfo.type == 'PSHELL':  # Nastran isotropic shell
                kcorr = propInfo.tst

                if 'T' in elemDict[propertyID]['dvs']:
                    thickness = elemDict[propertyID]['dvs']['T'].dvids_ref[0].xinit
                    tNum = elemDict[propertyID]['dvs']['T'].dvids[0] - 1
                    minThickness = elemDict[propertyID]['dvs']['T'].dvids_ref[0].xlb
                    maxThickness = elemDict[propertyID]['dvs']['T'].dvids_ref[0].xub
                    name = elemDict[propertyID]['dvs']['T'].dvids_ref[0].label
                    self.scaleList[tNum - 1] = elemDict[propertyID]['dvs']['T'].coeffs[0]
                else:
                    thickness = propInfo.t
                    tNum = -1
                    minThickness = 0.0
                    maxThickness = 1e20

                con = tacs.constitutive.IsoShellConstitutive(mat, t=thickness,
                                                             tlb=minThickness, tub=maxThickness, tNum=tNum)

            elif propInfo.type == 'PCOMP':  # Nastran composite shell
                numPlies = propInfo.nplies
                plyThicknesses = []
                plyAngles = []
                plyMats = []

                # if the laminate is symmetric, mirror the ply indices
                if propInfo.lam == 'SYM':
                    plyIndices = list(range(numPlies / 2))
                    plyIndices.extend(plyIndices[::-1])
                else:
                    plyIndices = range(numPlies)

                # Loop through plies and setup each entry in layup
                for ply_i in plyIndices:
                    plyThicknesses.append(propInfo.thicknesses[ply_i])
                    plyMat = tacs.constitutive.OrthotropicPly(plyThicknesses[ply_i], mat[ply_i])
                    plyMats.append(plyMat)
                    plyAngles.append(propInfo.thetas[ply_i] * DEG2RAD)

                # Convert thickness/angles to appropriate numpy array
                plyThicknesses = np.array(plyThicknesses, dtype=self.dtype)
                plyAngles = np.array(plyAngles, dtype=self.dtype)

                if propInfo.lam is None or propInfo.lam in ['SYM', 'MEM']:
                    # Discrete laminate class (not for optimization)
                    con = tacs.constitutive.CompositeShellConstitutive(plyMats, plyThicknesses, plyAngles)
                    # Need to add functionality to consider only membrane in TACS for type = MEM

                else:
                    raise self.TACSError("Unrecognized LAM type '%s' for PCOMP number %d." % (propInfo.lam, propertyID))

            elif propInfo.type == 'PSOLID':  # Nastran solid property
                if 'T' in elemDict[propertyID]['dvs']:
                    thickness = elemDict[propertyID]['dvs']['T'].dvids_ref[0].xinit
                    tNum = elemDict[propertyID]['dvs']['T'].dvids[0] - 1
                    minThickness = elemDict[propertyID]['dvs']['T'].dvids_ref[0].xlb
                    maxThickness = elemDict[propertyID]['dvs']['T'].dvids_ref[0].xub
                    name = elemDict[propertyID]['dvs']['T'].dvids_ref[0].label
                    self.scaleList[tNum - 1] = elemDict[propertyID]['dvs']['T'].coeffs[0]
                else:
                    thickness = 1.0
                    tNum = -1
                    minThickness = 0.0
                    maxThickness = 10.0

                con = tacs.constitutive.SolidConstitutive(mat, t=thickness,
                                                          tlb=minThickness, tub=maxThickness, tNum=tNum)

            elif propInfo.type == 'PBUSH':  # Nastran spring
                k = numpy.zeros(6)
                for j in range(len(k)):
                    if (propInfo.Ki[j]):
                        k[j] = propInfo.Ki[j]
                con = tacs.constitutive.DOFSpringConstitutive(k=k)

            else:
                raise self.TACSError("Unsupported property type '%s' for property number %d. " % (propInfo.type, propertyID))

            # Set up transform object which may be required for certain elements
            transform = None
            if propInfo.type in ['PSHELL', 'PCOMP']:
                mcid = elemDict[propertyID]['elements'][0].theta_mcid_ref
                if mcid:
                    if mcid.type == 'CORD2R':
                        refAxis = mcid.i
                        transform = tacs.elements.ShellRefAxisTransform(refAxis)
                    else:  # Don't support spherical/cylindrical yet
                        raise self.TACSError("Unsupported material coordinate system type "
                                    "'%s' for property number %d." % (mcid.type, propertyID))
            elif propInfo.type == 'PBUSH':
                if elemDict[propertyID]['elements'][0].cid_ref:
                    refAxis_i = elemDict[propertyID]['elements'][0].cid_ref.i
                    refAxis_j = elemDict[propertyID]['elements'][0].cid_ref.j
                    transform = tacs.elements.SpringRefFrameTransform(refAxis_i, refAxis_j)
                elif elemDict[propertyID]['elements'][0].x[0]:
                    refAxis = numpy.array(elemDict[propertyID]['elements'][0].x) \
                              - elemDict[propertyID]['elements'][0].nodes_ref[0].xyz
                    transform = tacs.elements.SpringRefAxisTransform(refAxis)
                elif elemDict[propertyID]['elements'][0].g0_ref:
                    refAxis = elemDict[propertyID]['elements'][0].g0_ref.xyz \
                              - elemDict[propertyID]['elements'][0].nodes_ref[0].xyz
                    transform = tacs.elements.SpringRefAxisTransform(refAxis)

            # Finally set up the element objects belonging to this component
            elemList = []
            for descript in elemDescripts:
                if descript in ['CQUAD4', 'CQUADR']:
                    elem = tacs.elements.Quad4Shell(transform, con)
                elif descript in ['CQUAD9', 'CQUAD']:
                    elem = tacs.elements.Quad9Shell(transform, con)
                elif descript in ['CTRIA3', 'CTRIAR']:
                    elem = tacs.elements.Tri3Shell(transform, con)
                elif 'CTETRA' in descript:
                    # May have variable number of nodes in card
                    nnodes = len(elemInfo.nodes)
                    if nnodes == 4:
                        basis = tacs.elements.LinearTetrahedralBasis()
                    elif nnodes == 10:
                        basis = tacs.elements.QuadraticTetrahedralBasis()
                    else:
                        raise self.TACSError("TACS does not currently support CTETRA elements with %d nodes." % nnodes)
                    model = tacs.elements.LinearElasticity3D(con)
                    elem = tacs.elements.Element3D(model, basis)
                elif descript in ['CHEXA8', 'CHEXA']:
                    basis = tacs.elements.LinearHexaBasis()
                    model = tacs.elements.LinearElasticity3D(con)
                    elem = tacs.elements.Element3D(model, basis)
                elif descript == 'CBUSH':
                    elem = tacs.elements.SpringElement(transform, con)
                else:
                    raise self.TACSError("Unsupported element type "
                                "'%s' specified for property number %d." % (descript, propertyID))
                elemList.append(elem)

            return elemList, scaleList

        return elemCallBack

    def getOrigDesignVars(self):
        """
        get the original design variables that were specified with
        during assembler creation.

        Returns
        ----------
        x : array
            The current design variable vector set in tacs.

        """
        if self.assembler is None:
            raise self.TACSError("TACS assembler has not been created. "
                        "Assembler must created first by running 'initalize' method.")

        return self.x0.getArray().copy()

    def createDesignVec(self, asBVec=False):
        """
        Create a new tacs distributed design vector.
        Values are initialized to zero.

        Parameters
        ----------
        asBVec : bool
            Flag that determines whether to return
            design vector as tacs BVec (True) or numpy array (False).
            Defaults to False.

        Returns
        ----------
        x : ndarray or BVec
            Distributed design variable vector
        """
        if self.assembler is None:
            raise self.TACSError("TACS assembler has not been created. "
                        "Assembler must created first by running 'initalize' method.")
        xVec = self.assembler.createDesignVec()
        if asBVec:
            return xVec
        else:
            return xVec.getArray()

    def getNumDesignVars(self):
        """
        Return the number of design variables on this processor.
        """
        if self.assembler is None:
            raise self.TACSError("TACS assembler has not been created. "
                        "Assembler must created first by running 'initalize' method.")

        return self.x0.getSize()

    def getTotalNumDesignVars(self):
        """
        Return the number of design variables across all processors.
        """
        if self.assembler is None:
            raise self.TACSError("TACS assembler has not been created. "
                        "Assembler must created first by running 'initalize' method.")

        return self.dvNum

    def getOrigNodes(self):
        """
        Return the original mesh coordiantes read in from the meshLoader.

        Returns
        -------
        coords : array
            Structural coordinate in array of size (N, 3) where N is
            the number of structural nodes on this processor.
        """
        if self.assembler is None:
            raise self.TACSError("TACS assembler has not been created. "
                        "Assembler must created first by running 'initalize' method.")

        return self.Xpts0.getArray().copy()

    def createNodeVec(self, asBVec=False):
        """
        Create a new tacs distributed node vector.
        Values are initialized to zero.

        Parameters
        ----------
        asBVec : bool
            Flag that determines whether to return
            node vector as tacs BVec (True) or numpy array (False).
            Defaults to False.

        Returns
        ----------
        xpts : ndarray or BVec
            Distributed node coordinate vector
        """
        if self.assembler is None:
            raise self.TACSError("TACS assembler has not been created. "
                                 "Assembler must created first by running 'initalize' method.")
        xptVec = self.assembler.createNodeVec()
        if asBVec:
            return xptVec
        else:
            return xptVec.getArray()

    def getNumOwnedNodes(self):
        """
        Get the number of nodes owned by this processor.
        """
        if self.assembler is None:
            raise self.TACSError("TACS assembler has not been created. "
                        "Assembler must created first by running 'initalize' method.")

        return self.assembler.getNumOwnedNodes()

    def createVec(self, asBVec=False):
        """
        Create a new tacs distributed state variable vector.
        Values are initialized to zero.

        Parameters
        ----------
        asBVec : bool
            Flag that determines whether to return
            state vector as tacs BVec (True) or numpy array (False).
            Defaults to False.

        Returns
        ----------
        vars : ndarray or BVec
            Distributed state variable vector
        """
        if self.assembler is None:
            raise self.TACSError("TACS assembler has not been created. "
                                 "Assembler must created first by running 'initalize' method.")
        vars = self.assembler.createVec()
        if asBVec:
            return vars
        else:
            return vars.getArray()

    def getVarsPerNode(self):
        """
        Get the number of variables per node for the model.
        """
        if self.assembler is None:
            raise self.TACSError("TACS assembler has not been created. "
                        "Assembler must created first by running 'initalize' method.")

        return self.assembler.getVarsPerNode()

    def applyBCsToVec(self, vec):
        """
        Applies zeros to boundary condition dofs in input vector.
        """
        if self.assembler is None:
            raise self.TACSError("TACS assembler has not been created. "
                        "Assembler must created first by running 'initalize' method.")

        varVec = self.assembler.createVec()
        varArray = varVec.getArray()

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

    def createStaticProblem(self, name, options={}):
        """
        Create a new staticProblem for modeling a static load cases.
        This object can be used to set loads, evalFunctions as well as perform
        solutions and sensitivities related to static problems

        Parameters
        ----------
        name : str
            Name to assign problem.
        options : dict
            Problem-specific options to pass to StaticProblem instance.

        Returns
        ----------
        problem : StaticProblem
            StaticProblem object used for modeling and solving static cases.
        """
        if self.assembler is None:
            raise self.TACSError("TACS assembler has not been created. "
                        "Assembler must created first by running 'initalize' method.")

        problem = tacs.problems.static.StaticProblem(name, self.assembler, self.comm,
                                                     self.outputViewer, self.meshLoader, options)
        # Set with original design vars and coordinates, in case they have changed
        problem.setDesignVars(self.x0)
        problem.setNodes(self.Xpts0)
        return problem

    def createTransientProblem(self, name, tInit, tFinal, numSteps, options={}):
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
            Problem-specific options to pass to TransientProblem instance.

        Returns
        ----------
        problem : TransientProblem
            TransientProblem object used for modeling and solving transient cases.
        """
        if self.assembler is None:
            raise self.TACSError("TACS assembler has not been created. "
                                 "Assembler must created first by running 'initalize' method.")

        problem = tacs.problems.transient.TransientProblem(name, tInit, tFinal, numSteps,
                                                           self.assembler, self.comm, self.outputViewer,
                                                           self.meshLoader,
                                                           options)
        # Set with original design vars and coordinates, in case they have changed
        problem.setDesignVars(self.x0)
        problem.setNodes(self.Xpts0)
        return problem

    def createModalProblem(self, name, sigma, numEigs, options={}):
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
            Problem-specific options to pass to ModalProblem instance.

        Returns
        ----------
        problem : ModalProblem
            ModalProblem object used for performing modal eigenvalue analysis.
        """
        if self.assembler is None:
            raise self.TACSError("TACS assembler has not been created. "
                        "Assembler must created first by running 'initalize' method.")

        problem = tacs.problems.modal.ModalProblem(name, sigma, numEigs,
                                                   self.assembler, self.comm, self.outputViewer, self.meshLoader,
                                                   options)
        # Set with original design vars and coordinates, in case they have changed
        problem.setDesignVars(self.x0)
        problem.setNodes(self.Xpts0)
        return problem

    def createTACSProbsFromBDF(self):
        """
        Automatically define tacs problem classes with loads using information contained in BDF file.
        This function assumes all loads are specified in the BDF and allows users to
        skip setting loads in Python.

        Returns
        ----------
        structProblems : dict[TACSProblems]
            Dictionary containing a predfined TACSProblem for every loadcase found int the BDF.
            The dictionary keys are the loadcase IDs from the BDF.

        Notes
        -----
        Currently only supports LOAD, FORCE, MOMENT, PLOAD2, and PLOAD4 cards.
        Currently only supports staticProblem (SOL 101) and modalProblems (SOL 103)
        """

        if self.assembler is None:
            raise self.TACSError("TACS assembler has not been created. "
                        "Assembler must created first by running 'initalize' method.")

        # Make sure cross-referencing is turned on in pynastran
        if self.bdfInfo.is_xrefed is False:
            self.bdfInfo.cross_reference()
            self.bdfInfo.is_xrefed = True

        vpn = self.varsPerNode
        loads = self.bdfInfo.loads
        nloads = len(loads)

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

            if 'SUBTITLE' in subCase.params:
                name = subCase.params['SUBTITLE'][0]
            else:
                name = 'load_set_%.3d' % (subCase.id)

            if self.bdfInfo.sol == 103:
                methodID = subCase.params['METHOD'][0]
                methodInfo = self.bdfInfo.methods[methodID]
                if methodInfo.v1 is not None:
                    sigma = (2 * np.pi * methodInfo.v1) ** 2
                elif methodInfo.v2 is not None:
                    sigma = (2 * np.pi * methodInfo.v2) ** 2
                else:
                    sigma = 1.0
                nEigs = methodInfo.maxset
                problem = self.createModalProblem(name, sigma, nEigs)

            else:
                problem = self.createStaticProblem(name)

                if 'LOAD' in subCase.params:
                    loadsID = subCase.params['LOAD'][0]
                    # Get loads and scalers for this load case ID
                    loadSet, loadScale, _ = self.bdfInfo.get_reduced_loads(loadsID)
                    # Loop through every load in set and add it to problem
                    for loadInfo, scale in zip(loadSet, loadScale):
                        # Add any point force or moment cards
                        if loadInfo.type == 'FORCE' or loadInfo.type == 'MOMENT':
                            nodeID = loadInfo.node_ref.nid

                            loadArray = numpy.zeros(vpn)
                            if loadInfo.type == 'FORCE' and vpn >= 3:
                                loadArray[:3] += scale * loadInfo.scaled_vector
                            elif loadInfo.type == 'MOMENT' and vpn >= 6:
                                loadArray[3:6] += scale * loadInfo.scaled_vector
                            problem.addLoadToNodes(nodeID, loadArray, nastranOrdering=True)

                        # Add any pressure loads
                        # Pressure load card specific to shell elements
                        elif loadInfo.type == 'PLOAD2':
                            elemIDs = loadInfo.eids
                            pressure = scale * loadInfo.pressure
                            problem.addPressureToElements(elemIDs, pressure, nastranOrdering=True)

                        # Alternate more general pressure load type
                        elif loadInfo.type == 'PLOAD4':
                            self._addPressureFromPLOAD4(problem, loadInfo, scale)

                        else:
                            self.TACSWarning("Unsupported load type "
                                        " '%s' specified for load set number %d, skipping load" %(loadInfo.type, loadInfo.sid))

            # append to list of structural problems
            structProblems[subCase.id] = problem

        return structProblems

    def _addPressureFromPLOAD4(self, staticProb, loadInfo, scale=1.0):
        """
        Add pressure to tacs static problem from pynastran PLOAD4 card.
        Should only be called by createTACSProbsFromBDF and not directly by user.
        """
        # Dictionary mapping nastran element face indices to TACS equivilent numbering
        nastranToTACSFaceIDDict = {'CTETRA4': {1: 1, 2: 3, 3: 2, 4: 0},
                                   'CTETRA': {2: 1, 4: 3, 3: 2, 1: 0},
                                   'CHEXA': {1: 4, 2: 2, 3: 0, 4: 3, 5: 0, 6: 5}}

        # We don't support pressure variation across elements, for now just average it
        pressure = scale * np.mean(loadInfo.pressures)
        for elemInfo in loadInfo.eids_ref:
            elemID = elemInfo.eid

            # Get the correct face index number based on element type
            if 'CTETRA' in elemInfo.type:
                for faceIndex in elemInfo.faces:
                    if loadInfo.g1 in elemInfo.faces[faceIndex] and \
                            loadInfo.g34 not in elemInfo.faces[faceIndex]:
                        # For some reason CTETRA4 is the only element that doesn't
                        # use ANSYS face numbering convention by default
                        if len(elemInfo.nodes) == 4:
                            faceIndex = nastranToTACSFaceIDDict['CTETRA4'][faceIndex]
                        else:
                            faceIndex = nastranToTACSFaceIDDict['CTETRA'][faceIndex]
                        # Positive pressure is inward for solid elements, flip pressure if necessary
                        # We don't flip it for face 0, because the normal for that face points inward by convention
                        # while the rest point outward
                        if faceIndex != 0:
                            pressure *= -1.0
                        break

            elif 'CHEXA' in elemInfo.type:
                for faceIndex in elemInfo.faces:
                    if loadInfo.g1 in elemInfo.faces[faceIndex] and \
                            loadInfo.g34 in elemInfo.faces[faceIndex]:
                        faceIndex = nastranToTACSFaceIDDict['CHEXA'][faceIndex]
                        # Pressure orientation is flipped for solid elements per Nastran convention
                        pressure *= -1.0
                        break

            elif 'CQUAD' in elemInfo.type or 'CTRIA' in elemInfo.type:
                # Face index doesn't matter for shells, just use 0
                faceIndex = 0

            else:
                raise self.TACSError("Unsupported element type "
                            "'%s' specified for PLOAD4 load set number %d." % (elemInfo.type, loadInfo.sid))

            # Figure out whether or not this is a traction based on if a vector is defined
            if np.linalg.norm(loadInfo.nvector) == 0.0:
                staticProb.addPressureToElements(elemID, pressure, faceIndex,
                                           nastranOrdering=True)
            else:
                trac = pressure * loadInfo.nvector
                staticProb.addTractionToElements(elemID, trac, faceIndex,
                                           nastranOrdering=True)

    def getNumComponents(self):
        """
        Return number of components (property) groups found in bdf.
        """
        return self.nComp

    def _createOutputGroups(self):
        """Automatically determine how to split out the output file
        for easier viewing"""

        self.fam = []
        for i in range(self.nComp):
            aux = self.compDescripts[i].split(self.getOption('familySeparator'))
            self.fam.append(aux[0])

        # Uniqify them and sort
        self.fam = sorted(numpy.unique(self.fam))

        self.compFam = numpy.zeros(self.nComp, dtype='intc')
        for i in range(self.nComp):
            aux = self.compDescripts[i].split(self.getOption('familySeparator'))
            self.compFam[i] = self.fam.index(aux[0])

    def _createOutputViewer(self):
        """
        Internal function to create the appropriate output viewer
        (TACSToFH5 object) for TACS.
        """

        # Depending on the user supplied options generate the
        # write_flag
        write_flag = 0
        if self.getOption('writeConnectivity'):
            write_flag |= tacs.TACS.OUTPUT_CONNECTIVITY
        if self.getOption('writeNodes'):
            write_flag |= tacs.TACS.OUTPUT_NODES
        if self.getOption('writeDisplacements'):
            write_flag |= tacs.TACS.OUTPUT_DISPLACEMENTS
        if self.getOption('writeStrains'):
            write_flag |= tacs.TACS.OUTPUT_STRAINS
        if self.getOption('writeStresses'):
            write_flag |= tacs.TACS.OUTPUT_STRESSES
        if self.getOption('writeExtras'):
            write_flag |= tacs.TACS.OUTPUT_EXTRAS
        if self.getOption('writeCoordinateFrame'):
            write_flag |= tacs.TACS.OUTPUT_COORDINATES

        # Create actual viewer
        if self.getOption('outputElement') is not None:
            elementType = self.getOption('outputElement')
        elif self.varsPerNode == 6 or self.varsPerNode == 7:
            elementType = tacs.TACS.BEAM_OR_SHELL_ELEMENT
        elif self.varsPerNode == 3:
            elementType = tacs.TACS.SOLID_ELEMENT
        else:
            self.TACSWarning("'outputElement' not specified in options. "
                             "No elements will be written out in f5 file.")
            elementType = tacs.TACS.ELEMENT_NONE

        self.outputViewer = tacs.TACS.ToFH5(
            self.assembler, elementType, write_flag)

        # Set the names of each of the output families
        for i in range(len(self.fam)):
            self.outputViewer.setComponentName(i, self.fam[i])

    def _getCompIDs(self, op, *inList):
        """ Internal function to return the component IDs mathing
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
                    self.TACSWarning('Trying to add component ID of %d, which\
                    is out of the range 0 <= compID < %d' % (item, self.nComp))

            elif isinstance(item, str):
                # This is a little inefficinet here; loop over
                # self.compDescripts and see if 'item' (a string) in
                # part of the description. if so add.
                item = item.upper()
                for i in range(self.nComp):
                    if item in self.compDescripts[i].upper():
                        compIDs[-1].append(i)
            else:
                self.TACSWarning('Unidentifiable information given for \'include\'\
                or \'exclude\'. Valid data are integers 0 <= i < %d, or \
                strings.' % self.nComp)

        if op == 'and':
            # First convert each entry to a set:
            for i in range(len(compIDs)):
                compIDs[i] = set(compIDs[i])

            # We want to go though and take only the intersection of
            # each of the sets we have found:
            tmp = copy.deepcopy(compIDs[0])

            for i in range(1, len(compIDs)):
                tmp = tmp.intersection(compIDs[i])
            compIDs = tmp

        # Finally convert to a list
        compIDs = self._flatten(list(compIDs))

        return compIDs

    def _createElements(self, elemCallBack):
        """
        Create all the constitutive objects by calling the
        userSupplied or default callback function
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
            result = elemCallBack(self.dvNum, compID, compDescript, self.elemDescripts[i], self.globalDVs,
                                  propID=propID)

            # For maximum flexibiliy, multiple pieces of information
            # can be returned. At a minimum, the element objects
            # must be returned!

            # Note: If two objects are returned, the
            # first one is used as the element object list and the
            # second one is treated as a scale list for the added dvs.

            # Check that result is an element object instance or .
            foundElem = False
            numFoundElements = 0
            scaleList = None

            if isinstance(result, tuple):
                elemObjects = result[0]
                if hasattr(result[1], '__iter__'):
                    # Iterable item, the scale list:
                    scaleList = result[1]
                elif isinstance(result[1], numbers.Number):
                    scaleList = [result[1]]
                else:
                    print(result[1])
                    # Don't know what it is:
                    self.TACSWarning("Could not identify objects returned \
                    from elemCallBack. Valid return objects are: \
                    A list of TACS element objects (required, first), \
                    an iterable object \
                    (eg, list or array) containing the scaling parameters \
                    for the added design variables (optional, second). The \
                    string representation of the offending object is: \
                    '%s'" % repr(result[1]))

            else:
                elemObjects = result

            if isinstance(elemObjects, tacs.TACS.Element) and numElements == 1:
                # There was only one element, recast it as a list and continue
                elemObjects = [elemObjects]
                numFoundElements += 1
            elif isinstance(elemObjects, list):
                # Multiple elements were returned, count how many
                for object in elemObjects:
                    if isinstance(object, tacs.TACS.Element):
                        numFoundElements += 1
                    else:
                        self.TACSError("Object of type %s returned in elemCallBack function "
                              "is not a valid TACS element object. The \
                               string representation of the offending object is: \
                               '%s'"%(type(object), repr(object)))

            if numFoundElements != numElements:
                raise self.TACSError("Could not find all required element objects in the "
                            "return arguments from user-supplied "
                            "elemCallBack function. {} element types ({}) are contained in Component {}, "
                            "but only {} were returned by elemCallback.".format(numElements, repr(self.elemDescripts[i]),
                                                                                i, numFoundElements))

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
            newVars = numpy.unique(newVars)
            newVars.sort()

            if len(newVars) > 0:
                # Now the length of newVars must the same as
                # newVars[-1]-newVars[0] + 1
                if not len(newVars) == newVars[-1] - newVars[0] + 1:
                    raise self.TACSError("Inconsistent design variables detected. "
                                "The added design variables are not continuous."
                                " The added design varibales are %s." %
                                repr(newVars))

            # Finally increment the dvcounter
            self.dvNum += len(newVars)

            if len(newVars) > 0:
                if scaleList is None:
                    self.scaleList.extend(numpy.ones(len(newVars)))
                else:
                    # Make sure that the scaleList is the correct length.
                    if len(scaleList) != len(newVars):
                        self.TACSWarning('An incorrect number of scale variables \
                        were returned. There were %d variables added, but only \
                        %d scale variables returned. The scale for these \
                        variables will be set to 1.0. The scale variables are %s.' % (
                            len(newVars), len(scaleList), repr(scaleList)))
                        self.scaleList.extend(numpy.ones(len(newVars)))
                    else:
                        self.scaleList.extend(scaleList)

            # Loop through every element type in this component,
            # there may be multiple (e.g CQUAD4 + CTRIA3)
            for j, elemObject in enumerate(elemObjects):
                # Set each of the elements for this component
                self.meshLoader.setElementObject(i, j, elemObject)
                # set varsPerNode
                elemVarsPerNode = elemObject.getVarsPerNode()
                if self.varsPerNode is None:
                    self.varsPerNode = elemVarsPerNode
                elif self.varsPerNode != elemVarsPerNode:
                    raise self.TACSError("Model references elements with differing numbers of variables per node (%d and %d). "
                                "All elements must use same number of variables to be compatible."%(self.varsPerNode,
                                                                                                    elemVarsPerNode))

        # If varsPerNode still hasn't been set (because there were no elements added in the callback)
        # Default to 6
        if self.varsPerNode is None:
            self.varsPerNode = 6
