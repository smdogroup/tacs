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
import tacs.TACS, tacs.constitutive, tacs.elements, tacs.functions, tacs.problems
from pyNastran.bdf.bdf import read_bdf

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

DEG2RAD = numpy.pi / 180.0


class pyTACS(object):

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

        startTime = time.time()

        # Default Option List
        defOpts = {

            'probname': [str, 'defaultName'],
            'outputdir': [str, './'],

            # Solution Options
            'solutionType': [str, 'linear'],
            'KSMSolver': [str, 'GMRES'],
            'orderingType': [str, 'ND'],
            'varsPerNode': [int, 6],
            'PCFillLevel': [int, 1000],
            'PCFillRatio': [float, 20.0],
            'subSpaceSize': [int, 10],
            'nRestarts': [int, 15],
            'flexible': [int, 1],
            'L2Convergence': [float, 1e-12],
            'L2ConvergenceRel': [float, 1e-12],
            'useMonitor': [bool, False],
            'monitorFrequency': [int, 10],
            'maxLoadCases': [int, 20],
            'MaxItersPCReset': [int, 25],
            'nNLIter': [int, 10],
            'NLSolver': [str, 'FNR'],
            'resNormUB': [float, 1e20],

            # selectCompID Options
            'projectVector': [list, [0.0, 1.0, 0.0]],

            # Output Options
            'outputElement': [int, None],
            'writeBDF': [bool, False],
            'writeSolution': [bool, True],
            'writeConnectivity': [bool, True],
            'writeNodes': [bool, True],
            'writeDisplacements': [bool, True],
            'writeStrains': [bool, True],
            'writeStresses': [bool, True],
            'writeExtras': [bool, True],
            'writeCoordinateFrame': [bool, False],
            'familySeparator': [str, '/'],
            'numberSolutions': [bool, True],
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
        self._scanBdfFile(fileName)
        self.bdfName = fileName

        meshLoadTime = time.time()

        # Retrieve the number of components. This is the maximum
        # number of unique constituive objects possible in this model.
        self.nComp = self.bdfInfo.nproperties

        # Set the starting dvNum and scaleList
        self.dvNum = dvNum
        self.scaleList = scaleList
        if scaleList is None:
            self.scaleList = []

        DVPreprocTime = time.time()

        # List of DV groups
        self.specialDVs = {}
        self.compIDBounds = {}
        self.addedCompIDs = set()
        self.varName = 'struct'
        self.coordName = 'Xpts'

        self._usedLCNums = 0
        self.curSP = None
        self.curLoadCase = None
        self.doDamp = False
        self._factorOnNext = True
        self._PCfactorOnNext = False
        # List of functions
        self.functionList = OrderedDict()
        self.adjointList = OrderedDict()
        self.dIduList = OrderedDict()
        self.dvSensList = OrderedDict()
        self.xptSensList = OrderedDict()

        # Element IDs (tacs ordering) partitioned to this proc
        self.globalToLocalNodeIDDict = None

        # List of initial coordinates
        self.coords0 = None

        # Norms
        self.initNorm = 0.0
        self.startNorm = 0.0
        self.finalNorm = 0.0
        # Flag for mat/vector creation
        self._variablesCreated = False

        # TACS assembler and creator objects
        self.creator = None
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

    ####### pyNastran BDF parsing methods ########

    def _scanBdfFile(self, fileName):
        """
        Scan nastran bdf file using pyNastran's bdf parser.
        We also set up arrays that will be require later to build tacs.
        """

        # Only print debug info on root
        if self.comm.rank == 0:
            debugPrint = True
        else:
            debugPrint = False

        # Read in bdf file as pynsatran object
        # By default we avoid cross-referencing unless we actually need it,
        # since its expensive for large models
        self.bdfInfo = read_bdf(fileName, validate=False, xref=False, debug=debugPrint)
        # Set flag letting us know model is not xrefed yet
        self.bdfInfo.is_xrefed = False
        self.bdfInfo.no_properties = False

        # If no property cards were found in bdf, we have to add dummy cards
        # so pynastran doesn't through errors when cross-referencing
        if self.bdfInfo.nproperties == 0:
            # Flag to warn the user that the property cards are placeholders
            # and should not be read in using pytacs _elemCallBackFromBDF method
            self.bdfInfo.no_properties = True
            # If no material properties were found,
            # add dummy properties and materials
            matID = 1
            E = 70.0
            G = 35.0
            nu = 0.3
            self.bdfInfo.add_mat1(matID, E, G, nu)

            # Loop through all elements and add dummy property, as necessary
            for element_id in self.bdfInfo.elements:
                element = self.bdfInfo.elements[element_id]
                if element.pid not in self.bdfInfo.property_ids:
                    self.bdfInfo.add_pbar(element.pid, matID)

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
        self.compDescript = []

        # Populate list entries with default values
        for pID in self.bdfInfo.property_ids:
            self.elemDescripts.append([])
            self.elemObjectNumByComp.append([])
            # Check if there is a Femap label for this component
            propComment = self.bdfInfo.properties[pID].comment
            if '$ Femap Property' in propComment:
                # Pick off last word from comment, this is the name
                propName = propComment.split()[-1]
                self.compDescript.append(propName)
            #  Default component name
            else:
                self.compDescript.append('untitled')

        # Element connectivity information
        self.elemConnectivity = [None] * self.bdfInfo.nelements
        self.elemConnectivityPointer = [None] * (self.bdfInfo.nelements + 1)
        self.elemConnectivityPointer[0] = 0
        elementTypeCounter = 0
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
                self.elemObjectNumByComp[componentID].append(elementTypeCounter)
                elementTypeCounter += 1

            # Find the index number corresponding to the element object number for this component
            componentTypeIndex = self.elemDescripts[componentID].index(elementType)
            self.elemObjectNumByElem[tacsElementID] = self.elemObjectNumByComp[componentID][componentTypeIndex]

            # We've identified a ICEM property label
            if 'Shell element data for family' in element.comment:
                componentName = element.comment.split()[-1]
                self.compDescript[componentID] = componentName

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
        self.elemObjects = [None] * elementTypeCounter

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

    def _elemCallBackFromBDF(self):
        """
        Automatically setup elemCallBack using information contained in BDF file.
        This function assumes all material properties are specified in the BDF.
        """
        # Make sure cross-referencing is turned on in pynastran
        if self.bdfInfo.is_xrefed is False:
            self.bdfInfo.cross_reference()
            self.bdfInfo.is_xrefed = True

        # Check if any properties are in the BDF
        if self.bdfInfo.no_properties:
            raise Error("BDF file '%s' has no properties included in it. "
                        "User must define own elemCallBack function." % (self.bdfName))

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
                raise Error("Unsupported material type '%s' for material number %d. " % (matInfo.type, matInfo.mid))
            
            return mat

        def elemCallBack(dvNum, compID, compDescript, elemDescripts, specialDVs, **kwargs):
            # Initialize scale list for design variables we will add
            scaleList = []
            
            # Get the Nastran property ID
            propertyID = kwargs['propID']
            propInfo = self.bdfInfo.properties[propertyID]
            elemInfo = elemDict[propertyID]['elements'][0]
            
            # First we define the material object
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
                    plyIndices = range(numPlies / 2)
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
                    raise Error("Unrecognized LAM type '%s' for PCOMP number %d." % (propInfo.lam, propertyID))

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

            else:
                raise Error("Unsupported property type '%s' for property number %d. " % (propInfo.type, propertyID))

            # Set up transform object which may be required for certain elements
            transform = None
            if hasattr(elemInfo, 'theta_mcid_ref'):
                mcid = elemDict[propertyID]['elements'][0].theta_mcid_ref
                if mcid:
                    if mcid.type == 'CORD2R':
                        refAxis = mcid.i
                        transform = tacs.elements.ShellRefAxisTransform(refAxis)
                    else:  # Don't support spherical/cylindrical yet
                        raise Error("Unsupported material coordinate system type "
                                    "'%s' for property number %d." % (mcid.type, propertyID))

            # Finally set up the element objects belonging to this component
            elemList = []
            for descript in elemDescripts:
                if descript in ['CQUAD4', 'CQUADR']:
                    elem = tacs.elements.Quad4Shell(transform, con)
                elif descript in ['CQUAD9', 'CQUAD']:
                    elem = tacs.elements.Quad9Shell(transform, con)
                elif descript in ['CTRIA3', 'CTRIAR']:
                    elem = tacs.elements.Tri3Shell(transform, con)
                elif descript in ['CTETRA10', 'CTETRA']:
                    basis = tacs.elements.QuadraticTetrahedralBasis()
                    model = tacs.elements.LinearElasticity3D(con)
                    elem = tacs.elements.Element3D(model, basis)
                elif descript in ['CHEXA8', 'CHEXA']:
                    basis = tacs.elements.LinearHexaBasis()
                    model = tacs.elements.LinearElasticity3D(con)
                    elem = tacs.elements.Element3D(model, basis)
                else:
                    raise Error("Unsupported element type "
                                "'%s' specified for property number %d." % (descript, propertyID))
                elemList.append(elem)

            return elemList, scaleList

        return elemCallBack

    def getBDFInfo(self):
        """
        Return pynastran bdf object.
        """
        return self.bdfInfo

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
        elementIDs = self.bdfInfo.get_element_ids_dict_with_pids()[propertyID + 1]
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

    def _getGlobalToLocalNodeIDDict(self):
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
        globalToLocalNodeIDDict = { gID : lID for lID, gID in enumerate(ownedElementIDs) }

        return globalToLocalNodeIDDict

    def addSpecialDV(self, descript, value,
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
        self.specialDVs[descript] = {'num': self.dvNum,
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

        if includeBounds is not None:
                includeBoundIDs = \
                    self._getCompIDsInBounds(includeBounds, projectVector)

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
                TACSWarning('nGroup=%d is larger than the number of\
                selected components=%d. nGroup will be clipped to %d' %
                            (nGroup, len(compIDs), nGroup), self.comm)
                nGroup = len(compIDs)

            # Pluck out the component descriptions again and we will
            # sort them
            compDescript = []
            for i in range(len(compIDs)):
                compDescript.append(self.compDescript[compIDs[i]])

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

    # TODO: Most options still not supported, most importantly safetyFactor
    def addFunction(self, funcName, funcHandle, include=None, exclude=None,
                    includeBound=None, compIDs=None, KSWeight=80.0, safetyFactor=1.0,
                    volume=1.0, axisIndex=0, axisOffset=0, funcType=0):
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

        include : varries
            Argument passed to selctCompIDs. See this function for
            more information

        exclude : varries
            Argument passed to selctCompIDs. See this function for
            more information

        includeBounds : varries
            Argument passed to selctCompIDs. See this function for
            more information

        compIDs: list
            List of compIDs to select. Alternative to selectCompIDs
            arguments.

        KSWeight : float
            Weight for KS function. Only used if KS function is being
            created

        safetyFactor : float
            Internal safety factor to use for failure and buckling functions.

        volume: float
            Volume value used in denominator to compute AverageTemperature.
            Use an area for 2D/shell elements.

        axisIndex : int
            Index of axis for MassMoment and AggregateDisplacement function.
            0 for x, 1 for y, and 2 for z.
        axisOffset : flaot
            Offset for axis in MassMoment calc.

        funcType : int
            A constant that determines the type of function. This
            argument is only used for specific functions. E.g., KSFailure
            may use the discrete (0) or continous (1) form.
            See the tacs documentation for the specific options.
        """

        # First we will get the required domain, but only if both
        # include is None and exclude is None. If so, just use the
        # entire domain:

        # Note nGroup is one since we only want exactly one domain
        if compIDs is None:
            compIDs = self.selectCompIDs(include, exclude, includeBound,
                                         nGroup=1)[0]

        # Flatten and get element numbers on each proc corresponding to specified compIDs
        compIDs = self._flatten(compIDs)
        elemIDs = self.getLocalElementIDsForComps(compIDs)

        # We have to call certain functions differently:
        f = tacs.functions
        if funcHandle in [f.KSFailure, f.KSTemperature]:
            self.functionList[funcName] = funcHandle(self.assembler, KSWeight)
        elif funcHandle in [f.AverageTemperature]:
            self.functionList[funcName] = funcHandle(self.assembler, volume)
        else:
            try:
                self.functionList[funcName] = funcHandle(self.assembler)
            except:
                TACSWarning("Function type %s is not currently supported "
                            "in pyTACS" % funcHandle, self.comm)
                return

        # Finally set the domain information
        self.functionList[funcName].setDomain(elemIDs)

        # Create additional tacs BVecs to hold adjoint and sens info
        self.adjointList[funcName] = self.assembler.createVec()
        self.dIduList[funcName] = self.assembler.createVec()
        self.dvSensList[funcName] = self.assembler.createDesignVec()
        self.xptSensList[funcName] = self.assembler.createNodeVec()

        return compIDs

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
            compDescripts.append(self.compDescript[compIDs[i]])

        return compDescripts

    def getElementDescripts(self):
        """
        Get nested list containing all element types owned by each component group
        example: [['CQUAD4', 'CTRIA3], ['CQUAD4'], ['CQUAD4', CQUAD9']]
        """
        return self.elemDescripts

    def getFunctionKeys(self):
        """Return a list of the current function key names"""
        return list(self.functionList.keys())

    def setStructProblem(self, structProblem):
        """Set the structProblem. This function can be called by the
        user but typically will be called automatically by functions
        that accept a structProblem object.

        Parameters
        ----------
        structProblem : instance of pyStruct_problem
            Description of the sturctural problem to solve
            """

        if structProblem is self.curSP:

            return

        if self.comm.rank == 0:
            print('+' + '-' * 70 + '+')
            print('|  Switching to Struct Problem: %-39s|' % structProblem.name)
            print('+' + '-' * 70 + '+')

        try:
            structProblem.tacsData
        except AttributeError:
            structProblem.tacsData = TACSLoadCase()
            structProblem.tacsData.F = self.assembler.createVec()
            structProblem.tacsData.u = self.assembler.createVec()
            structProblem.tacsData.auxElems = tacs.TACS.AuxElements()

            # We will assign this structProblem a 'load case number'
            i = self._usedLCNums
            if i < self.getOption('maxLoadCases'):
                structProblem.tacsData.lcnum = self._usedLCNums
                self._usedLCNums += 1
            else:
                raise Error("The maximum number of internal load cases "
                            "conditions, %d, has beeen exceeded. Use a larger "
                            "value for the 'maxLoadCases' option." %
                            self.getOption('maxLoadCases'))

        # We are now ready to associate self.curSP with the supplied SP
        self.curSP = structProblem
        self.curSP.adjointRHS = None
        # Force and displacement vectors for problem
        self.F = self.curSP.tacsData.F
        self.u = self.curSP.tacsData.u
        # Set auxiliary elements for adding tractions/pressures
        self.auxElems = self.curSP.tacsData.auxElems
        self.assembler.setAuxElements(self.auxElems)

        # Create numpy array representation for easier access to vector values
        vpn = self.getOption('varsPerNode')
        self.F_array = self.F.getArray()
        self.u_array = self.u.getArray()
        self.F_array = self.F_array.reshape(len(self.F_array) // vpn, vpn)
        self.u_array = self.u_array.reshape(len(self.u_array) // vpn, vpn)

        # Put a short-cut to the loadcase num for the current loadcase
        self.curLoadCase = self.curSP.tacsData.lcnum
        self.assembler.setVariables(self.u)

        # Reset the Aitken acceleration for multidisciplinary analyses
        self.doDamp = False

    def _setupTACSCreator(self):
        """
        Setup TACSCreator object responsible for creating TACSAssembler
        """
        varsPerNode = self.getOption('varsPerNode')
        self.creator = tacs.TACS.Creator(self.comm, varsPerNode)

        # Append RBE elements to element list, these are not setup by the user
        for rbe in self.bdfInfo.rigid_elements.values():
            if rbe.type == 'RBE2':
                self._addTACSRBE2(rbe)
            elif rbe.type == 'RBE3':
                self._addTACSRBE3(rbe)
            else:
                raise NotImplementedError("Rigid element of type '{}' is not supported".format(rbe.type))

        if self.comm.rank == 0:
            # Set connectivity for all elements
            ptr = np.array(self.elemConnectivityPointer, dtype=np.intc)
            conn = np.array(self._flatten(self.elemConnectivity), dtype=np.intc)
            objectNums = np.array(self.elemObjectNumByElem, dtype=np.intc)
            self.creator.setGlobalConnectivity(self.bdfInfo.nnodes, ptr, conn, objectNums)

            # Set up the boundary conditions
            bcnodes = []
            bcdofs = []
            bcptr = [0]
            bcvals = []
            for spc_id in self.bdfInfo.spcs:
                for spc in self.bdfInfo.spcs[spc_id]:
                    # Loop through every node specifed in this spc and record bc info
                    for j, nastranNode in enumerate(spc.nodes):
                        tacsNode = self.idMap(nastranNode, self.nastranToTACSNodeIDDict)
                        bcnodes.append(tacsNode)
                        bcptr.append(bcptr[-1])

                        for dof in range(varsPerNode):
                            dofAsString = '%d' % (dof + 1)
                            if spc.type == 'SPC':
                                # each node may have its own dofs uniquely constrained
                                constrained_dofs = spc.components[j]
                                # The boundary condition may be forced to a non-zero value
                                bcvals.append(spc.enforced[j])
                            else:  # SPC1?
                                # All nodes always have the same dofs constrained
                                constrained_dofs = spc.components
                                # This boundary condition is always 0
                                bcvals.append(0.0)
                            # Find if dof number occurs in constraint string
                            loc = constrained_dofs.find(dofAsString)
                            # if dof is constrained, add it to list
                            if loc > -1:
                                bcdofs.append(dof)
                                bcptr[-1] += 1

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

    def _addTACSRBE2(self, rbeInfo):
        """
        Method to automatically set up RBE2 element from bdf file for user.
        User should *NOT* set these up in their elemCallBack function.
        """
        varsPerNode = self.getOption('varsPerNode')
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

    def _addTACSRBE3(self, rbeInfo):
        """
        Method to automatically set up RBE3 element from bdf file for user.
        User should *NOT* set these up in their elemCallBack function.
        """
        varsPerNode = self.getOption('varsPerNode')
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

    def createTACSAssembler(self, elemCallBack=None):
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
                             specialDVs, **kwargs):

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

           specialDVs is a dictionary containing information about any
           special DVs that have been added.

           elemCallBack must return a list containing as many TACS element
           objects as there are element types in elemDescripts (one for each).

        """

        if elemCallBack is None:
            elemCallBack = self._elemCallBackFromBDF()
        self._createOutputGroups()
        self._createElements(elemCallBack)

        self._setupTACSCreator()
        self.assembler = self.creator.createTACS()

        self.globalToLocalNodeIDDict = self._getGlobalToLocalNodeIDDict()

        self._createVariables()
        self._createOutputViewer()

        # Initial set of nodes for geometry manipuation if necessary
        self.coords0 = self.getCoordinates()

    ####### Static load methods ########

    def addLoadToComponents(self, structProblem, compIDs, F, averageLoad=False):
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
            F = numpy.atleast_2d(F)

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
                self.addLoadToComponents(structProblem, compID, F[i], averageLoad=True)

        # Average one force vector over all components
        else:
            F = np.atleast_1d(F)

            self.setStructProblem(structProblem)

            # First determine the actual physical nodal location in the
            # original BDF ordering of the nodes we want to add forces
            # to. Only the root rank need do this:
            uniqueNodes = None
            if self.comm.rank == 0:
                allNodes = []
                compIDs = set(compIDs)
                for cID in compIDs:
                    tmp = self.getConnectivityForComp(cID, nastranOrdering=True)
                    allNodes.extend(self._flatten(tmp))

                # Now just unique all the nodes:
                uniqueNodes = numpy.unique(allNodes)

            uniqueNodes = self.comm.bcast(uniqueNodes, root=0)

            # Now generate the final average force vector
            Favg = F / len(uniqueNodes)

            self.addLoadToNodes(structProblem, uniqueNodes, Favg, nastranOrdering=True)

            # Write out a message of what we did:
            self._info("Added a fixed load of %s to %d components, "
                       "distributed over %d nodes." % (
                           repr(F), len(compIDs), len(uniqueNodes)),
                       maxLen=80, box=True)

    def addLoadToPoints(self, structProblem, points, F):
        """"
        The function is used to add a fixed point load of F to the
        selected physical locations, points. A closest point search is
        used to determine the FE nodes that are the closest to the
        requested nodes. It is most efficient if many point loads are
        necessary that points and F, contain many entries.

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
        """
        try:
            from scipy.spatial import cKDTree
        except:
            raise Error("scipy.spatial "
                        "must be available to use addLoadToPoints")

        points = numpy.atleast_2d(points)
        F = numpy.atleast_2d(F)

        # If the user only specified one force vector,
        # we assume the force should be the same for each node
        if F.shape[0] == 1:
            F = np.repeat(F, [len(points)], axis=0)
        # If the dimensions still don't match, raise an error
        elif F.shape[0] != len(points):
            raise Error("Number of forces must match number of points,"
                        " {} forces were specified for {} points".format(F.shape[0], len(points)))

        vpn = self.getOption('varsPerNode')
        if len(F[0]) != vpn:
            raise Error("Length of force vector must match varsPerNode specified "
                        "for problem, which is {}, "
                        "but length of vector provided was {}".format(vpn, len(F[0])))

        self.setStructProblem(structProblem)

        # Pull out the local nodes on the proc and search "points" in the tree
        self.assembler.getNodes(self.Xpts)
        localNodes = np.real(self.Xpts.getArray())
        nNodes = len(localNodes) // 3
        xNodes = localNodes.reshape((nNodes, 3)).copy()
        tree = cKDTree(xNodes)

        d, index = tree.query(points, k=1)

        # Now figure out which proc has the best distance for this

        for i in range(len(points)):
            proc = self.comm.allreduce((d[i], self.comm.rank), op=MPI.MINLOC)
            print((i, self.comm.rank, proc, d[i], index[i], F[i]))
            if proc[1] == self.comm.rank:
                # Add contribution to global force array
                self.F_array[index[i], :] += F[i]

    def addLoadToNodes(self, structProblem, nodeIDs, F, nastranOrdering=False):
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
        nodeIDs = numpy.atleast_1d(nodeIDs)
        F = numpy.atleast_2d(F)

        numNodes = len(nodeIDs)

        # If the user only specified one force vector,
        # we assume the force should be the same for each node
        if F.shape[0] == 1:
            F = np.repeat(F, [numNodes], axis=0)
        # If the dimensions still don't match, raise an error
        elif F.shape[0] != numNodes:
            raise Error("Number of forces must match number of nodes,"
                        " {} forces were specified for {} node IDs".format(F.shape[0], numNodes))

        vpn = self.getOption('varsPerNode')
        if len(F[0]) != vpn:
            raise Error("Length of force vector must match varsPerNode specified "
                        "for problem, which is {}, "
                        "but length of vector provided was {}".format(vpn, len(F[0])))

        # First find the cooresponding local node ID on each processor
        localNodeIDs = self.getLocalNodeIDsFromGlobal(nodeIDs, nastranOrdering)

        # Set the structural problem
        self.setStructProblem(structProblem)

        # Flag to make sure we find all user-specified nodes
        nodeFound = np.zeros(numNodes, dtype=int)

        # Loop through every node and if it's owned by this processor, add the load
        for i, nodeID in enumerate(localNodeIDs):
            # The node was found on this proc
            if nodeID >= 0:
                # Add contribution to global force array
                self.F_array[nodeID, :] += F[i]
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
                TACSWarning("Can't add load to node ID {} ({} ordering), node not found in model. "
                            "Double check BDF file.".format(nodeIDs[i], orderString), self.comm)

    def addTractionToComponents(self, structProblem, compIDs, tractions,
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
        elemIDs = self.getGlobalElementIDsForComps(compIDs, nastranOrdering=False)
        # Add tractions element by element
        self.addTractionToElements(structProblem, elemIDs, tractions, faceIndex, nastranOrdering=False)

        # Write out a message of what we did:
        self._info("Added a fixed traction of %s to %d components, "
                   "distributed over %d elements." % (
                       repr(tractions), len(compIDs), len(elemIDs)),
                   maxLen=80, box=True)

    def addTractionToElements(self, structProblem, elemIDs, tractions,
                              faceIndex=0, nastranOrdering=False):
        """
        The function is used to add a fixed traction to the
        selected element IDs. Tractions can be specified on an
        element by element basis (if tractions is a 2d array) or
        set to a uniform value (if tractions is a 1d array)

        Parameters
        ----------

        elemIDs : The components with added loads. Use selectCompIDs()
            to determine this.
        tractions : Numpy 1d or 2d array length varsPerNodes or (elemIDs, varsPerNodes)
            Array of traction vectors for each element
        faceIndex : int
            Indicates which face (side) of element to apply traction to.
            Note: not required for certain elements (i.e. shells)
        nastranOrdering : bool
            Flag signaling whether elemIDs are in TACS (default)
            or NASTRAN ordering
        """
        # Convert to tacs ordering for rest of procedure
        if nastranOrdering:
            globalElemIDs = self.idMap(elemIDs, self.nastranToTACSElemIDDict)
        else:
            globalElemIDs = elemIDs

        # Make sure the inputs are the correct shape
        globalElemIDs = numpy.atleast_1d(globalElemIDs)
        tractions = numpy.atleast_2d(tractions).astype(dtype=self.dtype)

        numElems = len(globalElemIDs)

        # If the user only specified one traction vector,
        # we assume the force should be the same for each element
        if tractions.shape[0] == 1:
            tractions = np.repeat(tractions, [numElems], axis=0)
        # If the dimensions still don't match, raise an error
        elif tractions.shape[0] != numElems:
            raise Error("Number of tractions must match number of elements,"
                        " {} tractions were specified for {} element IDs".format(tractions.shape[0], numElems))

        # First find the coresponding local element ID on each processor
        localElemIDs = self.getLocalElementIDsFromGlobal(globalElemIDs, nastranOrdering=False)

        # Set the structural problem
        self.setStructProblem(structProblem)

        # Flag to make sure we find all user-specified elements
        elemFound = np.zeros(numElems, dtype=int)

        # Loop through every element and if it's owned by this processor, add the traction
        for i, elemID in enumerate(localElemIDs):
            # The element was found on this proc
            if elemID >= 0:
                # Mark element as found
                elemFound[i] = 1
                # Get the pointer for the tacs element object for this element
                elemObjNum = self.elemObjectNumByElem[globalElemIDs[i]]
                elemObj = self.elemObjects[elemObjNum]
                # Create appropriate traction object for this element type
                tracObj = elemObj.createElementTraction(faceIndex, tractions[i])
                # Traction not implemented for element
                if tracObj is None:
                    TACSWarning("TACS element of type {} does not hav a traction implimentation. "
                                "Skipping element in addTractionToElement procedure.".format(elemObj.getObjectName()),
                                self.comm)
                # Traction implemented
                else:
                    # Add new traction to auxiliary element object
                    self.auxElems.addElement(elemID, tracObj)

        # Reduce the element flag and make sure that every element was found on exactly 1 proc
        elemFound = self.comm.allreduce(elemFound, op=MPI.SUM)

        # Warn the user if any elements weren't found
        if nastranOrdering:
            orderString = 'Nastran'
        else:
            orderString = 'TACS'

        for i in range(numElems):
            if not elemFound[i]:
                TACSWarning("Can't add traction to element ID {} ({} ordering), element not found in model. "
                            "Double check BDF file.".format(elemIDs[i], orderString), self.comm)

    def addPressureToComponents(self, structProblem, compIDs, pressures,
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
        elemIDs = self.getGlobalElementIDsForComps(compIDs, nastranOrdering=False)
        # Add pressure element by element
        self.addPressureToElements(structProblem, elemIDs, pressures, faceIndex, nastranOrdering=False)

        # Write out a message of what we did:
        self._info("Added a fixed pressure of %s to %d components, "
                   "distributed over %d elements." % (
                       repr(pressures), len(compIDs), len(elemIDs)),
                   maxLen=80, box=True)

    def addPressureToElements(self, structProblem, elemIDs, pressures,
                              faceIndex=0, nastranOrdering=False):
        """
        The function is used to add a fixed presure to the
        selected element IDs. Pressures can be specified on an
        element by element basis (if pressures is an array) or
        set to a uniform value (if pressures is a scalar)

        Parameters
        ----------

        elemIDs : The components with added loads. Use selectCompIDs()
            to determine this.
        pressures : Numpy array length 1 or elemIDs
            Array of pressure values for each element
        faceIndex : int
            Indicates which face (side) of element to apply pressure to.
            Note: not required for certain elements (i.e. shells)
        nastranOrdering : bool
            Flag signaling whether elemIDs are in TACS (default)
            or NASTRAN ordering
        """
        # Convert to tacs ordering for rest of procedure
        if nastranOrdering:
            globalElemIDs = self.idMap(elemIDs, self.nastranToTACSElemIDDict)
        else:
            globalElemIDs = elemIDs

        # Make sure the inputs are the correct shape
        globalElemIDs = numpy.atleast_1d(globalElemIDs)
        pressures = numpy.atleast_1d(pressures)

        numElems = len(globalElemIDs)

        # If the user only specified one pressure,
        # we assume the force should be the same for each element
        if pressures.shape[0] == 1:
            pressures = np.repeat(pressures, [numElems], axis=0)
        # If the dimensions still don't match, raise an error
        elif pressures.shape[0] != numElems:
            raise Error("Number of pressures must match number of elements,"
                        " {} pressures were specified for {} element IDs".format(pressures.shape[0], numElems))

        # First find the coresponding local element ID on each processor
        localElemIDs = self.getLocalElementIDsFromGlobal(globalElemIDs, nastranOrdering=False)

        # Set the structural problem
        self.setStructProblem(structProblem)

        # Flag to make sure we find all user-specified elements
        elemFound = np.zeros(numElems, dtype=int)

        # Loop through every element and if it's owned by this processor, add the pressure
        for i, elemID in enumerate(localElemIDs):
            # The element was found on this proc
            if elemID >= 0:
                elemFound[i] = 1
                # Get the pointer for the tacs element object for this element
                elemObjNum = self.elemObjectNumByElem[globalElemIDs[i]]
                elemObj = self.elemObjects[elemObjNum]
                # Get the pointer for the tacs element object for this element
                elemObjNum = self.elemObjectNumByElem[globalElemIDs[i]]
                elemObj = self.elemObjects[elemObjNum]
                # Create appropriate pressure object for this element type
                pressObj = elemObj.createElementPressure(faceIndex, pressures[i])
                # Pressure not implemented for element
                if pressObj is None:
                    TACSWarning("TACS element of type {} does not hav a pressure implimentation. "
                                "Skipping element in addPressureToElement procedure.".format(elemObj.getObjectName()),
                                self.comm)
                # Pressure implemented
                else:
                    # Add new pressure to auxiliary element object
                    self.auxElems.addElement(elemID, pressObj)


        # Reduce the element flag and make sure that every element was found on exactly 1 proc
        elemFound = self.comm.allreduce(elemFound, op=MPI.SUM)

        # Warn the user if any elements weren't found
        if nastranOrdering:
            orderString = 'Nastran'
        else:
            orderString = 'TACS'

        for i in range(numElems):
            if not elemFound[i]:
                TACSWarning("Can't add pressure to element ID {} ({} ordering), element not found in model. "
                            "Double check BDF file.".format(elemIDs[i], orderString), self.comm)

    def createTACSProbsFromBDF(self):
        """
        Automatically define tacs problem class using information contained in BDF file.
        This function assumes all loads are specified in the BDF and allows users to
        skip setting loads in Python.
        NOTE: Currently only supports FORCE, MOMENT, and PLOAD2 cards.
        NOTE: Currently only supports staticProblem (SOL 101)
        """
        # Make sure cross-referencing is turned on in pynastran
        if self.bdfInfo.is_xrefed is False:
            self.bdfInfo.cross_reference()
            self.bdfInfo.is_xrefed = True

        vpn = self.getOption('varsPerNode')
        loads = self.bdfInfo.loads
        nloads = len(loads)

        # Check if any loads are in the BDF
        if nloads == 0:
            raise Error("BDF file '%s' has no loads included in it. " % (self.bdfName))

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
            sp = tacs.problems.StaticProblem(name=name)

            if 'LOAD' in subCase.params:
                loadsID = subCase.params['LOAD'][0]
                # Check if there are any loads in this set
                if loadsID in loads:
                    loadSet = loads[loadsID]

                    loadDict = {}
                    for load in loadSet:
                        # Record any point force or moment cards
                        if load.type == 'FORCE' or load.type == 'MOMENT':
                            nodeID = load.node_ref.nid

                            if nodeID not in loadDict:
                                loadDict[nodeID] = numpy.zeros(vpn)

                            if load.type == 'FORCE' and vpn >= 3:
                                loadDict[nodeID][:3] += load.scaled_vector
                            elif load.type == 'MOMENT' and vpn >= 6:
                                loadDict[nodeID][3:6] += load.scaled_vector

                        # Add any pressure loads
                        if load.type == 'PLOAD2':
                            elemIDs = load.eids
                            pressure = load.pressure
                            self.addPressureToElements(sp, elemIDs, pressure, nastranOrdering=True)

                    # Add all point forces/moments for this problem at once
                    if len(loadDict) > 0:
                        self.addLoadToNodes(sp, list(loadDict.keys()), list(loadDict.values()), nastranOrdering=True)
            # append to list of structural problems
            structProblems[subCase.id] = sp

        return structProblems

    ####### Static solver methods ########

    def reset(self, SP):
        """ Reset each of the solution to last converged value."""
        self.setStructProblem(SP)
        self.u.copyValues(self.u_old)

    def _initializeSolve(self):
        """
        Initialze the solution of the structural system for the
        loadCase. The stiffness matrix is assembled and factored.
        """

        if self._factorOnNext:
            self.assembler.assembleJacobian(self.alpha, self.beta, self.gamma, self.res, self.K)
            self.PC.factor()
            self.old_update.zeroEntries()
            self._factorOnNext = False
            self._PCfactorOnNext = False

    def _updateStiffness(self):
        """
        The stiffness matrix is updated, assembled, and potentially factored.
        The residual is also computed for the current nonlinear iteration.
        """
        self.assembler.assembleMat(self.curLoadCase, self.K, self.res)
        self.PCjustFactored = False
        if self._PCfactorOnNext:
            self.PC.factor()
            self.PCjustFactored = True
            self._PCfactorOnNext = False

    def _updateResidual(self, Res, FullRes=None):
        """Compute the residual at the current loadScale and store it in Res and,
        if requested, compute the residual at the max loadScale and store it in FullRes.

        Everything except the aero forces are included in the residual assembly method, assembleRes,
        we then add the negative aero forces (stored in F) to compute the true residual.

        During a nonlinear analysis where load incrementation is employed, Res is the vector used
        to compute the next Newton step and to judge the convergence of each increment. FullRes,
        on the other hand, can be used as a measure of how close the current solution is to the
        solution at the max load scale.

        In other words:
        Res = `self.structure.assembleRes` + `self.loadScale`*`self.F`
            = -(Internal Forces + loadScale*Traction/Point/Inertial Loads etc) + loadScale*(-Aero Forces)
            = -(Internal Forces) - loadScale*((Traction/Point/Inertial Loads) + (Aero Forces))

        FullRes = `self.structure.assembleRes` + `self.maxLoadScale`*`self.F`
            = -(Internal Forces + maxLoadScale*Traction/Point/Inertial Loads etc) + maxLoadScale*(-Aero Forces)
            = -(Internal Forces) - maxLoadScale*((Traction/Point/Inertial Loads) + (Aero Forces))

        Parameters
        ----------
            Res : TACSVec
                The vector object in which to store the residual for the currently set load scale
            FullRes : TACSVec
                The vector object in which to store the residual computed wrt the max load
        """

        self.assembler.assembleRes(self.curLoadCase, Res)
        Res.axpy(-self.loadScale, self.F)

        if FullRes is not None:
            # Compute the residual at the max load scale
            self.assembler.setLoadFactor(self.curSP.tacsData.lcnum, self.maxLoadScale)
            self.assembler.assembleRes(self.curLoadCase, FullRes)
            FullRes.axpy(-self.maxLoadScale, self.F)
            # Reset the load scale
            self.assembler.setLoadFactor(self.curSP.tacsData.lcnum, self.loadScale)

    def _getForces(self, Fext=None, Fint=None):
        """Compute and store the internal and external force vectors.

        This method splits the residual and combines with the aero forces in order to
        obtain the full internal and external force vectors.

        These computations are based on the definition of the residual vector used in tacs:

        `self.res` = `self.structure.assembleRes` - `self.loadScale`*`self.F`
                   = -(Internal Forces + loadScale*Traction/Point/Inertial Loads etc) + loadScale*(-Aero Forces)
                   = -(Internal Forces) - loadScale*((Traction/Point/Inertial Loads) + (Aero Forces))
                   = -Fint - loadScale*Fext

        .. note:: This method works on the assumption that self.res contains the up to date
        residual at the current load scale

        Parameters
        ----------
            Fext : TACSVec
                The vector object in which to store the external nodal forces
            Fint : TACSVec
                The vector object in which to store the internal nodal forces
        """

        if Fext is not None:
            # assembleResLoadFactor returns the derivative of assembleRes with respect to the load factor.
            # Based on the definition above, this derivative is simply:

            # d(self.structure.assembleRes)/d(loadFactor) = d/d(loadFactor)(-(Internal Forces + loadScale*Traction/Point/Inertial Loads etc))
            # = -Traction/Point/Inertial Loads etc

            # We then combine with the negative aero loads and reverse to get the full external loads vector
            self.assembler.assembleResLoadFactor(self.curLoadCase, Fext)
            Fext.axpy(-1.0, self.F)
            Fext.scale(-1.0)

        if Fint is not None:
            # Compute internal loads by recomputing the residual with a load factor of zero:
            # `self.res` = -Fint - loadScale*Fext = -Fint - 0.0*Fext = -Fint
            LF = self.assembler.getLoadFactor(self.curLoadCase)
            self.assembler.setLoadFactor(self.curLoadCase, 0.0)
            self.assembler.assembleRes(self.curLoadCase, Fint)
            Fint.scale(-1.0)
            # Reset the loadFactor
            self.assembler.setLoadFactor(self.curLoadCase, LF)

    def __call__(self, structProblem, damp=1.0, useAitkenAcceleration=False,
                 dampLB=0.2, loadScale=1.0):
        """
        Solution of the structural system for loadCase. The
        forces must already be set.

        Parameters
        ----------
        structProblem

        Optional Arguments:

        damp, float: Value to use to damp the solution update. Default is 1.0

        useAitkenAcceleration, boolen: Flag to use
        aitkenAcceleration. Only applicable for aerostructural
        problems. Default is False.

        loadScale, float: value to scale external loads by. Only useful for
        load step approach on nonlinear problems.
        """
        startTime = time.time()

        self.setStructProblem(structProblem)
        self.curSP.tacsData.callCounter += 1

        # Set loadScale attributes, during load incrementation, self.loadScale is the current loadScale
        # while self.maxLoadScale is the target/final load scale.
        # For now, maxLoadScale is set equal to self.loadScale to make _updateResidual
        # and _getForces work, this will be addressed in future when new NL solver is merged
        self.loadScale = loadScale
        self.maxLoadScale = loadScale

        setupProblemTime = time.time()

        # Check if we need to initialize
        self._initializeSolve()

        initSolveTime = time.time()

        # Compute the RHS
        # TODO: Auxiliary forces still need to be load scaled
        # self.structure.setLoadFactor(self.curSP.tacsData.lcnum,loadScale)
        self.assembler.assembleRes(self.res)
        # Zero out bc terms in F
        self.assembler.applyBCs(self.F)
        # Add the -F
        self.res.axpy(-loadScale, self.F)

        # Set initnorm as the norm of F
        self.initNorm = numpy.real(self.F.norm()) * loadScale

        # Starting Norm for this compuation
        self.startNorm = numpy.real(self.res.norm())

        initNormTime = time.time()

        # TODO: Fix nonlinear solver
        # Check if user requested nonlinear analysis
        if self.getOption('solutionType').lower() != 'linear':
            # Copy current displacement vector to compute the update delta when
            # the solve is complete
            self.u_old.copyValues(self.u)
            nIter = self.getOption('nNLIter')
            for Iter in range(nIter):

                # If full Newton-Raphson requested, do that
                if self.getOption('NLSolver').upper() in ['FNR']:
                    self._FNRIter(loadScale)

                # We can add new nonlinear solvers here as well
                else:
                    TACSWarning('Unrecognized nonlinear solver, %s,requested by'
                                + 'user. Currently only FNR has been implemented. Using'
                                + 'FNR instead.' % (self.getOption('NLSolver')),
                                self.comm)
                    self.setOption('NLSolver', 'FNR')
                    self._FNRIter(loadScale)

                # Do Convergence Check
                self._convergenceCheck(Iter)

                # If we're converged, compute update delta relative to previous
                # displacement and exit
                if self.converged:
                    self.update.copyValues(self.u_old)
                    self.update.axpy(-1.0, self.u)
                    self.u.copyValues(self.u_old)
                    self.assembler.setVariables(self.curLoadCase, self.u)
                    break

                # If we failed, reset solution and exit
                if self.curSP.fatalFail:
                    self.reset(structProblem)
                    break
            self.finalIter = Iter

        # Otherwise, the system is linear, so perform one linear solve
        else:
            # Solve Linear System for the update
            self.KSM.solve(self.res, self.update)

        self.update.scale(-1.0)

        solveTime = time.time()

        # Apply Aitken Acceleration if necessary:
        if useAitkenAcceleration:
            if self.doDamp:
                # Comput: temp0 = update - old_update
                self.temp0.zeroEntries()
                self.temp0.axpy(1.0, self.update)
                self.temp0.axpy(-1.0, self.old_update)
                dnom = self.temp0.dot(self.temp0)
                damp = damp * (1.0 - self.temp0.dot(self.update) / dnom)

                # Clip to a reasonable range
                damp = numpy.clip(damp, dampLB, 1.0)
            self.doDamp = True

        # Update State Variables
        self.assembler.getVariables(self.u)
        self.u.axpy(damp, self.update)
        self.assembler.setVariables(self.u)

        # Set the old update
        self.old_update.copyValues(self.update)

        stateUpdateTime = time.time()

        # Compute final FEA Norm
        self.assembler.assembleRes(self.res)
        self.res.axpy(-loadScale, self.F)  # Add the -F
        self.finalNorm = numpy.real(self.res.norm())

        finalNormTime = time.time()

        # If timing was was requested print it, if the solution is nonlinear
        # print this information automatically if prinititerations was requested.
        if (self.getOption('printTiming') or (self.getOption('printIterations')
                                              and self.getOption('solutionType').lower() != 'linear')):
            self.pp('+--------------------------------------------------+')
            self.pp('|')
            self.pp('| TACS Solve Times:')
            self.pp('|')
            self.pp('| %-30s: %10.3f sec' % ('TACS Setup Time', setupProblemTime - startTime))
            self.pp('| %-30s: %10.3f sec' % ('TACS Solve Init Time', initSolveTime - setupProblemTime))
            self.pp('| %-30s: %10.3f sec' % ('TACS Init Norm Time', initNormTime - initSolveTime))
            self.pp('| %-30s: %10.3f sec' % ('TACS Solve Time', solveTime - initNormTime))
            self.pp('| %-30s: %10.3f sec' % ('TACS State Update Time', stateUpdateTime - solveTime))
            self.pp('| %-30s: %10.3f sec' % ('TACS Final Norm Time', finalNormTime - stateUpdateTime))
            self.pp('|')
            self.pp('| %-30s: %10.3f sec' % ('TACS Total Solution Time', finalNormTime - startTime))
            self.pp('+--------------------------------------------------+')

        return damp

    def loadStep(self, structProblem, nLoadSteps=10):
        """
        Perform a series of nonlinear solve while ramping up loads to
        desired level. This method can help convergence for otherwise difficult
        nonlinear structural problems.

        Parameters
        ----------
        structProblem

        Optional Arguments:

        nLoadSteps, int: Maximum number of load steps to take to reach final
                        load.
                        (Only useful for converging a nonlinear problem)
        """

        structProblem.tacsData.u.zeroEntries()
        dLS = 1.0 / (nLoadSteps)
        self.finalIter = 0
        loadScale = 0.0
        self.pp("+----------------------------+")
        self.pp("|Begining load step procedure|")
        self.pp("+----------------------------+")

        for j in range(nLoadSteps):
            # Don't check for anything on first iteration
            if self.finalIter == 0:
                pass
            # If the previous step did not converge lets take a step back and
            # try again with a smaller step size
            elif not self.converged:
                structProblem.tacsData.u.copyValues(self.u_old)
                self.assembler.setVariables(structProblem.tacsData.lcnum, self.u_old)
                loadScale -= dLS
                dLS *= 0.5
                self._PCfactorOnNext = True

            # If the previous nonlinear solve ran quickly, lets be a little
            # more aggresive with our step size for next iteration
            elif self.finalIter <= 6:
                dLS *= 1.25
            # If the previous nonlinear solve ran slowly, lets be a little
            # less aggresive with our step size for next iteration
            elif self.finalIter >= 10:
                dLS *= 0.75

            # Set next load step, clip it to 1.0 if it exceeds target value
            loadScale = min(loadScale + dLS, 1.0)

            self.pp("Load factor is %f" % (loadScale))

            self.__call__(structProblem, loadScale=loadScale)
            if loadScale == 1.0 and self.converged:
                self.pp("Load step converged!")
                break

    def _FNRIter(self, loadScale):
        """
        Function to conduct 1 iteration of the full Newton-Raphson
        """

        # Update stiffness matrix
        self._updateStiffness()
        # Compute the RHS
        self.assembler.setLoadFactor(self.curSP.tacsData.lcnum, loadScale)
        self.assembler.assembleRes(self.curLoadCase, self.res)
        self.res.axpy(-loadScale, self.F)  # Add the -F

        # Compute the current starting norm
        self.curNorm = numpy.real(self.res.norm())

        # Solve Linear System for the update
        self.KSM.solve(self.res, self.update)
        # Record number of subiterations to solve linear system
        nLinIters = self.KSM.getIterCount() + 1
        # If the solve takes too long, refactor preconitioner on next iter
        if nLinIters >= self.getOption('MaxItersPCReset'):
            self._PCfactorOnNext = True

        self.update.scale(-1.0)

        # Update State Variables
        self.assembler.getVariables(self.curLoadCase, self.u)
        self.u.axpy(1.0, self.update)
        self.assembler.setVariables(self.curLoadCase, self.u)

    def _convergenceCheck(self, Iter):
        """ Internal function to check nonlinear structural
        convergence"""

        # Get the requried data
        data = self._getMonitorData()

        self.curSP.fatalFail = False
        self.curSP.solveFailed = False
        self.converged = False

        # Only do check on root:
        if self.comm.rank == 0:
            if self.getOption('printIterations'):
                # Write header on first call
                if Iter == 0:
                    self._writeSolveHeader()

                # Write solution Data
                self._writeSolveData(Iter, data)

            # Check for Structural Convergence
            structChk = data['curResNorm'] <= (
                    self.getOption('L2ConvergenceRel') * data['startResNorm']) or \
                        data['curResNorm'] <= (self.getOption('L2Convergence'))

            # Check for divergence
            structDiv = data['curResNorm'] >= min(
                1. / self.getOption('L2ConvergenceRel') * data['startResNorm'],
                self.getOption('resNormUB'))

            if structChk:
                if self.getOption('printIterations'):
                    self.pp('--------------------------------------')
                    self.pp('      Nonlinear FEA CONVERGED!        ')
                    self.pp('--------------------------------------')
                self.converged = True

            if structDiv:
                self.curSP.fatalFail = True

            # Check for 'close' convergence
            if Iter == self.getOption('nNLIter') - 1 and not self.converged:
                # Check for "Close" Convergence at last iteration
                structChk = data['curResNorm'] <= (
                        self.getOption('L2ConvergenceRel') * 100 * data['startResNorm'])

                if structChk:
                    if self.getOption('printIterations'):
                        self.pp('--------------------------------------')
                        self.pp('      Nonlinear FEA SEMI-CONVERGED!   ')
                        self.pp('--------------------------------------')
                    self.converged = True
                else:
                    self.curSP.solveFailed = True
        # end if (rank == 0)

        # Broadcast flags from aero root proc to all
        self.curSP.fatalFail = self.comm.bcast(self.curSP.fatalFail, root=0)
        self.converged = self.comm.bcast(self.converged, root=0)

    def _getMonitorData(self):
        """This is an extra private function that retrieves and
        communicates the residual norm data, and a few other bits for
        the convergence check function and output"""
        data = {}  # struct data dictionary
        data['initResNorm'], data['startResNorm'], data['finalResNorm'] = (
            self.getResNorms())
        data['curResNorm'] = self.curNorm
        data['norm_u'] = self.u.norm()
        if self.PCjustFactored:
            data['PCFactored'] = '       *'
        else:
            data['PCFactored'] = ''
        return data

    def _writeSolveHeader(self):
        """Output the solver heading to the screen"""

        # We always have the "Iter" column, res, and u_norm.
        if self.dtype != 'D':
            headerStr = '%s  %17s  %17s %s' % (
                'Iter', 'Res', 'u_norm', 'Factored')
        else:
            headerStr = '%s  %17s %35s  %s' % (
                'Iter', 'Res', 'u_norm', 'Factored')
        print(headerStr)

    def _writeSolveData(self, Iter, data):
        """Write solve data to screen. Only called on one processor"""
        #  We always have the "Iter" column, res, and u_norm
        if self.dtype != 'D':
            dataStr = '%4.3d  %17.10e  %17.10e %s' % (
                Iter, data['curResNorm'], data['norm_u'], data['PCFactored'])

        else:
            dataStr = '%4.3d %17.10e% 17.10e   %17.10e% 17.10e %s' % (
                Iter,
                numpy.real(data['curResNorm']),
                numpy.imag(data['curResNorm']), numpy.real(data['norm_u']),
                numpy.imag(data['norm_u']), data['PCFactored'])

        print(dataStr)

    ####### Function eval/sensitivity methods ########

    def evalFunctions(self, structProblem, funcs, evalFuncs=None,
                      ignoreMissing=False):
        """
        This is the main routine for returning useful information from
        pytacs. The functions corresponding to the strings in
        EVAL_FUNCS are evaluated and updated into the provided
        dictionary.

        Parameters
        ----------
        structProblem : pyStructProblem class
            Structural problem to get the solution for
        funcs : dict
            Dictionary into which the functions are saved.
        evalFuncs : iterable object containing strings.
            If not none, use these functions to evaluate.
        ignoreMissing : bool
            Flag to supress checking for a valid function. Please use
            this option with caution.

        Examples
        --------
        >>> funcs = {}
        >>> FEAsolver(sp)
        >>> FEAsolver.evalFunctions(sp, funcs, ['mass'])
        >>> funcs
        >>> # Result will look like (if structProblem, sp, has name of 'c1'):
        >>> # {'cl_mass':12354.10}
        """
        startTime = time.time()
        # Set the structural problem
        self.setStructProblem(structProblem)
        if evalFuncs is None:
            evalFuncs = sorted(list(self.curSP.evalFuncs))
        else:
            evalFuncs = sorted(list(evalFuncs))

        if not ignoreMissing:
            for f in evalFuncs:
                if not f in self.functionList:
                    raise Error("Supplied function '%s' has not been added "
                                "using addFunction()." % f)

        setupProblemTime = time.time()

        # Fast parallel function evaluation of structural funcs:
        handles = [self.functionList[f] for f in evalFuncs if
                   f in self.functionList]
        funcVals = self.assembler.evalFunctions(handles)

        functionEvalTime = time.time()

        # Assign function values to appropriate dictionary
        i = 0
        for f in evalFuncs:
            if f in self.functionList:
                key = self.curSP.name + '_%s' % f
                self.curSP.funcNames[f] = key
                funcs[key] = funcVals[i]
                i += 1

        dictAssignTime = time.time()

        if self.getOption('printTiming'):
            self.pp('+--------------------------------------------------+')
            self.pp('|')
            self.pp('| TACS Function Times:')
            self.pp('|')
            self.pp('| %-30s: %10.3f sec' % ('TACS Function Setup Time', setupProblemTime - startTime))
            self.pp('| %-30s: %10.3f sec' % ('TACS Function Eval Time', functionEvalTime - setupProblemTime))
            self.pp('| %-30s: %10.3f sec' % ('TACS Dict Time', dictAssignTime - functionEvalTime))
            self.pp('|')
            self.pp('| %-30s: %10.3f sec' % ('TACS Function Time', dictAssignTime - startTime))
            self.pp('+--------------------------------------------------+')

    def evalFunctionsSens(self, structProblem, funcsSens, evalFuncs=None):
        """
        This is the main routine for returning useful (sensitivity)
        information from pytacs. The derivatives of the functions
        corresponding to the strings in EVAL_FUNCS are evaluated and
        updated into the provided dictionary.

        Parameters
        ----------
        structProblem : pyStructProblem class
            Structural problem to get the solution for
        funcsSens : dict
            Dictionary into which the derivatives are saved.
        evalFuncs : iterable object containing strings
            The functions the user wants returned

        Examples
        --------
        >>> funcsSens = {}
        >>> FEAsolver.evalFunctionsSens(sp, funcsSens, ['mass'])
        >>> funcs
        >>> # Result will look like (if structProblem, sp, has name of 'c1'):
        >>> # {'c1_mass':{'struct':[1.234, ..., 7.89]}
        """

        startTime = time.time()

        # Set the structural problem
        self.setStructProblem(structProblem)
        if evalFuncs is None:
            evalFuncs = sorted(list(self.curSP.evalFuncs))
        else:
            evalFuncs = sorted(list(evalFuncs))
        # Check that the functions are all ok.
        # and prepare tacs vecs for adjoint procedure
        dvSenses = []
        xptSenses = []
        dIdus = []
        adjoints = []
        for f in evalFuncs:
            if f not in self.functionList:
                raise Error("Supplied function has not beed added "
                            "using addFunction()")
            else:
                # Populate the lists with the tacs bvecs
                # we'll need for each adjoint/sens calculation
                dvSens = self.dvSensList[f]
                dvSens.zeroEntries()
                dvSenses.append(dvSens)

                xptSens = self.xptSensList[f]
                xptSens.zeroEntries()
                xptSenses.append(xptSens)

                dIdu = self.dIduList[f]
                dIdu.zeroEntries()
                dIdus.append(dIdu)

                adjoint = self.adjointList[f]
                adjoint.zeroEntries()
                adjoints.append(adjoint)

        setupProblemTime = time.time()

        adjointStartTime = {}
        adjointEndTime = {}

        # Next we will solve all the adjoints
        # Set adjoint rhs
        self.addSVSens(evalFuncs, dIdus)
        adjointRHSTime = time.time()
        for i, f in enumerate(evalFuncs):
            adjointStartTime[f] = time.time()
            self.solveAdjoint(dIdus[i], adjoints[i])
            adjointEndTime[f] = time.time()

        adjointFinishedTime = time.time()
        # Evaluate all the adoint res prooduct at the same time for
        # efficiency:
        self.addDVSens(evalFuncs, dvSenses)
        self.addAdjointResProducts(adjoints, dvSenses)
        self.addXptSens(evalFuncs, xptSenses)
        self.addAdjointResXptSensProducts(adjoints, xptSenses)

        # Recast sensititivities into dict for user
        for i, f in enumerate(evalFuncs):
            key = self.curSP.name + '_%s' % f
            # Finalize sensitivity arrays across all procs
            dvSenses[i].beginSetValues()
            dvSenses[i].endSetValues()
            xptSenses[i].beginSetValues()
            xptSenses[i].endSetValues()
            # Return sensitivities as array in sens dict
            funcsSens[key] = {self.varName: dvSenses[i].getArray().copy(),
                              self.coordName: xptSenses[i].getArray().copy()}

        totalSensitivityTime = time.time()

        if self.getOption('printTiming'):
            self.pp('+--------------------------------------------------+')
            self.pp('|')
            self.pp('| TACS Adjoint Times:')
            print('|')
            print('| %-30s: %10.3f sec' % ('TACS Sens Setup Problem Time', setupProblemTime - startTime))
            print('| %-30s: %10.3f sec' % (
                'TACS Adjoint RHS Time', adjointRHSTime - setupProblemTime))
            for f in evalFuncs:
                print('| %-30s: %10.3f sec' % (
                'TACS Adjoint Solve Time - %s' % (f), adjointEndTime[f] - adjointStartTime[f]))
            print('| %-30s: %10.3f sec' % ('Total Sensitivity Time', totalSensitivityTime - adjointFinishedTime))
            print('|')
            print('| %-30s: %10.3f sec' % ('Complete Sensitivity Time', totalSensitivityTime - startTime))
            print('+--------------------------------------------------+')

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
        else:
            raise ValueError("setDesignVars must be called with either a numpy array or dict as input.")

        # Set the variables in tacs, and the constriant objects
        self.assembler.setDesignVars(self.x)
        self._factorOnNext = True

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

        # Set the variables in tacs, and the constriant objects
        # Set the variables in tacs, and the constriant objects
        self.assembler.getDesignVars(self.x)
        return self.x.getArray().copy()

    def getNumDesignVars(self):
        """
        Return the number of design variables on this processor.
        """
        return self.x.getSize()

    def getTotalNumDesignVars(self):
        """
        Return the number of design variables across all processors.
        """
        return self.dvNum

    def getCoordinates(self):
        """
        Return the mesh coordiantes of the structure.

        Returns
        -------
        coords : array
            Structural coordinate in array of size (N, 3) where N is
            the number of structural nodes on this processor.
        """
        Xpts = self.assembler.createNodeVec()
        self.assembler.getNodes(Xpts)
        coords = Xpts.getArray()

        return coords

    def setCoordinates(self, coords):
        """
        Set the mesh coordinates of the structure.

        Returns
        -------
        coords : array
            Structural coordinate in array of size (N, 3) where N is
            the number of structural nodes on this processor.
        """
        XptsArray = self.Xpts.getArray()
        # Make sure input is raveled (1D) in case user changed shape
        XptsArray[:] = numpy.ravel(coords)
        self.assembler.setNodes(self.Xpts)
        self._factorOnNext = True

    def getNumCoordinates(self):
        """
        Return the number of mesh coordinates on this processor.
        """
        return self.Xpts.getSize()

    ####### Post processing methods ########

    def getVariablesAtPoints(self, structProblem, points):
        '''The function is used to get the state variables DOF's at the
        selected physical locations, points. A closest point search is
        used to determine the FE nodes that are the closest to the
        requested nodes.

        NOTE: The number and units of the entries of the state vector
        depends on the physics problem being solved and the dofs included
        in the model.

        A couple of examples of state vector components for common problem are listed below:

        In Elasticity with varsPerNode = 3,
        q = [u, v, w] # displacements
        In Elasticity with varsPerNode = 6,
        q = [u, v, w, tx, ty, tz] # displacements + rotations
        In Thermoelasticity with varsPerNode = 4,
        q = [u, v, w, T] # displacements + temperature
        In Thermoelasticity with varsPerNode = 7,
        q = [u, v, w, tx, ty, tz, T] # displacements + rotations + temperature
        '''
        try:
            from scipy.spatial import cKDTree
        except:
            raise Error("scipy.spatial "
                        "must be available to use getDisplacements")

        points = numpy.atleast_2d(points)
        self.setStructProblem(structProblem)

        # Pull out the local nodes on the proc and search "points" in the tree
        vpn = self.getOption('varsPerNode')
        Xpts = self.assembler.createNodeVec()
        self.assembler.getNodes(Xpts)
        localNodes = np.real(Xpts.getArray())
        nNodes = len(localNodes) // vpn
        xNodes = localNodes[0:nNodes * 3].reshape((nNodes, 3)).copy()
        tree = cKDTree(xNodes)

        d, index = tree.query(points, k=1)

        # Now figure out which proc has the best distance for this
        localu = np.real(structProblem.tacsData.u.getArray())
        uNodes = localu[0:nNodes * vpn].reshape((nNodes, vpn)).copy()
        u_req = numpy.zeros([len(points), vpn])
        for i in range(len(points)):
            proc = self.comm.allreduce((d[i], self.comm.rank), op=MPI.MINLOC)
            u_req[i, :] = uNodes[index[i], :]
            u_req[i, :] = self.comm.bcast(uNodes[index[i], :], root=proc[1])
        return u_req

    def writeForcesFile(self, structProblem, fileName, internal=False):

        """Write a solution file containing the nodal forces

        This routine can be used help debug the structural load
        vector. It works by writing out a typical .f5 results file
        with the nodal displacement values replaced by nodal forces.

        Parameters
        ----------
        structProblem : pyStructProblem class
            Structural problem to use.
        fileName : str
            Filename of output file. This will be a f5 file so should
            have a .f5 extension.
        internal : bool
            Flag to determine whether to write internal or external nodal forces.
            External forces are those applied to the structure to the user,
            Internal forces are those resulting from internal stresses caused
            by displacements in the structure, at equilibrium these are equal
            and opposite to the external forces.
        """

        self.setStructProblem(structProblem)

        self.u_tmp = self.assembler.createVec()
        self.assembler.getVariables(self.curLoadCase, self.u_tmp)

        self._updateResidual(self.res)
        if internal:
            self._getForces(Fext=None, Fint=self.temp0)
        else:
            self._getForces(Fext=self.temp0, Fint=None)

        # Write the forces as displacements
        self.assembler.setVariables(self.curLoadCase, self.temp0)
        self.writeOutputFile(fileName)

        # Reset the displacements
        self.assembler.setVariables(self.curLoadCase, self.u_tmp)

    def writeDVVisualization(self, fileName, n=17):
        """
        This function writes a standard f5 output file, but with
        design variables defined by x=mod(arange(nDV, n)), where n an
        integer supplied by the user. The idea is to use contouring in
        a post processing program to visualize the structural design
        variables.

        Parameters
        ----------
        fileName : str
            Filename to use. Since it is a f5 file, shoud have .f5 extension
        n : int
            Modulus value. 17 is the default which tends to work well.
            """

        nDVs = self.getNumDesignVars()

        # Save the current variables
        xSave =  self.getDesignVars()

        # Generate and set the 'mod' variables
        x = numpy.mod(numpy.arange(nDVs), n)
        self.setDesignVars(x)

        # Normal solution write
        self.writeOutputFile(fileName)

        # Reset the saved variables
        self.setDesignVars(xSave)

    def writeOutputFile(self, fileName):
        """Low-level command to write the current loadcase to a file

        Parameters
        ----------
        fileName : str
            Filename for output. Should have .f5 extension.
         """
        self.outputViewer.writeToFile(fileName)

    def writeSolution(self, outputDir=None, baseName=None, number=None):
        """This is a generic shell function that writes the output
        file(s).  The intent is that the user or calling program can
        call this function and pyTACS writes all the files that the
        user has defined. It is recommneded that this function is used
        along with the associated logical flags in the options to
        determine the desired writing procedure

        Parameters
        ----------

        outputDir : str or None
            Use the supplied output directory
        baseName : str or None
            Use this supplied string for the base filename. Typically
            only used from an external solver.
        number : int or None
            Use the user spplied number to index solution. Again, only
            typically used from an external solver
        """
        # Check input
        if outputDir is None:
            outputDir = self.getOption('outputDir')

        if baseName is None:
            baseName = self.curSP.name

        # If we are numbering solution, it saving the sequence of
        # calls, add the call number
        if number is not None:
            # We need number based on the provided number:
            baseName = baseName + '_%3.3d' % number
        else:
            # if number is none, i.e. standalone, but we need to
            # number solutions, use internal counter
            if self.getOption('numberSolutions'):
                baseName = baseName + '_%3.3d' % self.curSP.tacsData.callCounter

        # Unless the writeSolution option is off write actual file:
        if self.getOption('writeSolution'):
            base = os.path.join(outputDir, baseName) + '.f5'
            self.outputViewer.writeToFile(base)

        if self.getOption('writeBDF'):
            base = os.path.join(outputDir, baseName) + '.bdf'
            self.writeBDF(base)

    # =========================================================================
    #   The remainder of the routines should not be needed by a user
    #   using this class directly. However, many of the functions are
    #   still public since they are used by a solver that uses this
    #   class, i.e. an Aerostructural solver.
    # =========================================================================

    def getNumComponents(self):
        """
        Return number of components (property) groups found in bdf.
        """
        return self.nComp

    def solveAdjoint(self, rhs, phi, damp=1.0):
        """
        Solve the structural adjoint.

        Parameters
        ----------
        rhs : TACS BVec
            right hand side vector for adjoint solve
        phi : TACS BVec
            BVec into which the adjoint is saved
        damp : float
            A damping variable for adjoint update. Typically only used
            in multidisciplinary analysis
        """

        # First compute the residual
        self.K.mult(phi, self.res)
        self.res.axpy(-1.0, rhs)  # Add the -RHS

        # Starting Norm for this compuation
        self.startNorm = numpy.real(self.res.norm())

        # Solve Linear System
        zeroGuess = 0
        self.update.zeroEntries()
        self.KSM.solve(self.res, self.update, zeroGuess)

        # Update the adjoint vector with the (damped) update
        phi.axpy(-damp, self.update)

        # Compute actual final FEA Norm
        self.K.mult(phi, self.res)
        self.res.axpy(-1.0, rhs)  # Add the -RHS
        self.finalNorm = numpy.real(self.res.norm())

    def getNumVariables(self):
        """Return the number of degrees of freedom (states) that are
        on this processor

        Returns
        -------
        nstate : int
            number of states.
        """
        return self.u.getSize()

    def getVariables(self, structProblem, states=None):
        """Return the current state values for the current
        structProblem"""
        self.setStructProblem(structProblem)

        if states is None:
            states = self.u.getArray().copy()
        else:
            states[:] = self.u.getArray()

        return states

    def setVariables(self, structProblem, states):
        """ Set the structural states for current load case. Typically
        only used for aerostructural analysis

        Parameters
        ----------
        states : array
            Values to set. Must be the size of getNumVariables()
        """
        self.setStructProblem(structProblem)
        self.u.setValues(states)
        self.assembler.setVariables(self.u)

    def setAdjoint(self, adjoint, objective):
        """
        Set externally supplied adjoint values.

        Arguments:

        adjoint, float, tacs bvec: An array of size getNumVariables() to
        be copied to the structural adjoint variables
        """
        self.phi.copyValues(adjoint)
        if objective is not None:
            if not objective in self.curSP.tacsData.adjoints:
                self.curSP.tacsData.adjoints[objective] = (
                    self.assembler.createVec())
            self.curSP.tacsData.adjoints[objective].copyValues(adjoint)

    def getAdjoint(self, objective):
        """ Return the adjoint values for objective if they
        exist. Otherwise just return zeros"""

        if objective in self.curSP.tacsData.adjoints:
            vals = numpy.zeros(self.getNumVariables())
            self.curSP.tacsData.adjoints[objective].getValues(vals)
            return vals
        else:
            return numpy.zeros(self.getNumVariables(), self.dtype)

    def addSVSens(self, evalFuncs, dIduList):
        """ Add the state variable sensitivity to the ADjoint RHS for given evalFuncs"""
        funcHandles = [self.functionList[f] for f in evalFuncs if
                       f in self.functionList]
        self.assembler.addSVSens(funcHandles, dIduList, self.alpha, self.beta, self.gamma)

    def addDVSens(self, evalFuncs, dvSensList, scale=1.0):
        """ Add pratial sensitivity contribution due to design vars for evalFuncs"""
        funcHandles = [self.functionList[f] for f in evalFuncs if
                       f in self.functionList]
        self.assembler.addDVSens(funcHandles, dvSensList, scale)

    def addAdjointResProducts(self, adjointlist, dvSensList, scale=-1.0):
        """ Add the adjoint product contribution to the design variable sensitivity arrays"""
        self.assembler.addAdjointResProducts(adjointlist, dvSensList, scale)

    def addXptSens(self, evalFuncs, xptSensList, scale=1.0):
        """ Add pratial sensitivity contribution due to nodal coordinates for evalFuncs"""
        funcHandles = [self.functionList[f] for f in evalFuncs if
                       f in self.functionList]
        self.assembler.addXptSens(funcHandles, xptSensList, scale)

    def addAdjointResXptSensProducts(self, adjointlist, xptSensList, scale=-1.0):
        """ Add the adjoint product contribution to the nodal coordinates sensitivity arrays"""
        self.assembler.addAdjointResXptSensProducts(adjointlist, xptSensList, scale)

    def getResidual(self, structProblem, res=None, Fext=None):
        """
        This routine is used to evaluate directly the structural
        residual. Only typically used with aerostructural analysis.

        Parameters
        ----------
        structProblem : pyStructProblem class
            Structural problem to use
        res : numpy array
            If res is not None, place the residuals into this array.

        Returns
        -------
        res : array
            The same array if res was provided, (otherwise a new
            array) with evaluated residuals
        """
        self.setStructProblem(structProblem)
        self.assembler.assembleRes(self.res)
        self.res.axpy(1.0, self.curSP.tacsData.F)  # Add the -F

        if Fext is not None:
            resArray = self.res.getArray()
            resArray[:] -= Fext[:]

        if res is None:
            res = self.res.getArray().copy()
        else:
            res[:] = self.res.getArray()

        return res

    def getResNorms(self):
        """Return the initial, starting and final Res Norms. Note that
        the same norms are used for both solution and adjoint
        computations"""

        return (numpy.real(self.initNorm), numpy.real(self.startNorm),
                numpy.real(self.finalNorm))

    def zeroVectors(self):
        """Zero all the tacs b-vecs"""
        if not self._variablesCreated:
            self._createVariables()
        self.res.zeroEntries()
        self.u.zeroEntries()
        self.assembler.setVariables(self.u)
        self.update.zeroEntries()

    def setOption(self, name, value):
        """
        Set a solver option value. The name is not case sensitive.
        """
        name = name.lower()

        # Try to the option in the option dictionary
        defOptions = self.options['defaults']
        try:
            defOptions[name]
        except:
            TACSWarning('Option: \'%-30s\' is not a valid TACS option |' % name,
                        self.comm)
            return

        # Now we know the option exists, lets check if the type is ok:
        #        if type(value) == self.options[name][0]:
        if isinstance(value, self.options[name][0]):
            # Just set:
            self.options[name] = [type(value), value]
        else:
            raise Error("Datatype for Option %s was not valid. "
                        "Expected data type is %s. Received data type "
                        " is %s." % (name, self.options[name][0], type(value)))

    def getOption(self, name):

        # Redefine the getOption def from the base class so we can
        # mane sure the name is lowercase

        def_options = self.options['defaults']
        if name.lower() in def_options:
            return self.options[name.lower()][1]
        else:
            raise AttributeError(repr(name) + ' is not a valid option name')

    def _createOutputGroups(self):
        """Automatically determine how to split out the output file
        for easier viewing"""

        self.fam = []
        for i in range(self.nComp):
            aux = self.compDescript[i].split(self.getOption('familySeparator'))
            self.fam.append(aux[0])

        # Uniqify them and sort
        self.fam = sorted(numpy.unique(self.fam))

        self.compFam = numpy.zeros(self.nComp, dtype='intc')
        for i in range(self.nComp):
            aux = self.compDescript[i].split(self.getOption('familySeparator'))
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
        elif self.getOption('varsPerNode') == 6:
            elementType = tacs.TACS.BEAM_OR_SHELL_ELEMENT
        elif self.getOption('varsPerNode') == 3:
            elementType = tacs.TACS.SOLID_ELEMENT

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
                    TACSWarning('Trying to add component ID of %d, which\
                    is out of the range 0 <= compID < %d' % (item, self.nComp),
                                self.comm)

            elif isinstance(item, str):
                # This is a little inefficinet here; loop over
                # self.compDescript and see if 'item' (a string) in
                # part of the description. if so add.
                item = item.upper()
                for i in range(self.nComp):
                    if item in self.compDescript[i]:
                        compIDs[-1].append(i)
            else:
                TACSWarning('Unidentifiable information given for \'include\'\
                or \'exclude\'. Valid data are integers 0 <= i < %d, or \
                strings.' % self.nComp, self.comm)

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

    def _getCompIDsInBounds(self, includeBounds, projectVector=None):

        # This function Requires scipy spatital for doing the convex
        # hull in (either 2d or 3d)

        if self.rank == 0:
            try:
                from scipy.spatial import Delaunay
            except ImportError:
                raise Error("scipy.spatial could not be loaded. Selecting "
                            "components by bounds is not supported.")

            # Use default projectVector if not given:
            if projectVector is None:
                projectVector = numpy.array(self.getOption('projectVector'),
                                            dtype='d')
            else:
                projectVector = numpy.array(projectVector, dtype='d')

            # Individually select the compIDs for each entry in
            # includeBounds...this is recursive.
            compIDs = {}
            for item in includeBounds:
                compIDs[item] = self.selectCompIDs(include=item)[0]

            # Now we want to find the "intersection" between each of the
            # sets of componets in compIDs. I don't think there is a
            # better way than an n^2 algorithm for this.

            shared_edges = []
            keys = sorted(compIDs.keys())
            for i in range(len(keys)):  # |-> Double loop over each
                for j in range(i + 1, len(keys)):  # |   group of items
                    if i != j:
                        iComps = compIDs[keys[i]]  # | -> These are lists of the
                        jComps = compIDs[keys[j]]  # | -> compontIDs

                        # Now we want to find the edges shared between
                        # each possible set of components
                        for ii in iComps:  # | -> Double group over each
                            for jj in jComps:  # | pair of componets
                                if ii < jj:
                                    key = (ii, jj)
                                else:
                                    key = (jj, ii)

                                if key in self.dvEdges:
                                    # These are the shared edges.
                                    shared_edges.append(self.dvEdges[key])

            # That finishes the logical section. Now we proceed
            # geometrically. Get a unique set of nodes:
            uniqueNodes = numpy.array(numpy.unique(self._flatten(shared_edges)),
                                      dtype='intc')
            # Extract the coordinates for these
            nodes = self.getBDFNodes(uniqueNodes)

            # We first do a 3d bounding box (+ delta) elimination of
            # possible components. This is by far the fastest way to
            # eliminate components:
            xMin = numpy.min(nodes, axis=0)
            xMax = numpy.max(nodes, axis=0)
            # Expand the box by 10%
            d = numpy.linalg.norm(xMax - xMin) * 0.1
            xMin -= d
            xMax += d

            # Now we have to contruct a coordinate system in the plane
            # defined by projectVector.
            length = numpy.linalg.norm(projectVector)
            if length < 1e-12:
                raise Error("projectVector has almost zero length!")

            n = projectVector / length  # Plane (unit) normal
            pt = [1, 1, 1]
            if numpy.linalg.norm(numpy.cross(pt, n)) < 1e-12:
                pt = [2, 1, 1]

            v1 = pt - numpy.dot(pt, n) * n
            v1 /= numpy.linalg.norm(v1)
            v2 = numpy.cross(v1, projectVector)

            def projectNodesOntoPlane(nodes, projectVector):
                projNodes = numpy.zeros((len(nodes), 2))
                for i in range(len(nodes)):
                    # Project onto v1
                    projNodes[i, 0] = numpy.dot(nodes[i], v1)
                    projNodes[i, 1] = numpy.dot(nodes[i] -
                                                projNodes[i, 0] * v1, v2)
                return projNodes

            # Project to plane
            nodes = projectNodesOntoPlane(nodes, projectVector)

            # Create delaunay triagularization
            hull = Delaunay(nodes)

            # Simple loop over each component:
            compIDs = []
            for compID in range(self.nComp):
                cenPt = self.compCenters[compID]
                if (cenPt[0] > xMin[0] and cenPt[0] < xMax[0] and
                        cenPt[1] > xMin[1] and cenPt[1] < xMax[1] and
                        cenPt[2] > xMin[2] and cenPt[2] < xMax[2]):
                    # We are in the bounds, so check in hull:
                    proj = projectNodesOntoPlane([cenPt],
                                                 projectVector)[0]
                    res = hull.find_simplex(proj)
                    if res >= 0:
                        compIDs.append(compID)

        else:  # rank == 0
            compIDs = None

        # Broadcast and return
        return self.comm.bcast(compIDs, root=0)

    def _createElements(self, elemCallBack):
        """
        Create all the constitutive objects by calling the
        userSupplied or default callback function
        """

        for i in range(self.nComp):

            # Get a list of compDescripts to help the user
            compDescript = self.compDescript[i]
            numElements = len(self.elemDescripts[i])
            # TACS component ID
            compID = i
            # Nastran property ID
            propID = list(self.bdfInfo.property_ids)[i]

            # Call the user function
            result = elemCallBack(self.dvNum, compID, compDescript, self.elemDescripts[i], self.specialDVs,
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
                    TACSWarning("Could not identify objects returned \
                    from elemCallBack. Valid return objects are: \
                    A list of TACS element objects (required, first), \
                    an iterable object \
                    (eg, list or array) containing the scaling parameters \
                    for the added design variables (optional, second). The \
                    string representation of the offending object is: \
                    '%s'" % repr(result[1]), self.comm)

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
                        Error("Object of type %s returned in elemCallBack function "
                              "is not a valid TACS element object. The \
                               string representation of the offending object is: \
                               '%s'"%(type(object), repr(object)))

            if numFoundElements != numElements:
                raise Error("Could not find all required element objects in the "
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
                    raise Error("Inconsistent design variables detected. "
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
                        TACSWarning('An incorrect number of scale variables \
                        were returned. There were %d variables added, but only \
                        %d scale variables returned. The scale for these \
                        variables will be set to 1.0. The scale variables are %s.' % (
                            len(newVars), len(scaleList), repr(scaleList)),
                                    self.comm)
                        self.scaleList.extend(numpy.ones(len(newVars)))
                    else:
                        self.scaleList.extend(scaleList)

            # Loop through every element type in this component,
            # there may be multiple (e.g CQUAD4 + CTRIA3)
            for j, elemObject in enumerate(elemObjects):
                # Set each of the elements for this component
                pointer = self.elemObjectNumByComp[i][j]
                self.elemObjects[pointer] = elemObject

    def _createVariables(self):
        """Internal to create the variable required by TACS"""

        if not self._variablesCreated:

            self.x = self.assembler.createDesignVec()
            self.Xpts = self.assembler.createNodeVec()

            # Generic residual vector
            self.res = self.assembler.createVec()

            self.u_old = self.assembler.createVec()
            self.old_update = self.assembler.createVec()
            self.temp0 = self.assembler.createVec()

            # Current adjoint vector
            self.phi = self.assembler.createVec()

            # Current adjoint RHS
            self.adjRHS = self.assembler.createVec()

            # Current derivative of objective wrt states
            self.dIdu = self.assembler.createVec()

            # The 'pseudo load vector'. Contribution from aerodynamics
            # for coupled problems.
            self.pLoad = self.assembler.createVec()

            # Derivative of structural residuals wrt states, from
            # aerodynamics for coupled problems
            self.dSdu = self.assembler.createVec()
            opt = self.getOption

            # Tangent Stiffness --- process the ordering option here:
            tmp = opt('orderingType').lower()
            if tmp == 'natural':
                ordering = tacs.TACS.NATURAL_ORDER
            elif tmp == 'nd':
                ordering = tacs.TACS.ND_ORDER
            elif tmp == 'rcm':
                ordering = tacs.TACS.RCM_ORDER
            elif tmp == 'tacs_amd':
                ordering = tacs.TACS.TACS_AMD_ORDER
            elif tmp == 'multicolor':
                ordering = tacs.TACS.MULTICOLOR_ORDER
            else:
                raise Error("Unrecognized 'orderingType' option value: "
                            "'%s'. Valid values are: 'natural', 'nd', 'rcm', "
                            "'tacs_amd', or 'multicolor'." % tmp)

            self.K = self.assembler.createSchurMat(ordering)

            # Additional Vecs for updates
            self.update = self.assembler.createVec()

            # Setup PCScMat and KSM solver
            self.alpha = 1.0
            self.beta = 0.0
            self.gamma = 0.0
            self.assembler.assembleJacobian(self.alpha, self.beta, self.gamma, self.res, self.K)

            reorderSchur = 1
            self.PC = tacs.TACS.Pc(self.K, lev_fill=opt('PCFillLevel'),
                                   ratio_fill=opt('PCFillRatio'), reorder=reorderSchur)

            # Operator, fill level, fill ratio, msub, rtol, ataol
            if opt('KSMSolver').upper() == 'GMRES':
                self.KSM = tacs.TACS.KSM(
                    self.K, self.PC, opt('subSpaceSize'), opt('nRestarts'),
                    opt('flexible'))
            # TODO: Fix this
            # elif opt('KSMSolver').upper() == 'GCROT':
            #    self.KSM = tacs.TACS.GCROT(
            #        self.K, self.PC, opt('subSpaceSize'), opt('subSpaceSize'),
            #        opt('nRestarts'), opt('flexible'))
            else:
                raise Error("Unknown KSMSolver option. Valid options are "
                            "'GMRES' or 'GCROT'")

            self.KSM.setTolerances(self.getOption('L2ConvergenceRel'),
                                   self.getOption('L2Convergence'))

            if opt('useMonitor'):
                self.KSM.setMonitor(tacs.TACS.KSMPrintStdout(
                    opt('KSMSolver'), self.comm.rank, opt('monitorFrequency')))

            self._variablesCreated = True

    # ----------------------------------------------------------------------------
    #                      Utility Functions
    # ---------------------------------------------------------------------------
    def pp(self, printStr):
        """ Simple parallel print"""
        if self.comm.rank == 0:
            print(printStr)

    def _info(self, message, maxLen=80, box=False):
        """ Generic function for writing an info message. """

        if self.rank == 0:
            if not box:
                i = 9
                print('INFO: ', end="")
                aux = message.split()
                for word in aux:
                    if len(word) + i > 120:
                        print(' ')
                        print(' ' * 6, end="")
                        i = 0

                    print(word, end=" ")
                    i += len(word) + 1

                print()
            else:
                print('+' + '-' * (maxLen - 2) + '+')
                print('| INFO: ', end="")
                i = 9
                aux = message.split()
                for word in aux:
                    if len(word) + i > maxLen - 2:
                        print(' ' * (80 - i) + '|')
                        print('|', end="")
                        i = 2
                        print(word, end=" ")
                        i += len(word) + 1
                    else:
                        print(word, end=" ")
                        i += len(word) + 1

                print(' ' * (maxLen - i) + '|')
                print('+' + '-' * (maxLen - 2) + '+', )

    # Misc Functions
    def _flatten(self, l, ltypes=(list, tuple)):
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

    def print_scientific_8(self, value):
        """
        Prints a value in 8-character scientific notation.
        This is a sub-method and shouldnt typically be called

        See Also
        --------
        print_float_8 : for a better method
        """
        python_value = '%8.11e' % value
        (svalue, sExponent) = python_value.strip().split('e')
        exponent = int(sExponent)  # removes 0s

        sign = '-' if abs(value) < 0.01 else '+'

        sExp2 = str(exponent).strip('-+')  # the exponent will be added later...
        value2 = float(svalue)

        lenSExp = len(sExp2) + 1  # the plus 1 is for the sign
        leftover = 8 - lenSExp

        if value < 0:
            Format = "%%1.%sf" % (leftover - 3)
        else:
            Format = "%%1.%sf" % (leftover - 2)

        svalue3 = Format % value2
        svalue4 = svalue3.strip('0')
        field = "%8s" % (svalue4 + sign + sExp2)
        return field

    def print_float_8(self, value):
        """
        Prints a float in nastran 8-character width syntax using the
        highest precision possbile.
        """
        value = float(value)
        if value == 0.0:
            return '%8s' % '0.'
        elif value > 0.:  # positive, not perfect...
            if value < 5e-8:
                field = self.print_scientific_8(value)
                return field
            elif value < 0.001:
                field = self.print_scientific_8(value)
                field2 = "%8.7f" % value  # small value
                field2 = field2.strip('0 ')

                field1 = field.replace('-', 'e-')

                if field2 == '.':
                    return self.print_scientific_8(value)
                if len(field2) <= 8 and float(field1) == float(field2):
                    field = field2
                    field = field.strip(' 0')
            elif value < 0.1:
                field = "%8.7f" % value
            elif value < 1.:
                field = "%8.7f" % value  # same as before...
            elif value < 10.:
                field = "%8.6f" % value
            elif value < 100.:
                field = "%8.5f" % value
            elif value < 1000.:
                field = "%8.4f" % value
            elif value < 10000.:
                field = "%8.3f" % value
            elif value < 100000.:
                field = "%8.2f" % value
            elif value < 1000000.:
                field = "%8.1f" % value
            else:  # big value
                field = "%8.1f" % value
                if field.index('.') < 8:
                    field = '%8.1f' % round(value)
                    field = field[0:8]
                    assert '.' != field[0], field
                else:
                    field = self.print_scientific_8(value)
                return field
        else:
            if value > -5e-7:
                field = self.print_scientific_8(value)
                return field
            elif value > -0.01:  # -0.001
                field = self.print_scientific_8(value)
                field2 = "%8.6f" % value  # small value
                field2 = field2.strip('0 ')

                # get rid of the first minus sign, add it on afterwards
                field1 = '-' + field.strip(' 0-').replace('-', 'e-')

                if len(field2) <= 8 and float(field1) == float(field2):
                    field = field2.rstrip(' 0')
                    field = field.replace('-0.', '-.')

            elif value > -0.1:
                # -0.01 >x>-0.1...should be 5 (maybe scientific...)
                field = "%8.6f" % value
                field = field.replace('-0.', '-.')
            elif value > -1.:
                # -0.1  >x>-1.....should be 6, but the baseline 0 is kept...
                field = "%8.6f" % value
                field = field.replace('-0.', '-.')
            elif value > -10.:
                field = "%8.5f" % value  # -1    >x>-10
            elif value > -100.:
                field = "%8.4f" % value  # -10   >x>-100
            elif value > -1000.:
                field = "%8.3f" % value  # -100  >x>-1000
            elif value > -10000.:
                field = "%8.2f" % value  # -1000 >x>-10000
            elif value > -100000.:
                field = "%8.1f" % value  # -10000>x>-100000
            else:
                field = "%8.1f" % value
                if field.index('.') < 8:
                    field = '%7s.' % int(round(value, 0))
                    assert '.' != field[0], field
                else:
                    field = self.print_scientific_8(value)
                return field
        field = field.strip(' 0')
        field = '%8s' % field

        assert len(field) == 8, ('value=|%s| field=|%s| is not 8 characters '
                                 'long, its %s' % (value, field, len(field)))
        return field

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
        if isinstance(fromIDs, int):
            return tacsIDDict[fromIDs]

        toIDs = [None] * len(fromIDs)
        for i, id in enumerate(fromIDs):
            if id in tacsIDDict:
                toIDs[i] = tacsIDDict[id]
            else:
                toIDs[i] = -1

        return toIDs


class TACSLoadCase(object):
    """
    A container class for storing data related to a particular load case
    """

    def __init__(self):
        self.F = None
        self.u = None
        self.auxElems = None
        self.adjoints = {}
        self.callCounter = -1


class Error(Exception):
    """
    Format the error message in a box to make it clear this
    was a expliclty raised exception.
    """

    def __init__(self, message):
        msg = '\n+' + '-' * 78 + '+' + '\n' + '| pyTACS Error: '
        i = 15
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
            msg = '\n+' + '-' * 78 + '+' + '\n' + '| pyTACS Warning: '
            i = 17
            for word in message.split():
                if len(word) + i + 1 > 78:  # Finish line and start new one
                    msg += ' ' * (78 - i) + '|\n| ' + word + ' '
                    i = 1 + len(word) + 1
                else:
                    msg += word + ' '
                    i += len(word) + 1
            msg += ' ' * (78 - i) + '|\n' + '+' + '-' * 78 + '+' + '\n'
            print(msg)
