"""
This wingbox is a simplified version of the one of the University of Michigan uCRM-9.
We use a couple of pyTACS load generating methods to model various wing loads under cruise.
The script runs the structural analysis, evaluates the wing mass and von misses failure index
and computes sensitivities with respect to wingbox thicknesses and node xyz locations.
"""
# ==============================================================================
# Standard Python modules
# ==============================================================================
from __future__ import print_function
import os

# ==============================================================================
# External Python modules
# ==============================================================================
from pprint import pprint
import numpy as np
from mpi4py import MPI

# ==============================================================================
# Extension modules
# ==============================================================================
from tacs import functions, constitutive, elements, pyTACS, problems

comm = MPI.COMM_WORLD

# Instantiate FEASolver
structOptions = {
    'printtiming':True,
}

bdfFile = os.path.join(os.path.dirname(__file__), 'CRM_box_2nd.bdf')
FEASolver = pyTACS(bdfFile, options=structOptions, comm=comm)

# Material properties
rho = 2780.0        # density kg/m^3
E = 73.1e9          # Young's modulus (Pa)
nu = 0.33           # Poisson's ratio
kcorr = 5.0/6.0     # shear correction factor
ys = 324.0e6        # yield stress

# Shell thickness
t = 0.01            # m
tMin = 0.002        # m
tMax = 0.05         # m

# Callback function used to setup TACS element objects and DVs
def elemCallBack(dvNum, compID, compDescript, elemDescripts, globalDVs, **kwargs):
    # Setup (isotropic) property and constitutive objects
    prop = constitutive.MaterialProperties(rho=rho, E=E, nu=nu, ys=ys)
    # Set one thickness dv for every component
    con = constitutive.IsoShellConstitutive(prop, t=t, tNum=dvNum)

    # Define reference axis for local shell stresses
    if 'SKIN' in compDescript: # USKIN + LSKIN
        sweep = 35.0 / 180.0 * np.pi
        refAxis = np.array([np.sin(sweep), np.cos(sweep), 0])
    else: # RIBS + SPARS + ENGINE_MOUNT
        refAxis = np.array([0.0, 0.0, 1.0])

    # For each element type in this component,
    # pass back the appropriate tacs element object
    elemList = []
    transform = elements.ShellRefAxisTransform(refAxis)
    for elemDescript in elemDescripts:
        if elemDescript in ['CQUAD4', 'CQUADR']:
            elem = elements.Quad4Shell(transform, con)
        elif elemDescript in ['CTRIA3', 'CTRIAR']:
            elem = elements.Tri3Shell(transform, con)
        else:
            print("Uh oh, '%s' not recognized" % (elemDescript))
        elemList.append(elem)

    # Add scale for thickness dv
    scale = [100.0]
    return elemList, scale

# Set up elements and TACS assembler
FEASolver.createTACSAssembler(elemCallBack)

# Add TACS Functions
FEASolver.addFunction('wing_mass', functions.StructuralMass)
FEASolver.addFunction('ks_vmfailure', functions.KSFailure, safetyFactor=1.5,
                      ksWeight=100.0)

# ==============================================================================
# Setup static problem
# ==============================================================================
# Static problem
evalFuncs = ['wing_mass', 'ks_vmfailure']
sp = problems.StaticProblem(name='wing', evalFuncs=evalFuncs, loadFactor=2.5)

# Various methods for adding loads to structural problem:
# Let's model an engine load of 75kN weight, and 64kN thrust attached to spar
compIDs = FEASolver.selectCompIDs(["WING_SPARS/LE_SPAR/SEG.16", "WING_SPARS/LE_SPAR/SEG.17"])
We = 75.0e3 # N
Te = 64.0e3 # N
FEASolver.addLoadToComponents(sp, compIDs, [-Te, 0.0, -We, 0.0, 0.0, 0.0], averageLoad=True)
# Next we'll approximate aerodynamic loads on upper/lower skin with a uniform traction
L = 3e3 # N/m^2
D = 150 # N/m^2
tracVec = np.array([D, 0.0, L])
compIDs = FEASolver.selectCompIDs(include='SKIN')
FEASolver.addTractionToComponents(sp, compIDs, tracVec)
# Finally, we can approximate fuel load by adding a pressure load to the lower skin
P = 2e3 # N/m^2
compIDs = FEASolver.selectCompIDs(include='L_SKIN')
FEASolver.addPressureToComponents(sp, compIDs, P)

# Solve structural problem
FEASolver(sp)

# Evaluate functions
funcs = {}
FEASolver.evalFunctions(sp, funcs)
if comm.rank == 0:
    pprint(funcs)

# Solve adjoints and evaluate function sensitivities
funcsSens = {}
FEASolver.evalFunctionsSens(sp, funcsSens)
if comm.rank == 0:
    pprint(funcsSens)

FEASolver.writeSolution(outputDir=os.path.dirname(__file__))