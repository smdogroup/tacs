"""
This is similar to the crm.py example, but demonstrates how to solve problems through
the pytacs interface rather than using TACS directly.
This wingbox is a simplified version of the one of the University of Michigan uCRM-9.
We use a couple of pyTACS load generating methods to model various wing loads under cruise.
The script runs the structural analysis, evaluates the wing mass and von misses failure index
and computes sensitivities with respect to wingbox thicknesses and node xyz locations.
The sensitivities are then verified agains finite-difference or complex step approximations.
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
from tacs import functions, constitutive, elements, pyTACS, TACS

comm = MPI.COMM_WORLD

# Instantiate FEAAssembler
structOptions = {
    "printtiming": True,
}

bdfFile = os.path.join(os.path.dirname(__file__), "CRM_box_2nd.bdf")
FEAAssembler = pyTACS(bdfFile, options=structOptions, comm=comm)

# Material properties
rho = 2780.0  # density kg/m^3
E = 73.1e9  # Young's modulus (Pa)
nu = 0.33  # Poisson's ratio
ys = 324.0e6  # yield stress

# Shell thickness
t = 0.01  # m
tMin = 0.002  # m
tMax = 0.05  # m

# Callback function used to setup TACS element objects and DVs
def elemCallBack(dvNum, compID, compDescript, elemDescripts, globalDVs, **kwargs):
    # Setup (isotropic) property and constitutive objects
    prop = constitutive.MaterialProperties(rho=rho, E=E, nu=nu, ys=ys)
    # Set one thickness dv for every component
    con = constitutive.IsoShellConstitutive(prop, t=t, tNum=dvNum)

    # Define reference axis for local shell stresses
    if "SKIN" in compDescript:  # USKIN + LSKIN
        sweep = 35.0 / 180.0 * np.pi
        refAxis = np.array([np.sin(sweep), np.cos(sweep), 0])
    else:  # RIBS + SPARS + ENGINE_MOUNT
        refAxis = np.array([0.0, 0.0, 1.0])

    # For each element type in this component,
    # pass back the appropriate tacs element object
    elemList = []
    transform = elements.ShellRefAxisTransform(refAxis)
    for elemDescript in elemDescripts:
        if elemDescript in ["CQUAD4", "CQUADR"]:
            elem = elements.Quad4Shell(transform, con)
        elif elemDescript in ["CTRIA3", "CTRIAR"]:
            elem = elements.Tri3Shell(transform, con)
        else:
            print("Uh oh, '%s' not recognized" % (elemDescript))
        elemList.append(elem)

    # Add scale for thickness dv
    scale = [100.0]
    return elemList, scale


# Set up elements and TACS assembler
FEAAssembler.initialize(elemCallBack)

# ==============================================================================
# Setup static problem
# ==============================================================================
# Static problem
problem = FEAAssembler.createStaticProblem("cruise")

# Add TACS Functions
problem.addFunction("mass", functions.StructuralMass)
problem.addFunction(
    "ks_vmfailure", functions.KSFailure, safetyFactor=1.5, ksWeight=100.0
)

# Various methods for adding loads to structural problem:
# Let's model an engine load of 75kN weight, and 64kN thrust attached to spar
compIDs = FEAAssembler.selectCompIDs(
    ["WING_SPARS/LE_SPAR/SEG.16", "WING_SPARS/LE_SPAR/SEG.17"]
)
We = 75.0e3  # N
Te = 64.0e3  # N
problem.addLoadToComponents(compIDs, [-Te, 0.0, -We, 0.0, 0.0, 0.0], averageLoad=True)
# Next we'll approximate aerodynamic loads on upper/lower skin with a uniform traction
L = 3e3  # N/m^2
D = 150  # N/m^2
tracVec = np.array([D, 0.0, L])
compIDs = FEAAssembler.selectCompIDs(include="SKIN")
problem.addTractionToComponents(compIDs, tracVec)
# We can approximate fuel load by adding a pressure load to the lower skin
P = 2e3  # N/m^2
compIDs = FEAAssembler.selectCompIDs(include="L_SKIN")
problem.addPressureToComponents(compIDs, P)
# Add gravity load to all structural elements
g = np.array([0.0, 0.0, -9.81])  # m/s^2
problem.addInertialLoad(g)

# ==============================================================================
# Solve static problem
# ==============================================================================

# Solve structural problem
problem.solve()

# Evaluate functions
funcs = {}
problem.evalFunctions(funcs)
if comm.rank == 0:
    pprint(funcs)

# Solve adjoints and evaluate function sensitivities
funcsSens = {}
problem.evalFunctionsSens(funcsSens)
if comm.rank == 0:
    pprint(funcsSens)

# Write out solution
problem.writeSolution(outputDir=os.path.dirname(__file__))

# ==============================================================================
# Perform sensitivity check
# ==============================================================================
# Perform a fd/cs sensisitivity check on design  variable sensitivity
x_orig = problem.getDesignVars()

# Get number of design variables owned by this proc
ndvs = FEAAssembler.getNumDesignVars()
x_pert = np.random.rand(ndvs)

if TACS.dtype == complex:
    dh = 1e-50 * 1j
else:
    dh = 1e-6
x_new = x_orig + x_pert * dh

# Re-solve and evaluate function with new perturbed design variable
funcs_new = {}
problem.setDesignVars(x_new)
# Solve
problem.solve()
# Evaluate pertrubed functions
problem.evalFunctions(funcs_new)

# Loop through each function and compare sensitivities
for funcName in funcs:
    dfunc_approx = np.real((funcs_new[funcName] - funcs[funcName]) / dh)
    # Project sensitivity against perturbation vector
    dfunc_exact_local = np.real(np.dot(funcsSens[funcName]["struct"], x_pert))
    # The sens vector is distributed across multiple processors,
    # accumulate sensitivity contribution from each processor to get total sensitivity
    dfunc_exact = comm.allreduce(dfunc_exact_local)
    if comm.rank == 0:
        print("Func name:      ", funcName)
        print("FD (DVSens):      ", dfunc_approx)
        print("Result (DVSens):  ", dfunc_exact)
        print("Rel err (DVSens): ", (dfunc_exact - dfunc_approx) / dfunc_exact)

# Reset design variables
problem.setDesignVars(x_orig)

# Perform a fd/cs sensisitivity check on nodal coordinate sensitivity
xpts_orig = problem.getNodes()

# Get number of nodes owned by this proc
nnodes = FEAAssembler.getNumOwnedNodes()
xpts_pert = np.random.rand(3 * nnodes)

xpts_new = xpts_orig + xpts_pert * dh

# Re-solve and evaluate function with new perturbed node coordinates
problem.setNodes(xpts_new)
# Solve
problem.solve()
# Evaluate pertrubed functions
problem.evalFunctions(funcs_new)

# Loop through each function and compare sensitivities
for funcName in funcs:
    dfunc_approx = np.real((funcs_new[funcName] - funcs[funcName]) / dh)
    # Project sensitivity against perturbation vector
    dfunc_exact_local = np.real(np.dot(funcsSens[funcName]["Xpts"], xpts_pert))
    # The sens vector is distributed across multiple processors,
    # accumulate sensitivity contribution from each processor to get total sensitivity
    dfunc_exact = comm.allreduce(dfunc_exact_local)
    if comm.rank == 0:
        print("Func name:      ", funcName)
        print("FD (XptSens):      ", dfunc_approx)
        print("Result (XptSens):  ", dfunc_exact)
        print("Rel err (XptSens): ", (dfunc_exact - dfunc_approx) / dfunc_exact)
