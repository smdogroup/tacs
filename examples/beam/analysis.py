"""
6 noded beam model 1 meter long in x direction.
We apply a tip load in the z direction and clamp it at the root.
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
from mpi4py import MPI
import numpy as np

# ==============================================================================
# Extension modules
# ==============================================================================
from tacs import constitutive, elements, functions, pyTACS

comm = MPI.COMM_WORLD

# Instantiate FEAAssembler
structOptions = {}

bdfFile = os.path.join(os.path.dirname(__file__), "beam_model.bdf")
# Load BDF file
FEAAssembler = pyTACS(bdfFile, comm, options=structOptions)

# Material properties
rho = 2700.0  # density kg/m^3
E = 70.0e9  # Young's modulus (Pa)
nu = 0.3  # Poisson's ratio
ys = 270.0e6  # yield stress

# Shell thickness
A = 0.1  # m
Iz = 0.2  # m
Iy = 0.3  # m
J = 0.4

# Callback function used to setup TACS element objects and DVs
def elemCallBack(dvNum, compID, compDescript, elemDescripts, globalDVs, **kwargs):
    # Setup (isotropic) property and constitutive objects
    prop = constitutive.MaterialProperties(rho=rho, E=E, nu=nu, ys=ys)
    con = constitutive.BasicBeamConstitutive(
        prop, A=A, Iy=Iy, Iz=Iz, J=J, ky=1000, kz=1000
    )

    refAxis = np.array([0.0, 1.0, 0.0])

    # For each element type in this component,
    # pass back the appropriate tacs element object
    transform = elements.BeamRefAxisTransform(refAxis)
    elem = elements.Beam2(transform, con)
    return elem


# Set up elements and TACS assembler
FEAAssembler.initialize(elemCallBack)

# ==============================================================================
# Setup static problem
# ==============================================================================
# Static problem
evalFuncs = ["mass", "ks_vmfailure"]

# Create a static problem with a simple z shear load at tip node
problems = FEAAssembler.createTACSProbsFromBDF().values()
# Add some eval funcs
for problem in problems:
    problem.addFunction("mass", functions.StructuralMass)
    problem.addFunction(
        "ks_vmfailure", functions.KSFailure, safetyFactor=1.5, ksWeight=100.0
    )

for problem in problems:
    # Solve state
    problem.solve()

# Evaluate functions
funcs = {}
for problem in problems:
    problem.evalFunctions(funcs, evalFuncs=evalFuncs)

if comm.rank == 0:
    pprint(funcs)

# Evaluate function sensitivities
funcsSens = {}
for problem in problems:
    problem.evalFunctionsSens(funcsSens, evalFuncs=evalFuncs)
if comm.rank == 0:
    pprint(funcsSens)

# Write solution out
for problem in problems:
    problem.writeSolution(outputDir=os.path.dirname(__file__))
