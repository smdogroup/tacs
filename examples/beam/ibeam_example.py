"""
This example models a cantilevered I-beam, with a tip shear load applied.
The parameters for the problem are given below:
    h_web = 1.0
    t_web = 0.01
    A_cap = 0.005
    I_beam = 1/12 * h_web^3 * t_web + 2 * A_cap * (h_web/2)^2 = 0.00333
    L_beam = 10.0
    V = 1000000.
The I-beam is modeled through a combination of two element types:
    Shell elements: To model the webs
    Beam elements (only axial stiffness): To model the caps
The tip deflection of the model is given by the following formula
    v_tip = V * L^3 / (3 * E * I) = 1.4285714285714286
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

bdfFile = os.path.join(os.path.dirname(__file__), "I_beam.bdf")
# Load BDF file
FEAAssembler = pyTACS(bdfFile, comm, options=structOptions)

# Material properties
rho = 2700.0  # density kg/m^3
E = 70.0e9  # Young's modulus (Pa)
nu = 0.3  # Poisson's ratio
ys = 270.0e6  # yield stress

# Web thickness
t = 0.01  # m
# Flange area
A = 0.005  # m^2

# Tip shear
V = 1000000.0

# Callback function used to setup TACS element objects and DVs
def elemCallBack(dvNum, compID, compDescript, elemDescripts, globalDVs, **kwargs):
    # Setup (isotropic) property and constitutive objects
    prop = constitutive.MaterialProperties(rho=rho, E=E, nu=nu, ys=ys)
    # For each element type in this component,
    # pass back the appropriate tacs element object
    elemList = []
    for descript in elemDescripts:
        if descript == "CQUAD4":
            con = constitutive.IsoShellConstitutive(prop, t=t, tNum=dvNum)
            # TACS shells are sometimes a little overly-rigid in shear
            # We can reduce this effect by decreasing the drilling regularization
            con.setDrillingRegularization(0.1)
            refAxis = np.array([1.0, 0.0, 0.0])
            transform = elements.ShellRefAxisTransform(refAxis)
            elem = elements.Quad4Shell(transform, con)
        elif descript == "CROD":
            # Shear corrections and bending stiffness are zero for pure axial members
            con = constitutive.BasicBeamConstitutive(prop, A=A, ky=0.0, kz=0.0)
            refAxis = np.array([0.0, 0.0, 1.0])
            transform = elements.BeamRefAxisTransform(refAxis)
            elem = elements.Beam2(transform, con)
        else:
            raise ValueError(f'Element type "{descript}" not recognized.')
        elemList.append(elem)
    return elemList


# Set up elements and TACS assembler
FEAAssembler.initialize(elemCallBack)

# ==============================================================================
# Setup static problem
# ==============================================================================
# Static problem

# Create a static problem with a simple z shear load at tip node
problem = FEAAssembler.createStaticProblem("I_Beam")
# Apply load at centroid of RBE3 at tip of beam
problem.addLoadToNodes(89, [0.0, V, 0.0, 0.0, 0.0, 0.0], nastranOrdering=True)
# Add some eval funcs
problem.addFunction("mass", functions.StructuralMass)
problem.addFunction("compliance", functions.Compliance)
problem.addFunction(
    "ks_vmfailure", functions.KSFailure, safetyFactor=1.5, ksWeight=100.0
)
problem.addFunction(
    "y_disp", functions.KSDisplacement, ksWeight=100.0, direction=[0.0, 1.0, 0.0]
)
problem.addFunction(
    "z_disp", functions.KSDisplacement, ksWeight=100.0, direction=[0.0, 0.0, 1.0]
)

# Solve state
problem.solve()

# Evaluate functions
funcs = {}
problem.evalFunctions(funcs)

if comm.rank == 0:
    pprint(funcs)

# Evaluate function sensitivities
funcsSens = {}
problem.evalFunctionsSens(funcsSens)
if comm.rank == 0:
    pprint(funcsSens)

# Write solution out
problem.writeSolution(outputDir=os.path.dirname(__file__))
