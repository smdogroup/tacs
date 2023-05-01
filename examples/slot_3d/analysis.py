"""
Structural model of a slot made from 3D tetrahedral elements.
The slot is constrained at its bolt holes and a set of point
forces are applied at the top of the slot.
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

# ==============================================================================
# Extension modules
# ==============================================================================
from tacs import functions, pyTACS

comm = MPI.COMM_WORLD

# Instantiate FEAAssembler
structOptions = {}

bdfFile = os.path.join(os.path.dirname(__file__), "slot.bdf")
# Load BDF file
FEAAssembler = pyTACS(bdfFile, comm, options=structOptions)
# Set up TACS Assembler
# Don't need a elemCallBack since property info exists in bdf
FEAAssembler.initialize()

# ==============================================================================
# Setup static problem
# ==============================================================================
# Static problem
evalFuncs = ["mass", "ks_vmfailure"]
# Read in forces from BDF and create tacs static problems
SPs = FEAAssembler.createTACSProbsFromBDF()

# Set up eval functions
for spID in SPs:
    SPs[spID].addFunction("mass", functions.StructuralMass)
    SPs[spID].addFunction(
        "ks_vmfailure", functions.KSFailure, safetyFactor=1.5, ksWeight=100.0
    )

# Solve state
for spID in SPs:
    SPs[spID].solve()

# Evaluate functions
funcs = {}
for spID in SPs:
    SPs[spID].evalFunctions(funcs, evalFuncs=evalFuncs)

if comm.rank == 0:
    pprint(funcs)

# Evaluate function sensitivities
funcsSens = {}
for spID in SPs:
    SPs[spID].evalFunctionsSens(funcsSens, evalFuncs=evalFuncs)
if comm.rank == 0:
    pprint(funcsSens)

SPs[spID].writeSolution(outputDir=os.path.dirname(__file__))
