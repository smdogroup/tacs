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

# Instantiate FEASolver
structOptions = {}

bdfFile = os.path.join(os.path.dirname(__file__), 'slot.bdf')
# Load BDF file
FEASolver = pyTACS(bdfFile, comm, options=structOptions)
# Set up TACS Assembler
# Don't need a elemCallBack since property info exists in bdf
FEASolver.createTACSAssembler()

# Add Functions
FEASolver.addFunction('mass', functions.StructuralMass)
FEASolver.addFunction('ks_vmfailure', functions.KSFailure, safetyFactor=1.5,
                      ksWeight=100.0)

# ==============================================================================
# Setup static problem
# ==============================================================================
# Static problem
evalFuncs = ['mass', 'ks_vmfailure']
# Read in forces from BDF and create tacs static problems
SPs = FEASolver.createTACSProbsFromBDF()

# Solve state
for spID in SPs:
    FEASolver(SPs[spID])

# Evaluate functions
funcs = {}
for spID in SPs:
    FEASolver.evalFunctions(SPs[spID], funcs, evalFuncs=evalFuncs)

if comm.rank == 0:
    pprint(funcs)

funcsSens = {}
for spID in SPs:
    FEASolver.evalFunctionsSens(SPs[spID], funcsSens, evalFuncs=evalFuncs)
if comm.rank == 0:
    pprint(funcsSens)

FEASolver.writeSolution(outputDir=os.path.dirname(__file__))
