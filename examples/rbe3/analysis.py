"""
Cylindrical beam constructed from shell elements. The beam is cantilevered at
one end and loaded at the other using an RBE3.
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
from tacs import pyTACS, functions

comm = MPI.COMM_WORLD

bdfFile = os.path.join(os.path.dirname(__file__), "rbe3.bdf")
# Load BDF file
FEAAssembler = pyTACS(bdfFile, comm=comm)
# Set up TACS Assembler
# Don't need a elemCallBack since property info exists in bdf
FEAAssembler.initialize()

# Read in forces from BDF and create tacs struct problems
SPs = FEAAssembler.createTACSProbsFromBDF()

# Solve each structural problem and write solutions
funcs = {}
for problem in SPs.values():
    # Add eval functions to problem
    problem.addFunction("mass", functions.StructuralMass)
    problem.addFunction(
        "ks_disp", functions.KSDisplacement, ksWeight=100.0, direction=[1.0, 1.0, 1.0]
    )
    # Solve
    problem.solve()
    # Evaluate functions
    problem.evalFunctions(funcs)
    problem.writeSolution(outputDir=os.path.dirname(__file__))

if comm.rank == 0:
    pprint(funcs)
