"""
This problem performs a buckling analysis on a 4m x 3m flat plate.
The perimeter of the plate is pinned and loaded in compression (20kN/m) on its horizontal edges.
We use TACS buckling eigenvalue solver through the pyTACS BucklingProblem interface.
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
from tacs import pyTACS

comm = MPI.COMM_WORLD

# Instantiate FEAAssembler
bdfFile = os.path.join(os.path.dirname(__file__), "ss_plate_buckle.bdf")
FEAAssembler = pyTACS(bdfFile, comm=comm)
# Set up constitutive objects and elements
FEAAssembler.initialize()

# Setup static problem
bucklingProb = FEAAssembler.createBucklingProblem(name="buckle", sigma=10.0, numEigs=5)
bucklingProb.setOption("printLevel", 2)
# Add Loads
bucklingProb.addLoadFromBDF(loadID=1)

# solve and evaluate functions/sensitivities
funcs = {}
funcsSens = {}
bucklingProb.solve()
bucklingProb.evalFunctions(funcs)
bucklingProb.evalFunctionsSens(funcsSens)
bucklingProb.writeSolution(outputDir=os.path.dirname(__file__))


if comm.rank == 0:
    pprint(funcs)
    pprint(funcsSens)
