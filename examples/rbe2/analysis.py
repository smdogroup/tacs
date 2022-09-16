"""
RBE2 example where an RBE2 is used to connect two separate cantilevered plate sections.
A load is applied at the centroid of the RBE2.
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
import numpy

# ==============================================================================
# Extension modules
# ==============================================================================
from tacs import pyTACS

comm = MPI.COMM_WORLD

# Instantiate FEASolver
structOptions = {
    "printtiming": True,
}

bdfFile = os.path.join(os.path.dirname(__file__), "rbe_test.bdf")
FEAAssembler = pyTACS(bdfFile, comm, options=structOptions)
# Set up TACS Assembler
# Don't need a elemCallBack since property info exists in bdf
FEAAssembler.initialize()

# Read in forces from BDF and create tacs static problems
SPs = FEAAssembler.createTACSProbsFromBDF()

# Solve problems
for problem in SPs.values():
    problem.solve()
    problem.writeSolution(outputDir=os.path.dirname(__file__))
