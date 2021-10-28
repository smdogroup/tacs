"""
RBE3 example where an RBE3 is used to connect two separate cantilevered plate sections.
A load is applied at the centroid of the RBE3. NOTE: This is the same case as the RBE2
example, with the RBE2 replaced with an RBE3.
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
    'printtimings':True,
    'writeSolution':True,
}

bdfFile = os.path.join(os.path.dirname(__file__), 'rbe3.bdf')
# Set up TACS Assembler
# Don't need a elemCallBack since property info exists in bdf
FEASolver = pyTACS(bdfFile, options=structOptions)
# Set up constitutive objects and elements
FEASolver.createTACSAssembler()

# Read in forces from BDF and create tacs static problems
SPs = FEASolver.createTACSProbsFromBDF()

# Solve state
for caseID in SPs:
    FEASolver(SPs[caseID])
    FEASolver.writeSolution(outputDir=os.path.dirname(__file__))


