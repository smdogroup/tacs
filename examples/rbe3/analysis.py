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
# Extension modules
# ==============================================================================
from tacs import pyTACS

structOptions = {'writeSolution': True, }

bdfFile = os.path.join(os.path.dirname(__file__), 'rbe3.bdf')
# Load BDF file
FEAAssembler = pyTACS(bdfFile, options=structOptions)
# Set up TACS Assembler
# Don't need a elemCallBack since property info exists in bdf
FEAAssembler.initialize()

# Read in forces from BDF and create tacs struct problems
SPs = FEAAssembler.createTACSProbsFromBDF()

# Solve each structural problem and write solutions
for caseID in SPs:
    SPs[caseID].solve()
    SPs[caseID].writeSolution(outputDir=os.path.dirname(__file__))
