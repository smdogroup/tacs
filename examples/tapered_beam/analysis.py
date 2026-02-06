"""
==============================================================================
Modal analysis of a tapered beam using TACS
==============================================================================
@File    :   analysis.py
@Date    :   2025/09/29
@Author  :   Alasdair Christison Gray
@Description :
"""

# ==============================================================================
# Standard Python modules
# ==============================================================================
import argparse

# ==============================================================================
# External Python modules
# ==============================================================================

# ==============================================================================
# Extension modules
# ==============================================================================
from tacs import pyTACS

parser = argparse.ArgumentParser(description="Tapered Beam Analysis")
parser.add_argument(
    "-i", "--input", type=str, default="tapered_beam_pbeam.bdf", help="Input BDF file"
)
args = parser.parse_args()
baseName = args.input.split(".bdf")[0]

FEASolver = pyTACS(args.input)
FEASolver.initialize()
problems = FEASolver.createTACSProbsFromBDF()
for problem in problems.values():
    problem.setOption("printLevel", 2)
    problem.setOption("writeSolution", True)
    problem.solve()
    problem.writeSolution(baseName=baseName)
