"""
==============================================================================
Nonlinear cantilever beam analysis
==============================================================================
@File    :   analysis.py
@Date    :   2023/01/24
@Author  :   Alasdair Christison Gray
@Description : This code runs an analysis of a cantilever beam modeled with
shell elements subject to a vertical tip force. The problem is taken from
section 3.1 of "Popular benchmark problems for geometric nonlinear analysis of
shells" by Sze et al (https://doi.org/10.1016/j.finel.2003.11.001).
"""

# ==============================================================================
# Standard Python modules
# ==============================================================================
import os
import pickle
import argparse

# ==============================================================================
# External Python modules
# ==============================================================================
import numpy as np
from mpi4py import MPI

# ==============================================================================
# Extension modules
# ==============================================================================
from tacs import pyTACS, constitutive, elements

# ==============================================================================
# Parse command line arguments
# ==============================================================================
parser = argparse.ArgumentParser()
parser.add_argument(
    "--strainType", type=str, default="nonlinear", choices=["linear", "nonlinear"]
)
parser.add_argument(
    "--rotationType",
    type=str,
    default="quadratic",
    choices=["linear", "quadratic", "quaternion"],
)
args = parser.parse_args()

# ==============================================================================
# Constants
# ==============================================================================
COMM = MPI.COMM_WORLD
PWD = os.path.dirname(__file__)
BDF_FILE = os.path.join(PWD, "Beam.bdf")
E = 1.2e6  # Young's modulus
NU = 0.0  # Poisson's ratio
RHO = 1.0  # density
YIELD_STRESS = 1.0  # yield stress
THICKNESS = 0.1  # Shell thickness
FORCE_MULTIPLIER = 4.0  # Multiplier applied to the baseline force of EI/L^2
STRAIN_TYPE = args.strainType
ROTATION_TYPE = args.rotationType

elementType = None
if STRAIN_TYPE == "linear":
    if ROTATION_TYPE == "linear":
        elementType = elements.Quad4Shell
    elif ROTATION_TYPE == "quadratic":
        elementType = elements.Quad4ShellModRot
    elif ROTATION_TYPE == "quaternion":
        elementType = elements.Quad4ShellQuaternion
elif STRAIN_TYPE == "nonlinear":
    if ROTATION_TYPE == "linear":
        elementType = elements.Quad4NonlinearShell
    elif ROTATION_TYPE == "quadratic":
        elementType = elements.Quad4NonlinearShellModRot
    elif ROTATION_TYPE == "quaternion":
        elementType = elements.Quad4NonlinearShellQuaternion

if elementType is None:
    raise RuntimeError("Invalid element type, check STRAIN_TYPE and ROTATION_TYPE.")

# ==============================================================================
# Create pyTACS Assembler and problems
# ==============================================================================
structOptions = {
    "printtiming": True,
}
FEAAssembler = pyTACS(BDF_FILE, options=structOptions, comm=COMM)


def elemCallBack(dvNum, compID, compDescript, elemDescripts, specialDVs, **kwargs):
    matProps = constitutive.MaterialProperties(rho=RHO, E=E, nu=NU, ys=YIELD_STRESS)
    con = constitutive.IsoShellConstitutive(
        matProps, t=THICKNESS, tNum=dvNum, tlb=1e-2 * THICKNESS, tub=1e2 * THICKNESS
    )
    transform = None
    element = elementType(transform, con)
    tScale = [10.0]
    return element, tScale


FEAAssembler.initialize(elemCallBack)

probOptions = {
    "printTiming": True,
    "printLevel": 1,
}
newtonOptions = {"useEW": True, "MaxLinIters": 10}
continuationOptions = {
    "CoarseRelTol": 1e-3,
    "InitialStep": 0.2,
    "UsePredictor": True,
    "NumPredictorStates": 7,
}
forceProblem = FEAAssembler.createStaticProblem("TipForce", options=probOptions)
try:
    forceProblem.nonlinearSolver.innerSolver.setOptions(newtonOptions)
    forceProblem.nonlinearSolver.setOptions(continuationOptions)
except AttributeError:
    pass

# ==============================================================================
# Determine beam dimensions and other properties based on node coordinates
# ==============================================================================
bdfInfo = FEAAssembler.getBDFInfo()
# cross-reference bdf object to use some of pynastrans advanced features
bdfInfo.cross_reference()
nodeCoords = bdfInfo.get_xyz_in_coord()
beamLength = np.max(nodeCoords[:, 0]) - np.min(nodeCoords[:, 0])
beamWidth = np.max(nodeCoords[:, 1]) - np.min(nodeCoords[:, 1])
Izz = beamWidth * THICKNESS**3 / 12.0

# ==============================================================================
# Add tip loads for each case
# ==============================================================================
tipForce = FORCE_MULTIPLIER * E * Izz / beamLength**2

# In order to work for different mesh sizes, we need to find the tip node IDs
# ourselves, we do this by finding the indices of the nodes whose x coordinate
# is within a tolerance of the max X coordinate in the mesh
tipNodeInds = np.nonzero(np.abs(np.max(nodeCoords[:, 0]) - nodeCoords[:, 0]) <= 1e-6)[0]
nastranNodeNums = list(bdfInfo.node_ids)
tipNodeIDs = [nastranNodeNums[ii] for ii in tipNodeInds]
numTipNodes = len(tipNodeIDs)

forceProblem.addLoadToNodes(
    tipNodeIDs, [0, 0, tipForce / numTipNodes, 0, 0, 0], nastranOrdering=True
)

# ==============================================================================
# Run analysis with load scales in 5% increments from 5% to 100%
# ==============================================================================
forceFactor = np.arange(0.05, 1.01, 0.05)
ForceVec = np.copy(forceProblem.F_array)

results = {"zDisp": [0.0], "xDisp": [0.0], "yRot": [0.0], "tipForce": [0.0]}

for scale in forceFactor:
    Fext = (scale - 1.0) * ForceVec
    forceProblem.solve(Fext=Fext)

    fileName = f"{STRAIN_TYPE}_{ROTATION_TYPE}"
    forceProblem.writeSolution(outputDir=PWD, baseName=fileName)
    disps = forceProblem.u_array
    xDisps = disps[0::6]
    zDisps = disps[2::6]
    yRot = disps[4::6]
    results["tipForce"].append(scale)
    results["xDisp"].append(xDisps[-1])
    results["zDisp"].append(zDisps[-1])
    results["yRot"].append(yRot[-1])

for key in results:
    results[key] = np.array(results[key])

with open(os.path.join(PWD, f"TACS-Disps-{fileName}.pkl"), "wb") as f:
    pickle.dump(results, f)
