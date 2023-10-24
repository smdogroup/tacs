"""
==============================================================================
Box beam with skin buckling
==============================================================================
@File    :   analysis.py
@Date    :   2023/01/25
@Author  :   Alasdair Christison Gray
@Description : This code runs a nonlinear rectangular thin wall wingbox model
with a single rib made from plate elements. The model is clamped at one end,
with a smallout-of-plane force applied at the tip. Since the model only has one
rib in the middle, the ribs and spars buckle almost immediatly for a small
amount of load.
This case was originally created by Tim Brooks.
"""

# ==============================================================================
# Standard Python modules
# ==============================================================================
import os

# ==============================================================================
# External Python modules
# ==============================================================================
import numpy as np
from mpi4py import MPI
from pprint import pprint

# ==============================================================================
# Extension modules
# ==============================================================================
from tacs import pyTACS, constitutive, elements, functions


# ==============================================================================
# Constants
# ==============================================================================
COMM = MPI.COMM_WORLD
BDF_FILE = os.path.join(os.path.dirname(__file__), "Box.bdf")
E = 1.2e8  # Young's modulus
NU = 0.0  # Poisson's ratio
RHO = 1.0  # density
YIELD_STRESS = 1.0  # yield stress
THICKNESS = 0.02  # Shell thickness
STRAIN_TYPE = "nonlinear"
ROTATION_TYPE = "quadratic"

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
    tScale = [50.0]
    return element, tScale


FEAAssembler.initialize(elemCallBack)

probOptions = {
    "nRestarts": 3,
    "subSpaceSize": 20,
    "printLevel": 1,
}
newtonOptions = {
    "MaxLinIters": 10,
    "UseEW": True,
}
continuationOptions = {
    "InitialStep": 0.05,
    "TargetIter": 6,
    "RelTol": 1e-7,
    "UsePredictor": True,
    "NumPredictorStates": 4,
    "MaxIter": 60,
}

problem = FEAAssembler.createStaticProblem("TipForce", options=probOptions)
problem.nonlinearSolver.setOptions(continuationOptions)
problem.nonlinearSolver.innerSolver.setOptions(newtonOptions)

# ==============================================================================
# Find tip force points
# ==============================================================================
bdfInfo = FEAAssembler.getBDFInfo()
# cross-reference bdf object to use some of pynastran's advanced features
bdfInfo.cross_reference()
nodeCoords = bdfInfo.get_xyz_in_coord()
nastranNodeNums = list(bdfInfo.node_ids)
loadPoints = np.array(
    [[3.0, 40.0, 0.75], [3.0, 40.0, -0.75], [-3.0, 40.0, 0.75], [-3.0, 40.0, -0.75]]
)
loadPointNodeIDs = []
for ii in range(loadPoints.shape[0]):
    # find the closest node to the load point
    dists = np.linalg.norm(nodeCoords - loadPoints[ii, :], axis=1)
    closestNode = np.argmin(dists)
    loadPointNodeIDs.append(nastranNodeNums[closestNode])

# ==============================================================================
# Add tip loads
# ==============================================================================
tipForceTotal = 250.0
problem.addLoadToNodes(
    loadPointNodeIDs, [0.0, 0.0, tipForceTotal / 4, 0.0, 0.0, 0.0], nastranOrdering=True
)

# ==============================================================================
# Add functions
# ==============================================================================

# KS approximation of the maximum failure value
problem.addFunction("KSFailure", functions.KSFailure, ksWeight=80.0, ftype="discrete")

# Maximum displacement in the z-direction (KS with a very large weight to get a true max)
problem.addFunction(
    "MaxZDisp",
    functions.KSDisplacement,
    direction=np.array([0.0, 0.0, 1.0]),
    ksWeight=1e20,
    ftype="discrete",
)

# Compliance
problem.addFunction("Compliance", functions.Compliance)

# ==============================================================================
# Solve all problems and evaluate functions
# ==============================================================================
funcs = {}
funcsSens = {}
problem.solve()
problem.evalFunctions(funcs)
problem.evalFunctionsSens(funcsSens)
problem.writeSolution(outputDir=os.path.dirname(__file__))

if COMM.rank == 0:
    pprint(funcs)
