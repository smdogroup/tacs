"""
==============================================================================
Large deformation annular shell
==============================================================================
@File    :   analysis.py
@Date    :   2023/01/25
@Author  :   Alasdair Christison Gray
@Description : This code runs a geometrically nonlinear analysis of a
annular shell subject to vertical forces at one end. The problem
is taken from section 3.3 of "Popular benchmark problems for geometric
nonlinear analysis of shells" by Sze et al
(https://doi.org/10.1016/j.finel.2003.11.001).
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
BDF_FILE = os.path.join(os.path.dirname(__file__), "annulus.bdf")
E = 21e6  # Young's modulus
NU = 0.0  # Poisson's ratio
RHO = 1.0  # density
YIELD_STRESS = 1.0  # yield stress
THICKNESS = 0.03  # Shell thickness

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
    "isNonlinear": STRAIN_TYPE != "linear" or ROTATION_TYPE != "linear",
}
FEAAssembler = pyTACS(BDF_FILE, options=structOptions, comm=COMM)


def elemCallBack(dvNum, compID, compDescript, elemDescripts, specialDVs, **kwargs):
    matProps = constitutive.MaterialProperties(rho=RHO, E=E, nu=NU, YS=YIELD_STRESS)
    con = constitutive.IsoShellConstitutive(
        matProps, t=THICKNESS, tNum=dvNum, tlb=1e-2 * THICKNESS, tub=1e2 * THICKNESS
    )
    transform = None
    element = elementType(transform, con)
    tScale = [50.0]
    return element, tScale


FEAAssembler.initialize(elemCallBack)

probOptions = {
    "continuationInitialStep": 0.01,
    "continuationTargetIter": 6,
    "newtonSolverMaxIter": 50,
    "continuationRelTol": 1e-7,
    "newtonSolverMaxLinIters": 5,
    "continuationUsePredictor": True,
    "continuationNumPredictorStates": 8,
    # "newtonSolverUseEW": True,
    "nRestarts": 3,
    "subSpaceSize": 20,
    "nonlinearSolverMonitorVars": [
        "lambda",
        "linsolveriters",
        "linsolverres",
        "EWTol",
        "linesearchstep",
        "linesearchiters",
    ],
    "writeNLIterSolutions": True,
}
problem = FEAAssembler.createStaticProblem("Annulus", options=probOptions)


# ==============================================================================
# Find tip force points
# ==============================================================================
bdfInfo = FEAAssembler.getBDFInfo()
# cross-reference bdf object to use some of pynastran's advanced features
bdfInfo.cross_reference()
nodeCoords = bdfInfo.get_xyz_in_coord()
nastranNodeNums = list(bdfInfo.node_ids)
SPCNodeIDs = bdfInfo.get_SPCx_node_ids(1)

# The load points lie along the theta = 0 axis, however just searching for points along this axis will return
# both ends of the annulus as it covers a full 360 degrees. So we need to exclude any nodes that are subject to SPC's
nodeTheta = np.arctan2(nodeCoords[:, 1], nodeCoords[:, 0])
loadPointNodeIDs = np.argwhere(nodeTheta == 0.0).flatten() + 1
loadPointNodeIDs = [ii for ii in loadPointNodeIDs if ii not in SPCNodeIDs]

# ==============================================================================
# Add tip loads
# ==============================================================================
PMax = 0.8  # force per length
# Now we need to compute the nodal forces due to the line load, the nodes at the inner and outer radius will be subject to half the load of the other nodes
nodeRadius = np.linalg.norm(nodeCoords, axis=1)
innerRadius = np.min(nodeRadius)
outerRadius = np.max(nodeRadius)
totalForce = PMax * (outerRadius - innerRadius)
loadPointFactors = [
    1.0
    if nodeRadius[ii - 1] <= (innerRadius + 1e-4)
    or nodeRadius[ii - 1] >= (outerRadius - 1e-4)
    else 2.0
    for ii in loadPointNodeIDs
]
factorSum = np.sum(loadPointFactors)

nodalForces = np.zeros((len(loadPointNodeIDs), 6))
nodalForces[:, 2] = totalForce * np.array(loadPointFactors) / factorSum
problem.addLoadToNodes(loadPointNodeIDs, nodalForces, nastranOrdering=True)

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
problem.writeSolutionHistory(outputDir=os.path.dirname(__file__))

if COMM.rank == 0:
    pprint(funcs)
