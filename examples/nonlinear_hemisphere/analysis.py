"""
==============================================================================
Large deformation hemispherical shell
==============================================================================
@File    :   analysis.py
@Date    :   2023/01/25
@Author  :   Alasdair Christison Gray
@Description : This code runs a geometrically nonlinear analysis of a
hemispherical shell subject to radial point forces around its rim. The problem
is taken from section 3.4 of "Popular benchmark problems for geometric
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
BDF_FILE = os.path.join(os.path.dirname(__file__), "hemisphere.bdf")
E = 6.825e7  # Young's modulus
NU = 0.3  # Poisson's ratio
RHO = 1.0  # density
YIELD_STRESS = 1.0  # yield stress
THICKNESS = 0.04  # Shell thickness
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
    "continuationInitialStep": 1.0,
    "newtonSolverMaxIter": 50,
    "newtonSolverUseEW": True,
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
    "newtonSolverMaxLinIters": 10,
    "printTiming": True,
}
problem = FEAAssembler.createStaticProblem("RadialForces", options=probOptions)

# ==============================================================================
# Find tip force points
# ==============================================================================
bdfInfo = FEAAssembler.getBDFInfo()
# cross-reference bdf object to use some of pynastran's advanced features
bdfInfo.cross_reference()
nodeCoords = bdfInfo.get_xyz_in_coord()
nastranNodeNums = list(bdfInfo.node_ids)
loadPoints = np.array([[10.0, 0.0, 0.0], [0.0, 10.0, 0.0]])
loadPointNodeIDs = []
for ii in range(loadPoints.shape[0]):
    # find the closest node to the load point
    dists = np.linalg.norm(nodeCoords - loadPoints[ii, :], axis=1)
    closestNode = np.argmin(dists)
    loadPointNodeIDs.append(nastranNodeNums[closestNode])

# ==============================================================================
# Add tip loads
# ==============================================================================
PMax = 400.0
nodalForces = np.array(
    [[-PMax / 2, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, PMax / 2, 0.0, 0.0, 0.0, 0.0]]
)
problem.addLoadToNodes(loadPointNodeIDs, nodalForces, nastranOrdering=True)

# ==============================================================================
# Add functions
# ==============================================================================

# KS approximation of the maximum failure value
problem.addFunction("KSFailure", functions.KSFailure, ksWeight=80.0, ftype="discrete")

# Maximum displacement in the y and z-directions (KS with a very large weight to get a true max)
problem.addFunction(
    "MaxYDisp",
    functions.KSDisplacement,
    direction=np.array([0.0, 1.0, 0.0]),
    ksWeight=1e20,
    ftype="discrete",
)

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
    pprint(funcsSens)
