"""
==============================================================================
Nonlinear cantilever beam analysis
==============================================================================
@File    :   analysis.py
@Date    :   2023/01/24
@Author  :   Alasdair Christison Gray
@Description :
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
BDF_FILE = os.path.join(os.path.dirname(__file__), "Beam.bdf")
E = 1.2e6  # Young's modulus
NU = 0.0  # Poisson's ratio
RHO = 1.0  # density
YIELD_STRESS = 1.0  # yield stress
THICKNESS = 0.1  # Shell thickness
FORCE_MULTIPLIER = 1.0  # Multiplier applied to the baseline force of EI/L^2
MOMENT_MULTIPLIER = 0.1  # Multiplier applied to the baseline moment of 2pi * EI/L (which results in a full rotation)
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
    "isNonlinear": True,
}
FEAAssembler = pyTACS(BDF_FILE, options=structOptions, comm=COMM)


def elemCallBack(dvNum, compID, compDescript, elemDescripts, specialDVs, **kwargs):
    matProps = constitutive.MaterialProperties(rho=RHO, E=E, nu=NU, YS=YIELD_STRESS)
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
    "skipFirstNLineSearch": 0,
    "newtonSolverCoarseRelTol": 1e-3,
    "continuationInitialStep": 0.05,
    "continuationUsePredictor": True,
    "continuationNumPredictorStates": 7,
    "writeNLIterSolutions": True,
}
forceProblem = FEAAssembler.createStaticProblem("TipForce", options=probOptions)
momentProblem = FEAAssembler.createStaticProblem("TipMoment", options=probOptions)
problems = [forceProblem, momentProblem]


# ==============================================================================
# Determine beam dimensions and other properties
# ==============================================================================
bdfInfo = FEAAssembler.getBDFInfo()
# cross-reference bdf object to use some of pynastrans advanced features
bdfInfo.cross_reference()
nodeCoords = bdfInfo.get_xyz_in_coord()
beamLength = np.max(nodeCoords[:, 0]) - np.min(nodeCoords[:, 0])
beamWidth = np.max(nodeCoords[:, 1]) - np.min(nodeCoords[:, 1])
I = beamWidth * THICKNESS**3 / 12.0

# ==============================================================================
# Add tip loads for each case
# ==============================================================================
tipForce = FORCE_MULTIPLIER * E * I / beamLength**2
tipMoment = MOMENT_MULTIPLIER * -2 * np.pi * E * I / beamLength

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
momentProblem.addLoadToNodes(
    tipNodeIDs, [0, 0, 0, 0, tipMoment / numTipNodes, 0], nastranOrdering=True
)

# ==============================================================================
# Add functions for each problem
# ==============================================================================

for problem in problems:
    # KS approximation of the maximum failure value
    problem.addFunction(
        "KSFailure", functions.KSFailure, ksWeight=80.0, ftype="discrete"
    )

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
for problem in problems:
    problem.solve()
    problem.evalFunctions(funcs)
    problem.evalFunctionsSens(funcsSens)
    problem.writeSolution(outputDir=os.path.dirname(__file__))

if COMM.rank == 0:
    pprint(funcs)
