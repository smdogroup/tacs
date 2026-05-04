"""
==============================================================================
Analysis of the Simple Transonic Wing with blade stiffened constitutive model
==============================================================================
@Author  :   Alasdair Christison Gray
@Description : This file performs a basic structural analysis of the Simple Transonic
wing's wingbox using the blade stiffened shell constitutive model and
a composite material. This case is based on the Simple Transonic Wing Benchmark setup used in the paper.
"A Proposed Benchmark Model for Practical Aeroelastic Optimization of Aircraft Wings"
(https://doi.org/10.2514/6.2024-2775)
"""

# ==============================================================================
# Standard Python modules
# ==============================================================================
import os
from pprint import pprint

# ==============================================================================
# External Python modules
# ==============================================================================
import numpy as np
import pandas as pd
from tacs import elements, constitutive, functions, pyTACS, TACS
from mpi4py import MPI

import argparse

parser = argparse.ArgumentParser(description="Run the STW structural analysis example.")
parser.add_argument(
    "--output",
    type=str,
    default=os.path.join(os.path.dirname(__file__), "output"),
    help="Directory to write output files to.",
)
parser.add_argument(
    "--useDVFile",
    action="store_true",
    default=False,
    help="Whether to use the design variable values from the StructDVs.csv file or to use the default values defined in the script.",
)
args = parser.parse_args()
os.makedirs(args.output, exist_ok=True)
np.set_printoptions(linewidth=200, precision=3, suppress=True)

comm = MPI.COMM_WORLD
dtype = TACS.dtype

# ==============================================================================
# Composite properties
# ==============================================================================
# The material properties are taken from "Aerostructural tradeoffs for tow-steered
# composite wings" by Tim Brooks, https://doi.org/10.2514/1.C035699
# The ply fractions are taken from chapter 7 of the PhD thesis of Johannes Dillinger,
# available at https://repository.tudelft.nl/islandora/object/uuid%3A20484651-fd5d-49f2-9c56-355bc680f2b7
compositeProperties = {
    "E11": 117.9e9,  # Young's modulus in 11 direction (Pa)
    "E22": 9.7e9,  # Young's modulus in 22 direction (Pa)
    "G12": 4.8e9,  # in-plane 1-2 shear modulus (Pa)
    "G13": 4.8e9,  # Transverse 1-3 shear modulus (Pa)
    "G23": 4.8e9,  # Transverse 2-3 shear modulus (Pa)
    "nu12": 0.35,  # 1-2 poisson's ratio
    "rho": 1.55e3,  # density kg/m^3
    "T1": 1648e6,  # Tensile strength in 1 direction (Pa)
    "C1": 1034e6,  # Compressive strength in 1 direction (Pa)
    "T2": 64e6,  # Tensile strength in 2 direction (Pa)
    "C2": 228e6,  # Compressive strength in 2 direction (Pa)
    "S12": 71e6,  # Shear strength direction (Pa)
}
skinPlyAngles = np.deg2rad(np.array([0.0, -45.0, 45.0, 90.0])).astype(dtype)
skinPlyFracs = np.array([44.41, 22.2, 22.2, 11.19], dtype=dtype) / 100.0
sparRibPlyAngles = np.deg2rad(np.array([0.0, -45.0, 45.0, 90.0])).astype(dtype)
sparRibPlyFracs = np.array([10.0, 35.0, 35.0, 20.0], dtype=dtype) / 100.0

kcorr = 5.0 / 6.0  # shear correction factor

# ==============================================================================
# Default design variable values, bounds, and scaling factors
# ==============================================================================
# The design variable values defined here are those specified for use in the benchmark STW structural and aeroelastic analyses.
flangeFraction = 1.0
# Panel length
defaultPanelLength = 0.5
defaultPanelLengthMax = np.inf
defaultPanelLengthMin = 0.0
defaultPanelLengthScale = 1.0

# Stiffener pitch
defaultStiffenerPitch = 0.15  # m
defaultStiffenerPitchMax = 0.5  # m
defaultStiffenerPitchMin = 0.15  # m
defaultStiffenerPitchScale = 1.0

# Panel thickness
defaultPanelThickness = 0.0065  # m
defaultPanelThicknessMax = 0.1  # m
defaultPanelThicknessMin = 0.6e-3  # m
defaultPanelThicknessScale = 100.0

# ply fraction bounds
defaultPlyFractionMax = 1.0
defaultPlyFractionMin = 0.1
defaultPlyFractionScale = 1.0

# Stiffener height
defaultStiffenerHeight = 0.05  # m
defaultStiffenerHeightMax = 0.15 / flangeFraction  # m
defaultStiffenerHeightMin = max(25e-3 * flangeFraction, 18e-3)  # m
defaultStiffenerHeightScale = 10.0

# Stiffener thickness
defaultStiffenerThickness = 0.006  # m
defaultStiffenerThicknessMax = 0.1  # m
defaultStiffenerThicknessMin = 0.6e-3  # m
defaultStiffenerThicknessScale = 100.0

# --- Stiffener axis directions ---
TESparDirection = np.array([0.34968083, 0.93686889, 0.0])
SpanwiseDirection = np.array([0.0, 1.0, 0.0])
VerticalDirection = np.array([0.0, 0.0, 1.0])

# These sizing values are taken from the Case 1 aeroelastic optimization results presented in
# "Aerostructural Optimization of the Simple Transonic Wing Using MPhys, ADflow, and TACS"
# (https://doi.org/10.2514/6.2025-2813)
structDVs = pd.read_csv(
    os.path.join(os.path.dirname(__file__), "StructDVs.csv"), index_col=0
)


# ==============================================================================
# Element callback function
# ==============================================================================
def element_callback(
    dvNum,
    compID,
    compDescript,
    elemDescripts,
    globalDVs,
    useDVsFromFile=False,
    **kwargs,
):
    # Get the design variable values for this component from the CSV file/dataframe
    dvValues = structDVs.loc[compDescript]

    prop = constitutive.MaterialProperties(
        rho=compositeProperties["rho"],
        E1=compositeProperties["E11"],
        E2=compositeProperties["E22"],
        G12=compositeProperties["G12"],
        G13=compositeProperties["G13"],
        G23=compositeProperties["G23"],
        nu12=compositeProperties["nu12"],
        T1=compositeProperties["T1"],
        C1=compositeProperties["C1"],
        T2=compositeProperties["T2"],
        C2=compositeProperties["C2"],
        S12=compositeProperties["S12"],
    )
    ply = constitutive.OrthotropicPly(1.0, prop)

    # Use a 0-deg biased layup for the skin and a +-45-deg biased layup spars and ribs.
    # Align the stiffeners in the skins with the trailing edge spar, and the stiffeners
    # in the spars and ribs vertically.
    if "SKIN" in compDescript:
        isInboard = any(
            [name in compDescript for name in ["SKIN.000", "SKIN.001", "SKIN.002"]]
        )
        if isInboard:
            refAxis = SpanwiseDirection
        else:
            refAxis = TESparDirection
        plyAngles = skinPlyAngles
        panelPlyFractions = skinPlyFracs
    else:
        plyAngles = sparRibPlyAngles
        panelPlyFractions = sparRibPlyFracs
        refAxis = VerticalDirection
    panelLength = dvValues.panelLength

    # Always use the 0-deg biased layup for the stiffeners
    stiffenerPlyFractions = skinPlyFracs
    numPlies = len(plyAngles)

    # --- Setup DV numbering and scaling ---

    # The ordering of the DVs used by the BladeStiffenedShell model is:
    # - panel length
    # - stiffener pitch (not used asa DV in this case)
    # - panel thickness
    # - panel ply fractions (not used asa DV in this case)
    # - stiffener height
    # - stiffener thickness
    # - stiffener ply fractions (not used asa DV in this case)
    currentDVNum = dvNum
    DVScales = []

    panelLength = dvValues.panelLength
    panelLengthNum = currentDVNum
    DVScales.append(defaultPanelLengthScale)
    currentDVNum += 1

    stiffenerPitch = (
        dvValues.stiffenerPitch if useDVsFromFile else defaultStiffenerPitch
    )

    panelThickness = (
        dvValues.panelThickness if useDVsFromFile else defaultPanelThickness
    )
    panelThicknessNum = currentDVNum
    DVScales.append(defaultPanelThicknessScale)
    currentDVNum += 1

    stiffenerHeight = (
        dvValues.stiffenerHeight if useDVsFromFile else defaultStiffenerHeight
    )
    stiffenerHeightNum = currentDVNum
    DVScales.append(defaultStiffenerHeightScale)
    currentDVNum += 1

    stiffenerThickness = (
        dvValues.stiffenerThickness if useDVsFromFile else defaultStiffenerThickness
    )
    stiffenerThicknessNum = currentDVNum
    DVScales.append(defaultStiffenerThicknessScale)
    currentDVNum += 1

    con = constitutive.BladeStiffenedShellConstitutive(
        panelPly=ply,
        stiffenerPly=ply,
        kcorr=kcorr,
        panelLength=dtype(panelLength),
        panelLengthNum=panelLengthNum,
        stiffenerPitch=dtype(stiffenerPitch),
        panelThick=dtype(panelThickness),
        panelThickNum=panelThicknessNum,
        panelPlyAngles=plyAngles.astype(dtype),
        panelPlyFracs=panelPlyFractions.astype(dtype),
        stiffenerHeight=dtype(stiffenerHeight),
        stiffenerHeightNum=stiffenerHeightNum,
        stiffenerThick=dtype(stiffenerThickness),
        stiffenerThickNum=stiffenerThicknessNum,
        stiffenerPlyAngles=plyAngles.astype(dtype),
        stiffenerPlyFracs=stiffenerPlyFractions.astype(dtype),
        flangeFraction=flangeFraction,
    )
    con.setPanelThicknessBounds(defaultPanelThicknessMin, defaultPanelThicknessMax)
    con.setStiffenerHeightBounds(defaultStiffenerHeightMin, defaultStiffenerHeightMax)
    con.setStiffenerThicknessBounds(
        defaultStiffenerThicknessMin, defaultStiffenerThicknessMax
    )

    # --- Create reference axis transform to define the stiffener direction ---
    transform = elements.ShellRefAxisTransform(refAxis)

    # --- Create the element object ---
    elem = elements.Quad4Shell(transform, con)

    return elem, DVScales


# ==============================================================================
# Setup constraints
# ==============================================================================
def setupConstraints(FEAAssembler):
    # From the benchmark problem specification:
    """Adjacency constraints are enforced to avoid abrupt changes in panel sizing. The change in panel and stiffener
    thicknesses between adjacent skin and spar panels is limited to 2.5 mm
    and the change in stiffener height to 10 mm. Some basic structural sizing rules suggested by Kassapoglou should be used
    on all panels:

    - The skin and stiffener thicknesses should be at least 0.6 mm
    - The stiffener heights should be at least 18 mm
    - The stiffener flange widths should be at least 25.4 mm (This and the above are satisfied by the stiffener height lower bound in this case.)
    - The aspect-ratio of the stiffener web (ℎ/t) should be between 5 and 30.
    - The thickness of the stiffener flanges on a panel should be no more than 15 times the panel thickness.
    - The stiffener flange width should be less than the stiffener pitch to avoid overlapping flanges. (This is covered by the stiffener height upper bound in this case)
    """
    thicknessAdjCon = 2.5e-3
    heightAdjCon = 10e-3
    stiffAspectMax = 5.0
    stiffAspectMin = 30.0

    constraints = []
    compIDs = {}
    for group in ["SPAR", "U_SKIN", "L_SKIN", "RIB"]:
        compIDs[group] = FEAAssembler.selectCompIDs(include=group.upper())

    localPanelLengthInd = 0
    localPanelThicknessInd = 1
    localstiffenerHeightInd = 2
    localStiffenerThicknessInd = 3

    # --- Adjacency constraints ---
    adjCon = FEAAssembler.createAdjacencyConstraint("AdjCon")
    for group in ["SPAR", "U_SKIN", "L_SKIN"]:
        adjCon.addConstraint(
            conName=group + "_panelThicknessAdj",
            compIDs=compIDs[group],
            lower=-thicknessAdjCon,
            upper=thicknessAdjCon,
            dvIndex=localPanelThicknessInd,
        )
        adjCon.addConstraint(
            conName=group + "_stiffenerThicknessAdj",
            compIDs=compIDs[group],
            lower=-thicknessAdjCon,
            upper=thicknessAdjCon,
            dvIndex=localStiffenerThicknessInd,
        )
        adjCon.addConstraint(
            conName=group + "_stiffenerHeightAdj",
            compIDs=compIDs[group],
            lower=-heightAdjCon,
            upper=heightAdjCon,
            dvIndex=localstiffenerHeightInd,
        )
    constraints.append(adjCon)

    # --- Arbitrary linear DV constraints ---
    dvCon = FEAAssembler.createDVConstraint("DVCon")

    # Flange thickness should be no more than 15x skin thickness
    # stiffenerThickness - 15 * panelThickness <= 0
    dvCon.addConstraint(
        conName="flangeThicknessMax",
        upper=0.0,
        dvIndices=[localPanelThicknessInd, localStiffenerThicknessInd],
        dvWeights=[-15.0, 1.0],
    )

    # Limit the aspect ratio of the stiffener
    # stiffenerHeight - stiffAspectMax * stiffenerThickness <= 0
    dvCon.addConstraint(
        conName="stiffenerAspectMax",
        upper=0.0,
        dvIndices=[localstiffenerHeightInd, localStiffenerThicknessInd],
        dvWeights=[1.0, -stiffAspectMax],
    )
    # stiffenerHeight - stiffAspectMin * stiffenerThickness >= 0
    dvCon.addConstraint(
        conName="stiffenerAspectMin",
        lower=0.0,
        dvIndices=[localstiffenerHeightInd, localStiffenerThicknessInd],
        dvWeights=[1.0, -stiffAspectMin],
    )
    constraints.append(dvCon)

    # --- Panel length constraints ---
    # These are what keep the panel length DVs consistent with the true panel lengths
    panelLengthCon = FEAAssembler.createPanelLengthConstraint("PanelLengthCon")
    panelLengthCon.addConstraint("PanelLength", dvIndex=localPanelLengthInd)
    constraints.append(panelLengthCon)
    return constraints

    # ==============================================================================
    # Setup static problem
    # ==============================================================================


def setupProblem(FEAAssembler):
    problem = FEAAssembler.createStaticProblem("StructAnalysis")

    # Add TACS Functions
    problem.addFunction("mass", functions.StructuralMass)

    # Add failure and mass functions for each skin and spar/rib group
    failureGroups = ["l_skin", "u_skin", "spar", "rib"]
    if failureGroups is not None:
        for group in failureGroups:
            compIDs = FEAAssembler.selectCompIDs(include=group.upper())
            problem.addFunction(
                f"{group}_ksFailure",
                functions.KSFailure,
                compIDs=compIDs,
                safetyFactor=1.5,
                ksWeight=100.0,
            )
            problem.addFunction(
                f"{group}_mass",
                functions.StructuralMass,
                compIDs=compIDs,
            )
    problem.addFunction("compliance", functions.Compliance)

    # Add gravity load
    LOAD_FACTOR = 2.5
    g = np.zeros(3)
    g[-1] = -9.81
    problem.addInertialLoad(LOAD_FACTOR * g)

    # Apply a uniform pressure over the lower skin, 30 kPa is roughly the right value to get a vertical
    # load equivalent to 2.5g
    lSkinCompIDs = FEAAssembler.selectCompIDs("L_SKIN")
    problem.addPressureToComponents(lSkinCompIDs, -LOAD_FACTOR * 30e3 / 2.5)
    return problem


# ==============================================================================
# Create pyTACS assembler
# ==============================================================================
structOptions = {
    "printtiming": True,
    "writeCoordinateFrame": True,
}

bdfFile = os.path.join(os.path.dirname(__file__), "wingbox-L2-Order2.bdf")
FEAAssembler = pyTACS(bdfFile, options=structOptions, comm=comm)


def element_callback_wrapper(
    dvNum, compID, compDescript, elemDescripts, globalDVs, **kwargs
):
    return element_callback(
        dvNum,
        compID,
        compDescript,
        elemDescripts,
        globalDVs,
        **kwargs,
        useDVsFromFile=args.useDVFile,
    )


FEAAssembler.initialize(element_callback_wrapper)
problem = setupProblem(FEAAssembler)
problem.setOptions(
    {
        "outputDir": args.output,
        "useMonitor": True,
        "monitorFrequency": 1,
    }
)
constraints = setupConstraints(FEAAssembler)

# Solve structural problem
problem.solve()

# Evaluate functions
funcs = {}
problem.evalFunctions(funcs)
for constraint in constraints:
    constraint.evalConstraints(funcs)


# ==============================================================================
# Extract tip displacement and twist
# ==============================================================================
components = ["SPAR.00", "SPAR.01", "RIB.22", "U_SKIN", "L_SKIN"]
nodes = {}
for comp in components:
    compIDs = FEAAssembler.selectCompIDs(include=comp)
    nodes[comp] = set(
        FEAAssembler.getGlobalNodeIDsForComps(compIDs, nastranOrdering=False)
    )

# The node at the front upper corner of the tip rib is the one node that is common to the upper skin, the front spar and the tip rib
frontUpperNodeGlobalID = list(
    nodes["U_SKIN"].intersection(nodes["RIB.22"]).intersection(nodes["SPAR.00"])
)[0]

# Similarly, the node at the rear upper corner of the tip rib is the one node that is common to the upper skin, the rear spar and the tip rib
rearUpperNodeGlobalID = list(
    nodes["U_SKIN"].intersection(nodes["RIB.22"]).intersection(nodes["SPAR.01"])
)[0]

frontUpperNodeLocalID = FEAAssembler.meshLoader.getLocalNodeIDsFromGlobal(
    frontUpperNodeGlobalID, nastranOrdering=False
)[0]
rearUpperNodeLocalID = FEAAssembler.meshLoader.getLocalNodeIDsFromGlobal(
    rearUpperNodeGlobalID, nastranOrdering=False
)[0]

# To compute the tip rotation we need the node coordinates
frontUpperCoord = FEAAssembler.meshLoader.getBDFNodes(
    frontUpperNodeGlobalID, nastranOrdering=False
)
rearUpperCoord = FEAAssembler.meshLoader.getBDFNodes(
    rearUpperNodeGlobalID, nastranOrdering=False
)

# Now retrieve the displacements at these nodes and compute the overall tip displacement and rotation
disp = problem.getVariables()
frontUpperDisp = None
rearUpperDisp = None
if frontUpperNodeLocalID != -1:
    frontUpperDisp = disp[6 * frontUpperNodeLocalID : 6 * frontUpperNodeLocalID + 3]
if rearUpperNodeLocalID != -1:
    rearUpperDisp = disp[6 * rearUpperNodeLocalID : 6 * rearUpperNodeLocalID + 3]

# broadcast front and rear upper displacements to all procs
hasFrontDisp = comm.allgather(frontUpperDisp is not None)
hasRearDisp = comm.allgather(rearUpperDisp is not None)
frontUpperDisp = comm.bcast(frontUpperDisp, root=np.argmax(hasFrontDisp))
rearUpperDisp = comm.bcast(rearUpperDisp, root=np.argmax(hasRearDisp))

# Compute the tip twist as the change in the angle of the line in the XZ plane between the front and rear upper nodes
x1 = frontUpperCoord[0]
z1 = frontUpperCoord[2]
dx1 = frontUpperDisp[0]
dz1 = frontUpperDisp[2]
x2 = rearUpperCoord[0]
z2 = rearUpperCoord[2]
dx2 = rearUpperDisp[0]
dz2 = rearUpperDisp[2]

tipZDisp = (dz1 + dz2) / 2

tipTwist = np.rad2deg(
    np.arctan2((z2 + dz2) - (z1 + dz1), (x2 + dx2) - (x1 + dx1))
    - np.arctan2(z2 - z1, x2 - x1)
)
funcs["tipZDisp"] = tipZDisp
funcs["tipTwist"] = tipTwist

if comm.rank == 0:
    pprint(funcs)

# ==============================================================================
# Solve adjoints and evaluate function sensitivities
# ==============================================================================
funcsSens = {}
problem.evalFunctionsSens(funcsSens)
for constraint in constraints:
    constraint.evalConstraintsSens(funcsSens)
if comm.rank == 0:
    pprint(funcsSens)

problem.writeSolution()
