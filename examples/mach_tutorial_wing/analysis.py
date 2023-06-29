"""
==============================================================================
Analysis of the MACH tutorial wing with blade stiffened constitutive model
==============================================================================
@Date    :   2023/06/27
@Author  :   Alasdair Christison Gray
@Description : This file performs a basic structural analysis of the MACH
tutorial wing's wingbox using the blade stiffened shell constitutive model and
a composite material. This case is based on the setup used in the paper:
"Geometrically Nonlinear High-fidelity Aerostructural Optimization Including
Geometric Design Variables" (https://doi.org/10.2514/6.2023-3316)
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
from tacs import elements, constitutive, functions, pyTACS, TACS
from mpi4py import MPI

comm = MPI.COMM_WORLD

dtype = TACS.dtype

# Instantiate FEAAssembler
structOptions = {
    "printtiming": True,
}

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

kcorr = dtype(5.0 / 6.0)  # shear correction factor

# ==============================================================================
# Design variable values, bounds, and scaling factors
# ==============================================================================
# Panel length
panelLengthMax = np.inf
panelLengthMin = 0.0
panelLengthScale = 1.0

# Stiffener pitch
stiffenerPitch = dtype(0.2)  # m
stiffenerPitchMax = 0.5  # m
stiffenerPitchMin = 0.05  # m
stiffenerPitchScale = 1.0

# Panel thickness
panelThickness = dtype(0.02)  # m
panelThicknessMax = 0.1  # m
panelThicknessMin = 0.002  # m
panelThicknessScale = 100.0

# ply fraction bounds
plyFractionMax = 1.0
plyFractionMin = 0.1
plyFractionScale = 1.0

# Stiffener height
stiffenerHeight = dtype(0.05)  # m
stiffenerHeightMax = 0.1  # m
stiffenerHeightMin = 0.002  # m
stiffenerHeightScale = 10.0

# Stiffener thickness
stiffenerThickness = dtype(0.02)  # m
stiffenerThicknessMax = 0.1  # m
stiffenerThicknessMin = 0.002  # m
stiffenerThicknessScale = 100.0

# --- Stiffener axis directions ---
TESparDirection = np.array([0.34968083, 0.93686889, 0.0])
VerticalDirection = np.array([0.0, 0.0, 1.0])

# ==============================================================================
# Element callback function
# ==============================================================================


def elemCallBack(dvNum, compID, compDescript, elemDescripts, specialDVs, **kwargs):
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
    # The panel length values I set here are approximate, to get the real values, you'd
    # need to run an optimization with panel length design variables and constraints.
    if "SKIN" in compDescript:
        plyAngles = skinPlyAngles
        panelPlyFractions = skinPlyFracs
        refAxis = TESparDirection
        panelLength = 0.65
    else:
        plyAngles = sparRibPlyAngles
        panelPlyFractions = sparRibPlyFracs
        refAxis = VerticalDirection
        if "RIB" in compDescript:
            panelLength = 0.38
        elif "SPAR" in compDescript:
            panelLength = 0.36

    # Always use the 0-deg biased layup for the stiffeners
    stiffenerPlyFractions = skinPlyFracs
    numPlies = len(plyAngles)

    # --- Setup DV numbering and scaling ---

    # The ordering of the DVs used by the BladeStiffenedShell model is:
    # - panel length
    # - stiffener pitch
    # - panel thickness
    # - panel ply fractions (not used in this case)
    # - stiffener height
    # - stiffener thickness
    # - stiffener ply fractions (not used in this case)
    currDVNum = dvNum
    DVScales = []

    panelLengthNum = currDVNum
    DVScales.append(panelLengthScale)
    currDVNum += 1

    stiffenerPitchNum = currDVNum
    DVScales.append(stiffenerPitchScale)
    currDVNum += 1

    panelThicknessNum = currDVNum
    DVScales.append(panelThicknessScale)
    currDVNum += 1

    panelPlyFracNums = -np.ones(numPlies).astype(np.intc)

    stiffenerHeightNum = currDVNum
    DVScales.append(stiffenerHeightScale)
    currDVNum += 1

    stiffenerThicknessNum = currDVNum
    DVScales.append(stiffenerThicknessScale)
    currDVNum += 1

    stiffenerPlyFracNums = -np.ones(numPlies).astype(np.intc)

    con = constitutive.BladeStiffenedShellConstitutive(
        panelPly=ply,
        stiffenerPly=ply,
        kcorr=kcorr,
        panelLength=panelLength,
        panelLengthNum=panelLengthNum,
        stiffenerPitch=stiffenerPitch,
        stiffenerPitchNum=stiffenerPitchNum,
        panelThick=panelThickness,
        panelThickNum=panelThicknessNum,
        numPanelPlies=numPlies,
        panelPlyAngles=plyAngles,
        panelPlyFracs=panelPlyFractions,
        panelPlyFracNums=panelPlyFracNums,
        stiffenerHeight=stiffenerHeight,
        stiffenerHeightNum=stiffenerHeightNum,
        stiffenerThick=stiffenerThickness,
        stiffenerThickNum=stiffenerThicknessNum,
        numStiffenerPlies=numPlies,
        stiffenerPlyAngles=plyAngles,
        stiffenerPlyFracs=stiffenerPlyFractions,
        stiffenerPlyFracNums=stiffenerPlyFracNums,
    )
    con.setStiffenerPitchBounds(stiffenerPitchMin, stiffenerPitchMax)
    con.setPanelThicknessBounds(panelThicknessMin, panelThicknessMax)
    con.setStiffenerHeightBounds(stiffenerHeightMin, stiffenerHeightMax)
    con.setStiffenerThicknessBounds(stiffenerThicknessMin, stiffenerThicknessMax)

    # --- Create reference axis transform to define the stiffener direction ---
    transform = elements.ShellRefAxisTransform(refAxis)

    # --- Create the element object ---
    if elemDescripts[0] == "CQUAD4":
        elem = elements.Quad4Shell(transform, con)
    elif elemDescripts[0] == "CQUAD9":
        elem = elements.Quad9Shell(transform, con)
    elif elemDescripts[0] == "CQUAD16":
        elem = elements.Quad16Shell(transform, con)

    return elem, DVScales


# ==============================================================================
# Create pyTACS assembler
# ==============================================================================
structOptions = {
    "printtiming": True,
}

bdfFile = os.path.join(os.path.dirname(__file__), "wingbox-L2-Order2.bdf")
FEAAssembler = pyTACS(bdfFile, options=structOptions, comm=comm)

# Set up elements and TACS assembler
FEAAssembler.initialize(elemCallBack)

# ==============================================================================
# Setup static problem
# ==============================================================================
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
g = np.zeros(3)
g[-1] = -9.81
problem.addInertialLoad(2.5 * g)

# Add simple distributed vertical load over the skins
totalForce = 0.5 * 2.5 * 55e3 * 9.81
F = np.zeros(6)
F[2] = totalForce
skinCompIDs = FEAAssembler.selectCompIDs("SKIN")
problem.addLoadToComponents(skinCompIDs, F, averageLoad=True)

# ==============================================================================
# Setup constraints
# ==============================================================================

thicknessAdjCon = 2.5e-3  # 2.5mm
heightAdjCon = 1e-2  # 1cm
pitchAdjCon = 5e-2  # 5cm
thickDiffMax = (
    2.5e-3  # 2.5mm, Max allowable thickness difference between skin and stiffener
)
stiffAspectMax = 10.0  # Maximum allowable stiffener aspect ratio (height/thickness)
stiffAspectMin = 2.0  # Minimum allowable stiffener aspect ratio (height/thickness)

compIDs = {}
for group in ["SPAR", "U_SKIN", "L_SKIN", "RIB"]:
    compIDs[group] = FEAAssembler.selectCompIDs(include=group.upper())

constraints = []

# Add adjacency constraints on panel thickness, stiffener thickness and stiffener height to everything but the ribs
adjCon = FEAAssembler.createAdjacencyConstraint("AdjCon")
for group in ["SPAR", "U_SKIN", "L_SKIN"]:
    adjCon.addConstraint(
        conName=group + "_panelThicknessAdj",
        compIDs=compIDs[group],
        lower=-thicknessAdjCon,
        upper=thicknessAdjCon,
        dvIndex=2,
    )
    adjCon.addConstraint(
        conName=group + "_stiffenerThicknessAdj",
        compIDs=compIDs[group],
        lower=-thicknessAdjCon,
        upper=thicknessAdjCon,
        dvIndex=4,
    )
    adjCon.addConstraint(
        conName=group + "_stiffenerHeightAdj",
        compIDs=compIDs[group],
        lower=-heightAdjCon,
        upper=heightAdjCon,
        dvIndex=3,
    )
    adjCon.addConstraint(
        conName=group + "_stiffenerPitchAdj",
        compIDs=compIDs[group],
        lower=-pitchAdjCon,
        upper=pitchAdjCon,
        dvIndex=1,
    )
constraints.append(adjCon)

# Add constraints between the DV's on each panel
dvCon = FEAAssembler.createDVConstraint("DVCon")
# Limit the difference in thickness between the panel and stiffener
# -thickDiffMax <= (panelThickness - stiffenerThickness) <= thickDiffMax
dvCon.addConstraint(
    conName="thickDiffLimit",
    lower=-thickDiffMax,
    upper=thickDiffMax,
    dvIndices=[2, 4],
    dvWeights=[1.0, -1.0],
)
# Limit the aspect ratio of the stiffener
# stiffenerHeight - (stiffAspectMax * stiffenerThickness) <= 0
dvCon.addConstraint(
    conName="stiffenerAspectMax",
    upper=0.0,
    dvIndices=[3, 4],
    dvWeights=[1.0, -stiffAspectMax],
)
# stiffenerHeight - (stiffAspectMin * stiffenerThickness) >= 0
dvCon.addConstraint(
    conName="stiffenerAspectMin",
    lower=0.0,
    dvIndices=[3, 4],
    dvWeights=[1.0, -stiffAspectMin],
)
# Ensure there is space between the stiffeners
# 2*flangeFraction - stiffenerPitch <= 0
dvCon.addConstraint(
    conName="stiffSpacingMin",
    upper=0.0,
    dvIndices=[3, 2],
    dvWeights=[2.0, -1.0],
)
constraints.append(dvCon)

panelLengthCon = FEAAssembler.createPanelLengthConstraint("PanelLengthCon")
panelLengthCon.addConstraint("PanelLength", dvIndex=0)
constraints.append(panelLengthCon)

# ==============================================================================
# Solve static problem
# ==============================================================================
# Solve structural problem
problem.solve()

# Evaluate functions
funcs = {}
problem.evalFunctions(funcs)
for constraint in constraints:
    constraint.evalConstraints(funcs)

if comm.rank == 0:
    pprint(funcs)

# Solve adjoints and evaluate function sensitivities
funcsSens = {}
problem.evalFunctionsSens(funcsSens)
for constraint in constraints:
    constraint.evalConstraintsSens(funcsSens)
# if comm.rank == 0:
#     pprint(funcsSens)

# Write out solution
problem.writeSolution(outputDir=os.path.dirname(__file__))
