__all__ = ["blade_elemCallBack"]

"""
Adapted by Sean Engelstad
Source : Alasdair Christian Gray
"""

from tacs import elements, constitutive, TACS
import numpy as np

tacs_dtype = TACS.dtype

# define the element callback for TACS
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
skinPlyAngles = np.deg2rad(np.array([0.0, -45.0, 45.0, 90.0])).astype(tacs_dtype)
skinPlyFracs = np.array([44.41, 22.2, 22.2, 11.19], dtype=tacs_dtype) / 100.0
sparRibPlyAngles = np.deg2rad(np.array([0.0, -45.0, 45.0, 90.0])).astype(tacs_dtype)
sparRibPlyFracs = np.array([10.0, 35.0, 35.0, 20.0], dtype=tacs_dtype) / 100.0

# ==============================================================================
# Design variable values, bounds, and scaling factors
# ==============================================================================
# Panel length
panelLengthMax = np.inf
panelLengthMin = 0.0
panelLengthScale = 1.0

# Stiffener pitch
stiffenerPitch = tacs_dtype(0.2)  # m
stiffenerPitchMax = 0.5  # m
stiffenerPitchMin = 0.05  # m
stiffenerPitchScale = 1.0

# Panel thickness
panelThickness = tacs_dtype(0.02)  # m
panelThicknessMax = 0.1  # m
panelThicknessMin = 0.002  # m
panelThicknessScale = 100.0

# ply fraction bounds
plyFractionMax = 1.0
plyFractionMin = 0.1
plyFractionScale = 1.0

# Stiffener height
stiffenerHeight = tacs_dtype(0.05)  # m
stiffenerHeightMax = 0.1  # m
stiffenerHeightMin = 0.002  # m
stiffenerHeightScale = 10.0

# Stiffener thickness
stiffenerThickness = tacs_dtype(0.02)  # m
stiffenerThicknessMax = 0.1  # m
stiffenerThicknessMin = 0.002  # m
stiffenerThicknessScale = 100.0

# --- Stiffener axis directions ---
TESparDirection = np.array([0.34968083, 0.93686889, 0.0])
VerticalDirection = np.array([0.0, 0.0, 1.0])

# ==============================================================================
# Element callback function
# ==============================================================================


def blade_elemCallBack(
    dvNum, compID, compDescript, elemDescripts, specialDVs, **kwargs
):

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
    if "OML" in compDescript:
        plyAngles = skinPlyAngles
        panelPlyFractions = skinPlyFracs
        refAxis = TESparDirection
        panelLength = 0.65
    else:
        plyAngles = sparRibPlyAngles
        panelPlyFractions = sparRibPlyFracs
        refAxis = VerticalDirection
        if "rib" in compDescript:
            panelLength = 0.38
        elif "spar" in compDescript:
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

    stiffenerHeightNum = currDVNum
    DVScales.append(stiffenerHeightScale)
    currDVNum += 1

    stiffenerThicknessNum = currDVNum
    DVScales.append(stiffenerThicknessScale)
    currDVNum += 1

    # print(f"making blade stiffeners")
    con = constitutive.BladeStiffenedShellConstitutive(
        panelPly=ply,
        stiffenerPly=ply,
        panelLength=panelLength,
        stiffenerPitch=stiffenerPitch,
        panelThick=panelThickness,
        panelPlyAngles=plyAngles,
        panelPlyFracs=panelPlyFractions,
        stiffenerHeight=stiffenerHeight,
        stiffenerThick=stiffenerThickness,
        stiffenerPlyAngles=plyAngles,
        stiffenerPlyFracs=stiffenerPlyFractions,
        panelLengthNum=panelLengthNum,
        stiffenerPitchNum=stiffenerPitchNum,
        panelThickNum=panelThicknessNum,
        stiffenerHeightNum=stiffenerHeightNum,
        stiffenerThickNum=stiffenerThicknessNum,
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
