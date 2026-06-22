"""
==============================================================================
Nastran beam compatibility regression suite — BDF generator
==============================================================================
@File    :   generate_inputs.py
@Author  :   Alasdair Christison Gray
@Description : Emits 72 self-contained Nastran BDFs (36 configurations, each
with a SOL 101 static driver and a SOL 103 modal driver) for the beam
compatibility regression suite.

Each configuration is a 1 m aluminium cantilever discretised into 30 beam
elements. The configurations span all 8 supported element/property/section
combinations (CBAR/PBAR, CBAR/PBARL+{BAR,TUBE,T}, CBEAM/PBEAM,
CBEAM/PBEAML+{BAR,TUBE,T}), each extended along a compounding ladder of
feature "rungs" (taper, section rotation, shear-centre offset, neutral-axis
offset, non-structural mass, NSM offset). See ``DESIGN.md`` in this directory
for the authoritative spec.

The BDFs are self-contained (no INCLUDE) and are written end-to-end through
pyNastran's BDF writer. ``NASTRAN BARMASS=1`` is set via
``model.system_command_lines`` so CBAR torsional rotary inertia is included in
the mass matrix (matching TACS Beam2). The SOL 101 driver carries 6 SUBCASEs
(one per tip-load direction); the SOL 103 driver carries a single modal
SUBCASE plus an EIGRL card requesting the first 10 modes. The case matrix
(which element/property/feature configurations are generated) is defined once
in ``cases.py``; this script imports it.

Usage:
    python generate_inputs.py
"""

# ==============================================================================
# Standard Python modules
# ==============================================================================
from pathlib import Path

# ==============================================================================
# External Python modules
# ==============================================================================
import numpy as np
from pyNastran.bdf.bdf import BDF
from pyNastran.bdf.case_control_deck import CaseControlDeck
from cases import iterCases, STATIC_LOAD_NAMES

# ==============================================================================
# Material and geometry constants
# ==============================================================================
E = 71.7e9  # Young's modulus, Pa (Aluminium 7075)
NU = 0.33  # Poisson's ratio
G = E / (2.0 * (1.0 + NU))  # Shear modulus, Pa
RHO = 2810.0  # Density, kg/m^3

BEAM_LENGTH = 1.0  # m
NUM_ELEMENTS = 30

# Rectangular section (PBAR / PBARL+BAR / PBEAM / PBEAML+BAR).
# WIDTH is the local-z dimension; DEPTH (the through-thickness dimension) is the
# local-y dimension. This matches the TACS convention that thickness lies along
# local-y (see the figure in src/constitutive/TACSIsoRectangleBeamConstitutive.h)
# and the Nastran PBARL/PBEAML BAR ordering dims=[local-y, local-z] used below.
WIDTH1, DEPTH1 = 0.20, 0.05  # end-A (root) dimensions
WIDTH2, DEPTH2 = 0.10, 0.025  # end-B (tip) dimensions = half of end-A

# Tube sections are deliberately THIN-WALLED (wall = 10% of the outer radius).
# TACS's TACSIsoTubeBeamConstitutive hard-codes the *thin-walled* transverse
# shear-correction factor kappa = 2(1+nu)/(4+3nu), which only matches Nastran
# (and reality) in the thin-wall limit. A thick or near-solid tube would
# disagree with Nastran by ~1-2% in shear-affected responses. See the note in
# src/constitutive/TACSIsoTubeBeamConstitutive.cpp about generalising to the
# full Cowper hollow-circular coefficient, which would let us test thick tubes
# and solid rods against Nastran too.
#
# Tube section (PBARL+TUBE / PBEAML+TUBE): Nastran TUBE dims are
# [outer radius, inner radius].
TUBE_OUTER_R, TUBE_INNER_R = 0.10, 0.09  # end-A: outer radius, inner radius

# Tube2 section (PBARL+TUBE2 / PBEAML+TUBE2): Nastran TUBE2 dims are
# [outer radius, wall thickness]. Describes the same physical tube as TUBE above
# (wall = outer radius - inner radius = 0.10 - 0.09 = 0.01).
TUBE2_OUTER_R, TUBE2_WALL = 0.10, 0.01  # end-A: outer radius, wall thickness

# Rod section (PBARL+ROD / PBEAML+ROD): a solid circle, single dim [radius].
# NOTE: a ROD is solid (inner radius 0), so TACS's thin-wall shear factor is
# wrong for it (the solid-circle value is ~0.89, not ~0.53). ROD results will
# disagree with Nastran in shear-affected responses until the Cowper factor is
# implemented; see the .cpp note referenced above.
ROD_R1 = 0.05  # end-A radius

# T-section (PBARL+T / PBEAML+T): [flange width, total depth, web t, flange t].
T_DIMS1 = [0.20, 0.10, 0.02, 0.02]  # end-A
T_DIMS2 = [0.10, 0.05, 0.01, 0.01]  # end-B = half of end-A

# Shear correction factor for the rectangular explicit-property cards.
SHEAR_CORRECTION = 5.0 / 6.0

# ==============================================================================
# Feature rung values
# ==============================================================================
# Section orientation vector rotated 30 degrees about the beam axis. The
# default orientation is [0, 1, 0]; the rotation tilts it towards local-z.
_ROTATION_RAD = np.deg2rad(30.0)
_ROTATION_VEC = np.array([0.0, np.cos(_ROTATION_RAD), np.sin(_ROTATION_RAD)])
_ROTATION_VEC = _ROTATION_VEC / np.linalg.norm(_ROTATION_VEC)
_DEFAULT_ORIENTATION_VEC = np.array([0.0, 1.0, 0.0])

# Shear-centre (grid) offset applied via the CBAR/CBEAM WA/WB fields.
_WA_OFFSET = np.array([0.0, 2.0 * WIDTH1, 2.0 * DEPTH1])

# PBEAM neutral-axis offset (N1/N2), same magnitudes as the shear-centre offset.
_N1, _N2 = 2.0 * WIDTH1, 2.0 * DEPTH1

# Non-structural mass per unit length (comparable to the section's own mass).
_NSM = RHO * WIDTH1 * DEPTH1

# PBEAM NSM offset (M1/M2).
_M1, _M2 = WIDTH1, DEPTH1

# ==============================================================================
# Loads
# ==============================================================================
FORCE_MAG = 10000.0  # N
MOMENT_MAG = 1000.0  # N*m

# (subcase name, load type, direction vector, load set id). Load ids are
# unique so a single model can hold all 6 static load cases simultaneously.
_STATIC_LOAD_DEFINITIONS = [
    ("force", [1.0, 0.0, 0.0], 10),
    ("force", [0.0, 1.0, 0.0], 20),
    ("force", [0.0, 0.0, 1.0], 30),
    ("moment", [1.0, 0.0, 0.0], 40),
    ("moment", [0.0, 1.0, 0.0], 50),
    ("moment", [0.0, 0.0, 1.0], 60),
]
_STATIC_LOADS = [
    (name, *spec) for name, spec in zip(STATIC_LOAD_NAMES, _STATIC_LOAD_DEFINITIONS)
]

# ==============================================================================
# Fixed identifiers
# ==============================================================================
_MATERIAL_ID = 1
_SPC_ID = 100
_EIGRL_ID = 1
_NUM_MODES = 10
_NODE_ID_START = 1
_ELEM_ID_START = 1
_PROP_ID_START = 1

_SCRIPT_DIR = Path(__file__).parent
_BDF_OUT_DIR = _SCRIPT_DIR / "nastran_input_files"


def computeRectangleSectionProperties(width, depth):
    """Compute area, bending inertias, and torsion constant of a solid rectangle.

    Parameters
    ----------
    width : float
        Section dimension in the local-z direction.
    depth : float
        Section dimension in the local-y direction (the through-thickness
        dimension, per the TACS thickness-along-local-y convention).

    Returns
    -------
    tuple of float
        (area, i1, i2, j) where i1 = Izz (bending in local-y), i2 = Iyy
        (bending in local-z), and j is the torsion constant from Roark's
        formula for a solid rectangular bar.
    """
    area = width * depth
    i1 = width * depth**3 / 12.0  # Izz = int(y^2 dA), bending in the local-y direction
    i2 = depth * width**3 / 12.0  # Iyy = int(z^2 dA), bending in the local-z direction

    a = 0.5 * max(width, depth)
    b = 0.5 * min(width, depth)
    j = (a * b**3) * (16.0 / 3.0 - 3.36 * (b / a) * (1.0 - (b**4) / (12.0 * a**4)))
    return area, i1, i2, j


def addPropertyCard(
    model, element, prop, section, propertyId, widthA, depthA, widthB, depthB, features
):
    """Add a single per-element beam property card to the model.

    Parameters
    ----------
    model : BDF
        pyNastran model to add the card to.
    element : str
        Element type ("CBAR" or "CBEAM").
    prop : str
        Property card type ("PBAR", "PBARL", "PBEAM", or "PBEAML").
    section : str or None
        Section type for PBARL/PBEAML ("BAR", "ROD", "TUBE", "TUBE2", or "T");
        None otherwise. ("T" is used only when building negative-test BDFs.)
    propertyId : int
        ID for this property card.
    widthA, depthA : float
        Rectangular section dimensions at the element's end-A. For TUBE/T
        sections these act as scale factors via the depthA/DEPTH1 ratio.
    widthB, depthB : float
        Rectangular section dimensions at the element's end-B.
    features : set of str
        Active feature names for this configuration.
    """
    isNsmActive = "nsm" in features
    nsm = _NSM if isNsmActive else 0.0

    # Scale factor used to taper the non-rectangular sections. depthA/DEPTH1 is
    # 1.0 at the root and 0.5 at the tip when the taper feature is active.
    scaleA = depthA / DEPTH1
    scaleB = depthB / DEPTH1

    if prop == "PBAR":
        # CBAR/PBAR is uniform by design (no taper rung), so end-A dims apply.
        area, i1, i2, j = computeRectangleSectionProperties(widthA, depthA)
        model.add_pbar(
            propertyId,
            _MATERIAL_ID,
            A=area,
            i1=i1,
            i2=i2,
            i12=0.0,
            j=j,
            nsm=nsm,
            k1=SHEAR_CORRECTION,
            k2=SHEAR_CORRECTION,
        )

    elif prop == "PBARL":
        # PBARL is uniform (no taper rung), so end-A dims apply.
        if section == "BAR":
            # BAR dims = [depth (local-z), width (local-y)] per Nastran convention.
            dims = [depthA, widthA]
        elif section == "TUBE":
            dims = [TUBE_OUTER_R, TUBE_INNER_R]
        elif section == "TUBE2":
            dims = [TUBE2_OUTER_R, TUBE2_WALL]
        elif section == "ROD":
            dims = [ROD_R1]
        elif section == "T":
            dims = list(T_DIMS1)
        else:
            raise ValueError(f"Unsupported PBARL section {section!r}")
        model.add_pbarl(propertyId, _MATERIAL_ID, section, dims, nsm=nsm)

    elif prop == "PBEAM":
        # Explicit A/I/J at each of the element's two stations (supports taper).
        areas, i1s, i2s, js = [], [], [], []
        for width, depth in [(widthA, depthA), (widthB, depthB)]:
            area, i1, i2, j = computeRectangleSectionProperties(width, depth)
            areas.append(area)
            i1s.append(i1)
            i2s.append(i2)
            js.append(j)

        n1 = _N1 if "n" in features else 0.0
        n2 = _N2 if "n" in features else 0.0
        m1 = _M1 if "m" in features else 0.0
        m2 = _M2 if "m" in features else 0.0

        model.add_pbeam(
            propertyId,
            _MATERIAL_ID,
            xxb=[0.0, 1.0],
            so="NO",
            area=areas,
            i1=i1s,
            i2=i2s,
            i12=[0.0, 0.0],
            j=js,
            nsm=[nsm, nsm],
            n1a=n1,
            n2a=n2,
            n1b=n1,
            n2b=n2,
            m1a=m1,
            m2a=m2,
            m1b=m1,
            m2b=m2,
            k1=SHEAR_CORRECTION,
            k2=SHEAR_CORRECTION,
        )

    elif prop == "PBEAML":
        # Linear dimension taper between the element's two stations.
        if section == "BAR":
            dimsA = [depthA, widthA]
            dimsB = [depthB, widthB]
        elif section == "TUBE":
            dimsA = [TUBE_OUTER_R * scaleA, TUBE_INNER_R * scaleA]
            dimsB = [TUBE_OUTER_R * scaleB, TUBE_INNER_R * scaleB]
        elif section == "TUBE2":
            dimsA = [TUBE2_OUTER_R * scaleA, TUBE2_WALL * scaleA]
            dimsB = [TUBE2_OUTER_R * scaleB, TUBE2_WALL * scaleB]
        elif section == "ROD":
            dimsA = [ROD_R1 * scaleA]
            dimsB = [ROD_R1 * scaleB]
        elif section == "T":
            dimsA = [dim * scaleA for dim in T_DIMS1]
            dimsB = [dim * scaleB for dim in T_DIMS1]
        else:
            raise ValueError(f"Unsupported PBEAML section {section!r}")
        model.add_pbeaml(
            propertyId,
            _MATERIAL_ID,
            section,
            xxb=[0.0, 1.0],
            dims=[dimsA, dimsB],
            nsm=[nsm, nsm],
        )

    else:
        raise ValueError(f"Unsupported property type {prop!r}")


def addElementCard(
    model, element, elementId, propertyId, nodeA, nodeB, orientationVec, waOffset
):
    """Add a single CBAR or CBEAM element card to the model."""
    if element == "CBAR":
        model.add_cbar(
            elementId,
            propertyId,
            [nodeA, nodeB],
            list(orientationVec),
            g0=None,
            wa=list(waOffset),
            wb=list(waOffset),
        )
    elif element == "CBEAM":
        model.add_cbeam(
            elementId,
            propertyId,
            [nodeA, nodeB],
            list(orientationVec),
            g0=None,
            wa=list(waOffset),
            wb=list(waOffset),
        )
    else:
        raise ValueError(f"Unsupported element type {element!r}")


def buildModel(element, prop, section, features, addEigrl):
    """Build a complete pyNastran model for one configuration.

    Parameters
    ----------
    element, prop, section : str / str / (str or None)
        Element type, property card type, and section type.
    features : set of str
        Active feature names for this configuration.
    addEigrl : bool
        When True, add the EIGRL modal-extraction card (for the SOL 103 model).

    Returns
    -------
    BDF
        The populated model. Loads, SPC, and PARAMs are included; the EIGRL is
        added only when ``addEigrl`` is True.
    """
    model = BDF(debug=False)
    model.add_mat1(_MATERIAL_ID, E, G, NU, rho=RHO, comment="Aluminium 7075")

    # Orientation vector: rotated when the rotation feature is active.
    orientationVec = (
        _ROTATION_VEC if "rotation" in features else _DEFAULT_ORIENTATION_VEC
    )

    # Shear-centre offset: applied via WA/WB when the wa feature is active.
    waOffset = _WA_OFFSET if "wa" in features else np.zeros(3)

    # Nodes along the beam axis (basic-x).
    xNodes = np.linspace(0.0, BEAM_LENGTH, NUM_ELEMENTS + 1)
    nodeIds = []
    for ii, xNode in enumerate(xNodes):
        nodeId = _NODE_ID_START + ii
        model.add_grid(nodeId, [xNode, 0.0, 0.0])
        nodeIds.append(nodeId)

    # Per-element tapered rectangular dimensions. When taper is inactive the
    # end-B dimensions equal the end-A dimensions, giving a uniform beam.
    isTaperActive = "taper" in features
    widthTip = WIDTH2 if isTaperActive else WIDTH1
    depthTip = DEPTH2 if isTaperActive else DEPTH1
    widths = np.linspace(WIDTH1, widthTip, NUM_ELEMENTS + 1)
    depths = np.linspace(DEPTH1, depthTip, NUM_ELEMENTS + 1)

    # One property card per element so taper can vary station-to-station.
    for ii in range(NUM_ELEMENTS):
        elementId = _ELEM_ID_START + ii
        propertyId = _PROP_ID_START + ii
        addPropertyCard(
            model,
            element,
            prop,
            section,
            propertyId,
            widths[ii],
            depths[ii],
            widths[ii + 1],
            depths[ii + 1],
            features,
        )
        addElementCard(
            model,
            element,
            elementId,
            propertyId,
            nodeIds[ii],
            nodeIds[ii + 1],
            orientationVec,
            waOffset,
        )

    # Fix the root node in all 6 DOFs.
    model.add_spc1(_SPC_ID, "123456", [nodeIds[0]])

    # Tip loads (all 6 static load cases live in every SOL 101 model).
    tipNodeId = nodeIds[-1]
    for _name, loadType, direction, loadId in _STATIC_LOADS:
        if loadType == "force":
            model.add_force(loadId, tipNodeId, FORCE_MAG, np.array(direction))
        else:
            model.add_moment(loadId, tipNodeId, MOMENT_MAG, np.array(direction))

    if addEigrl:
        model.add_eigrl(_EIGRL_ID, nd=_NUM_MODES)

    model.add_param("COUPMASS", 1)
    model.add_param("POST", 1)

    return model


def buildStaticDeck(title, log):
    """Build the SOL 101 case-control deck: one SUBCASE per static load."""
    deck = CaseControlDeck([], log=log)
    deck.add_parameter_to_global_subcase(f"TITLE = {title}")
    deck.add_parameter_to_global_subcase(f"SPC = {_SPC_ID}")
    deck.add_parameter_to_global_subcase("DISP = ALL")
    deck.add_parameter_to_global_subcase("STRESS = ALL")
    for subcaseId, (name, _type, _direction, loadId) in enumerate(
        _STATIC_LOADS, start=1
    ):
        deck.create_new_subcase(subcaseId)
        deck.add_parameter_to_local_subcase(subcaseId, f"SUBTITLE = {name}")
        deck.add_parameter_to_local_subcase(subcaseId, f"LOAD = {loadId}")
    return deck


def buildModalDeck(title, log):
    """Build the SOL 103 case-control deck: a single modal SUBCASE."""
    deck = CaseControlDeck([], log=log)
    deck.add_parameter_to_global_subcase(f"TITLE = {title}")
    deck.add_parameter_to_global_subcase(f"SPC = {_SPC_ID}")
    deck.add_parameter_to_global_subcase("DISP = ALL")
    deck.create_new_subcase(1)
    deck.add_parameter_to_local_subcase(1, "SUBTITLE = modal")
    deck.add_parameter_to_local_subcase(1, f"METHOD = {_EIGRL_ID}")
    return deck


def writeBdf(path, sol, title, deckBuilder, model):
    """Write a self-contained BDF using pyNastran's writer.

    Parameters
    ----------
    path : Path
        Output file path.
    sol : int
        Solution sequence number (101 or 103).
    title : str
        TITLE line for the case-control deck.
    deckBuilder : callable
        ``(title, log) -> CaseControlDeck`` factory for this BDF's case-control.
    model : BDF
        Populated model whose bulk data is written.
    """
    # NASTRAN BARMASS=1 makes CBAR's mass matrix include torsional rotary inertia,
    # matching TACS Beam2's mass formulation. CBEAM ignores it (no-op there).
    model.system_command_lines = ["NASTRAN BARMASS=1"]
    model.sol = sol
    model.case_control_deck = deckBuilder(title, model.log)
    model.write_bdf(path, size=8, enddata=True, close=True, write_header=False)


def main():
    """Generate all 88 BDFs into the nastran_input_files/ subdirectory."""
    _BDF_OUT_DIR.mkdir(parents=True, exist_ok=True)
    numWritten = 0
    for element, prop, section, stem, features in iterCases():
        staticModel = buildModel(element, prop, section, features, addEigrl=False)
        writeBdf(
            _BDF_OUT_DIR / f"{stem}_sol101.bdf",
            sol=101,
            title=f"{stem} static analysis",
            deckBuilder=buildStaticDeck,
            model=staticModel,
        )
        numWritten += 1

        modalModel = buildModel(element, prop, section, features, addEigrl=True)
        writeBdf(
            _BDF_OUT_DIR / f"{stem}_sol103.bdf",
            sol=103,
            title=f"{stem} modal analysis",
            deckBuilder=buildModalDeck,
            model=modalModel,
        )
        numWritten += 1

    print(f"Generated {numWritten} BDF files in {_BDF_OUT_DIR}")


if __name__ == "__main__":
    main()
