"""
Helpers for translating Nastran beam property cards to TACS constitutive inputs.
"""

import numpy as np
import pyNastran.bdf as pn

import tacs.constitutive
from tacs.utilities import Error

# Ceiling used to truncate PBAR/PBEAM shear stiffness factors (k1/k2); see _truncateShearStiffnessFactor.
_SHEAR_STIFFNESS_FACTOR_CEILING = 1e3


def _truncateShearStiffnessFactor(k):
    """Truncate a PBAR/PBEAM shear stiffness factor (k1/k2) to ``_SHEAR_STIFFNESS_FACTOR_CEILING``.

    pynastran defaults PBAR/PBEAM shear stiffness factors (k1/k2) to 1e8, which
    can lead to scaling issues in the stiffness matrix. We truncate this value
    to prevent this.
    """
    return (
        _SHEAR_STIFFNESS_FACTOR_CEILING
        if k is None or k > _SHEAR_STIFFNESS_FACTOR_CEILING
        else k
    )


def isoTubeBeamDims(sectionType, dims):
    """Map a circular Nastran section's dims to IsoTubeBeamConstitutive inputs.

    Returns ``(innerDiameter, wallThickness)`` for ROD/TUBE/TUBE2. The math is
    elementwise, so ``dims`` may hold scalars (PBARL, a 1-D dim list) or
    per-station arrays (PBEAML, pass ``propInfo.dim.T`` so ``dims[i]`` is the
    i-th dimension across all stations).

    - ROD is a solid circle (a tube with inner diameter zero); ``dims[0]`` is the
      radius.
    - TUBE dims are ``[r_outer, r_inner]``.
    - TUBE2 dims are ``[r_outer, wall_thickness]``.
    """
    if sectionType == "ROD":
        wallThickness = dims[0]
        innerDiameter = 0.0 * dims[0]  # scalar 0.0 or zeros array matching dims[0]
    elif sectionType == "TUBE":
        innerDiameter = 2.0 * dims[1]
        wallThickness = dims[0] - dims[1]
    elif sectionType == "TUBE2":
        wallThickness = dims[1]
        innerDiameter = 2.0 * (dims[0] - dims[1])
    else:
        raise ValueError(f"isoTubeBeamDims: non-circular section type '{sectionType}'")
    return innerDiameter, wallThickness


def cowperHollowCircleShearFactor(m, nu):
    """Cowper's (1966) transverse-shear correction factor for a hollow circular section.

    ``m`` is the inner/outer diameter ratio and ``nu`` is Poisson's ratio.
    Reduces to the solid-circle value ``6(1+nu)/(7+6nu)`` as ``m -> 0`` and to
    the thin-wall value ``2(1+nu)/(4+3nu)`` as ``m -> 1``.
    """
    return (
        6.0
        * (1.0 + nu)
        * (1.0 + m**2) ** 2
        / ((7.0 + 6.0 * nu) * (1.0 + m**2) ** 2 + (20.0 + 12.0 * nu) * m**2)
    )


def shearCentreOffset(elem0, bdfInfo):
    """Extract the shear-centre (WA/WB) offset of a Nastran CBAR/CBEAM element.

    ``elem0`` is the pyNastran element card and ``bdfInfo`` is the pyNastran
    BDF object (needed by ``elem0.get_axes``).

    Returns
    -------
    shearCenterYOffset, shearCenterZOffset : float
        The averaged WA/WB offset vector projected onto the element's local
        y/z section axes.
    hasShearCenterOffset : bool
        True if either projected offset is non-zero.
    yElem, zElem : ndarray
        The element's local y/z section axes, returned so callers that need
        to project onto per-station section axes (PBARL/PBEAML BAR sections)
        don't have to call ``get_axes`` again.
    offset_vector : ndarray
        The averaged WA/WB offset vector itself, before projection.
    """
    _, (_, _, yElem, zElem, wa, wb) = elem0.get_axes(bdfInfo)
    # Take the average of the offset vectors at either end of bar
    offset_vector = (wa + wb) / 2.0
    # Project the offset vector onto the local section axes
    shearCenterYOffset = np.dot(yElem, offset_vector)
    shearCenterZOffset = np.dot(zElem, offset_vector)
    hasShearCenterOffset = shearCenterYOffset != 0.0 or shearCenterZOffset != 0.0
    return (
        shearCenterYOffset,
        shearCenterZOffset,
        hasShearCenterOffset,
        yElem,
        zElem,
        offset_vector,
    )


def averageStationProps(props, xxb):
    """Collapse a dict of per-station beam property arrays into single averaged values.

    ``props`` is a dict ``{name: array}`` of per-station values; ``xxb`` is
    ``propInfo.xxb``, the station fraction-of-length coordinates.

    ``xxb`` spans ``[0, 1]``, so ``np.trapezoid(value, xxb)`` integrates over a
    UNIT interval: it is therefore the length-average of ``value`` along the
    beam, NOT a raw integral. If there is only one station, there is nothing
    to average and the sole value is used directly.
    """
    if len(xxb) == 1:
        return {name: value[0] for name, value in props.items()}
    return {name: np.trapezoid(value, xxb) for name, value in props.items()}


def _translatePBAR(propInfo, mat, shearCenterYOffset, shearCenterZOffset):
    """Translate a PBAR card (section constants given directly, constant along the element)."""
    area = propInfo.A
    I1 = propInfo.i1
    I2 = propInfo.i2
    # Nastran uses negative convention for POI's
    I12 = -propInfo.i12
    J = propInfo.j
    k1 = propInfo.k1
    k2 = propInfo.k2
    nsm = propInfo.nsm

    k1 = _truncateShearStiffnessFactor(k1)
    k2 = _truncateShearStiffnessFactor(k2)

    # All section reference points (mass, centroid, shear, NSM) are assumed
    # coincident and placed at the WA/WB offset from the nodes.
    return tacs.constitutive.BasicBeamConstitutive(
        mat,
        A=area,
        I33=I1,
        I22=I2,
        I23=I12,
        J=J,
        k2=k1,
        k3=k2,
        nsm=nsm,
        xm2=shearCenterYOffset,
        xm3=shearCenterZOffset,
        xc2=shearCenterYOffset,
        xc3=shearCenterZOffset,
        xk2=shearCenterYOffset,
        xk3=shearCenterZOffset,
        xnsm2=shearCenterYOffset,
        xnsm3=shearCenterZOffset,
    )


def _translatePBARL(
    propInfo, mat, shearCenterYOffset, shearCenterZOffset, hasShearCenterOffset
):
    """Translate a PBARL card (a standard cross-section type, constant along the element)."""
    nsm = propInfo.nsm
    if propInfo.Type == "BAR":
        w = propInfo.dim[0]
        t = propInfo.dim[1]
        # Normalize the offsets by the section dimensions to get non-dimensional offsets for TACS
        wOffset = -shearCenterZOffset / w
        tOffset = -shearCenterYOffset / t
        return tacs.constitutive.IsoRectangleBeamConstitutive(
            mat, w=w, t=t, tOffset=tOffset, wOffset=wOffset, nsm=nsm
        )

    elif propInfo.Type in ("TUBE", "TUBE2") and (not hasShearCenterOffset):
        # Hollow circular sections without a shear-center offset go
        # to IsoTubeBeamConstitutive, which computes its own
        # (correct) J. TUBE/TUBE2 differ only in dim convention.
        innerDiameter, wallThickness = isoTubeBeamDims(propInfo.Type, propInfo.dim)
        return tacs.constitutive.IsoTubeBeamConstitutive(
            mat, d=innerDiameter, t=wallThickness, nsm=nsm
        )

    elif propInfo.Type in ("ROD", "TUBE", "TUBE2"):
        # Solid ROD (any case) and hollow tubes with a shear-center
        # offset use BasicBeamConstitutive: IsoTubeBeamConstitutive
        # cannot carry an offset, and its shear-correction factor is
        # the thin-walled limit 2(1+nu)/(4+3*nu), which is wrong for
        # a solid ROD. pyNastran does not expose a shear factor for
        # PBARL, so we compute Cowper's (1966) hollow-circular value
        # here from the inner/outer diameter ratio m. It reduces to
        # the solid-circle value 6(1+nu)/(7+6*nu) as m -> 0 and to
        # the thin-wall value above as m -> 1. pyNastran's PBARL.J()
        # is exact for circular sections (J = polar moment I1 + I2).
        A, I1, I2, I12 = pn.cards.properties.bars._bar_areaL(
            "PBARL", propInfo.Type, propInfo.dim, propInfo
        )
        J = propInfo.J()
        innerDiameter, wallThickness = isoTubeBeamDims(propInfo.Type, propInfo.dim)
        outerDiameter = innerDiameter + 2.0 * wallThickness
        m = innerDiameter / outerDiameter
        nu = propInfo.mid_ref.nu
        kShear = cowperHollowCircleShearFactor(m, nu)
        return tacs.constitutive.BasicBeamConstitutive(
            mat,
            A=A,
            J=J,
            I33=I1,
            I22=I2,
            I23=-I12,
            k2=kShear,
            k3=kShear,
            nsm=nsm,
            xm2=shearCenterYOffset,
            xm3=shearCenterZOffset,
            xc2=shearCenterYOffset,
            xc3=shearCenterZOffset,
            xk2=shearCenterYOffset,
            xk3=shearCenterZOffset,
            xnsm2=shearCenterYOffset,
            xnsm3=shearCenterZOffset,
        )

    else:
        raise Error(
            "pyTACS",
            f"Unsupported PBARL section type '{propInfo.Type}' for property number "
            f"{propInfo.pid}. TACS supports BAR, ROD, TUBE, and TUBE2. pyNastran does "
            "not compute a correct torsion constant J for other section types.",
        )


def _translatePBEAM(propInfo, mat, shearCenterYOffset, shearCenterZOffset):
    """Translate a PBEAM card (section constants given directly, tapered across stations)."""
    area = propInfo.A
    I1 = propInfo.i1
    I2 = propInfo.i2
    # Nastran uses negative convention for POI's
    I12 = -propInfo.i12
    J = propInfo.j
    k1 = propInfo.k1
    k2 = propInfo.k2
    nsm = propInfo.nsm
    # Offsets relative to shear center
    # Y coordinate of non-structural mass for end A
    nsmYOffsetA = propInfo.m1a
    # Z coordinate of non-structural mass for end A
    nsmZOffsetA = propInfo.m2a
    # Y coordinate of non-structural mass for end B
    nsmYOffsetB = propInfo.m1b
    # Z coordinate of non-structural mass for end B
    nsmZOffsetB = propInfo.m2b
    # Y coordinate of neutral axis for end A
    neutralAxisYOffsetA = propInfo.n1a
    # Z coordinate of neutral axis for end A
    neutralAxisZOffsetA = propInfo.n2a
    # Y coordinate of neutral axis for end B
    neutralAxisYOffsetB = propInfo.n1b
    # Z coordinate of neutral axis for end B
    neutralAxisZOffsetB = propInfo.n2b

    # Average offsets at two ends, then add the shear center offset because TACS offsets are not
    # relative to shear center like Nastran's are
    nsmYOffset = (nsmYOffsetA + nsmYOffsetB) / 2 + shearCenterYOffset
    nsmZOffset = (nsmZOffsetA + nsmZOffsetB) / 2 + shearCenterZOffset
    neutralAxisYOffset = (
        neutralAxisYOffsetA + neutralAxisYOffsetB
    ) / 2 + shearCenterYOffset
    neutralAxisZOffset = (
        neutralAxisZOffsetA + neutralAxisZOffsetB
    ) / 2 + shearCenterZOffset

    k1 = _truncateShearStiffnessFactor(k1)
    k2 = _truncateShearStiffnessFactor(k2)

    averaged = averageStationProps(
        {"A": area, "I1": I1, "I2": I2, "I12": I12, "J": J, "nsm": nsm}, propInfo.xxb
    )
    area = averaged["A"]
    I1 = averaged["I1"]
    I2 = averaged["I2"]
    I12 = averaged["I12"]
    J = averaged["J"]
    nsm = averaged["nsm"]

    return tacs.constitutive.BasicBeamConstitutive(
        mat,
        A=area,
        I22=I2,
        I33=I1,
        I23=I12,
        J=J,
        k2=k1,
        k3=k2,
        nsm=nsm,
        xk2=shearCenterYOffset,
        xk3=shearCenterZOffset,
        xc2=neutralAxisYOffset,
        xc3=neutralAxisZOffset,
        xm2=neutralAxisYOffset,
        xm3=neutralAxisZOffset,
        xnsm2=nsmYOffset,
        xnsm3=nsmZOffset,
    )


def _translatePBEAML(
    propInfo,
    mat,
    shearCenterYOffset,
    shearCenterZOffset,
    hasShearCenterOffset,
    yElem,
    zElem,
    offset_vector,
):
    """Translate a PBEAML card (a standard cross-section type, tapered across stations)."""
    sectionType = propInfo.beam_type
    sectionProps = {}
    if sectionType == "BAR":
        sectionProps["w"] = propInfo.dim[:, 0]
        sectionProps["t"] = propInfo.dim[:, 1]
        sectionProps["nsm"] = propInfo.nsm
        # Project the offset vector onto the width and thickness axes
        sectionProps["wOffset"] = -np.dot(zElem, offset_vector) / sectionProps["w"]
        sectionProps["tOffset"] = -np.dot(yElem, offset_vector) / sectionProps["t"]
        conType = tacs.constitutive.IsoRectangleBeamConstitutive
    elif sectionType in ("ROD", "TUBE", "TUBE2") and (not hasShearCenterOffset):
        # Circular sections without a shear-center offset go to
        # IsoTubeBeamConstitutive, which computes its own (correct)
        # J. dims are per-station, so transpose so the helper sees
        # one dimension across all stations.
        innerDiameter, wallThickness = isoTubeBeamDims(sectionType, propInfo.dim.T)
        sectionProps["d"] = innerDiameter
        sectionProps["t"] = wallThickness
        sectionProps["nsm"] = propInfo.nsm
        conType = tacs.constitutive.IsoTubeBeamConstitutive

    elif sectionType in ("ROD", "TUBE", "TUBE2"):
        # Circular section with a shear-center offset.
        # IsoTubeBeamConstitutive cannot carry the offset, and the
        # BasicBeamConstitutive fallback needs a J that PBEAML cannot
        # provide (pyNastran's PBEAML.J() always returns None), so
        # this combination is unsupported.
        raise Error(
            "pyTACS",
            f"PBEAML section type '{sectionType}' with a shear-center (WA/WB) offset "
            f"is unsupported for property number {propInfo.pid}: IsoTubeBeamConstitutive "
            "cannot carry an offset, and pyNastran does not compute a torsion constant J "
            "for PBEAML, so there is no J for the BasicBeamConstitutive fallback.",
        )

    else:
        raise Error(
            "pyTACS",
            f"Unsupported PBEAML section type '{sectionType}' for property number "
            f"{propInfo.pid}. TACS supports BAR, ROD, TUBE, and TUBE2 (circular types "
            "without WA/WB offsets). pyNastran does not compute a correct J for other "
            "section types.",
        )

    # Whatever properties we're going to pass to the TACS
    # constitutive model, average them along the element
    sectionProps = averageStationProps(sectionProps, propInfo.xxb)

    return conType(mat, **sectionProps)


def beamPropertyToConstitutive(propInfo, elem0, bdfInfo, mat):
    """Translate a Nastran beam property card to a TACS beam constitutive object.

    propInfo : the pyNastran property card (PBAR/PBARL/PBEAM/PBEAML)
    elem0    : the first pyNastran element card referencing this property (a CBAR/CBEAM),
               used for the WA/WB shear-centre offset
    bdfInfo  : the pyNastran BDF object (needed by elem0.get_axes)
    mat      : the TACS MaterialProperties object for this property's material
    returns  : a tacs.constitutive beam constitutive object (BasicBeam / IsoRectangleBeam / IsoTubeBeam)
    raises   : tacs.utilities.Error for unsupported card/section/offset combinations
    """
    # Get shear center offset from the associated element card
    (
        shearCenterYOffset,
        shearCenterZOffset,
        hasShearCenterOffset,
        yElem,
        zElem,
        offset_vector,
    ) = shearCentreOffset(elem0, bdfInfo)

    if propInfo.type == "PBAR":  # Nastran bar
        return _translatePBAR(propInfo, mat, shearCenterYOffset, shearCenterZOffset)

    elif propInfo.type == "PBARL":  # Nastran bar w/ cross-section
        return _translatePBARL(
            propInfo, mat, shearCenterYOffset, shearCenterZOffset, hasShearCenterOffset
        )

    elif propInfo.type == "PBEAM":
        return _translatePBEAM(propInfo, mat, shearCenterYOffset, shearCenterZOffset)

    elif propInfo.type == "PBEAML":
        return _translatePBEAML(
            propInfo,
            mat,
            shearCenterYOffset,
            shearCenterZOffset,
            hasShearCenterOffset,
            yElem,
            zElem,
            offset_vector,
        )

    else:
        # Should not happen: the caller only dispatches here for PBAR/PBARL/PBEAM/PBEAML.
        raise Error(
            "pyTACS",
            f"Unsupported beam property type '{propInfo.type}'. Expected one of "
            "PBAR, PBARL, PBEAM, PBEAML.",
        )
