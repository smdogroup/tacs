"""
Helpers for translating Nastran beam property cards to TACS constitutive inputs.
"""

import numpy as np
import pyNastran.bdf as pn

import tacs.constitutive
from tacs.utilities import Error

# Ceiling used to truncate PBAR/PBEAM shear stiffness factors (k1/k2); see _truncateShearStiffnessFactor.
SHEAR_STIFFNESS_FACTOR_CEILING = 1e3


def _truncateShearStiffnessFactor(k):
    """Truncate a PBAR/PBEAM shear stiffness factor (k1/k2) to ``SHEAR_STIFFNESS_FACTOR_CEILING``.

    pyNastran defaults PBAR/PBEAM shear stiffness factors (k1/k2) to 1e8, which
    can lead to scaling issues in the stiffness matrix. We truncate this value
    to prevent this.

    Parameters
    ----------
    k : float or None
        The raw shear stiffness factor from the property card. ``None`` is
        treated as an unset/defaulted value.

    Returns
    -------
    float
        ``k`` unchanged when it is finite and at or below
        ``SHEAR_STIFFNESS_FACTOR_CEILING``, otherwise the ceiling value.
    """
    return (
        SHEAR_STIFFNESS_FACTOR_CEILING
        if k is None or k > SHEAR_STIFFNESS_FACTOR_CEILING
        else k
    )


def tubeBeamDims(sectionType, dims):
    """Map a circular Nastran section's dims to IsoTubeBeamConstitutive inputs.

    The math is elementwise, so ``dims`` may hold scalars (PBARL, a 1-D dim
    list) or per-station arrays (PBEAML, pass ``propInfo.dim.T`` so ``dims[i]``
    is the i-th dimension across all stations).

    The per-section dim conventions are:

    - ROD is a solid circle (a tube with inner diameter zero); ``dims[0]`` is the radius.
    - TUBE dims are ``[r_outer, r_inner]``.
    - TUBE2 dims are ``[r_outer, wall_thickness]``.

    Parameters
    ----------
    sectionType : str
        The circular section type, one of ``"ROD"``, ``"TUBE"``, or ``"TUBE2"``.
    dims : sequence of float or ndarray
        The section dimensions, either as scalars or per-station arrays (see
        above).

    Returns
    -------
    innerDiameter : float or ndarray
        The inner diameter of the section (zero for a solid ROD).
    wallThickness : float or ndarray
        The wall thickness of the section (the full radius for a solid ROD).

    Raises
    ------
    ValueError
        If ``sectionType`` is not one of the supported circular types.
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
        raise ValueError(f"tubeBeamDims: non-circular section type '{sectionType}'")
    return innerDiameter, wallThickness


def cowperHollowCircleShearFactor(innerDiameter, outerDiameter, nu):
    """Cowper's (1966) transverse-shear correction factor for a hollow circular section.

    The factor depends on the section only through the inner/outer diameter
    ratio ``m = innerDiameter / outerDiameter``. It reduces to the solid-circle
    value ``6(1+nu)/(7+6nu)`` as ``m -> 0`` (solid rod) and to the thin-wall
    value ``2(1+nu)/(4+3nu)`` as ``m -> 1``.

    Parameters
    ----------
    innerDiameter : float or ndarray
        The inner diameter of the section (zero for a solid rod).
    outerDiameter : float or ndarray
        The outer diameter of the section.
    nu : float
        Poisson's ratio of the material.

    Returns
    -------
    float or ndarray
        The transverse-shear correction factor.
    """
    m = innerDiameter / outerDiameter
    return (
        6.0
        * (1.0 + nu)
        * (1.0 + m**2) ** 2
        / ((7.0 + 6.0 * nu) * (1.0 + m**2) ** 2 + (20.0 + 12.0 * nu) * m**2)
    )


def shearCentreOffset(elem, bdfInfo):
    """Extract the shear-centre (WA/WB) offset of a Nastran CBAR/CBEAM element.

    Parameters
    ----------
    elem : pyNastran element card
        The CBAR/CBEAM element card whose WA/WB offsets are read.
    bdfInfo : pyNastran.bdf.BDF
        The pyNastran BDF object (needed by ``elem.get_axes``).

    Returns
    -------
    shearCenterYOffset : float
        The averaged WA/WB offset vector projected onto the element's local y section axis.
    shearCenterZOffset : float
        The averaged WA/WB offset vector projected onto the element's local z section axis.
    """
    _, (_, _, yElem, zElem, wa, wb) = elem.get_axes(bdfInfo)
    # Take the average of the offset vectors at either end of bar
    offset_vector = (wa + wb) / 2.0
    # Project the offset vector onto the local section axes
    shearCenterYOffset = np.dot(yElem, offset_vector)
    shearCenterZOffset = np.dot(zElem, offset_vector)
    return shearCenterYOffset, shearCenterZOffset


def _hasShearCenterOffset(shearCenterYOffset, shearCenterZOffset):
    """Return whether a beam section has a non-zero shear-centre (WA/WB) offset.

    Single source of truth for the offset/no-offset branch decision shared by
    the PBARL and PBEAML translators, so the flag can never disagree with the
    offset values it is derived from.

    Parameters
    ----------
    shearCenterYOffset : float
        The shear-centre y offset projected onto the local section axis.
    shearCenterZOffset : float
        The shear-centre z offset projected onto the local section axis.

    Returns
    -------
    bool
        True if either projected offset is non-zero.
    """
    return shearCenterYOffset != 0.0 or shearCenterZOffset != 0.0


def averageStationProps(props, xxb):
    """Collapse a dict of per-station beam property arrays into single averaged values.

    ``xxb`` spans ``[0, 1]``, so ``np.trapezoid(value, xxb)`` naturally computes the length-average of each property
    along the element.

    Parameters
    ----------
    props : dict of {str: ndarray}
        A mapping from property name to its per-station values.
    xxb : ndarray
        ``propInfo.xxb``, the station fraction-of-length coordinates spanning ``[0, 1]``.

    Returns
    -------
    dict of {str: float}
        The same keys as ``props``, each mapped to its length-average along
        the beam.
    """
    if len(xxb) == 1:
        return {name: value[0] for name, value in props.items()}
    return {name: np.trapezoid(value, xxb) for name, value in props.items()}


def _translatePBAR(propInfo, mat, shearCenterYOffset, shearCenterZOffset):
    """Translate a PBAR card (section constants given directly, constant along the element).

    Parameters
    ----------
    propInfo : pyNastran PBAR card
        The property card to translate.
    mat : tacs.constitutive.MaterialProperties
        The TACS material properties object for this property's material.
    shearCenterYOffset : float
        The shear-centre y offset from the nodes, used for all section
        reference points.
    shearCenterZOffset : float
        The shear-centre z offset from the nodes, used for all section
        reference points.

    Returns
    -------
    tacs.constitutive.BasicBeamConstitutive
        The translated beam constitutive object.
    """
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


def _translatePBARL(propInfo, mat, shearCenterYOffset, shearCenterZOffset):
    """Translate a PBARL card (a standard cross-section type, constant along the element).

    The target constitutive class depends on the section type and whether a
    shear-centre offset is present:

    - BAR sections map to ``IsoRectangleBeamConstitutive``.
    - TUBE/TUBE2 without an offset map to ``IsoTubeBeamConstitutive``.
    - ROD (any case) and offset tubes map to ``BasicBeamConstitutive`` with a Cowper (1966) shear-correction factor computed here.

    Parameters
    ----------
    propInfo : pyNastran PBARL card
        The property card to translate.
    mat : tacs.constitutive.MaterialProperties
        The TACS material properties object for this property's material.
    shearCenterYOffset : float
        The shear-centre y offset from the nodes.
    shearCenterZOffset : float
        The shear-centre z offset from the nodes.

    Returns
    -------
    tacs.constitutive.IsoRectangleBeamConstitutive or \
tacs.constitutive.IsoTubeBeamConstitutive or \
tacs.constitutive.BasicBeamConstitutive
        The translated beam constitutive object, chosen per the rules above.

    Raises
    ------
    tacs.utilities.Error
        If the PBARL section type is not one of BAR, ROD, TUBE, or TUBE2.
    """
    nsm = propInfo.nsm
    hasShearCenterOffset = _hasShearCenterOffset(shearCenterYOffset, shearCenterZOffset)
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
        innerDiameter, wallThickness = tubeBeamDims(propInfo.Type, propInfo.dim)
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
        innerDiameter, wallThickness = tubeBeamDims(propInfo.Type, propInfo.dim)
        outerDiameter = innerDiameter + 2.0 * wallThickness
        nu = propInfo.mid_ref.nu
        kShear = cowperHollowCircleShearFactor(innerDiameter, outerDiameter, nu)
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
            f"PBARL property card {propInfo.pid} has '{propInfo.Type}' section type. "
            "This is unsupported: TACS only supports BAR, ROD, TUBE, and TUBE2, "
            "because pyNastran does not compute a correct torsion constant J for "
            "other section types.",
        )


def _translatePBEAM(propInfo, mat, shearCenterYOffset, shearCenterZOffset):
    """Translate a PBEAM card to a TACS constitutive object.

    Section constants and offsets may be given per station; they are averaged
    along the element. Nastran's neutral-axis and NSM offsets are relative to
    the shear centre, so the shear-centre offset is added back to express them
    in TACS' node-relative convention.

    Parameters
    ----------
    propInfo : pyNastran PBEAM card
        The property card to translate.
    mat : tacs.constitutive.MaterialProperties
        The TACS material properties object for this property's material.
    shearCenterYOffset : float
        The shear-centre y offset from the nodes.
    shearCenterZOffset : float
        The shear-centre z offset from the nodes.

    Returns
    -------
    tacs.constitutive.BasicBeamConstitutive
        The translated beam constitutive object.
    """
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


def _translatePBEAML(propInfo, mat, shearCenterYOffset, shearCenterZOffset):
    """Translate a PBEAML card (a standard cross-section type, tapered across stations).

    The per-station section properties are averaged along the element before
    being passed to the constitutive model. BAR sections map to
    ``IsoRectangleBeamConstitutive`` and circular sections without a
    shear-centre offset map to ``IsoTubeBeamConstitutive``; circular sections
    with an offset are unsupported (see Raises).

    Parameters
    ----------
    propInfo : pyNastran PBEAML card
        The property card to translate.
    mat : tacs.constitutive.MaterialProperties
        The TACS material properties object for this property's material.
    shearCenterYOffset : float
        The shear-centre y offset from the nodes.
    shearCenterZOffset : float
        The shear-centre z offset from the nodes.

    Returns
    -------
    tacs.constitutive.BasicBeamConstitutive or \
        tacs.constitutive.IsoRectangleBeamConstitutive or \
        tacs.constitutive.IsoTubeBeamConstitutive
        The translated beam constitutive object.

    Raises
    ------
    tacs.utilities.Error
        If the section is a circular type with a shear-centre offset, or if the PBEAML section type is not one
        of BAR, ROD, TUBE, or TUBE2. Both these scenarios would require computing the torsion constant J, to pass to
        BasicBeamConstitutive, which pyNastran does not do for PBEAML.
    """
    sectionType = propInfo.beam_type
    hasShearCenterOffset = _hasShearCenterOffset(shearCenterYOffset, shearCenterZOffset)
    sectionProps = {}
    if sectionType == "BAR":
        sectionProps["w"] = propInfo.dim[:, 0]
        sectionProps["t"] = propInfo.dim[:, 1]
        sectionProps["nsm"] = propInfo.nsm
        # Normalize the offsets by the per-station section dimensions to get
        # non-dimensional offsets for TACS (identical to the PBARL BAR branch;
        # shearCenterZ/YOffset are already the offset projected onto the local
        # section axes).
        sectionProps["wOffset"] = -shearCenterZOffset / sectionProps["w"]
        sectionProps["tOffset"] = -shearCenterYOffset / sectionProps["t"]
        conType = tacs.constitutive.IsoRectangleBeamConstitutive
    elif sectionType in ("ROD", "TUBE", "TUBE2") and (not hasShearCenterOffset):
        # Circular sections without a shear-center offset go to
        # IsoTubeBeamConstitutive, which computes its own (correct)
        # J. dims are per-station, so transpose so the helper sees
        # one dimension across all stations.
        innerDiameter, wallThickness = tubeBeamDims(sectionType, propInfo.dim.T)
        sectionProps["d"] = innerDiameter
        sectionProps["t"] = wallThickness
        sectionProps["nsm"] = propInfo.nsm
        conType = tacs.constitutive.IsoTubeBeamConstitutive

    elif sectionType in ("ROD", "TUBE", "TUBE2"):
        raise Error(
            "pyTACS",
            f"PBEAML property card {propInfo.pid} has '{sectionType}' section type and a "
            "non-zero shear-center offset (WA/WB). This is unsupported: "
            "IsoTubeBeamConstitutive cannot carry an offset, and pyNastran does not "
            "compute a torsion constant J for PBEAML, so there is no J for the "
            "BasicBeamConstitutive fallback.",
        )

    else:
        raise Error(
            "pyTACS",
            f"PBEAML property card {propInfo.pid} has '{sectionType}' section type. "
            "This is unsupported: TACS supports BAR, ROD, TUBE, and TUBE2 (circular "
            "types without WA/WB offsets), because pyNastran does not compute a correct "
            "J for other section types.",
        )

    # Whatever properties we're going to pass to the TACS
    # constitutive model, average them along the element
    sectionProps = averageStationProps(sectionProps, propInfo.xxb)

    return conType(mat, **sectionProps)


def beamPropertyToConstitutive(propInfo, elem0, bdfInfo, mat):
    """Translate a Nastran beam property card to a TACS beam constitutive object.

    Dispatches to the per-card translator (PBAR/PBARL/PBEAM/PBEAML) after
    extracting the shear-centre offset from the associated element card.

    Parameters
    ----------
    propInfo : pyNastran property card
        The beam property card to translate (PBAR, PBARL, PBEAM, or PBEAML).
    elem0 : pyNastran element card
        The first element card referencing this property (a CBAR/CBEAM), used
        for the WA/WB shear-centre offset.
    bdfInfo : pyNastran.bdf.BDF
        The pyNastran BDF object (needed by ``elem0.get_axes``).
    mat : tacs.constitutive.MaterialProperties
        The TACS material properties object for this property's material.

    Returns
    -------
    tacs.constitutive.BasicBeamConstitutive or \
        tacs.constitutive.IsoRectangleBeamConstitutive or tacs.constitutive.IsoTubeBeamConstitutive
        The translated beam constitutive object.

    Raises
    ------
    tacs.utilities.Error
        For unsupported card, section, or offset combinations.
    """
    # Get shear center offset from the associated element card
    shearCenterYOffset, shearCenterZOffset = shearCentreOffset(elem0, bdfInfo)

    if propInfo.type == "PBAR":  # Nastran bar
        return _translatePBAR(propInfo, mat, shearCenterYOffset, shearCenterZOffset)

    elif propInfo.type == "PBARL":  # Nastran bar w/ cross-section
        return _translatePBARL(propInfo, mat, shearCenterYOffset, shearCenterZOffset)

    elif propInfo.type == "PBEAM":
        return _translatePBEAM(propInfo, mat, shearCenterYOffset, shearCenterZOffset)

    elif propInfo.type == "PBEAML":
        return _translatePBEAML(propInfo, mat, shearCenterYOffset, shearCenterZOffset)

    else:
        # Should not happen: the caller only dispatches here for PBAR/PBARL/PBEAM/PBEAML.
        raise Error(
            "pyTACS",
            f"Property card {propInfo.pid} has '{propInfo.type}' type. This is "
            "unsupported: expected one of PBAR, PBARL, PBEAM, PBEAML.",
        )
