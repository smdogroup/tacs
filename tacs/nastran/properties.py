"""
Helpers for translating Nastran beam property cards to TACS constitutive inputs.
"""


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
