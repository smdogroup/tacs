"""
==============================================================================
Nastran beam compatibility regression suite — shared configuration
==============================================================================
Single source of truth for the element/property/feature combinations
exercised by the suite. Imported by ``generate_inputs.py``,
``extract_nastran_refs.py``, and ``tests/integration_tests/test_nastran_beam_compat.py``.

Each entry in ``CASE_MATRIX`` maps an (element, property, section) base
configuration to an ordered list of features that get applied cumulatively:
the first generated case has no features, the next has the first feature,
the next has the first two, and so on. So an entry with N features yields
N+1 cases (baseline plus N rungs of feature accumulation).
"""

CASE_MATRIX = {
    ("CBAR", "PBAR", None): ["rotation", "wa", "nsm"],
    ("CBAR", "PBARL", "BAR"): ["rotation", "wa", "nsm"],
    ("CBAR", "PBARL", "TUBE"): ["rotation", "nsm"],
    # ("CBAR", "PBARL", "T"): ["rotation", "wa", "nsm"],  # T not supported: pyNastran J() incorrect
    # ("CBEAM", "PBEAM", None): ["taper", "rotation", "wa", "n", "nsm", "m"],
    # ("CBEAM", "PBEAML", "BAR"): ["taper", "rotation", "wa", "nsm"],
    # ("CBEAM", "PBEAML", "TUBE"): ["taper", "rotation", "nsm"],
    # ("CBEAM", "PBEAML", "T"): ["taper", "rotation", "wa", "nsm"],  # T not supported: pyNastran does not compute J
}

STATIC_LOAD_NAMES = (
    "static_fx",
    "static_fy",
    "static_fz",
    "static_mx",
    "static_my",
    "static_mz",
)


def iterCases():
    """Yield (element, prop, section, stem, features) for each case.

    ``features`` is a tuple of strings (possibly empty for the baseline rung).
    ``stem`` is the canonical file/directory stem,
    e.g. ``"cbeam_pbeam_taper_rotation_wa"``.
    """
    for (element, prop, section), featureLadder in CASE_MATRIX.items():
        parts = [element.lower(), prop.lower()]
        if section is not None:
            parts.append(section.lower())
        stemPrefix = "_".join(parts)

        for numFeatures in range(len(featureLadder) + 1):
            activeFeatures = tuple(featureLadder[:numFeatures])
            rungName = "baseline" if not activeFeatures else "_".join(activeFeatures)
            stem = f"{stemPrefix}_{rungName}"
            yield element, prop, section, stem, activeFeatures
