"""
==============================================================================
Nastran beam compatibility regression suite — reference data extractor
==============================================================================
@File    :   extract_nastran_refs.py
@Description : Reads Nastran .op2 outputs from ``_nastran_outputs/`` and
    writes per-configuration CSV reference files consumed by the test harness.

    For every configuration the script writes:
      - static_fx.csv  …  static_mz.csv  (31 × 6 displacement arrays)
      - modal_frequencies.csv             (10 frequencies in Hz)
      - modal_mode01.csv  …  modal_mode10.csv  (31 × 6 max-normalised mode shapes)

Usage:
    cd tests/integration_tests/input_files/nastran_beam_compat
    python extract_nastran_refs.py
"""

# ==============================================================================
# Standard Python modules
# ==============================================================================
from pathlib import Path

# ==============================================================================
# External Python modules
# ==============================================================================
import numpy as np
from pyNastran.op2.op2 import OP2  # pyNastran.op2.op2 → OP2

# ==============================================================================
# Local modules
# ==============================================================================
from cases import (
    iterCases,
    STATIC_LOAD_NAMES,
)  # cases.iterCases, cases.STATIC_LOAD_NAMES

_SCRIPT_DIR = Path(__file__).parent
_NASTRAN_OUTPUTS_DIR = _SCRIPT_DIR / "nastran_outputs"


# ==============================================================================
# Helpers
# ==============================================================================


def _maxNormalize(shape: np.ndarray) -> np.ndarray:
    """Divide a mode shape by its largest-magnitude element so max(|shape|) == 1."""
    idx = np.argmax(np.abs(shape))
    return shape / shape.flatten()[idx]


# ==============================================================================
# Static extraction
# ==============================================================================


def extractStaticRefs(stem: str, refDir: Path) -> None:
    """Read the SOL 101 .op2 for *stem* and write per-load-case displacement CSVs.

    Parameters
    ----------
    stem:
        Configuration stem, e.g. ``cbar_pbar_baseline``.
    refDir:
        Output directory for this configuration's reference CSVs.
    """
    op2Path = _NASTRAN_OUTPUTS_DIR / f"{stem}_sol101.op2"
    if not op2Path.exists():
        print(f"  SKIP (not found): {op2Path}")
        return

    op2 = OP2(log=None, debug=False)
    op2.set_results("displacements")
    op2.read_op2(str(op2Path))

    # Map subcase results to static load names using the subtitle attribute.
    # Fall back to label if subtitle is absent.
    found = set()
    for subcaseId, dispObj in op2.displacements.items():
        # Attempt subtitle first, then label.
        subtitle = None
        if hasattr(dispObj, "subtitle") and dispObj.subtitle is not None:
            subtitle = dispObj.subtitle.strip().lower()
        if subtitle is None or subtitle == "":
            if hasattr(dispObj, "label") and dispObj.label is not None:
                subtitle = dispObj.label.strip().lower()

        if subtitle not in STATIC_LOAD_NAMES:
            print(
                f"  WARNING: subcase {subcaseId} subtitle {subtitle!r} not in expected names — skipping"
            )
            continue

        # data shape: (1, n_nodes, 6) for a single load step
        disps = dispObj.data[0]  # (n_nodes, 6)
        outPath = refDir / f"{subtitle}.csv"
        np.savetxt(
            str(outPath),
            disps,
            delimiter=",",
            header="u,v,w,rotx,roty,rotz",
            fmt="%.10e",
            comments="# ",
        )
        found.add(subtitle)

    missing = set(STATIC_LOAD_NAMES) - found
    if missing:
        print(f"  WARNING: {stem}: no subcase found for {sorted(missing)}")
    else:
        print(f"  OK static: {stem}")


# ==============================================================================
# Modal extraction
# ==============================================================================


def extractModalRefs(stem: str, refDir: Path) -> None:
    """Read the SOL 103 .op2 for *stem* and write frequency and mode-shape CSVs.

    Parameters
    ----------
    stem:
        Configuration stem, e.g. ``cbar_pbar_baseline``.
    refDir:
        Output directory for this configuration's reference CSVs.
    """
    op2Path = _NASTRAN_OUTPUTS_DIR / f"{stem}_sol103.op2"
    if not op2Path.exists():
        print(f"  SKIP (not found): {op2Path}")
        return

    op2 = OP2(log=None, debug=False)
    op2.set_results("eigenvectors")
    op2.read_op2(str(op2Path))

    eigvecObj = list(op2.eigenvectors.values())[0]
    shapes = eigvecObj.data  # (n_modes, n_nodes, 6)
    freqs = np.array(eigvecObj.mode_cycles, dtype=float)  # (n_modes,) in Hz

    # Write frequency CSV: 10 rows × 1 column
    freqPath = refDir / "modal_frequencies.csv"
    np.savetxt(
        str(freqPath),
        freqs.reshape(-1, 1),
        delimiter=",",
        header="Frequency (Hz)",
        fmt="%.10e",
        comments="# ",
    )

    # Write max-normalised mode shape CSVs: 1-indexed, zero-padded to 2 digits
    nModes = shapes.shape[0]
    for ii in range(nModes):
        shape = shapes[ii]  # (n_nodes, 6)
        shape = _maxNormalize(shape)
        modePath = refDir / f"modal_mode{ii + 1:02d}.csv"
        np.savetxt(
            str(modePath),
            shape,
            delimiter=",",
            header="u,v,w,rotx,roty,rotz",
            fmt="%.10e",
            comments="# ",
        )

    print(f"  OK modal:  {stem}  ({nModes} modes, {len(freqs)} freqs)")


# ==============================================================================
# Main
# ==============================================================================


def main() -> None:
    """Extract Nastran reference data for all 36 configurations."""
    print(f"Nastran outputs directory: {_NASTRAN_OUTPUTS_DIR}")
    print(f"Reference CSV output root: {_SCRIPT_DIR}")
    print()

    for _element, _prop, _section, stem, _features in iterCases():
        refDir = _SCRIPT_DIR / stem
        refDir.mkdir(parents=True, exist_ok=True)

        extractStaticRefs(stem, refDir)
        extractModalRefs(stem, refDir)

    print()
    print("Done.")


if __name__ == "__main__":
    main()
