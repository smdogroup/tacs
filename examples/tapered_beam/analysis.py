"""
==============================================================================
Modal analysis of a tapered beam using TACS
==============================================================================
@File    :   analysis.py
@Date    :   2025/09/29
@Author  :   Alasdair Christison Gray
@Description :
"""

# ==============================================================================
# Standard Python modules
# ==============================================================================
import argparse
import os
from pathlib import Path

# ==============================================================================
# External Python modules
# ==============================================================================
import numpy as np

# ==============================================================================
# Extension modules
# ==============================================================================
import tacs
from tacs import pyTACS


def solveBdfFile(bdfPath: Path, outputDir: Path) -> None:
    """Solve all TACS problems in a BDF file and write results."""
    baseName = bdfPath.stem
    feaSolver = pyTACS(str(bdfPath))
    feaSolver.initialize()

    problems = feaSolver.createTACSProbsFromBDF()
    for problem in problems.values():
        problem.setOption("printLevel", 2)
        problem.setOption("writeSolution", True)
        problem.solve()
        problem.writeSolution(baseName=baseName, outputDir=str(outputDir))

        if isinstance(problem, tacs.problems.modal.ModalProblem):
            # Write out all available frequencies to a csv file.
            numModes = problem.getNumEigs()
            freqs = np.zeros((numModes, 1))
            for ii in range(numModes):
                eigVal, _ = problem.getVariables(ii)
                freqs[ii] = np.sqrt(eigVal) / (2.0 * np.pi)

            freqFile = outputDir / f"{baseName}_frequencies.csv"
            np.savetxt(freqFile, freqs, delimiter=",", header="Frequency (Hz)")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Tapered Beam SOL101/SOL103 Analysis")
    parser.add_argument(
        "-i",
        "--input-base",
        type=str,
        default="tapered_beam_pbeam",
        help="Base name for generated input files (expects <base>_sol101.bdf and <base>_sol103.bdf).",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        default="output",
        help="Output directory for results",
    )
    args = parser.parse_args()

    outputDir = Path(args.output)
    outputDir.mkdir(parents=True, exist_ok=True)

    inputBase = Path(args.input_base)
    workingDir = Path(__file__).resolve().parent
    sol101Bdf = workingDir / f"{inputBase.name}_sol101.bdf"
    sol103Bdf = workingDir / f"{inputBase.name}_sol103.bdf"

    for bdfPath in [sol101Bdf, sol103Bdf]:
        if not bdfPath.exists():
            raise FileNotFoundError(
                f"Could not find expected analysis input file: {bdfPath}"
            )

    for bdfPath in [sol101Bdf, sol103Bdf]:
        solveBdfFile(bdfPath, outputDir)
