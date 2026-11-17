r"""
Shear / bending-moment / torque report for a wingbox under several load cases.

The script solves the coarse MDO tutorial wingbox with pyTACS for a 2.5 g
pull-up, a -1 g push-over and a lateral gust, writes one f5 file per case
and reduces the nodal loads to V, M and T along the quarter-chord line of the
wing box. The box is unswept out to the planform break at z = 1.5 m and swept
and tapered from there to the tip at z = 13.8 m, so the quarter-chord axis is
a three-point piecewise-linear line. It produces ``wingbox_vmt.pdf``, a
multi-page report with the envelope of all cases first and one page per
case, plus one CSV per case and an envelope CSV. The model is in SI units
(N, m); the report labels its axes accordingly and shows forces in kN, while
the CSV files and the printed table stay in N.

The same post-processing is available from the command line::

    python -m tacs.postprocess.vmt pullup_000.f5 pushover_000.f5 gust_000.f5 \
        --axis 1.9375 0 0  1.9375 0 1.5  7.95 0 13.8 --shear-dir 0 -1 0 \
        --num-stations 30 --output wingbox_vmt.pdf --csv wingbox_vmt.csv \
        --force-scale 1e-3 --force-unit kN --length-unit m
"""

import argparse
import os

from mpi4py import MPI

from tacs import pyTACS
from tacs.postprocess import (
    computeEnvelope,
    computeVMTCases,
    loadF5,
    plotVMTReport,
    writeVMTCsv,
    writeVMTEnvelopeCsv,
)

parser = argparse.ArgumentParser()
parser.add_argument(
    "--bdf",
    default=os.path.join(
        os.path.dirname(__file__),
        "../../tests/integration_tests/input_files/coarse_mdo_tutorial_wingbox.bdf",
    ),
)
parser.add_argument("--num-stations", type=int, default=30)
args = parser.parse_args()

comm = MPI.COMM_WORLD

# Load cases as (name, acceleration vector); the model's vertical axis is y
g = 9.81
loadCases = {
    "pullup": [0.0, -2.5 * g, 0.0],
    "pushover": [0.0, 1.0 * g, 0.0],
    "gust": [0.5 * g, -1.5 * g, 0.0],
}

# Solve every case and write its f5 file (loads and reactions are written by
# default)
FEAAssembler = pyTACS(args.bdf, comm)
FEAAssembler.initialize()
f5Files = []
for name, accel in loadCases.items():
    problem = FEAAssembler.createStaticProblem(name)
    problem.addInertialLoad(accel)
    problem.solve()
    problem.writeSolution(outputDir=".")
    f5Files.append(f"{name}_000.f5")

if comm.rank == 0:
    # Beam axis along the quarter chord of the wing box: the leading and
    # trailing spars sit at x = 1.25 and 4.0 out to the planform break at
    # z = 1.5, and at x = 7.725 and 8.625 at the tip (z = 13.8). The axis lies
    # in the box mid-plane y = 0 and the main load acts in -y.
    axisPts = [
        [1.9375, 0.0, 0.0],  # root, quarter chord of the unswept section
        [1.9375, 0.0, 1.5],  # planform break
        [7.95, 0.0, 13.8],  # tip
    ]
    shearDir = [0.0, -1.0, 0.0]

    results = computeVMTCases(
        f5Files, axisPts, shearDir=shearDir, numStations=args.num_stations
    )
    envelope = computeEnvelope(results)
    for name, result in results.items():
        print(
            f"{name:>9s}: root V = {result.V[0]:11.4e}  "
            f"M = {result.M[0]:11.4e}  T = {result.T[0]:11.4e}"
        )
        writeVMTCsv(result, f"wingbox_vmt_{name}.csv", metadata={"case": name})
    print(f"envelope: root Vmax = {envelope.Vmax[0]:11.4e} ({envelope.VmaxCase[0]})")
    writeVMTEnvelopeCsv(envelope, "wingbox_vmt_envelope.csv")

    # Label the plot axes with the model units. The scale factor converts the
    # plotted forces (and moments) from N to kN; the results themselves and
    # the CSV files above are left in model units.
    pages = plotVMTReport(
        results,
        "wingbox_vmt.pdf",
        data=loadF5(f5Files[0]),
        shearDir=shearDir,
        title="Wingbox",
        forceUnit="kN",
        lengthUnit="m",
        forceScale=1e-3,
    )
    print(f"wrote wingbox_vmt.pdf with {pages} pages")
