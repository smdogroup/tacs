"""
Benchmark of the BCSRMat block sparse matrix kernels (through ``TACS.BCSRMatrix``).

For each block size a structured problem is built (2-D grid of nodes, 9-point block stencil,
natural ordering, random block-diagonally dominant values) and each operation is timed as the
best of ``--repeats`` runs. Results are printed as a table (milliseconds per call) and optionally
written as JSON of the form ``{"meta": {...}, "times": {operation: {block size: seconds}}}``.

Example (pinned to one core)::

    taskset -c 2 python tests/bcsr_tests/benchmark_bcsr_mat.py --out bcsr-baseline.json
"""

from __future__ import annotations

import argparse
import json
import platform
import time
from collections.abc import Callable

import numpy as np

from tacs import TACS

SEED = 20261111
FULL_LEV_FILL = 1000
FILL_ESTIMATE = 10.0
STENCIL_REACH = 1  # 9-point stencil: neighbours within one node in each grid direction
NUM_VECTOR_CALLS = 20  # calls per timed sample for the cheap vector operations
SOR_OMEGA = 1.0
OPERATIONS = (
    "mult",
    "multAdd",
    "multTranspose",
    "factorFull",
    "applyFactorFull",
    "factorILU0",
    "applyFactorILU0",
    "applySOR",
    "applyLowerFactor",
    "applyUpperFactor",
    "matMultAdd",
)


def buildGridPattern(gridSize: int) -> tuple[np.ndarray, np.ndarray]:
    """Return (rowp, cols) of the 9-point stencil on a gridSize x gridSize grid in natural ordering."""
    rowp = [0]
    cols = []
    for jj in range(gridSize):
        for ii in range(gridSize):
            for jjNbr in range(
                max(jj - STENCIL_REACH, 0), min(jj + STENCIL_REACH, gridSize - 1) + 1
            ):
                for iiNbr in range(
                    max(ii - STENCIL_REACH, 0),
                    min(ii + STENCIL_REACH, gridSize - 1) + 1,
                ):
                    cols.append(jjNbr * gridSize + iiNbr)
            rowp.append(len(cols))
    return np.array(rowp), np.array(cols)


def buildGridMatrix(
    rng: np.random.Generator, bsize: int, gridSize: int
) -> TACS.BCSRMatrix:
    """Return the block-diagonally dominant random stencil matrix for one block size."""
    rowp, cols = buildGridPattern(gridSize)
    blocks = rng.uniform(-1.0, 1.0, (len(cols), bsize, bsize)).astype(TACS.dtype)
    maxRowSum = (2 * STENCIL_REACH + 1) ** 2 * bsize
    rows = np.repeat(np.arange(len(rowp) - 1), np.diff(rowp))
    blocks[rows == cols] += (maxRowSum + 1.0) * np.eye(bsize)
    return TACS.BCSRMatrix(bsize, rowp, cols, A=blocks)


def timeBest(
    func: Callable[[], None],
    repeats: int,
    numCalls: int = 1,
    setup: Callable[[], None] | None = None,
) -> float:
    """Return the best time per call of `func` in seconds; `setup` runs untimed before each sample."""
    best = float("inf")
    for _ in range(repeats):
        if setup is not None:
            setup()
        start = time.perf_counter()
        for _ in range(numCalls):
            func()
        best = min(best, (time.perf_counter() - start) / numCalls)
    return best


def timeFactorization(
    mat: TACS.BCSRMatrix, levFill: int, rhs: np.ndarray, repeats: int
) -> tuple[float, float]:
    """Return (numeric factor time, applyFactor time) for the ILU(levFill) pattern of `mat`."""
    factor = TACS.BCSRMatrix.createILU(mat, levFill, FILL_ESTIMATE)
    solution = np.zeros_like(rhs)
    factorTime = timeBest(factor.factor, repeats, setup=lambda: factor.copyValues(mat))
    applyTime = timeBest(
        lambda: factor.applyFactor(rhs, solution), repeats, NUM_VECTOR_CALLS
    )
    return factorTime, applyTime


def timeSchurSequence(mat: TACS.BCSRMatrix, repeats: int) -> dict[str, float]:
    """Time applyLowerFactor, applyUpperFactor and matMultAdd(-1) on a half/half partition with full fill."""
    parts = mat.partition(mat.getRowDim() // 2)
    matBpc, matEpc, matFpc, matSc = TACS.BCSRMatrix.createSchur(
        *parts, FULL_LEV_FILL, FILL_ESTIMATE, True
    )

    def resetValues() -> None:
        for target, source in zip((matBpc, matEpc, matFpc, matSc), parts, strict=True):
            target.copyValues(source)
        matBpc.factor()

    return {
        "applyLowerFactor": timeBest(
            lambda: matBpc.applyLowerFactor(matEpc), repeats, setup=resetValues
        ),
        "applyUpperFactor": timeBest(
            lambda: matBpc.applyUpperFactor(matFpc), repeats, setup=resetValues
        ),
        "matMultAdd": timeBest(
            lambda: matSc.matMultAdd(-1.0, matFpc, matEpc), repeats, setup=resetValues
        ),
    }


def benchmarkBlockSize(bsize: int, gridSize: int, repeats: int) -> dict[str, float]:
    """Return the time in seconds of every operation for one block size."""
    rng = np.random.default_rng(SEED + bsize)
    mat = buildGridMatrix(rng, bsize, gridSize)
    size = bsize * mat.getRowDim()
    x = rng.uniform(-1.0, 1.0, size).astype(TACS.dtype)
    y = np.zeros(size, dtype=TACS.dtype)
    times = {
        "mult": timeBest(lambda: mat.mult(x, y), repeats, NUM_VECTOR_CALLS),
        "multAdd": timeBest(lambda: mat.multAdd(x, y, y), repeats, NUM_VECTOR_CALLS),
        "multTranspose": timeBest(
            lambda: mat.multTranspose(x, y), repeats, NUM_VECTOR_CALLS
        ),
    }
    times["factorFull"], times["applyFactorFull"] = timeFactorization(
        mat, FULL_LEV_FILL, x, repeats
    )
    times["factorILU0"], times["applyFactorILU0"] = timeFactorization(
        mat, 0, x, repeats
    )
    mat.factorDiag()
    y[:] = 0.0
    times["applySOR"] = timeBest(
        lambda: mat.applySOR(x, y, SOR_OMEGA, 1), repeats, NUM_VECTOR_CALLS
    )
    times.update(timeSchurSequence(mat, repeats))
    return times


def printTable(times: dict[str, dict[str, float]], blockSizes: list[int]) -> None:
    """Print the times in milliseconds, one row per operation and one column per block size."""
    print(
        f"{'operation [ms]':<18}"
        + "".join(f"{'b=' + str(bsize):>11}" for bsize in blockSizes)
    )
    for operation in OPERATIONS:
        print(
            f"{operation:<18}"
            + "".join(
                f"{1e3 * times[operation][str(bsize)]:>11.4f}" for bsize in blockSizes
            )
        )


def main() -> None:
    """Parse the command line, run the benchmark and write the results."""
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--grid", type=int, default=40, help="nodes per grid side (default: 40)"
    )
    parser.add_argument(
        "--repeats",
        type=int,
        default=5,
        help="timed samples per operation; the best is kept",
    )
    parser.add_argument(
        "--block-sizes",
        type=int,
        nargs="+",
        default=list(range(1, 9)),
        help="block sizes to run",
    )
    parser.add_argument("--out", type=str, default=None, help="JSON output file")
    args = parser.parse_args()

    times: dict[str, dict[str, float]] = {operation: {} for operation in OPERATIONS}
    for bsize in args.block_sizes:
        for operation, seconds in benchmarkBlockSize(
            bsize, args.grid, args.repeats
        ).items():
            times[operation][str(bsize)] = seconds
    printTable(times, args.block_sizes)

    if args.out is not None:
        meta = {
            "grid": args.grid,
            "numBlockRows": args.grid**2,
            "repeats": args.repeats,
            "blockSizes": args.block_sizes,
            "dtype": np.dtype(TACS.dtype).name,
            "machine": platform.machine(),
            "date": time.strftime("%Y-%m-%dT%H:%M:%S"),
            "units": "seconds per call, best of repeats",
        }
        with open(args.out, "w") as outFile:
            json.dump({"meta": meta, "times": times}, outFile, indent=2)


if __name__ == "__main__":
    main()
