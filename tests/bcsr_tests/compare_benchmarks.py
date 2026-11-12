"""Compare BCSRMat benchmark results written by benchmark_bcsr_mat.py.

Usage: python compare_benchmarks.py <baseline.json> <new.json> [<new2.json> ...]

Prints a markdown table of the speedup (baseline time / new time) for every operation and block size.
When several new files are given, the fastest time of each entry is used.
"""

from __future__ import annotations

import json
import sys


def loadTimes(path: str) -> dict[str, dict[str, float]]:
    """Return the ``times`` table (operation -> block size -> seconds) of a benchmark file."""
    with open(path) as resultFile:
        return json.load(resultFile)["times"]


def main() -> None:
    """Print the speedup table."""
    if len(sys.argv) < 3:
        raise SystemExit(__doc__)
    baseline = loadTimes(sys.argv[1])
    newRuns = [loadTimes(path) for path in sys.argv[2:]]

    blockSizes = sorted(next(iter(baseline.values())), key=int)
    print("| operation | " + " | ".join(f"b={size}" for size in blockSizes) + " |")
    print("|---|" + "---|" * len(blockSizes))
    for operation, baseTimes in baseline.items():
        cells = []
        for size in blockSizes:
            newTime = min(run[operation][size] for run in newRuns)
            cells.append(f"{baseTimes[size] / newTime:.2f}x")
        print(f"| {operation} | " + " | ".join(cells) + " |")


if __name__ == "__main__":
    main()
