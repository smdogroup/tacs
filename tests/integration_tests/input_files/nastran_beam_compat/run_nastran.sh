#!/usr/bin/env bash
# Run Nastran on every BDF in this directory.
# Edit NASTRAN_CMD once for your local Nastran install path.
#
# Usage: ./run_nastran.sh
#
# Outputs land in _nastran_outputs/ (gitignored).

set -e

NASTRAN_CMD="${NASTRAN_CMD:-nastran}"   # override with: NASTRAN_CMD=/path/to/nastran ./run_nastran.sh
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
OUTPUT_DIR="${SCRIPT_DIR}/_nastran_outputs"

mkdir -p "${OUTPUT_DIR}"

for bdf in "${SCRIPT_DIR}"/*.bdf; do
    stem="$(basename "${bdf}" .bdf)"
    echo "Running Nastran on ${stem}..."
    "${NASTRAN_CMD}" "${bdf}" \
        out="${OUTPUT_DIR}/${stem}" \
        mem=2gb \
        old=no \
        news=no
done

echo "Done. Results in ${OUTPUT_DIR}/"
