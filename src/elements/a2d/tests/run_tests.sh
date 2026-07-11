#!/bin/sh
# Repeatable runner for the standalone A2D second-order op-level tests
# (feature-beam-element-methods, SPEC.md sec 6.6). These tests are
# header-only: they compile directly against src/elements/a2d/*.h +
# TACSObject.h via mpicxx, with no libtacs/Python-extension build required
# (per VALIDATION.md's E6 experiment) -- so this script, not pytest/ctest,
# is the regression net for this directory. Each test source's own header
# comment documents the same build recipe this script automates,
# including the TACSObject.h-before-a2d.h include-ordering workaround
# (a2dmatops.h/a2dmatcore.h use TacsScalar before a2dobjs.h -- the file
# that #includes TACSObject.h -- would otherwise be reached; this is a
# standalone-compilation-only issue, not a defect in a2dmatops.h itself).
#
# Usage: sh src/elements/a2d/tests/run_tests.sh
# Exits 0 if every test passes in both real and -DTACS_USE_COMPLEX builds,
# nonzero otherwise (per-test/per-mode failures are listed at the end).

set -u

SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO_ROOT=$(cd "$SCRIPT_DIR/../../../.." && pwd)
BUILD_DIR=$(mktemp -d)
trap 'rm -rf "$BUILD_DIR"' EXIT

INCLUDES="-I$REPO_ROOT/src -I$REPO_ROOT/src/elements/a2d"

fail_count=0
total_count=0
failed_list=""

for src in "$SCRIPT_DIR"/test_*.cpp; do
  name=$(basename "$src" .cpp)

  for mode in real complex; do
    total_count=$((total_count + 1))
    bin="$BUILD_DIR/${name}_${mode}"
    if [ "$mode" = "complex" ]; then
      defs="-DTACS_USE_COMPLEX"
    else
      defs=""
    fi

    if ! mpicxx -std=c++11 $defs $INCLUDES "$src" -o "$bin" 2>"$BUILD_DIR/${name}_${mode}.compile.log"; then
      echo "[COMPILE FAIL] $name ($mode)"
      cat "$BUILD_DIR/${name}_${mode}.compile.log"
      fail_count=$((fail_count + 1))
      failed_list="$failed_list $name($mode:compile)"
      continue
    fi

    if ! "$bin" > "$BUILD_DIR/${name}_${mode}.out" 2>&1; then
      echo "[RUN FAIL] $name ($mode)"
      cat "$BUILD_DIR/${name}_${mode}.out"
      fail_count=$((fail_count + 1))
      failed_list="$failed_list $name($mode:run)"
    else
      echo "[PASS] $name ($mode)"
    fi
  done
done

echo ""
echo "$((total_count - fail_count))/$total_count passed"
if [ "$fail_count" -ne 0 ]; then
  echo "FAILED:$failed_list"
  exit 1
fi
exit 0
