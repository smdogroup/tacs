# Nastran Beam Compatibility References

These files support the regression suite in `tests/integration_tests/test_nastran_beam_compat.py`.

## Regenerating reference CSVs (requires a local Nastran install)

```bash
cd tests/integration_tests/input_files/nastran_beam_compat

# Step 1 — Regenerate the 72 BDFs (already committed; only needed after model changes)
python generate_inputs.py

# Step 2 — Run Nastran on every BDF (edit run_nastran.sh once for your local install)
./run_nastran.sh

# Step 3 — Extract Nastran results to CSV references (committed)
python extract_nastran_refs.py
```

The `.op2` / `.f06` / `.log` Nastran outputs land in `_nastran_outputs/` (gitignored).
Only the CSV reference files and BDFs are committed.
