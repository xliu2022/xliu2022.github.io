# Turán density verification

Exact arithmetic checks for *On possible Turán densities under general uniformity conditions* by Heng Li and Xizhi Liu.

## Running the checks

Requires Python 3.9 or later and only the Python standard library. From the repository root, run:

```bash
python3 -B verify_certificates.py --main /path/to/manuscript.tex
```

Supply the manuscript's LaTeX source separately. The vertex-pair verifier reads the 16 rows of Table 2, identified by the label `vp-tab:extension-data`, directly from that file. If the source is saved as `main.tex` in the repository root, `--main` can be omitted.

Add `--verbose` to show the individual calculations. A failed check produces a nonzero exit status. Normal verification reads the data without changing the stored certificates.

## Contents

- `verify_certificates.py`: runs all verification groups.
- `verify_degree_intervals.py`: checks the minimum-degree interval calculations.
- `parts/ordinary/verify_all.py`: runs the rounding, 37-leaf star, and book checks.
- `parts/ordinary/certificates/`: exact polynomial arithmetic, the three ordinary-density verifiers, and the two JSON files containing the 132 book test points and their certificates.
- `parts/vertex_pair/scripts/verify_extension.py`: checks the interval-extension table and its explicit palette relations.

All arithmetic checks use exact rational values. They verify the finite calculations accompanying the analytic and combinatorial proofs in the manuscript.
