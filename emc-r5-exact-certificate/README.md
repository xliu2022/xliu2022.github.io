# Exact certificate for the five-uniform Erdős matching conjecture

This repository contains the computer-verifiable certificate used for the
projected quadratic majorant in Proposition 3.3 of the paper
*The Erdős matching conjecture for 5-uniform hypergraphs with large ground
sets*, by Jie Han, Xizhi Liu, and Jian Wang.

The repository contains no manuscript source. The analytic reductions,
the equality analysis, the transfer to hypergraph matchings, and the final
exactification remain part of the paper. This repository verifies the finite
computer-assisted step used in the zero-endpoint branch.

## What the certificate proves

Let $p_5$ be the positive solution of
$(1-p_5/5)^5=1-p_5^5$, and put $z=(1-p_5)/p_5$. The exact computation is
performed in the quartic quotient algebra determined by

$$
3125z^4+11250z^3+15250z^2+9225z-1024=0.
$$

The certificate checks all 5,164,740 degree-four Bernstein rows:

- 5,151,946 rows have a strictly positive outward-rounded lower endpoint;
- 12,794 rows are represented by 1,620 exact identities in the quotient
  basis $1,z,z^2,z^3$;
- the smallest strict lower bound is
  `5411340156129015928571 / 2^90`;
- the row-stream fingerprint is `49063aa32fa17a6f`;
- the closed top-face corrections are included.

Floating-point computations are used only to propose candidate active rows
and convenient pivots. Exact rational algebra, quotient-ring identities, and
outward-rounded interval bounds are the acceptance tests.

## Requirements

The archived verification was run with:

- Python 3.12.13;
- NumPy 2.5.2;
- SciPy 1.18.1;
- python-flint 0.9.0;
- an Apple Clang compiler supporting C++20.

The Python package versions are pinned in `requirements.txt`. The verifier
compiles its C++ scanner in a temporary directory; no compiled binary needs
to be committed.

## Quick verification

From the repository root, run:

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python work/emc_audit/scripts/r5_zs5_d4_orbit_geometry_exact.py --verify
.venv/bin/python work/emc_audit/scripts/r5_zs5_d4_quadratic_exact_manifest.py
```

The geometry replay must end with:

```text
CERTIFIED
```

The complete two-stage replay must end with:

```text
EXACT TWO-STAGE MANIFEST PASS
```

The second command first reconstructs and verifies all exact active
identities, then regenerates the dyadic interval table and scans the complete
Bernstein row stream.

## Independent red-team check

The independent certificate bundle can be checked from the repository root:

```bash
shasum -a 256 -c work/emc_audit/computer/r5_zs5_d4_quadratic_redteam_SHA256SUMS.txt
```

This checks the red-team report, checker, output, and final exact manifest.

## Repository layout

- `work/emc_audit/scripts/`: exact generators, exhaustive verifier, and
  discovery programs;
- `work/emc_audit/computer/`: frozen geometry, coefficient data, active
  system, interval certificate, manifests, and verification log;
- `work/emc_audit/agents/r5_zs5_d4_local_quadratic_redteam.md`: the
  independent local certificate audit;
- `FILE_MANIFEST.md`: the role of every published artifact;
- `SHA256SUMS.txt`: repository-level checksums, excluding the checksum file
  itself.

The relative path `work/emc_audit/` is intentional. The verification
programs resolve their frozen inputs from this layout.

## Proof-critical and discovery-only files

The proof-critical entry points are
`r5_zs5_d4_orbit_geometry_exact.py` and
`r5_zs5_d4_quadratic_exact_manifest.py`. The latter orchestrates the exact
active-system reconstruction, exact face repair, interval regeneration, and
the C++ exhaustive scan.

The cutting-plane, floating scan, and candidate-audit programs are retained
for provenance. Their numerical output does not by itself establish an
identity, sign, rank, or feasible majorant.

## Integrity

To verify the published tree, run:

```bash
shasum -a 256 -c SHA256SUMS.txt
```

The checksum manifest does not include itself. A tagged release should be
identified by both its Git commit and the checksum of its release archive.

## Citation and license

Citation metadata are provided in `CITATION.cff`. Repository and publication
identifiers should be added there after the GitHub release and any archival
DOI are created.

The public license is still awaiting agreement among the authors. Until
`LICENSE` is replaced, the repository is distributed without an
open-source license.
