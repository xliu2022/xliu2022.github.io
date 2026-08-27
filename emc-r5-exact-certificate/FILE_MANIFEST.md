# Published artifact manifest

This file separates formal verification inputs from generated outputs,
discovery-only programs, and independent audit material.

## Formal verification programs

| File | Role |
| --- | --- |
| `r5_zs5_d4_orbit_geometry_exact.py` | Reconstructs the rational cell geometry and pulling triangulations and compares them with the frozen geometry manifest. |
| `r5_zs5_d4_quadratic_active_exact_system.py` | Reconstructs the 1,620 exact equations representing all 12,794 candidate active rows. |
| `r5_zs5_d4_quadratic_face_repair_exact.py` | Reconstructs the coefficient table using exact quotient-ring and FLINT calculations. |
| `r5_zs5_d4_quadratic_exact_flat.py` | Regenerates the outward-rounded dyadic interval certificate. |
| `r5_zs5_d4_projected_quadratic_exact_verify.cpp` | Scans the complete Bernstein row stream with checked integer and interval arithmetic. |
| `r5_zs5_d4_quadratic_exact_manifest.py` | Runs the two exact stages as one verification command and writes the final manifest and log. |

All paths in this table are relative to `work/emc_audit/scripts/`.

## Frozen certificate inputs and outputs

| File | Role |
| --- | --- |
| `r5_zs5_d4_relevant_cell_orbits_exact_flat.txt` | Exact upstream records for the relevant cell orbits. |
| `r5_zs5_d4_orbit_geometry_exact.json` | Frozen rational geometry and pulling-triangulation manifest. |
| `r5_zs5_d4_quadratic_geometry_exact_flat.txt` | Flattened exact geometry used by the row verifier. |
| `r5_zs5_d4_projected_quadratic_symmetric_zero_highs_probe.json` | Floating candidate retained only as input to exact reconstruction. |
| `r5_zs5_d4_projected_quadratic_active_audit.json` | Candidate active-row identifiers and multiplicities. |
| `r5_zs5_d4_projected_quadratic_active_exact_system.json` | Frozen exact active system. |
| `r5_zs5_d4_projected_quadratic_face_repair_exact.json` | Frozen exactly reconstructed coefficient table. |
| `r5_zs5_d4_projected_quadratic_exact_flat.txt` | Frozen outward-rounded interval table. |
| `r5_zs5_d4_projected_quadratic_exact_manifest.json` | Final machine-readable acceptance manifest. |
| `r5_zs5_d4_projected_quadratic_exact_verify.log` | Frozen output of the exhaustive row verifier. |

All paths in this table are relative to `work/emc_audit/computer/`.

## Discovery and provenance programs

| File | Role |
| --- | --- |
| `r5_zs5_d4_local_arrangement_counts.cpp` | Enumerates the upstream arrangement and orbit data. |
| `r5_zs5_d4_quadratic_geometry_flat.py` | Produces the flattened geometry representation. |
| `r5_zs5_d4_projected_quadratic_cutting_probe.py` | Searches numerically for a candidate projected quadratic table. |
| `r5_zs5_d4_projected_quadratic_scan.cpp` | Performs the exploratory row scan used during candidate discovery. |
| `r5_zs5_d4_quadratic_active_audit.py` | Audits the proposed active-row list before exact reconstruction. |

These programs document how the candidate was found. The proof does not
accept a row because of a floating residual or discovery-stage classification.

## Independent red-team material

| File | Role |
| --- | --- |
| `work/emc_audit/scripts/r5_zs5_d4_quadratic_redteam_exact.py` | Independent exact checker. |
| `work/emc_audit/computer/r5_zs5_d4_quadratic_redteam_exact.json` | Frozen red-team result. |
| `work/emc_audit/computer/r5_zs5_d4_quadratic_redteam_SHA256SUMS.txt` | Checks the red-team report, checker, output, and exact manifest. |
| `work/emc_audit/agents/r5_zs5_d4_local_quadratic_redteam.md` | Human-readable independent audit report. |

