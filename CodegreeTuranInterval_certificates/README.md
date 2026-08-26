# Exact certificates

The three data files contain the exact finite routing tables used in the manuscript:

- `weighted_path_chain_certificate.txt`: 113 homogeneous rows;
- `heterogeneous_satellite_chain_certificate.txt`: 57 heterogeneous rows;
- `pi2_4_high_satellite_rows.txt`: 42 high-density rows.

## Main verifier

Run from the manuscript directory:

```bash
python3 CodegreeTuranInterval_certificates/verify_pi2_4_exact.py
```

The verifier uses only the Python standard library and exact `Fraction` arithmetic. It checks the 113 homogeneous rows, 57 heterogeneous rows, and 42 high-density rows. For the `D25` bridge it verifies the support-coverage identity, `det(C)=4`, the inverse and `K`-matrix patterns, the three principal blocks, and the endpoint margins. Use `--show-witnesses` to print the labeled graph witness selected for every heterogeneous row.

Expected final line:

```text
ALL PI_2(4) EXACT CERTIFICATE CHECKS PASSED
```

## Optional cross-check

An independent cross-check of the high-density module requires SymPy and is not part of the main verification dependency:

```bash
python3 CodegreeTuranInterval_certificates/pi2_4_high_module_certificate.py
```

Its expected final line is:

```text
HIGH MODULE EXACT CERTIFICATES PASSED
```

## Higher-uniformity checks

The rational Bernstein coefficients, endpoint gaps, and safety values in the
$k=5$ and $k=6$ compilers are checked independently by
`verify_higher_uniformities_exact.py`.  Run

```bash
python3 CodegreeTuranInterval_certificates/verify_higher_uniformities_exact.py
```

from the manuscript directory.  Its expected final line is:

```text
ALL HIGHER-UNIFORMITY EXACT CERTIFICATES PASSED
```

## File integrity

The SHA-256 digests of the three data files and the three verification scripts are recorded in `SHA256SUMS`. From the certificate directory, verify them with

```bash
shasum -a 256 -c SHA256SUMS
```

or, on a system with GNU coreutils, with

```bash
sha256sum -c SHA256SUMS
```
