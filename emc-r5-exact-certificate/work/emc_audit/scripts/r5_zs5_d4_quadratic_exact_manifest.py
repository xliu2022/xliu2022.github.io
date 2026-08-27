#!/usr/bin/env python3
"""One-command exact audit for the d=4 projected quadratic ZS5 dual.

This wrapper deliberately makes the two logically different verification
stages inseparable:

1. replay all 1,620 distinct active-face equations in Q[z]/(P), whose
   multiplicities sum to all 12,794 active Bernstein rows;
2. regenerate the outward dyadic enclosure byte-for-byte, compile the exact
   C++ scanner, and require strict positivity on every one of the remaining
   5,151,946 rows among the full 5,164,740-row stream.

The C++ stage may classify a row as active only if stage 1 has already passed.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
import sys
import tempfile
from fractions import Fraction as F
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
COMPUTER = ROOT / "work/emc_audit/computer"

FILES = {
    "geometry_json": COMPUTER / "r5_zs5_d4_orbit_geometry_exact.json",
    "geometry_flat": COMPUTER / "r5_zs5_d4_quadratic_geometry_exact_flat.txt",
    "floating_candidate": COMPUTER / "r5_zs5_d4_projected_quadratic_symmetric_zero_highs_probe.json",
    "active_audit": COMPUTER / "r5_zs5_d4_projected_quadratic_active_audit.json",
    "active_system": COMPUTER / "r5_zs5_d4_projected_quadratic_active_exact_system.json",
    "face_repair": COMPUTER / "r5_zs5_d4_projected_quadratic_face_repair_exact.json",
    "interval_certificate": COMPUTER / "r5_zs5_d4_projected_quadratic_exact_flat.txt",
    "active_system_generator": HERE / "r5_zs5_d4_quadratic_active_exact_system.py",
    "face_repair_generator": HERE / "r5_zs5_d4_quadratic_face_repair_exact.py",
    "interval_generator": HERE / "r5_zs5_d4_quadratic_exact_flat.py",
    "cpp_verifier": HERE / "r5_zs5_d4_projected_quadratic_exact_verify.cpp",
    "manifest_wrapper": HERE / "r5_zs5_d4_quadratic_exact_manifest.py",
}

EXPECTED_SHA256 = {
    "geometry_json": "681ab4c16db61bdad98d121c6c265c18918169e9b0081f5d47f4591eebf361c2",
    "geometry_flat": "567faf32d329ebb1a3a653edcd4918136d1053a0e93ce4ab6f2faa86d8307fd4",
    "floating_candidate": "322348e1f8451a6e7757c0ca8ba184230cd0d8e5d98652e7b977f5d941571897",
    "active_audit": "d08c9b2dccb89c0ae3a1d10abbc48b8028a4ecf9ee6b433e72419e211a18c7c0",
    "active_system": "aaf2b6fdf682386c93e2cdaca0c958071d8e9d253c2b40d112ad3a7ccd102a2e",
    "face_repair": "4d537d614540337e4b55226fb7edfe91109fd38a00b6bc0e6b817a37b785b5f1",
    "interval_certificate": "4dc19e34b4e13ce3a849ac286e8b7c8dc4c2ef4bcf0c60db4e047bc7be565239",
}

DEFAULT_OUTPUT = COMPUTER / "r5_zs5_d4_projected_quadratic_exact_manifest.json"
DEFAULT_LOG = COMPUTER / "r5_zs5_d4_projected_quadratic_exact_verify.log"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def exact_active_replay(system: dict, repair: dict) -> dict:
    coefficients = [[F(value) for value in polynomial]
                    for polynomial in repair["collapsed_coefficients_qz_mod_cusp"]]
    assert len(coefficients) == system["stabilizer_collapsed_variable_count"] == 6759
    equations = system["equations"]
    assert len(equations) == system["unique_exact_equation_count"] == 1620
    multiplicity = 0
    maximum_support = 0
    for equation_index, equation in enumerate(equations):
        residual = [F(value) for value in equation["rhs_mod_cusp"]]
        support = equation["row"]
        maximum_support = max(maximum_support, len(support))
        for column, value in support:
            column, value = int(column), F(value)
            for degree in range(4):
                residual[degree] += value * coefficients[column][degree]
        if any(residual):
            raise AssertionError(("nonzero active residual", equation_index, residual))
        multiplicity += int(equation["multiplicity"])
    assert multiplicity == system["active_row_count"] == 12794
    assert repair["exact_active_failure_count"] == 0
    return {
        "unique_equations": len(equations),
        "represented_active_rows": multiplicity,
        "maximum_sparse_support": maximum_support,
        "status": "EXACT_QZ_MOD_CUSP_PASS",
    }


def parse_scanner(stdout: str) -> dict:
    required = {
        "rows": r"rows (\d+)",
        "strict_rows": r"strict_rows (\d+)",
        "active_rows": r"exact_active_rows (\d+)",
        "minimum_strict_lower_scaled": r"minimum_strict_lower_scaled (\d+)",
        "scale_bits": r"scale_bits (\d+)",
        "largest_active_interval_width_scaled": r"largest_active_interval_width_scaled (\d+)",
        "row_stream_fnv1a64": r"row_stream_fnv1a64 ([0-9a-f]{16})",
    }
    parsed = {}
    for key, pattern in required.items():
        match = re.search(pattern, stdout)
        if not match:
            raise AssertionError(("missing scanner field", key))
        parsed[key] = match.group(1)
    for key in required:
        if key != "row_stream_fnv1a64":
            parsed[key] = int(parsed[key])
    assert "R5_ZS5_D4_QUADRATIC_EXACT_INTERVAL_PASS" in stdout
    assert "closed_top_face_corrections INCLUDED" in stdout
    assert parsed["rows"] == 5_164_740
    assert parsed["strict_rows"] == 5_151_946
    assert parsed["active_rows"] == 12_794
    assert parsed["minimum_strict_lower_scaled"] > 0
    assert parsed["scale_bits"] == 90
    return parsed


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cxx", default="clang++")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--log", type=Path, default=DEFAULT_LOG)
    args = parser.parse_args()

    hashes = {name: sha256(path) for name, path in FILES.items()}
    for name, expected in EXPECTED_SHA256.items():
        if hashes[name] != expected:
            raise AssertionError(("frozen artifact hash mismatch", name,
                                  hashes[name], expected))

    active = json.loads(FILES["active_audit"].read_text())
    assert active["active_count"] == 12794

    with tempfile.TemporaryDirectory(prefix="r5_zs5_d4_exact_") as temp_name:
        temp = Path(temp_name)
        regenerated_system = temp / "active_system.json"
        subprocess.run([
            sys.executable, str(FILES["active_system_generator"]),
            "--output", str(regenerated_system),
        ], cwd=ROOT, check=True)
        if regenerated_system.read_bytes() != FILES["active_system"].read_bytes():
            raise AssertionError("active system reconstruction differs byte-for-byte")

        regenerated_repair = temp / "face_repair.json"
        subprocess.run([
            sys.executable, str(FILES["face_repair_generator"]),
            "--system", str(regenerated_system),
            "--floating", str(FILES["floating_candidate"]),
            "--output", str(regenerated_repair),
        ], cwd=ROOT, check=True)
        if regenerated_repair.read_bytes() != FILES["face_repair"].read_bytes():
            raise AssertionError("exact face repair differs byte-for-byte")

        system = json.loads(regenerated_system.read_text())
        repair = json.loads(regenerated_repair.read_text())
        assert repair["system_sha256"] == hashes["active_system"]
        assert repair["floating_source_sha256"] == hashes["floating_candidate"]
        assert system["geometry_sha256"] == hashes["geometry_json"]
        assert system["active_source_sha256"] == hashes["active_audit"]

        # This exact replay is intentionally completed before the interval
        # scanner is compiled, so no active row can be silently waived.  The
        # equations have just been reconstructed from all 12,794 row metadata
        # items and the exact geometry, rather than trusted as an intermediate.
        active_result = exact_active_replay(system, repair)
        print("active exact replay PASS", active_result, flush=True)

        regenerated = temp / "interval_flat.txt"
        subprocess.run([
            sys.executable, str(FILES["interval_generator"]),
            "--output", str(regenerated),
        ], cwd=ROOT, check=True)
        if regenerated.read_bytes() != FILES["interval_certificate"].read_bytes():
            raise AssertionError("interval certificate regeneration differs byte-for-byte")

        executable = temp / "exact_verify"
        compile_command = [
            args.cxx, "-std=c++20", "-O3",
            str(FILES["cpp_verifier"]), "-o", str(executable),
        ]
        subprocess.run(compile_command, cwd=ROOT, check=True)
        verify_command = [
            str(executable), str(FILES["geometry_flat"]),
            str(FILES["interval_certificate"]),
        ]
        completed = subprocess.run(verify_command, cwd=ROOT, check=True,
                                   text=True, capture_output=True)
        sys.stderr.write(completed.stderr)
        print(completed.stdout, end="", flush=True)

    scanner_result = parse_scanner(completed.stdout)
    args.log.parent.mkdir(parents=True, exist_ok=True)
    args.log.write_text(completed.stdout)
    payload = {
        "schema": "r5_zs5_d4_projected_quadratic_exact_manifest_v1",
        "proof_status": "EXACT D4 PROJECTED-DUAL REPLAY PASS",
        "logical_stages": {
            "active_face": active_result,
            "all_rows_interval_replay": scanner_result,
            "closed_top_face_corrections": "INCLUDED",
            "interval_certificate_regeneration": "BYTE_IDENTICAL_PASS",
        },
        "hashes_sha256": hashes,
        "expected_frozen_hashes_sha256": EXPECTED_SHA256,
        "scanner_stdout_sha256": hashlib.sha256(completed.stdout.encode()).hexdigest(),
        "reproduction_command": (
            ".venv/bin/python work/emc_audit/scripts/"
            "r5_zs5_d4_quadratic_exact_manifest.py"
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print("EXACT TWO-STAGE MANIFEST PASS")
    print("wrote", args.output)
    print("manifest sha256", sha256(args.output))


if __name__ == "__main__":
    main()
