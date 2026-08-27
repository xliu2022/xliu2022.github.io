#!/usr/bin/env python3
"""Floating cutting-plane search for the d=4 projected-cell quadratic dual.

The C++ separator streams all 5,164,740 degree-four Bernstein rows on the
exact five-dimensional arrangement.  This driver keeps only the worst sparse
rows, solves a bounded min-max LP, and repeats.  A pass is still diagnostic;
the candidate must later be rationalized and replayed with Fraction/Q(z).
"""
from __future__ import annotations

import argparse
import itertools
import json
import subprocess
from pathlib import Path

import numpy as np
from scipy.optimize import linprog
from scipy.sparse import coo_matrix, csr_matrix, eye, hstack, vstack


BASIS = 15
BASIS_EXPONENTS = (
    (0, 0, 0, 0),
    (1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1),
    (2, 0, 0, 0), (1, 1, 0, 0), (1, 0, 1, 0), (1, 0, 0, 1),
    (0, 2, 0, 0), (0, 1, 1, 0), (0, 1, 0, 1),
    (0, 0, 2, 0), (0, 0, 1, 1), (0, 0, 0, 2),
)
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
CPP = HERE / "r5_zs5_d4_projected_quadratic_scan.cpp"
GEOMETRY = (ROOT / "work/emc_audit/computer/"
            "r5_zs5_d4_quadratic_geometry_exact_flat.txt")
MANIFEST = (ROOT / "work/emc_audit/computer/"
            "r5_zs5_d4_orbit_geometry_exact.json")


def write_candidate(path, coefficients):
    path.write_text("\n".join(f"{q:.19g}" for q in coefficients) + "\n")


def scan(executable, coefficient_path, top, elevation=4):
    command = [str(executable), str(GEOMETRY), str(coefficient_path), str(top),
               str(elevation)]
    result = subprocess.run(command, text=True, check=True,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if result.stderr:
        print(result.stderr, end="", flush=True)
    summary = None
    parsed = []
    for line in result.stdout.splitlines():
        tokens = line.split()
        if not tokens:
            continue
        if tokens[0] == "SUMMARY":
            summary = {
                "row_count": int(tokens[2]),
                "minimum": float(tokens[4]),
                "returned": int(tokens[6]),
            }
        if tokens[0] != "ROW":
            continue
        raw_value = float(tokens[1])
        base = float(tokens[2])
        simplex = int(tokens[3])
        assert tokens[4] == "A"
        alpha = tuple(int(q) for q in tokens[5:11])
        cursor = 11
        coefficients = {}
        for _ in range(4):
            assert tokens[cursor] == "P"
            key = int(tokens[cursor + 1])
            values = [float(q) for q in tokens[cursor + 2:cursor + 2 + BASIS]]
            cursor += 2 + BASIS
            for basis, value in enumerate(values):
                column = BASIS * key + basis
                coefficients[column] = coefficients.get(column, 0.0) + value
        assert cursor == len(tokens)
        parsed.append((raw_value, base, coefficients, (simplex, alpha)))
    assert summary is not None and len(parsed) == summary["returned"]
    return summary, parsed


def permute_mask(mask, permutation):
    result = mask & 8
    for old in range(3):
        if mask & (1 << old):
            result |= 1 << permutation[old]
    return result


def stabilizers(key):
    b_box = key[0]
    boxes = tuple(key[1:4])
    signs = key[4:]
    answer = []
    for permutation in itertools.permutations(range(3)):
        permuted_boxes = [None] * 3
        for old in range(3):
            permuted_boxes[permutation[old]] = boxes[old]
        permuted_signs = [None] * 15
        for mask in range(1, 16):
            permuted_signs[permute_mask(mask, permutation) - 1] = signs[mask - 1]
        descriptor = b_box + "".join(permuted_boxes) + "".join(permuted_signs)
        if descriptor == key:
            answer.append(permutation)
    assert answer
    return tuple(answer)


def symmetry_equalities(projected_keys, variable_count):
    exponent_where = {exponent: index for index, exponent in enumerate(BASIS_EXPONENTS)}
    pairs = set()
    for key_index, key in enumerate(projected_keys):
        for permutation in stabilizers(key):
            for basis, exponent in enumerate(BASIS_EXPONENTS):
                transformed = [0, 0, 0, exponent[3]]
                for old in range(3):
                    transformed[permutation[old]] = exponent[old]
                other = exponent_where[tuple(transformed)]
                left = BASIS * key_index + basis
                right = BASIS * key_index + other
                if left != right:
                    pairs.add(tuple(sorted((left, right))))
    rows, columns, values = [], [], []
    for row, (left, right) in enumerate(sorted(pairs)):
        rows.extend((row, row)); columns.extend((left, right)); values.extend((1.0, -1.0))
    return coo_matrix((values, (rows, columns)),
                      shape=(len(pairs), variable_count)).tocsr()


def solve(rows, variable_count, projected_keys, target, bound):
    rr, cc, vv, rhs = [], [], [], []
    for coefficients, base in rows.values():
        row = len(rhs)
        nonstructural = abs(base) > 1e-13 or any(abs(q) > 1e-13
                                                for q in coefficients.values())
        rhs.append(base - (target if nonstructural else 0.0))
        # base + coefficients*x >= target becomes -coefficients*x <= base-target.
        for column, value in coefficients.items():
            if abs(value) > 1e-15:
                rr.append(row); cc.append(column); vv.append(-value)
    matrix = coo_matrix((vv, (rr, cc)), shape=(len(rhs), variable_count)).tocsr()
    identity = eye(variable_count, format="csr")
    t_column = csr_matrix(-np.ones((variable_count, 1)))
    bound_block = vstack((hstack((identity, t_column)),
                          hstack((-identity, t_column))), format="csr")
    augmented = vstack((hstack((matrix, csr_matrix((len(rhs), 1)))),
                        bound_block), format="csr")
    objective = np.zeros(variable_count + 1)
    objective[-1] = 1
    equality = symmetry_equalities(projected_keys, variable_count)
    augmented_equality = hstack((equality, csr_matrix((equality.shape[0], 1))),
                                format="csr")
    result = linprog(objective, A_ub=augmented,
                     b_ub=np.concatenate((rhs, np.zeros(2 * variable_count))),
                     A_eq=augmented_equality, b_eq=np.zeros(equality.shape[0]),
                     bounds=[(-bound, bound)] * variable_count + [(0, bound)],
                     method="highs")
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--executable", type=Path,
                        default=Path("/private/tmp/r5_zs5_d4_projected_quadratic_scan"))
    parser.add_argument("--candidate", type=Path,
                        default=Path("/private/tmp/r5_zs5_d4_projected_quadratic_candidate.txt"))
    parser.add_argument("--top", type=int, default=20_000)
    parser.add_argument("--iterations", type=int, default=50)
    parser.add_argument("--elevation", type=int, default=4)
    parser.add_argument("--target-slack", type=float, default=0.001)
    parser.add_argument("--bound", type=float, default=100)
    parser.add_argument("--json-out", type=Path)
    args = parser.parse_args()
    subprocess.run(["c++", "-O3", "-std=c++17", str(CPP), "-o",
                    str(args.executable)], check=True)
    geometry = json.loads(MANIFEST.read_text())
    projected_keys = geometry["projected_keys"]
    variable_count = BASIS * len(projected_keys)
    coefficients = np.zeros(variable_count)
    rows = {}
    passed = False
    adjusted_minimum = None
    summary = None

    for iteration in range(args.iterations):
        write_candidate(args.candidate, coefficients)
        summary, cuts = scan(args.executable, args.candidate, args.top,
                             args.elevation)
        adjusted = []
        added = 0
        for raw_value, base, row, meta in cuts:
            nonstructural = abs(base) > 1e-13 or any(abs(q) > 1e-13
                                                    for q in row.values())
            value = raw_value - (args.target_slack if nonstructural else 0.0)
            adjusted.append(value)
            if value < -1e-10 and meta not in rows:
                rows[meta] = (row, base)
                added += 1
        adjusted_minimum = min(adjusted)
        print("iteration", iteration, "LP rows", len(rows), "added", added,
              "raw minimum", summary["minimum"], "adjusted minimum",
              adjusted_minimum, flush=True)
        if adjusted_minimum >= -2e-10:
            passed = True
            break
        if not added:
            raise RuntimeError("negative separator but no new sparse rows")
        result = solve(rows, variable_count, projected_keys,
                       args.target_slack, args.bound)
        print("LP status", result.message, flush=True)
        if not result.success:
            break
        coefficients = result.x[:-1]
        print("T", result.x[-1], flush=True)

    if args.json_out:
        payload = {
            "proof_status": ("FLOATING FULL-CONTINUOUS QUADRATIC CANDIDATE ONLY"
                             if passed else
                             "FLOATING QUADRATIC CUTTING DIAGNOSTIC ONLY"),
            "full_continuous_scan_passed": passed,
            "basis": ["1", "h1", "h2", "h3", "hb",
                      "h1^2", "h1h2", "h1h3", "h1hb", "h2^2",
                      "h2h3", "h2hb", "h3^2", "h3hb", "hb^2"],
            "projected_key_count": len(projected_keys),
            "variable_count": variable_count,
            "LP_row_count": len(rows),
            "target_slack": args.target_slack,
            "bernstein_elevation_degree": args.elevation,
            "last_raw_minimum": summary["minimum"] if summary else None,
            "last_adjusted_minimum": adjusted_minimum,
            "closed_top_corrections_included": True,
            "g": [{"key": key,
                   "coefficients": [float(coefficients[BASIS * index + j])
                                    for j in range(BASIS)]}
                  for index, key in enumerate(projected_keys)],
        }
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(json.dumps(payload, indent=2) + "\n")
        print("wrote", args.json_out)


if __name__ == "__main__":
    main()
