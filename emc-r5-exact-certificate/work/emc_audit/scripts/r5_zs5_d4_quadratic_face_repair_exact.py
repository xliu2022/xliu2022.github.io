#!/usr/bin/env python3
"""Exact Q[z]/(cusp) repair on the d=4 quadratic candidate's zero face.

The 12,794 numerical zero rows collapse to 1,620 exact equations.  After S3
stabilizer identification the system splits by the four b-boxes and is highly
underdetermined.  We round free coordinates to a common rational denominator,
choose a numerically well-conditioned exact pivot minor in each block, and
solve its four cusp-basis right-hand sides with FLINT fmpq matrices.  Every
remaining active equation is then replayed exactly.

Passing this script proves only exact membership in the active face.  The full
5,164,740-row positivity replay is a separate required stage.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from collections import defaultdict
from fractions import Fraction as F
from pathlib import Path

import numpy as np
import scipy.linalg as la
from flint import fmpq, fmpq_mat


Z_FLOAT = 0.09500753725790037273535596965
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SYSTEM = (ROOT / "work/emc_audit/computer/"
          "r5_zs5_d4_projected_quadratic_active_exact_system.json")
FLOATING = (ROOT / "work/emc_audit/computer/"
            "r5_zs5_d4_projected_quadratic_symmetric_zero_highs_probe.json")
OUTPUT = (ROOT / "work/emc_audit/computer/"
          "r5_zs5_d4_projected_quadratic_face_repair_exact.json")


def as_fmpq(value):
    value = F(value)
    return fmpq(value.numerator, value.denominator)


def as_fraction(value):
    return F(str(value))


def eval_poly(poly, z=Z_FLOAT):
    answer = 0.0
    for coefficient in reversed(poly):
        answer = answer * z + float(coefficient)
    return answer


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--free-denominator", type=int, default=10**15)
    parser.add_argument("--system", type=Path, default=SYSTEM)
    parser.add_argument("--floating", type=Path, default=FLOATING)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    system = json.loads(args.system.read_text())
    floating = json.loads(args.floating.read_text())
    original = np.asarray([q for item in floating["g"] for q in item["coefficients"]])
    mapping = np.asarray(system["original_to_collapsed_variable"], dtype=int)
    collapsed_count = int(system["stabilizer_collapsed_variable_count"])
    collapsed_float = np.zeros(collapsed_count)
    multiplicity = np.zeros(collapsed_count)
    np.add.at(collapsed_float, mapping, original)
    np.add.at(multiplicity, mapping, 1)
    collapsed_float /= multiplicity
    coefficient = [[F(round(value * args.free_denominator), args.free_denominator),
                    F(0), F(0), F(0)] for value in collapsed_float]

    equations_by_box = defaultdict(list)
    for item in system["equations"]:
        row = {int(column): F(value) for column, value in item["row"]}
        rhs = tuple(F(q) for q in item["rhs_mod_cusp"])
        equations_by_box[int(item["b_box"])].append((row, rhs))

    block_summary = []
    pivot_columns_all = set()
    for b_box in range(1, 5):
        equations = equations_by_box[b_box]
        columns = sorted({column for row, _ in equations for column in row})
        where = {column: index for index, column in enumerate(columns)}
        matrix = np.zeros((len(equations), len(columns)))
        for r, (row, _) in enumerate(equations):
            for column, value in row.items():
                matrix[r, where[column]] = float(value)
        _, triangular, column_order = la.qr(matrix, mode="economic", pivoting=True)
        diagonal = np.abs(np.diag(triangular))
        rank = int(np.sum(diagonal > diagonal[0] * 1e-10))
        pivot_local = list(map(int, column_order[:rank]))
        pivot_columns = [columns[q] for q in pivot_local]
        pivot_columns_all.update(pivot_columns)
        pivot_set = set(pivot_columns)
        free_columns = [column for column in columns if column not in pivot_set]
        square_float = matrix[:, pivot_local]
        _, row_triangular, row_order = la.qr(square_float.T, mode="economic",
                                             pivoting=True)
        row_diagonal = np.abs(np.diag(row_triangular))
        assert np.sum(row_diagonal > row_diagonal[0] * 1e-10) == rank
        pivot_rows = list(map(int, row_order[:rank]))

        square_entries = []
        rhs_entries = []
        for row_index in pivot_rows:
            row, rhs = equations[row_index]
            square_entries.extend(as_fmpq(row.get(column, F(0)))
                                  for column in pivot_columns)
            free_constant = sum(row.get(column, F(0)) * coefficient[column][0]
                                for column in free_columns)
            rhs_entries.append(as_fmpq(-rhs[0] - free_constant))
            rhs_entries.extend(as_fmpq(-rhs[degree]) for degree in range(1, 4))
        square = fmpq_mat(rank, rank, square_entries)
        right = fmpq_mat(rank, 4, rhs_entries)
        solution = square.solve(right, algorithm="dixon")
        for local, column in enumerate(pivot_columns):
            coefficient[column] = [as_fraction(solution[local, degree])
                                   for degree in range(4)]

        condition = float(np.linalg.cond(square_float[pivot_rows]))
        block_summary.append({
            "b_box": b_box,
            "equation_count": len(equations),
            "touched_column_count": len(columns),
            "exact_rank": rank,
            "pivot_condition_estimate": condition,
            "pivot_column_count": len(pivot_columns),
            "free_column_count": len(free_columns),
        })
        print("box", b_box, "rank", rank, "condition", condition, flush=True)

    # Exact replay of every unique active equation in the quotient basis.
    failures = []
    for b_box, equations in equations_by_box.items():
        for index, (row, rhs) in enumerate(equations):
            residual = list(rhs)
            for column, value in row.items():
                for degree in range(4):
                    residual[degree] += value * coefficient[column][degree]
            if any(residual):
                failures.append((b_box, index, tuple(residual)))
                if len(failures) >= 10:
                    break
        if failures:
            break
    if failures:
        raise AssertionError(f"active exact replay failed: {failures[:2]}")

    repaired_float = np.asarray([eval_poly(poly) for poly in coefficient])
    deviation = repaired_float - collapsed_float
    expanded = [coefficient[q] for q in mapping]
    payload = {
        "schema": "r5_zs5_d4_projected_quadratic_face_repair_exact_v1",
        "proof_status": "EXACT ACTIVE-FACE REPAIR ONLY; FULL POSITIVITY NOT YET VERIFIED",
        "system_sha256": hashlib.sha256(args.system.read_bytes()).hexdigest(),
        "floating_source_sha256": hashlib.sha256(args.floating.read_bytes()).hexdigest(),
        "free_denominator": args.free_denominator,
        "block_summary": block_summary,
        "unique_active_equation_count": sum(len(q) for q in equations_by_box.values()),
        "exact_active_failure_count": len(failures),
        "maximum_collapsed_coefficient_deviation_at_float_cusp": float(np.max(np.abs(deviation))),
        "rms_collapsed_coefficient_deviation_at_float_cusp": float(np.sqrt(np.mean(deviation**2))),
        "collapsed_coefficients_qz_mod_cusp": [[str(q) for q in poly]
                                                for poly in coefficient],
        "expanded_coefficients_qz_mod_cusp": [[str(q) for q in poly]
                                               for poly in expanded],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print("exact active replay PASS", len(equations_by_box), "blocks")
    print("maximum deviation", payload["maximum_collapsed_coefficient_deviation_at_float_cusp"])
    print("wrote", args.output)


if __name__ == "__main__":
    main()
