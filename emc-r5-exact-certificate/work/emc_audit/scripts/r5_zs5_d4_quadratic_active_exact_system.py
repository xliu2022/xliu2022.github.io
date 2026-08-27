#!/usr/bin/env python3
"""Reconstruct the frozen quadratic candidate's active system over Q[z]/(cusp).

This stage does not yet solve the system.  It independently rebuilds every
near-active degree-four Bernstein row using Fraction, collapses the exact S3
stabilizer equalities by union-find, reduces the RHS modulo the cusp quartic,
and deduplicates identical equations.  The output is a compact audit manifest
for the subsequent exact nullspace/rational repair.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import sys
from collections import Counter, defaultdict
from fractions import Fraction as F
from pathlib import Path


M = F(1, 5)
CUSP = (F(-1024), F(9225), F(15250), F(11250), F(3125))
BASIS_EXPONENTS = (
    (0, 0, 0, 0),
    (1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1),
    (2, 0, 0, 0), (1, 1, 0, 0), (1, 0, 1, 0), (1, 0, 0, 1),
    (0, 2, 0, 0), (0, 1, 1, 0), (0, 1, 0, 1),
    (0, 0, 2, 0), (0, 0, 1, 1), (0, 0, 0, 2),
)
BASIS = len(BASIS_EXPONENTS)
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(HERE))
import r5_zs5_d4_projected_quadratic_cutting_probe as common

GEOMETRY = common.MANIFEST
ACTIVE = (ROOT / "work/emc_audit/computer/"
          "r5_zs5_d4_projected_quadratic_active_audit.json")
OUTPUT = (ROOT / "work/emc_audit/computer/"
          "r5_zs5_d4_projected_quadratic_active_exact_system.json")


class UnionFind:
    def __init__(self, size):
        self.parent = list(range(size))

    def find(self, value):
        while self.parent[value] != value:
            self.parent[value] = self.parent[self.parent[value]]
            value = self.parent[value]
        return value

    def union(self, left, right):
        left, right = self.find(left), self.find(right)
        if left != right:
            if right < left:
                left, right = right, left
            self.parent[right] = left


def symmetry_representatives(projected_keys):
    union = UnionFind(BASIS * len(projected_keys))
    exponent_where = {q: i for i, q in enumerate(BASIS_EXPONENTS)}
    for key_index, key in enumerate(projected_keys):
        for permutation in common.stabilizers(key):
            for basis, exponent in enumerate(BASIS_EXPONENTS):
                transformed = [0, 0, 0, exponent[3]]
                for old in range(3):
                    transformed[permutation[old]] = exponent[old]
                other = exponent_where[tuple(transformed)]
                union.union(BASIS * key_index + basis,
                            BASIS * key_index + other)
    roots = sorted({union.find(q) for q in range(len(union.parent))})
    where = {root: index for index, root in enumerate(roots)}
    return tuple(where[union.find(q)] for q in range(len(union.parent))), roots


def decode_key(key):
    return int(key[0]), tuple(int(q) for q in key[1:5])


def global_point(local, box, b_box):
    return tuple((F(box[j]) + local[j]) / 5 for j in range(4)) + (
        (F(b_box) + local[4]) / 5,)


def add_poly(target, source, scale=F(1), shift=0):
    for degree, coefficient in enumerate(source):
        target[degree + shift] += scale * coefficient


def reduce_cusp(poly):
    poly = list(poly)
    while len(poly) > 4:
        degree = len(poly) - 1
        lead = poly.pop()
        if not lead:
            continue
        shift = degree - 4
        # z^4 = -(c0+c1*z+c2*z^2+c3*z^3)/c4.
        for j in range(4):
            poly[shift + j] -= lead * CUSP[j] / CUSP[4]
    poly.extend([F(0)] * (4 - len(poly)))
    return tuple(poly)


def kernel_polynomials(point):
    sum_zero = [F(0)] * 5
    sum_one = [F(0)] * 5
    for state in itertools.product((0, 1), repeat=4):
        subtotal = sum(state[j] * point[j] for j in range(4))
        degree = 4 - sum(state)
        if subtotal > 1:
            sum_zero[degree] += 1
        if subtotal + point[4] > 1:
            sum_one[degree] += 1
    return sum_zero, sum_one


def basis_g(basis, h, r, s):
    if basis == 0:
        return F(1)
    if basis <= 4:
        j = basis - 1
        return (h[j][r] + h[j][s]) / 2
    at = 5
    for j in range(4):
        for k in range(j, 4):
            if at == basis:
                if j == k:
                    return h[j][r] * h[j][s]
                return (h[j][r] * h[k][s] + h[j][s] * h[k][r]) / 2
            at += 1
    raise AssertionError("bad basis")


def product_blossom(b, u, h, basis, slots):
    answer = F(0)
    for a in range(4):
        for c in range(4):
            if a == c:
                continue
            remaining = [j for j in range(4) if j not in (a, c)]
            answer += (b[slots[a]] * u[slots[c]]
                       * basis_g(basis, h, slots[remaining[0]], slots[remaining[1]]))
    return answer / 12


def reconstruct(cell, simplex, alpha, vertices, key_where, representative):
    b_box, box = decode_key(cell["key"])
    cell_points = [global_point(vertices[q], box, b_box) for q in cell["vertices"]]
    center = tuple(sum(point[j] for point in cell_points) / len(cell_points)
                   for j in range(5))
    sum_zero, sum_one = kernel_polynomials(center)
    points = [global_point(vertices[q], box, b_box) for q in simplex]
    slots = tuple(q for q, multiplicity in enumerate(alpha)
                  for _ in range(multiplicity))
    assert len(slots) == 4
    b = tuple(point[4] for point in points)
    b_average = sum(b[q] for q in slots) / 4
    rhs = [F(0)] * 7
    rhs[0] += b_average
    add_poly(rhs, sum_zero, M - b_average)
    add_poly(rhs, sum_zero, -b_average, shift=1)
    add_poly(rhs, sum_one, -M)
    for axis in range(4):
        if all(points[q][axis] == 1 for q in slots):
            rhs[3] += M - b_average
            rhs[4] -= b_average
    if all(points[q][4] == 1 for q in slots):
        rhs[4] -= M

    row = {}
    for omitted in range(4):
        projection = cell["projections"][omitted]
        axes = projection["axes"]
        key_index = key_where[projection["key"]]
        u = tuple(point[omitted] - M for point in points)
        h = (*(tuple(point[axis] - M for point in points) for axis in axes),
             tuple(point[4] - M for point in points))
        for basis in range(BASIS):
            value = product_blossom(b, u, h, basis, slots)
            column = representative[BASIS * key_index + basis]
            row[column] = row.get(column, F(0)) + value
    row = tuple(sorted((column, value) for column, value in row.items() if value))
    return row, reduce_cusp(rhs), b_box


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--print-json", action="store_true")
    args = parser.parse_args()
    geometry = json.loads(GEOMETRY.read_text())
    active = json.loads(ACTIVE.read_text())
    target = {(int(item["simplex"]), tuple(item["alpha"]))
              for item in active["active_rows"]}
    assert len(target) == active["active_count"]
    target_by_sid = defaultdict(list)
    for sid, alpha in target:
        target_by_sid[sid].append(alpha)
    vertices = [tuple(F(q) for q in point) for point in geometry["vertices"]]
    key_where = {key: index for index, key in enumerate(geometry["projected_keys"])}
    representative, roots = symmetry_representatives(geometry["projected_keys"])
    equations = Counter()
    by_box = Counter()
    sid = 0
    for cell in geometry["cells"]:
        for simplex in cell["simplices"]:
            for alpha in target_by_sid.get(sid, ()):
                row, rhs, b_box = reconstruct(cell, simplex, alpha, vertices,
                                              key_where, representative)
                equations[(row, rhs, b_box)] += 1
                by_box[b_box] += 1
            sid += 1
    assert sum(equations.values()) == len(target)
    zero_row_bad = [(rhs, count) for (row, rhs, _), count in equations.items()
                    if not row and any(rhs)]
    assert not zero_row_bad
    unique_by_box = Counter(b_box for (_, _, b_box) in equations)
    payload = {
        "schema": "r5_zs5_d4_projected_quadratic_active_exact_system_v1",
        "proof_status": "EXACT ACTIVE-SYSTEM RECONSTRUCTION; NOT YET A PROOF",
        "geometry_sha256": hashlib.sha256(GEOMETRY.read_bytes()).hexdigest(),
        "active_source_sha256": hashlib.sha256(ACTIVE.read_bytes()).hexdigest(),
        "original_variable_count": BASIS * len(geometry["projected_keys"]),
        "stabilizer_collapsed_variable_count": len(roots),
        "active_row_count": len(target),
        "unique_exact_equation_count": len(equations),
        "active_rows_by_b_box": dict(sorted(by_box.items())),
        "unique_equations_by_b_box": dict(sorted(unique_by_box.items())),
        "duplicate_multiplicity_distribution": dict(sorted(Counter(equations.values()).items())),
        "zero_row_bad_count": len(zero_row_bad),
        "original_to_collapsed_variable": list(representative),
        "collapsed_roots": roots,
        "equations": [{
            "b_box": b_box,
            "row": [[column, str(value)] for column, value in row],
            "rhs_mod_cusp": [str(q) for q in rhs],
            "multiplicity": multiplicity,
        } for (row, rhs, b_box), multiplicity in sorted(
            equations.items(), key=lambda item: (item[0][2], item[0][0], item[0][1]))],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    if args.print_json:
        print(json.dumps(payload, indent=2))
    print("active exact reconstruction rows", len(target),
          "unique equations", len(equations))
    print("wrote", args.output)


if __name__ == "__main__":
    main()
