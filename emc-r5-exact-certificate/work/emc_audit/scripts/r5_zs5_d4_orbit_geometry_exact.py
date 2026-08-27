#!/usr/bin/env python3
"""Exact geometry replay for S4 orbit representatives of d=4 relevant cells.

Consumes the exact-rational flat discovery output of
``r5_zs5_d4_local_arrangement_counts.cpp``.  For every canonical global-box
cell it reconstructs the H-polytope, filters its true vertices by exact active
rank, and builds the canonical pulling triangulation.  It also canonicalizes
the four projected (three random axes plus b) cells under S3, producing the
variable keys needed by the projected-dual LP.

The upstream C++ enumerator is responsible for exhaustive vertex/intersection
and S4-orbit discovery.  This script independently rechecks every reported
cell and all downstream face/triangulation algebra using Fraction.
"""
from __future__ import annotations

import argparse
import itertools
import json
from fractions import Fraction as F
from pathlib import Path


DIM = 5
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
DEFAULT_FLAT = (ROOT / "work/emc_audit/computer/"
                "r5_zs5_d4_relevant_cell_orbits_exact_flat.txt")
DEFAULT_JSON = (ROOT / "work/emc_audit/computer/"
                "r5_zs5_d4_orbit_geometry_exact.json")


def dot(a, x):
    return sum((q * r for q, r in zip(a, x)), F(0))


def matrix_rank(rows):
    a = [list(map(F, row)) for row in rows if any(row)]
    if not a:
        return 0
    rank = 0
    for column in range(len(a[0])):
        pivot = next((r for r in range(rank, len(a)) if a[r][column]), None)
        if pivot is None:
            continue
        a[rank], a[pivot] = a[pivot], a[rank]
        value = a[rank][column]
        a[rank] = [q / value for q in a[rank]]
        for r in range(len(a)):
            if r != rank and a[r][column]:
                value = a[r][column]
                a[r] = [q - value * p for q, p in zip(a[r], a[rank])]
        rank += 1
    return rank


def affine_dim(points):
    if len(points) <= 1:
        return len(points) - 1
    origin = points[0]
    return matrix_rank([[q - r for q, r in zip(point, origin)]
                        for point in points[1:]])


def parse_flat(path):
    lines = path.read_text().splitlines()
    assert lines[0] == "VERTICES 587"
    vertices = [tuple(F(q) for q in line.split()) for line in lines[1:588]]
    assert lines[588].startswith("ORBITS ") and lines[-1] == "END"
    orbit_count = int(lines[588].split()[1])
    keys = lines[589:-1]
    assert len(keys) == orbit_count and len(set(keys)) == len(keys)
    return vertices, keys


def decode_key(key):
    assert len(key) == 36
    b_box = int(key[0])
    box = tuple(int(q) for q in key[1:5])
    signs = tuple(key[5 + mask - 1] for mask in range(1, 32))
    return b_box, box, signs


def relevant_level(mask, box, b_box, denominator=5, dimension=5):
    size = mask.bit_count()
    if size < 2:
        return None
    base = tuple(box) + (b_box,)
    level = denominator - sum(base[i] for i in range(dimension) if mask & (1 << i))
    return level if 1 <= level < size else None


def constraints_for(key, denominator=5):
    b_box, box, signs = decode_key(key)
    constraints = []  # normal dot x <= rhs
    for axis in range(DIM):
        e = tuple(int(i == axis) for i in range(DIM))
        constraints.append((tuple(-q for q in e), F(0)))
        constraints.append((e, F(1)))
    for mask in range(1, 32):
        level = relevant_level(mask, box, b_box, denominator)
        sign = signs[mask - 1]
        assert (level is None) == (sign == "0")
        if level is None:
            continue
        normal = tuple(int(mask & (1 << i) != 0) for i in range(DIM))
        if sign == "-":
            constraints.append((normal, F(level)))
        else:
            assert sign == "+"
            constraints.append((tuple(-q for q in normal), F(-level)))
    return tuple(constraints)


def true_vertex_ids(vertices, constraints):
    answer = []
    for vertex_id, point in enumerate(vertices):
        if any(dot(normal, point) > rhs for normal, rhs in constraints):
            continue
        active = [normal for normal, rhs in constraints if dot(normal, point) == rhs]
        if matrix_rank(active) == DIM:
            answer.append(vertex_id)
    return tuple(answer)


def face_facets(vertex_ids, vertices, constraints):
    dimension = affine_dim([vertices[i] for i in vertex_ids])
    facets = set()
    for normal, rhs in constraints:
        face = tuple(i for i in vertex_ids if dot(normal, vertices[i]) == rhs)
        if face != vertex_ids and len(face) >= dimension:
            if affine_dim([vertices[i] for i in face]) == dimension - 1:
                facets.add(face)
    return tuple(sorted(facets))


def pulling(vertex_ids, vertices, constraints, memo):
    key = tuple(vertex_ids)
    if key in memo:
        return memo[key]
    dimension = affine_dim([vertices[i] for i in vertex_ids])
    if len(vertex_ids) == dimension + 1:
        answer = (key,)
    else:
        apex = min(vertex_ids)
        simplices = []
        for facet in face_facets(vertex_ids, vertices, constraints):
            if apex in facet:
                continue
            for simplex in pulling(facet, vertices, constraints, memo):
                simplices.append(tuple(sorted((apex,) + simplex)))
        answer = tuple(sorted(set(simplices)))
        if not answer:
            raise AssertionError("empty nonsimplex pulling triangulation")
    memo[key] = answer
    return answer


def permute_local_mask(mask, permutation):
    result = mask & (1 << 3)  # b is fixed local axis 3
    for old in range(3):
        if mask & (1 << old):
            result |= 1 << permutation[old]
    return result


def projected_descriptor(full_key, omitted):
    b_box, box, full_signs = decode_key(full_key)
    remaining = tuple(i for i in range(4) if i != omitted)
    local_signs = ["0"] * 15
    for local_mask in range(1, 16):
        global_mask = 0
        for local_axis in range(3):
            if local_mask & (1 << local_axis):
                global_mask |= 1 << remaining[local_axis]
        if local_mask & 8:
            global_mask |= 1 << 4
        local_signs[local_mask - 1] = full_signs[global_mask - 1]

    best = None
    best_axes = None
    for permutation in itertools.permutations(range(3)):
        permuted_box = [0] * 3
        permuted_axes = [None] * 3
        for old in range(3):
            permuted_box[permutation[old]] = box[remaining[old]]
            permuted_axes[permutation[old]] = remaining[old]
        permuted_signs = ["0"] * 15
        for mask in range(1, 16):
            permuted_signs[permute_local_mask(mask, permutation) - 1] = local_signs[mask - 1]
        descriptor = str(b_box) + "".join(map(str, permuted_box)) + "".join(permuted_signs)
        if best is None or descriptor < best:
            best = descriptor
            best_axes = tuple(permuted_axes)
    return best, best_axes


def build(flat_path, denominator=5):
    vertices, keys = parse_flat(flat_path)
    projected_keys = set()
    cells = []
    simplex_distribution = {}
    vertex_distribution = {}
    total_simplices = 0
    for index, key in enumerate(keys):
        constraints = constraints_for(key, denominator)
        vertex_ids = true_vertex_ids(vertices, constraints)
        assert affine_dim([vertices[i] for i in vertex_ids]) == DIM
        memo = {}
        simplices = pulling(vertex_ids, vertices, constraints, memo)
        for simplex in simplices:
            assert len(simplex) == DIM + 1
            assert affine_dim([vertices[i] for i in simplex]) == DIM
        projections = []
        for omitted in range(4):
            projected_key, axes = projected_descriptor(key, omitted)
            projected_keys.add(projected_key)
            projections.append({"key": projected_key, "axes": list(axes)})
        cells.append({
            "key": key,
            "vertices": list(vertex_ids),
            "simplices": [list(q) for q in simplices],
            "projections": projections,
        })
        total_simplices += len(simplices)
        vertex_distribution[len(vertex_ids)] = vertex_distribution.get(len(vertex_ids), 0) + 1
        simplex_distribution[len(simplices)] = simplex_distribution.get(len(simplices), 0) + 1
        if (index + 1) % 500 == 0:
            print("processed", index + 1, "cells", flush=True)
    payload = {
        "schema": "r5_zs5_d4_orbit_geometry_exact_v1",
        "vertices": [[str(q) for q in point] for point in vertices],
        "projected_keys": sorted(projected_keys),
        "cells": cells,
        "summary": {
            "universal_vertex_count": len(vertices),
            "s4_cell_orbit_count": len(cells),
            "projected_s3_cell_count": len(projected_keys),
            "pulling_simplex_orbit_count": total_simplices,
            "cell_vertex_distribution": {str(k): v for k, v in sorted(vertex_distribution.items())},
            "cell_simplex_distribution": {str(k): v for k, v in sorted(simplex_distribution.items())},
        },
    }
    if denominator != 5:
        payload["denominator"] = denominator
    return payload


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--flat", type=Path, default=DEFAULT_FLAT)
    parser.add_argument("--output", type=Path, default=DEFAULT_JSON)
    parser.add_argument("--denominator", type=int, default=5)
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    assert args.denominator in (5, 10)
    rebuilt = build(args.flat, args.denominator)
    if args.verify:
        recorded = json.loads(args.output.read_text())
        if json.dumps(recorded, sort_keys=True) != json.dumps(rebuilt, sort_keys=True):
            raise AssertionError("geometry manifest mismatch")
        print(json.dumps(rebuilt["summary"], indent=2))
        print("CERTIFIED")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(rebuilt, indent=2) + "\n")
        print(json.dumps(rebuilt["summary"], indent=2))
        print("wrote", args.output)


if __name__ == "__main__":
    main()
