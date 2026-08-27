#!/usr/bin/env python3
"""Independent exact red-team audit for the d=4 local-quadratic certificate.

This checker deliberately does not solve for a dual.  It verifies the glue
between the floating active-row discovery, the exact Q[z]/(cusp) active-face
repair, the S3 stabilizer quotient, and the outward dyadic interval flat file.

The optional floating scan is only a reproducibility check for the claim that
the recorded 12,794 rows are *all* rows below the stated 1e-9 threshold among
the 5,164,740 degree-four rows.  Mathematical validity does not rely on those
rows having been exact zeros before repair: this script reconstructs their
exact equations and checks that the repaired point lies in their intersection.
The remaining rows must still be checked by the separate full interval replay.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
import re
import subprocess
import tempfile
from collections import Counter, defaultdict
from fractions import Fraction as F
from pathlib import Path


M = F(1, 5)
BASIS = 15
CUSP = (F(-1024), F(9225), F(15250), F(11250), F(3125))
BASIS_EXPONENTS = (
    (0, 0, 0, 0),
    (1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1),
    (2, 0, 0, 0), (1, 1, 0, 0), (1, 0, 1, 0), (1, 0, 0, 1),
    (0, 2, 0, 0), (0, 1, 1, 0), (0, 1, 0, 1),
    (0, 0, 2, 0), (0, 0, 1, 1), (0, 0, 0, 2),
)
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
COMPUTER = ROOT / "work/emc_audit/computer"
GEOMETRY = COMPUTER / "r5_zs5_d4_orbit_geometry_exact.json"
GEOMETRY_FLAT = COMPUTER / "r5_zs5_d4_quadratic_geometry_exact_flat.txt"
FLOATING = COMPUTER / "r5_zs5_d4_projected_quadratic_symmetric_zero_highs_probe.json"
ACTIVE = COMPUTER / "r5_zs5_d4_projected_quadratic_active_audit.json"
SYSTEM = COMPUTER / "r5_zs5_d4_projected_quadratic_active_exact_system.json"
REPAIR = COMPUTER / "r5_zs5_d4_projected_quadratic_face_repair_exact.json"
INTERVAL_FLAT = COMPUTER / "r5_zs5_d4_projected_quadratic_exact_flat.txt"
SCANNER = HERE / "r5_zs5_d4_projected_quadratic_scan.cpp"
EXACT_VERIFIER = HERE / "r5_zs5_d4_projected_quadratic_exact_verify.cpp"
REPAIR_GENERATOR = HERE / "r5_zs5_d4_quadratic_face_repair_exact.py"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class UnionFind:
    def __init__(self, size: int):
        self.parent = list(range(size))

    def find(self, value: int) -> int:
        while self.parent[value] != value:
            self.parent[value] = self.parent[self.parent[value]]
            value = self.parent[value]
        return value

    def union(self, left: int, right: int) -> None:
        left, right = self.find(left), self.find(right)
        if left == right:
            return
        if right < left:
            left, right = right, left
        self.parent[right] = left


def permute_mask(mask: int, permutation: tuple[int, ...]) -> int:
    result = mask & 8  # projected b-axis is fixed
    for old in range(3):
        if mask & (1 << old):
            result |= 1 << permutation[old]
    return result


def permuted_projected_key(key: str, permutation: tuple[int, ...]) -> str:
    assert len(key) == 19
    boxes = key[1:4]
    signs = key[4:]
    new_boxes = [""] * 3
    new_signs = [""] * 15
    for old in range(3):
        new_boxes[permutation[old]] = boxes[old]
    for mask in range(1, 16):
        new_signs[permute_mask(mask, permutation) - 1] = signs[mask - 1]
    return key[0] + "".join(new_boxes) + "".join(new_signs)


def stabilizer_mapping(projected_keys: list[str]):
    union = UnionFind(BASIS * len(projected_keys))
    exponent_where = {value: index for index, value in enumerate(BASIS_EXPONENTS)}
    stabilizer_distribution = Counter()
    for key_index, key in enumerate(projected_keys):
        stabilizers = [p for p in itertools.permutations(range(3))
                       if permuted_projected_key(key, p) == key]
        assert stabilizers
        stabilizer_distribution[len(stabilizers)] += 1
        for permutation in stabilizers:
            for basis, exponent in enumerate(BASIS_EXPONENTS):
                transformed = [0, 0, 0, exponent[3]]
                for old in range(3):
                    transformed[permutation[old]] = exponent[old]
                other = exponent_where[tuple(transformed)]
                union.union(BASIS * key_index + basis,
                            BASIS * key_index + other)
    roots = sorted({union.find(q) for q in range(len(union.parent))})
    root_where = {root: index for index, root in enumerate(roots)}
    mapping = [root_where[union.find(q)] for q in range(len(union.parent))]
    return mapping, roots, stabilizer_distribution


def decode_full_key(key: str):
    assert len(key) == 36
    return int(key[0]), tuple(int(q) for q in key[1:5]), key[5:]


def projected_descriptor_candidates(full_key: str, omitted: int):
    b_box, box, full_signs = decode_full_key(full_key)
    remaining = tuple(q for q in range(4) if q != omitted)
    local_signs = ["0"] * 15
    for local_mask in range(1, 16):
        global_mask = 0
        for local_axis in range(3):
            if local_mask & (1 << local_axis):
                global_mask |= 1 << remaining[local_axis]
        if local_mask & 8:
            global_mask |= 1 << 4
        local_signs[local_mask - 1] = full_signs[global_mask - 1]
    candidates = []
    for permutation in itertools.permutations(range(3)):
        new_box = [None] * 3
        new_axes = [None] * 3
        new_signs = [None] * 15
        for old in range(3):
            new_box[permutation[old]] = box[remaining[old]]
            new_axes[permutation[old]] = remaining[old]
        for mask in range(1, 16):
            new_signs[permute_mask(mask, permutation) - 1] = local_signs[mask - 1]
        descriptor = str(b_box) + "".join(map(str, new_box)) + "".join(new_signs)
        candidates.append((descriptor, tuple(new_axes)))
    best = min(descriptor for descriptor, _ in candidates)
    return {(descriptor, axes) for descriptor, axes in candidates if descriptor == best}


def rebuilt_flat_geometry(data) -> bytes:
    key_where = {key: index for index, key in enumerate(data["projected_keys"])}
    lines = [f"VERTICES {len(data['vertices'])}"]
    lines.extend(" ".join(point) for point in data["vertices"])
    lines.append(f"PROJECTED_KEYS {len(key_where)}")
    lines.extend(data["projected_keys"])
    simplex_count = sum(len(cell["simplices"]) for cell in data["cells"])
    lines.append(f"SIMPLICES {simplex_count}")
    for cell in data["cells"]:
        b_box, box, _ = decode_full_key(cell["key"])
        prefix = [str(b_box), *(str(q) for q in box)]
        for projection in cell["projections"]:
            prefix.extend((str(key_where[projection["key"]]),
                           *(str(q) for q in projection["axes"])))
        for simplex in cell["simplices"]:
            lines.append(" ".join((*prefix, *(str(q) for q in simplex))))
    lines.append("END")
    return ("\n".join(lines) + "\n").encode()


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
        for j in range(4):
            poly[shift + j] -= lead * CUSP[j] / CUSP[4]
    poly.extend([F(0)] * (4 - len(poly)))
    return tuple(poly)


def global_point(local, box, b_box):
    return tuple((F(box[j]) + local[j]) / 5 for j in range(4)) + (
        (F(b_box) + local[4]) / 5,)


def kernel_signature(point):
    signature = []
    for state in range(16):
        subtotal = sum((point[j] for j in range(4) if state & (1 << j)), F(0))
        signature.append((subtotal > 1, subtotal + point[4] > 1))
    return tuple(signature)


def kernel_polynomials(point):
    sum_zero = [F(0)] * 5
    sum_one = [F(0)] * 5
    for state in range(16):
        subtotal = sum((point[j] for j in range(4) if state & (1 << j)), F(0))
        degree = 4 - state.bit_count()
        if subtotal > 1:
            sum_zero[degree] += 1
        if subtotal + point[4] > 1:
            sum_one[degree] += 1
    return sum_zero, sum_one


def basis_blossom(basis, h, r, s):
    if basis == 0:
        return F(1)
    if basis <= 4:
        axis = basis - 1
        return (h[axis][r] + h[axis][s]) / 2
    at = 5
    for left in range(4):
        for right in range(left, 4):
            if at == basis:
                if left == right:
                    return h[left][r] * h[left][s]
                return (h[left][r] * h[right][s]
                        + h[left][s] * h[right][r]) / 2
            at += 1
    raise AssertionError("bad quadratic basis")


def product_blossom(b_values, u_values, h, basis, slots):
    answer = F(0)
    for b_position in range(4):
        for u_position in range(4):
            if b_position == u_position:
                continue
            remaining = [q for q in range(4)
                         if q not in (b_position, u_position)]
            answer += (b_values[slots[b_position]] * u_values[slots[u_position]]
                       * basis_blossom(basis, h,
                                       slots[remaining[0]], slots[remaining[1]]))
    return answer / 12


def reconstruct_row(cell, simplex, alpha, vertices, key_where, mapping):
    b_box, box, _ = decode_full_key(cell["key"])
    cell_points = [global_point(vertices[q], box, b_box) for q in cell["vertices"]]
    cell_center = tuple(sum((point[j] for point in cell_points), F(0)) / len(cell_points)
                        for j in range(5))
    sum_zero, sum_one = kernel_polynomials(cell_center)
    points = [global_point(vertices[q], box, b_box) for q in simplex]
    slots = tuple(vertex for vertex, multiplicity in enumerate(alpha)
                  for _ in range(multiplicity))
    assert len(slots) == 4
    b_values = tuple(point[4] for point in points)
    b_average = sum((b_values[q] for q in slots), F(0)) / 4
    rhs = [F(0)] * 6
    rhs[0] += b_average
    add_poly(rhs, sum_zero, M - b_average)
    add_poly(rhs, sum_zero, -b_average, 1)
    add_poly(rhs, sum_one, -M)
    for axis in range(4):
        if all(points[q][axis] == 1 for q in slots):
            rhs[3] += M - b_average
            rhs[4] -= b_average
    if all(points[q][4] == 1 for q in slots):
        rhs[4] -= M

    row = defaultdict(F)
    for omitted in range(4):
        projection = cell["projections"][omitted]
        axes = tuple(projection["axes"])
        assert omitted not in axes and set(axes) == set(range(4)) - {omitted}
        key_index = key_where[projection["key"]]
        u_values = tuple(point[omitted] - M for point in points)
        h = tuple(tuple(point[axis] - M for point in points) for axis in axes) + (
            tuple(point[4] - M for point in points),)
        for basis in range(BASIS):
            value = product_blossom(b_values, u_values, h, basis, slots)
            row[mapping[BASIS * key_index + basis]] += value
    exact_row = tuple(sorted((column, value) for column, value in row.items() if value))
    return exact_row, reduce_cusp(rhs), b_box


def parse_recorded_equations(system):
    answer = Counter()
    for item in system["equations"]:
        row = tuple((int(column), F(value)) for column, value in item["row"])
        rhs = tuple(F(q) for q in item["rhs_mod_cusp"])
        answer[(row, rhs, int(item["b_box"]))] += int(item["multiplicity"])
    return answer


def eval_poly(poly, x):
    answer = F(0)
    for coefficient in reversed(poly):
        answer = answer * x + coefficient
    return answer


def interval_eval(poly, lower, upper):
    low = high = poly[-1]
    for coefficient in reversed(poly[:-1]):
        products = (low * lower, low * upper, high * lower, high * upper)
        low = min(products) + coefficient
        high = max(products) + coefficient
    return low, high


def cusp_interval(steps=64):
    lower, upper = F(19, 200), F(12, 125)
    assert eval_poly(CUSP, lower) < 0 < eval_poly(CUSP, upper)
    # P'(z)>0 termwise on this positive interval, so exactly one root is enclosed.
    assert F(9225) + F(30500) * lower + F(33750) * lower**2 + F(12500) * lower**3 > 0
    for _ in range(steps):
        middle = (lower + upper) / 2
        if eval_poly(CUSP, middle) < 0:
            lower = middle
        else:
            upper = middle
    return lower, upper


def floor_fraction(value):
    return value.numerator // value.denominator


def ceil_fraction(value):
    return -floor_fraction(-value)


def parse_interval_flat(path: Path):
    tokens = path.read_text().split()
    at = 0

    def take(expected=None):
        nonlocal at
        value = tokens[at]
        at += 1
        if expected is not None:
            assert value == expected, (value, expected)
        return value

    take("R5_ZS5_D4_QUADRATIC_EXACT_INTERVAL_V1")
    take("SCALE_BITS"); bits = int(take())
    hashes = {}
    for label in ("REPAIR_SHA256", "ACTIVE_SHA256", "GEOMETRY_SHA256"):
        take(label); hashes[label] = take()
    take("CUSP_INTERVAL"); root_interval = (F(take()), F(take()))
    take("POWERS"); power_count = int(take())
    powers = [(int(take()), int(take())) for _ in range(power_count)]
    take("COEFFICIENT_INTERVALS"); coefficient_count = int(take())
    coefficients = [(int(take()), int(take())) for _ in range(coefficient_count)]
    take("ORIGINAL_TO_COLLAPSED"); mapping_count = int(take())
    mapping = [int(take()) for _ in range(mapping_count)]
    take("ACTIVE_ROWS"); active_count = int(take())
    active_rows = []
    for _ in range(active_count):
        active_rows.append((int(take()), tuple(int(take()) for _ in range(6))))
    take("END")
    assert at == len(tokens)
    return bits, hashes, root_interval, powers, coefficients, mapping, active_rows


def audit_interval_flat(path, repair, active, geometry_hash, mapping):
    bits, hashes, root_interval, powers, intervals, flat_mapping, flat_active = \
        parse_interval_flat(path)
    assert bits == 90
    assert hashes == {
        "REPAIR_SHA256": sha256(REPAIR),
        "ACTIVE_SHA256": sha256(ACTIVE),
        "GEOMETRY_SHA256": sha256(GEOMETRY_FLAT),
    }
    lower, upper = cusp_interval()
    assert root_interval == (lower, upper)
    scale = 1 << bits
    expected_powers = [(floor_fraction(lower**degree * scale),
                        ceil_fraction(upper**degree * scale))
                       for degree in range(4)]
    assert powers == expected_powers
    polynomials = [[F(q) for q in poly]
                   for poly in repair["collapsed_coefficients_qz_mod_cusp"]]
    expected_intervals = []
    maximum_width = 0
    for poly in polynomials:
        low, high = interval_eval(poly, lower, upper)
        enclosure = (floor_fraction(low * scale), ceil_fraction(high * scale))
        expected_intervals.append(enclosure)
        maximum_width = max(maximum_width, enclosure[1] - enclosure[0])
    assert intervals == expected_intervals
    assert flat_mapping == mapping
    recorded_active = sorted((int(item["simplex"]), tuple(item["alpha"]))
                             for item in active["active_rows"])
    assert flat_active == recorded_active
    return {
        "scale_bits": bits,
        "coefficient_interval_count": len(intervals),
        "maximum_scaled_coefficient_width": maximum_width,
        "active_row_count": len(flat_active),
        "flat_sha256": sha256(path),
    }


def floating_scan(active, floating, tolerance, compiler):
    with tempfile.TemporaryDirectory(prefix="r5_d4_redteam_") as directory:
        directory = Path(directory)
        executable = directory / "scan"
        coefficient_path = directory / "coefficients.txt"
        coefficients = [value for item in floating["g"] for value in item["coefficients"]]
        coefficient_path.write_text("\n".join(f"{value:.19g}" for value in coefficients) + "\n")
        subprocess.run([compiler, "-O3", "-std=c++17", str(SCANNER), "-o", str(executable)],
                       check=True)
        process = subprocess.Popen(
            [str(executable), str(GEOMETRY_FLAT), str(coefficient_path), "20000", "4"],
            text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        assert process.stdout is not None
        summary = None
        discovered = []
        returned_values = []
        for line in process.stdout:
            tokens = line.split()
            if not tokens:
                continue
            if tokens[0] == "SUMMARY":
                summary = (int(tokens[2]), float(tokens[4]), int(tokens[6]))
            elif tokens[0] == "ROW":
                value = float(tokens[1])
                sid = int(tokens[3])
                assert tokens[4] == "A"
                alpha = tuple(int(q) for q in tokens[5:11])
                returned_values.append(value)
                if abs(value) < tolerance:
                    discovered.append((sid, alpha))
        stderr = process.stderr.read() if process.stderr is not None else ""
        return_code = process.wait()
        if return_code:
            raise RuntimeError(f"floating scanner failed ({return_code}): {stderr[-2000:]}")
    assert summary == (5_164_740, summary[1], 20_000)
    assert len(returned_values) == 20_000
    assert returned_values == sorted(returned_values)
    expected = sorted((int(item["simplex"]), tuple(item["alpha"]))
                      for item in active["active_rows"])
    assert sorted(discovered) == expected
    nonactive_values = [value for value in returned_values if abs(value) >= tolerance]
    assert nonactive_values and min(nonactive_values) >= tolerance
    return {
        "all_row_count": summary[0],
        "returned_worst_count": summary[2],
        "active_below_tolerance_count": len(discovered),
        "maximum_active_absolute_gap": max(abs(value) for value in returned_values
                                           if abs(value) < tolerance),
        "minimum_returned_nonactive_gap": min(nonactive_values),
        "largest_returned_gap": returned_values[-1],
        "scanner_sha256": sha256(SCANNER),
    }


def rerun_repair_byte_exact(python_executable):
    with tempfile.TemporaryDirectory(prefix="r5_d4_repair_redteam_") as directory:
        output = Path(directory) / "repair.json"
        subprocess.run([
            python_executable, str(REPAIR_GENERATOR), "--output", str(output),
        ], cwd=ROOT, check=True)
        assert output.read_bytes() == REPAIR.read_bytes()
    return {
        "byte_identical": True,
        "repair_sha256": sha256(REPAIR),
        "generator_sha256": sha256(REPAIR_GENERATOR),
    }


def full_interval_replay(compiler, use_ubsan):
    with tempfile.TemporaryDirectory(prefix="r5_d4_full_redteam_") as directory:
        executable = Path(directory) / "exact_verify"
        flags = ["-O1", "-g", "-fsanitize=undefined",
                 "-fno-sanitize-recover=undefined"] if use_ubsan else ["-O3"]
        subprocess.run([compiler, *flags, "-std=c++17", str(EXACT_VERIFIER),
                        "-o", str(executable)], cwd=ROOT, check=True)
        completed = subprocess.run([
            str(executable), str(GEOMETRY_FLAT), str(INTERVAL_FLAT),
        ], cwd=ROOT, check=True, text=True, capture_output=True)
    stdout = completed.stdout
    required = {
        "rows": r"rows (\d+)",
        "strict_rows": r"strict_rows (\d+)",
        "active_rows": r"exact_active_rows (\d+)",
        "minimum_strict_lower_scaled": r"minimum_strict_lower_scaled (\d+)",
        "scale_bits": r"scale_bits (\d+)",
        "largest_active_interval_width_scaled":
            r"largest_active_interval_width_scaled (\d+)",
        "row_stream_fnv1a64": r"row_stream_fnv1a64 ([0-9a-f]{16})",
    }
    parsed = {}
    for key, pattern in required.items():
        match = re.search(pattern, stdout)
        assert match, (key, stdout)
        parsed[key] = match.group(1)
    for key in required:
        if key != "row_stream_fnv1a64":
            parsed[key] = int(parsed[key])
    assert "R5_ZS5_D4_QUADRATIC_EXACT_INTERVAL_PASS" in stdout
    assert "closed_top_face_corrections INCLUDED" in stdout
    assert parsed == {
        "rows": 5_164_740,
        "strict_rows": 5_151_946,
        "active_rows": 12_794,
        "minimum_strict_lower_scaled": 5_411_340_156_129_015_928_571,
        "scale_bits": 90,
        "largest_active_interval_width_scaled": 267_468_378,
        "row_stream_fnv1a64": "49063aa32fa17a6f",
    }
    parsed["undefined_behavior_sanitizer_enabled"] = use_ubsan
    parsed["verifier_sha256"] = sha256(EXACT_VERIFIER)
    parsed["stdout_sha256"] = hashlib.sha256(stdout.encode()).hexdigest()
    return parsed


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--interval-flat", type=Path, default=INTERVAL_FLAT)
    parser.add_argument("--skip-floating-scan", action="store_true")
    parser.add_argument("--compiler", default="c++")
    parser.add_argument("--rerun-repair", action="store_true")
    parser.add_argument("--full-replay", action="store_true")
    parser.add_argument("--ubsan", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    geometry = json.loads(GEOMETRY.read_text())
    floating = json.loads(FLOATING.read_text())
    active = json.loads(ACTIVE.read_text())
    system = json.loads(SYSTEM.read_text())
    repair = json.loads(REPAIR.read_text())

    assert geometry["summary"] == {
        **geometry["summary"],
        "universal_vertex_count": 587,
        "s4_cell_orbit_count": 3493,
        "projected_s3_cell_count": 545,
        "pulling_simplex_orbit_count": 40990,
    }
    assert len(geometry["vertices"]) == 587
    assert len(geometry["projected_keys"]) == 545
    assert len(geometry["cells"]) == 3493
    assert sum(len(cell["simplices"]) for cell in geometry["cells"]) == 40990
    assert rebuilt_flat_geometry(geometry) == GEOMETRY_FLAT.read_bytes()

    projected_keys = geometry["projected_keys"]
    key_where = {key: index for index, key in enumerate(projected_keys)}
    assert len(key_where) == 545
    projection_count = 0
    ambiguous_projection_count = 0
    for cell in geometry["cells"]:
        for omitted, projection in enumerate(cell["projections"]):
            candidates = projected_descriptor_candidates(cell["key"], omitted)
            actual = (projection["key"], tuple(projection["axes"]))
            assert actual in candidates
            assert omitted not in projection["axes"]
            assert set(projection["axes"]) == set(range(4)) - {omitted}
            projection_count += 1
            ambiguous_projection_count += len(candidates) > 1

    mapping, roots, stabilizer_distribution = stabilizer_mapping(projected_keys)
    assert len(mapping) == 8175 and len(roots) == 6759
    assert mapping == system["original_to_collapsed_variable"]
    assert roots == system["collapsed_roots"]

    assert system["geometry_sha256"] == sha256(GEOMETRY)
    assert system["active_source_sha256"] == sha256(ACTIVE)
    assert repair["system_sha256"] == sha256(SYSTEM)
    assert repair["floating_source_sha256"] == sha256(FLOATING)
    assert active["active_count"] == len(active["active_rows"]) == 12794
    targets = {(int(item["simplex"]), tuple(item["alpha"]))
               for item in active["active_rows"]}
    assert len(targets) == 12794
    by_sid = defaultdict(list)
    for sid, alpha in targets:
        assert 0 <= sid < 40990 and sum(alpha) == 4 and len(alpha) == 6
        by_sid[sid].append(alpha)

    vertices = [tuple(F(q) for q in point) for point in geometry["vertices"]]
    reconstructed = Counter()
    signature_mismatch = 0
    sid = 0
    for cell in geometry["cells"]:
        b_box, box, _ = decode_full_key(cell["key"])
        cell_points = [global_point(vertices[q], box, b_box) for q in cell["vertices"]]
        cell_center = tuple(sum((point[j] for point in cell_points), F(0)) / len(cell_points)
                            for j in range(5))
        cell_signature = kernel_signature(cell_center)
        for simplex in cell["simplices"]:
            simplex_points = [global_point(vertices[q], box, b_box) for q in simplex]
            simplex_center = tuple(sum((point[j] for point in simplex_points), F(0)) / 6
                                   for j in range(5))
            if kernel_signature(simplex_center) != cell_signature:
                signature_mismatch += 1
            for alpha in by_sid.get(sid, ()):
                reconstructed[reconstruct_row(cell, simplex, alpha, vertices,
                                              key_where, mapping)] += 1
            sid += 1
    assert sid == 40990 and signature_mismatch == 0
    assert sum(reconstructed.values()) == 12794
    recorded = parse_recorded_equations(system)
    assert reconstructed == recorded
    assert len(reconstructed) == system["unique_exact_equation_count"] == 1620

    coefficients = [[F(q) for q in poly]
                    for poly in repair["collapsed_coefficients_qz_mod_cusp"]]
    expanded = [[F(q) for q in poly]
                for poly in repair["expanded_coefficients_qz_mod_cusp"]]
    assert len(coefficients) == 6759 and all(len(poly) == 4 for poly in coefficients)
    assert len(expanded) == 8175
    assert expanded == [coefficients[column] for column in mapping]
    replay_failures = []
    for row, rhs, _ in reconstructed:
        residual = list(rhs)
        for column, value in row:
            for degree in range(4):
                residual[degree] += value * coefficients[column][degree]
        if any(residual):
            replay_failures.append(tuple(residual))
    assert not replay_failures
    assert repair["exact_active_failure_count"] == 0

    interval_summary = audit_interval_flat(
        args.interval_flat, repair, active, sha256(GEOMETRY_FLAT), mapping)
    scan_summary = None
    if not args.skip_floating_scan:
        scan_summary = floating_scan(
            active, floating, F(str(active["active_tolerance"])), args.compiler)
    repair_replay = rerun_repair_byte_exact(str(Path(__import__("sys").executable))) \
        if args.rerun_repair else None
    full_replay = full_interval_replay(args.compiler, args.ubsan) \
        if args.full_replay else None

    payload = {
        "schema": "r5_zs5_d4_quadratic_redteam_exact_v1",
        "proof_status": (
            "PASS: independent active-face, interval-flat, and full-row replay"
            if full_replay else
            "PASS: active-face glue and interval-flat generation; "
            "full nonactive positivity remains the separate C++ replay"),
        "artifact_sha256": {
            "geometry_json": sha256(GEOMETRY),
            "geometry_flat": sha256(GEOMETRY_FLAT),
            "floating_candidate": sha256(FLOATING),
            "active_audit": sha256(ACTIVE),
            "active_system": sha256(SYSTEM),
            "face_repair": sha256(REPAIR),
            "interval_flat": sha256(args.interval_flat),
        },
        "geometry": {
            "vertices": 587,
            "full_cell_orbits": 3493,
            "simplex_orbits": 40990,
            "projected_keys": 545,
            "projection_records_checked": projection_count,
            "ambiguous_canonical_projection_records": ambiguous_projection_count,
            "simplex_cell_kernel_signature_mismatch_count": signature_mismatch,
        },
        "stabilizer": {
            "original_variables": 8175,
            "collapsed_variables": 6759,
            "key_stabilizer_size_distribution": dict(sorted(stabilizer_distribution.items())),
            "mapping_exact_match": True,
        },
        "active_face": {
            "recorded_rows": len(targets),
            "reconstructed_rows_with_multiplicity": sum(reconstructed.values()),
            "unique_exact_equations": len(reconstructed),
            "counter_exact_match": True,
            "exact_replay_failure_count": len(replay_failures),
            "selected_near_zero_rows_need_not_have_been_symbolic_zeros": True,
        },
        "interval_flat": interval_summary,
        "floating_exhaustion_replay": scan_summary,
        "deterministic_face_repair_replay": repair_replay,
        "full_interval_replay": full_replay,
    }
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))
    print("REDTEAM ACTIVE-FACE/INTERVAL-FLAT PASS")


if __name__ == "__main__":
    main()
