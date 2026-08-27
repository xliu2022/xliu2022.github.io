#!/usr/bin/env python3
"""Generate the flat exact interval certificate for the repaired d=4 dual."""
from __future__ import annotations

import argparse
import hashlib
import json
from fractions import Fraction as F
from pathlib import Path


CUSP = (F(-1024), F(9225), F(15250), F(11250), F(3125))
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
REPAIR = (ROOT / "work/emc_audit/computer/"
          "r5_zs5_d4_projected_quadratic_face_repair_exact.json")
ACTIVE = (ROOT / "work/emc_audit/computer/"
          "r5_zs5_d4_projected_quadratic_active_audit.json")
GEOMETRY = (ROOT / "work/emc_audit/computer/"
            "r5_zs5_d4_quadratic_geometry_exact_flat.txt")
OUTPUT = (ROOT / "work/emc_audit/computer/"
          "r5_zs5_d4_projected_quadratic_exact_flat.txt")


def eval_poly(poly, x):
    answer = F(0)
    for coefficient in reversed(poly):
        answer = answer * x + coefficient
    return answer


def interval_eval(poly, lower, upper):
    low = high = poly[-1]
    for coefficient in reversed(poly[:-1]):
        products = (low * lower, low * upper, high * lower, high * upper)
        low, high = min(products) + coefficient, max(products) + coefficient
    return low, high


def cusp_interval(steps=64):
    lower, upper = F(19, 200), F(12, 125)
    # P'(x)=9225+30500x+33750x^2+12500x^3 is positive for x>0.
    # Thus the endpoint sign change isolates exactly one real root, and every
    # bisection update below retains that root.
    assert lower > 0 and all(coefficient > 0 for coefficient in CUSP[1:])
    assert eval_poly(CUSP, lower) < 0 < eval_poly(CUSP, upper)
    for _ in range(steps):
        middle = (lower + upper) / 2
        if eval_poly(CUSP, middle) < 0:
            lower = middle
        else:
            upper = middle
    assert eval_poly(CUSP, lower) < 0 < eval_poly(CUSP, upper)
    return lower, upper


def floor_fraction(value):
    return value.numerator // value.denominator


def ceil_fraction(value):
    return -floor_fraction(-value)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--scale-bits", type=int, default=90)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    repair = json.loads(REPAIR.read_text())
    active = json.loads(ACTIVE.read_text())
    assert repair["exact_active_failure_count"] == 0
    system_path = (ROOT / "work/emc_audit/computer/"
                   "r5_zs5_d4_projected_quadratic_active_exact_system.json")
    assert repair["system_sha256"] == hashlib.sha256(system_path.read_bytes()).hexdigest()
    assert active["active_count"] == 12794
    polynomials = [[F(q) for q in poly]
                   for poly in repair["collapsed_coefficients_qz_mod_cusp"]]
    mapping = json.loads(system_path.read_text())["original_to_collapsed_variable"]
    lower, upper = cusp_interval()
    scale = 1 << args.scale_bits
    intervals = []
    maximum_width = 0
    for poly in polynomials:
        low, high = interval_eval(poly, lower, upper)
        low_scaled = floor_fraction(low * scale)
        high_scaled = ceil_fraction(high * scale)
        intervals.append((low_scaled, high_scaled))
        maximum_width = max(maximum_width, high_scaled - low_scaled)
    powers = []
    for degree in range(4):
        low = lower ** degree
        high = upper ** degree
        powers.append((floor_fraction(low * scale), ceil_fraction(high * scale)))

    lines = [
        "R5_ZS5_D4_QUADRATIC_EXACT_INTERVAL_V1",
        f"SCALE_BITS {args.scale_bits}",
        f"REPAIR_SHA256 {hashlib.sha256(REPAIR.read_bytes()).hexdigest()}",
        f"ACTIVE_SHA256 {hashlib.sha256(ACTIVE.read_bytes()).hexdigest()}",
        f"GEOMETRY_SHA256 {hashlib.sha256(GEOMETRY.read_bytes()).hexdigest()}",
        f"CUSP_INTERVAL {lower} {upper}",
        "POWERS 4",
    ]
    lines.extend(f"{low} {high}" for low, high in powers)
    lines.append(f"COEFFICIENT_INTERVALS {len(intervals)}")
    lines.extend(f"{low} {high}" for low, high in intervals)
    lines.append(f"ORIGINAL_TO_COLLAPSED {len(mapping)}")
    # Wrap mapping lines to keep the file readable and parsing simple.
    for start in range(0, len(mapping), 30):
        lines.append(" ".join(str(q) for q in mapping[start:start + 30]))
    active_rows = [(int(item["simplex"]), tuple(item["alpha"]))
                   for item in active["active_rows"]]
    lines.append(f"ACTIVE_ROWS {len(active_rows)}")
    lines.extend(" ".join((str(sid), *(str(q) for q in alpha)))
                 for sid, alpha in sorted(active_rows))
    lines.append("END")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text("\n".join(lines) + "\n")
    print("coefficient intervals", len(intervals), "maximum scaled width", maximum_width)
    print("active rows", len(active_rows), "scale bits", args.scale_bits)
    print("wrote", args.output)


if __name__ == "__main__":
    main()
