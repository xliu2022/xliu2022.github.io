#!/usr/bin/env python3
"""Extract and classify near-active rows of the frozen d=4 quadratic candidate."""
from __future__ import annotations

import argparse
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(HERE))
import r5_zs5_d4_projected_quadratic_cutting_probe as common


DEFAULT_CANDIDATE = (ROOT / "work/emc_audit/computer/"
                     "r5_zs5_d4_projected_quadratic_symmetric_zero_highs_probe.json")
DEFAULT_OUTPUT = (ROOT / "work/emc_audit/computer/"
                  "r5_zs5_d4_projected_quadratic_active_audit.json")


def simplex_boxes():
    geometry = json.loads(common.MANIFEST.read_text())
    answer = []
    for cell in geometry["cells"]:
        b_box = int(cell["key"][0])
        box = tuple(int(q) for q in cell["key"][1:5])
        answer.extend((b_box, box) for _ in cell["simplices"])
    return answer


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidate", type=Path, default=DEFAULT_CANDIDATE)
    parser.add_argument("--executable", type=Path,
                        default=Path("/private/tmp/r5_zs5_d4_projected_quadratic_scan_elevated"))
    parser.add_argument("--top", type=int, default=20_000)
    parser.add_argument("--tolerance", type=float, default=1e-9)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    data = json.loads(args.candidate.read_text())
    coefficients = np.asarray([q for item in data["g"] for q in item["coefficients"]])
    coefficient_path = Path("/private/tmp/r5_zs5_d4_quadratic_active_candidate.txt")
    common.write_candidate(coefficient_path, coefficients)
    summary, rows = common.scan(args.executable, coefficient_path, args.top, 4)
    boxes = simplex_boxes()
    active = [row for row in rows if abs(row[0]) < args.tolerance]
    by_b_box = Counter(boxes[meta[0]][0] for _, _, _, meta in active)
    by_alpha = Counter(str(meta[1]) for _, _, _, meta in active)
    by_full_box = Counter(str(boxes[meta[0]]) for _, _, _, meta in active)
    touched_keys = defaultdict(set)
    for _, _, row, meta in active:
        b_box = boxes[meta[0]][0]
        touched_keys[b_box].update(column // common.BASIS for column in row)
    payload = {
        "proof_status": "FLOATING ACTIVE-ROW AUDIT ONLY",
        "candidate": str(args.candidate),
        "scan_summary": summary,
        "active_tolerance": args.tolerance,
        "active_count": len(active),
        "active_by_b_box": dict(sorted(by_b_box.items())),
        "distinct_active_alpha_count": len(by_alpha),
        "active_by_alpha": dict(by_alpha.most_common()),
        "distinct_active_full_box_count": len(by_full_box),
        "most_common_active_full_boxes": by_full_box.most_common(30),
        "touched_projected_keys_by_b_box": {
            str(key): len(value) for key, value in sorted(touched_keys.items())},
        "active_rows": [{
            "gap": raw,
            "base_at_float_cusp": base,
            "simplex": meta[0],
            "alpha": list(meta[1]),
        } for raw, base, _, meta in active],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps({key: payload[key] for key in (
        "active_count", "active_by_b_box", "distinct_active_alpha_count",
        "distinct_active_full_box_count", "touched_projected_keys_by_b_box")},
        indent=2))
    print("wrote", args.output)


if __name__ == "__main__":
    main()
