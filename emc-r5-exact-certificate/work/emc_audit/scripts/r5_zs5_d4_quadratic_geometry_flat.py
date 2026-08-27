#!/usr/bin/env python3
"""Flatten the exact d=4 orbit geometry for the quadratic C++ scanner."""
from __future__ import annotations

import argparse
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
DEFAULT_GEOMETRY = (ROOT / "work/emc_audit/computer/"
                    "r5_zs5_d4_orbit_geometry_exact.json")
DEFAULT_OUTPUT = (ROOT / "work/emc_audit/computer/"
                  "r5_zs5_d4_quadratic_geometry_exact_flat.txt")


def decode_key(key):
    return int(key[0]), tuple(int(q) for q in key[1:5])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--geometry", type=Path, default=DEFAULT_GEOMETRY)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    data = json.loads(args.geometry.read_text())
    key_where = {key: index for index, key in enumerate(data["projected_keys"])}
    lines = [f"VERTICES {len(data['vertices'])}"]
    lines.extend(" ".join(point) for point in data["vertices"])
    simplex_count = sum(len(cell["simplices"]) for cell in data["cells"])
    lines.append(f"PROJECTED_KEYS {len(key_where)}")
    lines.extend(data["projected_keys"])
    lines.append(f"SIMPLICES {simplex_count}")
    for cell in data["cells"]:
        b_box, box = decode_key(cell["key"])
        prefix = [str(b_box), *(str(q) for q in box)]
        for projection in cell["projections"]:
            prefix.extend((str(key_where[projection["key"]]),
                           *(str(q) for q in projection["axes"])))
        for simplex in cell["simplices"]:
            lines.append(" ".join((*prefix, *(str(q) for q in simplex))))
    lines.append("END")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text("\n".join(lines) + "\n")
    print("vertices", len(data["vertices"]), "projected keys", len(key_where),
          "simplices", simplex_count)
    print("wrote", args.output)


if __name__ == "__main__":
    main()
