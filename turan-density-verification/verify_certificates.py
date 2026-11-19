#!/usr/bin/env python3
"""Run the exact finite arithmetic checks used in the manuscript.

Requires only the Python standard library. Reads inputs without rewriting
certificates or saving reports. Use --verbose for the individual calculations.
These checks supplement the analytic and combinatorial proofs in the text.

Run from the repository root with the manuscript's LaTeX source:
    python3 -B verify_certificates.py --main /path/to/manuscript.tex
The default source is main.tex in the repository root. Relative paths supplied
to --main are resolved from the current working directory.
"""
import argparse
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
DEFAULT_MAIN = ROOT / "main.tex"
CHECKS = (
    ("ordinary", "parts/ordinary/verify_all.py"),
    ("minimum_degree", "verify_degree_intervals.py"),
    ("vertex_pair", "parts/vertex_pair/scripts/verify_extension.py"),
)

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--main", type=Path, default=DEFAULT_MAIN,
                        help="LaTeX source (default: main.tex in the repository root)")
    parser.add_argument("--verbose", action="store_true", help="show each verifier's output")
    args = parser.parse_args()
    main_path = args.main.expanduser().resolve()
    if not main_path.is_file():
        parser.error(f"LaTeX source not found: {main_path}; supply --main PATH")
    failures = []
    command = [sys.executable, "-B"]
    if sys.flags.optimize:
        command.append("-" + "O" * sys.flags.optimize)
    for name, relative in CHECKS:
        script = ROOT / relative
        check_command = command + [str(script)]
        if name == "vertex_pair":
            check_command.extend(["--main", str(main_path)])
        result = subprocess.run(check_command, cwd=ROOT,
                                text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        print(name + (": PASS" if result.returncode == 0 else ": FAIL"), flush=True)
        if args.verbose or result.returncode:
            print(result.stdout, end="", flush=True)
        if result.returncode:
            failures.append(name)
    if not failures:
        print(f"All finite arithmetic checks for {main_path.name} passed.", flush=True)
    return int(bool(failures))

if __name__ == "__main__":
    raise SystemExit(main())
