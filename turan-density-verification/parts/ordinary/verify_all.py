"""Run the ordinary-density arithmetic checks used in the manuscript.

Requires Python 3.9+ and its standard library. No numerical optimizer is used.
Analytic reductions and the probabilistic argument are proved in the manuscript;
these scripts check finite arithmetic and read the stored book certificate.
"""
from pathlib import Path
import subprocess
import sys

def main():
    root=Path(__file__).resolve().parent
    command = [sys.executable, "-B"]
    if sys.flags.optimize:
        command.append("-" + "O" * sys.flags.optimize)
    for name in ['verify_rounding.py','verify_star37.py','verify_books.py']:
        print(name, flush=True)
        result = subprocess.run(command + [str(root/'certificates'/name)], cwd=root)
        if result.returncode:
            return result.returncode
    print('Ordinary-density arithmetic checks: PASS', flush=True)
    return 0
if __name__=='__main__': raise SystemExit(main())
