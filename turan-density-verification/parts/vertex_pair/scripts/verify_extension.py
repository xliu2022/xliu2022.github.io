#!/usr/bin/env python3
"""Exact certificate for the finite vertex-pair palette interval bridge.

Reads the 16 rows of the manuscript table labeled vp-tab:extension-data
and checks the extension from [0,1/4] to [0,4/9]. The palette realization
theorem and the construction of the initial interval are proved separately
in the manuscript. This script does not change the source or data files.

Each outer palette has two color classes of sizes a,b. The four entries
of k specify the numbers (0,1,2) of cyclic orbits on distinct-color triples
of types AAA, AAB, ABB, BBB, respectively. There are no repeated-color
outer triples. Each orbit contributes three ordered triples.

All calculations use fractions.Fraction. No floating-point computation,
optimization package, or third-party dependency is used.

Run from the repository root with the manuscript's LaTeX source:
    python3 -B parts/vertex_pair/scripts/verify_extension.py --main /path/to/manuscript.tex
The default source is main.tex in the repository root. Relative paths supplied
to --main are resolved from the current working directory.
"""
import argparse
import re
from dataclasses import dataclass
from fractions import Fraction as F
from math import comb
from pathlib import Path
from typing import Tuple

ROOT = Path(__file__).resolve().parents[3]
DEFAULT_MAIN = ROOT / "main.tex"


@dataclass(frozen=True)
class Step:
    a: int
    b: int
    k: Tuple[int, int, int, int]
    split: F
    target: F
    subdivisions: int
    certified_upper: F
    expected_margin: F


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ArithmeticError(message)


def read_steps(main: Path):
    """Read this table's numeric cells; reject unsupported LaTeX syntax."""
    source = re.sub(r"(?<!\\)%[^\n]*", "", main.read_text(encoding="utf-8"))
    tables = [table for table in re.findall(r"\\begin\{table\}(.*?)\\end\{table\}",
                                           source, flags=re.S)
              if r"\label{vp-tab:extension-data}" in table]
    require(len(tables) == 1, f"{main}: expected one table labeled vp-tab:extension-data")
    blocks = re.findall(r"\\midrule(.*?)\\bottomrule", tables[0], flags=re.S)
    require(len(blocks) == 1, "Expected one data block in the interval-extension table")
    rows = [row.strip() for row in blocks[0].split(r"\\") if row.strip()]
    require(len(rows) == 16, f"Expected 16 table rows, found {len(rows)}")

    def rational(cell):
        # Accept both inline a/b and the usual TeX fraction commands.
        fraction = re.fullmatch(r"\\(?:frac|tfrac|dfrac)\{(-?\d+)\}\{(\d+)\}", cell)
        if fraction:
            cell = "/".join(fraction.groups())
        require(re.fullmatch(r"-?\d+(?:/[1-9]\d*)?", cell) is not None,
                f"Expected an exact integer or fraction, found {cell!r}")
        return F(cell)

    def integer(cell):
        require(re.fullmatch(r"\d+", cell) is not None,
                f"Expected a nonnegative integer, found {cell!r}")
        return int(cell)

    def integer_tuple(cell, length):
        require(cell.startswith("(") and cell.endswith(")"),
                f"Expected a tuple, found {cell!r}")
        entries = cell[1:-1].split(",")
        require(len(entries) == length, f"Expected {length} tuple entries in {cell!r}")
        return tuple(integer(entry) for entry in entries)

    steps = []
    for index, row in enumerate(rows, 1):
        cells = [re.sub(r"\s+", "", cell).strip("$") for cell in row.split("&")]
        require(len(cells) == 8, f"Table row {index}: expected eight cells")
        require(integer(cells[0]) == index, f"Table row {index}: row index mismatch")
        a, b = integer_tuple(cells[1], 2)
        k = integer_tuple(cells[2], 4)
        steps.append(Step(a, b, k, rational(cells[3]), rational(cells[4]),
                          integer(cells[5]), rational(cells[6]), rational(cells[7])))
    return tuple(steps)


def polynomial(step: Step) -> Tuple[F, F, F, F]:
    """Power coefficients of f(x), in increasing order."""
    a, b = step.a, step.b
    k0, k1, k2, k3 = step.k
    # f(x) = A*x^3+B*x^2*(1-x)+C*x*(1-x)^2+D*(1-x)^3.
    A = F(k0*(a-1)*(a-2), 2*a*a)
    B = F(3*k1*(a-1), 2*a)
    C = F(3*k2*(b-1), 2*b)
    D = F(k3*(b-1)*(b-2), 2*b*b)
    return (D, C-3*D, B-2*C+3*D, A-B+C-D)


def value(coefficients: Tuple[F, F, F, F], x: F) -> F:
    result = F(0)
    for coefficient in reversed(coefficients):
        result = result*x + coefficient
    return result


def bernstein_on_interval(coefficients: Tuple[F, F, F, F],
                          left: F, right: F) -> Tuple[F, ...]:
    """Degree-three Bernstein coefficients after x=left+(right-left)*z.

    If all returned coefficients are <= U, then f(x)<=U on [left,right],
    because the Bernstein basis functions are nonnegative and sum to one.
    """
    length = right-left
    powers = tuple(
        sum((coefficients[j]*comb(j,k)*left**(j-k)*length**k
             for j in range(k,4)), F(0))
        for k in range(4)
    )
    return tuple(
        sum((powers[k]*F(comb(i,k),comb(3,k)) for k in range(i+1)), F(0))
        for i in range(4)
    )


def verify(steps) -> None:
    previous = F(1,4)
    print("All quantities below are exact rational numbers.\n")
    print("step | input endpoint | certified max f | output endpoint | witness margin")
    for i, step in enumerate(steps, 1):
        require(step.a > 0 and step.b > 0, f"Step {i}: invalid class size")
        require(all(k in (0,1,2) for k in step.k), f"Step {i}: invalid orbit count")
        require(F(0) < step.split < F(1), f"Step {i}: split outside (0,1)")
        require(previous < step.target < 1, f"Step {i}: endpoint does not increase")
        coeff = polynomial(step)
        n = step.subdivisions
        require(n > 0, f"Step {i}: subdivision count must be positive")
        certified = max(
            coefficient
            for j in range(n)
            for coefficient in bernstein_on_interval(coeff, F(j,n), F(j+1,n))
        )
        require(certified == step.certified_upper,
                f"Step {i}: stored Bernstein certificate is incorrect")
        require(certified <= previous,
                f"Step {i}: outer Lagrangian upper bound exceeds input endpoint")
        x = step.split
        diagonal = x**3 / step.a**2 + (1-x)**3 / step.b**2
        require(F(0) < diagonal < F(1), f"Step {i}: invalid diagonal probability")
        margin = value(coeff, x) - (1-diagonal)*step.target
        require(margin == step.expected_margin, f"Step {i}: stored witness margin is incorrect")
        require(margin >= 0, f"Step {i}: witness does not reach target fixed point")
        print(f"{i:>4} | {str(previous):>14} | {str(certified):>15} | "
              f"{str(step.target):>15} | {margin}")
        previous = step.target
    require(previous == F(4,9), "The bridge does not end at 4/9")
    print("\nPASS: all 16 finite interval-extension certificates verified exactly.")




def explicit_patterns(step: Step):
    """Build the finite relation from its cyclic tuple orbits."""
    from itertools import combinations
    patterns=set()
    for e in combinations(range(step.a+step.b), 3):
        nA=sum(v<step.a for v in e)
        orbits=step.k[3-nA]
        u,v,w=e
        if orbits>=1:
            patterns.update(((u,v,w),(v,w,u),(w,u,v)))
        if orbits==2:
            patterns.update(((u,w,v),(w,v,u),(v,u,w)))
    return patterns


def verify_relations(steps):
    """Independently reconcile each explicit relation and its polynomial."""
    from collections import Counter
    for index, step in enumerate(steps, 1):
        rel=explicit_patterns(step)
        size=step.a+step.b
        require(all(len(set(e))==3 for e in rel),
                f"Step {index}: repeated-color outer tuple")
        require(all((v,w,u) in rel for u,v,w in rel),
                f"Step {index}: relation is not cyclically invariant")
        # Construct polynomial coefficients directly by summing tuple monomials
        # for weights x/a in A and (1-x)/b in B.
        coefficients=[F(0) for _ in range(4)]
        for e in rel:
            nA=sum(v<step.a for v in e)
            nB=3-nA
            factor=F(1, step.a**nA*step.b**nB)
            for j in range(nB+1):
                coefficients[nA+j] += factor*comb(nB,j)*((-1)**j)
        require(tuple(coefficients)==polynomial(step),
                f"Step {index}: explicit relation has a different polynomial")
        # Exact witness evaluation directly from the relation.
        x=step.split
        weights=[x/step.a]*step.a+[(1-x)/step.b]*step.b
        raw=sum((weights[u]*weights[v]*weights[w] for u,v,w in rel), F(0))
        require(raw==value(polynomial(step),x),
                f"Step {index}: explicit witness mismatch")
        # Cyclic invariance gives the coordinate derivative identity used in
        # the substitution lemma, independently for every color.
        degree=[Counter(e[i] for e in rel) for i in range(3)]
        require(degree[0]==degree[1]==degree[2],
                f"Step {index}: coordinate marginal counts disagree")
        print(f"relation {index:2}: {size:2} colors, {len(rel):4} tuples; PASS")
    print("\nPASS: all 16 explicit relations agree with the certified polynomials.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--main", type=Path, default=DEFAULT_MAIN,
                        help="LaTeX source containing the interval-extension table (default: main.tex)")
    args = parser.parse_args()
    args.main = args.main.expanduser().resolve()
    if not args.main.is_file():
        parser.error(f"LaTeX source not found: {args.main}; supply --main PATH")
    steps = read_steps(args.main)
    print("Reading interval-extension data from:", args.main)
    verify(steps)
    print()
    verify_relations(steps)
