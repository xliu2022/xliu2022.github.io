"""Exact checks for all 132 book test points and interval bounds in the manuscript.

At each test point x, check the leading minors of -Hessian(P_D)(x)-(2/3)I.
The stationary-value lemma in the manuscript converts these pointwise bounds
into the upper bound U = P_D(x) + ||gradient(P_D)(x)||^2.

The default run recomputes and compares book_certificate.json without writing
files. Use --write-certificate explicitly to regenerate it. Standard library
only; the analytic reduction is not proved by this finite arithmetic check.
"""
import argparse
import json
import sys
from fractions import Fraction as Q
from pathlib import Path

sys.dont_write_bytecode = True
from exact_poly import book_polynomial, det3

ROOT = Path(__file__).resolve().parent
c = Q(4, 9)


def require(condition, message):
    if not condition:
        raise ArithmeticError(message)


def beta(m):
    require(m is None or (type(m) is int and m >= 1),
            "A star index must be a positive integer or null (infinity)")
    return c if m is None else c*(1-Q(1, m))


def verify(windows):
    rows = json.loads(windows.read_text(encoding="utf-8"))["rows"]
    require(isinstance(rows, list) and len(rows) == 132, "Expected 132 book rows")
    minima = [None] * 3
    intervals, records = [], []
    thresholds = (Q(6), Q(2), Q(2))
    h = Q(2, 3)
    require(Q(2, 3) / h == 1, "Stationary-value error coefficient")

    for index, row in enumerate(rows):
        context = f"Book row {index}"
        require(type(row["index"]) is int and row["index"] == index,
                f"{context}: row indices must be 0,...,131")
        D, den = row["D"], row["denominator"]
        numerators = row["numerators"]
        require(type(D) is int and D >= 10, f"{context}: D must be at least 10")
        require(type(den) is int and den == 10**8, f"{context}: denominator must be 10^8")
        require(isinstance(numerators, list) and len(numerators) == 4
                and all(type(v) is int for v in numerators),
                f"{context}: expected four integer numerators")
        a, u, v, z = (Q(v, den) for v in numerators)
        require(a+u+v+z == 1 and min(a, u, v, z) > 0,
                f"{context}: weights must be positive and sum to one")
        x = (a, u+v, u-v)
        q1, q2 = beta(row["m"]), beta(row["n"])
        P = book_polynomial(D, q1, q2)
        derivatives = [P.derivative(i) for i in range(3)]
        H = [[derivatives[i].derivative(j).at(x) for j in range(3)]
             for i in range(3)]
        require(all(H[i][j] == H[j][i] for i in range(3) for j in range(3)),
                f"{context}: Hessian is not symmetric")
        matrix = [[-H[i][j] - (h if i == j else 0) for j in range(3)]
                  for i in range(3)]
        minors = [matrix[0][0],
                  matrix[0][0]*matrix[1][1] - matrix[0][1]*matrix[1][0],
                  det3(matrix)]
        for j, (minor, threshold) in enumerate(zip(minors, thresholds)):
            require(minor > threshold,
                    f"{context}: leading minor {j+1} is {minor}, not > {threshold}")
            minima[j] = minor if minima[j] is None else min(minima[j], minor)

        lower = P.at(x)
        # Independently evaluate the original four-weight expression.
        raw = (3*a*((1-a)**2-u*u-v*v-z*z/D) + 6*u*v*z
               + q1*u**3 + q2*v**3 + c*z**3/D**2)
        require(lower == raw, f"{context}: coordinate substitution mismatch")
        gradient = [derivative.at(x) for derivative in derivatives]
        upper = lower + sum(g*g for g in gradient)
        width = Q(3*(D-7), 250*D**3)
        require(width == (D-7)*c*Q(3, 10*D)**3, f"{context}: interval width")
        require(lower > c, f"{context}: L must exceed 4/9")
        require(upper-width <= lower, f"{context}: empty certified interval")
        intervals.append((upper-width, lower, index))
        records.append({
            "index": index, "L": str(lower), "U": str(upper),
            "width": str(width), "left": str(upper-width), "right": str(lower),
            "gradient": [str(g) for g in gradient],
            "hessian_leading_minors": [str(v) for v in minors],
        })

    intervals.sort()
    require(intervals[0][0] < c, "The cover must start below 4/9")
    right, overlaps = c, []
    for left, end, index in intervals:
        require(left <= right, f"Book row {index}: gap in the interval cover")
        if end > right:
            overlaps.append(right-left)
            right = end
    require(min(overlaps) > Q(2, 10**7), "Extending overlap must exceed 2*10^(-7)")
    require(right-Q(162, 361) > Q(1, 100000), "Final excess must exceed 10^(-5)")
    return {
        "hessian_point_checks": len(rows), "hessian_lower_bound": str(h),
        "leading_minor_thresholds": [str(v) for v in thresholds],
        "leading_minor_minima": [str(v) for v in minima],
        "interval_count": len(rows), "minimum_extending_overlap": str(min(overlaps)),
        "final_right_endpoint": str(right),
        "final_excess_over_162_361": str(right-Q(162, 361)), "records": records,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--windows", type=Path, default=ROOT / "book_windows.json")
    parser.add_argument("--certificate", type=Path, default=ROOT / "book_certificate.json")
    parser.add_argument("--write-certificate", action="store_true",
                        help="replace the certificate only after all mathematical checks pass")
    args = parser.parse_args()
    require(args.windows.resolve() != args.certificate.resolve(),
            "Input rows and output certificate must be different files")
    result = verify(args.windows)
    if args.write_certificate:
        args.certificate.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
        print("Wrote verified certificate:", args.certificate)
    else:
        saved = json.loads(args.certificate.read_text(encoding="utf-8"))
        require(saved == result,
                "Stored book certificate differs from the exact recomputation; "
                "inspect the data before using --write-certificate")
        print("Stored book certificate matches the exact recomputation: PASS")
    print("132 pointwise Hessian checks with lower bound 2/3: PASS")
    print("Leading principal minor minima:", ", ".join(result["leading_minor_minima"]))
    print("132 interval enclosures and exact cover of [4/9,162/361]: PASS")
    print("Smallest extending overlap:", result["minimum_extending_overlap"])
    print("Final excess:", result["final_excess_over_162_361"])


if __name__ == "__main__":
    main()
