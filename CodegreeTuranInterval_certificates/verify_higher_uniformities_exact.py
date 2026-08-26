#!/usr/bin/env python3
"""Exact verifier for the rational certificates in the k=5 and k=6 proofs.

The script uses only the Python standard library.  It checks the Bernstein
coefficients and endpoint gaps in the four k=5 compiler rows, as well as the
safety values, image endpoints, and overlaps in the two k=6 compiler branches.
"""

from __future__ import annotations

from fractions import Fraction as F
from math import comb


class CertificateError(RuntimeError):
    """Raised when an exact rational check fails."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise CertificateError(message)


def bernstein_coefficients(theta: F, g: F, left: F, right: F) -> tuple[F, ...]:
    """Return the degree-three Bernstein coefficients on [left,right]."""
    power = (3 * g, 8 * theta - 9 * g, 3 * g - 9 * theta, 3 * theta)
    width = right - left
    shifted = (
        power[0] + power[1] * left + power[2] * left**2 + power[3] * left**3,
        width * (power[1] + 2 * power[2] * left + 3 * power[3] * left**2),
        width**2 * (power[2] + 3 * power[3] * left),
        width**3 * power[3],
    )
    return tuple(
        sum(
            (shifted[j] * F(comb(i, j), comb(3, j)) for j in range(i + 1)),
            F(0),
        )
        for i in range(4)
    )


def verify_k5() -> None:
    rows = (
        (
            F(8, 27),
            F(3, 32),
            (
                (F(0), F(3, 2), (F(9, 32), F(1805, 1728), F(65, 3456), F(235, 1152))),
                (F(3, 2), F(3), (F(235, 1152), F(1345, 3456), F(3085, 1728), F(2129, 288))),
            ),
            F(2, 5),
            F(3, 2),
            F(1, 9),
            F(3, 20),
            (F(449, 10800), F(29, 9261), F(1, 200)),
        ),
        (
            F(12, 25),
            F(3, 32),
            (
                (F(0), F(3, 2), (F(9, 32), F(2847, 1600), F(159, 640), F(351, 640))),
                (F(3, 2), F(3), (F(351, 640), F(543, 640), F(4767, 1600), F(9441, 800))),
            ),
            F(12, 25),
            F(5, 3),
            F(3, 20),
            F(47, 200),
            (F(1653, 125000), F(1959, 1013060), F(269, 409600)),
        ),
        (
            F(240, 289),
            F(3, 32),
            (
                (F(0), F(3, 2), (F(9, 32), F(58839, 18496), F(25395, 36992), F(44595, 36992))),
                (F(3, 2), F(3), (F(44595, 36992), F(63795, 36992), F(97239, 18496), F(186921, 9248))),
            ),
            F(89, 170),
            F(9, 5),
            F(47, 200),
            F(3, 8),
            (F(9903, 26726720), F(577013, 3474795800), F(400233, 25376512)),
        ),
        (
            F(1),
            F(120, 343),
            (
                (F(0), F(3, 4), (F(360, 343), F(776, 343), F(10891, 5488), F(32523, 21952))),
                (F(3, 4), F(3, 2), (F(32523, 21952), F(10741, 10976), F(1375, 5488), F(1545, 2744))),
                (F(3, 2), F(3), (F(1545, 2744), F(815, 686), F(4099, 686), F(8592, 343))),
            ),
            F(2, 5),
            F(41, 43),
            F(3, 8),
            F(1, 2),
            (F(859, 8575), F(29, 2744), F(1, 148176)),
        ),
    )

    for number, (theta, g, intervals, r_minus, r_plus, left, right, gaps) in enumerate(rows, 1):
        for u, v, expected in intervals:
            actual = bernstein_coefficients(theta, g, u, v)
            require(actual == expected, f"k=5 row {number}: Bernstein coefficients fail on [{u},{v}]")
            require(all(value > 0 for value in actual), f"k=5 row {number}: a Bernstein coefficient is not positive")

        def phi(x: F, r: F) -> F:
            return (theta + 3 * x * r**2) / (1 + r) ** 3

        def psi(r: F) -> F:
            return (3 * g * r + 3 * theta * r**2) / (1 + r) ** 3

        actual_gaps = (
            theta - (3 * g * r_minus + 3 * theta * r_minus**2),
            left - phi(F(0), r_minus),
            min(phi(F(1), r_plus), psi(r_plus)) - right,
        )
        require(actual_gaps == gaps, f"k=5 row {number}: endpoint-gap certificate fails")
        require(all(value > 0 for value in actual_gaps), f"k=5 row {number}: an endpoint gap is not positive")

    print("[PASS] k=5: four compiler rows, Bernstein coefficients, and endpoint gaps")


def verify_k6() -> None:
    def target(theta: F, r: F) -> F:
        return theta / ((1 - r) * (1 + r) ** 3)

    def mixed_one(theta: F, r: F) -> F:
        return 8 * r**2 * theta - 9 * r * theta + 3 * r + 3 * theta - 3

    def mixed_zero(theta: F, r: F) -> F:
        return 17 * r**3 * theta - 18 * r**2 * theta + 12 * r**2 + 3 * r * theta - 12 * r + 3 * theta

    branches = (
        (
            F(3, 32),
            F(6, 7),
            F(18, 19),
            F(7203, 70304),
            F(390963, 1620896),
            F(-33, 11552),
            F(-46107, 219488),
        ),
        (
            F(24, 125),
            F(6, 7),
            F(9, 10),
            F(57624, 274625),
            F(1920, 6859),
            F(-219, 6250),
            F(-6336, 15625),
        ),
    )
    for number, (theta, left, right, image_left, image_right, safety_one, safety_zero) in enumerate(branches, 1):
        require(target(theta, left) == image_left, f"k=6 branch {number}: left image endpoint fails")
        require(target(theta, right) == image_right, f"k=6 branch {number}: right image endpoint fails")
        require(mixed_one(theta, right) == safety_one <= 0, f"k=6 branch {number}: first safety check fails")
        require(mixed_zero(theta, right) == safety_zero <= 0, f"k=6 branch {number}: second safety check fails")

    require(F(3, 16) - F(7203, 70304) > 0, "k=6: first bridge does not meet the source interval")
    require(
        F(390963, 1620896) - F(57624, 274625) == F(13965702771, 445138564000) > 0,
        "k=6: the two bridge branches do not overlap",
    )
    require(F(1920, 6859) - F(5, 18) == F(265, 123462) > 0, "k=6: second bridge misses the terminal threshold")
    print("[PASS] k=6: two compiler branches, safety values, endpoints, and overlaps")


def main() -> int:
    try:
        verify_k5()
        verify_k6()
    except CertificateError as exc:
        print(f"[FAIL] {exc}")
        return 1
    print("ALL HIGHER-UNIFORMITY EXACT CERTIFICATES PASSED")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
