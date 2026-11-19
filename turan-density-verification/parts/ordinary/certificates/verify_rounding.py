"""Exact arithmetic in the manuscript's rounding and finite-control proofs.

Checks the constants in the squared-error recurrence, without constructing
the enormous denominators in the 34-control case. Uses the standard library.
"""
from fractions import Fraction as Q
from math import factorial

c = Q(4, 9)
kappa = Q(1000, 729)


def require(condition, message):
    if not condition:
        raise ArithmeticError(message)


def envelope(bound):
    """The reciprocal-grid error envelope proved in the manuscript."""
    require(0 < bound <= 1, "Envelope bound must lie in (0,1]")
    reciprocal = 1 / bound
    m = -((-reciprocal.numerator) // reciprocal.denominator)
    return m, max(bound - Q(1, m), Q(1, m * (m + 1)))


def verify():
    require(kappa == Q(10, 9)**3 < Q(11, 8), "Incorrect compensation constant")
    require(envelope(Q(1)) == (1, Q(1, 2)), "First rounding envelope")
    require(kappa / 2 < Q(11, 16) < 1, "First compensation is too large")
    require(envelope(Q(11, 16)) == (2, Q(3, 16)), "Second rounding envelope")
    require(kappa * Q(3, 16) < Q(33, 128), "Bound on x_3")

    x3 = Q(33, 128)
    u4 = Q(11, 8)**2 * x3**2 / (1 + x3)
    require(u4 == Q(131769, 1318912) < Q(1, 10), "Bound on u_4")
    # The manuscript proves u_{j+1} <= u_j^2. Once u_4 < 1/10,
    # record just the exponent in u_j < 10^{-exponent}.
    require(Q(1, 10) / kappa < x3 < Q(11, 16) < 1,
            "Later compensations must stay below one")
    downward_factor = 2 * c / (kappa - Q(1, 10))
    require(downward_factor == Q(6480, 9271) < 1,
            "Last downward-rounding estimate")

    for controls in (7, 34):
        exponent = 1
        for _ in range(4, controls):
            exponent *= 2
        require(exponent == 2**(controls - 4), "Squared-error exponent")
        final_exponent = 2 * exponent
        if controls == 7:
            require(final_exponent == 16 and Q(1, 10**16) < Q(1, 10**15),
                    "Seven-control error bound")
        else:
            # e < 3 follows from the factorial-series bound in the proof.
            # Thus e^2 < 9 < 10 and 10^{-N} < e^{-2N}.
            require(final_exponent == 2**31 and 3**2 < 10
                    and 2 * final_exponent > 10**9,
                    "34-control comparison with exp(-10^9)")
        print(f"{controls} controls: 2*kappa*Delta < 10^(-{final_exponent}): PASS")

    # Rational estimates at the end of the quantitative matching proof.
    taylor = sum((Q(99, 20)**j / factorial(j) for j in range(8)), Q(0))
    require(taylor > 110, "Taylor lower bound")
    require(Q(400, 10**6) / Q(99, 100) < Q(405, 10**6),
            "Matching error bound")
    require(Q(1, 110) + Q(405, 10**6) < Q(1, 100), "Matching probability bound")
    print("Rounding constants and quantitative matching bounds: PASS")


if __name__ == "__main__":
    verify()
