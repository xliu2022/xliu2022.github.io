"""Exact constants for the 37-leaf star in the manuscript.

Checks the curvature at the rational test point and the resulting enclosure.
The analytic stationary-value estimate is proved in the manuscript.
"""
from fractions import Fraction as Q
import sys

sys.dont_write_bytecode = True
from exact_poly import Poly

c = Q(4, 9)


def require(condition, message):
    if not condition:
        raise ArithmeticError(message)


def verify():
    a = Poly.variable(0, 1)
    q, t, r = Q(8, 27), Q(36, 37), c / 37**2
    f = q*a**3 + 3*t*a*(1-a)**2 + r*(1-a)**3
    fp = f.derivative(0)
    fpp = fp.derivative(0)
    a0 = Q(3527805417, 10**10)
    require(Q(1, 4) < a0 < Q(1, 2), "Test point outside the analytic domain")
    curvature = -fpp.at((a0,))
    require(curvature == Q(6248829722647, 1283437500000) > 4,
            "Curvature at a_0 must exceed 4")

    lower = f.at((a0,))
    upper = lower + fp.at((a0,))**2 / 4
    ell = Q(73, 5550)
    require(Q(2, 10**6) < lower-c <= upper-c < Q(25, 10**7),
            "Star maximum enclosure")
    require(upper-c < 3*c*ell**3 == Q(389017, 128215406250),
            "Three-control crossing")
    require(lower-c == Q(208915215870873322169027313,
                         85562500000000000000000000000000),
            "Exact lower excess in the manuscript")
    require(upper-c == Q(6435110936862607623201324021954481082401,
                         2635538906250000000000000000000000000000000000),
            "Exact upper excess in the manuscript")
    require(Q(53, 54) > Q(73, 75) and Q(53, 3996) > ell,
            "Leaf-weight bounds")
    require(Q(75, 73)**3 < Q(1000, 729), "Leaf-weight ratio")
    require(q/64 + Q(27, 64) + r < c, "First boundary bound")
    require(Q(89, 216) + r < c, "Second boundary bound")
    print("37-leaf star curvature, enclosure and three-control crossing: PASS")
    print("Exact lower excess:", lower-c)
    print("Exact upper excess:", upper-c)


if __name__ == "__main__":
    verify()
