"""Finite exact-arithmetic checks for the growth proof in the article's Appendix A.

Run with Python 3.10 or later. The geometric counts below are independent of
the finite network catalogue. The all-order bounds are proved in the paper.
"""
from collections import Counter
from fractions import Fraction
from math import ceil, floor, isqrt


def check_rectangle_intervals(T, width, height):
    """Count inverse points by their missing interval and its residual anchors."""
    D = T * (T - 1) // 2
    rises = [b * (b - 1) // 2 for b in range(T + 1)]
    starts = Counter()
    ends = Counter()
    full = points = 0
    for x in range(-T, width):
        for y in range(-D, height):
            missing = [b for b, rise in enumerate(rises)
                       if 0 <= x + b < width and 0 <= y + rise < height]
            if not missing:
                continue
            points += 1
            first, last = missing[0], missing[-1]
            assert missing == list(range(first, last + 1))
            if first == 0 and last < T:
                ends[last] += 1
            elif first > 0 and last == T:
                starts[first] += 1
            else:
                full += 1
                if first > 0:
                    ends[first - 1] -= 1
                    starts[last + 1] -= 1
    for j in range(1, T + 1):
        assert ends[j - 1] == starts[j] == height - D + (width - T + 1) * (j - 1)
    assert full == (width - T) * (height - D) + (T - 1) * D
    assert points == (width + T) * height + (width - 1) * D


def check_rectangular_factor(k, T, m):
    t, D = T + 1, T * (T - 1) // 2
    A, B = 4 * k * T**2, 6 * k * T**2
    H, rho = divmod(m - 1, A * B)
    Q = B - 2 * T + 1
    S = Fraction(A * H, H - 3 * D)
    z = ((T * B * H + A * B * D - D * (T - 1) * (S + (A * B - S) / Q))
         / ((B - 2 * T) * (H - 3 * D - Fraction(D * (T - 1), Q))))
    x = S - z
    y = (A * B - x) / Q - z
    assert x > 1 and y > 2 and T - 1 < z <= Fraction(11, 8) * T + 2

    horizontal = (T, H - 2 * D)
    vertical = (B - T, D)
    whole = (B - T, H - 2 * D)
    real = [(x, *horizontal), (y, *vertical), (z, *whole)]
    nh, nv, nf = floor(x), floor(y) - 1, floor(z)
    L = floor((z - nf) * (H - 2 * D))
    rounded = [(nh, *horizontal), (nv, *vertical), (nf, *whole)]
    if L >= D:
        rounded.append((1, B - T, L))

    def demands(rectangles):
        def proper(delta):
            return sum(n * (h - D + (w - T + 1) * delta) for n, w, h in rectangles)
        full = sum(n * ((w - T) * (h - D) + (T - 1) * D) for n, w, h in rectangles)
        return proper(0), proper(T - 1), full

    supply = (A * H, A * H + A * B * (T - 1), T * B * H + A * B * D)
    assert demands(real) == supply
    assert all(required <= available for required, available in zip(demands(rounded), supply))
    tag_box = (A + T) * (B + T) * (H + 2 * D + 1)
    R = tag_box - m - sum(n * w * h for n, w, h in rounded)
    Rstar = tag_box - m - sum(n * w * h for n, w, h in real)
    assert Rstar == (A * B + 2 * A * D * T + A * T + 2 * B * D * T
                     + B * T + 2 * D * T**2 + H * T**2 + T**2 - rho - 1)
    assert 0 < Rstar <= R <= Fraction(46, 3) * k * T**5 + 80 * k**2 * T**4
    levels = (ceil(Fraction(nh, B // (2 * T))) + ceil(Fraction(nv, (H - D) // (2 * D)))
              + nf + (L >= D))
    assert levels <= Fraction(83, 24) * T + 6
    capacity = ((A - 1) // t - levels) * ((B - 1) // t) * (H - 2 * D)
    assert capacity - R >= k**3 * T**4 * (2 * T - 421)
    primary_owners = t * sum(n * ((w + T) * h + (w - 1) * D) for n, w, h in rounded)
    assert primary_owners <= t * levels * B * H
    assert k * (1 + primary_owners + t**2 * R) <= m
    if k >= 2:
        assert k * (1 + t * levels * B * H + t**2 * R) <= m
    return L >= D


def check_planar_factor(k, T, m):
    D, B = T * (T - 1) // 2, (T - 1)**2 + 1
    p, t = 32 * k, T + 1
    w = p * t
    h, rho = divmod(m - 1, w)
    R0 = (T**3 - T) // 2 + w - 1
    h0 = max(4 * D + B, 3 * D + (R0 + D + 2 * B - 1 + p - 3) // (p - 2))
    assert h >= h0
    R = ((w + T) * (h + 1) + (w - 1) * D - m
         - (T * h + (w - 3 * T) * D + T * B))
    assert 0 < R == R0 - rho <= (p - 2) * (h - 3 * D) - (D + 2 * B - 1)
    support_bound = 2 * T * h + 2 * D * (w - T) + 2 * T * B + t * R + 1
    assert k * support_bound <= m


def main():
    for T in range(2, 9):
        D = T * (T - 1) // 2
        for width in (T, T + 1, 2 * T + 3):
            for height in (D, D + 1, 2 * D + 5):
                check_rectangle_intervals(T, width, height)
    rectangle_count = planar_count = 0
    rounding_cases = set()
    for k in (1, 2, 3, 8, 32, 225, 1024):
        for T in (*range(224, 256), 299, 503, 997, 4096):
            lower, upper = 128 * k**3 * T**7, 128 * k**3 * (T + 1)**7
            base = 24 * k**2 * T**4
            first_layer = ((lower + base - 1) // base) * base
            for m in sorted({lower, lower + 1, first_layer, first_layer + 1,
                             (lower + upper) // 2, upper - 1}):
                rounding_cases.add(check_rectangular_factor(k, T, m))
                rectangle_count += 1
    # At these consecutive heights z crosses 308. The partial rectangle
    # is respectively omitted and retained.
    for height in (60_190_275, 60_190_276):
        m = 24 * 224**4 * height + 1
        rounding_cases.add(check_rectangular_factor(1, 224, m))
        rectangle_count += 1
    assert rounding_cases == {False, True}
    for k in (1, 2, 8):
        for T in range(2, 224**2):
            lower = isqrt((128 * k**3)**2 * T**7) + 1
            w = 32 * k * (T + 1)
            next_layer = ((lower + w - 1) // w) * w
            for m in sorted({lower, next_layer, next_layer + 1}):
                check_planar_factor(k, T, m)
                planar_count += 1
    print('63 exhaustive rectangle interval counts passed.')
    print(f'{rectangle_count} rectangular and {planar_count} planar parameter checks passed.')


if __name__ == '__main__':
    main()
