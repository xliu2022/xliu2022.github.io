#!/usr/bin/env python3
"""Exact, standalone verifier for the pi_2(4) certificate tables.

The verifier uses only the Python standard library and the three accompanying
text tables:

* weighted_path_chain_certificate.txt;
* heterogeneous_satellite_chain_certificate.txt;
* pi2_4_high_satellite_rows.txt.

Put this file beside those tables and run it with no arguments, or pass their
directory explicitly with ``--data-dir``.  Every numerical operation is over
``fractions.Fraction``.  For each low heterogeneous row, the verifier enumerates
every labeled copy in the graph family named by the table, selects a canonical
witness, and reports when the witness is not unique.
"""

from __future__ import annotations

import argparse
import ast
import itertools
import re
import sys
from dataclasses import dataclass
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path
from typing import Iterable, Sequence


class CertificateError(RuntimeError):
    """Raised for a failed exact certificate check."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise CertificateError(message)


def Q(token: str) -> F:
    token = token.strip()
    if token.lower() in {"infinity", "inf", "none"}:
        raise CertificateError(f"expected a rational number, found {token!r}")
    try:
        return F(token)
    except (ValueError, ZeroDivisionError) as exc:
        raise CertificateError(f"invalid rational number {token!r}") from exc


def alpha(s: int | None) -> F:
    return F(1) if s is None else F(s, s - 1)


def kappa(rho: F) -> F:
    return 4 * (1 - rho) / (2 - rho) ** 2


Edge = tuple[int, int]
Edges = tuple[Edge, ...]


def normal_edges(edges: Iterable[Sequence[int]]) -> Edges:
    result = tuple(sorted((min(int(i), int(j)), max(int(i), int(j))) for i, j in edges))
    require(all(i != j for i, j in result), "a graph contains a loop")
    require(len(set(result)) == len(result), "a graph contains a repeated edge")
    return result


def graph_stats(m: int, edges: Edges, context: str) -> tuple[int, ...]:
    require(m >= 3, f"{context}: m must be at least 3")
    degree = [0] * m
    edge_set = set(edges)
    require(len(edge_set) == len(edges), f"{context}: duplicate edge")
    for i, j in edges:
        require(0 <= i < j < m, f"{context}: invalid edge {(i, j)}")
        degree[i] += 1
        degree[j] += 1
    require(min(degree) >= 1, f"{context}: isolated vertex")
    require(max(degree) <= 2, f"{context}: degree exceeds 2")
    for i, j, h in itertools.combinations(range(m), 3):
        triangle = {(i, j), (i, h), (j, h)}
        require(not triangle <= edge_set, f"{context}: triangle on {(i, j, h)}")
    return tuple(degree)


def satellite_data(s_values: Sequence[int | None], edges: Edges) -> tuple[F, F, F, F, F, F, F]:
    weights = tuple(alpha(s) for s in s_values)
    S = sum(weights, F(0))
    K = sum((weights[i] * weights[j] for i, j in edges), F(0))
    neighbour_sums = [F(0)] * len(weights)
    for i, j in edges:
        neighbour_sums[i] += weights[j]
        neighbour_sums[j] += weights[i]
    L = max(neighbour_sums)
    rho = 1 - K / S
    T = 1 / (2 * S)
    B = T * kappa(rho)
    xstar = 4 * rho * T
    return S, K, L, rho, T, B, xstar


# ---------------------------------------------------------------------------
# Weighted path table
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class WeightedRow:
    step: int | None
    left: F
    target: F
    m: int
    s: int | None
    c: int
    rho: F
    B: F
    xstar: F


def parse_weighted(path: Path) -> tuple[list[WeightedRow], list[WeightedRow]]:
    special: list[WeightedRow] = []
    chain: list[WeightedRow] = []
    for line_number, raw in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        if "|" not in raw:
            continue
        parts = [part.strip() for part in raw.split("|")]
        try:
            if len(parts) == 8 and "/" in parts[0]:
                left, target = Q(parts[0]), Q(parts[1])
                s = None if parts[3].lower() == "infinity" else int(parts[3])
                special.append(
                    WeightedRow(None, left, target, int(parts[2]), s, int(parts[4]),
                                Q(parts[5]), Q(parts[6]), Q(parts[7]))
                )
            elif len(parts) == 9 and parts[0].isdigit():
                s = None if parts[4].lower() == "infinity" else int(parts[4])
                chain.append(
                    WeightedRow(int(parts[0]), Q(parts[1]), Q(parts[2]), int(parts[3]),
                                s, int(parts[5]), Q(parts[6]), Q(parts[7]), Q(parts[8]))
                )
        except (ValueError, CertificateError) as exc:
            raise CertificateError(f"{path.name}:{line_number}: cannot parse row") from exc
    require(len(special) == 4, f"{path.name}: expected 4 special rows, found {len(special)}")
    require(len(chain) == 109, f"{path.name}: expected 109 chain rows, found {len(chain)}")
    return special, chain


def verify_weighted(path: Path) -> None:
    special, chain = parse_weighted(path)
    all_rows = special + chain
    for position, row in enumerate(all_rows, 1):
        context = f"weighted row {position}"
        require(8 <= row.m <= 55, f"{context}: m outside the homogeneous box")
        require(1 <= row.c <= row.m // 2,
                f"{context}: no union of {row.c} nontrivial paths on {row.m} vertices")
        if row.s is not None:
            require(row.s >= 2, f"{context}: s<2")
        gamma = F(1) if row.s is None else F(row.s - 1, row.s)
        edges = row.m - row.c
        rho = (gamma * row.m - edges) / (gamma * row.m)
        target = gamma / (2 * row.m)
        B = target * kappa(rho)
        xstar = 4 * rho * target
        require(rho == row.rho, f"{context}: rho mismatch: {rho} != {row.rho}")
        require(target == row.target, f"{context}: T mismatch: {target} != {row.target}")
        require(B == row.B, f"{context}: B mismatch: {B} != {row.B}")
        require(xstar == row.xstar, f"{context}: x_* mismatch: {xstar} != {row.xstar}")
        require(F(4, 5) <= gamma <= 1, f"{context}: gamma outside the stated box")
        require(F(1, 7) <= rho <= F(1, 4), f"{context}: rho outside the stated box")
        require(B <= row.left and xstar <= row.left < target,
                f"{context}: routing inequalities fail")
    for previous, current in zip(all_rows, all_rows[1:]):
        require(previous.target == current.left,
                f"weighted route breaks at {previous.target} -> {current.left}")
    require(all(row.step == i for i, row in enumerate(chain, 1)),
            "weighted chain step labels are not 1,...,109")
    require(all_rows[0].left == F(1, 114), "weighted route has wrong initial input")
    require(all_rows[-1].target == F(1, 20), "weighted route has wrong endpoint")

    # Startup reciprocal certificate.  The integer loop computes the ceiling
    # without floating point.  The parity margins and their monotonicity are
    # checked by exact polynomial evaluations.
    for h in range(57, 64):
        q = 0
        while (q + 2) ** 2 < 4 * (h + 1):
            q += 1
        require((q + 1) ** 2 < 4 * (h + 1) <= (q + 2) ** 2,
                f"startup h={h}: ceiling computation failed")
        require(q * q + 4 * q >= 4 * h, f"startup h={h}: first inequality fails")
        require(4 * q * (h + 1) <= h * h, f"startup h={h}: second inequality fails")

    odd = lambda r: r**4 - 4*r**3 - 14*r**2 - 20*r - 7
    odd_d1 = lambda r: 4*r**3 - 12*r**2 - 28*r - 20
    odd_d2 = lambda r: 12*r**2 - 24*r - 28
    require(odd(7) > 0 and odd_d1(7) > 0 and odd_d2(7) > 0,
            "odd startup margin/monotonicity fails at r=7")
    # odd_d2 is increasing for r>=7 because its derivative is 24r-24>0.
    require(24*7 - 24 > 0, "odd startup second derivative is not increasing")

    even = lambda r: r * (r**3 - 6*r**2 - 7*r - 8)
    even_d1 = lambda r: 4*r**3 - 18*r**2 - 14*r - 8
    even_d2 = lambda r: 12*r**2 - 36*r - 14
    require(even(8) > 0 and even_d1(8) > 0 and even_d2(8) > 0,
            "even startup margin/monotonicity fails at r=8")
    require(24*8 - 36 > 0, "even startup second derivative is not increasing")
    print("[PASS] weighted path: 4 special + 109 chain rows; exact route 1/114 -> 1/20")


# ---------------------------------------------------------------------------
# Low heterogeneous table, including finite reconstruction of missing graphs
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class HeterogeneousRow:
    step: int
    family: str
    left: F
    target: F
    s_values: tuple[int | None, ...]
    rho: F
    B: F
    xstar: F


def parse_s_vector(text: str) -> tuple[int | None, ...]:
    values: list[int | None] = []
    for token in text.split(","):
        token = token.strip()
        values.append(None if token.lower() in {"none", "inf", "infinity"} else int(token))
    return tuple(values)


def parse_heterogeneous(path: Path) -> list[HeterogeneousRow]:
    rows: list[HeterogeneousRow] = []
    for line_number, raw in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        if "|" not in raw:
            continue
        parts = [part.strip() for part in raw.split("|")]
        if len(parts) != 8 or not parts[0].isdigit():
            continue
        try:
            rows.append(
                HeterogeneousRow(int(parts[0]), parts[1], Q(parts[2]), Q(parts[3]),
                                 parse_s_vector(parts[4]), Q(parts[5]), Q(parts[6]), Q(parts[7]))
            )
        except (ValueError, CertificateError) as exc:
            raise CertificateError(f"{path.name}:{line_number}: cannot parse row") from exc
    require(len(rows) == 57, f"{path.name}: expected 57 rows, found {len(rows)}")
    require(all(row.step == i for i, row in enumerate(rows, 1)),
            "heterogeneous step labels are not 1,...,57")
    return rows


def path_forest_base(component_sizes: tuple[int, ...]) -> tuple[int, Edges]:
    cursor = 0
    edges: list[Edge] = []
    for size in component_sizes:
        require(size >= 2, "path component must be nontrivial")
        edges.extend((cursor + i, cursor + i + 1) for i in range(size - 1))
        cursor += size
    return cursor, normal_edges(edges)


@lru_cache(maxsize=None)
def labeled_copies(m: int, base: Edges) -> tuple[Edges, ...]:
    copies: set[Edges] = set()
    for permutation in itertools.permutations(range(m)):
        copies.add(normal_edges((permutation[i], permutation[j]) for i, j in base))
    return tuple(sorted(copies))


@lru_cache(maxsize=None)
def family_graphs(family: str) -> tuple[Edges, ...]:
    match = re.search(r"m=(\d+)", family)
    require(match is not None, f"cannot read m from family {family!r}")
    m = int(match.group(1))
    bases: list[Edges] = []
    if "P4+P2+P2" in family:
        m0, base = path_forest_base((4, 2, 2))
        require(m0 == m, f"family {family!r} has inconsistent m")
        bases.append(base)
    elif "P6+P2" in family:
        m0, base = path_forest_base((6, 2))
        require(m0 == m, f"family {family!r} has inconsistent m")
        bases.append(base)
    elif "final C4+P2" in family:
        require(m == 6, "C4+P2 family must have m=6")
        bases.append(normal_edges(((0, 1), (1, 2), (2, 3), (0, 3), (4, 5))))
    elif "path forest/full path" in family and m == 7:
        for signature in ((7,), (5, 2), (4, 3), (3, 2, 2)):
            _, base = path_forest_base(signature)
            bases.append(base)
    elif "path forest/full path" in family and m == 6:
        for signature in ((6,), (4, 2), (3, 3), (2, 2, 2)):
            _, base = path_forest_base(signature)
            bases.append(base)
    else:
        raise CertificateError(f"unsupported heterogeneous family {family!r}")
    graphs: set[Edges] = set()
    for base in bases:
        graphs.update(labeled_copies(m, base))
    return tuple(sorted(graphs))


@dataclass(frozen=True)
class RecoveredWitness:
    edges: Edges
    S: F
    K: F
    L: F
    matches: int
    distinct_L: int


def preferred_low_template(step: int) -> Edges:
    """Natural fixed templates recovered from the rowwise state changes.

    These seven canonical templates reproduce all 57 rows.  The exhaustive
    family enumeration below is the authoritative existence check and also
    measures label-level non-uniqueness.
    """
    if 1 <= step <= 8:
        return normal_edges(((0, 1), (1, 2), (2, 3), (4, 5), (6, 7)))
    if 9 <= step <= 17:
        return normal_edges(((0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (6, 7)))
    if 18 <= step <= 24:
        return normal_edges(((0, 1), (1, 2), (3, 4), (4, 5), (5, 6)))
    if 25 <= step <= 36:
        return normal_edges((i, i + 1) for i in range(6))
    if 37 <= step <= 41:
        return normal_edges(((0, 1), (1, 2), (3, 4), (4, 5)))
    if 42 <= step <= 56:
        return normal_edges((i, i + 1) for i in range(5))
    if step == 57:
        return normal_edges(((1, 4), (0, 2), (2, 3), (3, 5), (0, 5)))
    raise CertificateError(f"no preferred low template for row {step}")


def recover_heterogeneous_witness(row: HeterogeneousRow) -> RecoveredWitness:
    m = len(row.s_values)
    weights = tuple(alpha(s) for s in row.s_values)
    S = sum(weights, F(0))
    expected_K = S * (1 - row.rho)
    matches: list[tuple[Edges, F]] = []
    for edges in family_graphs(row.family):
        K = sum((weights[i] * weights[j] for i, j in edges), F(0))
        if K != expected_K:
            continue
        neighbour_sums = [F(0)] * m
        for i, j in edges:
            neighbour_sums[i] += weights[j]
            neighbour_sums[j] += weights[i]
        L = max(neighbour_sums)
        if L <= F(5, 2):
            matches.append((edges, L))
    require(matches, f"heterogeneous row {row.step}: no graph witness in {row.family}")
    matches.sort()
    edges = preferred_low_template(row.step)
    match_map = dict(matches)
    require(edges in match_map,
            f"heterogeneous row {row.step}: recovered natural template is not a valid witness")
    L = match_map[edges]
    graph_stats(m, edges, f"heterogeneous row {row.step}")
    return RecoveredWitness(edges, S, expected_K, L, len(matches), len({item[1] for item in matches}))


def check_low_tails(rows: Sequence[HeterogeneousRow]) -> None:
    # P6+P2 homogeneous tail: n>=12, with the first row attached to finite row 17.
    n = 12
    T = F(n - 1, 16 * n)
    B = F(3 * (n - 1) ** 2, (7*n - 4) ** 2)
    xstar = F(n - 4, 16*n)
    M = rows[16].target
    require(B <= M and xstar <= M < T, "P6+P2 tail does not attach at n=12")
    # For n>=12, T_{n-1}-B_n has positive numerator
    # (n+2)(n^2-12n+8); the quadratic is positive and increasing there.
    p6_left = poly_add(
        poly_mul((-2, 1), poly_mul((-4, 7), (-4, 7))),
        poly_scale(poly_mul(poly_mul((-1, 1), (-1, 1)), (-1, 1)), -48),
    )
    p6_right = poly_mul((2, 1), (8, -12, 1))
    require(p6_left == p6_right, "P6+P2 B-overlap rational identity fails")
    q = lambda z: z*z - 12*z + 8
    require(q(12) > 0 and 2*12 - 12 > 0,
            "P6+P2 B-overlap polynomial is not positive for n>=12")
    # T_{n-1}-x_n reduces to 3n-4>=0.
    p6_x_left = poly_add(poly_mul((0, 1), (-2, 1)),
                         poly_scale(poly_mul((-4, 1), (-1, 1)), -1))
    require(p6_x_left == (-4, 3), "P6+P2 x_*-overlap rational identity fails")
    require(3*12 - 4 > 0, "P6+P2 x_*-overlap polynomial fails")
    require(rows[17].left == F(1, 16), "P6+P2 tail has wrong limiting endpoint")

    # The m=7,P7 finite block ends at T_1000.  For every n>=21,
    # T_{n-1}-B_n has positive numerator (n+5)(n^2-21n+14).
    require(rows[35].target == F(999, 14000), "P7 finite block does not end at T_1000")
    p7_left = poly_add(
        poly_mul((-2, 1), poly_mul((-7, 13), (-7, 13))),
        poly_scale(poly_mul(poly_mul((-1, 1), (-1, 1)), (-1, 1)), -168),
    )
    p7_right = poly_mul((5, 1), (14, -21, 1))
    require(p7_left == p7_right, "P7 B-overlap rational identity fails")
    q7 = lambda z: z*z - 21*z + 14
    require(q7(21) > 0 and 2*21 - 21 > 0,
            "P7 B-overlap polynomial is not positive for n>=21")
    # With denominator 98n(n-1), the x_* comparison has numerator
    # 3n^2+18n-28, positive and increasing for n>=21.
    p7_x_left = poly_add(poly_scale(poly_mul((0, 1), (-2, 1)), 7),
                         poly_scale(poly_mul((-7, 1), (-1, 1)), -4))
    require(p7_x_left == (-28, 18, 3), "P7 x_*-overlap rational identity fails")
    x7 = lambda z: 3*z*z + 18*z - 28
    require(x7(21) > 0 and 6*21 + 18 > 0,
            "P7 x_*-overlap polynomial fails")
    require(rows[36].left == F(1, 14), "P7 tail has wrong limiting endpoint")
    require(F(1, 16) == F(1, 2) * F(2*1, 4**2), "1/16 wrapper identity fails")
    require(F(1, 14) == F(27, 28) * F(3*2, 9**2), "1/14 wrapper identity fails")


def verify_heterogeneous(path: Path, show_witnesses: bool) -> None:
    rows = parse_heterogeneous(path)
    recovered: list[RecoveredWitness] = []
    for row in rows:
        match = re.search(r"m=(\d+)", row.family)
        require(match is not None and int(match.group(1)) == len(row.s_values),
                f"heterogeneous row {row.step}: family/s-vector size mismatch")
        require(all(s is None or s >= 2 for s in row.s_values),
                f"heterogeneous row {row.step}: invalid s value")
        witness = recover_heterogeneous_witness(row)
        recovered.append(witness)
        S, K, L, rho, T, B, xstar = satellite_data(row.s_values, witness.edges)
        require((S, K, L) == (witness.S, witness.K, witness.L),
                f"heterogeneous row {row.step}: witness statistics changed")
        require(rho == row.rho, f"heterogeneous row {row.step}: rho mismatch")
        require(T == row.target, f"heterogeneous row {row.step}: T mismatch")
        require(B == row.B, f"heterogeneous row {row.step}: B mismatch")
        require(xstar == row.xstar, f"heterogeneous row {row.step}: x_* mismatch")
        require(S <= 10 and L <= F(5, 2) and 0 < rho <= F(1, 4),
                f"heterogeneous row {row.step}: analytic-box bound fails")
        require(B <= row.left and xstar <= row.left < T,
                f"heterogeneous row {row.step}: routing inequalities fail")
        if show_witnesses:
            print(f"  low row {row.step:02d}: edges={witness.edges}; S={S}; K={K}; "
                  f"L={L}; labeled matches={witness.matches}")

    for first, last in ((0, 16), (17, 35), (36, 56)):
        for previous, current in zip(rows[first:last+1], rows[first+1:last+1]):
            require(previous.target == current.left,
                    f"heterogeneous route breaks after row {previous.step}")
    require(rows[0].left == F(1, 20), "heterogeneous route has wrong initial input")
    require(rows[56].target == F(3996, 50315), "heterogeneous route has wrong finite endpoint")
    require(rows[56].target > F(3, 38), "heterogeneous finite endpoint does not cross 3/38")
    check_low_tails(rows)
    ambiguous = sum(witness.matches > 1 for witness in recovered)
    max_matches = max(witness.matches for witness in recovered)
    print(f"[PASS] heterogeneous low module: 57 rows; exact S,K,L,rho,T,B,x_* and routes")
    print(f"       reconstructed witnesses: {ambiguous}/57 rows non-unique; "
          f"largest labeled witness class={max_matches}")


# ---------------------------------------------------------------------------
# High 42-row table and its analytic/ordinary-link companions
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class HighRow:
    row: int
    m: int
    s_values: tuple[int | None, ...]
    edges: Edges
    input_M: F
    B: F
    xstar: F
    T: F
    rho: F
    S: F
    L: F


def parse_high(path: Path) -> list[HighRow]:
    rows: list[HighRow] = []
    for line_number, raw in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        if "|" not in raw:
            continue
        parts = [part.strip() for part in raw.split("|")]
        if len(parts) != 11 or not parts[0].isdigit():
            continue
        try:
            s_values = tuple(ast.literal_eval(parts[2]))
            edges = normal_edges(ast.literal_eval(parts[3]))
            rows.append(HighRow(int(parts[0]), int(parts[1]), s_values, edges,
                                Q(parts[4]), Q(parts[5]), Q(parts[6]), Q(parts[7]),
                                Q(parts[8]), Q(parts[9]), Q(parts[10])))
        except (ValueError, SyntaxError, TypeError, CertificateError) as exc:
            raise CertificateError(f"{path.name}:{line_number}: cannot parse row") from exc
    require(len(rows) == 42, f"{path.name}: expected 42 rows, found {len(rows)}")
    require(all(row.row == i for i, row in enumerate(rows, 1)),
            "high-module row labels are not 1,...,42")
    return rows


def poly_eval(coefficients: Sequence[F], x: F) -> F:
    value = F(0)
    for coefficient in reversed(coefficients):
        value = value * x + coefficient
    return value


def poly_derivative(coefficients: Sequence[F]) -> tuple[F, ...]:
    return tuple(i * coefficients[i] for i in range(1, len(coefficients)))


def poly_add(left: Sequence[F | int], right: Sequence[F | int]) -> tuple[F, ...]:
    size = max(len(left), len(right))
    result = [F(0)] * size
    for i, value in enumerate(left):
        result[i] += F(value)
    for i, value in enumerate(right):
        result[i] += F(value)
    while len(result) > 1 and result[-1] == 0:
        result.pop()
    return tuple(result)


def poly_scale(coefficients: Sequence[F | int], scalar: F | int) -> tuple[F, ...]:
    return tuple(F(scalar) * F(value) for value in coefficients)


def poly_mul(left: Sequence[F | int], right: Sequence[F | int]) -> tuple[F, ...]:
    result = [F(0)] * (len(left) + len(right) - 1)
    for i, x in enumerate(left):
        for j, y in enumerate(right):
            result[i+j] += F(x) * F(y)
    while len(result) > 1 and result[-1] == 0:
        result.pop()
    return tuple(result)


def check_high_analytic_box() -> None:
    L, S = F(3), F(31, 5)
    # Coefficients in ascending powers of rho after substituting L=3,S=31/5.
    R1 = (
        F(384),
        -192*L - 960,
        16*L*L + 320*L + 784,
        -16*L*L - 128*L - 8*S - 208,
        L*S + 5*S,
    )
    R2 = (
        F(192),
        -96*L - 480,
        4*L*L + 152*L + 388,
        -4*L*L - 56*L - 2*S - 100,
        S,
    )
    require(poly_eval(R1, F(1, 4)) == F(16947, 160), "high-box R1 corner mismatch")
    require(poly_eval(R2, F(1, 4)) == F(64103, 1280), "high-box R2 corner mismatch")
    expected_dR1 = tuple(F(8, 5)*v for v in (-960, 2360, -1473, 124))
    expected_dR2 = tuple(F(2, 5)*v for v in (-1920, 4400, -2373, 62))
    require(poly_derivative(R1) == expected_dR1, "high-box R1 derivative identity fails")
    require(poly_derivative(R2) == expected_dR2, "high-box R2 derivative identity fails")
    inequalities = (
        F(9) + F(31, 320) + F(25) - F(64),
        F(10, 4) - 8,
        F(18) + F(31, 320) + F(80) - F(192),
        F(9, 2) + F(38) - F(96),
        F(124, 64) + F(590) - F(960),
        F(62, 64) + F(1100) - F(1920),
    )
    require(all(value < 0 for value in inequalities), "high-box monotonicity inequality fails")


def check_high_tail_and_links() -> None:
    # P3 tail, with exact rational evaluation at n=8 and polynomial certificates.
    def tail_values(n: int) -> tuple[F, F, F, F]:
        T = F(5*(n-1), 4*(8*n-3))
        B = F(30*n*(n-1), (14*n-3)**2)
        xstar = F(5*(n-1)*(2*n-3), (8*n-3)**2)
        previous = F(5*(n-2), 4*(8*n-11))
        return T, B, xstar, previous

    T8, B8, x8, _ = tail_values(8)
    require(B8 <= F(36, 251) and x8 <= F(36, 251) < T8,
            "high P3 tail does not attach at n=8")
    n_minus_2 = (-2, 1)
    n_minus_1 = (-1, 1)
    two_n_minus_3 = (-3, 2)
    eight_n_minus_3 = (-3, 8)
    eight_n_minus_11 = (-11, 8)
    fourteen_n_minus_3 = (-3, 14)
    # Clearing the positive denominators verifies the two rational-function
    # identities used in the original SymPy certificate, coefficient by coefficient.
    first_left = poly_add(
        poly_scale(poly_mul(n_minus_2, poly_mul(eight_n_minus_3, eight_n_minus_3)), 5),
        poly_scale(poly_mul(poly_mul(n_minus_1, two_n_minus_3), eight_n_minus_11), -20),
    )
    first_right = poly_scale((114, -211, 72), 5)
    require(first_left == first_right, "high P3 first difference identity fails")
    second_left = poly_add(
        poly_scale(poly_mul(n_minus_2, poly_mul(fourteen_n_minus_3, fourteen_n_minus_3)), 5),
        poly_scale(poly_mul(poly_mul((0, 1), n_minus_1), eight_n_minus_11), -120),
    )
    second_right = poly_scale((-18, -87, -20, 4), 5)
    require(second_left == second_right, "high P3 second difference identity fails")
    p1 = lambda n: 72*n*n - 211*n + 114
    p2 = lambda n: 4*n**3 - 20*n*n - 87*n - 18
    require(p1(8) > 0 and 144*8 - 211 > 0,
            "high P3 first overlap polynomial fails")
    require(p2(8) > 0 and (12*8*8 - 40*8 - 87) > 0 and 24*8 - 40 > 0,
            "high P3 second overlap polynomial fails")
    require(F(5, 4*8) == F(5, 32), "high P3 leading-coefficient limit fails")

    overlaps = (
        (F(2, 9), F(8, 27), F(7, 12)),
        (F(8, 27), F(3, 8), F(3, 5)),
        (F(3, 8), F(12, 25), F(8, 13)),
        (F(12, 25), F(5, 9), F(7, 11)),
        (F(5, 9), F(21, 32), F(13, 20)),
        (F(21, 32), F(55, 72), F(43, 64)),
        (F(55, 72), F(552, 625), F(20, 29)),
    )
    for q0, q1, s in overlaps:
        threshold = (2*s - 1) / (2*s*(1-s))
        t = s*s
        pvalue = t**3 + (1 + 4*q0)*t**2 + (3 - 4*q0)*t - 1
        require(q1 < threshold and pvalue < 0,
                f"ordinary-link overlap fails at q={q0}->{q1}")
    t = F(1, 3)
    require(t**3 + (1 + 4*F(2, 9))*t**2 + (3 - 4*F(2, 9))*t - 1 < 0,
            "ordinary-link initial crossing fails")
    t = F(1, 2)
    pfinal = t**3 + (1 + 4*F(552, 625))*t**2 + (3 - 4*F(552, 625))*t - 1
    require(pfinal == -F(41, 5000), "ordinary-link final crossing mismatch")

    # d/dr (2r^3-4r^2+7r-2)=6r^2-8r+7 has discriminant -104,
    # hence is positive.  At r=1/sqrt(2) the cubic is 4sqrt(2)-4>0.
    require((-8)**2 - 4*6*7 == -104 and 2 > 1,
            "mixed-pair cubic certificate fails")

    endpoints = (
        (F(1, 2), 8, 5, F(5, 32)),
        (F(1, 2), 6, 4, F(1, 6)),
        (F(2, 3), 6, 4, F(2, 9)),
        (F(1, 2), 12, 9, F(1, 4)),
        (F(1, 2), 30, 25, F(1, 3)),
    )
    for beta, a, b, target in endpoints:
        require(beta * F(b*(b-1), a*a) == target,
                f"wrapper endpoint {target} fails")
    for m, expected in ((3, F(2, 9)), (4, F(3, 8)), (5, F(12, 25)),
                        (6, F(5, 9)), (8, F(21, 32)), (12, F(55, 72)),
                        (25, F(552, 625))):
        require(F((m-1)*(m-2), m*m) == expected,
                f"ordinary multipartite value for m={m} fails")


def verify_high(path: Path) -> None:
    rows = parse_high(path)
    M = F(2, 25)
    for row in rows:
        context = f"high row {row.row}"
        require(row.m == len(row.s_values), f"{context}: m/s-vector mismatch")
        require(row.input_M == M, f"{context}: recorded input does not continue route")
        require(all(s is None or s >= 2 for s in row.s_values), f"{context}: invalid s value")
        graph_stats(row.m, row.edges, context)
        S, K, L, rho, T, B, xstar = satellite_data(row.s_values, row.edges)
        require(S == row.S, f"{context}: S mismatch")
        require(L == row.L, f"{context}: L mismatch")
        require(rho == row.rho, f"{context}: rho mismatch")
        require(T == row.T, f"{context}: T mismatch")
        require(B == row.B, f"{context}: B mismatch")
        require(xstar == row.xstar, f"{context}: x_* mismatch")
        require(K == S*(1-rho), f"{context}: K/rho identity mismatch")
        require(S <= F(31, 5) and L <= 3 and 0 < rho <= F(1, 4),
                f"{context}: high analytic-box bound fails")
        require(B <= M and xstar <= M < T, f"{context}: routing inequalities fail")
        M = T
    require(M == F(36, 251), "high 42-row route has wrong endpoint")
    check_high_analytic_box()
    check_high_tail_and_links()
    print("[PASS] high module: 42 explicit rows + analytic box, P3 tail, overlaps, endpoints")


# ---------------------------------------------------------------------------
# D25 terminal bridge, implemented over Fraction matrices (no SymPy)
# ---------------------------------------------------------------------------


Matrix = list[list[F]]


def matrix(rows: int, columns: int, value: F = F(0)) -> Matrix:
    return [[value for _ in range(columns)] for _ in range(rows)]


def transpose(A: Matrix) -> Matrix:
    return [list(column) for column in zip(*A)]


def matmul(A: Matrix, B: Matrix) -> Matrix:
    require(bool(A) and bool(B) and len(A[0]) == len(B), "matrix dimension mismatch")
    BT = transpose(B)
    return [[sum((a*b for a, b in zip(row, column)), F(0)) for column in BT] for row in A]


def inverse_and_determinant(A: Matrix) -> tuple[Matrix, F]:
    n = len(A)
    require(n > 0 and all(len(row) == n for row in A), "matrix must be square")
    augmented = [list(row) + [F(int(i == j)) for j in range(n)] for i, row in enumerate(A)]
    determinant = F(1)
    for column in range(n):
        pivot_row = next((i for i in range(column, n) if augmented[i][column]), None)
        require(pivot_row is not None, "matrix is singular")
        if pivot_row != column:
            augmented[column], augmented[pivot_row] = augmented[pivot_row], augmented[column]
            determinant *= -1
        pivot = augmented[column][column]
        determinant *= pivot
        augmented[column] = [entry / pivot for entry in augmented[column]]
        for i in range(n):
            if i == column:
                continue
            factor = augmented[i][column]
            if factor:
                augmented[i] = [x - factor*y for x, y in zip(augmented[i], augmented[column])]
    return [row[n:] for row in augmented], determinant


def verify_d25() -> None:
    N = 25
    support = {1, 2, 3, 5, 6, 7, 9, 10, 11, 13, 14, 15, 17, 18, 19, 21, 22, 23, 24}
    coverage = {j: 0 for j in range(1, N)}
    for i in support:
        for d in range(1, 5):
            j = (i + d) % N
            if j:
                coverage[j] += 1
    require(len(support) == 19 and set(coverage.values()) == {3},
            "D25 support dual does not cover each nonzero color three times")

    C = matrix(N, N)
    for i in range(N):
        for d in range(1, 5):
            C[i][(i+d) % N] = F(1)
    Cinv, determinant = inverse_and_determinant(C)
    require(determinant == 4, f"D25 cyclic-window determinant is {determinant}, not 4")
    Bmat = matrix(N, N)
    for i in range(N):
        for j in range(N):
            Bmat[i][j] = 25*Cinv[i][j] - F(1, 4)
            expected = F(-19 if (j-i) % 25 in {0, 4, 8, 12, 16, 20} else 6)
            require(Bmat[i][j] == expected, f"D25 inverse pattern fails at {(i, j)}")

    AH = matrix(N, N)
    for i in range(N):
        for d in (3, 8, 10):
            j = (i+d) % N
            AH[i][j] = AH[j][i] = F(1)
    AH[0][1] = AH[1][0] = F(1)
    require(sum((sum(row, F(0)) for row in AH), F(0))/2 == 76,
            "D25 auxiliary graph does not have 76 edges")
    for i, j, h in itertools.combinations(range(N), 3):
        require(not (AH[i][j] and AH[i][h] and AH[j][h]),
                f"D25 auxiliary graph contains triangle {(i, j, h)}")

    R = matmul(matmul(transpose(Bmat), AH), Bmat)
    Lset = {0, 1, 4, 5, 8, 9, 12, 13, 16, 17, 20, 21}
    Kmat = matrix(N, N)
    for i in range(N):
        for j in range(N):
            Kmat[i][j] = 297 + 475*(int(i in Lset) + int(j in Lset)) - R[i][j]
    paths = ((6, 3, 18, 15), (10, 7, 22, 19), (2, 24, 14, 11))
    expected_negative: set[Edge] = set()
    for a, b, c, d in paths:
        expected_negative.update({tuple(sorted((a, b))), tuple(sorted((b, c))), tuple(sorted((c, d)))})
    negative = {(i, j) for i in range(N) for j in range(i, N) if Kmat[i][j] < 0}
    require(negative == expected_negative, "D25 K-matrix has an unexpected negative-entry pattern")
    target = (
        (1, -1, 1, 2),
        (-1, 1, -1, 1),
        (1, -1, 1, -1),
        (2, 1, -1, 1),
    )
    for vertices in paths:
        block = tuple(tuple(Kmat[i][j] / 625 for j in vertices) for i in vertices)
        require(block == target, f"D25 principal block mismatch on {vertices}")

    require(F(19, 242) < F(3, 38), "D25 endpoint 19/242 comparison fails")
    require(F(48, 625) < F(3, 38), "D25 endpoint 48/625 comparison fails")
    require(F(3996, 50315) > F(3, 38), "D25 finite endpoint comparison fails")
    # For 0<=u<=1/19, discard the positive cubic and quadratic terms and use
    # P(u)>=18304-23028/19>0.
    require(F(18304) - F(23028, 19) > 0, "D25 cubic endpoint margin fails")
    print("[PASS] D25: support dual, det(C)=4, exact inverse/K blocks, endpoint margins")


def resolve_tables(data_dir: Path) -> tuple[Path, Path, Path]:
    names = (
        "weighted_path_chain_certificate.txt",
        "heterogeneous_satellite_chain_certificate.txt",
        "pi2_4_high_satellite_rows.txt",
    )
    paths = tuple(data_dir / name for name in names)
    missing = [str(path) for path in paths if not path.is_file()]
    require(not missing, "missing certificate table(s): " + ", ".join(missing))
    return paths  # type: ignore[return-value]


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=Path(__file__).resolve().parent,
        help="directory containing the three certificate .txt tables (default: script directory)",
    )
    parser.add_argument(
        "--show-witnesses",
        action="store_true",
        help="print the deterministically reconstructed edge set and S,K,L for all 57 low rows",
    )
    arguments = parser.parse_args(argv)
    try:
        weighted, heterogeneous, high = resolve_tables(arguments.data_dir.resolve())
        verify_weighted(weighted)
        verify_heterogeneous(heterogeneous, arguments.show_witnesses)
        verify_high(high)
        verify_d25()
    except (CertificateError, OSError) as exc:
        print(f"[FAIL] {exc}", file=sys.stderr)
        return 1
    print("ALL PI_2(4) EXACT CERTIFICATE CHECKS PASSED")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
