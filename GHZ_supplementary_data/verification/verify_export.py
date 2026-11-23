#!/usr/bin/env python3
"""Independently enumerate all perfect matchings in an exported small example.

No solver or nonstandard Python package is required. This is an exhaustive
certificate checker, not a heuristic search for a new example.
"""
from pathlib import Path
from collections import Counter
from itertools import combinations
import argparse
import json


def verify(path):
    data = json.loads(Path(path).read_text())
    planted = [frozenset(tuple(sorted(e)) for e in M) for M in data['matchings']]
    assert planted and len(set(planted)) == len(planted)
    edges = set().union(*planted)
    assert sum(map(len, planted)) == len(edges)
    vertices = set().union(*(set(e) for e in edges))
    n = len(vertices)
    assert n % 3 == 0
    for M in planted:
        assert len(M) == n // 3
        counts = Counter(v for e in M for v in e)
        assert set(counts) == vertices and set(counts.values()) == {1}
        assert all(len(e) == len(set(e)) == 3 for e in M)
    incidence = {v: sorted(e for e in edges if v in e) for v in vertices}
    actual = []
    nodes = 0

    def visit(uncovered, chosen):
        nonlocal nodes
        nodes += 1
        if not uncovered:
            actual.append(frozenset(chosen))
            return
        choices = [tuple(e for e in incidence[v] if set(e) <= uncovered)
                   for v in sorted(uncovered)]
        candidates = min(choices, key=len)
        for e in candidates:
            visit(uncovered.difference(e), chosen + [e])

    visit(vertices, [])
    assert len(actual) == len(planted) and set(actual) == set(planted)
    pairs = Counter(pair for e in edges for pair in combinations(e, 2))
    return {
        'file': str(path), 'vertices': n, 'edges': len(edges),
        'prescribed_matchings': len(planted), 'actual_matchings': len(actual),
        'search_nodes': nodes, 'maximum_codegree': max(pairs.values()),
        'matches_prescribed_exactly': True,
    }


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('export', type=Path)
    print(json.dumps(verify(p.parse_args().export), indent=2))
