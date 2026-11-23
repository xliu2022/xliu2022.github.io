"""Universal grafting, with global labels 0,...,3m-1.

The proof is analytic. Structural checks do not replace that proof.
"""
from collections import Counter


def globalize(seed):
    q = len(seed[0])
    return ([[tuple(g*q + e[g] for g in range(3)) for e in row]
             for row in seed], [set(range(g*q,(g+1)*q)) for g in range(3)])


def check(matchings, parts=None):
    t = len(matchings)
    n = 3*len(matchings[0])
    edges = [tuple(sorted(e)) for row in matchings for e in row]
    assert len(set(edges)) == len(edges) == t*n//3
    for row in matchings:
        assert all(len(e) == len(set(e)) == 3 for e in row)
        assert sorted(v for e in row for v in e) == list(range(n))
    assert Counter(v for e in edges for v in e) == Counter({v:t for v in range(n)})
    if parts is not None:
        assert len(parts) == 3 and all(len(P) == n//3 for P in parts)
        assert set.union(*parts) == set(range(n))
        assert all(len(set(e) & P) == 1 for e in edges for P in parts)
    return {'m':n//3, 'vertices':n, 'states':t, 'edges':len(edges)}


def graft(matchings, seed, parts=None, selected=None):
    """Replace any three old vertices by one balanced tripartite seed.

    seed is tripartite: its triples use part-relative coordinates.
    The seed edge through first-part vertex 0 is deleted in each color.
    """
    n = 3*len(matchings[0])
    t, q = len(matchings), len(seed[0])
    assert len(seed) == t
    if selected is None:
        selected = [min(P) for P in parts] if parts is not None else [0,1,2]
    assert len(selected) == len(set(selected)) == 3
    assert set(selected) <= set(range(n))
    if parts is not None:
        assert all(selected[g] in parts[g] for g in range(3))
    slots = {v:g for g,v in enumerate(selected)}
    old = [v for v in range(n) if v not in slots]
    old_index = {v:k for k,v in enumerate(old)}
    offset = len(old)
    result = []
    for i, row in enumerate(matchings):
        removed = next(e for e in seed[i] if e[0] == 0)
        new_row = []
        for e in row:
            new_row.append(tuple(sorted(
                old_index[v] if v in old_index else
                offset + slots[v]*q + removed[slots[v]] for v in e)))
        new_row.extend(tuple(offset + g*q + e[g] for g in range(3))
                       for e in seed[i] if e[0] != 0)
        result.append(new_row)
    new_parts = None
    if parts is not None:
        new_parts = [{old_index[v] for v in parts[g] if v in old_index} |
                     set(range(offset+g*q,offset+(g+1)*q)) for g in range(3)]
    check(result, new_parts)
    return result, new_parts
