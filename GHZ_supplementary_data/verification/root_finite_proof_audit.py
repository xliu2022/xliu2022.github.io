"""Independent orbit coverage and positive-control checks for the finite proof."""
from itertools import combinations, permutations, product
from pathlib import Path
from math import comb
import json

ROOT = Path(__file__).resolve().parent
pairs = ([(a, 3+b) for a in range(3) for b in range(3)]
         + [(3+b, 6+c) for b in range(3) for c in range(3)]
         + [(6+c, a) for c in range(3) for a in range(3)])
lookup = {frozenset(p): i for i, p in enumerate(pairs)}
ps = list(permutations(range(3)))
group = []
for blocks in ps:
    for inside in product(ps, repeat=3):
        vp = [3*blocks[b]+inside[b][i] for b in range(3) for i in range(3)]
        group.append([lookup[frozenset((vp[u], vp[v]))] for u, v in pairs])
assert len({tuple(g) for g in group}) == 1296
orbit_data = json.loads((ROOT/'universal_nine_missing_orbits.json').read_text())
seen = {}
for o in orbit_data['orbits']:
    mask = o['representative_mask']
    bits = [i for i in range(27) if mask >> i & 1]
    orbit = {sum(1 << g[i] for i in bits) for g in group}
    assert len(orbit) == o['orbit_size']
    for x in orbit:
        assert x not in seen
        seen[x] = mask
assert len(seen) == sum(comb(27, h) for h in range(4)) == 3304

V = frozenset(range(9))
coords = {v: divmod(v, 3) for v in V}
sums = {(0, 0), (1, 0), (0, 1), (1, 1)}
edges = {frozenset(e) for e in combinations(V, 3)
         if tuple(sum(coords[v][k] for v in e) % 3 for k in range(2)) in sums}

def partitions(vs):
    if not vs:
        yield ()
        return
    a = min(vs)
    for bc in combinations(sorted(vs-{a}), 2):
        e = frozenset((a, *bc))
        for rest in partitions(vs-e):
            yield (e,)+rest

all_partitions = list(partitions(V))
assert len(all_partitions) == 280
pms = [p for p in all_partitions if set(p) <= edges]
assert len(pms) == 13 and len(edges) == 39
assert all(sum(e in p for p in pms) == 1 for e in edges)
controls = []
for M in pms:
    blocks = list(M)
    owner = {v: i for i, X in enumerate(blocks) for v in X}
    arc_edges = {}
    for x in V:
        X = blocks[owner[x]]
        for y in V-X:
            e = (X-{x}) | {y}
            if e in edges:
                assert (y, x) not in arc_edges
                arc_edges[x, y] = e
    cycles = []
    for x, y, z in product(*blocks):
        for cyc in [((x, y), (y, z), (z, x)), ((y, x), (x, z), (z, y))]:
            if all(a in arc_edges for a in cyc):
                cycles.append(frozenset(cyc))
    used = set().union(*cycles) if cycles else set()
    assert len(used) == 3*len(cycles)
    residual = set(arc_edges)-used
    covered = set()
    free = []
    for vs in product(*blocks):
        T = frozenset(vs)
        completions = []
        for ar, br in combinations(arc_edges, 2):
            e, f = arc_edges[ar], arc_edges[br]
            if not e & f and (e | f) == V-T:
                completions.append({ar, br})
        if not completions:
            free.append(T)
        elif len(completions) == 1 and not completions[0] & used:
            covered |= completions[0]
    assert covered == residual
    a, c = len(arc_edges), len(cycles)
    assert (a-3*c) % 2 == 0
    f = (a-3*c)//2
    z = 12-c-f
    assert z >= 0 and len(free) >= 3*z
    for i in range(3):
        outgoing = sum(owner[x] == i for x, y in arc_edges)
        other_pairs = sum(owner[x] != i and owner[y] != i for x, y in arc_edges)
        assert outgoing == other_pairs
    controls.append({'arcs': a, 'triangles': c, 'fork_matchings': f,
                     'pure_matchings': z, 'free_triples': len(free)})

out = {'status': 'PASS', 'direct_group_elements': 1296,
       'all_missing_pair_sets': len(seen), 'orbit_counts': orbit_data['orbit_counts'],
       'all_partitions': 280, 'positive_control_perfect_matchings': 13,
       'positive_control_all_local_conditions_passed': True,
       'positive_control_local_counts': controls}
(ROOT/'root_finite_proof_audit.json').write_text(json.dumps(out, indent=2)+'\n')
print(json.dumps({k: v for k, v in out.items() if k != 'positive_control_local_counts'}))
