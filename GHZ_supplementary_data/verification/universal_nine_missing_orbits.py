#!/usr/bin/env python3
"""All missing-edge orbits of K333 with at most three absent pairs.

Integer breadth-first orbit generation under adjacent transpositions
inside each part and adjacent transpositions of the three parts.
"""
from itertools import combinations
from pathlib import Path
from collections import deque, Counter
from math import comb
import json

pairs = [(a, 3+b) for a in range(3) for b in range(3)]
pairs += [(3+b, 6+c) for b in range(3) for c in range(3)]
pairs += [(6+c, a) for c in range(3) for a in range(3)]
lookup = {frozenset(e): i for i, e in enumerate(pairs)}
generators = []
for part in range(3):
    for j in range(2):
        p = list(range(9))
        u, v = 3*part+j, 3*part+j+1
        p[u], p[v] = p[v], p[u]
        generators.append(p)
for part in range(2):
    p = list(range(9))
    for j in range(3):
        u, v = 3*part+j, 3*(part+1)+j
        p[u], p[v] = p[v], p[u]
    generators.append(p)
maps = [[lookup[frozenset([p[u], p[v]])] for u, v in pairs] for p in generators]

orbits = []
for missing in range(4):
    unseen = {sum(1 << i for i in combo) for combo in combinations(range(27), missing)}
    count = 0
    while unseen:
        start = min(unseen)
        seen, queue = {start}, deque([start])
        while queue:
            mask = queue.popleft()
            bits = [i for i in range(27) if mask >> i & 1]
            for g in maps:
                image = sum(1 << g[i] for i in bits)
                if image not in seen:
                    seen.add(image); queue.append(image)
        unseen -= seen
        count += len(seen)
        orbits.append({"missing_count": missing, "representative_mask": start,
                       "orbit_size": len(seen),
                       "missing_pairs": [pairs[i] for i in range(27) if start >> i & 1]})
    assert count == comb(27, missing)
out = {"status": "all missing-edge subsets covered exactly once",
       "orbit_counts": dict(Counter(o['missing_count'] for o in orbits)),
       "orbits": orbits}
base = Path(__file__).with_suffix('')
base.with_suffix('.json').write_text(json.dumps(out, indent=2)+'\n')
base.with_suffix('.txt').write_text(''.join(f"{o['representative_mask']} {o['orbit_size']}\n" for o in orbits))
print(json.dumps(out, indent=2))
