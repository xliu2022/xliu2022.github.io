#!/usr/bin/env python3
"""Independent exhaustive checks of explicit compact rigid 3-graphs.

The recursion always covers the least uncovered vertex. Every edge containing
that vertex and disjoint from the previously chosen edges is explored. Thus
each perfect matching occurs once, with no heuristic assumption or solver.
"""
import json
from pathlib import Path
from collections import Counter

HERE=Path(__file__).resolve().parent
def verify_file(path):
    data=json.loads(path.read_text())
    n=data['n']
    d=len(data['matchings'])
    triples=[tuple(e) for e in data['edges']]
    assert len(triples)==d*n//3 and len(set(triples))==len(triples)
    assert all(len(e)==3 and len(set(e))==3 and min(e)>=1 and max(e)<=n for e in triples)
    masks=[sum(1<<(v-1) for v in e) for e in triples]
    expected={tuple(sorted(M)) for M in data['matchings']}
    assert len(expected)==d
    assert Counter(i for M in expected for i in M)==Counter(range(len(triples)))
    full=(1<<n)-1
    for M in expected:
        covered=0
        for i in M:
            assert not covered&masks[i]
            covered |=masks[i]
        assert covered==full
    incident=[tuple(i for i,e in enumerate(masks) if e>>v&1) for v in range(n)]
    assert all(len(row)==d for row in incident)
    found=[]
    nodes=0
    def visit(covered,chosen):
        nonlocal nodes
        nodes+=1
        if covered==full:
            found.append(tuple(sorted(chosen)))
            return
        v=((full^covered)&-(full^covered)).bit_length()-1
        for i in incident[v]:
            if not masks[i]&covered:
                visit(covered|masks[i],chosen+(i,))
    visit(0,())
    assert len(found)==d
    assert set(found)==expected
    print(f'n={n}: {len(triples)} distinct triples, {d}-regular; '
          f'exhausted {nodes} recursive nodes; exactly the {d} planted perfect matchings.')
    return nodes
if __name__=='__main__':
    for filename in ('nine_on_6.json','nine_on_9.json','nine_on_12.json','12_on_12.json','nine_on_15.json','10_on_15.json','11_on_18.json','nine_on_21.json'):
        verify_file(HERE/filename)
