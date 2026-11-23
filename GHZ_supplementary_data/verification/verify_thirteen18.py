"""Independent exhaustive verification of a13-state rigid 3-graph on18vertices.
Standard Python only. Run from any directory.
"""
from pathlib import Path
import json

def verify():
    path=Path(__file__).with_name('13_on_18.json')
    data=json.loads(path.read_text())
    n=data['n'];triples=[tuple(e) for e in data['edges']]
    assert n==18 and len(triples)==78
    assert len(set(triples))==78
    assert all(len(e)==len(set(e))==3 and min(e)>=1 and max(e)<=n for e in triples)
    masks=[sum(1<<(v-1) for v in e) for e in triples]
    full=(1<<n)-1
    specified={frozenset(m) for m in data['matchings']}
    assert len(specified)==13
    assert sorted(e for m in specified for e in m)==list(range(78))
    for m in specified:
        covered=0
        for e in m:
            assert not covered&masks[e]
            covered|=masks[e]
        assert covered==full
    incident=[[i for i,e in enumerate(masks) if e>>v&1] for v in range(n)]
    assert all(len(a)==13 for a in incident)
    actual=[];nodes=0
    def search(used,chosen):
        nonlocal nodes
        nodes+=1
        if used==full:
            actual.append(frozenset(chosen));return
        v=next(v for v in range(n) if not used>>v&1)
        for i in incident[v]:
            if not used&masks[i]:search(used|masks[i],chosen+[i])
    search(0,[])
    assert len(actual)==13 and set(actual)==specified
    return {'vertices':n,'edges':len(triples),'degree':13,'perfect_matchings':len(actual),'search_nodes':nodes,'matches_prescribed_exactly':True}
if __name__=='__main__':print(json.dumps(verify(),indent=2))
