#!/usr/bin/env python3
"""Exact standard-library verifier: eleven and only eleven perfect matchings
on 9q vertices for every positive integer q.

The all-q proof is the integer height certificate, not the finite exact-cover
cross-checks. Run `python verify.py --direct`.
"""
from __future__ import annotations
import argparse
import json
from collections import Counter
from itertools import combinations
from pathlib import Path

# Each row is one planted color; entries are (triple, outgoing subset).
ROWS = (
    (('123','3'),('456',''),('789','')),
    (('128','8'),('369',''),('457','4')),
    (('129',''),('356',''),('478','')),
    (('134','34'),('279',''),('568','')),
    (('137','3'),('246','6'),('589','')),
    (('145','4'),('238',''),('679','67')),
    (('147','4'),('269','26'),('358','3')),
    (('149','4'),('268','6'),('357','3')),
    (('159','9'),('248',''),('367','')),
    (('167','6'),('235','3'),('489','8')),
    (('178','8'),('259',''),('346','46')),
)
LOOPS={0,4,12,36,44,46,104,136,164,168,256}
HEIGHT_LEVELS=(
    (2,6,10,14,64,68,70,72,76,78,96,110,130,134,138,174,192,194,
     196,198,200,202,226,230,236,238,258,260,264,266,270,292,294,
     298,302,320,324,328,332,334,356,360,364,386,390,394,428,448,
     452,454,456,458,484,490,494),
    (0,136,384,424),
    (36,38,100,132,162,172,224),
    (40,46,104,160,164,170,232,288,416),
    (12,42,108,140,256,290,352),
    (8,32,168),
    (4,44,166,228,296),
    (128,),
)
FULL=511
BOUNDARY=494  # {2,3,4,6,7,8,9}; vertices 1 and 5 are never outgoing.


def require(ok: bool, message: str) -> None:
    if not ok:
        raise RuntimeError(message)


def mask(text: str) -> int:
    return sum(1 << (int(v)-1) for v in text)


def check_template():
    base=[]
    planted_states=[]
    for row in ROWS:
        require(Counter(v for e,_ in row for v in e)==Counter('123456789'),
                'A row is not a perfect matching')
        state=0
        for e,o in row:
            ee,oo=mask(e),mask(o)
            require(ee.bit_count()==3 and oo & ~ee == 0 and oo!=ee,
                    'Invalid split edge')
            base.append((ee,oo))
            state |= oo
        planted_states.append(state)
    require(len(base)==len({e for e,_ in base})==33, 'Repeated base edge')
    require(len(set(planted_states))==11 and set(planted_states)==LOOPS,
            'Incorrect planted loop states')
    require(Counter(v for e,_ in base for v in range(9) if e>>v&1)
            ==Counter({v:11 for v in range(9)}), 'Base is not 11-regular')
    outgoing=0
    for _,o in base:
        outgoing |= o
    require(outgoing==BOUNDARY, 'Wrong boundary coordinates')
    return base


def enumerate_transfer(base):
    """Completely enumerate all local configurations, without a search heuristic.

    There are nine local slots and seven possible outgoing slots. Each selected
    edge uses three distinct slots. Hence at most floor(16/3)=5 edges can be
    chosen. All sum_{k=0}^5 binom(33,k)=284274 candidate subsets are tested.
    """
    transitions=Counter()
    records=[]
    by_size=Counter()
    candidates=0
    for k in range(6):
        for selected in combinations(range(33),k):
            candidates+=1
            local=outgoing=0
            for i in selected:
                e,o=base[i]
                loc=e^o
                if local & loc or outgoing & o:
                    break
                local |= loc
                outgoing |= o
            else:
                missing=FULL ^ local
                if missing & ~BOUNDARY:
                    continue
                transitions[missing,outgoing]+=1
                by_size[k]+=1
                records.append({'edges':list(selected),
                                'source':missing,'target':outgoing})
    require(candidates==284274, 'Incorrect candidate count')
    require(dict(by_size)=={1:2,2:80,3:209,4:42}, 'Incorrect size counts')
    require(len(records)==333 and len(transitions)==250, 'Incorrect transfer size')
    print('Exhaustive enumeration: 284274 candidate subsets; 333 admissible configurations.')
    print('By selected edge count: k=1:2, k=2:80, k=3:209, k=4:42; no others.')
    return transitions,records


def check_height(transitions):
    height={s:h for h,states in enumerate(HEIGHT_LEVELS) for s in states}
    states={s for arc in transitions for s in arc}
    require(len(states)==91 and set(height)==states, 'Incomplete height certificate')
    require({s for s,t in transitions if s==t}==LOOPS, 'Incorrect loop support')
    for (s,t),w in transitions.items():
        if s==t:
            require(w==1, 'A planted self-loop has multiplicity greater than one')
        else:
            require(height[s]<height[t], f'Height does not strictly increase: {s}->{t}')
    print('Integer certificate: 91 incident states, 250 weighted arcs, 11 unit self-loops.')
    print('Every other arc strictly increases a height in {0,...,7}.')
    print('Therefore every positive-length closed walk is constant at one of the 11 loops.')
    return height


def exact_cover(n,edges):
    incident=[[i for i,e in enumerate(edges) if e>>v&1] for v in range(n)]
    require({len(x) for x in incident}=={11}, 'Lift is not 11-regular')
    full=(1<<n)-1
    solutions=[]
    nodes=0
    def visit(used,chosen):
        nonlocal nodes
        nodes+=1
        if used==full:
            solutions.append(chosen)
            return
        remaining=full ^ used
        best=None
        while remaining:
            bit=remaining & -remaining
            remaining ^= bit
            v=bit.bit_length()-1
            choices=[i for i in incident[v] if not(edges[i]&used)]
            if not choices:
                return
            if best is None or len(choices)<len(best):
                best=choices
                if len(best)==1:
                    break
        for i in best:
            visit(used|edges[i],chosen+(i,))
    visit(0,())
    return solutions,nodes


def direct_check(base,q):
    edges=[sum(1 << (9*((j+((o>>v)&1))%q)+v)
               for v in range(9) if e>>v&1)
           for j in range(q) for e,o in base]
    require(len(edges)==len(set(edges))==33*q, 'Repeated lifted edges')
    require(all(e.bit_count()==3 for e in edges), 'Nonuniform lifted edge')
    solutions,nodes=exact_cover(9*q,edges)
    require(len(solutions)==11, 'Incorrect number of perfect matchings')
    require(Counter(i for M in solutions for i in M)==Counter(range(33*q)),
            'The perfect matchings do not partition the edges')
    expected={frozenset(33*j+3*c+k for j in range(q) for k in range(3))
              for c in range(11)}
    require({frozenset(M) for M in solutions}==expected,
            'Perfect matchings do not coincide with the planted colors')
    print(f'Independent exact cover q={q}: n={9*q}, edges={33*q}, '
          f'perfect matchings=11; search nodes={nodes}.')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--direct',action='store_true')
    parser.add_argument('--export',type=Path)
    args=parser.parse_args()
    base=check_template()
    transitions,records=enumerate_transfer(base)
    height=check_height(transitions)
    if args.direct:
        for q in range(1,9):
            direct_check(base,q)
    payload={'base':base,'boundary':BOUNDARY,'configurations':records,
             'arcs':[[s,t,w] for (s,t),w in sorted(transitions.items())],
             'height_levels':HEIGHT_LEVELS,'loops':sorted(LOOPS)}
    stored=Path(__file__).with_name('eleven_certificate.json')
    if stored.exists():
        require(json.loads(stored.read_text())==json.loads(json.dumps(payload)),
                'Exported certificate disagrees with recomputation')
        print('Stored eleven_certificate.json agrees exactly with recomputation.')
    if args.export:
        args.export.write_text(json.dumps(payload,indent=2)+'\n')
    print('VERIFIED: f_3(9q) >= 11 for EVERY positive integer q.')
    print('This is a constant lower bound on an infinite arithmetic progression, '
          'not a proof that f_3(n) is unbounded.')

if __name__=='__main__':
    main()
