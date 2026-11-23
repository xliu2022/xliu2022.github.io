#!/usr/bin/env python3
"""Exact verifier for two periodic perfect-matching constructions.

Standard library only. No optimization solver, floating point computation,
or unproved extrapolation from small examples is used for the all-q results.
Run: python verify.py
Optional: python verify.py --direct --export certificate.json
"""
from __future__ import annotations
import argparse
import json
from collections import Counter, defaultdict
from itertools import combinations, permutations
from pathlib import Path

# Each pair is (base triple, subset of its vertices sent to the next layer).
TEN_TEMPLATE = (
    ('123','13'), ('124','14'), ('125','5'), ('126',''),
    ('134','34'), ('135','5'), ('136','13'), ('145',''),
    ('146','46'), ('156','5'), ('234','24'), ('235',''),
    ('236','36'), ('245',''), ('246','46'), ('256','5'),
    ('345','3'), ('346',''), ('356','56'), ('456','6'),
)
TEN_LOOPS = {4,5,16,26,28,36,37,40,56,57}
TEN_CYCLE = {(0,63):2, (63,0):1}
TEN_HEIGHT_LEVELS = (
    (1,2,3,6,7,8,17,18,19,24,27,28,33,34,35,38,39,59),
    (4,10,31,37,47,56,61,62),
    (13,25,32,36,41,44),
    (0,9,29,30,40,57,58,63),
    (5,12,14,15,26,42,46,52,53,60),
    (16,20,21,48),
    (45,),
)
LATIN = (
    'hacebdfg','ahdcebgf','gfebdach','fgbdaehc',
    'dbfghcae','begfchda','echafgbd','cdahgfeb',
)
TRANSVERSALS = (
    (0,2,3,1,5,7,4,6), (1,3,2,0,4,6,5,7),
    (2,0,1,3,7,5,6,4), (3,1,0,2,6,4,7,5),
    (4,6,7,5,0,3,1,2), (5,7,6,4,1,2,0,3),
    (6,4,5,7,3,0,2,1), (7,5,4,6,2,1,3,0),
)
# Zero-based (row,column) -> (B-shift,C-shift). Unlisted shifts are (0,0).
EIGHT_SHIFTS = {
    (0,0):(1,0), (0,3):(1,0), (1,0):(1,0), (1,2):(1,0),
    (1,4):(0,1), (2,0):(1,0), (2,2):(0,1), (2,3):(1,0),
    (3,0):(1,0), (3,2):(1,0), (3,3):(1,0), (4,2):(1,0),
    (4,3):(1,0), (5,0):(1,0), (5,1):(0,1), (5,2):(1,0),
    (5,3):(1,0), (7,0):(1,0), (7,6):(0,1),
}
EIGHT_LOOPS = {2,4,5,7,9,11,13,15}
EIGHT_CYCLE = {(3,6):5, (6,3):1}
# Complete transfer matrix, using bit labels B_1,B_3,B_4,C_e in that order.
EIGHT_ARCS = (
    (2,1,5),(2,2,1),(2,4,6),(3,5,3),(3,6,5),
    (4,1,7),(4,4,1),(5,5,1),(6,3,1),(6,5,9),
    (7,7,1),(9,9,1),(9,10,13),(9,12,2),(11,11,1),
    (11,14,4),(13,11,12),(13,13,1),(13,14,9),(15,15,1),
)


def require(ok: bool, message: str) -> None:
    if not ok:
        raise RuntimeError(message)


def bits(text: str) -> int:
    return sum(1 << (int(v)-1) for v in text)


def ten_transfer():
    base = [(bits(e),bits(o)) for e,o in TEN_TEMPLATE]
    require({e for e,_ in base} == {sum(1<<v for v in e)
            for e in combinations(range(6),3)}, 'Incorrect complete base')
    require(all(o & ~e == 0 and o != e for e,o in base), 'Invalid split')
    counts, records, by_size = Counter(), [], Counter()
    # A local configuration uses at most four split edges: each uses three
    # of the twelve distinct slots (local v / outgoing v, v in [6]).
    # Enumerate all 1+20+190+1140+4845 = 6196 candidate subsets.
    for k in range(5):
        for selected in combinations(range(20), k):
            local = outgoing = 0
            for i in selected:
                e, o = base[i]
                loc = e ^ o
                if local & loc or outgoing & o:
                    break
                local |= loc
                outgoing |= o
            else:
                s, t = 63 ^ local, outgoing
                counts[s,t] += 1
                by_size[k] += 1
                records.append({'edges':list(selected),'source':s,'target':t})
    require(dict(by_size) == {0:1,1:20,2:73,3:44,4:2}, 'Bad size counts')
    require(len(records)==140 and len(counts)==123, 'Bad transfer counts')
    height = {s:h for h,ss in enumerate(TEN_HEIGHT_LEVELS) for s in ss}
    check_height(counts, TEN_LOOPS, TEN_CYCLE, height)
    require(len(height)==55, 'Bad state count')
    by_edge = dict(base)
    pure_states = {o | by_edge[63 ^ e] for e,o in base if e & 1}
    require(pure_states == TEN_LOOPS, 'Loops do not identify the ten colors')
    print('TEN: all 6196 candidate local subsets checked; 140 configurations, '
          '123 weighted arcs, 55 incident states.')
    print('TEN: ten unit loops; only other recurrent arcs are 0 -> 63 '
          '(weight 2), 63 -> 0 (weight 1).')
    return base, counts, records, height


def check_height(counts, loops, cycle, height):
    states = {s for arc in counts for s in arc}
    require(set(height)==states, 'Height certificate does not cover states')
    require({s for (s,t) in counts if s==t}==loops, 'Wrong loop states')
    require(all(counts[s,s]==1 for s in loops), 'Nonunit planted loop')
    for arc,w in cycle.items():
        require(counts[arc]==w, 'Incorrect recurrent arc weight')
    for (s,t),w in counts.items():
        if s==t:
            continue
        if (s,t) in cycle:
            require(height[s]==height[t], 'Bad recurrent height')
        else:
            require(height[s]<height[t], f'Height fails on {s}->{t}')


def eight_transfer():
    alphabet = set('abcdefgh')
    require(all(set(r)==alphabet for r in LATIN), 'Non-Latin row')
    require(all({LATIN[r][c] for r in range(8)}==alphabet
                for c in range(8)), 'Non-Latin column')
    actual = {p for p in permutations(range(8))
              if len({LATIN[r][p[r]] for r in range(8)})==8}
    require(actual==set(TRANSVERSALS), 'Wrong transversal list')
    require(Counter((r,p[r]) for p in actual for r in range(8))
            == Counter((r,c) for r in range(8) for c in range(8)),
            'Transversals do not partition the cells')
    options=[[] for _ in range(8)]
    allowed=0
    for r in range(8):
        for c in range(8):
            sym=ord(LATIN[r][c])-97
            db,dc=EIGHT_SHIFTS.get((r,c),(0,0))
            loc=(0 if db else 1<<c)|(0 if dc else 1<<(8+sym))
            out=(1<<c if db else 0)|(1<<(8+sym) if dc else 0)
            allowed |= out
            options[r].append((c,loc,out))
    positions=(0,2,3,12)  # B_1,B_3,B_4,C_e in 16-coordinate encoding.
    require(allowed==sum(1<<p for p in positions), 'Wrong boundary')
    def pack(s):
        require(s & ~allowed == 0, 'Unencodable state')
        return sum(1<<i for i,p in enumerate(positions) if s>>p&1)
    counts, records = Counter(),[]
    def visit(r,local,outgoing,cols):
        if r==8:
            missing=65535 ^ local
            if missing & ~allowed:
                return
            s,t=pack(missing),pack(outgoing)
            counts[s,t]+=1
            records.append({'columns':list(cols),'source':s,'target':t})
            return
        for c,loc,out in options[r]:
            if not(local & loc or outgoing & out):
                visit(r+1,local|loc,outgoing|out,cols+(c,))
    visit(0,0,0,())
    require(counts==Counter({(s,t):w for s,t,w in EIGHT_ARCS}),
            'Incorrect complete eight-transfer matrix')
    require(len(records)==84, 'Incorrect eight-configuration count')
    height={2:0,7:0,9:0,13:0,15:0,4:1,3:0,6:0,
            10:1,11:1,12:1,1:2,5:1,14:2}
    check_height(counts,EIGHT_LOOPS,EIGHT_CYCLE,height)
    print('EIGHT: exactly eight seed transversals; 84 local configurations, '
          '20 weighted arcs, 14 incident states.')
    print('EIGHT: eight unit loops; only other recurrent arcs are 3 -> 6 '
          '(weight 5), 6 -> 3 (weight 1).')
    return counts,records,height


def trace_checks(counts,degree,cycle_product):
    states={s for arc in counts for s in arc}
    adjacency=defaultdict(list)
    for (s,t),w in counts.items():
        adjacency[s].append((t,w))
    power={(s,s):1 for s in states}
    for q in range(1,13):
        nxt=defaultdict(int)
        for (s,u),a in power.items():
            for t,b in adjacency[u]:
                nxt[s,t]+=a*b
        power=nxt
        val=sum(a for (s,t),a in power.items() if s==t)
        expected=degree+(2*cycle_product**(q//2) if q%2==0 else 0)
        require(val==expected, f'Bad trace at q={q}')
    print(f'Trace cross-checks q=1,...,12 agree for degree {degree}.')


def exact_cover(n,edges):
    require(len(set(edges))==len(edges), 'Repeated lifted edge')
    require(all(e.bit_count()==3 for e in edges), 'Invalid lifted triple')
    incident=[[i for i,e in enumerate(edges) if e>>v&1] for v in range(n)]
    full=(1<<n)-1
    matchings=[]
    nodes=0
    def visit(used,selected):
        nonlocal nodes
        nodes+=1
        if used==full:
            matchings.append(selected)
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
            visit(used|edges[i],selected+(i,))
    visit(0,())
    return matchings,nodes,{len(a) for a in incident}


def direct_ten(base,q):
    edges=[sum(1<<(6*((j+((o>>v)&1))%q)+v)
               for v in range(6) if e>>v&1)
           for j in range(q) for e,o in base]
    matchings,nodes,degrees=exact_cover(6*q,edges)
    expected=10+(2*2**(q//2) if q%2==0 else 0)
    require(len(matchings)==expected and degrees=={10}, 'Direct ten test fails')
    if q%2:
        require(Counter(i for M in matchings for i in M)==Counter(range(len(edges))),
                'Odd-q matchings not edge-disjoint')
    print(f'Direct TEN q={q}: n={6*q}, |E|={len(edges)}, '
          f'PM={len(matchings)}, search nodes={nodes}.')


def direct_eight(q):
    m=8*q
    edges=[]
    for j in range(q):
        for r in range(8):
            for c in range(8):
                db,dc=EIGHT_SHIFTS.get((r,c),(0,0))
                sym=ord(LATIN[r][c])-97
                vv=(8*j+r,m+8*((j+db)%q)+c,2*m+8*((j+dc)%q)+sym)
                edges.append(sum(1<<v for v in vv))
    matchings,nodes,degrees=exact_cover(3*m,edges)
    expected=8+(2*5**(q//2) if q%2==0 else 0)
    require(len(matchings)==expected and degrees=={8}, 'Direct eight test fails')
    pairs=Counter(pair for e in edges
                  for pair in combinations([v for v in range(3*m) if e>>v&1],2))
    require(max(pairs.values())==1, 'Lift is not linear')
    if q%2:
        require(Counter(i for M in matchings for i in M)==Counter(range(len(edges))),
                'Eight matchings not disjoint')
    print(f'Direct EIGHT q={q}: n={3*m}, |E|={len(edges)}, '
          f'PM={len(matchings)}, search nodes={nodes}.')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--direct',action='store_true')
    parser.add_argument('--export',type=Path)
    args=parser.parse_args()
    base,d10,rec10,h10=ten_transfer()
    d8,rec8,h8=eight_transfer()
    trace_checks(d10,10,2)
    trace_checks(d8,8,5)
    if args.direct:
        for q in range(1,10):
            direct_ten(base,q)
        for q in (1,2,3):
            direct_eight(q)
    if args.export:
        payload={
            'ten':{'template':TEN_TEMPLATE,'configurations':rec10,
                   'arcs':[[s,t,w] for (s,t),w in sorted(d10.items())],
                   'height':h10,'loops':sorted(TEN_LOOPS),
                   'extra_cycle':[[s,t,w] for (s,t),w in TEN_CYCLE.items()]},
            'eight':{'latin_square':LATIN,'transversals':TRANSVERSALS,
                     'shifts':[[r,c,b,h] for (r,c),(b,h) in sorted(EIGHT_SHIFTS.items())],
                     'configurations':rec8,
                     'arcs':[[s,t,w] for (s,t),w in sorted(d8.items())],
                     'height':h8,'loops':sorted(EIGHT_LOOPS),
                     'extra_cycle':[[s,t,w] for (s,t),w in EIGHT_CYCLE.items()]}}
        args.export.write_text(json.dumps(payload,indent=2)+'\n')
    print('VERIFIED: f_3(6q)>=10 for every positive odd q.')
    print('VERIFIED: eight pairwise edge-disjoint and no other perfect matchings '
          'in a linear tripartite 3-graph on 24q vertices, for every positive odd q.')
    print('The all-q statements follow from the integer height certificates; '
          'the small exact-cover tests are independent cross-checks only.')

if __name__=='__main__':
    main()
