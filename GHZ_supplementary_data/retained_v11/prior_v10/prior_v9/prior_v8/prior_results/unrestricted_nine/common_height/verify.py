#!/usr/bin/env python3
"""Exact verifier for nine disjoint perfect matchings at every order n >= 24,
3 | n. Standard library only. The all-order proof is the common-height
certificate; --direct performs optional independent exact-cover checks.
"""
from __future__ import annotations
import argparse
import json
from collections import Counter
from itertools import combinations
from math import comb
from pathlib import Path

# Vertices 1,2,3,5,6 form the common boundary. Vertex 4 never crosses blocks.
# Entries: (base triple, vertices sent to the next block's boundary).
A_ROWS = (
    (('178','1'),('239',''),('456','')),
    (('578','5'),('136',''),('249','2')),
    (('678',''),('139',''),('245','')),
    (('127','12'),('468',''),('359','')),
    (('147','1'),('238','3'),('569','')),
    (('267','2'),('358','3'),('149','1')),
    (('679','6'),('258',''),('134','')),
    (('347','3'),('189','1'),('256','5')),
    (('457','5'),('689',''),('123','23')),
)
B_ROWS = (
    (('124','1'),('356','')),
    (('135','5'),('246','2')),
    (('234',''),('156','')),
    (('145','1'),('236','2')),
    (('235','3'),('146','1')),
    (('245','2'),('136','13')),
    (('125',''),('346','6')),
    (('345','35'),('126','1')),
    (('123','23'),('456','5')),
)
PHASES = (1,18,0,3,5,7,32,21,22)
HEIGHT_LEVELS = {
    1: (0,18,23),
    2: (2,17,21,52),
    3: (5,19,32),
    4: (3,4,20,22,36),
    5: (1,38),
    6: (7,16),
}
BOUNDARY = 55  # {1,2,3,5,6}
DIM = 64

def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)

def mask(text: str) -> int:
    return sum(1 << (int(v)-1) for v in text)

def check_rows(rows, size):
    base=[]
    for c,row in enumerate(rows):
        require(Counter(v for e,_ in row for v in e)==Counter(str(i)for i in range(1,size+1)),
                'A color does not partition the base vertices')
        state=0
        for e,o in row:
            em,om=mask(e),mask(o)
            require(em.bit_count()==3 and om & ~em == 0 and em!=om,
                    'Invalid split triple')
            require(om & ~BOUNDARY == 0,'An outgoing vertex is not on the boundary')
            base.append((em,om))
            state |= om
        require(state==PHASES[c], 'Incompatible planted boundary states')
    require(len(base)==len({e for e,_ in base}), 'Repeated base triple')
    return base

def enumerate_transfer(base,size):
    """Check EVERY candidate subset of size <= floor((size+5)/3)."""
    full=(1<<size)-1
    limit=(size+BOUNDARY.bit_count())//3
    arcs=Counter()
    configs=[]
    counts=Counter()
    candidates=0
    for k in range(limit+1):
        for chosen in combinations(range(len(base)),k):
            candidates+=1
            local=outgoing=0
            for i in chosen:
                e,o=base[i]
                l=e^o
                if local & l or outgoing & o:
                    break
                local |= l
                outgoing |= o
            else:
                missing=full^local
                if missing & ~BOUNDARY:
                    continue
                arcs[missing,outgoing]+=1
                counts[k]+=1
                configs.append([list(chosen),missing,outgoing])
    require(candidates==sum(comb(len(base),k)for k in range(limit+1)),
            'Incomplete exhaustive subset enumeration')
    D=[[0]*DIM for _ in range(DIM)]
    for (s,t),w in arcs.items(): D[s][t]=w
    print(f'{size}-vertex block: {len(base)} edges; {candidates} candidate subsets; '
          f'{len(configs)} admissible configurations; {len(arcs)} weighted arcs; '
          f'counts by size {dict(sorted(counts.items()))}.')
    return D,arcs,configs,candidates

def matmul(X,Y):
    Z=[[0]*DIM for _ in range(DIM)]
    for s in range(DIM):
        for t in range(DIM):
            if X[s][t]:
                for u in range(DIM):
                    if Y[t][u]: Z[s][u]+=X[s][t]*Y[t][u]
    return Z

def check_common_height(A,C):
    height=[0]*DIM
    for h,ss in HEIGHT_LEVELS.items():
        for s in ss:height[s]=h
    expected=[set(PHASES),set(PHASES)|{2}]
    for name,D,loops in zip(('A','C=BA'),(A,C),expected):
        require({s for s in range(DIM) if D[s][s]}==loops,
                f'Wrong diagonal support for {name}')
        require(all(D[s][s]==1 for s in loops), f'Nonunit loop in {name}')
        for s in range(DIM):
            for t in range(DIM):
                if s!=t and D[s][t]:
                    require(height[s]<height[t], f'{name} does not raise height: {s}->{t}')
        print(f'{name}: {sum(D[s][t]>0 for s in range(DIM)for t in range(DIM))} weighted arcs; '
              f'{len(loops)} unit loops; every nonloop strictly raises common height.')
    return height

def word_for_n(n):
    require(n>=24 and n%3==0,'The explicit all-order rule requires n>=24 and 3|n')
    if n in (24,27,30):return {24:'ABA',27:'AAA',30:'BBBBB'}[n]
    # n=9a+15b, a>=1. b=0,1,2 is fixed by n modulo 9.
    b={0:0,6:1,3:2}[n%9]
    a=(n-15*b)//9
    require(a>=1 and 9*a+15*b==n,'Failed arithmetic decomposition')
    return 'A'*a+'BA'*b

def build_hypergraph(word):
    sizes=[9 if t=='A' else 6 for t in word]
    offsets=[0]
    for size in sizes:offsets.append(offsets[-1]+size)
    edges=[]
    planted=[set()for _ in PHASES]
    for j,letter in enumerate(word):
        nxt=(j+1)%len(word)
        rows=A_ROWS if letter=='A' else B_ROWS
        for c,row in enumerate(rows):
            for e,o in row:
                vs=[(offsets[nxt] if v in o else offsets[j])+int(v)-1 for v in e]
                require(len(set(vs))==3,'A constructed edge has repeated vertices')
                planted[c].add(len(edges))
                edges.append(sum(1<<v for v in vs))
    n=offsets[-1]
    require(len(edges)==3*n and len(set(edges))==len(edges),'Wrong edge count or repeated edge')
    full=(1<<n)-1
    for M in planted:
        covered=0
        for i in M:
            require(not covered&edges[i],'A planted color is not a matching')
            covered |= edges[i]
        require(covered==full,'A planted matching is not perfect')
    return n,edges,{frozenset(M)for M in planted}

def exact_cover(n,edges):
    incident=[[i for i,e in enumerate(edges)if e>>v&1]for v in range(n)]
    require(all(len(x)==9 for x in incident),'Constructed hypergraph is not 9-regular')
    full=(1<<n)-1
    solutions=[]
    nodes=0
    def dfs(covered,chosen):
        nonlocal nodes
        nodes+=1
        if covered==full:
            solutions.append(frozenset(chosen))
            return
        unused=full^covered
        best=None
        while unused:
            bit=unused&-unused
            unused^=bit
            v=bit.bit_length()-1
            choices=[i for i in incident[v] if not edges[i]&covered]
            if not choices:return
            if best is None or len(choices)<len(best):
                best=choices
                if len(best)==1:break
        for i in best:dfs(covered|edges[i],chosen+(i,))
    dfs(0,())
    return solutions,nodes

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--direct',action='store_true',help='Independent exact-cover checks')
    parser.add_argument('--export',action='store_true',help='Write all finite certificate data')
    args=parser.parse_args()
    aa=check_rows(A_ROWS,9);bb=check_rows(B_ROWS,6)
    A,arcsA,configsA,cntA=enumerate_transfer(aa,9)
    B,arcsB,configsB,cntB=enumerate_transfer(bb,6)
    C=matmul(B,A)
    height=check_common_height(A,C)
    B5=B
    for _ in range(4):B5=matmul(B5,B)
    require(sum(B5[s][s] for s in range(DIM))==9, 'The 30-vertex exceptional word failed')
    print('Exceptional word B^5: trace(B^5)=9, certifying the 30-vertex example.')
    for n in range(24,3001,3):
        w=word_for_n(n)
        require(sum(9 if x=='A'else 6 for x in w)==n,'Bad word length')
    print('Arithmetic check: n=9a+15b, a>=1, b in {0,1,2}, for all 3|n, n>=33.')
    if args.direct:
        for n in (24,27,30,33,36,39,42):
            w=word_for_n(n)
            nn,edges,planted=build_hypergraph(w)
            sol,nodes=exact_cover(nn,edges)
            require(len(sol)==9 and set(sol)==planted,'Exact-cover cross-check failed')
            require(Counter(i for M in sol for i in M)==Counter(range(len(edges))),
                    'Perfect matchings do not partition edges')
            print(f'Independent exact cover: word={w}; n={nn}; edges={len(edges)}; '
                  f'perfect matchings=9; search nodes={nodes}.')
    payload={
        'A_rows':A_ROWS,'B_rows':B_ROWS,'phase_states':PHASES,
        'A_configurations':configsA,'B_configurations':configsB,
        'A_arcs':[[s,t,w]for(s,t),w in sorted(arcsA.items())],
        'B_arcs':[[s,t,w]for(s,t),w in sorted(arcsB.items())],
        'C_arcs':[[s,t,C[s][t]]for s in range(DIM)for t in range(DIM)if C[s][t]],
        'height':height,'trace_B5':sum(B5[s][s]for s in range(DIM)),
    }
    payload=json.loads(json.dumps(payload))
    file=Path(__file__).with_name('certificate.json')
    if file.exists():require(json.loads(file.read_text())==payload,'Stored certificate mismatch')
    if args.export:file.write_text(json.dumps(payload,indent=2)+'\n')
    print('VERIFIED: f_3(n)>=9 for every n>=24 divisible by 3.')
    print('This improves the all-large-order constant lower bound, not unboundedness or the asymptotic order.')

if __name__=='__main__':main()
