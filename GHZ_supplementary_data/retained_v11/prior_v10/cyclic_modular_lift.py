#!/usr/bin/env python3
"""Analytic cyclic lift for doubled modular rigid seeds.

Input: odd m, a single-sum rigid seed with shifts S and slope lam,
and a sheet count L coprime to m. For k shifts, put d=ceil(log2(k))
and mark {0,1,2,4,...,2**(d-2)}, using the empty set if k=1.
The theorem, not large-output enumeration, certifies rigidity.
Each part is labelled 0,...,L*m-1; (a,x) is encoded as a+m*x.
"""
from pathlib import Path
from collections import Counter
import argparse,itertools,json,math,time
from two_sheet_lift import SEEDS, structural, exact

def marking_set(k,m):
    assert 1<=k<=m and m%2
    d=(k-1).bit_length()
    D=[] if d==0 else [0]+[1<<j for j in range(d-1)]
    assert len(set(D))==d and all(0<=x<m for x in D) and sum(D)<m
    return D

def check_marking(D,m):
    count=0;d=len(D)
    for n in itertools.product(range(3),repeat=d):
        if sum(n)==d and (sum(a*x for a,x in zip(D,n))-sum(D))%m==0:
            assert all(x==1 for x in n),(D,m,n)
            count+=1
    return {'base_m':m,'marked_labels':D,'bounded_mean_solutions':count}

def construct(m,lam,S,L=2):
    S=list(S);k=len(S)
    assert m>=3 and m%2 and L>=1 and math.gcd(L,m)==1
    assert math.gcd(m,lam)==math.gcd(m,lam+1)==1
    assert len(set(S))==k and all(0<=s<m for s in S)
    D=marking_set(k,m);where={a:j for j,a in enumerate(D)}
    M=[]
    for h in [0,1]:
        for i,s in enumerate(S):
            row=[]
            for a in range(m):
                alpha=0 if a not in where else (i>>where[a])&1
                b=(lam*a+s)%m;c=(h-(lam+1)*a-s)%m
                for x in range(L):
                    row.append((a+m*x,b+m*((x+alpha)%L),c+m*((x+alpha+h)%L)))
            M.append(row)
    return M

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--full-small',action='store_true');a=ap.parse_args()
    root=Path(__file__).resolve().parent;reports=[]
    jobs=[(3,1,[0,1,2],2),(5,1,[0,1,2,3],2),(5,1,[0,1,2,3],3),
          (5,1,[0,1,2,3],4),(11,2,[0,1,2,3,5],2),
          (13,2,[0,1,2,4,5,8],2),(13,2,[0,1,2,4,5,8],3),
          (13,2,[0,1,2,4,5,8],4)]
    jobs += [(m,lam,S,2)for m,lam,S in SEEDS if m>=17]
    for m,lam,S,L in jobs:
        M=construct(m,lam,S,L)
        r={**structural(M),'base_m':m,'lambda':lam,'shifts':S,'sheets':L,'marking':check_marking(marking_set(len(S),m),m)}
        if a.full_small and m*L<=15:r['exact']=exact(M)
        reports.append(r)
        (root/'examples'/f'cyclic_seed_m{m*L}_base{m}_sheets{L}.json').write_text(json.dumps({'parameters':r,'matchings':M},indent=2)+'\n')
        print(json.dumps(r),flush=True)
    (root/'verification/cyclic_lift_checks.json').write_text(json.dumps(reports,indent=2)+'\n')
