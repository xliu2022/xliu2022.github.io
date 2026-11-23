#!/usr/bin/env python3
"""Analytic parity-changing two-sheet lift of doubled modular rigid seeds.

The input's single-layer rigidity is a hypothesis, checked separately.
Only Python's standard library is required. Vertices in each of three
parts are labelled 0,...,2m-1; (a,x) is encoded as a+m*x.
"""
from pathlib import Path
from collections import Counter
import argparse,itertools,json,math,time

SEEDS=[(3,1,[0,1,2]),(5,1,[0,1,2,3]),
       (11,2,[0,1,2,3,5]),(13,2,[0,1,2,4,5,8]),
       (17,2,[0,1,2,4,5]),(19,3,[0,1,3,12,14,15]),
       (23,2,[0,1,2,9,10,11]),(29,8,[0,1,2,7,20,21]),
       (31,4,[0,1,3,16,21,28])]

def marking_set(k,m):
    assert 1<=k<=m and m%2==1
    L=(k-1).bit_length()
    if L==0:return []
    D=[0]+[1<<j for j in range(L-1)]
    assert len(D)==L and sum(D)<m and len(set(D))==L
    return D

def construct(m,lam,S):
    S=list(S);k=len(S)
    assert m%2 and math.gcd(m,lam)==math.gcd(m,lam+1)==1
    assert len(set(S))==k and all(0<=s<m for s in S)
    D=marking_set(k,m);where={a:j for j,a in enumerate(D)}
    M=[]
    for h in [0,1]:
        for i,s in enumerate(S):
            row=[]
            for a in range(m):
                alpha=0 if a not in where else (i>>where[a])&1
                b=(lam*a+s)%m;c=(h-(lam+1)*a-s)%m
                for x in [0,1]:
                    row.append((a+m*x,b+m*(x^alpha),c+m*(x^alpha^h)))
            M.append(row)
    return M

def structural(M):
    t=len(M);q=len(M[0]);E=[e for row in M for e in row]
    assert len(set(E))==len(E)==q*t
    for row in M:
        assert len(row)==q
        for j in range(3):assert sorted(e[j]for e in row)==list(range(q))
    for j in range(3):assert Counter(e[j]for e in E)==Counter({a:t for a in range(q)})
    codeg=Counter((i,j,e[i],e[j])for e in E for i,j in [(0,1),(0,2),(1,2)])
    return {'m':q,'states':t,'vertices':3*q,'edges':q*t,'maximum_codegree':max(codeg.values())}

def check_marking(k,m):
    D=marking_set(k,m);count=0
    for delta in itertools.product([-1,0,1],repeat=len(D)):
        if sum(delta)==0 and sum(a*x for a,x in zip(D,delta))%m==0:
            assert not any(delta)
            count+=1
    return {'base_m':m,'shifts':k,'marked_labels':D,'balanced_zero_relations':count}

def exact(M):
    q=len(M[0]);E=sorted(set(e for row in M for e in row))
    masks=[];inc=[[]for _ in range(3*q)]
    for i,e in enumerate(E):
        vv=[e[j]+q*j for j in range(3)]
        masks.append(sum(1<<v for v in vv))
        for v in vv:inc[v].append(i)
    solutions=[];nodes=0
    def dfs(rem,sel):
        nonlocal nodes
        nodes+=1
        if rem==0:
            solutions.append(frozenset(E[i]for i in sel));return
        rr=rem;best=None
        while rr:
            bit=rr&-rr;rr-=bit;v=bit.bit_length()-1
            options=[i for i in inc[v]if masks[i]&rem==masks[i]]
            if not options:return
            if best is None or len(options)<len(best):best=options
            if len(best)==1:break
        for i in best:dfs(rem^masks[i],sel+[i])
    start=time.monotonic();dfs((1<<(3*q))-1,[])
    assert len(solutions)==len(M)
    assert set(solutions)=={frozenset(row)for row in M}
    return {'perfect_matchings':len(solutions),'search_nodes':nodes,'complete':True,'seconds':time.monotonic()-start}

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--full-small',action='store_true');a=ap.parse_args()
    reports=[];root=Path(__file__).resolve().parent
    for m,lam,S in SEEDS:
        M=construct(m,lam,S)
        r={**structural(M),'base_m':m,'lambda':lam,'shifts':S,'marking':check_marking(len(S),m)}
        if a.full_small and m<=5:r['exact']=exact(M)
        reports.append(r)
        (root/f'doubled_seed_m{2*m}.json').write_text(json.dumps({'parameters':r,'matchings':M},indent=2)+'\n')
        print(json.dumps(r),flush=True)
    (root/'two_sheet_checks.json').write_text(json.dumps(reports,indent=2)+'\n')
