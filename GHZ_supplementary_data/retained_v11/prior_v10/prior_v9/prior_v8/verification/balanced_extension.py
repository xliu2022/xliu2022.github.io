#!/usr/bin/env python3
"""Balanced binary-port extensions; pure Python analytic-construction checks.

The mathematical proof is in balanced_extension.tex.  This program verifies
binary-gadget inequalities, small seeds, degree/edge-disjointness invariants,
and small complete examples.  No extrapolation from enumeration is needed.

python balanced_extension.py
python balanced_extension.py --m 270 --export example_m270.json
"""
from collections import Counter,deque
from itertools import permutations
from pathlib import Path
import argparse,json,math

LATIN=('hacebdfg','ahdcebgf','gfebdach','fgbdaehc',
       'dbfghcae','begfchda','echafgbd','cdahgfeb')
TRANSVERSALS=((0,2,3,1,5,7,4,6),(1,3,2,0,4,6,5,7),
 (2,0,1,3,7,5,6,4),(3,1,0,2,6,4,7,5),
 (4,6,7,5,0,3,1,2),(5,7,6,4,1,2,0,3),
 (6,4,5,7,3,0,2,1),(7,5,4,6,2,1,3,0))

def eight_coprime(m):
    assert m>=5 and math.gcd(m,6)==1
    return [[(a,(a+s)%m,(h-2*a-s)%m)for a in range(m)]
            for h in range(2)for s in range(4)]

def eight_latin():
    return [[(a,p[a],ord(LATIN[a][p[a]])-97)for a in range(8)]
            for p in TRANSVERSALS]

EIGHT_SHIFTS={
 (0,0):(1,0),(0,3):(1,0),(1,0):(1,0),(1,2):(1,0),
 (1,4):(0,1),(2,0):(1,0),(2,2):(0,1),(2,3):(1,0),
 (3,0):(1,0),(3,2):(1,0),(3,3):(1,0),(4,2):(1,0),
 (4,3):(1,0),(5,0):(1,0),(5,1):(0,1),(5,2):(1,0),
 (5,3):(1,0),(7,0):(1,0),(7,6):(0,1)}

def eight_latin_lift(q):
    """Prior separately certified family, valid only for odd q."""
    assert q>=1 and q%2==1
    out=[]
    for p in TRANSVERSALS:
        row=[]
        for layer in range(q):
            for a in range(8):
                b=p[a];c=ord(LATIN[a][b])-97
                db,dc=EIGHT_SHIFTS.get((a,b),(0,0))
                row.append((8*layer+a,8*((layer+db)%q)+b,
                            8*((layer+dc)%q)+c))
        out.append(row)
    return out

def six_cyclic(m):
    assert m>=3
    cycle=list(range(m)) if m%2 else [0,2,1]+list(range(3,m))
    tau={cycle[j]:cycle[(j+1)%m]for j in range(m)}
    out=[]
    for step,offset in ((-1,0),(0,-1),(1,1)):
        for c in (0,1):
            out.append([(x,(x+step)%m,
                         tau[(x+offset)%m]if c else (x+offset)%m)
                        for x in range(m)])
    return out

def binary_gadget(t):
    X=[x for x in range(1<<t)if x.bit_count()%2==0]
    Y=[x for x in range(1<<t)if x.bit_count()%2==1]
    iy={y:j for j,y in enumerate(Y)}
    P=[[iy[x^(1<<i)]for x in X]for i in range(t)]
    B=[[(x&((1<<i)-1)).bit_count()%2 for x in X]for i in range(t)]
    return X,P,B

def rank(t,i,x):
    s=[(x&((1<<j)-1)).bit_count()%2 for j in range(t)]
    return sum((s[(i+a)%t]^s[i])<<(a-1)for a in range(1,t))

def check_gadget(t):
    X,P,B=binary_gadget(t);L=len(X)
    assert len({(x,B[i][x],P[i][x])for i in range(t)for x in range(L)})==t*L
    count=0
    for i in range(t):
        inv={y:x for x,y in enumerate(P[i])}
        ranks=[rank(t,i,x)for x in X]
        assert sorted(ranks)==list(range(L))
        for j in range(t):
            if i==j:continue
            for x in range(L):
                if B[i][x]==B[j][x]:
                    y=inv[P[j][x]]
                    assert X[y]==X[x]^(1<<i)^(1<<j)
                    assert ranks[x]<ranks[y]
                    count+=1
    return {'states':t,'rows':L,'strict_rank_inequalities':count}

def extend(matchings,substitutions=None):
    """substitutions maps distinct gadget row indices to t-rigid seeds."""
    substitutions={}if substitutions is None else substitutions
    t=len(matchings);m=len(matchings[0]);X,P,B=binary_gadget(t);L=len(X)
    assert m>=L
    out=[]
    for i,M in enumerate(matchings):
        N=[]
        for a,b,c in M:
            aa=a if a>=L or 1-B[i][a]==0 else m+a
            N.append((aa,b,c))
        new_offset=m+L
        for x in range(L):
            port=x if B[i][x]==0 else m+x
            oldB=m+x;oldC=m+P[i][x]
            if x not in substitutions:
                N.append((port,oldB,oldC));continue
            seed=substitutions[x];q=len(seed[0]);assert len(seed)==t
            for a,b,c in seed[i]:
                N.append((port if a==0 else new_offset+a-1,
                          oldB if b==0 else new_offset+b-1,
                          oldC if c==0 else new_offset+c-1))
            new_offset+=q-1
        out.append(N)
    return out

def structural_check(M):
    t=len(M);m=len(M[0]);all_edges=[]
    for row in M:
        assert len(row)==m
        for k in range(3):assert sorted(e[k]for e in row)==list(range(m))
        all_edges.extend(row)
    assert len(set(all_edges))==len(all_edges)==t*m
    for k in range(3):assert Counter(e[k]for e in all_edges)==Counter({a:t for a in range(m)})
    return {'m':m,'states':t,'vertices':3*m,'edges':t*m}

def exact_count(M,cap=10000):
    m=len(M[0]);options=[set()for _ in range(m)]
    for row in M:
        for a,b,c in row:options[a].add((b,c))
    count=0;nodes=0
    def visit(a,usedB,usedC):
        nonlocal count,nodes
        nodes+=1
        if a==m:count+=1;return
        for b,c in options[a]:
            if (usedB>>b)&1 or (usedC>>c)&1:continue
            visit(a+1,usedB|(1<<b),usedC|(1<<c))
            assert count<=cap
    visit(0,0,0)
    return count,nodes

def construct_states(m,t):
    """t=7 or8; includes the separately certified old Latin periodic family."""
    assert t in (7,8)
    if m>=5 and math.gcd(m,6)==1:
        return eight_coprime(m)[:t],{'method':'coprime_six','base_m':m}
    if m>=8 and m%8==0 and (m//8)%2==1:
        return eight_latin_lift(m//8)[:t],{'method':'prior_latin_periodic','q':m//8}
    L=1<<(t-1);residue=m%6
    if (t==8 and residue==3 and m>=261) or (t==7 and residue==3 and m>=129):
        delta=L;subs={};method='binary'
    elif ((t==8 and residue in (2,4) and m>=(266 if residue==2 else 268))
          or (t==7 and residue in (0,4) and m>=(138 if residue==0 else 136))):
        delta=L+7;subs={0:eight_latin()[:t]};method='binary_latin'
    elif (t==8 and residue==0 and m>=270) or (t==7 and residue==2 and m>=140):
        delta=L+11;subs={0:eight_latin()[:t],1:eight_coprime(5)[:t]};method='binary_two_seeds'
    else:
        raise ValueError('Not covered by this particular family')
    base=m-delta
    assert base>=L and math.gcd(base,6)==1
    out=extend(eight_coprime(base)[:t],subs)
    assert len(out[0])==m
    return out,{'method':method,'base_m':base,'increment':delta}

def construct_eight(m):return construct_states(m,8)
def construct_seven(m):return construct_states(m,7)

def construct_best(m):
    for t in (8,7):
        try:return construct_states(m,t)
        except ValueError:pass
    if m>=3:return six_cyclic(m),{'method':'prior_cyclic_six'}
    if m==2:
        return [[(0,b,c),(1,1-b,1-c)]for b in (0,1)for c in (0,1)],{'method':'complete_two'}
    if m==1:return [[(0,0,0)]],{'method':'single_edge'}
    raise ValueError('m must be positive')

def main():
    p=argparse.ArgumentParser();p.add_argument('--m',type=int);p.add_argument('--export',type=Path)
    a=p.parse_args()
    if a.m is not None:
        M,method=construct_best(a.m);result={**structural_check(M),**method}
        if a.export:a.export.write_text(json.dumps({'parameters':result,'matchings':M},indent=2)+'\n')
        print(json.dumps(result));return
    reports={'gadgets':[check_gadget(t)for t in range(2,10)],'seeds':[],
             'small_complete_checks':[],'large_structural_checks':[]}
    for seed in (eight_coprime(5),eight_latin()):
        s=structural_check(seed);cnt,nodes=exact_count(seed);assert cnt==8
        reports['seeds'].append({**s,'perfect_matchings':cnt,'search_nodes':nodes})
    for subs in ({},{0:six_cyclic(3)[:3]}):
        M=extend(six_cyclic(4)[:3],subs)
        s=structural_check(M);cnt,nodes=exact_count(M);assert cnt==3
        reports['small_complete_checks'].append({**s,'perfect_matchings':cnt,'search_nodes':nodes})
    for m in (261,266,268,270,265,391,396,999,1000):
        M,method=construct_eight(m)
        reports['large_structural_checks'].append({**structural_check(M),**method})
    # This checks the arithmetic covering, not perfect-matchings by extrapolation.
    for m in range(265,1001):
        M,method=construct_eight(m);assert len(M[0])==m
    reports['arithmetic_coverage']={'eight_from':263,'seven_from':135,'through':1000,'all_covered':True}
    for m in range(135,1001):
        M,method=construct_seven(m);assert len(M[0])==m
    for m in range(263,1001):
        M,method=construct_eight(m);assert len(M[0])==m
    records=[]
    for m in range(1,1001):
        M,method=construct_best(m)
        records.append({**structural_check(M),**method})
    Path(__file__).with_name('small_m_1_to_1000_constructions.json').write_text(json.dumps(records,indent=2)+'\n')
    dest=Path(__file__).with_name('balanced_extension_checks.json')
    dest.write_text(json.dumps(reports,indent=2)+'\n');print(json.dumps(reports,indent=2))
if __name__=='__main__':main()
