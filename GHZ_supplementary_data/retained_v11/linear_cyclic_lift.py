#!/usr/bin/env python3
"""Linear variants of the proved cyclic lift, with optional color deletion.

The single-sum rigidity of each finite affine seed is a separate hypothesis.
The cyclic lift proves rigidity; direct pair counting verifies linearity.
Each part uses labels a+m*x in 0,...,m*L-1.
"""
from pathlib import Path
from collections import Counter
from math import gcd
import argparse,json

SEEDS={
 5:(1,[0,1,2,3]),7:(1,[0,1,2,3]),11:(2,[0,1,2,3,5]),
 13:(2,[0,1,2,4,5,8]),17:(2,[0,1,2,4,5]),
 19:(3,[0,1,3,12,14,15]),23:(2,[0,1,2,9,10,11]),
 29:(8,[0,1,2,7,20,21]),31:(4,[0,1,3,16,21,28]),
 37:(1,[0,1,2,5,21]),41:(1,[0,1,2,5,16]),
 43:(1,[0,1,2,5,38]),47:(1,[0,1,2,11,23])}


# Retained eight-state Latin input and its previously certified odd lifts.
LATIN=('hacebdfg','ahdcebgf','gfebdach','fgbdaehc',
       'dbfghcae','begfchda','echafgbd','cdahgfeb')
TRANSVERSALS=((0,2,3,1,5,7,4,6),(1,3,2,0,4,6,5,7),
 (2,0,1,3,7,5,6,4),(3,1,0,2,6,4,7,5),
 (4,6,7,5,0,3,1,2),(5,7,6,4,1,2,0,3),
 (6,4,5,7,3,0,2,1),(7,5,4,6,2,1,3,0))

def eight_latin():
    return [[(a,p[a],ord(LATIN[a][p[a]])-97)for a in range(8)]
            for p in TRANSVERSALS]

EIGHT_SHIFTS={
 (0,0):(1,0),(0,3):(1,0),(1,0):(1,0),(1,2):(1,0),
 (1,4):(0,1),(2,0):(1,0),(2,2):(0,1),(2,3):(1,0),
 (3,0):(1,0),(3,2):(1,0),(3,3):(1,0),(4,2):(1,0),
 (4,3):(1,0),(5,0):(1,0),(5,1):(0,1),(5,2):(1,0),
 (5,3):(1,0),(7,0):(1,0),(7,6):(0,1)}

def latin_eight_lift(q):
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

def marking(d,m):
 D=[] if not d else [0]+[1<<j for j in range(d-1)]
 assert sum(D)<m and all(0<=a<m for a in D)
 return D

def monotone_codes(m,S,delta=1):
 S=list(S);ss=set(S);assert gcd(m,delta)==1
 starts=[s for s in S if (s-delta)%m not in ss]
 paths=[]
 for s in starts:
  p=[]
  while s in ss:p.append(s);s=(s+delta)%m
  paths.append(p)
 assert sum(map(len,paths))==len(S), 'The shift graph contains a full cycle'
 paths.sort(key=lambda p:-len(p));order=[s for p in paths for s in p]
 for d in range((len(S)-1).bit_length(),m):
  if d and (1<<(d-1))-1>=m:break
  if any(len(p)>d+1 for p in paths):continue
  vals=sorted(range(1<<d),key=lambda x:(-x.bit_count(),-x));assigned={};used=set()
  def back(j):
   if j==len(order):return True
   s=order[j];pred=(s-delta)%m
   for c in vals:
    if c in used:continue
    if pred in assigned and c&assigned[pred]!=c:continue
    assigned[s]=c;used.add(c)
    if back(j+1):return True
    used.remove(c);del assigned[s]
   return False
  if back(0):
   out=[assigned[s]for s in S]
   assert len(set(out))==len(S)
   for s in S:
    if (s+delta)%m in ss:assert assigned[s]&assigned[(s+delta)%m]==assigned[(s+delta)%m]
   marking(d,m)
   return out
 raise ValueError('No compatible binary coding found within the marking bound')

def construct(m,lam,S,L,delta=1,codes=None,delete_h1=()):
 S=list(S);codes=list(range(len(S)))if codes is None else list(codes)
 assert m>=3 and m%2 and L>=1 and gcd(m,L)==gcd(m,delta)==1
 assert gcd(m,lam)==gcd(m,lam+1)==1
 assert len(set(S))==len(S) and all(0<=s<m for s in S)
 assert len(codes)==len(S) and len(set(codes))==len(codes) and min(codes)>=0
 d=max(codes).bit_length();D=marking(d,m);where={a:j for j,a in enumerate(D)}
 delete_h1=set(delete_h1);assert delete_h1<=set(S)
 out=[]
 for h in [0,1]:
  for s,code in zip(S,codes):
   if h and s in delete_h1:continue
   row=[]
   for a in range(m):
    alpha=0 if a not in where else (code>>where[a])&1
    b=(lam*a+s)%m;c=(delta*h-(lam+1)*a-s)%m
    for x in range(L):row.append((a+m*x,b+m*((x+alpha+h)%L),c+m*((x+alpha+2*h)%L)))
   out.append(row)
 return out

def verify(M):
 m=len(M[0]);t=len(M);E=[e for row in M for e in row]
 assert len(set(E))==m*t
 for row in M:
  for j in range(3):assert sorted(e[j]for e in row)==list(range(m))
 pairs=Counter((i,j,e[i],e[j])for e in E for i,j in [(0,1),(0,2),(1,2)])
 assert max(pairs.values())==1
 return {'m':m,'states':t,'vertices':3*m,'edges':m*t,'maximum_codegree':1,'simple':True,'balanced_tripartite':True}

def best_two_delta(m,S):
 return min((len(set(S)&{(s+d)%m for s in S}),d,sorted(set(S)&{(s+d)%m for s in S}))for d in range(1,m)if gcd(d,m)==1)

if __name__=='__main__':
 root=Path(__file__).resolve().parent;reports=[]
 for m,(lam,S) in SEEDS.items():
  for L in [2,3,4]:
   if gcd(m,L)!=1:continue
   drop=[];delta=1;codes=list(range(len(S)))
   if L==2:_,delta,drop=best_two_delta(m,S)
   if L==3:codes=monotone_codes(m,S,delta)
   M=construct(m,lam,S,L,delta,codes,drop)
   r={**verify(M),'base_m':m,'lambda':lam,'shifts':S,'sheets':L,'delta':delta,'codes':codes,'delete_h1':drop}
   reports.append(r);print(json.dumps(r),flush=True)
   (root/f'linear_seed_m{m*L}_base{m}_L{L}.json').write_text(json.dumps({'parameters':r,'matchings':M},indent=2)+'\n')
 (root/'linear_lift_checks.json').write_text(json.dumps(reports,indent=2)+'\n')
