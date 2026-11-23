#!/usr/bin/env python3
"""Complete independent standard-library exact-cover checks of small linear lifts."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parent.parent))
import json,time
from linear_cyclic_lift import construct,verify,monotone_codes,SEEDS,best_two_delta

def exact(M):
 m=len(M[0]);E=sorted(set(e for row in M for e in row));n=3*m
 inc=[[]for _ in range(n)];masks=[]
 for k,e in enumerate(E):
  vv=[e[j]+m*j for j in range(3)];masks.append(sum(1<<v for v in vv))
  for v in vv:inc[v].append(k)
 nodes=0;solutions=[]
 def rec(rem,selected):
  nonlocal nodes
  nodes+=1
  if not rem:solutions.append(frozenset(E[k]for k in selected));return
  rr=rem;best=None
  while rr:
   bit=rr&-rr;rr-=bit;v=bit.bit_length()-1
   options=[k for k in inc[v]if masks[k]&rem==masks[k]]
   if not options:return
   if best is None or len(options)<len(best):best=options
   if len(best)==1:break
  for k in best:rec(rem^masks[k],selected+[k])
 st=time.monotonic();rec((1<<n)-1,[])
 assert len(solutions)==len(M) and set(solutions)=={frozenset(row)for row in M}
 return {'perfect_matchings':len(solutions),'nodes':nodes,'complete':True,'seconds':time.monotonic()-st}

if __name__=='__main__':
 records=[]
 for m,L in [(5,2),(7,2),(5,3),(5,4)]:
  lam,S=SEEDS[m];delta=1;drop=[];codes=list(range(len(S)))
  if L==2:_,delta,drop=best_two_delta(m,S)
  if L==3:codes=monotone_codes(m,S,delta)
  M=construct(m,lam,S,L,delta,codes,drop)
  r={**verify(M),'base_m':m,'sheets':L,'delta':delta,'codes':codes,'delete_h1':drop,'enumeration':exact(M)}
  records.append(r);print(json.dumps(r),flush=True)
 Path(__file__).with_name('linear_small_exact_checks.json').write_text(json.dumps(records,indent=2)+'\n')
