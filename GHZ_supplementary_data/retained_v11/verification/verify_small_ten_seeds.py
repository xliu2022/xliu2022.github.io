#!/usr/bin/env python3
"""Exhaustive verification of ten rigid states at each m=2,...,11."""
from pathlib import Path
import json
from collections import Counter

def verify(record):
 m=record['m'];n=3*m;M=record['matchings'];assert len(M)==10
 E=[tuple(sorted(e))for row in M for e in row]
 assert len(E)==len(set(E))==10*m
 assert all(len(e)==len(set(e))==3 for e in E)
 assert all(Counter(v for e in row for v in e)==Counter(range(n))for row in M)
 if 'parts'in record:
  P=[set(z)for z in record['parts']];assert len(P)==3 and all(len(z)==m for z in P)
  assert set.union(*P)==set(range(n)) and sum(map(len,P))==n
  assert all(len(set(e)&z)==1 for e in E for z in P)
 masks=[sum(1<<v for v in e)for e in E];inc=[[i for i,e in enumerate(masks)if e>>v&1]for v in range(n)]
 full=(1<<n)-1;nodes=0;sols=[]
 def visit(used,picked):
  nonlocal nodes
  nodes+=1
  if used==full:sols.append(frozenset(picked));return
  rest=full^used;v=(rest&-rest).bit_length()-1
  for i in inc[v]:
   if not used&masks[i]:visit(used|masks[i],picked+[i])
 visit(0,[])
 expected={frozenset(range(c*m,(c+1)*m))for c in range(10)}
 assert len(sols)==10 and set(sols)==expected
 return {'m':m,'vertices':n,'states':10,'edges':10*m,'nodes':nodes,'complete':True}
if __name__=='__main__':
 root=Path(__file__).resolve().parent
 R=json.loads((root/'ten_states_orders_two_to_eleven.json').read_text())
 assert [r['m']for r in R]==list(range(2,12))
 checks=[verify(r)for r in R]
 (root/'small_ten_seed_checks.json').write_text(json.dumps(checks,indent=2)+'\n')
 print(json.dumps(checks,indent=2))
