#!/usr/bin/env python3
"""Check the finite seed bank and every residue witness for linear all-order bounds.

This checks explicit algebraic constructions and integer witnesses. Rigidity
of a cyclic output follows from the analytic lift and retained base certificates.
No finite search is being extrapolated to an infinite family.
"""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parent.parent))
import json
from linear_cyclic_lift import SEEDS,construct,verify,monotone_codes,latin_eight_lift
ROOT=Path(__file__).resolve().parent

def build_seed(q,meta):
 if meta['method']=='old_linear_latin_lift':return latin_eight_lift(meta['sheets'])
 m=meta['base_m'];lam,S=SEEDS[m];L=meta['sheets'];delta=meta['delta']
 codes=monotone_codes(m,S,delta)if L==3 else list(range(len(S)))
 return construct(m,lam,S,L,delta,codes,meta.get('delete_h1',()))

if __name__=='__main__':
 records=json.loads((ROOT/'linear_thresholds_refined.json').read_text());out=[]
 done={}
 for row in records:
  t=row['states'];a=row['min_increment'];N=row['all_orders_from']
  used={int(q):meta for q,meta in row['used_seeds'].items()}
  assert a+1 in used
  for q,meta in used.items():
   key=json.dumps(meta,sort_keys=True)
   if key not in done:
    check=verify(build_seed(q,meta));assert check['m']==q and check['states']==meta['states']
    done[key]=check
   assert meta['states']>=t
  assert len(row['apery'])==a
  for r,x in enumerate(row['apery']):
   assert x['residue']==r and x['value']%a==r
   assert sum(x['increments'])==x['value']
   assert all(z+1 in used for z in x['increments'])
  assert max(x['value']for x in row['apery'])-a+2==N
  reps=row['consecutive_representatives']
  assert [r['m']for r in reps]==list(range(N,N+a))
  for r in reps:
   assert r['seed_orders'] and all(q in used for q in r['seed_orders'])
   assert 1+sum(q-1 for q in r['seed_orders'])==r['m']
  out.append({'states':t,'all_orders_from':N,'increment':a,'consecutive_orders_verified':len(reps),'used_seed_count':len(used),'valid':True})
 print(json.dumps(out,indent=2))
 (ROOT/'linear_arithmetic_checks.json').write_text(json.dumps({'all_order_checks':out,'distinct_structural_seed_checks':len(done),'seed_checks':list(done.values())},indent=2)+'\n')
