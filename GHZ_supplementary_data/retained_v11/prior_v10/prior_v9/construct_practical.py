#!/usr/bin/env python3
"""Construct the strongest bound recorded in the version-9 table.
Default: unrestricted f3 construction. --tripartite: balanced g3 construction.
Every exported vertex label is global, in 0,...,3m-1.
"""
from pathlib import Path
import argparse, importlib.util, json, sys
from functools import lru_cache
from grafting import globalize, graft, check
ROOT = Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'verification'))
from finite_slope_constructions import construct as modular
from verify_twenty_one import build as twenty_one
spec=importlib.util.spec_from_file_location('prior_v8',ROOT/'prior_v8/construct_small_order.py')
prior=importlib.util.module_from_spec(spec);spec.loader.exec_module(prior)
OLD_TABLE=json.loads((ROOT/'prior_v8/small_m_1_to_1000.json').read_text())
NEW_SEEDS={11:(2,(0,1,2,3,5)),13:(2,(0,1,2,4,5,8)),17:(2,(0,1,2,4,5)),
 19:(3,(0,1,3,12,14,15)),23:(2,(0,1,2,9,10,11)),
 29:(8,(0,1,2,7,20,21)),31:(4,(0,1,3,16,21,28))}
SEED_ORDERS={12:(13,19,23,29,31),11:(13,19,23,29,31),
 10:tuple(sorted(set(NEW_SEEDS)|set(prior.TEN))),
 8:tuple(sorted({5,7,8}|set(NEW_SEEDS)|set(prior.TEN)))}

def seed(q,t):
 if q in NEW_SEEDS:
  slope,S=NEW_SEEDS[q];M=modular(q,slope,S)
 else:M,_=prior.construct(q)
 assert len(M)>=t
 return M[:t]

def old_record(m):
 if m<=1000:
  r=OLD_TABLE[m-1]
  return r['tripartite_g3_lower_bound'],r['unrestricted_f3_at_3m_lower_bound']
 g=len(prior.construct(m)[0]);return g,prior.unrestricted_bound(m,g)[0]

@lru_cache(maxsize=20)
def plans(limit,t,partite):
 """Shortest graft plans; finite arithmetic only, no PM search."""
 paths={}
 for m in range(1,limit+1):
  old_g,old_f=old_record(m)
  if (old_g if partite else old_f)>=t:
   paths[m]=(m,(),'prior_partite' if partite else 'prior_unrestricted')
  if m in NEW_SEEDS and 2*len(NEW_SEEDS[m][1])>=t:
   paths[m]=(m,(),'new_partite_seed')
  if not partite and m==7 and t<=10:paths[m]=(m,(),'twenty_one')
 for m in range(1,limit+1):
  if m not in paths:continue
  base,steps,kind=paths[m]
  for q in SEED_ORDERS[t]:
   nxt=m+q-1
   if nxt<=limit and (nxt not in paths or len(paths[nxt][1])>len(steps)+1):
    paths[nxt]=(base,steps+(q,),kind)
 return paths

def bounds(m,limit=1000):
 if m<1:raise ValueError('m must be positive')
 limit=max(limit,m);g,f=old_record(m)
 for t in (8,10,12):
  if m in plans(limit,t,True):g=max(g,t)
 for t in (10,11,12):
  if m in plans(limit,t,False):f=max(f,t)
 old_g,old_f=old_record(m)
 return {'m':m,'vertices':3*m,'tripartite_g3_lower_bound':g,
         'tripartite_edges':m*g,'unrestricted_f3_at_3m_lower_bound':f,
         'unrestricted_edges':m*f,'previous_v8_tripartite_lower_bound':old_g,
         'previous_v8_unrestricted_lower_bound':old_f}

def construct(m,partite=False):
 r=bounds(m);t=r['tripartite_g3_lower_bound' if partite else 'unrestricted_f3_at_3m_lower_bound']
 if t in SEED_ORDERS and m in plans(max(m,1000),t,partite):
  base,steps,kind=plans(max(m,1000),t,partite)[m]
  if kind=='new_partite_seed':H,parts=globalize(seed(base,t))
  elif kind=='prior_partite':
   M,_=prior.construct(base);H,parts=globalize(M[:t])
  elif kind=='twenty_one':H,parts=twenty_one()[:t],None
  else:
   H,_=prior.construct_unrestricted(base);H=H[:t];parts=None
  for q in steps:H,parts=graft(H,seed(q,t),parts)
  return H,parts,{**check(H,parts),'setting':'balanced_tripartite' if partite else 'unrestricted',
         'method':'universal_grafting' if steps else kind,'base_m':base,'base_kind':kind,
         'seed_orders':list(steps),'graft_steps':len(steps),'tripartition_available':parts is not None}
 if partite:
  M,meta=prior.construct(m);H,parts=globalize(M)
 else:H,meta=prior.construct_unrestricted(m);parts=None
 assert len(H)==t
 return H,parts,{**check(H,parts),**meta,'setting':'balanced_tripartite' if partite else 'unrestricted'}

def export(m,path,partite=False):
 H,parts,meta=construct(m,partite)
 data={'parameters':meta,'vertex_labels':'0,...,3m-1','matchings':H}
 if parts is not None:data['parts']=[sorted(P)for P in parts]
 path.parent.mkdir(parents=True,exist_ok=True);path.write_text(json.dumps(data,indent=2)+'\n')
 return meta

def main():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('--m',type=int);p.add_argument('--tripartite',action='store_true')
 p.add_argument('--export',type=Path);p.add_argument('--table',action='store_true');args=p.parse_args()
 if args.table:
  R=[bounds(m)for m in range(1,1001)]
  for r in R:
   m=r['m'];g=r['tripartite_g3_lower_bound'];f=r['unrestricted_f3_at_3m_lower_bound']
   assert g>=r['previous_v8_tripartite_lower_bound'] and f>=r['previous_v8_unrestricted_lower_bound']
   if m>=11:assert g>=8
   if m%2 and m>=17:assert g>=10
   if m%2 and m>=41:assert g>=12
   for threshold,t in [(9,10),(21,11),(43,12)]:
    if m>=threshold:assert f>=t
   if m%2 and m>=19:assert f>=12
  (ROOT/'small_m_1_to_1000_v9.json').write_text(json.dumps(R,indent=2)+'\n')
  changes={'tripartite':sum(r['tripartite_g3_lower_bound']>r['previous_v8_tripartite_lower_bound']for r in R),
           'unrestricted':sum(r['unrestricted_f3_at_3m_lower_bound']>r['previous_v8_unrestricted_lower_bound']for r in R)}
  print(json.dumps({'rows':len(R),'improved_entries':changes}))
 if args.m is not None:
  meta=export(args.m,args.export,args.tripartite)if args.export else construct(args.m,args.tripartite)[2]
  print(json.dumps(meta))
 if not args.table and args.m is None:p.print_help()
if __name__=='__main__':main()
