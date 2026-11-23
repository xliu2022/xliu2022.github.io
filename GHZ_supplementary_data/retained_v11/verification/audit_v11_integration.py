#!/usr/bin/env python3
"""Independent structural and metadata audit of the integrated version-11 constructor.

No complete perfect-matching search is needed: the finite seeds and analytic
operations have separate certificates. This audit targets integration risks.
"""
from pathlib import Path
from collections import Counter
from itertools import combinations
import importlib.util,sys,json,math,hashlib,time
ROOT=Path(__file__).resolve().parent.parent
OUT=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('audit_target_v11',ROOT/'construct_practical.py')
mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
KNOWN={3:(1,[0,1,2]),5:(1,[0,1,2,3]),7:(1,[0,1,2,3]),
 11:(2,[0,1,2,3,5]),13:(2,[0,1,2,4,5,8]),17:(2,[0,1,2,4,5]),
 19:(3,[0,1,3,12,14,15]),23:(2,[0,1,2,9,10,11]),
 29:(8,[0,1,2,7,20,21]),31:(4,[0,1,3,16,21,28]),
 37:(1,[0,1,2,5,21]),41:(1,[0,1,2,5,16]),
 43:(1,[0,1,2,5,38]),47:(1,[0,1,2,11,23])}

def certify_seed_meta(rec,t):
 q,L=rec['base_m'],rec['sheets'];method=rec['method']
 assert isinstance(q,int)and isinstance(L,int)and q>=1 and L>=1
 if method=='retained_latin_odd_lift':
  assert q==8 and L%2==1 and rec['states']==8 and t<=8
  return 8*L
 assert q in KNOWN
 lam,S=KNOWN[q]
 assert rec['slope']==lam and rec['shifts']==S
 assert q%2 and math.gcd(lam,q)==math.gcd(lam+1,q)==1
 if method=='single_sum_seed':
  assert L==1 and rec['states']==len(S) and t<=len(S)
  return q
 assert L>=2 and math.gcd(q,L)==1
 delta=rec['delta'];codes=rec['codes'];drop=set(rec['delete_h1'])
 assert math.gcd(q,delta)==1 and len(codes)==len(S) and len(set(codes))==len(S)
 assert all(isinstance(c,int)and c>=0 for c in codes)
 d=max(codes).bit_length();assert d==0 or (1<<(d-1))-1<q
 assert drop<=set(S) and rec['states']==2*len(S)-len(drop) and rec['states']>=t
 if L==2:
  assert method=='linear_two_sheet_deletion'
  assert set(S)&{(s+delta)%q for s in S}<=drop
 elif L==3:
  assert method=='linear_monotone_three_sheet'
  cd=dict(zip(S,codes))
  for s in S:
   if (s+delta)%q in cd:assert cd[s]&cd[(s+delta)%q]==cd[(s+delta)%q]
 else:
  assert method=='linear_cyclic_lift'
 return q*L

def structural(H,parts,m,t,linear):
 n=3*m;E=[]
 assert len(H)==t and all(len(row)==m for row in H)
 for row in H:
  assert all(len(e)==3 and len(set(e))==3 for e in row)
  assert sorted(v for e in row for v in e)==list(range(n))
  E.extend(tuple(sorted(e))for e in row)
 assert len(E)==m*t and len(set(E))==len(E)
 deg=Counter(v for e in E for v in e)
 assert deg==Counter({v:t for v in range(n)})
 pairs=Counter(p for e in E for p in combinations(e,2));codeg=max(pairs.values())
 if parts is not None:
  assert len(parts)==3 and all(len(P)==m for P in parts)
  flat=[v for P in parts for v in P]
  assert sorted(flat)==list(range(n))
  which={v:i for i,P in enumerate(parts)for v in P}
  assert all(sorted(which[v]for v in e)==[0,1,2]for e in E)
 if linear:assert parts is not None and codeg==1
 return {'m':m,'states':t,'vertices':n,'edges':len(E),'maximum_codegree':codeg,'valid':True}

st=time.monotonic()
old=json.loads((ROOT/'prior_v10/small_m_1_to_1000_v10.json').read_text())
new=json.loads((ROOT/'small_m_1_to_1000_v11.json').read_text())
assert len(old)==len(new)==1000
fd=[];gd=[]
for m,(a,b)in enumerate(zip(old,new),1):
 assert a['m']==b['m']==m
 assert b==mod.bounds(m)
 f=a['unrestricted_f3_at_3m_lower_bound'];g=a['tripartite_g3_lower_bound']
 ff=b['unrestricted_f3_at_3m_lower_bound'];gg=b['tripartite_g3_lower_bound'];ll=b['linear_tripartite_g3_lower_bound']
 assert ff>=f and gg>=g and 1<=ll<=gg<=ff
 assert b['previous_v10_unrestricted_lower_bound']==f
 assert b['previous_v10_tripartite_lower_bound']==g
 assert b['vertices']==3*m and b['unrestricted_edges']==m*ff and b['tripartite_edges']==m*gg and b['linear_tripartite_edges']==m*ll
 if ff>f:fd.append({'m':m,'old':f,'new':ff})
 if gg>g:gd.append({'m':m,'old':g,'new':gg})
 for t,N in [(8,33),(9,62),(10,74),(11,112),(12,113)]:
  if m>=N:assert ll>=t
assert fd==[{'m':8,'old':9,'new':10}],fd
assert gd==[],gd
selected=sorted(set(range(1,61))|{61,62,63,73,74,75,111,112,113,114,150,169,500,999,1000})
checks=[]
for m in selected:
 for setting in ['unrestricted','balanced','linear']:
  H,parts,meta=mod.construct(m,partite=setting=='balanced',linear=setting=='linear')
  row=new[m-1]
  key={'unrestricted':'unrestricted_f3_at_3m_lower_bound','balanced':'tripartite_g3_lower_bound','linear':'linear_tripartite_g3_lower_bound'}[setting]
  t=row[key];r=structural(H,parts,m,t,setting=='linear')
  assert meta['m']==m and meta['states']==t and meta['vertices']==3*m and meta['edges']==m*t
  if setting=='balanced':assert parts is not None
  if setting=='linear':
   if t==1:assert meta['method']=='disjoint_triples'
   elif t==2:assert meta['method']=='linear_two_color_cycle'and m>=3
   else:
    base=certify_seed_meta(meta['base_seed'],t)
    assert base==meta['base_m']
    seeds=meta['seed_orders'];infos=meta['graft_seeds']
    assert len(seeds)==len(infos)==meta['graft_steps']
    assert all(q==certify_seed_meta(rec,t)for q,rec in zip(seeds,infos))
    assert base+sum(q-1 for q in seeds)==m
   r['seed_metadata_valid']=True
  checks.append({'setting':setting,**r})
# The color-deletion interface is checked independently at a large linear output.
H,parts,meta=mod.construct(1000,linear=True,states=10)
color_delete=structural(H,parts,1000,10,True)
assert meta['whole_color_deletion']and meta['maximum_codegree']==1
report={
 'status':'passed','scope':'Independent structural, table-retention, parameter and global-label audit; no complete PM enumeration',
 'table_rows_checked':1000,'unrestricted_changes':fd,'balanced_changes':gd,
 'selected_orders':selected,'constructor_checks':len(checks),'linear_constructor_checks':len(selected),
 'whole_color_deletion_check':color_delete,'checks':checks,
 'source_sha256':{name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest()for name in ['construct_practical.py','linear_cyclic_lift.py','small_m_1_to_1000_v11.json']},
 'seconds':time.monotonic()-st}
(OUT/'v11_integration_independent_audit.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:v for k,v in report.items()if k not in ['checks','source_sha256','selected_orders']},indent=2))
