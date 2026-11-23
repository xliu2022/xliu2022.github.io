from math import gcd
from pathlib import Path
import json,heapq
ROOT=Path(__file__).resolve().parent
MAX=5000
from linear_cyclic_lift import SEEDS
BASES={q:(2*len(S),S)for q,(lam,S)in SEEDS.items()}

def bank(t):
 B={}
 if t<=8:
  for l in range(1,MAX//8+1,2):B[8*l]={'base_m':8,'sheets':l,'method':'old_linear_latin_lift','states':8}
 for q,(states,S) in BASES.items():
  if states<t:continue
  overlap,delta,drop=min((len(set(S)&{(s+d)%q for s in S}),d,sorted(set(S)&{(s+d)%q for s in S}))for d in range(1,q)if gcd(d,q)==1)
  for L in range(2,MAX//q+1):
   if gcd(q,L)!=1:continue
   if L==2:
    if states-overlap<t:continue
    rec={'base_m':q,'sheets':L,'method':'color_deleted_linear_two_lift','states':states-overlap,'delta':delta,'delete_h1':drop}
   elif L==3:rec={'base_m':q,'sheets':L,'method':'monotone_code_linear_three_lift','states':states,'delta':1}
   else:rec={'base_m':q,'sheets':L,'method':'linear_cyclic_lift','states':states,'delta':1}
   if q*L not in B or rec['states']>B[q*L]['states']:B[q*L]=rec
 return B

def compute(t):
 B=bank(t); inc=sorted(q-1 for q in B);a=inc[0]
 dist=[10**9]*a;path=[None]*a;dist[0]=0;path[0]=[];pq=[(0,0)]
 while pq:
  val,r=heapq.heappop(pq)
  if val!=dist[r]:continue
  for x in inc:
   rr=(r+x)%a;v=val+x
   if v<dist[rr]:dist[rr]=v;path[rr]=path[r]+[x];heapq.heappush(pq,(v,rr))
 F=max(dist)-a;threshold=F+2
 records=[{'residue':r,'value':dist[r],'increments':path[r]}for r in range(a)]
 used=sorted({a+1}|{x+1 for pp in path for x in pp})
 # Direct representatives in the first full block; their displays alone prove coverage.
 reps=[]
 for m in range(threshold,threshold+a):
  s=m-1;r=s%a;steps=path[r]+[a]*((s-dist[r])//a)
  assert sum(steps)==s and steps
  reps.append({'m':m,'seed_orders':[x+1 for x in steps]})
 out={'states':t,'min_increment':a,'frobenius':F,'all_orders_from':threshold,'apery':records,'used_seeds':{q:B[q]for q in used},'consecutive_representatives':reps}
 print('RESULT',t,'threshold',threshold,'used',used)
 print('APERY',dist)
 print('REPRESENTATIVES',reps)
 return out
if __name__=='__main__':
 out=[compute(t)for t in [8,9,10,11,12]]
 (ROOT/'linear_thresholds_refined.json').write_text(json.dumps(out,indent=2)+'\n')
