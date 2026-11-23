"""Independent exact checks; the extension theorem has a separate analytic proof."""
from pathlib import Path
from itertools import product
from collections import Counter
import json,time

def gadget(t):
 X=[x for x in product((0,1),repeat=t) if sum(x)%2==0]
 def prefix(x,k):return sum(x[:k])%2
 def rank(x,i):return sum((prefix(x,(i+a)%t)^prefix(x,i))*2**(a-1) for a in range(1,t))
 checks=0
 for i in range(t):
  for x in X:
   for j in range(t):
    if i==j or prefix(x,i)!=prefix(x,j):continue
    y=list(x);y[i]^=1;y[j]^=1;y=tuple(y)
    assert rank(y,i)>rank(x,i)
    checks+=1
 return {'t':t,'rows':len(X),'rank_inequalities':checks}

def single_sum(m,S):
 edges=[];incident=[[] for _ in range(3*m)]
 for a in range(m):
  for s in S:
   e=(a,m+(a+s)%m,2*m+(-2*a-s)%m);mask=sum(1<<v for v in e)
   edges.append((e,mask))
   for v in e:incident[v].append(len(edges)-1)
 full=(1<<(3*m))-1;count=0;nodes=0;solutions=[]
 def rec(used,selected):
  nonlocal count,nodes
  nodes+=1
  if used==full:
   count+=1;solutions.append(selected);return
  remaining=full^used;best=None
  while remaining:
   bit=remaining&-remaining;remaining-=bit;v=bit.bit_length()-1
   opts=[j for j in incident[v] if not used&edges[j][1]]
   if not opts:return
   if best is None or len(opts)<len(best):best=opts
   if len(best)==1:break
  for j in best:rec(used|edges[j][1],selected+[j])
 start=time.monotonic();rec(0,[])
 uniform=[]
 for sol in solutions:
  shifts={(edges[j][0][1]-m-edges[j][0][0])%m for j in sol}
  assert len(shifts)==1
  uniform.append(next(iter(shifts)))
 assert sorted(uniform)==sorted(S)
 return {'m':m,'S':S,'exact_matchings':count,'nodes':nodes,'seconds':time.monotonic()-start,'only_constant_shifts':True}

if __name__=='__main__':
 result={'gadgets':[gadget(t) for t in range(2,11)],'single_sum_seeds':[]}
 for m,S in [(23,[0,1,2,5,11]),(29,[0,1,2,5,17]),(29,[0,1,3,7,20])]:
  r=single_sum(m,S);result['single_sum_seeds'].append(r);print(r,flush=True)
 for t,threshold,L in [(7,135,64),(8,265,128)]:
  records=[]
  for m in range(threshold,1001):
   if m>=5 and m%6 in (1,5):continue
   poss=[(inc,m-inc) for inc in (L,L+7,L+11) if m-inc>=L and (m-inc)%6 in (1,5)]
   assert poss,(t,m)
   records.append({'m':m,'increment':poss[0][0],'base':poss[0][1]})
  result[f'coverage_t{t}']={'threshold':threshold,'through':1000,'noncoprime_orders':len(records),'verified':True}
 Path(__file__).with_name('independent_audit.json').write_text(json.dumps(result,indent=2)+'\n')
 print(json.dumps(result,indent=2))
