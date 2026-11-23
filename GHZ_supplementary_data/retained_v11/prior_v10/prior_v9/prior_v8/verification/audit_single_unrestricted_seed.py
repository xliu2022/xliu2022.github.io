"""Independent checks for one arbitrary nonpartite seed substituted into a port row."""
from pathlib import Path
from itertools import combinations
import sys,json
sys.path.insert(0,str(Path(__file__).parent.resolve()))
from audit_unrestricted_extension import exact
from balanced_extension import binary_gadget

def complementary_seed(t):
 rows=[]
 for pair in combinations(range(1,6),2):
  e=(0,)+pair;f=tuple(sorted(set(range(6))-set(e)));rows.append([e,f])
 return rows[:t]

def substitute(M,S,x0=0,distinguished=(0,1)):
 t=len(M);n=3*len(M[0]);q=len(S[0]);X,P,B=binary_gadget(t);L=len(X)
 assert n>=L and len(S)==t
 a0,c0=distinguished;seedV=set(range(3*q));others=sorted(seedV-{a0,c0})
 base_n=n+3*L;oldBx=n+L+x0
 fresh_labels=[oldBx]+list(range(base_n,base_n+3*(q-1)))
 assert len(fresh_labels)==3*q-2
 fresh=dict(zip(others,fresh_labels))
 out=[]
 for i,row in enumerate(M):
  new=[tuple(sorted(v if v>=L or B[i][v]==1 else n+v for v in e)) for e in row]
  for x in range(L):
   p=x if B[i][x]==0 else n+x;bb=n+L+x;cc=n+2*L+P[i][x]
   if x!=x0:new.append(tuple(sorted((p,bb,cc))))
   else:
    mapping={**fresh,a0:p,c0:cc}
    new.extend(tuple(sorted(mapping[v] for v in e)) for e in S[i])
  out.append(new)
 return out

if __name__=='__main__':
 results=[]
 for t in [2,3]:
  M=complementary_seed(t);S=complementary_seed(t)
  for pair in [(0,1),(0,3),(1,5)]:
   H=substitute(M,S,distinguished=pair)
   result={'states':t,'distinguished':pair,**exact(H)};results.append(result)
 prior=lambda m:m>=2 and (m%3==0 or m%4==2)
 def covered(m):return prior(m) or any(m-d>=171 and (m-d)%3==0 for d in (512,514))
 missing=[m for m in range(171,10001) if not covered(m)]
 assert max(missing)==680 and all(covered(m) for m in range(681,10001))
 report={'exact_checks':results,'all_order_threshold':681,'last_uncovered_by_these_rules':680,'arithmetic_checked_through':10000,'unrestricted_seed_substitutions':1}
 Path(__file__).with_name('single_unrestricted_seed_audit.json').write_text(json.dumps(report,indent=2)+'\n')
 print(json.dumps(report,indent=2))
