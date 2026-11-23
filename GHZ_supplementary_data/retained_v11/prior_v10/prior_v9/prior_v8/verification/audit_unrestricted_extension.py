"""Independent check of the nonpartite extension, including multiply hit old edges.

The arithmetic at the end checks the preliminary threshold 722. The stronger
threshold 681 is checked in audit_single_unrestricted_seed.py.
"""
from itertools import combinations
from pathlib import Path
import sys,json
sys.path.insert(0,str(Path(__file__).resolve().parent))
from balanced_extension import binary_gadget,six_cyclic

def extend_unrestricted(M,subs=None):
 subs={} if subs is None else subs
 t=len(M);V=set().union(*(set(e) for row in M for e in row));n=len(V)
 assert V==set(range(n))
 X,P,B=binary_gadget(t);L=len(X);assert n>=L
 out=[]
 for i,row in enumerate(M):
  new=[]
  for e in row:
   new.append(tuple(sorted(v if v>=L or (1-B[i][v])==0 else n+v for v in e)))
  fresh=n+3*L
  for x in range(L):
   port=x if B[i][x]==0 else n+x;b0=n+L+x;c0=n+2*L+P[i][x]
   if x not in subs:new.append(tuple(sorted((port,b0,c0))));continue
   seed=subs[x];q=len(seed[0])
   for a,b,c in seed[i]:
    a=port if a==0 else fresh+a-1
    b=b0 if b==0 else fresh+(q-1)+b-1
    c=c0 if c==0 else fresh+2*(q-1)+c-1
    new.append(tuple(sorted((a,b,c))))
   fresh+=3*(q-1)
  out.append(new)
 return out

def exact(M):
 edges=[e for row in M for e in row];assert len(edges)==len(set(edges))
 V=set().union(*map(set,edges));n=len(V);assert V==set(range(n))
 for row in M:
  assert len(row)*3==n and len(set().union(*map(set,row)))==n
 masks=[sum(1<<v for v in e) for e in edges];inc=[[] for _ in range(n)]
 for j,e in enumerate(edges):
  for v in e:inc[v].append(j)
 full=(1<<n)-1;actual=set();nodes=0
 def rec(used,sol):
  nonlocal nodes
  nodes+=1
  if used==full:actual.add(frozenset(edges[j] for j in sol));return
  rem=full^used;opts=None
  while rem:
   bit=rem&-rem;rem-=bit;v=bit.bit_length()-1
   choices=[j for j in inc[v] if not used&masks[j]]
   if not choices:return
   if opts is None or len(choices)<len(opts):opts=choices
   if len(opts)==1:break
  for j in opts:rec(used|masks[j],sol+[j])
 rec(0,[])
 assert actual=={frozenset(row) for row in M}
 return {'vertices':n,'part_size_parameter':n//3,'states':len(M),'exact_matchings':len(actual),'nodes':nodes}

if __name__=='__main__':
 M=[]
 for e in combinations(range(1,6),2):
  edge=tuple(sorted((0,)+e));comp=tuple(sorted(set(range(6))-set(edge)));M.append([edge,comp])
 M=M[:3];assert any(len(set(e)&set(range(4)))>=2 for row in M for e in row)
 results={'small_complete_examples':[]}
 for subs in ({},{0:six_cyclic(3)[:3]},{0:six_cyclic(3)[:3],1:six_cyclic(4)[:3]}):
  H=extend_unrestricted(M,subs);results['small_complete_examples'].append(exact(H))
 prior=lambda m:m>=2 and (m%3==0 or m%4==2)
 def covered(m):
  return prior(m) or any(m-d>=171 and prior(m-d) for d in (512,534,556))
 uncovered=[m for m in range(171,10001) if not covered(m)]
 assert max(uncovered)==721
 assert all(covered(m) for m in range(722,10001))
 starts={r:next(m for m in range(171,1000) if m%12==r and covered(m)) for r in [1,4,5,7,8,11]}
 results['residue_first_covered']=starts
 results['all_order_threshold']=722;results['last_uncovered_by_these_rules']=721
 results['finite_arithmetic_checked_through']=10000
 Path(__file__).with_name('unrestricted_extension_audit.json').write_text(json.dumps(results,indent=2)+'\n')
 print(json.dumps(results,indent=2))
