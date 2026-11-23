"""General-slope modular rigid matching seeds with complete exact-cover certificates.

An entry (m,q,S) means M_{h,s} = {(a,q*a+s,h-(q+1)*a-s):a in Z_m},
h=0,1 and s in S. All arithmetic is modulo m. Requires odd m and q,q+1 units.
The h-sum argument reduces absence of additional PMs to count_single.
"""
from pathlib import Path
from collections import Counter
import argparse,json,math,time
SEEDS=[(11,2,[0,1,2,3,5]),(13,2,[0,1,2,4,5,8]),(17,2,[0,1,2,4,5]),
       (19,3,[0,1,3,12,14,15]),(23,2,[0,1,2,9,10,11]),
       (29,8,[0,1,2,7,20,21]),(31,4,[0,1,3,16,21,28])]

def construct(m,q,S):
 assert m%2 and math.gcd(m,q)==math.gcd(m,q+1)==1
 return [[(a,(q*a+s)%m,(h-(q+1)*a-s)%m)for a in range(m)]for h in (0,1)for s in S]

def verify_structure(M):
 t=len(M);m=len(M[0]);E=[e for row in M for e in row]
 assert len(set(E))==m*t
 for row in M:
  for part in range(3):assert sorted(e[part]for e in row)==list(range(m))
 pairs=Counter((i,j,e[i],e[j])for e in E for i,j in [(0,1),(0,2),(1,2)])
 assert max(pairs.values())<=2
 return dict(m=m,states=t,vertices=3*m,edges=m*t,maximum_codegree=max(pairs.values()))

def exact(M):
 m=len(M[0]);E=sorted(set(e for row in M for e in row));n=3*m
 inc=[[]for _ in range(n)];masks=[]
 for k,e in enumerate(E):
  vv=[e[0],m+e[1],2*m+e[2]];masks.append(sum(1<<v for v in vv))
  for v in vv:inc[v].append(k)
 nodes=0;solutions=[]
 def rec(rem,selected):
  nonlocal nodes
  nodes+=1
  if not rem:solutions.append(frozenset(E[k]for k in selected));return
  r=rem;best=None
  while r:
   bit=r&-r;r-=bit;v=bit.bit_length()-1
   options=[k for k in inc[v]if masks[k]&rem==masks[k]]
   if not options:return
   if best is None or len(options)<len(best):best=options
   if len(best)==1:break
  for k in best:rec(rem^masks[k],selected+[k])
 start=time.monotonic();rec((1<<n)-1,[])
 assert set(solutions)=={frozenset(row)for row in M},'Extra or missing perfect matching'
 return dict(perfect_matchings=len(solutions),search_nodes=nodes,seconds=time.monotonic()-start,complete=True)

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--full',action='store_true');args=ap.parse_args()
 reports=[]
 for m,q,S in SEEDS:
  M=construct(m,q,S);r={**verify_structure(M),'q':q,'S':S,'single_sum':exact(M[:len(S)])}
  if args.full and m<=13:r['full_hypergraph']=exact(M)
  reports.append(r);print(r,flush=True)
  Path(__file__).with_name(f'example_slope_m{m}.json').write_text(json.dumps({'parameters':r,'matchings':M},indent=2)+'\n')
 Path(__file__).with_name('finite_slope_checks.json').write_text(json.dumps(reports,indent=2)+'\n')
