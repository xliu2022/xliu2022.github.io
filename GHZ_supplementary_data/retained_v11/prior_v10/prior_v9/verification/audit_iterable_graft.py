"""Independent exact-transversal grafting implementation and finite checks."""
from pathlib import Path
from itertools import combinations
from collections import Counter
import json,random,sys
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parent/'prior_v8/verification'))
from balanced_extension import six_cyclic,eight_coprime

def global_seed(S):
 q=len(S[0]);return [[tuple(g*q+e[g]for g in range(3))for e in row]for row in S]

def graft(M,Z,S,selected):
 n=3*len(M[0]);m=n//3;q=len(S[0]);t=len(M)
 assert len(S)==t and len(selected)==len(set(selected))==3
 assert len(Z)==m and not set(selected)&set(Z)
 assert all(len(set(e)&set(Z))==1 for row in M for e in row)
 old=[v for v in range(n)if v not in selected];oi={v:k for k,v in enumerate(old)};offset=len(old)
 slot={v:j for j,v in enumerate(selected)};out=[]
 for i,row in enumerate(M):
  special=next(e for e in S[i]if e[0]==0)
  N=[]
  for e in row:
   N.append(tuple(sorted(oi[v]if v in oi else offset+slot[v]*q+special[slot[v]]for v in e)))
  for e in S[i]:
   if e[0]!=0:N.append(tuple(offset+g*q+e[g]for g in range(3)))
  out.append(N)
 ZZ={oi[z]for z in Z}|{offset+a for a in range(1,q)}
 return out,ZZ

def verify(M,Z,enumerate_all=True):
 n=3*len(M[0]);t=len(M);E=[tuple(sorted(e))for row in M for e in row]
 assert len(E)==len(set(E))==t*n//3
 assert len(Z)==n//3
 assert all(len(set(e)&set(Z))==1 for e in E)
 for row in M:assert sorted(v for e in row for v in e)==list(range(n))
 result={'vertices':n,'m':n//3,'states':t,'edges':len(E),'transversal_size':len(Z),'exact_transversal':True}
 if not enumerate_all:return result
 masks=[sum(1<<v for v in e)for e in E];inc=[[]for _ in range(n)]
 for k,e in enumerate(E):
  for v in e:inc[v].append(k)
 found=set();nodes=0
 def rec(remaining,chosen):
  nonlocal nodes
  nodes+=1
  if not remaining:found.add(frozenset(E[k]for k in chosen));return
  options=None;rem=remaining
  while rem:
   bit=rem&-rem;rem-=bit;v=bit.bit_length()-1
   current=[k for k in inc[v]if masks[k]&remaining==masks[k]]
   if not current:return
   if options is None or len(current)<len(options):options=current
   if len(options)==1:break
  for k in options:rec(remaining^masks[k],chosen+[k])
 rec((1<<n)-1,[])
 expected={frozenset(tuple(sorted(e))for e in row)for row in M}
 assert found==expected
 return {**result,'exact_perfect_matchings':len(found),'nodes':nodes,'complete':True}

if __name__=='__main__':
 random.seed(617351);report={'iterated_tests':[]}
 cube=[[(0,b,c),(1,1-b,1-c)]for b in range(2)for c in range(2)]
 for S,rounds in [(cube,4),(six_cyclic(3),3),(eight_coprime(5),2)]:
  M=global_seed(S);q=len(S[0]);Z=set(range(q))
  for step in range(rounds):
   choices=sorted(set(range(3*len(M[0])))-Z)
   selected=random.sample(choices,3)
   M,Z=graft(M,Z,S,selected)
   r=verify(M,Z);r['iteration']=step+1;r['seed_part_size']=q;r['selected_input_vertices']=selected
   report['iterated_tests'].append(r);print(r,flush=True)
 generators=[6,9,11,14];reachable={0}
 for n in range(1,501):
  if any(n-d in reachable for d in generators):reachable.add(n)
 assert 19 not in reachable and all(n in reachable for n in range(20,501))
 report['numerical_semigroup']={'generators':generators,'conductor':20,'last_gap':19,'six_consecutive_witnesses':{'20':[6,14],'21':[6,6,9],'22':[11,11],'23':[9,14],'24':[6,6,6,6],'25':[11,14]},'analytically_implies_all_n_at_least_20':True}
 (ROOT/'iterable_graft_independent_checks.json').write_text(json.dumps(report,indent=2)+'\n')
