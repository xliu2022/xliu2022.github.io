"""Independent unrestricted-input graft: no exact-transversal hypothesis."""
from pathlib import Path
from itertools import combinations
import json,sys,time
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT))
sys.path.insert(0,str(ROOT.parent/'prior_v8/verification'))
from finite_slope_constructions import construct
from balanced_extension import six_cyclic,eight_coprime

def graft(M,S,selected):
 n=3*len(M[0]);m=n//3;q=len(S[0]);t=len(M)
 assert len(S)==t and len(selected)==len(set(selected))==3
 old=[v for v in range(n)if v not in selected];oi={v:k for k,v in enumerate(old)};offset=len(old)
 slot={v:j for j,v in enumerate(selected)};out=[]
 for i,row in enumerate(M):
  special=next(e for e in S[i]if e[0]==0)
  N=[tuple(sorted(oi[v]if v in oi else offset+slot[v]*q+special[slot[v]]for v in e))for e in row]
  N.extend(tuple(offset+g*q+e[g]for g in range(3))for e in S[i]if e[0]!=0)
  out.append(N)
 return out

def verify(M,enumerate_all=True):
 n=3*len(M[0]);t=len(M);E=[tuple(sorted(e))for row in M for e in row]
 assert len(E)==len(set(E))==t*n//3
 for row in M:assert sorted(v for e in row for v in e)==list(range(n))
 r={'vertices':n,'m':n//3,'states':t,'edges':len(E)}
 if not enumerate_all:return r
 masks=[sum(1<<v for v in e)for e in E];inc=[[]for _ in range(n)]
 for k,e in enumerate(E):
  for v in e:inc[v].append(k)
 found=set();nodes=0
 def rec(rem,sol):
  nonlocal nodes
  nodes+=1
  if not rem:found.add(frozenset(E[k]for k in sol));return
  opts=None;rr=rem
  while rr:
   bit=rr&-rr;rr-=bit;v=bit.bit_length()-1
   current=[k for k in inc[v]if masks[k]&rem==masks[k]]
   if not current:return
   if opts is None or len(current)<len(opts):opts=current
   if len(opts)==1:break
  for k in opts:rec(rem^masks[k],sol+[k])
 start=time.monotonic();rec((1<<n)-1,[])
 expected={frozenset(tuple(sorted(e))for e in row)for row in M};assert found==expected
 return {**r,'perfect_matchings':len(found),'nodes':nodes,'seconds':time.monotonic()-start,'complete':True}

if __name__=='__main__':
 data=json.loads((ROOT.parent/'prior_v8/prior_results/unrestricted_finite_examples.json').read_text())
 seeds={z['n']:[[tuple(v-1 for v in e)for e in row]for row in z['matchings']]for z in data['constructions']}
 cube=[[(0,b,c),(1,1-b,1-c)]for b in range(2)for c in range(2)]
 specs=[(6,4,cube),(6,6,six_cyclic(3)),(6,8,eight_coprime(5)),(6,10,construct(11,2,[0,1,2,3,5])),(9,12,construct(13,2,[0,1,2,4,5,8]))]
 reports=[]
 for n,t,S in specs:
  M=seeds[n][:t];E={tuple(sorted(e))for row in M for e in row}
  selected_cases=[tuple(M[0][0])]
  nonedge=next((z for z in combinations(range(n),3)if z not in E),None)
  if nonedge is not None:selected_cases.append(nonedge)
  else:selected_cases.append(tuple(M[-1][-1]))
  for selected in selected_cases:
   H=graft(M,S,selected);r=verify(H);r.update(input_n=n,seed_q=len(S[0]),selected_vertices=selected,selected_triple_is_input_edge=tuple(sorted(selected))in E)
   reports.append(r);print(r,flush=True)
 # Repeat unrestricted grafts, choosing all three vertices of an existing edge.
 M=seeds[6][:4]
 for step in range(3):
  selected=tuple(M[0][0]);M=graft(M,cube,selected);r=verify(M);r.update(iteration=step+1,repeated_whole_edge_selection=True);reports.append(r);print(r,flush=True)
 (ROOT/'unrestricted_graft_independent_checks.json').write_text(json.dumps(reports,indent=2)+'\n')
