"""Complete independent exact-cover verification of a 21-vertex ten-state example."""
from pathlib import Path
import json
from collections import Counter

A_ROWS=(
(('123','3'),('456',''),('789','')),
(('128','8'),('369',''),('457','4')),
(('129',''),('356',''),('478','')),
(('134','34'),('279',''),('568','')),
(('137','3'),('246','6'),('589','')),
(('145','4'),('238',''),('679','67')),
(('147','4'),('269','26'),('358','3')),
(('149','4'),('268','6'),('357','3')),
(('159','9'),('248',''),('367','')),
(('167','6'),('235','3'),('489','8')),
(('178','8'),('259',''),('346','46')))
B_PARTITIONS=[
[10,9,8,6,4,7,3,5,1,11,2,0],
[2,1,11,4,9,6,7,8,5,10,3,0],
[1,0,4,8,6,2,11,3,10,9,7,5],
[1,8,2,4,10,11,7,6,9,0,5,3],
[10,0,7,3,4,8,1,9,6,2,5,11],
[6,11,1,2,4,9,0,5,10,3,7,8],
[9,10,3,6,2,11,1,7,0,5,4,8],
[10,0,6,2,1,4,11,5,8,7,9,3],
[1,6,9,0,7,8,5,3,11,2,10,4],
[6,4,3,8,5,0,2,1,9,10,11,7]]

def build():
 rows=[row for i,row in enumerate(A_ROWS)if i!=4];M=[]
 for arow,brow in zip(rows,B_PARTITIONS):
  S={int(v)-1 for e,o in arow for v in o};R=[]
  for e,o in arow:R.append(tuple(sorted((9 if v in o else 0)+int(v)-1 for v in e)))
  for j in range(0,12,3):R.append(tuple(sorted(v if v in S else 9+v for v in brow[j:j+3])))
  assert Counter(v for e in R for v in e)==Counter(range(21));M.append(R)
 E=[e for row in M for e in row];assert len(E)==len(set(E))==70
 assert Counter(v for e in E for v in e)==Counter({v:10 for v in range(21)})
 return M

def verify(M):
 E=[e for row in M for e in row];masks=[sum(1<<v for v in e)for e in E]
 by_vertex=[[j for j,e in enumerate(E)if v in e]for v in range(21)]
 actual=[];nodes=0;branches=[]
 def rec(rem,chosen):
  nonlocal nodes
  nodes+=1
  if not rem:actual.append(frozenset(chosen));return
  v=(rem&-rem).bit_length()-1
  for j in by_vertex[v]:
   e=masks[j]
   if rem&e==e:
    before_nodes=nodes;before_solutions=len(actual)
    rec(rem^e,chosen+[j])
    if not chosen:branches.append({'first_edge':E[j],'nodes':nodes-before_nodes,'solutions':len(actual)-before_solutions})
 rec((1<<21)-1,[])
 expected={frozenset(range(7*i,7*(i+1)))for i in range(10)}
 assert set(actual)==expected and len(actual)==10
 return {'vertices':21,'part_size_parameter':7,'edges':70,'states':10,'perfect_matchings':10,'complete':True,'search_nodes':nodes,'first_vertex_branches':branches}

if __name__=='__main__':
 M=build();r=verify(M);path=Path(__file__).with_name('twenty_one_ten_certificate.json')
 path.write_text(json.dumps({'verification':r,'matchings':M},indent=2)+'\n');print(json.dumps(r,indent=2))
