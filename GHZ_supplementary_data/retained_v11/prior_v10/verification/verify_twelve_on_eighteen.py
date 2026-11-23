#!/usr/bin/env python3
"""Complete finite verification of a twelve-state 3-graph on 18 vertices."""
import json
from pathlib import Path
from collections import Counter

SPLIT_EDGES = [tuple(map(int,s.split())) for s in '''
7 4
56 0
448 0
131 128
292 32
88 80
259 0
52 0
200 0
13 12
322 0
176 0
69 4
42 32
400 0
25 8
134 0
352 256
73 8
290 256
148 128
265 8
162 32
84 80
273 256
138 0
100 0
97 0
22 4
392 128
193 128
274 0
44 8
161 0
82 0
268 8
'''.splitlines() if s.strip()]

def construct(q=2):
 return [[tuple(9*((j+((o>>v)&1))%q)+v for v in range(9) if e>>v&1)
          for j in range(q) for e,o in SPLIT_EDGES[3*i:3*i+3]] for i in range(12)]

def verify():
 M=construct(); n=18; full=(1<<n)-1
 edges=[tuple(sorted(e)) for row in M for e in row]
 assert len(edges)==len(set(edges))==72
 assert all(Counter(v for e in row for v in e)==Counter(range(n)) for row in M)
 masks=[sum(1<<v for v in e) for e in edges]
 incident=[[i for i,e in enumerate(masks) if e>>v&1] for v in range(n)]
 sols=[]; nodes=0; branches=[]
 def go(used,picked):
  nonlocal nodes
  nodes+=1
  if used==full:sols.append(tuple(picked));return
  rest=full^used;v=(rest&-rest).bit_length()-1
  for i in incident[v]:
   if not masks[i]&used:go(used|masks[i],picked+[i])
 nodes=1
 for i in incident[0]:
  oldn=nodes; olds=len(sols);go(masks[i],[i]);branches.append({'first_edge':edges[i],'nodes':nodes-oldn,'solutions':len(sols)-olds})
 expected={frozenset(range(6*i,6*i+6)) for i in range(12)}
 assert {frozenset(s) for s in sols}==expected and len(sols)==12
 return {'n':18,'m':6,'states':12,'edges':72,'complete':True,'nodes':nodes,'branches':branches,'matchings':M,'split_edges':SPLIT_EDGES}
if __name__=='__main__':
 result=verify();Path(__file__).with_name('twelve_on_eighteen_certificate.json').write_text(json.dumps(result,indent=2)+'\n')
 print(json.dumps({k:v for k,v in result.items() if k not in ('matchings','split_edges')},indent=2))
