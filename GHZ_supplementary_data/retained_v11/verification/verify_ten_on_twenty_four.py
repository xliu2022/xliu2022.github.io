#!/usr/bin/env python3
"""Complete finite certificate: exactly ten perfect matchings on 24 vertices."""
import json
from pathlib import Path
from collections import Counter
SPLIT_EDGES=[(2084, 2080), (522, 8), (208, 128), (1281, 256), (13, 0), (642, 640), (1312, 1024), (2128, 2048), (22, 4), (1096, 0), (161, 32), (2816, 512), (37, 32), (2184, 2048), (1600, 0), (274, 256), (772, 512), (1042, 0), (2113, 2112), (168, 0), (560, 32), (3200, 1024), (321, 64), (14, 0), (1544, 0), (292, 0), (2114, 2048), (145, 16), (392, 384), (112, 0), (2054, 2048), (1537, 1536), (1282, 1024), (580, 576), (2177, 0), (56, 0), (265, 8), (1092, 1088), (162, 160), (2576, 0)]

def construct():
 return [[tuple(12*((j+((o>>v)&1))%2)+v for v in range(12) if e>>v&1)
          for j in range(2) for e,o in SPLIT_EDGES[4*c:4*c+4]]for c in range(10)]

def verify():
 M=construct();E=[tuple(sorted(e))for row in M for e in row];assert len(E)==len(set(E))==80
 assert all(Counter(v for e in row for v in e)==Counter(range(24))for row in M)
 masks=[sum(1<<v for v in e)for e in E];inc=[[i for i,e in enumerate(masks)if e>>v&1]for v in range(24)]
 assert all(len(x)==10 for x in inc)
 full=(1<<24)-1;sols=[];nodes=1;branches=[]
 def visit(used,chosen):
  nonlocal nodes
  nodes+=1
  if used==full:sols.append(chosen);return
  r=full^used;v=(r&-r).bit_length()-1
  for i in inc[v]:
   if not used&masks[i]:visit(used|masks[i],chosen+[i])
 for i in inc[0]:
  n0=nodes;s0=len(sols);visit(masks[i],[i]);branches.append({'first_edge':E[i],'nodes':nodes-n0,'solutions':len(sols)-s0})
 assert len(sols)==10
 assert {frozenset(s)for s in sols}=={frozenset(range(8*i,8*i+8))for i in range(10)}
 return {'n':24,'m':8,'states':10,'edges':80,'nodes':nodes,'complete':True,'branches':branches,'matchings':M,'split_edges':SPLIT_EDGES}
if __name__=='__main__':
 r=verify();Path(__file__).with_name('ten_on_twenty_four_certificate.json').write_text(json.dumps(r,indent=2)+'\n')
 print(json.dumps({k:v for k,v in r.items()if k not in ('matchings','split_edges')},indent=2))
