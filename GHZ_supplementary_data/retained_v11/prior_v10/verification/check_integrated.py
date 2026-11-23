#!/usr/bin/env python3
"""Structural integration checks; rigidity follows from the proved constructions."""
import json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parent.parent
sys.path.insert(0,str(ROOT))
from construct_practical import construct,bounds,export
records=[]
orders=list(range(1,61))+[74,100,135,263,500,681,999,1000]
for m in orders:
 for partite in(False,True):
  H,parts,meta=construct(m,partite)
  r=bounds(m)
  assert meta['states']==r['tripartite_g3_lower_bound'if partite else'unrestricted_f3_at_3m_lower_bound']
  assert meta['vertices']==3*m
  if partite and m>=7:assert parts is not None and meta['maximum_codegree']<=2
  records.append(meta)
(ROOT/'verification/constructor_integration_checks.json').write_text(json.dumps(records,indent=2)+'\n')
for m in[6,10,18,21,22,24,26,30,31,35,38,39,40,42,43,44,52,100,999,1000]:
 export(m,ROOT/'examples'/f'unrestricted_m{m}.json',False)
for m in[10,22,26,31,35,38,39,40,42,43,44,52,100,999,1000]:
 export(m,ROOT/'examples'/f'tripartite_m{m}.json',True)
print(json.dumps({'checked_input_setting_pairs':len(records),'examples_exported':35}))
