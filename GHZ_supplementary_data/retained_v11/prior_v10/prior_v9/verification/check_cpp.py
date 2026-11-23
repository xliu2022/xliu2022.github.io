"""Optional independent C++ exact-cover certification of all seven finite seeds.

Compile dlx_slope.cpp to dlx_slope.so in this directory before running.
"""
from pathlib import Path
import json
from search_slopes import count
from finite_slope_constructions import SEEDS

records = []
for m,slope,S in SEEDS:
    r = count(m,S,slope,limit=100_000_000)
    assert r['complete'] and r['count'] == len(S), r
    records.append(r)
    print(r,flush=True)
Path(__file__).with_name('finite_slope_cpp_checks.json').write_text(
    json.dumps(records,indent=2)+'\n')
