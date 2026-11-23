#!/usr/bin/env python3
"""Version 10 practical rigid-matching constructors (standard library only).

Use --tripartite for the balanced function g3(m). Default realizes f3(3m),
prefering a balanced output when its best recorded state count is the same.
"""
from pathlib import Path
from functools import lru_cache
from collections import Counter
from itertools import combinations
import argparse, importlib.util, json, math, sys
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'prior_v9'))
sys.path.insert(0,str(ROOT/'verification'))
from grafting import globalize,graft,check
from cyclic_modular_lift import construct as cyclic
from verify_twelve_on_eighteen import construct as eighteen
spec=importlib.util.spec_from_file_location('retained_v9',ROOT/'prior_v9/construct_practical.py')
v9=importlib.util.module_from_spec(spec);spec.loader.exec_module(v9)
OLD=json.loads((ROOT/'prior_v9/small_m_1_to_1000_v9.json').read_text())
BASES={5:(1,(0,1,2,3)),7:(1,(0,1,2,3)),**v9.NEW_SEEDS}
for m,S in v9.prior.TEN.items():
 if m not in BASES:BASES[m]=(1,S)

def prior_record(m):
 return OLD[m-1]if m<=1000 else v9.bounds(m,m)

@lru_cache(maxsize=4)
def bank(limit):
    """Finite collection of analytic seeds, indexed by their exact part size."""
    B={}
    for m,(lam,S)in sorted(BASES.items()):
        for L in range(1,limit//m+1):
            if math.gcd(L,m)!=1:continue
            q=m*L;rec={'states':2*len(S),'method':'cyclic_modular_seed','base_m':m,
                'slope':lam,'shifts':list(S),'sheets':L}
            if q not in B or (rec['states'],-L)>(B[q]['states'],-B[q].get('sheets',1)):
                B[q]=rec
    if limit>=8:B[8]={'states':8,'method':'retained_latin_seed','base_m':8,'sheets':1}
    return B

def seed(q,t,limit):
    info=bank(limit)[q]
    assert info['states']>=t
    if info['method']=='retained_latin_seed':M,_=v9.prior.construct(8)
    else:M=cyclic(info['base_m'],info['slope'],info['shifts'],info['sheets'])
    return M[:t]

@lru_cache(maxsize=12)
def plans(limit,t,partite):
    B=bank(limit);paths={}
    if not partite:
        for m in range(1,limit+1):
            if prior_record(m)['unrestricted_f3_at_3m_lower_bound']>=t:
                paths[m]=(m,(),'prior_unrestricted')
        if t<=12 and limit>=6:paths[6]=(6,(),'eighteen_vertex_seed')
    for q,meta in B.items():
        if meta['states']>=t:paths[q]=(q,(),'balanced_seed')
    qs=sorted(q for q,meta in B.items()if meta['states']>=t)
    for m in range(1,limit+1):
        if m not in paths:continue
        base,steps,kind=paths[m]
        for q in qs:
            nxt=m+q-1
            if nxt>limit:break
            if nxt not in paths or len(paths[nxt][1])>len(steps)+1:
                paths[nxt]=(base,steps+(q,),kind)
    return paths

def bounds(m,limit=1000):
    if m<1:raise ValueError('m must be positive')
    limit=max(m,limit);old=prior_record(m)
    g=old['tripartite_g3_lower_bound'];f=old['unrestricted_f3_at_3m_lower_bound']
    for t in(8,10,12):
        if m in plans(limit,t,True):g=max(g,t)
    for t in(10,11,12):
        if m in plans(limit,t,False):f=max(f,t)
    f=max(f,g)
    return {'m':m,'vertices':3*m,'tripartite_g3_lower_bound':g,'tripartite_edges':m*g,
        'unrestricted_f3_at_3m_lower_bound':f,'unrestricted_edges':m*f,
        'previous_v9_tripartite_lower_bound':old['tripartite_g3_lower_bound'],
        'previous_v9_unrestricted_lower_bound':old['unrestricted_f3_at_3m_lower_bound']}

def max_codegree(H):
    return max(Counter(tuple(sorted(p))for row in H for e in row for p in combinations(e,2)).values())

def construct(m,partite=False):
    limit=max(m,1000);r=bounds(m,limit)
    t=r['tripartite_g3_lower_bound'if partite else'unrestricted_f3_at_3m_lower_bound']
    use_balanced=partite or r['tripartite_g3_lower_bound']==t
    if t in(8,10,11,12)and m in plans(limit,t,use_balanced):
        base,steps,kind=plans(limit,t,use_balanced)[m]
        if kind=='balanced_seed':H,parts=globalize(seed(base,t,limit))
        elif kind=='eighteen_vertex_seed':H,parts=eighteen()[:t],None
        else:H,parts,_=v9.construct(base,False);H=H[:t]
        for q in steps:H,parts=graft(H,seed(q,t,limit),parts)
        meta={**check(H,parts),'method':'universal_grafting'if steps else kind,
              'base_m':base,'base_kind':kind,'seed_orders':list(steps),
              'graft_steps':len(steps),'tripartition_available':parts is not None}
        if kind=='balanced_seed':meta['base_seed']=bank(limit)[base]
        meta['graft_seeds']=[bank(limit)[q]for q in steps]
    else:H,parts,meta=v9.construct(m,partite)
    assert len(H)==t
    meta={**meta,'setting':'balanced_tripartite'if partite else'unrestricted',
          'maximum_codegree':max_codegree(H)}
    if use_balanced and t>=8:
        assert parts is not None and meta['maximum_codegree']<=2
    return H,parts,meta

def export(m,path,partite=False):
    H,parts,meta=construct(m,partite)
    data={'parameters':meta,'vertex_labels':'0,...,3m-1','matchings':H}
    if parts is not None:data['parts']=[sorted(P)for P in parts]
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(data,indent=2)+'\n');return meta

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--m',type=int);p.add_argument('--tripartite',action='store_true')
    p.add_argument('--export',type=Path);p.add_argument('--table',action='store_true')
    a=p.parse_args()
    if a.table:
        R=[bounds(m)for m in range(1,1001)]
        for r in R:
            m=r['m'];g=r['tripartite_g3_lower_bound'];f=r['unrestricted_f3_at_3m_lower_bound']
            for threshold,t in [(7,8),(31,10),(43,12)]:
                if m>=threshold:assert g>=t
            if m%2 and m>=17:assert g>=10
            if m%2 and m>=35:assert g>=12
            if m>=9:assert f>=10
            if m>=21 or (m%2 and m>=19):assert f>=12
            assert g>=r['previous_v9_tripartite_lower_bound']and f>=r['previous_v9_unrestricted_lower_bound']
        (ROOT/'small_m_1_to_1000_v10.json').write_text(json.dumps(R,indent=2)+'\n')
        print(json.dumps({'rows':len(R),'improved_balanced':sum(r['tripartite_g3_lower_bound']>r['previous_v9_tripartite_lower_bound']for r in R),
            'improved_unrestricted':sum(r['unrestricted_f3_at_3m_lower_bound']>r['previous_v9_unrestricted_lower_bound']for r in R)}))
    if a.m is not None:
        print(json.dumps(export(a.m,a.export,a.tripartite)if a.export else construct(a.m,a.tripartite)[2]))
    if not a.table and a.m is None:p.print_help()
if __name__=='__main__':main()
