#!/usr/bin/env python3
"""Version 11: unrestricted, balanced, and linear balanced rigid 3-graphs.

Use --linear for a balanced construction with maximum pair codegree one.
Use --tripartite for the strongest retained balanced bound without linearity.
The default realizes the unrestricted bound on exactly 3*m vertices.
"""
from pathlib import Path
from functools import lru_cache
import argparse, importlib.util, json, math, sys

ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'prior_v10'))
spec=importlib.util.spec_from_file_location('retained_v10',ROOT/'prior_v10/construct_practical.py')
v10=importlib.util.module_from_spec(spec);spec.loader.exec_module(v10)
sys.path.insert(0,str(ROOT))
sys.path.insert(0,str(ROOT/'verification'))
from grafting import globalize,graft,check
from verify_ten_on_twenty_four import construct as twenty_four
from linear_cyclic_lift import SEEDS,construct as linear_lift,monotone_codes,best_two_delta
from balanced_extension import eight_latin_lift

OLD=json.loads((ROOT/'prior_v10/small_m_1_to_1000_v10.json').read_text())

@lru_cache(maxsize=4)
def linear_bank(limit):
    B={}
    def add(q,rec):
        if q not in B or rec['states']>B[q]['states']:B[q]=rec
    for L in range(1,limit//8+1,2):
        add(8*L,{'states':8,'method':'retained_latin_odd_lift','base_m':8,'sheets':L})
    singles={3:(1,[0,1,2]),**SEEDS}
    for q,(lam,S) in singles.items():
        if q<=limit:add(q,{'states':len(S),'method':'single_sum_seed','base_m':q,
                          'slope':lam,'shifts':list(S),'sheets':1})
    for q,(lam,S) in SEEDS.items():
        overlap,delta,drop=best_two_delta(q,S)
        for L in range(2,limit//q+1):
            if math.gcd(q,L)!=1:continue
            rec={'states':2*len(S),'method':'linear_cyclic_lift','base_m':q,
                 'slope':lam,'shifts':list(S),'sheets':L,'delta':1,
                 'codes':list(range(len(S))),'delete_h1':[]}
            if L==2:
                rec.update(states=2*len(S)-overlap,delta=delta,delete_h1=drop,
                           method='linear_two_sheet_deletion')
            elif L==3:
                rec.update(codes=monotone_codes(q,S),method='linear_monotone_three_sheet')
            add(q*L,rec)
    return B

def linear_seed(q,t,limit):
    rec=linear_bank(limit)[q];assert rec['states']>=t
    if rec['method']=='retained_latin_odd_lift':M=eight_latin_lift(rec['sheets'])
    elif rec['method']=='single_sum_seed':
        m,lam,S=rec['base_m'],rec['slope'],rec['shifts']
        M=[[(a,(lam*a+s)%m,(-(lam+1)*a-s)%m)for a in range(m)]for s in S]
    else:
        M=linear_lift(rec['base_m'],rec['slope'],rec['shifts'],rec['sheets'],
                      rec['delta'],rec['codes'],rec['delete_h1'])
    return M[:t]

@lru_cache(maxsize=16)
def linear_plans(limit,t):
    B=linear_bank(limit)
    qs=sorted(q for q,rec in B.items()if rec['states']>=t)
    paths={q:(q,())for q in qs}
    for m in range(1,limit+1):
        if m not in paths:continue
        base,steps=paths[m]
        for q in qs:
            nxt=m+q-1
            if nxt>limit:break
            if nxt not in paths or len(paths[nxt][1])>len(steps)+1:
                paths[nxt]=(base,steps+(q,))
    return paths

def bounds(m,limit=1000):
    if m<1:raise ValueError('m must be positive')
    limit=max(m,limit)
    old=OLD[m-1]if m<=1000 else v10.bounds(m,m)
    ell=1 if m<=2 else 2
    for t in range(3,13):
        if m in linear_plans(limit,t):ell=max(ell,t)
    g=max(old['tripartite_g3_lower_bound'],ell)
    f=max(old['unrestricted_f3_at_3m_lower_bound'],g,10 if m>=2 else 1)
    return {'m':m,'vertices':3*m,
            'unrestricted_f3_at_3m_lower_bound':f,'unrestricted_edges':m*f,
            'tripartite_g3_lower_bound':g,'tripartite_edges':m*g,
            'linear_tripartite_g3_lower_bound':ell,'linear_tripartite_edges':m*ell,
            'previous_v10_unrestricted_lower_bound':old['unrestricted_f3_at_3m_lower_bound'],
            'previous_v10_tripartite_lower_bound':old['tripartite_g3_lower_bound']}

def construct_linear(m):
    limit=max(1000,m);t=bounds(m,limit)['linear_tripartite_g3_lower_bound']
    if t==1:
        M=[[(a,a,a)for a in range(m)]];H,parts=globalize(M)
        meta={'method':'disjoint_triples'}
    elif t==2:
        assert m>=3
        M=[[(a,a,a)for a in range(m)],[(a,(a+1)%m,(a+2)%m)for a in range(m)]]
        H,parts=globalize(M);meta={'method':'linear_two_color_cycle'}
    else:
        base,steps=linear_plans(limit,t)[m]
        H,parts=globalize(linear_seed(base,t,limit))
        for q in steps:H,parts=graft(H,linear_seed(q,t,limit),parts)
        meta={'method':'linear_universal_grafting'if steps else'linear_seed',
              'base_m':base,'base_seed':linear_bank(limit)[base],
              'seed_orders':list(steps),'graft_seeds':[linear_bank(limit)[q]for q in steps],
              'graft_steps':len(steps)}
    assert v10.max_codegree(H)==1
    return H,parts,{**meta,**check(H,parts),'maximum_codegree':1,'setting':'linear_balanced_tripartite'}

def construct(m,partite=False,linear=False,states=None):
    if linear:H,parts,meta=construct_linear(m)
    elif not partite and m==8:
        H,parts=twenty_four(),None
        meta={**check(H),'method':'finite_twenty_four_vertex_seed','setting':'unrestricted',
              'maximum_codegree':v10.max_codegree(H)}
    else:
        r=bounds(m);key='tripartite_g3_lower_bound'if partite else'unrestricted_f3_at_3m_lower_bound'
        H,parts,meta=v10.construct(m,partite)
        if len(H)<r[key]:H,parts,meta=construct_linear(m)
    if states is not None:
        if not 1<=states<=len(H):raise ValueError('Requested state count is not supplied by the proved bound')
        H=H[:states]
        meta={**meta,**check(H,parts),'maximum_codegree':v10.max_codegree(H),
              'whole_color_deletion':True}
    return H,parts,meta

def export(m,path,partite=False,linear=False,states=None):
    H,parts,meta=construct(m,partite,linear,states)
    data={'parameters':meta,'vertex_labels':'0,...,3m-1','matchings':H}
    if parts is not None:data['parts']=[sorted(P)for P in parts]
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(data,indent=2)+'\n')
    return meta

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--m',type=int);p.add_argument('--states',type=int)
    setting=p.add_mutually_exclusive_group()
    setting.add_argument('--tripartite',action='store_true');setting.add_argument('--linear',action='store_true')
    p.add_argument('--export',type=Path);p.add_argument('--table',action='store_true')
    a=p.parse_args()
    if a.table:
        R=[bounds(m)for m in range(1,1001)]
        for r in R:
            m=r['m'];ell=r['linear_tripartite_g3_lower_bound']
            assert r['unrestricted_f3_at_3m_lower_bound']>=(10 if m>=2 else 1)
            assert ell<=r['tripartite_g3_lower_bound']<=r['unrestricted_f3_at_3m_lower_bound']
            for t,start in [(8,33),(9,62),(10,74),(11,112),(12,113)]:
                if m>=start:assert ell>=t
        (ROOT/'small_m_1_to_1000_v11.json').write_text(json.dumps(R,indent=2)+'\n')
        print(json.dumps({'rows':len(R),'new_unrestricted_orders':[r['m']for r in R if r['unrestricted_f3_at_3m_lower_bound']>r['previous_v10_unrestricted_lower_bound']],
                          'new_balanced_orders':[r['m']for r in R if r['tripartite_g3_lower_bound']>r['previous_v10_tripartite_lower_bound']],
                          'new_linear_setting':True}))
    if a.m is not None:
        print(json.dumps(export(a.m,a.export,a.tripartite,a.linear,a.states)if a.export else
                         construct(a.m,a.tripartite,a.linear,a.states)[2]))
    if not a.table and a.m is None:p.print_help()

if __name__=='__main__':main()
