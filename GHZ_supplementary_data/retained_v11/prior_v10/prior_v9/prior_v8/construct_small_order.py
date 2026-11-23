#!/usr/bin/env python3
"""Construct a certified rigid tripartite hypergraph at the requested part size.

python construct_small_order.py --m 23 --export examples/custom_m23.json
python construct_small_order.py --table --verify

The finite ten-state examples use a complete exact-cover certificate.
The other constructions follow the analytic and previously certified lemmas
in the accompanying report. No assertion of optimality is made.
"""
from pathlib import Path
import argparse, json, sys
import importlib.util

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / 'verification'))
from balanced_extension import construct_best, structural_check
from exact_cover_shifts import count_shifts
from unrestricted_extension import construct_ten,eleven_lift,ten_lift,check_structure

TEN = {23:(0,1,2,5,11), 29:(0,1,2,5,17), 31:(0,1,2,5,11),
       37:(0,1,2,5,21), 41:(0,1,2,5,16), 43:(0,1,2,5,38),
       47:(0,1,2,11,23)}

def unrestricted_bound(m,tripartite_bound):
    """Combine inherited unrestricted bounds with this tripartite example."""
    candidates=[(tripartite_bound,'tripartite_construction')]
    if m>=2:candidates.append((9,'prior_uniform_nine_state_certificate'))
    if m>=681:candidates.append((10,'new_unrestricted_binary_extension'))
    if m%3==0:candidates.append((11,'prior_eleven_state_cyclic_lift'))
    if m%4==2:candidates.append((10,'prior_ten_state_odd_period_lift'))
    if m in {2:10,3:13,4:12,5:10,6:11,7:9}:
        candidates.append(({2:10,3:13,4:12,5:10,6:11,7:9}[m],'prior_finite_example'))
    best=max(x[0] for x in candidates)
    return best,sorted({method for value,method in candidates if value==best})

def construct_unrestricted(m):
    """Realize the recorded unrestricted bound, including inherited families."""
    if m<1:raise ValueError('m must be positive')
    if m==1:return [[(0,1,2)]],{'method':'single_edge'}
    if m in (3,4,5):
        data=json.loads((ROOT/'prior_results/unrestricted_finite_examples.json').read_text())
        candidates=[r for r in data['constructions']if r['n']==3*m]
        best=max(candidates,key=lambda r:r['d'])
        M=[[tuple(v-1 for v in e)for e in row]for row in best['matchings']]
        return M,{'method':'prior_finite_example'}
    if m==7:
        data=json.loads((ROOT/'prior_results/unrestricted_nine/nine_on_21.json').read_text())
        return [[tuple(v-1 for v in data['edges'][j])for j in row]for row in data['matchings']],{'method':'prior_finite_example'}
    if m%3==0:return eleven_lift(m//3),{'method':'prior_eleven_lift','layers':m//3}
    if m%4==2:return ten_lift(m//2),{'method':'prior_ten_odd_lift','layers':m//2}
    if m in TEN:
        M,meta=construct(m)
        return [[(a,m+b,2*m+c)for a,b,c in row]for row in M],meta
    if m>=681:return construct_ten(m)
    path=ROOT/'prior_results/unrestricted_nine/common_height/verify.py'
    spec=importlib.util.spec_from_file_location('_prior_common_height',path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    n,edges,colors=module.build_hypergraph(module.word_for_n(3*m))
    def unpack(mask):
        vv=[]
        while mask:
            bit=mask&-mask;mask-=bit;vv.append(bit.bit_length()-1)
        return tuple(vv)
    return [[unpack(edges[j])for j in sorted(row)]for row in sorted(colors,key=lambda x:tuple(sorted(x)))],{'method':'prior_uniform_nine_state_certificate'}

def construct(m):
    if m in TEN:
        S=TEN[m]
        matchings=[[(a,(a+s)%m,(h-2*a-s)%m) for a in range(m)]
                   for h in (0,1) for s in S]
        method={'method':'modular_ten_finite_certificate','shift_set':list(S),
                'proof_status':'analytic reduction plus complete finite enumeration'}
    else:
        matchings,method=construct_best(m)
    return matchings,{**structural_check(matchings),**method}

def write_example(m,path):
    M,data=construct(m)
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps({'parameters':data,'matchings':M},indent=2)+'\n')
    return data

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--m',type=int)
    p.add_argument('--export',type=Path)
    p.add_argument('--table',action='store_true')
    p.add_argument('--verify',action='store_true')
    p.add_argument('--unrestricted',action='store_true',help='Export f3(3m) construction with global vertex labels')
    args=p.parse_args()
    if args.verify:
        records=[]
        for m,S in TEN.items():
            r=count_shifts(m,S,cap=6)
            assert r['complete'] and r['count']==5
            assert {tuple(sorted({s for a,s in sol})) for sol in r['solutions']}=={(s,)for s in S}
            records.append({k:v for k,v in r.items()if k!='solutions'})
        (ROOT/'verification'/'packaged_ten_verification.json').write_text(json.dumps(records,indent=2)+'\n')
        print(json.dumps({'verified_ten_state_orders':list(TEN)}))
    if args.table:
        records=[construct(m)[1]for m in range(1,1001)]
        for r in records:
            r['tripartite_g3_lower_bound']=r['states']
            val,methods=unrestricted_bound(r['m'],r['states'])
            r['unrestricted_f3_at_3m_lower_bound']=val
            r['unrestricted_methods']=methods
        assert all(x['states']>=7 for x in records if x['m']>=135)
        assert all(x['states']>=8 for x in records if x['m']>=263)
        (ROOT/'small_m_1_to_1000.json').write_text(json.dumps(records,indent=2)+'\n')
        print(json.dumps({'table_rows':len(records),'m_min':1,'m_max':1000}))
    if args.m is not None:
        if args.unrestricted:
            M,meta=construct_unrestricted(args.m)
            result={**check_structure(M),**meta,'setting':'unrestricted','vertex_labels':'0,...,3m-1'}
            if args.export:
                args.export.parent.mkdir(parents=True,exist_ok=True)
                args.export.write_text(json.dumps({'parameters':result,'matchings':M},indent=2)+'\n')
        elif args.export:result=write_example(args.m,args.export)
        else:result=construct(args.m)[1]
        print(json.dumps(result))
    if not(args.verify or args.table or args.m is not None):p.print_help()

if __name__=='__main__':main()
