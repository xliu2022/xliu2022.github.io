#!/usr/bin/env python3
"""Unrestricted binary-port extensions. Requires balanced_extension.py nearby.

The extension itself has a direct analytic proof. Prior 11- and10-state
cyclic templates retain their finite transfer certificates. This script
checks the new implementation, including old edges meeting several selected
vertices. It does not claim exhaustive enumeration of large hypergraphs.
"""
from collections import Counter
from itertools import combinations
from pathlib import Path
import argparse,json
from balanced_extension import binary_gadget,check_gadget

ELEVEN_ROWS=(
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
TEN_TEMPLATE=(
 ('123','13'),('124','14'),('125','5'),('126',''),
 ('134','34'),('135','5'),('136','13'),('145',''),
 ('146','46'),('156','5'),('234','24'),('235',''),
 ('236','36'),('245',''),('246','46'),('256','5'),
 ('345','3'),('346',''),('356','56'),('456','6'))

def eleven_lift(q):
    assert q>=1
    return [[tuple(sorted(9*((layer+(v in outgoing))%q)+int(v)-1 for v in edge))
             for layer in range(q)for edge,outgoing in row]
            for row in ELEVEN_ROWS]

def ten_lift(q):
    assert q>=1 and q%2==1
    d=dict(TEN_TEMPLATE);rows=[]
    for edge in sorted(e for e in d if '1'in e):
        comp=''.join(c for c in '123456'if c not in edge)
        row=[(edge,d[edge]),(comp,d[comp])]
        rows.append([tuple(sorted(6*((layer+(v in outgoing))%q)+int(v)-1 for v in e))
                     for layer in range(q)for e,outgoing in row])
    return rows

def check_structure(M):
    t=len(M);n=3*len(M[0]);E=[]
    for row in M:
        assert len(row)*3==n
        assert all(len(e)==3 and len(set(e))==3 for e in row)
        assert Counter(v for e in row for v in e)==Counter(range(n))
        E.extend(tuple(sorted(e))for e in row)
    assert len(set(E))==len(E)==t*n//3
    assert Counter(v for e in E for v in e)==Counter({v:t for v in range(n)})
    return {'n':n,'m':n//3,'states':t,'edges':len(E)}

def extend_unrestricted(M,seed=None,row_index=0):
    t=len(M);n=3*len(M[0]);X,P,B=binary_gadget(t);L=len(X)
    assert n>=L
    q=0 if seed is None else len(seed[0]);extra=0 if seed is None else 3*(q-1)
    if seed is not None:assert len(seed)==t and q>=1
    out=[]
    for i,matching in enumerate(M):
        N=[]
        for e in matching:
            N.append(tuple(sorted(v if v>=L or B[i][v]==1 else n+v for v in e)))
        for x in range(L):
            port=x if B[i][x]==0 else n+x
            cv=n+2*L+P[i][x]
            if seed is None or x!=row_index:
                N.append(tuple(sorted((port,n+L+x,cv))));continue
            # Reuse the removed B-row label for the first private vertex.
            private=[n+L+x]+list(range(n+3*L,n+3*L+extra))
            assert len(private)==3*q-2
            for e in seed[i]:
                N.append(tuple(sorted(port if v==0 else cv if v==1 else private[v-2]for v in e)))
        out.append(N)
    return out

def base_ten(m):
    if m>=3 and m%3==0:return eleven_lift(m//3)[:10],{'method':'prior_eleven_lift','layers':m//3}
    if m>=2 and m%4==2:return ten_lift(m//2),{'method':'prior_ten_odd_lift','layers':m//2}
    raise ValueError('No prior10-state family at this order')

def construct_ten(m):
    try:return base_ten(m)
    except ValueError:pass
    if m%3==2 and m>=683:
        base=m-512;seed=None;delta=512
    elif m%3==1 and m>=685:
        base=m-514;seed=eleven_lift(1)[:10];delta=514
    else:raise ValueError('Not covered by this theorem')
    assert base>=171 and base%3==0
    M=extend_unrestricted(eleven_lift(base//3)[:10],seed)
    assert len(M[0])==m
    return M,{'method':'unrestricted_binary','base_m':base,'increment':delta,
              'seed_q':None if seed is None else 3}

def exact_count(M):
    n=3*len(M[0]);edges=sorted({tuple(sorted(e))for row in M for e in row})
    masks=[sum(1<<v for v in e)for e in edges];inc=[[]for _ in range(n)]
    for k,e in enumerate(edges):
        for v in e:inc[v].append(k)
    count=nodes=0
    def rec(rem):
        nonlocal count,nodes
        nodes+=1
        if not rem:count+=1;return
        b=rem&-rem;v=b.bit_length()-1
        for j in inc[v]:
            if masks[j]&rem==masks[j]:rec(rem^masks[j])
    rec((1<<n)-1)
    return count,nodes

def main():
    p=argparse.ArgumentParser();p.add_argument('--m',type=int);p.add_argument('--export',type=Path)
    a=p.parse_args()
    if a.m is not None:
        M,method=construct_ten(a.m);record={**check_structure(M),**method}
        if a.export:a.export.write_text(json.dumps({'parameters':record,'matchings':M},indent=2)+'\n')
        print(json.dumps(record));return
    reports={'binary_gadget':check_gadget(10),'small_exact':[],'large_structure':[]}
    small=ten_lift(1)[:3]
    # Four selected vertices out of six ensure some old edges meet multiple selected vertices.
    for seed in (None,small,eleven_lift(1)[:3]):
        M=extend_unrestricted(small,seed);r=check_structure(M);count,nodes=exact_count(M)
        assert count==3
        reports['small_exact'].append({**r,'perfect_matchings':count,'search_nodes':nodes})
    for m in (681,682,683,684,685,688,692,721,722,727,999,1000):
        M,method=construct_ten(m);reports['large_structure'].append({**check_structure(M),**method})
    records=[]
    for m in range(1,1001):
        try:M,method=construct_ten(m)
        except ValueError:continue
        records.append({**check_structure(M),**method})
    for m in range(681,1001):assert any(r['m']==m for r in records)
    reports['coverage']={'from':681,'through':1000,'all_covered':True,
                         'latest_uncovered_order_by_this_family':680}
    here=Path(__file__).resolve().parent
    (here/'unrestricted_extension_checks.json').write_text(json.dumps(reports,indent=2)+'\n')
    (here/'unrestricted_ten_orders_to_1000.json').write_text(json.dumps(records,indent=2)+'\n')
    print(json.dumps(reports,indent=2))
if __name__=='__main__':main()
