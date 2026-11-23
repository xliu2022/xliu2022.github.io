#!/usr/bin/env python3
"""Independent exhaustive reconstruction of the 13-state nine-vertex example.

Vertices are labelled 3*x+y for (x,y) in F_3^2.  An edge is an unordered
triple whose vector sum belongs to {(0,0),(1,0),(0,1),(1,1)}.  All 280
unlabelled partitions of the vertex set into three triples are generated
directly; this file imports no previous search or verification code.

Outputs a complete edge/matching JSON, one CSV row per matching edge,
and the complete 280-row partition audit.  The data are combinatorial;
they make no claim about a physical implementation or success probability.
"""
from collections import Counter
from itertools import combinations, product
from pathlib import Path
from fractions import Fraction
import csv
import json
import math

ROOT=Path(__file__).resolve().parent
VERTICES=tuple(range(9))
SUM_ORDER=((0,0),(1,0),(0,1),(1,1))


def coordinate(v):
    return divmod(v,3)


def edge_sum(edge):
    points=[coordinate(v) for v in edge]
    return tuple(sum(p[i] for p in points)%3 for i in range(2))


def partitions(remaining):
    """Canonical recursion: the first block contains the least vertex."""
    if not remaining:
        yield ()
        return
    first=remaining[0]
    for other in combinations(remaining[1:],2):
        block=(first,)+other
        tail=tuple(v for v in remaining if v not in block)
        for rest in partitions(tail):
            yield (block,)+rest


def exchange_statistics(matching,edge_set):
    arcs=set()
    for block in matching:
        for x in block:
            for y in VERTICES:
                if y not in block:
                    replacement=tuple(sorted(set(block)-{x}|{y}))
                    if replacement in edge_set:
                        arcs.add((x,y))
    assert not any((y,x) in arcs for x,y in arcs)
    triangles=[]
    for x,y,z in product(*matching):
        forward=((x,y),(y,z),(z,x))
        backward=((x,z),(z,y),(y,x))
        if all(e in arcs for e in forward):
            triangles.append(forward)
        if all(e in arcs for e in backward):
            triangles.append(backward)
    usage=Counter(e for triangle in triangles for e in triangle)
    assert max(usage.values(),default=0)<=1
    return {'exchange_arcs':len(arcs),'directed_triangles':len(triangles)}


def main():
    edges=tuple(edge for edge in combinations(VERTICES,3) if edge_sum(edge) in SUM_ORDER)
    edge_set=set(edges)
    edge_id={e:f'E{i+1:02d}' for i,e in enumerate(edges)}
    all_partitions=list(partitions(VERTICES))
    expected=math.factorial(9)//(math.factorial(3)**3*math.factorial(3))
    assert len(all_partitions)==len(set(all_partitions))==expected==280
    valid=[p for p in all_partitions if all(e in edge_set for e in p)]
    for p in valid:
        sums={edge_sum(e) for e in p}
        assert len(sums)==1
    valid.sort(key=lambda p:(SUM_ORDER.index(edge_sum(p[0])),p))
    assert len(edges)==39 and len(valid)==13
    appearances=Counter(e for p in valid for e in p)
    assert set(appearances)==edge_set and set(appearances.values())=={1}
    assert Counter(edge_sum(p[0]) for p in valid)=={(0,0):4,(1,0):3,(0,1):3,(1,1):3}
    degrees=Counter(v for e in edges for v in e)
    assert [degrees[v] for v in VERTICES]==[13]*9

    owner={e:i for i,p in enumerate(valid) for e in p}
    matching_records=[]
    for i,p in enumerate(valid):
        matching_records.append({
            'state_index':i,'matching_label':f'M_{i:02d}',
            'vector_sum_class':edge_sum(p[0]),
            'edge_ids':[edge_id[e] for e in p],
            'edges':p,
            'edge_amplitude_product':' * '.join('a_'+edge_id[e] for e in p),
            **exchange_statistics(p,edge_set)
        })

    # Independent exact rational check of the local weight table and its line.
    weights=[Fraction(0)]
    for d in range(1,8):
        local=[]
        for k1 in range(4):
            for k2 in range(4):
                if k1+k2==d-1:
                    expectation=1-Fraction(k1+k2,6)+Fraction(k1*k2,27)
                    local.append(Fraction(d,3)*expectation)
        weights.append(max(local))
    assert [81*w for w in weights]==[0,27,45,57,62,65,63,63]
    assert all(81*w<=50+3*d for d,w in enumerate(weights))
    assert [d for d,w in enumerate(weights) if 81*w==50+3*d]==[4,5]
    pair_degrees=Counter(pair for e in edges for pair in combinations(e,2))
    assert sum(pair_degrees.values())==9*13

    record={
        'schema_version':1,
        'description':'Simple 3-uniform 13-rigid hypergraph on nine vertices; no partite restriction.',
        'vertex_numbering':'vertex_id = 3*x + y, with x,y in {0,1,2}; arithmetic modulo 3',
        'vertices':[{'vertex_id':v,'x':coordinate(v)[0],'y':coordinate(v)[1]} for v in VERTICES],
        'allowed_vector_sums':SUM_ORDER,
        'hyperedges':[{'edge_id':edge_id[e],'vertices':e,'vector_sum':edge_sum(e),'owner_state_index':owner[e]} for e in edges],
        'perfect_matchings':matching_records,
        'verification':{
            'all_unlabelled_three_triple_partitions_examined':280,
            'partitions_whose_three_edges_all_belong_to_H':13,
            'hyperedges':39,'pairwise_edge_disjoint':True,
            'each_edge_occurs_in_exactly_one_perfect_matching':True,
            'vertex_degrees':[degrees[v] for v in VERTICES],
            'pair_codegree_histogram':dict(sorted(Counter(pair_degrees.values()).items())),
            'c_M_histogram':dict(sorted(Counter(p['directed_triangles'] for p in matching_records).items())),
            'weighted_bound_table_81w':[int(81*w) for w in weights],
            'affine_weight_bound_equality_codegrees':[4,5],
        },
        'algebraic_identity':{
            'variables':'x_(v,j) for vertex v and state index j; a_E for each hyperedge E',
            'F':'sum over states j and edges E in M_j of a_E * product_{v in E} x_(v,j)',
            'operation':'Give every x_(v,j) multidegree e_v in Z^9; retain the component of multidegree (1,1,1,1,1,1,1,1,1) in F^3/3!.',
            'result':'sum_{j=0}^{12} (product_{E in M_j} a_E) * product_{v=0}^{8} x_(v,j)',
            'coefficient_of_each_listed_matching_term':1,
            'physical_claim':'None: this is a finite polynomial identity, independent of physical realizability.'
        }
    }
    json_path=ROOT/'nine_vertex_13_state_matchings.json'
    csv_path=ROOT/'nine_vertex_13_state_edges.csv'
    audit_path=ROOT/'nine_vertex_13_state_partition_audit.csv'
    json_path.write_text(json.dumps(record,indent=2)+'\n')
    with csv_path.open('w',newline='') as f:
        writer=csv.writer(f)
        writer.writerow(['state_index','matching_label','sum_x','sum_y','edge_position','edge_id','v1','v2','v3'])
        for p in matching_records:
            for j,(eid,edge) in enumerate(zip(p['edge_ids'],p['edges']),1):
                writer.writerow([p['state_index'],p['matching_label'],*p['vector_sum_class'],j,eid,*edge])
    valid_index={p:i for i,p in enumerate(valid)}
    with audit_path.open('w',newline='') as f:
        writer=csv.writer(f)
        writer.writerow(['partition_index','triple_1','triple_2','triple_3','all_edges_in_H','state_index'])
        for i,p in enumerate(all_partitions):
            writer.writerow([i,*(' '.join(map(str,e)) for e in p),int(p in valid_index),valid_index.get(p,'')])
    print(json.dumps(record['verification'],indent=2))
    print('Complete perfect matchings:')
    for p in matching_records:
        print(p['matching_label'],p['vector_sum_class'],p['edges'])
    print('Saved:',json_path,csv_path,audit_path,sep='\n')


if __name__=='__main__':
    main()
