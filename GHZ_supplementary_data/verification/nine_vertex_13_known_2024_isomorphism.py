#!/usr/bin/env python3
"""Exact incidence-graph isomorphism check against the supplied 2024 example.

No graph library is needed.  We backtrack through vertex bijections,
pruning by degrees, pair codegrees and fully mapped triples, and then
explicitly verify the induced bijection of the colored incidence graphs.
The published matching list is copied verbatim from the task.
"""
from collections import Counter
from itertools import combinations, permutations, product
from pathlib import Path
import json

ROOT=Path(__file__).resolve().parent
SOURCE_URL='https://math.stackexchange.com/questions/4994589/3-uniform-hypergraphs-with-only-disjoint-perfect-matchings'
PUBLISHED=[
    ['123','456','789'],['128','369','457'],['129','356','478'],
    ['134','279','568'],['137','246','589'],['145','238','679'],
    ['147','269','358'],['149','268','357'],['156','247','389'],
    ['159','248','367'],['167','235','489'],['168','257','349'],
    ['178','259','346'],
]
KNOWN_PM=[frozenset(tuple(int(c) for c in e) for e in M) for M in PUBLISHED]
KNOWN_EDGES=set().union(*KNOWN_PM)
OUR_DATA=json.loads((ROOT/'nine_vertex_13_state_matchings.json').read_text())
OUR_EDGES={tuple(e['vertices']) for e in OUR_DATA['hyperedges']}
OUR_PM={frozenset(tuple(e) for e in M['edges']):M['matching_label']
        for M in OUR_DATA['perfect_matchings']}


def partitions(vertices):
    if not vertices:
        yield frozenset()
        return
    first=vertices[0]
    for tail in combinations(vertices[1:],2):
        edge=(first,)+tail
        rest=tuple(v for v in vertices if v not in edge)
        for matching in partitions(rest):
            yield matching|{edge}


def main():
    assert len(KNOWN_EDGES)==len(OUR_EDGES)==39
    assert len(KNOWN_PM)==len(set(KNOWN_PM))==13
    for M in KNOWN_PM:
        assert sorted(v for e in M for v in e)==list(range(1,10))
    all_known_partitions=list(partitions(tuple(range(1,10))))
    assert len(all_known_partitions)==280
    valid=[M for M in all_known_partitions if M<=KNOWN_EDGES]
    assert len(valid)==13 and set(valid)==set(KNOWN_PM)

    known_degree=Counter(v for e in KNOWN_EDGES for v in e)
    our_degree=Counter(v for e in OUR_EDGES for v in e)
    known_pair=Counter(p for e in KNOWN_EDGES for p in combinations(e,2))
    our_pair=Counter(p for e in OUR_EDGES for p in combinations(e,2))
    mapping={}
    used=set()
    search_nodes=0

    def search():
        nonlocal search_nodes
        search_nodes+=1
        if len(mapping)==9:
            return True
        candidates_by_vertex=[]
        for x in range(1,10):
            if x in mapping:
                continue
            candidates=[]
            for y in range(9):
                if y in used or known_degree[x]!=our_degree[y]:
                    continue
                if any(known_pair[tuple(sorted((x,u)))]!=our_pair[tuple(sorted((y,v)))]
                       for u,v in mapping.items()):
                    continue
                compatible=True
                for u,v in combinations(mapping,2):
                    known_edge=tuple(sorted((x,u,v)))
                    our_edge=tuple(sorted((y,mapping[u],mapping[v])))
                    if (known_edge in KNOWN_EDGES)!=(our_edge in OUR_EDGES):
                        compatible=False
                        break
                if compatible:
                    candidates.append(y)
            candidates_by_vertex.append((len(candidates),x,candidates))
        count,x,candidates=min(candidates_by_vertex)
        if not count:
            return False
        for y in candidates:
            mapping[x]=y
            used.add(y)
            if search():
                return True
            used.remove(y)
            del mapping[x]
        return False

    isomorphic=search()

    def codegree_four_certificate(edge_set):
        pairs=Counter(p for e in edge_set for p in combinations(e,2))
        vertices={v for e in edge_set for v in e}
        adjacency={v:set() for v in vertices}
        for (u,v),d in pairs.items():
            if d==4:
                adjacency[u].add(v)
                adjacency[v].add(u)
        assert all(len(n)==2 for n in adjacency.values())
        unseen=set(vertices)
        components=[]
        while unseen:
            todo=[min(unseen)]
            found=set()
            while todo:
                v=todo.pop()
                if v not in found:
                    found.add(v)
                    todo.extend(adjacency[v]-found)
            unseen-=found
            components.append(tuple(sorted(found)))
        assert len(components)==3 and all(len(c)==3 for c in components)
        return {'codegree_four_graph_components':components,
                'components_that_are_hyperedges':[c for c in components if c in edge_set],
                'components_not_hyperedges':[c for c in components if c not in edge_set]}

    known_certificate=codegree_four_certificate(KNOWN_EDGES)
    our_certificate=codegree_four_certificate(OUR_EDGES)
    assert len(known_certificate['components_that_are_hyperedges'])==2
    assert len(our_certificate['components_that_are_hyperedges'])==3
    assert not isomorphic

    # Independent exhaustive check: every hypergraph isomorphism must preserve
    # the codegree-four graph, which is 3K3 on both sides.  There are precisely
    # 3!*(3!)^3 = 1296 vertex bijections preserving those components.
    candidate_maps=0
    successful_maps=[]
    known_components=known_certificate['codegree_four_graph_components']
    our_components=our_certificate['codegree_four_graph_components']
    for target_components in permutations(our_components):
        for within in product(*(permutations(c) for c in target_components)):
            candidate={u:v for source,target in zip(known_components,within)
                       for u,v in zip(source,target)}
            candidate_maps+=1
            image={tuple(sorted(candidate[v] for v in e)) for e in KNOWN_EDGES}
            if image==OUR_EDGES:
                successful_maps.append(candidate)
    assert candidate_maps==1296 and successful_maps==[]

    def cycle_histogram(edge_set,matchings):
        vertices={v for e in edge_set for v in e}
        counts=[]
        for M in matchings:
            arcs={(x,y) for B in M for x in B for y in vertices-set(B)
                  if tuple(sorted(set(B)-{x}|{y})) in edge_set}
            assert not any((y,x) in arcs for x,y in arcs)
            cycles=0
            for x,y,z in product(*M):
                cycles+=all(e in arcs for e in ((x,y),(y,z),(z,x)))
                cycles+=all(e in arcs for e in ((x,z),(z,y),(y,x)))
            counts.append(cycles)
        return dict(sorted(Counter(counts).items()))

    known_cycles=cycle_histogram(KNOWN_EDGES,KNOWN_PM)
    our_cycles=cycle_histogram(OUR_EDGES,OUR_PM)
    assert known_cycles!=our_cycles
    result={
        'status':'PASS','isomorphic':False,
        'primary_source_url':SOURCE_URL,
        'primary_source_author':'discrete_things',
        'primary_source_comment_timestamp_as_displayed':'2024-11-05 18:37',
        'source_list_provenance':'The 13 matchings appear in the comment by discrete_things on November 5, 2024.',
        'dimension_13_source_doi':'10.1103/PhysRevA.101.033816',
        'method':'Exact vertex backtracking, independent enumeration of all 1296 codegree-preserving maps, and a hand-verifiable nonisomorphism invariant. Equivalent to color-preserving incidence-graph isomorphism.',
        'incidence_graph_vertices':48,'incidence_graph_edges':117,
        'backtracking_nodes':search_nodes,
        'all_codegree_preserving_vertex_bijections_examined':candidate_maps,
        'isomorphisms_found':0,
        'published_to_our_vertex_map':None,
        'published_example_certificate':known_certificate,
        'vector_sum_example_certificate':our_certificate,
        'published_c_M_histogram':known_cycles,
        'vector_sum_c_M_histogram':our_cycles,
        'independent_published_partition_check':{'partitions_examined':280,'perfect_matchings':13},
        'conclusion':'Gu, Chen, and Krenn reported in 2020 that dimension 13 is attainable on nine paths. The finite-field hypergraph and the hypergraph listed by discrete_things in 2024 are nonisomorphic.'
    }
    path=ROOT/'nine_vertex_13_known_2024_isomorphism.json'
    path.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
    print(path)


if __name__=='__main__':
    main()
