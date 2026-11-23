#!/usr/bin/env python3
"""Export triple-source networks at the matching size and dimension requested.

Default: unrestricted 3-graph on 3*m vertices.
--tripartite: balanced tripartite construction.
--linear: balanced tripartite construction with maximum pair codegree one.
The table entries are lower bounds, not optimality claims.
"""
from pathlib import Path
import argparse
import importlib.util
import json
import sys

ROOT = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location(
    'cumulative_retained_v11', ROOT / 'retained_v11' / 'construct_practical.py')
v11 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(v11)


def bounds(m):
    row = v11.bounds(m)
    row.pop('previous_v10_unrestricted_lower_bound', None)
    row.pop('previous_v10_tripartite_lower_bound', None)
    if m == 6:
        row.update(unrestricted_f3_at_3m_lower_bound=13, unrestricted_edges=78)
        row['restored_historical_certificate'] = (
            'verification/13_on_18.json; first recorded in v2')
    return row


def construct(m, partite=False, linear=False, states=None):
    if m < 1:
        raise ValueError('The matching size m must be positive')
    if m == 3 and not partite and not linear:
        name = 'verification/nine_vertex_13_state_matchings.json'
        data = json.loads((ROOT / name).read_text())
        classes = [row['edges'] for row in data['perfect_matchings']]
        if states is not None and not 1 <= states <= len(classes):
            raise ValueError('The nine-path network supplies one to thirteen states')
        matchings = [[tuple(e) for e in row] for row in classes[:states]]
        metadata = {
            **v11.check(matchings),
            'method': 'affine_plane_order_three',
            'setting': 'unrestricted',
            'maximum_codegree': v11.v10.max_codegree(matchings),
            'finite_certificate': name,
        }
        return matchings, None, metadata
    if m == 6 and not partite and not linear:
        data = json.loads((ROOT / 'verification' / '13_on_18.json').read_text())
        edges = [tuple(v - 1 for v in e) for e in data['edges']]
        matchings = [[edges[e] for e in M] for M in data['matchings']]
        if states is not None:
            if not 1 <= states <= 13:
                raise ValueError('This certificate supplies from one to thirteen states')
            matchings = matchings[:states]
        metadata = {
            **v11.check(matchings),
            'method': 'restored_v2_thirteen_on_eighteen',
            'setting': 'unrestricted',
            'maximum_codegree': v11.v10.max_codegree(matchings),
            'finite_certificate': 'verification/13_on_18.json',
        }
        return matchings, None, metadata
    return v11.construct(m, partite=partite, linear=linear, states=states)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--m', type=int)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--tripartite', action='store_true')
    group.add_argument('--linear', action='store_true')
    parser.add_argument('--states', type=int)
    parser.add_argument('--export', type=Path)
    parser.add_argument('--table', action='store_true')
    args = parser.parse_args()
    if args.table:
        rows = [bounds(m) for m in range(1, 1001)]
        path = ROOT / 'small_m_1_to_1000_cumulative.json'
        path.write_text(json.dumps(rows, indent=2) + '\n')
        print(json.dumps({'table': str(path), 'rows': len(rows), 'restored_order': 6}))
    if args.m is not None:
        M, parts, meta = construct(args.m, args.tripartite, args.linear, args.states)
        if args.export:
            data = {'parameters': meta, 'vertex_labels': '0,...,3m-1', 'matchings': M}
            if parts is not None:
                data['parts'] = [sorted(P) for P in parts]
            args.export.parent.mkdir(parents=True, exist_ok=True)
            args.export.write_text(json.dumps(data, indent=2) + '\n')
        print(json.dumps(meta, indent=2))
    if not args.table and args.m is None:
        parser.print_help()


if __name__ == '__main__':
    main()
