"""Check the finite examples and the 1000-row table used in the paper."""
from collections import Counter
from itertools import combinations
from pathlib import Path
import json

from construct_cumulative import bounds
from verification.verify_export import verify

ROOT = Path(__file__).resolve().parent


def read_example(name):
    return json.loads((ROOT / 'examples' / (name + '.json')).read_text())


def colors(rows):
    return {frozenset(tuple(sorted(edge)) for edge in row) for row in rows}


def check_affine_base(data):
    # The sum of all labels in a matching is 0 modulo 13. All 13
    # edge sums therefore agree whenever each edge sum is either 0 or 1.
    classes = [[], []]
    for row in data['matchings']:
        for a, b, c in row:
            b -= 13
            c -= 26
            h = (a + b + c) % 13
            assert h in (0, 1)
            classes[h].append((a, b, c))
    for h, edges in enumerate(classes):
        by_a = [[edge for edge in edges if edge[0] == a] for a in range(13)]
        covers = []

        def visit(a, used_b, used_c, chosen):
            if a == 13:
                covers.append(chosen)
                return
            for edge in by_a[a]:
                _, b, c = edge
                if not (used_b & (1 << b) or used_c & (1 << c)):
                    visit(a + 1, used_b | (1 << b), used_c | (1 << c),
                          chosen + [(a, b + 13, c + 26)])

        visit(0, 0, 0, [])
        planted = [row for row in data['matchings'] if sum(row[0]) % 13 == h]
        assert len(covers) == 6 and colors(covers) == colors(planted)


def check_three_sheet_lift(data):
    # Parameters of the three-sheet lift in the section on linear refinements.
    code = dict(zip((0, 1, 2, 4, 5, 8), (7, 6, 4, 5, 1, 3)))
    markers = {0: 0, 1: 1, 2: 2}
    expected = []
    for h in (0, 1):
        for s, word in code.items():
            if (s + 1) % 13 in code:
                successor = code[(s + 1) % 13]
                assert successor & word == successor
            row = []
            for a in range(13):
                bit = (word >> markers[a]) & 1 if a in markers else 0
                b = (2 * a + s) % 13
                c = (h - 3 * a - s) % 13
                for x in range(3):
                    row.append((a + 13*x, 39 + b + 13*((x + bit + h) % 3),
                                78 + c + 13*((x + bit + 2*h) % 3)))
            expected.append(row)
    assert colors(data['matchings']) == colors(expected)
    edges = {tuple(e) for row in data['matchings'] for e in row}
    pairs = Counter(pair for edge in edges for pair in combinations(sorted(edge), 2))
    assert len(edges) == 468 and max(pairs.values()) == 1
    for row in data['matchings']:
        assert sorted(v for e in row for v in e) == list(range(117))


def main():
    cases = [
        ('twelve_on_9', 36, 12),
        ('thirteen_on_9', 39, 13),
        ('twelve_on_12', 48, 12),
        ('eight_on_15', 40, 8),
        ('ten_on_15', 50, 10),
        ('thirteen_on_18', 78, 13),
        ('ten_on_21', 70, 10),
        ('ten_on_24', 80, 10),
    ]
    for name, sources, dimension in cases:
        record = verify(ROOT / 'examples' / (name + '.json'))
        assert (record['edges'], record['actual_matchings']) == (sources, dimension)
        print(f"{name}: {sources} sources, {dimension} perfect matchings (exhaustive)")
    paper = json.loads((ROOT / 'verification/nine_vertex_13_state_matchings.json').read_text())
    classes = [row['edges'] for row in paper['perfect_matchings']]
    for name, dimension in [('twelve_on_9', 12), ('thirteen_on_9', 13)]:
        assert read_example(name)['matchings'] == classes[:dimension]
    five = [[(a, 5+(a+s)%5, 10+(h-2*a-s)%5) for a in range(5)]
            for h in (0, 1) for s in range(4)]
    assert colors(five) == colors(read_example('eight_on_15')['matchings'])
    check_affine_base(read_example('twelve_on_39'))
    print('twelve_on_39: 156 sources, 6 + 6 perfect matchings (exhaustive by edge sum)')
    check_three_sheet_lift(read_example('twelve_linear_on_117'))
    print('twelve_linear_on_117: 468 sources, 12 lifted classes, maximum pair codegree 1')
    rows = json.loads((ROOT / 'small_m_1_to_1000_cumulative.json').read_text())
    assert rows == [bounds(m) for m in range(1, 1001)]
    for key, first, threshold in [('unrestricted_f3_at_3m_lower_bound', 3, 21),
                                   ('tripartite_g3_lower_bound', 13, 43),
                                   ('linear_tripartite_g3_lower_bound', 39, 113)]:
        assert min(row['m'] for row in rows if row[key] >= 12) == first
        assert all(row[key] >= 12 for row in rows if row['m'] >= threshold)
    print('All 1000 table rows reproduced.')


if __name__ == '__main__':
    main()
