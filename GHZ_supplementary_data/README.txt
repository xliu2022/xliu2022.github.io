Triple-source networks for “Unbounded GHZ dimension from fixed-order multiphoton sources”
Heng Li and Xizhi Liu

A row of small_m_1_to_1000_cumulative.json specifies an exact matching size m. The network has n = 3m paths. The dimension columns give certified constructions for unrestricted triples, triples meeting three fixed groups once, and triples satisfying that grouping with no repeated pair of paths. The edge columns give the corresponding source counts. For a smaller target dimension, delete entire matching classes. The columns give lower bounds on the attainable dimensions.

To export a network, run these commands in the extracted supplementary_data directory. Python 3.10 or later suffices.

    python3 construct_cumulative.py --m 3 --states 12 --export my_9_paths.json
    python3 construct_cumulative.py --m 13 --tripartite --states 12 --export my_39_paths.json
    python3 construct_cumulative.py --m 39 --linear --states 12 --export my_117_paths.json

Every row of an exported matching list is one GHZ branch. Each triple in the row gives the path labels for one source with that internal label. Paths are numbered from 0 through 3m-1. The export includes the path groups whenever it uses a grouped construction. Source count is dm for dimension d.

The examples directory contains these networks:

    twelve_on_9.json           d = 12, 36 sources, first 12 classes of Table 1
    thirteen_on_9.json         d = 13, 39 sources, all classes of Table 1
    eight_on_15.json           d = 8, 40 sources, the modulo-5 construction
    thirteen_on_18.json        d = 13, 78 sources, unrestricted
    twelve_on_39.json          d = 12, 156 sources, three fixed groups
    twelve_linear_on_117.json  d = 12, 468 sources, three groups and no repeated pair

Gu, Chen, and Krenn reported in 2020 that dimension 13 is attainable on 9 paths. The 9-path example in this archive is the finite-field construction in Table 1, also used as the base of the 18-path example. https://doi.org/10.1103/PhysRevA.101.033816

The periodic family with 11 matchings on 9q vertices uses 11 rows of the 9-vertex list posted by discrete_things in 2024. https://math.stackexchange.com/posts/comments/10703741

The periodic Latin-square construction starts with entry 16 of Ian Wanless's file LS8alphan. Egan and Wanless give the order-8 enumeration. https://users.monash.edu.au/~iwanless/data/transversals/LS8alphan  https://doi.org/10.1002/jcd.20297

Run the paper-specific checks with

    python3 verify_paper_examples.py

This checks the source lists against Table 1 and the modulo-5 formula, enumerates the perfect matchings of the 9-, 15-, and 18-path examples, and enumerates the two affine sum classes on 39 paths. On 117 paths it reconstructs the three-sheet lift of the verified 13-element base and checks each path pair. Completeness of the lifted list follows from the cyclic signature-lift theorem in Section S4 and its linear refinement in Section S5 of the Supplemental Material. The script also regenerates every dimension and source count in the 1000-row table.

To check the estimates for the growth construction in Appendix A of the article, run

    python3 verification/verify_unified_bounds.py

It checks owner counts on small rectangles and the size and active-support inequalities for rectangular and planar parameters.

Two further checks for the finite networks are

    python3 verification/verify_export.py examples/thirteen_on_18.json
    python3 retained_v11/verification/verify_small_ten_seeds.py

The first accepts any small export and compares its entire family of perfect matchings with the supplied list. The second exhaustively checks the ten initial matching sizes used for the unrestricted dimension-10 construction. Larger instances are proved by grafting and lifting the finite inputs.

The retained_v11 directory holds the original construction routines, finite seed certificates, and their dependencies. The numerical table reproduces the cumulative PM3 construction record of 22 November 2026. For nine unrestricted paths the exporter uses the algebraic labeling of the paper. The optional command python3 verification/reproduce_finite.py additionally requires a C++17 compiler and runs two independent checks of the nine-path upper bound.
