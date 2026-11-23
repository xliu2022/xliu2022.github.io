# Construction certificates

The data accompany *Unbounded GHZ dimension from fixed-order multiphoton sources* by Heng Li and Xizhi Liu. A network of matching size $m$ has $3m$ paths. Each prescribed perfect matching gives one GHZ branch.

## Periodic constructions (Section S2)

The [periodic certificates](retained_v11/prior_v10/prior_v9/prior_v8/verification/prior_periodic_certificate.json) give the eight-transversal Latin square, the six-path network with ten complementary pairs, and the shifts of both constructions. Each record lists the local configurations, the weighted transitions between boundary states, their heights, the unit loops, and the exceptional two-cycle. The [independent enumeration](retained_v11/prior_v10/prior_v9/prior_v8/verification/prior_periodic_verifier.py) reconstructs these data.

The Latin square is entry 16 in Wanless’s `LS8alphan` catalogue. Its eight transversals partition the cells. The six-path construction uses the twenty triples of the network in Gu, Chen, and Krenn, Fig. 8. For even period $q$, walks on the two-cycle contribute an additional $2\cdot5^{q/2}$ matchings in the Latin construction and $2\cdot2^{q/2}$ in the six-path construction. For odd $q$ only the loops contribute. [Latin-square catalogue](https://users.monash.edu.au/~iwanless/data/transversals/LS8alphan) · [Egan and Wanless](https://doi.org/10.1002/jcd.20297) · [Gu, Chen, and Krenn](https://doi.org/10.1103/PhysRevA.101.033816)

The nine-vertex base retains eleven matchings from the list posted by `discrete_things`. Its certificate contains 333 configurations, 250 distinct transitions, and eleven loops of weight one. The height increases on all 239 other arcs. [Base and transfer certificate](retained_v11/prior_v10/prior_v9/prior_v8/prior_results/eleven_certificate.json) · [Enumeration](retained_v11/prior_v10/prior_v9/prior_v8/prior_results/eleven_verify.py) · [Original matching list](https://math.stackexchange.com/posts/comments/10703741)

## Affine seeds and cyclic lifts (Sections S4 and S5)

The single-sum matching indexed by $s$ consists of $(a,\lambda a+s,-(\lambda+1)a-s)$ for $a$ modulo $m$. The last column gives the three-sheet codes used by the linear threshold constructions. In each row the codes follow the displayed order of $S$. With marked coordinates $0,1,2$, the code for $s$ is $\alpha_s(0)+2\alpha_s(1)+4\alpha_s(2)$.

| Base size $m$ | Slope $\lambda$ | Shifts $S$ | Three-sheet codes used in S7 |
| ---: | ---: | --- | --- |
| 3 | 1 | 0, 1, 2 | |
| 5 | 1 | 0, 1, 2, 3 | |
| 7 | 1 | 0, 1, 2, 3 | 7, 6, 4, 0 |
| 11 | 2 | 0, 1, 2, 3, 5 | 7, 6, 4, 0, 5 |
| 13 | 2 | 0, 1, 2, 4, 5, 8 | 7, 6, 4, 5, 1, 3 |
| 17 | 2 | 0, 1, 2, 4, 5 | 7, 6, 4, 5, 1 |
| 19 | 3 | 0, 1, 3, 12, 14, 15 | 7, 6, 3, 2, 5, 4 |
| 23 | 2 | 0, 1, 2, 9, 10, 11 | 7, 6, 4, 5, 1, 0 |
| 29 | 8 | 0, 1, 2, 7, 20, 21 | 7, 6, 4, 3, 5, 1 |
| 31 | 4 | 0, 1, 3, 16, 21, 28 | 7, 6, 5, 3, 4, 2 |
| 37 | 1 | 0, 1, 2, 5, 21 | |
| 41 | 1 | 0, 1, 2, 5, 16 | |
| 43 | 1 | 0, 1, 2, 5, 38 | |
| 47 | 1 | 0, 1, 2, 11, 23 | |

The cases $m=3,5,7$ are proved in the supplement. For $m=11,13,17,19,23,29,31$, the Python and C++ records give the counts of perfect matchings and search nodes. For $m=37,41,43,47$ the slope is one. [Python enumeration records](retained_v11/prior_v10/prior_v9/verification/finite_slope_checks.json) · [C++ enumeration records](retained_v11/prior_v10/prior_v9/verification/finite_slope_cpp_checks.json) · [Larger slope-one seeds](retained_v11/prior_v10/prior_v9/prior_v8/verification/all_verified_ten_orders.json)

The base with $m=13$ gives the 39-path network. Its three-sheet lift gives the 117-path linear network. The lift records give the slopes, shifts, bit vectors, and two-sheet deletions. [39-path source list](examples/twelve_on_39.json) · [117-path source list](examples/twelve_linear_on_117.json) · [Lift parameters and pair checks](retained_v11/verification/linear_lift_checks.json) · [Construction routines](retained_v11/linear_cyclic_lift.py)

A linear seed with two matchings is available at every $m\ge3$. Take the triples $(a,a,a)$ and $(a,a+1,a+2)$ for $a$ modulo $m$. Triples from the two matchings share at most one vertex. In any perfect matching the sum of the second coordinates equals the sum of the first modulo $m$, so the number of triples of the second type is either zero or $m$.

## Unrestricted examples (Section S6)

The JSON exports use path labels from 0 through $n-1$. Every row of `matchings` lists the triples of one matching.

| Paths | Perfect matchings | Sources | Source list |
| ---: | ---: | ---: | --- |
| 12 | 12 | 48 | [12 on 12](examples/twelve_on_12.json) |
| 15 | 10 | 50 | [10 on 15](examples/ten_on_15.json) |
| 18 | 13 | 78 | [13 on 18](examples/thirteen_on_18.json) |
| 21 | 10 | 70 | [10 on 21](examples/ten_on_21.json) |
| 24 | 10 | 80 | [10 on 24](examples/ten_on_24.json) |

For the two-layer 18-path network, the vertex $(v,j)$ of the paper has JSON label $2v+j$. Each file in the table can be checked with `python3 verification/verify_export.py` followed by its path. The [ten initial matching lists](retained_v11/verification/ten_states_orders_two_to_eleven.json) and their [exact-cover check](retained_v11/verification/verify_small_ten_seeds.py) certify ten branches at each matching size from 2 through 11, including a balanced tripartite example at size 11.

## Residue representatives (Section S7)

The [threshold certificates](retained_v11/verification/linear_thresholds_refined.json) contain one record for each target dimension $t=8,9,10,11,12$. In each record `used_seeds` gives the base size, sheet count, dimension, and any deleted matchings. The entries in `apery` give every residue, its representative, and its decomposition into increments. Grafting a seed of part size $q$ adds $q-1$ to the matching size.

| Target $t$ | Graft increment $a$ | Largest representative $W$ | All part sizes from $W-a+2$ |
| ---: | ---: | ---: | ---: |
| 8 | 7 | 38 | 33 |
| 9 | 21 | 81 | 62 |
| 10 | 25 | 97 | 74 |
| 11 | 37 | 147 | 112 |
| 12 | 38 | 149 | 113 |

Run `python3 retained_v11/verification/verify_linear_arithmetic.py` from the extracted data directory to check the 128 representatives, their decompositions, and the seed parameters. The proof in Section S7 extends each representative by multiples of $a$. [Verification program](retained_v11/verification/verify_linear_arithmetic.py) · [Verification results](retained_v11/verification/linear_arithmetic_checks.json)
