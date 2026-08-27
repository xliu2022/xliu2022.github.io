# ZS5 (d=4) local-quadratic certificate: independent red-team audit

Date: 2026-08-25

## Outcome

**[PASS: finite (d=4) certificate, conditional only on the already stated
ZS5 reduction]**.  I independently checked the geometry, the omitted-axis
legality of the dual, the (S_3/S_4) symmetry glue, the exact repair of every
selected active Bernstein row, the outward algebraic intervals, the closed-top
corrections, and the full (5{,}164{,}740)-row replay.  I found no mathematical
or implementation defect in this certificate.

This is deliberately a component-level conclusion.  By itself it proves the
remaining (d=4) stratum of the zero-atom ZS5 statement; it does not by itself
re-audit the preceding reductions, the already certified (d\leq3) strata, or
the separate finite-support/threshold interface used in the final (r=5)
argument.

## 1. Dual legality and omitted-axis check

The certified pointwise inequality is

\[
 K_b(w_1,w_2,w_3,w_4)
 \leq 1+\sum_{i=1}^4 (w_i-m)g(w_{-i},b),\qquad m=1/5.
\]

On each canonical projected threshold cell, (g) is quadratic in the three
remaining centered random coordinates and centered (b).  I independently
reconstructed all (4\cdot3493=13{,}972) projection records.  In every record

\[
 i\notin\operatorname{axes}(g_i),\qquad
 \operatorname{axes}(g_i)=\{0,1,2,3\}\setminus\{i\}.
\]

Thus no coefficient or canonicalization step reintroduces the deleted random
axis.  Since (b) is a fixed parameter and the other (W_j)'s are independent
of (W_i),

\[
 \mathbb E[(W_i-m)g(W_{-i},b)]
 =\mathbb E(W_i-m)\,\mathbb E g(W_{-i},b)=0.
\]

The projected descriptor has multiple minimizing orientations in 8,045 of the
13,972 records.  These are not an ambiguity: every two such orientations
differ by a stabilizer of the same canonical projected key.  Reconstructing
all stabilizers independently gave the distribution

\[
 234\times1,\qquad268\times2,\qquad43\times6,
\]

and exactly the recorded quotient map

\[
 545\cdot15=8175\longrightarrow6759.
\]

The expanded repaired coefficients are exactly constant on every quotient
class.  Consequently a permutation of the four random axes only permutes the
four centered summands; checking one representative of every full (S_4)
cell orbit is legitimate.

## 2. Exhaustive exact geometry

I recompiled and reran

    work/emc_audit/scripts/r5_zs5_d4_local_arrangement_counts.cpp

from scratch, both normally and with Clang undefined-behavior sanitization.
Both runs returned

\[
 587\text{ universal vertices},\quad
 8410\text{ universal full cells},\quad
 48898\text{ cells in the 2500 global boxes},\quad
 3493\ S_4\text{-orbits}.
\]

The regenerated orbit flat was byte-identical to the frozen input, with
SHA-256

    07f674231616b548b340cb6c783f1736ae6070391931f68d290285094efd598f

I then reran the Fraction-only H-polytope reconstruction and pulling
triangulation verifier.  It rebuilt all 3,493 orbit representatives and all
40,990 five-dimensional pulling simplices and ended in `CERTIFIED`; the
recorded JSON and flattened geometry have hashes

    681ab4c16db61bdad98d121c6c265c18918169e9b0081f5d47f4591eebf361c2
    567faf32d329ebb1a3a653edcd4918136d1053a0e93ce4ab6f2faa86d8307fd4

respectively.  As a further cross-check, the exact kernel signature at the
centroid of every one of the 40,990 simplices agrees with the signature at the
centroid of its parent full cell.  Hence no simplex crosses an unrecorded
threshold wall.

## 3. The 12,794 active rows are an exact face, not a floating assumption

I reran the floating scanner on all 5,164,740 degree-four Bernstein rows and
retained the 20,000 smallest values.  Exactly 12,794 rows satisfy

\[
 |\text{gap}|<10^{-9};
\]

their largest absolute floating residual is
(1.39727\cdot10^{-14}).  The smallest returned row outside that set is already

\[
 4.371245767256782\cdot10^{-6}.
\]

More importantly, the proof does **not** infer from this numerical observation
that these rows must be symbolic zeros.  Starting only from their 12,794
simplex/multiindex metadata records and the exact geometry, my checker rebuilt
every row over 
(\mathbb Q[z]/(3125z^4+11250z^3+15250z^2+9225z-1024)).  The resulting Counter,
including all multiplicities, agrees exactly with the frozen system:

\[
 12794\text{ rows}\quad\longleftrightarrow\quad1620
 \text{ distinct exact equations}.
\]

Every one of those 1,620 equations has zero residual after inserting the
repaired coefficients.  Thus the selected numerical face is an actual
nonempty exact intersection.  Selecting a row that was not forced to be zero
by the original floating optimizer would still be harmless: the repaired
point is what matters.

## 4. FLINT repair and (\mathbb Q(z)) arithmetic

The four blocks choose exact invertible pivot minors of sizes

\[
 393,\quad381,\quad327,\quad393.
\]

The initial QR calculation is used only to choose a well-conditioned square
minor.  No proof claim depends on the floating rank threshold: FLINT solves
the chosen rational square systems, and all 1,620 equations are then replayed
with `Fraction`.  I reran the repair in a temporary file; it was byte-identical
to the frozen repair JSON, SHA-256

    4d537d614540337e4b55226fb7edfe91109fd38a00b6bc0e6b817a37b785b5f1

The maximum change from the stabilizer-averaged floating candidate at the
floating cusp approximation is (2.714\cdot10^{-12}).  This deviation is only
a diagnostic; exact feasibility follows from the quotient-ring replay, not
from its smallness.

The cusp polynomial is strictly increasing on the positive interval because

\[
 P'(z)=9225+30500z+33750z^2+12500z^3>0.
\]

Together with the opposite endpoint signs, this proves that the bisection
interval inside ((19/200,12/125)) contains exactly the intended root.  I
independently regenerated, with `Fraction`, all 6,759 coefficient enclosures
and the four enclosures of (1,z,z^2,z^3).  Every lower endpoint is an exact
floor and every upper endpoint an exact ceiling after multiplication by
(2^{90}); the regenerated flat is byte-identical to the certificate, SHA-256

    4dc19e34b4e13ce3a849ac286e8b7c8dc4c2ef4bcf0c60db4e047bc7be565239

Its maximum coefficient-interval width is
(162581340/2^{90}\approx1.31\cdot10^{-19}).  Treating occurrences of (z)
and the coefficients as independent intervals only enlarges the enclosure, so
it is conservative despite losing correlations.

## 5. Base polynomial, closed walls, and top corrections

Within a full threshold chamber the kernel is constant in the random axes,
and multiplication by (b>0) turns the dual gap into a polynomial of total
degree at most four.  The verifier constructs its symmetric blossom and checks
all

\[
 {4+5\choose5}=126
\]

degree-four Bernstein coefficients on every five-simplex.

On non-top closed walls, increasing every nonsaturated coordinate selects one
common northeast full chamber.  Its projections simultaneously select the
right-continuous value of every (g_i), so there is no incompatibility among
the four separately projected terms.  Every active threshold other than a
top singleton moves to the successful side.

The two exceptional singleton contributions after multiplying by (b) are
exactly

\[
 z^3\{b(1+z)-m\}\quad(w_i=1),\qquad mz^4\quad(b=1).
\]

I checked their signs and coefficient placement directly against the kernel:
the verifier subtracts the first correction for each random top face and the
second on the (b)-top face.  At intersections the loop adds them term by
term.  No mixed extra term is missing, since outcomes containing two positive
top variables were already successful in the inward chamber.

## 6. Full exact replay and active-row guard

I independently compiled and ran the final C++ interval verifier twice: once
at `-O3`, and once with

    -fsanitize=undefined -fno-sanitize-recover=undefined

The sanitizer run found no signed-`i128`, rational, indexing, or other checked
undefined behavior.  Both runs produced the identical final stream:

\[
\begin{array}{r|r}
\text{all rows}&5{,}164{,}740\\
\text{strictly positive rows}&5{,}151{,}946\\
\text{exact active rows}&12{,}794
\end{array}
\]

with

\[
 \min \text{strict lower endpoint}
 =\frac{5411340156129015928571}{2^{90}}
 \approx4.37124576668\cdot10^{-6}>0,
\]

and row-stream fingerprint

    49063aa32fa17a6f

An active row is not silently skipped.  Before the C++ pass, the wrapper
regenerates the 12,794-row exact system byte-for-byte and performs the 1,620
exact zero identities.  During the C++ pass, each listed active row must also
have an interval containing zero.  Every unlisted row must have a strictly
positive lower endpoint.  Thus the union of the two stages covers every row
exactly once.

The standalone C++ executable intentionally consumes rather than rederives
the algebraic coefficient intervals and active metadata.  It should therefore
not be cited alone.  The proof entry point is the wrapper

    .venv/bin/python work/emc_audit/scripts/r5_zs5_d4_quadratic_exact_manifest.py

which regenerates the active system, exact repair, and interval flat before
the complete scan.  At the time of this audit, the final manifest has SHA-256

    1e101b13ddec6f8b4496fe6f74ec3df0e5382fd12019a5569ea039fae6f1b208

and records the wrapper hash

    c42bb6d57f7c4361866ca8b8f1545b51782e353b80d6b05adcf963efbf91a5f0

## 7. Independent red-team artifacts

The independent checker and its machine-readable output are

    work/emc_audit/scripts/r5_zs5_d4_quadratic_redteam_exact.py
    work/emc_audit/computer/r5_zs5_d4_quadratic_redteam_exact.json

The checker separately verifies:

* exhaustive extraction of the 12,794 numerical near-zero rows from the
  complete 5,164,740-row floating scan;
* exact reconstruction and multiplicity comparison of all 1,620 equations;
* the complete stabilizer quotient and every omitted-axis projection;
* all exact active residuals;
* byte-level reconstruction of the flattened exact geometry;
* exact outward regeneration of every entry in the 90-bit interval flat.

No search output is used as a proof without a subsequent exact replay.
