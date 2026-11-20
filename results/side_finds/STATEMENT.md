# The partition tower: the exact statements, and what a proof needs (12 Nov 2026)

**In plain words.** For each prime ℓ, the partition numbers p(ℓᵇ·n + δ_b) form a tower: one level for each b. Going up two
levels multiplies the deep levels by a fixed ℓ-adic number, or (from 37 on) by a few of them at once. The statements below
say exactly which numbers. They come from a handful of classical modular forms, followed in ℓ-adic families to weight −2.
The same numbers govern the steps along every other prime Q, and they decide which primes Q give a law and how often.
A3 and E1 were first marked "proven"; they are short arguments from published lemmas that no referee has read, so
they are marked **derived** (changed 20 Nov 2026).
Each statement is marked **proven** (in print, QUOTED in `reports/Partition tower literature check.md`), **derived** (by an
argument here), **counted** (by our computations, with the file named) or **conjectured**.

## Setting

- ℓ ≥ 5 prime; for b ≥ 1, δ_b is the integer in [0, ℓᵇ) with 24·δ_b ≡ 1 (mod ℓᵇ). The tower: H_b(n) = p(ℓᵇ·n + δ_b).
  For 5, 7, 11, H_b is first divided by the power of ℓ that Watson's and Atkin's laws guarantee.
- **The class** S_ℓ: the newforms f of weight ℓ − 3 on Γ0(6) with U₂f = (2|ℓ)·2^((ℓ−5)/2)·f and U₃f = (3|ℓ)·3^((ℓ−5)/2)·f,
  that is Atkin–Lehner signs w₂ = −(2|ℓ), w₃ = −(3|ℓ) (Yang's G_m; proven for 13 ≤ ℓ ≤ 31, Guo–Ono, Yang).
- R_ℓ = ⌊(ℓ + 12)/24⌋, the rank bound of Folsom–Kent–Ono and Boylan–Webb (proven).
- For f in S_ℓ with a_ℓ(f) a unit on lane ℓ (ordinary), F_f is f's Hida family. Its members sit at the weights
  k = ℓ − 3 + (ℓ − 1)x, with U₂, U₃ values whose signs change by (2|ℓ), (3|ℓ) per step. u_f(κ) is the U_ℓ-number of the
  family at the weight κ, an ℓ-adic analytic function on the weight disc that contains ℓ − 3 and −2.
- The twist (3|ℓ) is Calegari's χ₁₂ at ℓ.

## Statement A: the constants (13 ≤ ℓ)

Let α_f = (3|ℓ)·u_f(−2), and P_ℓ(x) = ∏_{f ∈ S_ℓ} (x − α_f), a polynomial with ℓ-adic integer coefficients.

**A1 (conjectured; the frame is proven).** Σ_i [xⁱ]P_ℓ · H_(b+2i)(n) ≡ 0 modulo a power of ℓ that grows without bound
with b, for every n. That is, the operator that moves the tower's limit two levels up satisfies P_ℓ. This recurrence is
the main claim. P_ℓ is that operator's characteristic polynomial exactly where the tower has full rank R_ℓ: every
ℓ ≤ 1297 except 607, by Boylan–Webb's computation (NOTE.md §5). *Reworded 20 Nov 2026: an earlier version called the
recurrence and the characteristic-polynomial statement equivalent. They are not; the recurrence alone gives less.*

**A2 (counted).** dim S_ℓ = R_ℓ for every ℓ from 13 to 89 (exactly, over the integers) and to 1297 (on lane ℓ, the exact
class where forms of lower level interfere: 17, 29, 499). Every form of S_ℓ is ordinary for all ℓ ≤ 1297
(`partitions_three_forms.txt`, `partitions_flint_exceptions_new_37_1297.txt`).

**Evidence for A1 (counted):**
- The digits of α_f from the families agree with every digit that exact partition numbers settle: 13 digits for 13; 11 for
  17, 19, 23; 5 for 29, 31 (`partitions_flint_deeper.txt`).
- The law Σ [xⁱ]P_ℓ·H_(b+2i) ≡ 0 holds on lanes ℓ, ℓ², ℓ³ at levels (1,3,5), (3,5,7), (5,7,9) for 37, 41, 43, 47, 53 and 59,
  and on lanes 61, 61² for 61. Each time, P_ℓ came from the forms alone (`partitions_flint_arbiter.txt`,
  `partitions_hida_more_check.txt`, `partitions_flint_family_check_53_59.txt`).
- 73 and 97 (13 Nov 2026): the classes have 3 and 4 forms, found exactly over the integers. 73's law holds on lane 73 at levels (1,3,5,7) for n = 0 … 9; the laws along Q hold on 12 of 12 chains at Q = 5, 7, 11, 13 for both. P₇₃ and P₉₇ came from the forms alone, and no shorter law fits (`partitions_73_97_check.txt`).
- Boylan–Webb's printed relations for 53 (lanes 53 and 53²) and their periods c₃₇ = 36, c₄₁ = 10 are exactly what P₅₃,
  P₃₇, P₄₁ give (`partitions_bw_check.txt`).
- The digits themselves: 30 for 13 to 31 and 20 for 37 to 61 (`partitions_flint_family_ℓ*.txt`; the README, the eighteenth
  time).

**A3 (derived, mod ℓ only, from published results: AAT (2.4), (2.6) and (2.8), with AAD (1.4); written out in PROOF.md, "The first digit, written out". Counting the notes with multiplicity also uses Yang 2014 as a one-to-one correspondence.)** The first digit: α_f ≡ (3|ℓ)·a_ℓ(f) (mod ℓ).

**A4 (counted at 13 and 17; from 14 Nov also at every prime 29 … 613, below). What fades is the next family.**
- *14 Nov (README "The share condition, tested on partition numbers"):* at every prime 29 … 613 the part of the tower two
  digits deep at level 3 is the first-step families (Hecke matrices of partition numbers, Smith form). They were named by
  their Hecke numbers at 29 … 199, and every one of them is present at 29 … 397. Every family starts as a unit, except
  one slope-4 family at 59.
- Under its ordinary components, the tower carries a component of a slope-2 family g of level-6 forms with the class's Atkin–Lehner signs, in another system on lane ℓ.
- That component is multiplied by (3|ℓ)·u_g(−2) every two levels: 13 → 13²·(10, 11, 6, 2, 10, 10, …) and 17 → 17²·(3, 9, 1, 13, 4, 3, …).
- Measured on partition numbers as the fading part of the tower, these agree in 6 of 6 digits at each prime (`tower_fading.txt`, `fading_carry.txt`).
- Its valuation (2) is the towers' speed: one digit per level.
- At 19 and 23 the fading part is two such notes, two slope-2 families in two systems (4 of 4 digits each, counted).
- At 29 and 31 it is one slope-1 family, found to the two digits two weights give (23, 3 and 2, 28; the third weight was
  stopped by low memory).

## Statement B: the steps along another prime Q

**B (conjectured; the frame is proven).** Write the limit as a sum of components c_f, one for each f in S_ℓ. For each
prime Q ≥ 5, Q ≠ ℓ, and M = 24N − 1 with ℓ | M and Q ∤ M:
c_f(Q²M) = (λ_f(Q) + ε·δ)·c_f(M),  with  λ_f(Q) = (3|Q)·a_Q(F_f at weight −2),  ε = (M|Q),  δ = −(−3|Q)·Q⁻².
For two forms this gives the two-step law c(Q⁴M) ≡ (t + εδ)·c(Q²M) − (Q⁻³ + δ² + tεδ + n)·c(M), with t and n the sum and
product of the two λ_f(Q).

**Evidence (counted):** on lane ℓ, 48 of 48 chains for each of 37, 41, 43, 47, and for 61 the three-step law on 48 of 48
(`partitions_split.txt`, `partitions_three_check.txt`). For 13, 17, 19 at Q = 5, 7, 11, the multiplier measured at level
b shares about b digits with λ_f(Q) of the weight −2 member, but one digit with the classical f
(`partitions_hida_steps.txt`).

**B2 (derived from B; counted at 97). The components can be taken apart.** If h is a polynomial whose roots are the
λ_f(Q′) of some forms of the class, then h applied to the step at Q′ cancels those forms. What is left is the tower of
the other forms alone, with their own laws.
- At 97, the class is one form g and a conjugate triple on lane 97. With h the triple's cubic at Q′ = 5,
  G(M) = Σ κ_j·c(5^(2j)·M) is g's tower alone.
- Along 97: G(97²M) ≡ 28·G(M), with 28 = (3|97)·a₉₇(g), the first digit of the 97-adic constant. Holds 20 of 20,
  levels 1 and 3, numbers below 3·10¹¹.
- Along every other Q: G(Q²M) ≡ (λ_g(Q) + εδ)·G(M). Holds 60 of 60 for Q = 7 … 23.
- About 3 primes in 97 are lucky for g alone (53 is the first); their zeros and exact controls held 18 of 18
  (`partitions_97_split.txt`).

**B3 (derived from B; counted at 37 to 97). Every form can be heard alone.**
- **Steps with ℓ-adic numbers:** at the twins (37, 41, 43, 53), cancelling the other form leaves one form's tower.
  Its multiplier between levels b and b + 2 agrees with its own constant in at least (b + 1)/2 digits and with the
  other's in none (147 of 147; `partitions_twins_split.txt`).
- **Conjugate groups** (pairs at 47, 59, 61; triples at 73, 97): the same with steps whose numbers come from the
  group's own field (degree 2 or 3). λ₁(5) is paired with the constant through T₅ = R(T_ℓ) on lane ℓ.
- The same rule holds there, 108 of 108, and no member ever agrees with a conjugate (79 of 79;
  `partitions_braid_split.txt`).
- 97's constant shows 28, 26 at levels 3 → 5 (`partitions_97_split_digit2.txt`).

## Statement C: 5, 7 and 11

**C (conjectured; for 5 the slope is printed by Calegari).** For ℓ = 5, 7, 11, S_ℓ is empty and the tower follows the
smallest positive slope: a Coleman family of slope s = 2, 1, 2 in the tower's system on lane ℓ, through the weights
k = −2 + (ℓ − 1)·ℓ·(x + 1), with α_ℓ = (3|ℓ)·u(−2)/ℓ^s. The laws' exponents (Watson's, Atkin's) are this slope.

**Evidence (counted):** from forms alone on FLINT, 33 digits for 5, 20 for 7 and 12 for 11, with exactly one root of
slope s at every weight used (`partitions_flint_slopes_*.txt`). Exact partition numbers confirm 33, 20 and 7 of them
(`partitions_constants.txt`, `partitions_flint_lawprimes_check.txt`).

## Statement D: the rank, and the exception 607

**D1 (derived).** On lane ℓ the levels act on the class as (3|ℓ)·T_ℓ does (A3). So the space they fill has at most as many
directions as T_ℓ's minimal polynomial on S_ℓ (on lane ℓ) has degree.

**D2 (counted).** For every ℓ from 37 to 1297 except 607, that minimal polynomial has degree R_ℓ: no zero root, no repeated
root. At 607 it has degree 24 = R₆₀₇ − 1, since two forms share their number at 607 and T₅ tells them apart
(`partitions_flint_607_span_*.txt`, `partitions_flint_exceptions_new_37_1297.txt`).

**D3 (conjectured).** r_ℓ(1), the rank on lane ℓ, equals that degree. This holds for every ℓ ≤ 1297, against Boylan–Webb's
computed ranks: R_ℓ everywhere except r₆₀₇ = 24 (quoted). What happens at 607 on lanes 607², 607³, … is not settled here.

## Statement E: the lucky primes

**E1 (derived, once a_Q is known on lane ℓ, from AAT's Lemma 5.1, (2.4) and (2.6); written out in PROOF.md, "E1, the laws at a lucky prime, written out".)** If (3|Q)·a_Q(f) ≡ ε·(−3|Q)·Q⁻² (mod ℓ) for every f in
S_ℓ, then p(N) ≡ 0 (mod ℓ) for every N with 24N − 1 = Q²·M, ℓ | M, (M|Q) = ε. If (3|Q)·a_Q(f) ≡ 0 for every f, then
p(N) ≡ 0 (mod ℓ) whenever 24N − 1 = ℓ·Q³·m with Q ∤ m. So every lucky Q we list is a congruence, by a finite
computation: 62 to 174 of them below 10⁶ for each of 37 to 59 (`partitions_flint_lucky_*.txt`). Checked on exact partition
numbers as well: 2999, 4049, 8941, and 36 of 36 new ones, two of each kind for each of 37 to 59
(`partitions_flint_lucky_check.txt`).

**E2 (derived, from Chebotarev's density theorem, if each piece's mod-ℓ Galois image is full and the pieces independent; the
counts agree).** The share of primes Q with a law of a given kind with target c is the average over Q mod ℓ of
∏ over the pieces of S_ℓ on lane ℓ of (q + χ_q(c² − 4Q⁻³))/(q² − 1), with q = ℓ^(degree of the piece). That is about
3/ℓ^R_ℓ for all three kinds together. Counted to 10⁶ for 13 primes ℓ: χ² = 37.7 on 36 counts, and the ±1 correction is
confirmed (Σ(O − E)²/E = 55.0 on 48 cells, against 239.2 without it).

*Added 20 Nov 2026 from a reader's review, recomputed (`recheck/scripts/week2/e2_chi.py`, `e2_q5.py`):*
- *The 55.0 replaces 58.0, which the saved counts do not reproduce with this statistic (dividing by E(1 − p) gives 57.9).*
- *At 37 both "ifs" follow from print and one check. Ahlgren–Allen–Tang, Proposition 3.3 (arXiv:2112.09481): each image
  contains SL₂(𝔽_ℓ). Their Lemma 3.7(2): the joint image over ℚ(ζ_ℓ) is the full product if the projective
  representations differ. The check: a₅ ≡ 5 and 36 (mod 37), so a₅² ≡ 25 and 1 with the same determinant. The share
  for 37 is then 24,635/11,228,544, one in 455.8: 172.2 expected to 10⁶, 174 counted. Not refereed.*
- *At all 20 primes 13 … 97 the values a₅²/5^(ℓ−4) are distinct modulo ℓ and generate each piece's full field (computed).
  The density argument for pieces of degree above 1 is not written out, so E2 stays conditional there.*

**E3 (derived from B; counted). At a lucky prime nothing is left to chance.**
- Class law with ε: the other class has c(Q²M′) ≡ 2ε·(−3|Q)·Q⁻²·c(M′).
- Ono's kind: c(Q²M) ≡ −(M|Q)·(−3|Q)·Q⁻²·c(M) in both classes.
- Checked:
  - 14 of 14 controls at the laws of 37, 43 and 53;
  - 6 of 6 at the new laws of 73 and 61;
  - 16 of 16 Ono-type values at 73 and 61 (`partitions_flint_lucky_check_73.txt`, `partitions_flint_lucky_check_ono.txt`).

**E4 (counted): three forms.** Every prime Q < 5·10⁷, with the lighter test of `partitions_flint_lucky_dual.py` (a_n(f) =
a_1(T_n·f), so T_Q is one number on the class exactly when a_Q(S_j) ≡ c·a_1(S_j) for its basis; it agrees with the
whole-matrix test, 39 of 39):
- **73** (one piece of degree 3): 25 laws against 23.1 expected.
- **61** (pieces 1 and 2): 35 against 39.6 (χ² = 2.3 on the 6 counts).
- The first congruence modulo 73: p(107239327321929851108449·n + 9040523964878168) ≡ 0 (mod 73) for every n
  (107239327321929851108449 = 73·11367817³).
- Checked on exact partition numbers: 2 class laws for 73 and 4 for 61, 12 of 12 members and 6 of 6 controls
  (`partitions_flint_lucky_dual_73_61_max50000000.txt`, `partitions_lucky_progressions.txt`).
- **67, 71, 79, 83:** 92 laws against 89.3 expected (χ² = 7.2 on 12 counts). All six primes with three forms: 152
  against 151.9.
  - Checked: 6 class laws, 12 of 12 members and 6 of 6 controls (`partitions_flint_lucky_check_67_83.txt`).
  - Example: p(19728144816847346813·n + 6177904563614) ≡ 0 (mod 67).

## Statement F: the loops

**F1 (derived from A1).** On lane ℓ^m the tower comes back to itself after c_m two-level steps, c_m the order of the
constants (their companion matrix) modulo ℓ^m.

**F2 (counted, 5 to 61).** c_m = c₁·ℓ^(m−1) for every m the constants reach: no tower comes back early at a deeper digit.
The shortest loops close on exact partition numbers: 5 and 31 after 4 levels, 11 after 10, 7 after 12
(`partitions_loops.txt`).

**F3 (counted to 700; the expectations derived).** Across primes the loop lengths are distributed as for random units:
44 primitive roots among 119 constants on lane ℓ against 45.7 expected (`partitions_loops_stats.txt`).

## Statement G: other numbers of colours (r-colored partitions, ∏(1 − qⁿ)^(−r))

**G1 (counted for r = 1, 3, 5; conjectured for odd r).** The tower's rank bound R_ℓ(r) (Belmont, Lee, Musat and Trebat-Leder,
QUOTED) equals the number of newforms of weight ℓ − r − 2 in one Atkin–Lehner class. That class is on Γ0(6) when 3 ∤ r,
with U₂ = ε·(2|ℓ)·2^((k−2)/2) and U₃ = ε·(3|ℓ)·3^((k−2)/2) (ε = +1 for r = 1, −1 for r = 5). It is on Γ0(2) when 3 | r, with
U₂ = (2|ℓ)·2^((k−2)/2) (r = 3). Checked for r = 1 to 1297, r = 3 to 251 and r = 5 to 157 (`colors3_census.txt`,
`colors5_census.txt`).

**G2 (counted for r = 3, 5; conjectured).** The constants are ψ_r(ℓ)·u_f(−(r + 1)), with u_f the family through each class
form and ψ_r the character of η^r: (3|ℓ) for r = 1 and 5, (−1|ℓ) for r = 3, the character of η³'s theta series. Checked
digit by digit against exact r-colored partition numbers: for r = 3 at 13, 19, 23, 29, 41 and the three-term law at 31, 37;
for r = 5 at 13, 17, 19, 23, 29, 31 and the three-term law at 37, 41 (`colors_check_*.txt`, `colors3_digits_check.txt`).

**G3 (counted).** Where the class is empty, laws come from families of positive slope:
- r = 3, ℓ = 17: one form of slope 1 (the law 17^⌈b/2⌉);
- r = 3, ℓ = 11: two forms of slope 2 whose constants agree on lane 11, so a three-term law holds instead of a constant. It
  is checked to 11⁶ on exact values;
- r = 5, ℓ = 11: the law 11^(⌊b/2⌋+1).
  - It follows one family of slope 1, with constant (3|11)·u(−6)/11 = 2, 8, 5, 2, 10, 5, 2, 0, 1, 1, 3, … (11-adic).
  - Checked to every digit the exact numbers hold (`colors5_slope11.txt`, 12 Nov).

The pair (5, 23) has a fixed factor 23 at every level and one constant beyond it (Boylan's "superexceptional"
case).

**G4 (derived; the inputs counted exactly).** The classes of 1, 3 and 5 colours are ordinary at every prime where they
exist: ℓ ≥ 7, 11 and 11.
- The theta cycle sends a non-ordinary class form of weight ℓ − r − 2 to a form of weight r + 5 with fixed U₂ (and U₃)
  values, and none exists, except modulo 2, 3 or 5.
- So the colour laws (at 11 and 17 for 3 colours, at 11 for 5) come from families of positive slope at higher weights.
- Sources: `second_kind.py`, `second_kind_colors.txt`, NOTE.md §3.

## Statement H: other forms (24 colours, spt, j, 1/E₄, 1/E₆, 1/E₄², a pole at j = 2)

**H1 (counted).** Let F have level 1 and weight w, with its only pole at the cusp (1/Δ: w = −12; j: w = 0). Then the
tower F | U_ℓᵇ satisfies Σ [xⁱ]P·H_(b+i) ≡ 0, ever deeper. P is the polynomial of the U_ℓ-values at weight w of the ordinary
families of level 1: the Eisenstein family's 1 and the cusp forms'.
- At w = 0 the Eisenstein member is the constant 1, so for n ≥ 1 the cusp forms alone suffice. At 5, 7 and 11 that is
  Lehner's congruences.
- Checked at every prime 5…43 the data reach (`colors24_check.txt`, `level1_check_j_….txt`, `level1_j_cusp_….txt`).

**H2 (quoted).** For spt, two levels multiply the tower by (3|ℓ)·r, where r is a root of x² − (1 + ℓ)x + ℓ, the Eisenstein
series of weight 2.
- r = 1 for ℓ ≥ 11: Ahlgren–Bringmann–Lovejoy (2011).
- r = ℓ at 5, 7 and 13: Garvan (2012).
- The measured depths equal their moduli (`spt_tower_….txt`, `spt_check.txt`).

**H3 (counted; conjectured in general).** A pole inside: 1/E₄, with w = −4 and the pole at ρ.
- ℓ ≡ 2 mod 3 (ρ supersingular): the forms' P alone.
- ℓ ≡ 1 mod 3 (ρ ordinary): P(x)·(x − u^(w−2)), where u is the unit root of Frobenius of y² = x³ + 1 at ℓ. u⁻⁶ is
  measured to 9 digits at 7, 5 at 13, 3 at 19, and 1 at 31 and 43. 37 cannot tell (`level1_pole_check_invE4_….txt`).
- A second curve: 1/E₆, with w = −6 and the pole at i on y² = x³ − x. The prediction u⁻⁸ at ℓ ≡ 1 mod 4 was written
  before the run. Measured: u⁻⁸ to 10 digits at 5, 4 at 13, 3 at 17, 2 at 29, 1 at 37 and 41. At ℓ ≡ 3 mod 4 the forms'
  P alone (`level1_pole_check2_invE6_….txt`).
- A point with no complex multiplication: 1/(E₄³ − 2Δ), w = −12, the pole at j = 2 on y² = x³ + 10356x + 11916304.
  Predicted before the run: u_C⁻¹⁴, u_C the unit root of that curve at ℓ. Measured at 5, 7, 11, 13, 17, 19, 23; P alone at
  29 and 31, where the curve is supersingular. Read blind, the tower gives j = 2 back: it is the smallest of the two integers
  |j| ≤ 10⁶ that fit (`pole_j0_check_….txt`).
- Repeating and drifting (`repeat_check_A.txt`). U_ℓ moves a pole toward the canonical lift of its curve (a CM curve).
  A pole put on that lift (j = −32768 at 5, 54000 at 7 and 19, −3375 at 11) gives u⁻¹⁴ exactly. With the slope-1 families
  in the law, its depth grows 2 digits per level, while j = 2's grows 1. (Without them, 7, 11 and 19 grow 1 per level either
  way.)
- Two CM poles that Frobenius swaps (`repeat_check_B.txt`): 1/(Δ²·H₋₁₅(j)), weight −24. At 17 and 23 the two curves of
  discriminant −15 swap. The remainders' ratio is constant only two levels apart (U⁻²⁶, U the unit root over 𝔽_ℓ²), and only
  the pair law P(x)·(x² − U⁻²⁶) holds (with the slope-1 families, depth 6 of 7 at 17 to n < 10⁸). At 19 and 31 each curve
  is fixed: P(x)·(x − u⁻²⁶) (with the families, the whole lane at 19). The supersingular primes need P alone.
- Three CM poles in a 3-cycle (`repeat_check_C.txt`): 1/(Δ³·H₋₂₃(j)), weight −36. At 13 the remainders' ratio settles three
  levels apart, at U⁻³⁸ (U the unit root over 𝔽₁₃³), and P(x)·(x³ − U⁻³⁸)·S(x) holds (depth 5 where the data allow one
  check). At 29 the ratio settles neither one nor two levels apart. At 59 each curve is fixed (U⁻³⁸ one level apart). The
  cycle length is the order of ℓ's prime in the class group.
- The same at small primes (`repeat_check_C_D31_first.txt`): 1/(Δ³·H₋₃₁(j)), whose three curves cycle at 5, 7, 19 and 41. At
  5 the ratio settles three levels apart at all six levels (U⁻³⁸), and P(x)·(x³ − U⁻³⁸)·S(x) reaches the whole lane 5¹¹ from
  level 2. At 7, P(x)·(x³ − U⁻³⁸) gains two digits a level to the whole lane 7¹⁰. A single cube root fails everywhere.
- A double pole at ρ: 1/E₄², w = −8. The law P(x)·(x − ℓ·u^(w−4))·S(x), with S the slope-1 families at weight −8 computed
  from the forms (Lauder's algorithm, `slope1_families_…txt`), has no free constant. It holds to the whole lane at 7, 13 and
  19 from level 1 or 2 (`double_pole_full_law_….txt`). At weight 4 (E₄/j) the same formula gives Zhang's (−3|ℓ)·ℓ.
- Conjecture: a simple pole at a point whose curve is ordinary at ℓ adds the constant u^(w−2), u the unit root of that
  curve's Frobenius. It holds exactly at CM points, and digit by digit elsewhere. A double pole at ρ adds ℓ·u^(w−4). The derived
  heuristic is that T_ℓ multiplies the coefficient of (τ − τ₀)^(−j) by π^(j−1)·π̄^(w−1−j) + π̄^(j−1)·π^(w−1−j).

**H4 (counted): the boundary.** The data-only search finds a short law on lane ℓ in every modular function tried. That
includes three never computed with forms: distinct parts, overpartitions and the moonshine function T_2B, at 5, 7 and 11.
It finds none in partitions into squares or cubes, which have no modular origin, at 5, 7, 11 or 13: up to degree 6 where
the levels allow, 4 or 5 where they run short (`group_members_….txt`, `boundary_powers_….txt`).

## What a proof of A, B and C needs (derived outline)

*Superseded by `PROOF.md` (12 Nov 2026), which follows the published statements step by step (sources in
`research_notes/Proof ingredients/ingredients.md`), names the open links and gives a route for each. The outline below is
kept as it was first written.*

Everything below exists in print as a framework, and specialists would read A1 as what it predicts. What is missing is the
assembly, with exact normalizations:
1. **The two-level map as U(ℓ²).** Folsom–Kent–Ono and Calegari write the levels L_ℓ(b) as weight −1/2 overconvergent
   forms, and the step from b to b + 2 as an operator whose limit is the slope-0 part (FKO's appendix, Ramsey–Calegari).
   Needed: the exact operator and its normalization against H_b.
2. **Ramsey's ℓ-adic Shimura lift** carries the weight −1/2 eigenforms to weight −2 on Γ0(6), twisted by χ₁₂ (Calegari,
   Lemma 4.6.4), and is Hecke-equivariant (T(Q²) ↔ T_Q). Needed: that it carries U(ℓ²) to U_ℓ with the factor (3|ℓ), and the
   class's Atkin–Lehner signs.
3. **Hida's control theorem** at the weight −2. The slope-0 part of overconvergent forms of weight −2 (level 6, the class's
   signs) is the specialization at −2 of the ordinary Λ-adic forms through the ordinary stabilizations of S_ℓ. So its
   U_ℓ-numbers are the u_f(−2). For 5, 7, 11 the same with Coleman families of slope 2, 1, 2 (Coleman, Gouvêa–Mazur; the
   slopes are checked, not proven, beyond 5).
4. **The rank.** Boylan–Webb's computation (r_ℓ(1) = R_ℓ for ℓ ≤ 1297 except 607) supplies it where needed.
With 1–4, A1 and B follow for every ℓ ≤ 1297 except possibly 607, and C for 5 (7 and 11 need their slopes).

## What would make the statements fail (and did not)

- A form of lower level sharing the class's numbers on lane ℓ: it happens (17, 29, 37, 499), and is set aside exactly by its
  U₂, U₃ values (Deligne's bound forbids exact equality). It never changes a constant.
- A non-ordinary form in S_ℓ: none up to 1297. **Never, for any ℓ ≥ 7** (12 Nov; derived, `second_kind.py`):
  - Theta cycles turn a non-ordinary class form of weight ℓ − 3 into a weight-6 form on Γ0(6) with U₂ ≡ 4, U₃ ≡ 9, and
    none exists.
  - The mechanism's prediction for the neighbouring class came true: h₆ (the weight-6 newform, U₃ = −9) is non-ordinary
    at 859, and at 859 exactly that class has a non-ordinary form.
- Two forms with the same a_ℓ on lane ℓ: once up to 1297 (607). It lowers the rank by one, as Boylan–Webb found.
