# What is left, and where to dig next (12 Nov 2026)

What we found is in FINDINGS.md; the proof route is in PROOF.md. This page holds what remains open.

## 1. Finishing the proof (the classical route of PROOF.md)

| link | what it needs | status |
|---|---|---|
| 1. Each level is classical | Chen–Du–Hou–Sun, Thm 3.1. F_{ℓ,b,r} ≡ η(24z)^γ·φ(24z) (mod ℓᵇ), with φ ∈ M_λ(SL₂(ℤ)) and γ = (24β − r)/ℓᵇ, β ≡ r/24 (mod ℓᵇ), so 0 < γ < 24 for odd r < ℓᵇ. λ = ((ℓᵇ − ℓᵇ⁻¹)/2)·r − (γ + r)/2 for odd b, and (ℓᵇ − ℓᵇ⁻¹)·r − (γ + r)/2 for even b | **read and confirmed 12 Nov** (arXiv 1206.6642 v1, the only version). **Checked on numbers 12 Nov:** at 8 pairs (ℓ, b) up to (13, 3), 2,942 of 2,942 coefficients agree (`cdhs_check`). K_b = (ℓ − 1)ℓ^(b−1) − 2 (odd b) or 2(ℓ − 1)ℓ^(b−1) − 2 (even b) |
| 2. The step is exact | P(b + 2) = P(b) \| U(ℓ²) | classical |
| 3. U(ℓ²) ≡ T(ℓ²) | Yang's formula T(p²): a(p²n) + (12/p)((−1)ᵏn/p)pᵏ⁻¹a(n) + p²ᵏ⁻¹a(n/p²). The difference from U(p²) is divisible by pᵏ⁻¹, and k grows with b | follows from the quoted formula |
| 4. The eigenvalues | Ahlgren–Andersen–Dicks, Adv. Math. 447 (2024), Thm 1.1: the η-multiplier lift S_t maps into the new space on Γ0(6) with signs (−(8/r), −(12/r)), and S_t(T(p²)F) = (12/p)·T_p S_t(F) for EVERY prime p ≥ 5. With T(ℓ²) semisimple, the tower's minimal polynomial divides the class polynomial (NOTE.md §4.3) | **closed 12 Nov** (published theorem, every ℓ). Yang's Thm 1 (whose proof omits six classes mod 24) is no longer needed. Boylan–Swati (arXiv 2505.01018) also reprove Yang's result |
| 5. Only the ordinary part survives, and its values move continuously | Hida's theory at classical weights (Hida 1986, Wiles 1988). The ordinary part has bounded rank along the weights. Its Λ-adic characteristic polynomial makes A(k) ≡ A(k′) (mod ℓᵐ) when k ≡ k′ (mod (ℓ − 1)ℓᵐ⁻¹) | **cited (12 Nov):** Hida, Ann. Sci. ENS 19 (1986) Thm 3.1, Cor 3.2; Annals 128 (1988) Thm 3.2; via Ghate–Kumar. The θ-filtration fact is pinned via de Shalit–Goren's restatement |
| 6. The statement | Modulo every ℓᵐ, the two-level step's characteristic polynomial on the tower divides ∏_f (x − (3\|ℓ)·u_f(−2)), read off the Λ-adic polynomial at weight −2. Equality holds exactly at full rank (Boylan–Webb: every ℓ ≤ 1297 except 607). This form also covers congruent forms and non-ordinary ones | **drafted 12 Nov: NOTE.md** (Theorems A and B, every ingredient cited, Theorem A proved in full) |
| 7. A specialist's check | of a complete written note | the user decides whom to ask |

**Link 4, checked further (12 Nov).** The published text is Nagoya Math. J. 215 (2014) 1–66,
https://doi.org/10.1215/00277630-2684452 (Cambridge Core; the PDF downloaded). It has the same gap:
- p. 58: "Here we work out only the cases n ≡ 1 mod 24 and n ≡ 11 mod 24 and omit the proof for the other cases";
- p. 64: "We skip the proof of the other cases";
- Theorem 2's proof is "omitted".

What our route needs from Yang is the trace identity at n = ℓ, ℓ², ℓ³, … only. That fixes the characteristic polynomial
of T(ℓ²), through the Hecke recursion; derived. When ℓ ≡ 1 or 11 (mod 24), every power of ℓ is ≡ 1 or 11 (mod 24), so
**for those primes (59, 73, 83, 97, 107, 131, …) link 4 rests on the cases Yang wrote out.** For the other six classes,
someone has to do the omitted cases, following his lemmas (they include values for every class of n mod 24).

On this route (our reading) the overconvergent chain is no longer needed:
- Calegari's Lemma 4.6.4, Ramsey's lift, half-integral multiplicity one, and cuspidality all drop out.
- Multiplicity one is Yang's, at each classical weight.
- Each level is a holomorphic cusp form, since γ > 0.
- 607 is part of the statement, not a gap.

For 3 colours (3 | r) the same route needs Yang's Thm 2, whose proof the paper omits entirely. For 5 colours, Thm 1 applies.

**What we can do ourselves:**
- the numerical check of links 1 and 3 at small b;
- precise citations for link 5;
- the complete note (link 6).

**What needs outside help:** only link 7, the specialist's check.
- Link 4 was closed by Ahlgren–Andersen–Dicks (12 Nov).
- The routine points (Katz on Γ0(6), Hida on one Atkin–Lehner class, the continuity modulus, semisimplicity) are written
  out in NOTE.md §8.
- **The draft proof has no known gap.**

## 2. Where to dig next

Sizes: **quick** (minutes), **run** (hours on this PC), **WSL** (needs the Linux FLINT for long series), **outside**
(needs papers or people).

**The partition towers**
1. The next 607-type prime after 1297, found from the forms alone. **run**
   Realistic odds (derived from the same 1-in-ℓ chance per prime): about 0.7 such primes were expected between 37 and
   1297, and one (607) was seen. The chance of another is about 10% by 3,000, 22% by 10,000, and 38% by 100,000.
2. The first class form that is not ordinary: the "second kind" of exception. There is none up to 1297. **Reduced
   (12 Nov), see `second_kind.py`.**
   - Theta cycles turn such a form into one of weight 6 on Γ0(6) with U₂ ≡ 4 and U₃ ≡ 9.
   - No form of weight 6 has those values; the only newform, h₆, has U₃ = −9.
   - So the partition class is predicted to be ordinary at every prime.
   - h₆ is non-ordinary at 859. There the class (−(2|ℓ), +(3|ℓ)) is predicted to hold a non-ordinary form.
   - **Done 12 Nov: confirmed exactly.** At 859 that class, and only that class, has a non-ordinary form. Nothing at
     13–300 or at the controls 853, 857, 863.
   - **The second kind never happens for p(n).** Left: the argument written out carefully, for the specialist's check.
3. Why a loop's return depth is always 1 (STATEMENT F). **Answered (already in README, 20th dig): chance.** One early
   return was expected among 5–61 and none was seen, as for the primitive roots. Not a law.
4. The leftover extra multiples: 20.64% of p(N) divisible by 5 and 14.46% by 7 where no law reaches, and 8.94% for 13
   in its class. Not class numbers. Delaunay's Tate–Shafarevich odds come close (0.2067, 0.1454): a resemblance only, so
   test it. **run**
5. Why the constants settle at 5, 3, 1 and 2 digits per two levels (Dig 1, fourth time). **Partly answered (12 Nov),**
   from the slopes already computed (`partitions_flint_slopes_*.txt`):
   - At 5 the law's family has slope 2 and the next family, near weight −2, slope 7: a gap of 5, the measured speed.
   - At 7 the slopes are 1 and 7 (gap 6), and at 11 they are 2 and 4 (gap 2). The measured speeds 3 and 1 are exactly
     half the gap. Why half: open. **Answered (14 Nov):** the quiet forms at 7 and 11 sit in other systems, with slopes 4
     and 3. Their gaps to the law's slopes, 3 and 1, are the measured speeds. The "half" came from taking the next slope in
     the law's own system. And from 13 on, the rate is fast exactly when the class gains no form at the first weight step
     (README "What sets the settling rate").
   - 13's next slope is not yet computed.
   **Measured (13 Nov):** 13, 17, 19 and 23 add one digit per level (13's is Atkin–O'Brien's congruence, their theorem
   of 1967). 29, 31 and every note from 37 to 97 add one digit per two levels, Boylan–Webb's rate (README "The
   speed of the towers"). The guess "non-ordinary room in weight ℓ − 1" was tested and is wrong: every form there is
   ordinary for 13 ≤ ℓ ≤ 400. Why the four primes 13 to 23 are fast: open. **Answered (14 Nov, FINDINGS 22):** the
   first weight step adds exactly ⌊ℓ/24⌋ forms, the level-6 newforms of weight ℓ − 5 with signs ((2|ℓ), (3|ℓ)). They
   alone give slope-1 families, which make a tower slow. So it is fast exactly for ℓ < 24 (checked at every prime to
   613).
   **Narrowed (13 Nov, `tower_fading`):** the speed is the size of the fading part (ℓ² per two levels at 13–23, ℓ at 29
   and 31). It is one note at 13, 17, 29 and 31, with measured constants (13: 169 × (10, 11, 6, 2, 10, 10, …)), and a
   mix at 19 and 23. **Found at 13 and 17 (13 Nov):** a slope-2 family with the class's signs in another system,
   6 of 6 digits (STATEMENT A4). 19 and 23 (13 Nov): two slope-2 families each, 4 of 4 digits. 29 and 31 (13 Nov): one slope-1
   family each, 2 of 3 digits from two weights. Left: their third weight (2,434 and 2,788), which was stopped by low
   memory. A single run at a time would use about half the memory of the two together.
6. 5 colours at 11: the family behind the law 11^(⌊b/2⌋+1). **Done (12 Nov):** one slope-1 family; constant
   (3|11)·u(−6)/11 = 2, 8, 5, 2, 10, 5, 2, 0, 1, 1, 3, …; the exact numbers agree to every digit they hold (`colors5_slope11.txt`).
7. The lucky-prime density as a theorem: Chebotarev plus a large Galois image for the class forms (Ribet). **outside**
8. j at 13: read Atkin–O'Brien (1967), which the AMS site blocks for automated reads. **outside** (user's access)

21. Lucky primes at the other primes with three forms (67, 71, 79, 83) to 5·10⁷, with `partitions_flint_lucky_dual.py` (4 minutes
    each): more counts for the 3-in-ℓ³ law and explicit congruences. **Done (13 Nov):** 92 laws against 89.3 expected; 6 checked on exact partition numbers (README).
22. 97's third digit (81) and the twins' fifth: the next level pair needs p(N) near 10¹⁷ to 10¹⁸. **out of reach here**
23. A conjugate group (73's triple, 47's, 59's and 61's pairs, 97's triple) cannot be split by steps with ℓ-adic
    numbers. **Answered (13 Nov):** steps with numbers from the group's own field split it; every member shows its own
    constant and none of its conjugates (README "The braided forms taken apart", STATEMENT B3).
24. Why the forms of level 6ℓ grow by one or two per step in each system, and at which step a form reaches level 6.
    Together these would give every slope in advance (FINDINGS 18). Both are dimension counts. A start: the count per
    system at level 6ℓ against the count at level 6 over many steps (`level6l_counts.py`). The second forms of a system
    (the ghost's next zeros) are the next test. **quick** to **run**
    **Answered (14 Nov, FINDINGS 19):** one seed per system (a Serre weight); 2 per step for a seed of weight 4 … ℓ, 1
    for weight ℓ + 1 or an Eisenstein-type seed of weight ℓ − 1. Level 6 is the same seeds over Sym^(k−2). 132 of 132
    counts. Left: the full ghost series from the seeds against every slope, second forms included. **Done (FINDINGS
    20):** 225 of 225. Left: whether the proofs' genericity condition covers these systems (**outside**: the papers).
    **Done (FINDINGS 21):** 21 of 39 families are covered (very generic, ℓ ≥ 11). 18 are not (8 Eisenstein, 8 Steinberg,
    2 with a = ℓ − 3), and all 39 fit. Left for a specialist: that the class with fixed signs at 2 and 3 gives a
    primitive projective module, so that the local theorem applies to it directly. **Not needed (14 Nov):** the class
    is recovered from the whole families at levels 6ℓ, 3ℓ, 2ℓ, ℓ by subtraction, each covered by Theorem 1.3; the
    signs are fixed within a family (`sign_check.py`). And the 8 Eisenstein families are proved Eisenstein
    (`eisenstein_proof.py`).
25. The share condition: does the tower carry its slowest quiet notes? **Answered (14 Nov, FINDINGS 23):** yes, at
    every prime 29 … 613 (Hecke matrices of partition numbers, `share_test.py`). Every family starts as a unit, except
    one slope-4 family at 59. Each step then costs ⌈s/2⌉ or ⌊s/2⌋ digits. Left:
    - **A1 for the quiet notes** (14 Nov): read each tower's
      fading part from partition numbers at 29, 31, 37, … and match its multiplier, digit by digit, with the quiet
      family's number at weight −2 from the forms. So far it agrees to 6 digits at 13 and 17, 4 at 19 and 23, and 2 at
      29 and 31; the third weight at 29 is running. **run**
    - Whether the starting shares have a formula (an L-value?) that says when one is divisible by ℓ; 59's slope-4 family
      is the test case. **outside** or **run**
    - 23's shifted grid: one slot at level 1 is a digit deeper, unresolved (second families in play). **quick**

**The general law**
9. Constants from forms for distinct parts, overpartitions and T_2B, whose laws the data alone found. That means
   level-2 and level-16 forms. **run**
   - **T_2B done on lane ℓ (12 Nov):** its constants are the Eisenstein 1 and the level-2 cusp families. At 11 the
     constant 7 = a₁₁ of the weight-10 newform on Γ0(2) mod 11; at 13, Δ's family, as for j.
   - **Overpartitions done to 19 (12 Nov):** the constants are the Eisenstein 1, the level-1 cusp forms of weight ℓ − 3,
     and the level-2 newforms with w₂ = (2|ℓ) (Niwa: Γ0(4) ↔ Γ0(2)). Every prime 5–19 fits. Not yet tested out of
     sample: the next testable prime needs far longer series.
   - **Distinct parts done on lane ℓ (12 Nov):** the forms (η(2z)/η(z))·M_k⁰(Γ0(2)) give the data's laws at 5, 7, 11; at 7 the
     constant is confirmed to 7 digits. Left: 13–23 need towers near ℓ⁷–ℓ⁸; the next family at 7 and 11 (the rate).
10. More moonshine: other McKay–Thompson series, and Lehner-type laws wherever no ordinary cusp form exists. **run**
    - **T_2A and T_3C done on lane ℓ (12 Nov):** Lehner-type deaths where predicted (T_2A at 3, 5, 7; T_3C at 5, 7, 11),
      and the forms' laws everywhere the data test them. Left: levels 3, 4 and 6 (3A, 3B, 4A–4C, 6A–…) need engines of
      those levels; deep digits by Hida.
11. Mock theta functions, Hurwitz class numbers, and Mathieu moonshine's H(τ). **run**
    - **f(q) begun (12 Nov):** laws at 5, 7, 11, 13; the unit constants are the shadow's unit root (−3|ℓ) and (2|ℓ) from the
      level-6 class family of weight ℓ − 1 (its value at weight 0 is exactly 1, to 9 digits). Left: the slope-1 piece at
      5 and 7 (z ≡ 502 mod 5⁵, z ≡ −3 mod 7²; not level 6 at weight ℓ − 1); why the family is exactly 1 at weight 0; the
      further unit constants at 17, 19, 23 (towers to 10⁸ done via the Appell file; more levels need n near ℓ⁸); ω(q),
      fifth order, H(τ).
12. Non-congruence forms: Scholl's laws in our tower format, on a classical Atkin–Swinnerton-Dyer example. **run**
13. Higher rank: 1/χ₁₀ (Siegel; black-hole state counts), through its Fourier–Jacobi coefficients. **run**, after
    reading.
14. The modularity detector on many sequences (OEIS, physics counts), with controls. **run**

**The pole rule u^(w−2)**
15. Make the residue argument a proof. **outside** or careful work.
16. Double poles: 1/E₄², and E₄/j, where Zhang proved the constant (−3|p)·p. **Done in part (12 Nov):**
    - E₄/j reproduces Zhang at every prime 5–43.
    - For 1/E₄², the extended residue formula's ℓ·u⁻¹² is a root of the remainders' two-term law: all 6 digits at 7,
      one digit beyond chance at 13 and 19. The second root is a slope-1 holomorphic family.
    - WSL to 10⁸ (12 Nov): at 7, ℓ·u⁻¹² is confirmed to 4 digits (one slope-1 holomorphic form there, so the two-term fit is
      exact). At 13 there are three slope-1 forms near weight −8, too many for the fit; 19 can't be separated.
    - **Done (12 Nov):** the slope-1 families at weight −8 computed directly (Lauder's algorithm, `slope1_families.py`;
      checked against the classical forms at weight 148). With them the whole law has no free constant and holds to the
      whole lane at 7, 13 and 19 (`double_pole_full_law_….txt`). The earlier fit's family root at 7 was wrong in its 4th
      digit (3; the forms and the law give 6).
    - **Zhang's conjectures (12 Nov), `zhang_check.py`:** Theorem 1.3 and Conjectures 1.4 and 1.5 hold at j = 2 for 13
      primes, and 1.5 holds sharply. At supersingular primes X² − p² reaches 3l − r. **Next:** the sextic
      L₄·(X − p³)·(X − p³λ) at 19 and 23 (λ from Lauder at weight −2); a second curve (j = 4, supersingular at 31 and
      251); level 6 at 19. **WSL**
17. Poles at other CM points (discriminants −7, −8, −11, −19, …; class number 2, such as −15), and at points whose curve
    has no CM. **The no-CM case done (12 Nov):** 1/(E₄³ − 2Δ) adds u_C⁻¹⁴, reached one digit per level; the tower alone
    finds j = 2 (README "A pole at a point with no symmetry"). **Class number 2 done (12 Nov):** discriminant −15's two
    poles swap at 17 and 23 (period 2 seen in the tower), stay at 19 and 31; a pole on a CM lift is exact, a drifting one is
    one digit per level (README "Curves that repeat themselves"). **A 3-cycle seen (12 Nov):** −23 at 13 (settles three
    levels apart). **And in full at small primes (12 Nov):** discriminant −31 cycles at 5 and 7, where the cubic law reaches
    the whole lane. Left: −23 at 29, 31, 41 needs n beyond 10⁸ (or a leaner program); 37–43 with more levels. **WSL**
18. 1/E₄ at 37, where all four constants are ≡ 1: deeper data to separate them. **WSL**
