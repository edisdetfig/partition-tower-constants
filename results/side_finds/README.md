# Working notebook of the partition-tower study (11 to 15 Nov 2026)

This is the day-by-day record: each section names its scripts and kept outputs (`.py` beside `.txt`/`.json`).
Labels: counted (exact), measured (floating point), derived (reasoning).
In the surveys' four labels, counted and measured are both **computed**, and derived is **derived**: reasoning of the
project's own that no expert has checked. This notebook also writes "derived" beside a number calculated from counted
data. Where it says "proved" or "proven" of the project's own work, read "derived"; of a published theorem, read
"published". The wording of the dated entries below is left as written (note added 20 Nov 2026).
`FINDINGS.md` lists every finding with its status; `STATEMENT.md` states the main claims (A to H); `PROOF.md` and
`NOTE.md` hold the proof route and the draft proof; `NEXT.md` holds what is open.

Sections on an unrelated question (sums of three cubes) and the opening survey of another study were removed from
this public copy. A few later sections still refer to them by name ("Dig 2", "the side study").

## Dig 1: every 125th partition number from the 74th and the 124th (`partitions.py`, `partitions.txt`)

- **Known, and Ramanujan's own.** His unpublished manuscript (written about 1919 to 1920, published with proofs by
  Berndt and Ono, Séminaire Lotharingien de Combinatoire 42, 1999) states as eq. (21.3) that p(125n − 1),
  p(125n − 26) and p(125n − 51) are multiples of 125: the 124th, 99th and 74th. His proof is the one a reader rebuilt
  for the study. It uses his identity for p(25n + 24), and Jacobi's series for E(q)³ has nothing at those places on
  lane 5. MathWorld lists the three as ones that "seem also to hold". The study's finder found them blind.
- **The rule in one line** (derived here; counted below). Among the numbers N with 24N − 1 a multiple of 25, p(N) is
  a multiple of 125 exactly when (24N − 1)/25 is a square on lane 5 (0, 1 or 4). That is 3 of the 5 places, so there
  are 3 laws. For 49 it is the same with lane 7: (24N − 1)/7 a square there (0, 1, 2 or 4), which gives the 4 laws
  Ramanujan also wrote (eq. 7.8).
- **It reaches one level further, as forecast** (counted, p(n) to 400,000):
  - every 3125th from the 1849th, 2474th and 3099th is a multiple of 3125 (128 of 128 each);
  - for 7, every 2401st from the 929th, 1615th, 1958th and 2301st is a multiple of 343 (166 of 166 each).
  The 3125 laws are on MathWorld's list ("seem also to hold"). Of the 7 laws, only 2301 is there (Watson's); 929,
  1615 and 1958 I found in no source. Why these hold is the same reason as at 125, if the identity at that level has
  the same shape (derived, not checked in print).
- **Where the rule does not apply** (the controls, counted): at 625 it forecasts 3 laws and only Watson's 599 holds;
  at 343 it forecasts 4 and none holds. It applies only at the levels where the identity comes down to Jacobi's series.
- The full list found, power by power, for 5 (to 5⁵), 7 (to 7⁴) and 11 (to 11²), is in `partitions.txt`.

## Dig 1, deeper: the laws for 7 that no source states (`partitions_deeper.py`, `partitions_deeper.txt`)

p(n) for every n up to 2,400,000, worked out with whole numbers. All counts are counted.

- **The laws hold with far more cases.** Every 2401st partition number from the 929th, 1615th and 1958th is a multiple
  of 343: 1,000 of 1,000, 999 of 999 and 999 of 999.
- **They sort Ramanujan's one wrong guess.** He guessed that if 7³ divides 24N − 1, then 7³ divides p(N). Chowla
  found it false in 1934 at N = 243, and Watson (1938) proved only 7². Write 24N − 1 = 343q:
  - where q is a square on lane 7 (classes 929, 1615, 1958) the guess always holds: 2,998 of 2,998;
  - where q is 0 (class 2301) Watson's theorem makes it hold;
  - where q is not a square (classes 243, 586, 1272) it holds in only 153, 152 and 156 of 1,000. That is the 1 in 7
    that chance gives.
- **OEIS has the numbers but not the law.** A340957 (where the guess holds) and A340757 (where it fails), both by
  W. Bomfim, 2021, are rebuilt here from scratch, the first 39 and 40 terms exactly. Of the 39 terms shown in
  A340957, 38 come from the three laws and one, 25,282, is chance. Neither entry, nor MathWorld, states the law.
- **The reason, checked coefficient by coefficient.** At the level just below each set of laws, the generating
  function, divided by the power of ℓ that divides it, is on lane ℓ a constant times E(q^ℓ)^a · E(q)³. That is
  Jacobi's series, which is silent on whole classes. For 5: 3·E(q⁵)⁴E(q)³ at 25, 2·E(q⁵)⁴E(q)³ at 625,
  3·E(q⁵)⁴E(q)³ at 15,625. For 7: 1·E(q⁷)²E(q)³ at 7, 5·E(q⁷)²E(q)³ at 343, 4·E(q⁷)²E(q)³ at 16,807. This held for
  all n up to 3,000 (152 and 142 at the deepest levels). The constants step by −1 on lane 5 and by 5 on lane 7.
- **The pattern goes on, as forecast.** Every 78,125th from the 46,224th, 61,849th and 77,474th is a multiple of
  78,125 (30 each). Every 117,649th from the 45,519th, 79,133rd, 95,940th and 112,747th is a multiple of 2,401
  (20 each). 61,849 and 112,747 are Watson's classes; I found the others in no source either.
- **11 is different.** Only Atkin's classes appear (6, 116, 721), with no extra ones; 11's identity is not Jacobi's
  series.
- **Status.** The laws follow from the known theory (Watson 1938; Garvan's simpler proof for powers of 7, 1984), and a
  specialist could read them off those formulas. I found them stated nowhere. Reason: derived; laws: counted.

## Dig 1, again: the extra multiples (`partitions_fractal`, `partitions_bias`, `partitions_cells`, each .py and .txt)

- **The fractal, with the laws divided out.** Folsom, Kent and Ono (2012) proved that two levels apart the partition
  numbers repeat up to a constant. For 5, 7 and 11 their constant is 0, because the laws make every term divisible.
  Divide out what the laws force and a constant appears: p(625n + 599)/625 ≡ 4·p(25n + 24)/25 (mod 5) for every
  n ≤ 3,839, and the same 4 (= −1) at every level up to 5⁷. It is 5 (= −2) for 7 at every level up to 7⁶, and 9 (= −2)
  for 11 (levels 1 and 3, every n ≤ 1,802). Counted. This is the same thing as the constants 3, 2, 3 and 1, 5, 4 of
  "Dig 1, deeper". I found it stated nowhere.
- **Ramanujan's guesses for how often.** His manuscript says p(n) is a multiple of 5 "about 1/3" of the time and of 7
  "about 1/4" (from n ≤ 200). Up to 2,400,000 the shares are 0.3642 and 0.2728. Laws plus chance would give
  9/25 = 0.36 and 13/49 = 0.2653. His tally for 7 has slips: the true counts of the remainders 0 to 6 for n = 1 to 200
  are 51, 34, 23, 27, 23, 23, 19, not 50, 33, 22, 28, 23, 23, 21 (exact p(n)). Counted.
- **The extra multiples are of 5, 7 and 13, not 11.** Outside the laws, 20.52% of p(N) are multiples of 5 (chance:
  20%), 15.16% of 7 (14.29%) and 8.01% of 13 (7.69%), steady over 8 blocks of 300,000. For 2, 3, 11, 17, 19, 23 and
  31 it is chance. Counted. *This corrects what I said last time:* Berndt and Ono's 0.174 for 11 is exactly the laws
  plus chance, and their 0.080 for 13 is the extra.
- **Where they come from: Atkin's law of 1968, found again from the counts.** Write 24N − 1 = Q²·m, with Q a prime not
  dividing m, and let N* be the smaller number with 24N* − 1 = m. Radu proved in 2013 that a law p(An + B) ≡ 0 (mod ℓ)
  can only live where ℓ divides 24N − 1 or where −(24N − 1) is not a square on lane ℓ. Call the second "Radu's class".
  There, for 5, 7 and 13, p(N) ≡ μ·p(N*) (mod ℓ) with one μ per class, for every Q from 5 to 139. In the other half
  this never happens. For 11, 17, 19 and 23 it never happens in either half (one exception: 19 with Q = 103). Counted.
- **A formula for μ:** μ ≡ Q⁻²·[(3|Q)(1 + Q) − (−3|Q)(m|Q)] (mod ℓ), the same for 5, 7 and 13 ((a|Q) is the sign of a
  on lane Q). Fitted on Q < 140, then checked on every member for every prime Q up to 7,589: 0 wrong of 184,493 inside
  the fit and 0 wrong of 3,496 outside it. Fitted, then counted.
- **μ is 0 exactly when Q ≡ −2 (mod ℓ)** and (−m|Q) = −1. Then p(N) is always a multiple: 4,300 of 4,300 for 5,
  18,402 of 18,402 for 7 and 4,677 of 4,677 for 13. No other class is all multiples. As progressions:
  p(875n + B) ≡ 0 (mod 7) for B = 24, 274, 424, 549, 774, 799 (16,456 members); p(10985n + B) ≡ 0 (mod 5) for 12 values
  of B; p(17303n + B) ≡ 0 (mod 13) for 30 values, among them 237, Atkin's famous example. Counted. *Status:* this is
  Atkin's theorem (Proc. London Math. Soc. 18, 1968, Theorems 1 and 2, as quoted by Ahlgren, Allen and Tang 2022, eq.
  (1.7)). I did not see the formula for μ written down; Atkin's paper is behind a paywall.
- **Where Q ≡ −1 (mod ℓ)** the formula's first part is 0. By the Hecke rule p(N) is then a multiple whenever Q³
  divides 24N − 1 exactly. Derived, then counted: 179 of 179 for 5, 451 of 451 for 7, 1 of 1 for 13.
- **How much the law explains.** For 13: all of the extra outside the class 13 | 24N − 1 (7.683% left, against
  7.692%). For 7: about 90%. For 5: about 36%. The predicted share of all p(n) that are multiples, from the laws plus
  each prime's own odds, is 0.3615 for 5 (counted 0.3642), 0.2720 for 7 (0.2728) and 0.0800 for 13 (0.0801). Derived
  and counted.

## Still mysterious after this turn

- **A leftover the laws do not explain.** Take the numbers in Radu's class that no law touches and whose 24N − 1 has
  no square factor at all. Among them, 20.64% of p(N) are multiples of 5 (z = +15 against chance) and 14.46% of 7
  (z = +4.9). The share is flat across 177 (for 7: 186) residue classes and across sizes, and the extra comes evenly
  from the other remainders. In the class 13 | 24N − 1, p(13n + 6) ≡ 11 times the coefficient of qⁿ in ∏(1 − qᵏ)¹¹
  (mod 13), checked for every n ≤ 184,614. There 9.33% are multiples of 13. For 17, 19 and 23 in their classes the
  shares are 6.13%, 5.56% and 4.99% (chance: 5.88%, 5.26%, 4.35%). Counted. No law I know explains these. One
  resemblance only, not an explanation: for 5 the leftover is close to 1 − (1 − 1/5)(1 − 1/5³)(1 − 1/5⁵)⋯ = 0.2067,
  Delaunay's heuristic chance that 5 divides a Tate–Shafarevich group. For 7 the same formula gives 0.1454, two
  standard errors above the 0.1446 counted (both formula values derived).

## Dig 1, a third time: where the fractal constant leads (`partitions_tower.py`, `.txt`)

- **The constant settles into one ℓ-adic number.** Two levels apart, the laws divided out, one constant holds to more
  and more digits the deeper the levels: for 5 to 3, 6, 8, 11, 12 and 13 digits, all agreeing as far as they reach.
  Lowest digit first, the number is 4, 2, 3, 1, 0, 2, 0, 3, 0, 1, 1, 0, …. For 7 it is 5, 6, 1, 1, 2, …; for 11 it is
  9, 3, …; for 13 (no laws to divide) 6, 3, 1, …. Counted. For 13 this reproduces Folsom, Kent and Ono's published
  A₁₃(2, 4, 2) = 45 = 6 + 3·13, a check on the method. The number for 5 is not a simple fraction (the smallest fraction
  matching its 12 settled digits is −1437/592). I found these limits stated nowhere.
- **Each law, divided out, is again a Hecke eigenform.** Inside the laws of 5, 7 and 11, p(N)/ℓᵛ on lane ℓ obeys
  c(N) ≡ μ·c(N*) when 24N − 1 = Q²·(24N* − 1), one μ per sign, for every prime Q tested (25, 21 and 18 primes at the
  first level). In Folsom–Kent–Ono's terms, the "dimension 0" of these primes becomes 1 once the laws are divided out
  (my reading of the counts). The
  half-difference of the two μ is (−3|Q)·Q⁻² every time. That is the Hecke operator of weight −1/2 that Atkin wrote
  down, the weight of p(n) itself. Counted.
- **For 5 the multipliers are Atkin's own**: λ = (3|Q)·Q⁻²·(1 + Q) at levels 1, 2 and 3, as in Atkin's class. So the
  same primes Q ≡ −2 (mod 5) give zeros there too. For 7 and 11 the multipliers follow no formula of that kind, and none
  of 17 classic forms tried fits (the curves 11a to 64a, Ramanujan's Δ, and others). Their zeros fall at scattered
  primes: 13, 23, 41, 53, 61, 71, 89, … for 7, and 5, 67 for 11. *Update (next turn): the 17 forms were eta-products;
  Jameson's curves have conductors 42 and 726 and were not among them.*
- **Where μ is 0: one more factor of ℓ than the laws force.** Checked on every member, none missed:
  - p(10985n + B) ≡ 0 (mod 25) for 6 values of B (1,310 members);
  - p(60835n + B) ≡ 0 (mod 25) for 11 values (434);
  - p(54925n + B) ≡ 0 (mod 125) for 6 values (263);
  - p(15379n + B) ≡ 0 (mod 49) for 6 values (936);
  - p(85169n + B) ≡ 0 (mod 49) for 11 values (309);
  - p(107653n + B) ≡ 0 (mod 343) for 6 values (133);
  - and the shortest: **p(1375n + 149) and p(1375n + 424) ≡ 0 (mod 121)**, where Atkin's law forces only 11 (3,492
    members).

  *Status:* for 5 this is Lovejoy and Ono (J. reine angew. Math. 542, 2002); their example, p(10985N + 1514) ≡ 0
  (mod 25), is our first row. For 7 and 11 the analogue is Marie Jameson's "A refinement of Ramanujan's congruences
  modulo powers of 7 and 11", which describes the multipliers with two elliptic curves. Her paper would not download,
  so I cannot say which curves, or whether the 1375 law is printed there. *Update (next turn): it downloaded from the
  user's link; the curves are named in "Dig 1, a fourth time" below, and the 1375 law is her Corollary 1.3 with Q = 5.*
- **For 5, Atkin's class is level 0 of the tower.** For every N in Atkin's class, the partition number two levels deeper
  (24N′ − 1 = 25(24N − 1)), divided by 25, is 3·p(N) on lane 5: 38,400 pairs, no exception. In the other class the deeper
  value always has one more 5 (the square rule). So Ramanujan's laws, Watson's, the square rule, Atkin's class and the
  fractal are one structure for 5. For 7, 11 and 13 there is no such link (counted). I found this stated nowhere.
  *Correction and update (next turn): this bullet first said "25 levels deeper"; it is two levels (the factor 25 is
  5²). The three multiples of 125 behind the square rule, and the identity two levels down, are in Ramanujan's
  unpublished manuscript (Berndt–Ono, eqs. 21.3 and 21.4). The link itself is now one formula with a proof: see
  "Dig 1, a fourth time".*

## The leftover, dug into (`partitions_leftover.py`, `.txt`)

- **Where it shows** (24N − 1 with no square factor at all, so no Hecke rule reaches it; counted):
  - Atkin's class for 5 (20.64%, z = +15) and for 7 (14.46%, z = +4.9);
  - the class where 13 divides 24N − 1 (8.94%, z = +18);
  - inside the divided-out laws for 5 (z = +4.1 and +7.1), for 7 (+6.0), and for 11 (+7.5 and +5.8), although the
    plain values of 11 show nothing.
- **Where it does not show:** Atkin's class for 13 (z = +0.3, though it has Atkin's rule), every other class for every
  prime, and 17, 19 and 23 anywhere.
- **What it does not depend on.**
  - It does not depend on how many prime factors 24N − 1 has, or on which classes those primes fall in.
  - It holds in 177 (for 7: 186) residue classes, and at every size.
  - Deeper multiples get the same boost: of the multiples of 5, 20.09% hold 25 (chance: 20%). Derived: it behaves as
    if about 1 in 125 of these values were forced to be a multiple, with the rest random.
- **It is not the class numbers.** For N ≤ 3,000 I tested four fields (Q(√−M), Q(√−5M), Q(√−3M), Q(√−15M) with
  M = 24N − 1). The ratio p(N)/h on lane 5 takes all five values about equally, and 5 dividing h says nothing about 5
  dividing p(N). Counted. Still unexplained.

## Dig 1, a fourth time: the curves, the extra factors piling up, the one-line link, the constants

Files: `partitions_jameson`, `partitions_stack`, `partitions_multipliers`, `partitions_curve49`, `partitions_curve343`,
`partitions_tower_link`, `partitions_constants`, `partitions_simple` (each .py and .txt). Exact values of single p(N)
far beyond the recurrence come from the Hardy–Ramanujan–Rademacher formula (sympy). They agree with the recurrence on
201 sampled n, and every deep value is divisible by its law's power: p(146,230,061,849), with 425,986 digits, is a
multiple of 5¹⁶.

### What a zero multiplier means

The step rule: multiply 24N − 1 by Q² (Q a prime other than ℓ), and p(N) with its law divided out gets multiplied by a
fixed number on lane ℓ. That number depends only on Q and on which half of the remainders mod Q the number falls in. If
it is 0, every number reached by that step is divisible by one more ℓ than the law says. So one such prime gives a
whole family of laws, at every level of the tower. They are predictable, with no search.

### Jameson's paper (`partitions_jameson`)

- It downloaded from the link the user gave (12 pages). QUOTED: inside the laws, the multiplier at Q is
  Q^e·(3|Q)·(1 + Q − #E(F_Q) − sign term). The curve is E7: y² + xy + y = x³ + x² − 4x + 5 (conductor 42) for 7 and
  E11: y² + xy + y = x³ − 14x + 20 (conductor 726) for 11. Here #E(F_Q) is the number of points of the curve on the
  lanes mod Q. Lovejoy–Ono's theorem for 5 has the same shape, with 1 + Q in place of the curve. Her proof goes
  through the Shimura lift and Ramanujan's theta operator.
- Her theorem and theirs, on every value we have: 78,207 checks over levels 1–5 (5 and 7) and 1–2 (11), 0 wrong.
- Her table of primes is reproduced exactly from point counts (she also lists 11 itself for 11, which our count leaves
  out). Last turn's μ = 0 primes (13, 23, 41, 53, 61, 71, 89 for
  7; 5, 67 for 11) are exactly her ±1 primes. The 1375 law is her Corollary 1.3 with Q = 5.
- What share of primes gives a law one factor higher, counted to 30,000:
  - 7: 44.2% (limit 43.1%), with her counts reproduced to the 4th decimal;
  - 11: 26.6% (limit 27.3%);
  - 5: 50%, where only Q mod 5 matters.

  The limits are DERIVED by counting the matrices of GL2(F₇) and GL2(F₁₁) with each trace: 294, 287 and 287 of 2,016,
  and 1,210, 1,199 and 1,199 of 13,200, as she states. (One phrase in the .txt, "Q ≡ 4 zeroes both", should read:
  for Q ≡ 4 it is the multiplier when Q divides 24N − 1 exactly once that is 0, the Q³ rule.)
- Laws at seven primes we had not tested (Q = 103, 149, 31 for 7; 101, 109, 19 for 11; 113 for 5). Every member lies
  above 2.4 million and was checked with exact values: 10 of 10 each, N up to 161,111,732. Controls in the other class
  hit 1 to 5 times in 10. For the 11 cases each control hit is a member whose starting value was already a multiple
  (counted in a side run, not kept).

### The extra factors pile up (`partitions_stack`)

- The step rule holds beyond one digit (counted: one multiplier per class, every prime 5–37):
  - for 5, on lanes 25 and 125 at levels 1, 2 and 3, with the same multipliers at each level;
  - for 7 and 11, on lane ℓ² from level 2 on, but not at level 1, and not on lane ℓ³ at level 2.
- So zero steps stack where the rule holds deep enough. Each zero prime adds one factor (exact values; the smallest
  members; controls with every class flipped):

  | prime, level, zero primes | law | stacked | controls |
  |---|---|---|---|
  | 5, level 1, Q = 13, 23 | 5¹ | 5³ or more, 12 of 12 | 5¹, twice 5² |
  | 5, level 2, Q = 13, 23 | 5² | 5⁴ or more, 12 of 12 | 5², once 5³ |
  | 7, level 2, Q = 13, 23 | 7² | 7⁴, 12 of 12 | 7² |
  | 11, level 2, Q = 5, 67 | 11² | 11⁴ or more, 12 of 12 | 11², once 11³ |
  | 5, level 1, Q = 13, 23, 43 | 5¹ | 5⁴, 4 of 4 | 5¹, once 5² |

  The largest member here has 269,455 digits.
- **Contrast, where the rule fails on lane ℓ² (level 1 for 7 and 11):** only one extra factor is sure (7: 7² in 11 of
  12; 11: 11² to 11⁴, scattered). This is just what the step rule predicts.
- Lovejoy–Ono and Jameson prove one extra factor. Stacking I found in neither (nor in FKO).

### What the step rule is beyond one digit (`partitions_multipliers`, `partitions_curve49`)

- **The middle term is exactly Atkin's.** It is −(−3|Q)·Q⁻² on every lane tested: 5, 25, 125 for 5; 7, 49 for 7;
  11, 121 for 11. That is the weight −½ Hecke operator, not an approximation of it (counted).
- **For 5, the 1 + Q shape holds on lane 5 only.** The eigenvalue λ is (3|Q)·Q²·(1 + Q) on lane 5. On lanes 25 and 125
  no shape χ(Q)·Qᵃ·(1 + Q^c) fits. The object behind the 5-tower only looks like that simple form on one digit.
- **For 7, Jameson's curve holds on lane 49 too**, with Q¹⁹ where her lane-7 formula has Q¹:
  λ ≡ (3|Q)·Q¹⁹·(1 + Q − #E7(F_Q)) (mod 49).
  - At levels 2, 3 and 4 it holds for 74 of 74, 42 of 42 and 21 of 21 primes up to 397, 199 and 97, all with exact
    values (counted).
  - The exponent fits the weight. Both 1 and 19 make 2 + 2a ≡ −2 on the lane: mod 6 for lane 7 and mod 42 for lane 49.
    −2 is the weight that p(n)'s own series has under the Shimura lift (my reading).
  - For 11 the curve does not hold on lane 121 (counted).
- **What that predicts, and what was checked:** one prime alone gives two extra factors of 7 at level 2, where Watson
  gives 7² and Jameson's corollaries 7³.
  - Where 1 + Q − #E7 = 0 (Q = 31, 47): p(N) has 7⁴ or more in 24 of 24 members (24N − 1 = 49·Q³·m). The control
    Q = 41 (count −6) mostly stays at 7². (The .txt's "(any power, in fact)" assumes the formula holds on every lane;
    it is shown on lanes 7 and 49 only.)
  - Where it is ±48 or ±50, which is ±1 on lane 49 (Q = 607, 751, 1009, 1039, 1151, 1201): 7⁴ or more in 48 of 48
    members of the predicted class, and 7² in all 48 of the other class.
  - About 2% of primes fall in each of the three cases (counted to 30,000: 2.07%, 2.25%, 1.70%). The matrix share in
    GL2(Z/49) is 2.08% and 2.03% (derived).
  - I found this lane-49 form in no source.
- **How far the curve reaches (`partitions_curve343`).** From level 3 on, the step rule holds on lane 343 as well: one
  multiplier per class for every prime tried at levels 3 to 6. But no formula ±(3|Q)·Q^a·(1 + Q − #E7) with a = 0 to
  293 fits there. So E7 gets exactly two digits of the 7-tower right, and the third belongs to something else (counted).
  (That .txt's last line says "the same search finds a = 19 at levels 2, 3, 4". Precisely: the search found a = 19 at
  level 2, in partitions_multipliers, and partitions_curve49 confirmed that a at levels 2 to 4.)

### The tower link for 5, in one line (`partitions_tower_link`)

- **p(25N − 1) ≡ 25·N·(N − 3)·p(N) (mod 125), for every N:** 96,000 checked, 0 wrong. By class of N mod 5:
  - 0, 3, 4 give p(125n − 1), p(125n − 51), p(125n − 26) ≡ 0 (mod 125), Ramanujan's three (his eq. 21.3);
  - 1 and 2 give the ×3 link.

  With the laws divided out it holds at every depth: p(25N − 1)/5^(v+2) ≡ N(N − 3)·p(N)/5^v (mod 5), 96,000 checks over
  levels 0–7, 0 wrong. Below level 0, N(N − 3) ≡ −1, which is the first digit of the constant.
- **Why it holds at level 0 (derived; both halves also checked on all 96,000 values).** Two pieces:
  - (a) QUOTED, Ramanujan's eq. 21.4: p(25n + 24) ≡ 75·[qⁿ](q;q)²³ (mod 125).
  - (b) N(N − 3)·p(N) ≡ 3·[q^(N−1)](q;q)²³ (mod 5). With D = q·d/dq, E₂ = P, E₄ = Q, E₆ = R:
    - D(Σp(n)qⁿ) = (1 − P)/24 · Σp(n)qⁿ, and DP = (P² − Q)/12;
    - on lane 5, Q ≡ 1 and R ≡ P (his eq. 1.2), and Q³ − R² = 1728·q·(q;q)²⁴ (his eq. 1.1);
    - so (D² − 3D)(Σp(n)qⁿ) ≡ (1 − P²)/(q;q) ≡ (Q³ − R²)/(q;q) ≡ 3·q·(q;q)²³.

  In words: multiply each coefficient of the partition series by (24N − 1)² − 1 and, on lane 5, you get the series two
  levels down.
- **Only one digit at level 0:** on lane 25 the ratio takes all 5 possible values in every class of N mod 25 and mod 125.
- **Related primes: none.**
  - For 7 to 31, level 0 to 2 shows every possible ratio in every residue class.
  - From level 1 down, every prime has one constant: 5 (7), 9 (11), 6 (13), 6 (17), 10 (19), 5 (23), 3 (29), 30 (31).
  - For 7, two levels down is again (q;q)²³ (Jameson Cor. 2.2, checked on 48,979 values). But on lane 7 no factor fixed
    by N carries p(N) there.
  - The multipliers fit this (my reading). At level 0 they follow 1 + Q for 5, 7 and 13 alike (partitions_cells).
    Deeper, they follow 1 + Q for 5 but E7 for 7, so only 5's levels can be joined.
- **Predictions with exact values:**
  - the formula at N = 10⁸ to 10⁸ + 9 (p(25N − 1) of 55,691 digits): 10 of 10;
  - k steps down (24N_k − 1 = 25^k(24N − 1)), p(N_k)/5^(v+2k) ≡ (−1)^(k−1)·N(N − 3)·p(N)/5^v: 74 of 74, the largest
    with 121,782 digits.

### One number per prime (`partitions_constants`, `partitions_simple`)

- **Predictions from the 2.4-million data alone, against exact values:** every settled digit is right at every deeper
  level.

  | prime | digits predicted | levels | largest p(N) |
  |---|---|---|---|
  | 5 | 12 | 10–16 | 425,986 digits |
  | 7 | 5 | 8–13 | 291,829 digits |
  | 11 | 2 | 7–10 | 175,624 digits |
  | 13 | 2 | 6–9 | 77,655 digits |

- **The constants, lowest digit first (counted):**
  - 5 (33 digits): 4, 2, 3, 1, 0, 2, 0, 3, 0, 1, 1, 0, 1, 4, 0, 1, 3, 2, 3, 0, 3, 1, 0, 4, 3, 0, 4, 1, 2, 4, 2, 3, 0
  - 7 (14 digits): 5, 6, 1, 1, 2, 1, 3, 1, 0, 6, 5, 3, 1, 1
  - 11: 9, 3, 2, 0
  - 13: 6, 3, 1, 2, 5, 8
  - 17: 6, 12, 12, 7, 7
  - 19: 10, 2, 10, 8, 10
  - 23: 5, 20, 8, 4
  - 29: 3, 14
  - 31: 30, 6
- **They settle at fixed speeds.** Digits gained per two levels: 5 for 5, 3 for 7, 1 for 11, 2 for 13, and about 2 for
  17 to 23. My reading: this is the gap to the next-largest eigenvalue of the two-level step.
- **What is known:** the first digits are classical.
  - 5 for 7 and 4·5 ≡ 9 for 11: Watson's and Atkin's factors per level, in Jameson's Cor. 2.2 and 2.3.
  - 6 and 45 for 13: FKO.
  - −1 for 5: the one-line link above.

  The other digits I found nowhere. *Update ("a twelfth time" below):* every digit is the number of the tower form's
  ℓ-adic family at weight −2. Computed from modular forms alone, it agrees with every settled digit for 13 to 31.
- **Not simple numbers (`partitions_simple`).**
  - LLL (lattice reduction) looks for small whole-number relations. The constants' smallest relations of degree 1 to 4
    sit inside the range that 100 random ℓ-adic numbers give: ranks 5 to 99 of 100. The one low rank, 5, is 17 at
    degree 4, where 5 digits allow little; for 5, 7 and 13 the ranks are 21 to 99.
  - No unit root of a Hecke polynomial x² − a·x ± ℓ^k with |a| ≤ 2ℓ^(k/2) fits any of them. (*Explained below:* the
    constant belongs to weight −2, where no classical form lives.)
  - The section-4 run of `partitions_constants` stopped on an internal assertion in sympy's LLL. It was rerun on its
    own with an exact LLL written out in `partitions_simple.py`.
- **What kind of number (my reading; confirmed and made exact in "a twelfth time").** Ramsey's appendix to FKO puts the deep tower into ℓ-adic (overconvergent)
  modular forms of weight −½. In that frame each constant is the eigenvalue of one ℓ-adic eigenform under the two-level
  step. That fits all of the following:
  - one number per prime;
  - multipliers that match a classical object on one or two digits only (1 + Q for 5; E7 on lanes 7 and 49 for 7);
  - digits that look random.

### Still open after this turn

- What governs the multipliers beyond the classical object: 5 on lane 25, 7 on lane 343, 11 on lane 121. (*Answered for
  13, 17, 19 and 7 in "a thirteenth time" below, and for 5 and 11 in "a fourteenth time": the family's member at
  weight −2.*)
- Why the constants settle at 5, 3, 1 and 2 digits per two levels.
- The leftover (unchanged from "The leftover, dug into").

## Dig 1, a fifth time: why only 5's tower is joined to p(N), and a curve for 13

Files: `partitions_bolted`, `partitions_curves`, `partitions_curve13` (each .py and .txt). The question: every prime has
its own easily found constant, so why does it matter that only 5 has the level-0 link?

### What the level-0 link carries (`partitions_bolted`)

The constant joins a tower to itself: two levels down is a scaled copy. The level-0 link (5 only) joins the tower to
p(N) itself. Measured, for N with ℓ not dividing 24N − 1 (counted):

- **How much p(N) tells about the value two levels down.**
  - For 5: everything in Atkin's class. When p(N) is a multiple of 5, so is the deeper value, 100%. The information is
    1.16 bits (the shuffled noise floor is 0.0002); the other half of the N lies in the square-rule class, where the
    deeper value is always a multiple whatever p(N) is.
  - For 7 to 31: nothing measurable. The shares are equal within noise (7: 14.7% when p(N) is a multiple, 14.9% when
    not). The information equals the shuffled noise floor for every prime (7: 0.0037 against 0.0038 bits).
- **The step rules.** Level 0 (Atkin's class) and two levels down have the same multipliers at 14 of 14 primes Q for 5,
  at 0 of 13 for 7, and at 1 of 13 for 13. A link can only join two levels that follow the same rule.
- **What it means (my reading):**
  - With the link, everything about p(N) carries down 5's tower for free: its laws, its zeros, its leftover. And the
    tower's laws can be proved from p(N)'s own series (the one-line link of the fourth dig).
  - Without it, knowing p(N) tells nothing about the tower, and each tower has to be understood through its own rule.

### Each prime's rule (`partitions_curves`)

- **Every prime from 5 to 31 has a step rule.** At levels 1 and 2 there is one multiplier per class, with Atkin's middle
  term −(−3|Q)·Q⁻² every time: at 26 or 27 primes Q each at level 1, and 5 to 23 at level 2. Counted.
- **Only 5 follows the simple 1 + Q rule.** No shape ±(d|Q)·Q^a·(1 + Q^c) fits any other prime at either level.
- **The curves.** I tested every curve in J. Cremona's tables with conductor below 70,000 and no bad prime other than 2,
  3 and ℓ: 1,281 to 3,008 isogeny classes per prime (tables downloaded from github.com/JohnCremona/ecdata). The test
  was λ(Q) ≡ ±(d|Q)·Q^a·a_Q(E) (mod ℓ) at every prime Q with data, with the twists d ∈ {±1, ±2, ±3, ±6}. Counted:

  | ℓ | the curve | power Q^a | note |
  |---|---|---|---|
  | 5 | 48 classes, all of the kind whose counts mimic 1 + Q on lane 5 (50a, 50b, 75a, …) | — | the degenerate case |
  | 7 | 42a, Jameson's E7, and its twists | Q¹ | rediscovered |
  | 11 | 726d, Jameson's E11, and its twists | Q³ | rediscovered |
  | 13 | **1014f: y² + xy = x³ − 62x + 324**, conductor 2·3·13² | Q¹⁰ | also 1014b with Q⁴ (its twist by 13) |
  | 17, 19, 23, 29, 31 | none | — | |

  For 13 the fit holds at both levels: at 27 primes at level 1 and 13 at level 2. A chance fit would have odds near
  13⁻²⁷ against about 360,000 tries.
- **One recipe.** The three curves have conductors 42 = 2·3·7, 726 = 2·3·11² and 1014 = 2·3·13², the same twist (3|Q),
  and a power a with 2 + 2a ≡ −2 (mod ℓ − 1), the weight −2 of p(n)'s own series (my reading of the exponent).

### The curve behind 13 (`partitions_curve13`)

- **Far beyond the data, with exact values,** the formula λ = (3|Q)·Q¹⁰·a_Q holds at:
  - level 1: 106 of 106 primes up to 600;
  - level 2: 59 of 59 primes up to 300;
  - level 3: 32 of 32 primes up to 150.
- **What it predicts:** 13 has no Ramanujan law, but primes with a_Q ≡ 0 or ±1 on lane 13 give laws p(N) ≡ 0 (mod 13).
  - With a_Q ≡ ±1 (Q = 103, 191, 241, 251): one class of progressions modulo 13·Q³, 40 of 40 members.
  - With a_Q ≡ 0 (Q = 59, 73, 131): 24N − 1 = 13·Q³·m, 30 of 30 members.
  - Controls in the other class: 1 to 3 of 10, about chance.
- **A famous congruence is its first case.** QUOTED: p(59⁴·13n + 111247) ≡ 0 (mod 13). It is stated in Ono's
  "Distribution of the partition function modulo m" (Annals of Math. 151, 2000, eq. (2)), and a later thesis credits it
  to Atkin and O'Brien (as found by a web search; I read neither paper here). DERIVED: 24·111247 − 1 = 13·59³, and 59 is the first prime
  with a_Q = 0, where the curve has exactly 60 points. The same mechanism gives 58 progressions modulo 13·59⁴, one for
  each class of Y mod 59 in 24N − 1 = 13·59³·Y. Two members of each, with exact values: 116 of 116 are multiples of 13,
  and 111247 is among the 58.
- **How many primes give such laws:** 23.0% of the primes up to 30,000, against 23.1% from the matrices of GL2(F₁₃)
  with trace 0 or ±1 (derived).
- **Status:** Jameson's curves for 7 and 11 are hers. The curve for 13's tower I found in no source I read (Jameson,
  Lovejoy–Ono, Folsom–Kent–Ono, Ahlgren–Allen–Tang, Berndt–Ono's edition of Ramanujan, and the search result for Ono's
  example). *Correction (next dig): Ahlgren, Allen and Tang already name the rule of each prime from 13 on, as a
  classical newform of weight ℓ − 3 and level 6. What is new is only that, for 13, that form behaves like the curve 1014f
  on lane 13. See "Dig 1, a sixth time".*

### Still open after this dig

- What rules 17, 19, 23, 29 and 31. No curve over the rationals fits. Perhaps forms whose coefficients are not whole
  numbers (my reading, not tested). *Answered in the sixth dig: classical forms of level 6 and weight ℓ − 3; none of
  them behaves like a curve.*
- The three items of the fourth dig's list.

## Dig 1, a sixth time: the known rule behind each tower, and how new the 13 curve is (`partitions_newforms`)

**What the literature already has (QUOTED).** Ahlgren, Allen and Tang ("Congruences like Atkin's for the partition
function", arXiv 2112.09481) work with the newforms of weight ℓ − 3 and level 6 for congruences of this kind:
- S₂(6) = {0};
- S₄^new(6) is spanned by q − 2q² − 3q³ + 4q⁴ + 6q⁵ + …;
- S₁₀^new(6) by q − 16q² + 81q³ + 256q⁴ + 2694q⁵ + …;
- S₁₆^new(6) holds three newforms.

They prove that the mod-ℓ image of each of these contains SL₂(F_ℓ).

**Here (counted):**
- **The forms, computed from scratch.** I built the forms from products of E₂(z) − d·E₂(dz) and picked out the newforms
  by their U₂ and U₃ eigenvalues. This reproduces their printed coefficients exactly, for weights 4, 8, 10, 14, 16, 20,
  26 and 28 (one class at weight 28 is two-dimensional and was not split).
- **They are the towers' rules.** For every prime from 7 to 31, the measured eigenvalue is λ(Q) ≡ (3|Q)·a_Q(f) (mod ℓ)
  at every one of 26 or 27 primes Q, for one of these forms f:
  - weight ℓ − 3 for 13 to 31;
  - weight 8 for 11 (with an extra Q⁵);
  - for 7, a weight-16 form (the weight of Jameson's construction). On lane 7 the weights 4, 10, 16, 22, 28 count
    alike; the weight-4 form does not fit.
- **Which of the forms behave like a curve on lane ℓ.** I searched all 302,863 isogeny classes of Cremona's tables
  (conductor below 70,000; extra bad primes allowed and skipped). Every survivor of the filter at the primes 5 to 97
  was then checked at all primes up to 389:
  - 7: 42a (Jameson's E7) and its twists;
  - 11: 726d (E11) and its twists;
  - 13: 1014f, with power Q¹⁰ (and its twist 1014b, Q⁴), at all primes up to 389. This holds between two classical
    forms, with no partition numbers involved.
  - 17, 19, 23, 29, 31: no curve at all.
- **Why 5 is the odd one out, in these terms (derived).** Its rule would have weight 5 − 3 = 2, and S₂(6) = {0}. With no
  cusp form there, the rule can only be the 1 + Q kind, the same as p(N)'s own rule at level 0. That is why 5's tower is
  joined to p(N) (fifth dig).
- **A tie between primes (counted).** The weight-28 form that is 31's rule (a₂ = +8192, a₃ = −1594323) also reproduces
  7's rule when read on lane 7. The weight-16 form that is 7's rule sits beside 19's own weight-16 form (they differ in
  the signs of a₂ and a₃).

**How new.**
- The towers' rules as level-6 forms: Ahlgren–Allen–Tang's framework. Our data agree with it at every prime tested.
- The curves for 7 and 11: Jameson's.
- That the weight-10 form behaves like the curve 1014f on lane 13: I found it in no source I read. It is the 13-analogue
  of Jameson's paper. With her method it would turn into a proof by a finite check (a Sturm bound), which was not run
  here.
- That form's image contains SL₂(F₁₃) (their proposition). If the congruence holds, as all 73 primes up to 389 say,
  the curve's image is all of GL₂(F₁₃) (derived). Then the 23.1% share of law-giving primes for 13 is a real limit, not
  a guess. *Update: the congruence is proved in the seventh dig, so this holds.*

## Dig 1, a seventh time: the curve behind 13, proved the way Jameson proved hers (`partitions_curve13_proof`)

**The statement**, the 13-analogue of Jameson's Theorem 1.1. For every prime Q ≥ 5, Q ≠ 13, and every n ≥ 0:

  p(13Q²n + β) ≡ (3|Q)·Q⁴·[a_Q − ((−24n − 11)|Q)]·p(13n + 6) − Q⁹·p(13n/Q² + β′)   (mod 13),

where:
- β = (143Q² + 1)/24;
- a_Q = 1 + Q − #E(F_Q) for E = 1014b1: y² + xy + y = x³ − 10482x + 722308;
- the last term counts only when Q² divides 24n + 11.

Her powers of Q were Q¹ for 7 and Q³ for 11; this one is Q⁴. All three are Q^((ℓ−5)/2). 1014b1 is the twist by 13 of
the 1014f1 found in the fifth dig.

**The chain, every finite part of it checked (counted):**

| link | what is compared | range needed | result |
|---|---|---|---|
| (0) partitions to g = η(24z)¹¹ | p(13n + 6) ≡ 11·[qⁿ]∏(1 − qᵏ)¹¹ (mod 13), classical | every n with 13n + 6 ≤ 2.4 million | 184,615 values, 0 wrong |
| (1) g's eigenvalues to the level-6 form f | λ(Q) ≡ (3\|Q)·a_Q(f) (mod 13) | every prime ≤ 1,440, the Sturm bound for S₁₀(Γ0(864)), as in Jameson's proof | 225 primes, 0 wrong |
| (2) f to the curve | a_p(f) ≡ p⁴·a_p(E) (mod 13), with the 2-, 3- and 13-parts removed, both forms of weight 58 on Γ0(6084) modulo 13 | every prime ≤ 63,336, the Sturm bound 58·13104/12 | 6,342 primes, 0 wrong |
| (3) the statement itself | the formula above, on partition numbers | every n with N ≤ 2.4 million | 15,673 cases at 75 primes, 0 wrong |

For links (1) and (2), primes are enough because both sides are eigenforms with the same rule at prime powers. In (1)
the side not prime to 6 vanishes on both sides, and in (2) the weights fit: p⁹ = p⁴·p⁴·p.

**Standard theorems used and not re-proved:**
- Shimura's lift and its level (Niwa);
- Sturm's bound;
- the eigenform property of g (the spaces are one-dimensional for ℓ ≤ 31, as Ahlgren–Allen–Tang use);
- Serre's Θ operator on forms mod ℓ, through the Serre derivative and E₁₄ ≡ E₂ (mod 13);
- the modularity of E.

This is the same kind of proof Jameson gave for 7 and 11 (finite checks plus cited theorems). It is mine and has not
been checked by anyone else.

**What it gives (her Corollaries 1.2 and 1.3, for 13).** Write i = a_Q mod 13:
- i = ±1: p(13Q²(QN + r) + β) ≡ 0 (mod 13) for every r with ((−24r − 11)|Q) = i;
- i = 0: p(13Q²(Q²N + Qs + r) + β) ≡ 0 (mod 13) for r with Q | 24r + 11 and Q² ∤ 24(r + Qs) + 11.

The primes below 200 with i = 0 are 59, 73, 131 and 167; with i = +1 they are 97, 103, 109 and 191. The case Q = 59,
r = 2, s = 0 is exactly p(13·59⁴·N + 111247), the Atkin–O'Brien example in Ono's paper. These laws come from 23.1% of
all primes (derived from GL₂(F₁₃), whose full image now follows from Ahlgren–Allen–Tang's proposition and the
congruence).

**Not covered:** the deeper levels (p at 13ᵇ·Q²·n + …, b ≥ 2). There the same eigenvalues were counted (levels 2 and
3, at 59 and 32 primes, partitions_curve13.txt), but no Sturm check was run. *Update (eighth dig): no check is needed;
the theorem already covers every level.*

## Dig 1, an eighth time: every level of 13, and why no curve can exist for 17 to 31

Files: `partitions_curve13_levels`, `partitions_nocurve` (each .py and .txt).

### The 13 theorem at every level (`partitions_curve13_levels`)

- **Why no new Sturm check is needed (derived).**
  - 13 has no Ramanujan law, so nothing is divided out on the way down. A partition number at level b is also a
    level-1 number p(13n + 6), with 13^(b−1) dividing 24n + 11.
  - The seventh dig's theorem holds for every n, so it covers every level.
  - In level-b terms only the sign changes, by (13|Q)^(b−1).
  - This is unlike 7 and 11, where each level divides out a different power and Jameson needed a separate check for
    odd and even levels.
- **Shown directly (counted):**

  | level | cases | wrong | note |
  |---|---|---|---|
  | 2 | 1,089 (27 primes Q) | 0 | every case up to 2,400,000 |
  | 3 | 86 (12 primes Q) | 0 | every case up to 2,400,000 |
  | 4 to 7 | 12 each | 0 | exact values; the largest N at level 7 is 48,407,866,344 |

### No curve can exist for 17, 19, 23, 29, 31 (`partitions_nocurve`)

**The extrapolation.** The power of Q is not a free choice.

1. Suppose a curve E over Q carried the tower's rule up to any twist, that is, its mod-ℓ representation were that of
   the level-6 form f (weight k) times a character ψ.
2. Comparing determinants (ω for E, ω^(k−1) for f, ω the cyclotomic character) forces ψ² = ω^(2−k). So ψ = ω^((2−k)/2)
   times a quadratic character.
3. Twisting E removes the quadratic part. That is why the power came out as Q^((ℓ−5)/2) for 7, 11 and 13: 42a1, 726i1
   and 1014b1 have exactly the predicted traces at every prime up to 397.
4. So a curve would need trace u_Q = Q^(−(k−2)/2)·a_Q(f) (mod ℓ) at every prime Q not dividing 6ℓ.
5. There the representation is unramified, so E has good reduction (trace = 1 + Q − #E(F_Q), at most 2√Q in size, by
   Hasse) or multiplicative reduction (trace ±(Q + 1)). Additive reduction would make it ramified.

**The result (counted on the newforms, derived as above).** For each of these primes, at some Q the required trace
is none of the possible values:

| ℓ | impossible at Q = |
|---|---|
| 17 | 5 and 7 |
| 19 | 5 |
| 23 | 7, 11, 17 |
| 29 | 5, 7, 11, 31 |
| 31 | 7 |

For example, for 19 a curve would need a number of points mod 5 that Hasse's bound forbids.

So **no elliptic curve over Q, of any conductor, carries the rule of 17, 19, 23, 29 or 31**. The empty search below
70,000 was not bad luck. For 7, 11 and 13 the same test finds no obstruction up to 397, as it must. The .txt lists each
prime's "fingerprint": the only traces a curve could have at the first primes.

**What could carry the rule instead (my reading).** A weight-2 newform whose coefficients are not whole numbers, that
is, an abelian variety of dimension 2 or more. Its trace sizes are bounded only conjugate by conjugate, so Hasse's bound
does not block it. Finding it needs modular-symbol software (PARI or Sage).
- PARI has no Windows wheel.
- The user offered WSL (Ubuntu 22.04). I fetched the official Linux wheels of cypari2 2.2.0 and cysignals 1.12.5 into
  the scratchpad and verified their SHA-256 sums.
- Running them was blocked by the session's safety check, and was left for the user to decide. *Update: the user
  approved; see the ninth dig.*

## Dig 1, a ninth time: what takes the curve's place for 17 and up (`partitions_abelian`, PARI in WSL)

**How it ran.**
- PARI 2.15.4, through the official cypari2 2.2.0 wheel (SHA-256 checked), unpacked in the scratchpad and run in the
  user's WSL Ubuntu with the user's approval. The script runs in WSL; its header gives the command.
- Outputs: partitions_abelian_7_11_13_17.txt, _19.txt, _23.txt, _29.txt and _31.txt.
- Cross-check: PARI's level-6 newforms are ours. The weight-10 form has the same coefficients, and every tower's form
  was found with the a₂ and a₃ of the sixth dig.

**The search.**
- Every newform of weight 2 with trivial character at the levels ℓ, 2ℓ, 3ℓ, 6ℓ, ℓ², 2ℓ², 3ℓ² and 6ℓ², which is where
  Serre's conjecture places such a form here.
- Each was reduced at every prime above ℓ of its coefficient field and compared with the target of the eighth dig,
  u_Q = Q^(−(k−2)/2)·a_Q(f).
- Twists d = ±1, ±2, ±3, ±6, and those times ℓ, were included. The comparison ran at every prime Q ≤ 300 not dividing
  6ℓ, about 60 primes, so a chance match is out of the question.
- A form with a coefficient field of degree n corresponds to an abelian variety of dimension n over Q: a curve when
  n = 1, a surface when n = 2.

**Results (counted):**

| ℓ | the smallest match | its level | coefficient field | larger matches |
|---|---|---|---|---|
| 7 | the curve 42a (Jameson's) | 42 = 2·3·7 | Q | 294 (its twist by −7) |
| 11 | the curves 726i, 726d | 726 = 2·3·11² | Q | a pair with field Q(√5), at 726 |
| 13 | the curves 1014b, 1014f | 1014 = 2·3·13² | Q | a pair with field Q(ζ₇)⁺, at 1014 |
| 17 | **a surface** | 1734 = 2·3·17² | Q(√2) | Q(ζ₉)⁺ at 1734; degrees 3 and 4 at 867 |
| 19 | **a surface** | 2166 = 2·3·19² | Q(√5) | Q(ζ₉)⁺ at 2166 |
| 23 | **a surface** | 3174 = 2·3·23² | Q(√3) | Q(ζ₁₁)⁺ at 3174 (degree 5) |
| 29 | **a 3-dimensional one** (no surface) | 1682 = 2·29² and 5046 = 2·3·29² | Q(ζ₇)⁺ = Q(ζ₁₄)⁺ | Q(ζ₁₅)⁺ (degree 4) at both levels |
| 31 | **a 4-dimensional one** (predicted beforehand) | 5766 = 2·3·31² | Q(ζ₁₆)⁺ = y⁴ − 4y² + 2, and Q(ζ₁₅)⁺ | none smaller at any level |

Q(ζ_m)⁺ is the real part of the m-th cyclotomic field: Q(ζ₇)⁺ is y³ − y² − 2y + 1, Q(ζ₉)⁺ is y³ − 3y − 1, Q(ζ₁₁)⁺ is
y⁵ − y⁴ − 4y³ + 3y² + 3y − 1, and Q(ζ₁₅)⁺ is y⁴ − y³ − 4y² + 4y + 1. Q(√2), Q(√5) and Q(√3) are Q(ζ₈)⁺, Q(ζ₁₀)⁺ and Q(ζ₁₂)⁺.

**The rule behind the table (counted for 7 to 29; the reason is my reading).** Every matching form found has coefficient
field Q(ζ_m)⁺ with m = (ℓ − 1)/2 or m = (ℓ + 1)/2, of degree φ(m)/2:

| ℓ | m = (ℓ − 1)/2: field, dimension | m = (ℓ + 1)/2: field, dimension |
|---|---|---|
| 7 | 3: Q, 1 | 4: Q, 1 |
| 11 | 5: Q(√5), 2 | 6: Q, 1 |
| 13 | 6: Q, 1 | 7: Q(ζ₇)⁺, 3 |
| 17 | 8: Q(√2), 2 | 9: Q(ζ₉)⁺, 3 |
| 19 | 9: Q(ζ₉)⁺, 3 | 10: Q(√5), 2 |
| 23 | 11: Q(ζ₁₁)⁺, 5 | 12: Q(√3), 2 |
| 29 | 14: Q(ζ₇)⁺, 3 | 15: Q(ζ₁₅)⁺, 4 |
| 31 | 15: Q(ζ₁₅)⁺, 4 | 16: Q(ζ₁₆)⁺, 4 |

- **Why curves exist only for 7, 11 and 13.** φ(m)/2 = 1 only for m = 3, 4, 6. So one of (ℓ ± 1)/2 must be 3, 4 or 6,
  which happens for ℓ = 5, 7, 11 and 13 only. 5 is the degenerate case, whose rule is 1 + Q.
- **The prediction it made, then tested:** after 29, the rule said that for 31 the smallest carrier would have
  dimension 4 (m = 15 or 16), with nothing smaller. The 31 search, still running then, found exactly the two fields
  Q(ζ₁₅)⁺ and Q(ζ₁₆)⁺ at level 5766, both of degree 4, and nothing of degree 1, 2 or 3 at any of the eight levels.
- **What it predicts next (not run):** the smallest dimension is 3 for 37 (m = 18; 19 gives 9), 4 for 41 (m = 20; 21
  gives 6), 5 for 43 (m = 22; 21 gives 6), and 4 for 47 (m = 24; 23 gives 11).
- **My reading of the reason:** these forms have level ℓ² and a tame type at ℓ of order m. That type fixes the field of
  its traces, ζ_m + ζ_m⁻¹, inside the coefficient field.

**An example (17).**
- A curve would need 1 + 5 − #E(F₅) ≡ 8 (mod 17). Every curve over F₅ has between 2 and 10 points, so that is
  impossible (eighth dig).
- The surface's "count" at 5 is a₅ = 2 + √2, with its twin 2 − √2, both within Hasse's range ±2√5. On lane 17, √2 is 6
  (6² = 36 ≡ 2), and 2 + 6 = 8, as needed.
- For 19 likewise, a₅ = 2 + φ (φ the golden ratio), and on lane 19 φ is 5, giving the 7 that the tower needs.

**Other patterns (my reading).** The level is 2·3·ℓ² from 11 on (for 29 also 2·ℓ²), the same shape as Jameson's 726 and
our 1014.

**Status.**
- Counted at about 60 primes each. The three surfaces were also checked at every prime up to 1,000: 165 of 165 each
  (`partitions_surfaces`).
- Not proved. A Sturm check over the coefficient field would prove them, as for 13.
- I found none of them in a source I read.

## Dig 1, a tenth time: testing the rule at 37, 41, 43, 47 (`partitions_tower37`, `partitions_abelian_more`)

**The predictions, set down before the runs** (smallest dimension = min φ((ℓ ± 1)/2)/2):

| ℓ | predicted smallest dimension | field | the other side |
|---|---|---|---|
| 37 | 3 | Q(ζ₁₈)⁺ | m = 19 gives 9 |
| 41 | 4 | Q(ζ₂₀)⁺ | m = 21 gives 6 |
| 43 | 5 | Q(ζ₂₂)⁺ | m = 21 gives 6 |
| 47 | 4 | Q(ζ₂₄)⁺ | m = 23 gives 11 |

**From 37 on, the tower's rule splits (`partitions_tower37`, counted).**
- At level 1 there is no longer one multiplier per class. The ratios take almost every value in each class:
  - 37, 41 and 43: one multiplier at 0 of 11 primes Q;
  - 47: at 1 of 12.
- Up to 31 every prime Q had exactly one.
- This is Ahlgren–Allen–Tang's remark (QUOTED: "once ℓ ≥ 37 the spaces … are no longer one-dimensional"), seen in the
  partition numbers.
- In PARI the tower's class holds, for each of these primes, one pair of conjugate newforms over a real quadratic field
  (for 37: y² − y − 7187088385710). The class is the one with a₂ and a₃ of the signs (2|ℓ) and (3|ℓ), read off 13 to 31.
  - For 37, 41 and 43 the field splits at ℓ, so there are two separate rules on lane ℓ.
  - For 47 it stays whole (inert), so the rule lives in the field with 47² elements. The first 47 run looked only on lane
    47, found no target, and was stopped and its partial output removed. The script now compares in the right finite
    field; for split primes it gives the same results as before (checked again on 23).

**Results (counted; PARI, every prime Q ≤ 300 not dividing 6ℓ):**

| ℓ | predicted | found for the tower's own pair | field of the carrier |
|---|---|---|---|
| 37 | 3, Q(ζ₁₈)⁺ | **3**, at level 8214 = 6·37² | Q(ζ₁₈)⁺ exactly: **as predicted** |
| 41 | 4, Q(ζ₂₀)⁺ | **8**, at 10086 = 6·41² | contains Q(ζ₂₀)⁺, as a quadratic extension of it (branching at the prime 6481) |
| 43 | 5, Q(ζ₂₂)⁺ | **6**, at 11094 = 6·43² | Q(ζ₂₁)⁺ exactly (the m = (ℓ − 1)/2 side) |
| 47 | 4, Q(ζ₂₄)⁺ | **none** of dimension ≤ 6 (13254 = 6·47² alone took 79 minutes); then **8**, at 13254 (`partitions_abelian_more_47_8.txt`) | contains Q(ζ₂₄)⁺, as a quadratic extension of it (branching at the prime 73); 47 meets it in primes with 47² elements, as derived below |

- **Why 47 has no carrier that small (derived, then counted).** In Q(ζ_m)⁺ with m = (ℓ ± 1)/2 the prime ℓ is ±1
  modulo m, so ℓ splits there into primes with ℓ elements each. 47's two rules need the field with 47² elements: their
  numbers leave the field with 47 elements at 56 of the 59 primes Q from 5 to 300 (`partitions_split.txt`). So a carrier
  needs a prime above 47 with 47² elements. If its field contains Q(ζ₂₄)⁺ it has dimension at least 8; with Q(ζ₂₃)⁺, at
  least 22. The search up to 6 could not find one. **The search up to 8 then found it at exactly 8**, at 13254, with
  twists d = 1 and −47: field y⁸ − 24y⁶ + 170y⁴ − 408y² + 169 (discriminant 2¹⁶·3⁴·73²), whose quartic subfields include
  Q(ζ₂₄)⁺ = y⁴ − 4y² + 1 (PARI's nfsubfields). 47 splits fully in Q(ζ₂₄)⁺ and meets the octic in two primes with 47²
  elements (counted). The bound was set down before that run.

- **41, every level-6 form tried** (`partitions_abelian_more_41_all_8.txt`). The three single forms of the other sign
  classes do follow the rule: their carriers are Q(ζ₂₀)⁺ (4-dimensional) and Q(ζ₂₁)⁺ (6-dimensional), at 10086. So the
  predicted 4-dimensional objects exist for 41; they carry other forms' rules, not the tower's pair.
- **What survives, and what does not:**
  - Every carrier found from 7 to 43 has a field containing Q(ζ_m)⁺ with m = (ℓ ± 1)/2.
  - The prediction of the smallest dimension, min φ(m)/2, fails for the split towers of 41, 43 and 47. That prediction was
    drawn from 7 to 31, where the tower was a single form with whole-number coefficients. For 47 the failure is forced
    (derived above: its rules need the field with 47² elements).
  - From 37 on only one side, m = (ℓ − 1)/2, appeared for the tower's pair (37: 18; 41: 20; 43: 21); up to 31 both did.
- **For 37 the m = 19 side (dimension 9) did not appear,** although four 9-dimensional orbits exist at 8214.

## Dig 1, an eleventh time: what the split tells (`partitions_split`)

**The split, in plain words.** From 37 on, the tower is a blend of two rules, one from each twin form f, f′ (a conjugate
pair of level-6 newforms of weight ℓ − 3 over a real quadratic field). Each prime Q gives two numbers on lane ℓ,
(3|Q)·a_Q(f) and (3|Q)·a_Q(f′). Their sum t and product n are whole numbers; the two numbers alone are not.

**What it gives (counted, on exact partition numbers; PARI):**

- **A two-step law along Q.** Derived from the step rule, since two rules give T² − t·T + n = 0:
  c(Q⁴M) ≡ (t + εδ)·c(Q²M) − (Q⁻³ + δ² + t·ε·δ + n)·c(M) (mod ℓ), for M = 24N − 1 with ℓ | M and Q ∤ M, c(M) = p(N),
  ε = (M|Q), δ = −(−3|Q)·Q⁻².
  - Right on 48 of 48 chains for each of 37, 41, 43 and 47 (Q = 5, 7, 11, 13; 12 chains each).
  - One rule alone fits 9, 6 and 0 of the 48 (37, 41, 43).
- **Laws now need both numbers at once.** Primes Q below 10,000:

  | ℓ | p(N) ≡ 0 on a class, 24N − 1 = Q²·ℓ·n | both numbers 0, 24N − 1 = ℓ·Q³·m |
  |---|---|---|
  | 37 | 6599 ((n\|Q) = −1), 7541 (+1), 9547 (−1) | 2999 |
  | 41 | none | none |
  | 43 | 4049 (−1), 8941 (−1) | none |
  | 47 | none | none |

  - For 37 these are exactly Ahlgren–Allen–Tang's three, in their classes (QUOTED: "we have a congruence (1.6) when
    Q = 6599, 7541, and 9547"; 6599 with (n|6599) = −1, 7541 with (n|7541) = 1).
  - For each class law the two smallest members were computed exactly: all multiples of ℓ. One member of the other class
    each is not (14, 6, 35 for 37; 33 and 22 for 43). The largest was p(41,965,946,630).
  - 2999 for 37 was derived only: its smallest member, p(540,584,055,355), is past sympy's limit (149,832,547,102) and
    PARI is not trusted at that size. **Its cause was then checked** (`partitions_split_2999`): both numbers are 0 at 2999
    exactly when c(2999²·M) ≡ ε·δ·c(M) on lane 37 for every M; this held for 8 of 8 values of M (both classes; exact
    values from sympy, up to p(2,481,969,526)). So the law with 2999³ holds (derived from the step rule), and FLINT then
    computed the member itself: p(540,584,055,355) is a multiple of 37 (`partitions_flint_arbiter`, counted). In the
    literature's naming this is an Ono-type congruence (Yang, IMRN 2011, Remark 3.2, Theorem 7); Ahlgren–Allen–Tang's
    "Type II" is a different family, so "the second kind" below means Ono-type.
  - Rarer, as expected: up to 31 a law needed one number in place, about 2 primes in ℓ; now two, about 2 in ℓ² (derived
    estimate below 10,000: 1.8 class laws for 37, 1.5 for 41, 1.3 for 43, 1.1 for 47; laws of the second kind about half
    as many). None for 41 and 47 is within chance.
  - 4049 and 8941 for 43, and 2999 for 37, are not in the paper; I have not seen them listed elsewhere. *The literature
    check below found them nowhere either (the published lists stop at modulus 31), nor the class −1 for 9547, which
    the paper does not print.*
- **Along ℓ the constant splits too.** H_b(n) = p(ℓᵇn + δ_b), n = 0 to 9, levels 1, 3, 5:
  - No single constant links levels 1 and 3, or 3 and 5 (each of 37, 41, 43, 47).
  - Exactly one pair (T, N) gives H₅ ≡ T·H₃ − N·H₁, and it is the pair PREDICTED from the twins before the run:
    T = (3|ℓ)·(a_ℓ(f) + a_ℓ(f′)), N = a_ℓ(f)·a_ℓ(f′) (the Shimura correspondence with the twist (3|Q) found at every
    other prime).

    | ℓ | (T, N) from the partition numbers | predicted | the two numbers |
    |---|---|---|---|
    | 37 | (12, 29) | (12, 29) | 15 and 34 |
    | 41 | (13, 16) | (13, 16) | 23 and 31 |
    | 43 | (37, 14) | (37, 14) | 31 and 6 |
    | 47 | (29, 19) | (29, 19) | a conjugate pair in the field with 47² elements |

  - Neither number is 0, so neither rule fades out deeper down: at levels ≥ 3 the steps along Q still give several
    multipliers per class (Q = 5, 7, 11, 13).
- **The control, 13 to 31 (one form):** the constant from level 1 to level 3 equals (3|ℓ)·a_ℓ(f), six of six: 6, 6, 10,
  5, 3, 30 for 13, 17, 19, 23, 29, 31. These are the first digits of the constants in "One number per prime". FKO had 6
  for 13; the link to a_ℓ(f) I have not seen stated.

## Dig 1, a twelfth time: every digit of the constants, from modular forms (`partitions_hida`, `partitions_hida_digits`, `partitions_hida_check`)

**The idea (derived).** The partitions live in weight −½, whose partner in whole weights is −2. Suppose each constant is
(3|ℓ)·u(−2), where u(k) is the U_ℓ-number of the ℓ-adic family (Hida) of ordinary forms through f. Then:
- its first digit is (3|ℓ)·a_ℓ(f) (the control above: six of six);
- its first two digits are those of the family's member at weight (ℓ − 1)·ℓ − 2;
- all its digits follow from the members at weights ℓ − 3 + (ℓ − 1)·x, x = 0, 1, 2, …: Hida's theory makes u a power
  series with ℓ-adic integer coefficients, so its n-th difference in x is divisible by ℓⁿ, and Newton's series at x = −1
  (weight −2) is right to about as many digits as weights used.

**Results (counted). No partition number goes into the forms' digits.**
- **The second digit** (`partitions_hida`): one ordinary form of f's family at weight 154 (13), 270 (17) and 340 (19).
  Its (3|ℓ)·a_ℓ on lane ℓ² is 45, 210, 48: the constant each time. f itself gives 136, 57, 352, right on lane ℓ only, and
  no other form of those weights gives the constant. At weight 504 (23) PARI was hit by a signal (a backtrace, no
  message), hung, and was stopped; its file ends there. The next point covers 23.
- **All digits** (`partitions_hida_digits`, 11 weights per prime):

  | ℓ | from the forms (lowest first) | settled from the partition numbers | agree |
  |---|---|---|---|
  | 13 | 6, 3, 1, 2, 5, 8, 0, 9, 5, 9, 9 | 6, 3, 1, 2, 5, 8 | 6 of 6 |
  | 17 | 6, 12, 12, 7, 7, 5, 16, 3, 10, 14, 0 | 6, 12, 12, 7, 7 | 5 of 5 |
  | 19 | 10, 2, 10, 8, 10, 1, 16, 5, 7, 10, 6 | 10, 2, 10, 8, 10 | 5 of 5 |
  | 23 | 5, 20, 8, 4, 19, 3, 6, 6, 20, 20, 3 | 5, 20, 8, 4 | 4 of 4 |
  | 29 | 3, 14, 17, 22, 28, 6, 14, 19, 1, 14, 0 | 3, 14 | 2 of 2 |
  | 31 | 30, 6, 2, 28, 1, 8, 3, 15, 15, 13, 16 | 30, 6 | 2 of 2 |

  For every prime the n-th difference is divisible by exactly ℓⁿ (n = 0 to 10), and the series from 10 weights foretells
  the 11th to 10 digits.
- **The new digits, checked on partition numbers** (`partitions_hida_check`). The ratio of exact values p(δ_b) two levels
  apart agrees with the forms' digits in one more place per level: 4, 5, 6, 7, 8 places for the pairs (4, 6) to (8, 10).
  The pair (8, 10), with p(132,114,388,022) (404,903 digits), reads 6, 3, 1, 2, 5, 8, 0, 9. Digits 7 and 8 were foretold
  by the forms before any partition number had reached them.

**What this means.** The "one number per prime" is the U_ℓ-number, at weight −2, of the ordinary ℓ-adic family through
the tower's own form, twisted by (3|ℓ). Its digits look random because −2 is not a classical weight: no single classical
form carries the number, only the family does. This fits Ramsey's frame (overconvergent forms of weight −½) and makes it
exact and computable. I have not seen the constants computed this way, nor any digit beyond 13's first two.
- **An older find, read again (my reading, derived):** E7's law on lane 49 for 7, λ ≡ (3|Q)·Q¹⁹·a_Q(E7). The twist by Q¹⁹
  moves weight 2 to 40, and 40 ≡ −2 modulo 6·7, so E7 twisted agrees with the weight −2 point on lane 49. On lane 343 it
  cannot: 40 ≢ −2 modulo 6·49. That is where it was seen to fail.

## Dig 1, a thirteenth time: the same weight −2 member in every direction (`partitions_hida_steps`, `partitions_twins`, `partitions_hida_seven`)

**The steps along other primes Q, beyond lane ℓ (`partitions_hida_steps`).** If the deep tower is the family's member at
weight −2, its steps along every prime Q follow that member's a_Q, not only on lane ℓ. A_Q = a_Q of the family at weight −2
comes from Newton's series over 11 weights (forms only); the predicted multiplier is μ_ε = (3|Q)·A_Q − ε·(−3|Q)·Q⁻².
- PREDICTED before the run: at level b the measured ratio c(Q²M)/c(M) shares about b digits with μ_ε; f itself, 1 digit.
- Counted on exact partition numbers (Q = 5, 7, 11, both classes, the first suitable N at each level): at level b the ratio
  shares **at least b digits** with μ_ε in all 114 comparisons (13: levels 1 to 7; 17 and 19: levels 1 to 6), a few
  even more. With f itself: 1 digit in all 114.
- So the open question "what governs the multipliers beyond the classical object" has an answer for 13, 17, 19: the
  weight −2 member of the tower form's family. One object runs both directions, along ℓ (the constant) and along Q.

**The twins' two constants (`partitions_twins`).** From 37 on, each twin form is ordinary on lane ℓ and has its own family,
so each prime has two constants. From the forms (9 weights, digits lowest first):

| ℓ | twin 1 | twin 2 | on lane ℓ: sum, product |
|---|---|---|---|
| 37 | 15, 4, 4, 3, 1, 10, 35, 5, 36 | 34, 30, 35, 25, 2, 25, 28, 28, 1 | 12, 29 |
| 41 | 23, 14, 7, 38, 12, 26, 7, 1, 6 | 31, 11, 10, 38, 14, 13, 39, 40, 7 | 13, 16 |
| 43 | 31, 13, 8, 5, 40 | 6, 20, 10, 25, 3 | 37, 14 |

- The sums and products on lane ℓ are the (T, N) found in the partition numbers ("an eleventh time").
- On exact partition numbers H_b(0) = p(δ_b) the two-number law H_(b+4) ≡ T·H_(b+2) − N·H_b holds, with T and N from the
  forms to all their digits: on lane ℓ at levels 1, 3, 5 and 2, 4, 6; on lane **ℓ²** at levels 3, 5, 7. That is what was
  PREDICTED before the run (37 and 41; the deepest value used, p(56,803,329,882), has 265,496 digits).
- One constant alone fails: H_(b+2) − αᵢ·H_b is not divisible by ℓ for 41. For 37, twin 1 alone fits on lane 37 at n = 0
  (twin 2's part happens to vanish there on lane 37) but never on lane 37².
- 43 first crashed inside PARI at a high weight (a segmentation fault, then the 8 GB PARI stack ran out); it was rerun
  with 5 weights, which is enough for the lanes ℓ and ℓ² (`partitions_twins_43.txt`). That run reported the law failing
  even on lane 43 at levels 3, 5, 7, which the theory forbids. **The cause was a wrong partition number from PARI:**
  numbpart(56,628,877,314) is 23 on lane 43 (computed twice), but sympy's partition() gives 36, the value the law needs;
  the smaller values agree in both programs. With sympy's values (`partitions_twins_43_sympy`) 43 behaves like 37 and 41:
  lane 43 at levels 1, 3, 5 and 2, 4, 6, lane **43²** at levels 3, 5, 7.
- **What the PARI error touches.** Every other PARI partition number in these digs went into a check that agreed with an
  independent prediction in 2 to 12 digits; a wrong value would fake that only by a chance of 1 in ℓ² or less. Only a
  failure can come from such an error, and this was the one failure. The 13 check was redone with sympy
  (`partitions_hida_check_sympy`): the same ratio digits at every pair of levels, 4 to 8 places agreeing with the forms.
- **Settled with FLINT (`partitions_flint_arbiter`).** python-flint 0.9.0 (added to `requirements.txt`) computes p(n) with
  Arb's ball arithmetic, which FLINT documents as exact. p(56,628,877,314) is 36 on lane 43 and 63,332 on lane 43³:
  **sympy was right, PARI 2.15.4 was wrong** (counted). PARI's bug tracker has no such report (the literature check:
  its one wrong-value numbpart bug, #557, was off by one and fixed in 2007; this one is off by at least 357). With
  FLINT's values: the 15 law members and controls of `partitions_split` agree with PARI's at 15 of 15; the twins' law
  gives ℓ¹, ℓ¹, ℓ² for 37, 41 and 43; 13's deep check gives 4 to 8 digits as before. Partition numbers come from FLINT
  from here on.

**7, a prime with laws (`partitions_hida_seven`).** Watson's laws add one factor of 7 per two levels, so the two-level step
is 7·α₇ and the family must have slope 1 (its U_7-number divisible by 7 exactly once). By my reading of the Gouvêa–Mazur
picture such a family is one analytic piece only over weights that agree modulo 6·7 = 42; the weight that agrees with −2
there is 40, the same weight E7's lane-49 law pointed to (2 + 2·19).
- At every weight from 16 to 292 the forms of 7's tower system include exactly one of slope 1 (counted). The others
  have slopes 7, 14, 19, 24, 31, one more appearing every weight step of 42.
- Its (3|7)·a_7/7 on lane 7 climbs by one per weight step of 6: 1, 2, 3, 4, 5, 6 at weights 16 to 46. At weight 40 it is
  5, the first digit of α₇.
- From the weights 40, 82, …, 292 (forms only), Newton's series at weight −2 gives (3|7)·u(−2)/7 =
  **5, 6, 1, 1, 2, 1, 3, 1, 0, 6, 5, 3, 1**, 2. The partition numbers gave 5, 6, 1, 1, 2, 1, 3, 1, 0, 6, 5, 3, 1, 1:
  **13 of 14 digits agree** (counted). The differences shrink faster than for 13 (divisible by 7 to the powers 0, 1, 4,
  5, 7, 9, 11), so the series reaches about 13 digits; the 14th is past it. The first 7 were PREDICTED before the run
  (`partitions_hida_seven.txt`); the longer comparison is a second run (`partitions_hida_seven_14.txt`), measured.
- In the trials PARI crashed above weight 292 (a segmentation fault), and roots computed to 25 digits lost the members of
  weight 334 to the huge powers of 7 in their coefficients; the kept runs use roots to 400 digits and stop at 292.
- **7's steps along Q (`partitions_hida_seven_steps`).** The same slope-1 member's a_Q at weight −2 (forms only; its
  differences shrink like 7^(2n)) against exact partition numbers with Watson's factor divided out, Q = 5, 11, 13, both
  classes. PREDICTED before the run: more digits at deeper levels, past lane 49, while E7's rule stays at 2. Counted:
  digits shared at levels 1 to 9: 1, 2, 4, 5, 7, 8, 10, 11, 12 (the last at the 12-digit cap), every Q and class, a few
  one more; E7's rule λ ≡ (3|Q)·Q¹⁹·a_Q(E7): 2 from level 2 on, never 3. So the old question "7 on lane 343" is answered:
  E7 is what the weight −2 member looks like on lane 49, and the member itself holds on every lane.
- So the picture holds for a prime with laws too: the slope of the family is the law, and the constant is the family's
  number at weight −2.

## Dig 1, a fourteenth time: 5 and 11, the other primes with laws (`partitions_hida_laws`, `partitions_hida_law_steps`)

**The idea.** Watson's laws for 5 and Atkin's for 11 add two factors of ℓ per two levels, so each tower should be the weight
−2 member of a family of **slope 2** (its U_ℓ-number divisible by ℓ exactly twice), with α_ℓ = (3|ℓ)·u(−2)/ℓ². The tower's
system on lane ℓ: for 5 the "1 + Q" rule, a_Q ≡ Q²(1 + Q); for 11 the weight-28 form matched in `partitions_newforms`
(a₂ = −8192, a₃ = +3¹³), which sits in the sign class the rule for 13 to 31 gives.

**The survey (counted; `partitions_hida_laws.txt`, part 1).** Members of the system at the weights ≡ −2 modulo ℓ − 1:
- 5: slope-2 members only at the weights ≡ 18, 2, 6 (mod 20), with first digit 4, 2, 4; slope 3 at ≡ 10, 14 (mod 20).
  −2 ≡ 18, and α₅ begins with 4.
- 11: slope-2 members at 28 to 108 with first digits 9, 5, 10, 2, 3, 2, 10, 5, 9 (a mirror around 68), none at 118
  and 128 (slope 3), and 138, 148 repeat 28, 38. The pattern repeats every 110. −2 ≡ 108, and α₁₁ begins with 9.
- So the weights that agree with −2 repeat every (ℓ − 1)·ℓ: 20 for 5, 42 for 7, 110 for 11.

**The digits (counted, forms only; part 2).**

| ℓ | weights | from the forms (lowest first) | from the partition numbers | agree |
|---|---|---|---|---|
| 5 | 18, 38, …, 258 | 4, 2, 3, 1, 0, 2, 0, 3, 0, 1, 1, 0, 1, 4, 0, 1, 3, 2, 3, 4, … | 4, 2, 3, 1, 0, 2, 0, 3, 0, 1, 1, 0, 1, 4, 0, 1, 3, 2, 3, 0, … | **19** |
| 11 | 108, 218, 328 | 9, 3, 2, 0, 5, 8 | 9, 3, 2, 0 | **4 of 4** |

- 5: exactly one slope-2 member at each of the 13 weights; the differences are divisible by 5 to the powers 0, 1, 2, 4, 6,
  10, 10, 12, …, 18, and the series foretold the last weight to 18 digits. The 20th digit is past the series' reach.
- 11: three weights are what PARI managed here, so the series is good to about 3 digits; the 4th agreeing may be luck.

**The steps along Q (counted; `partitions_hida_law_steps_5.txt`, `_11.txt`; partition numbers from sympy).** The same
member's a_Q at weight −2 against the measured multipliers, the laws' factor ℓᵇ divided out, Q = 7, 11, 13 for 5:
- 5: digits shared at levels 1 to 5: 4, 6, 9, 11, 12 (the 12-digit cap), every Q and class, then 12 at every deeper
  level. The lane-5 rule (3|Q)·Q²·(1 + Q): 1 digit, always. So "on lanes 25 and 125 no shape fits" has its answer: the
  object is the weight −2 member, which no simple shape is.
- 11 (Q = 5, 7, 13): digits shared at levels 1 to 7: 1 to 2, 2, 2 to 3, 3, 3 to 4, 4, 4 to 5, growing about one per two
  levels, as 11's constant settles; the lane-11 rule: 1, always. So the rule that failed on lane 121 is again the weight
  −2 member, here followed to lane 11⁵.

**Where this leaves the picture.** Every prime from 5 to 43 now fits one statement: the partition tower at ℓ is the weight
−2 member of an ℓ-adic family of modular forms (two twin families from 37 on), whose slope is the law (0 for 13 to 43, 1 for
7, 2 for 5 and 11). Its U_ℓ-number gives the constant and its a_Q give the steps, on every lane checked.

## Dig 1, a fifteenth time: how deep the elliptic curves go (`partitions_curves_depth`)

**The question.** On lane ℓ each small tower's steps follow a weight-2 object twisted by a power of Q: the Eisenstein shape
1 + Q for 5 (with Q²), Jameson's E7 = 42a1 for 7 (Q¹), her E11 = 726d1 for 11 (Q³), our 1014b1 for 13 (Q⁴). A weight-2
object twisted by Q^j acts like weight 2 + 2j, so to agree with the weight −2 member it needs 2 + 2j ≡ −2: modulo ℓ − 1 on
lane ℓ (which gives exactly the exponents 2, 1, 3, 4), and modulo (ℓ − 1)·ℓ on lane ℓ² (8, 19, 53, 76). Here each family's
a_Q at weight −2 (forms only) is compared with every sign, quadratic twist (±1, ±2, ±3, ±6) and exponent j modulo (ℓ − 1)·ℓ²,
at the 22 primes Q from 5 to 100. PREDICTED before the run: E7 2 digits, the Eisenstein shape and E11 1, and for 1014b1 2
digits with j = 76 "if the weight rule is the whole story".

| ℓ | the weight-2 object | digits right at all 22 primes | lane ℓ² with the weight rule's exponent |
|---|---|---|---|
| 5 | 1 + Q | 1 | Q⁸: 2 of 22 primes (chance about 4) |
| 7 | E7 = 42a1 | **2**, with Q¹⁹ (also Q⁶¹, Q¹⁰³, …: j ≡ 19 mod 42) | Q¹⁹: **22 of 22** |
| 11 | E11 = 726d1 | 1 | Q⁵³: 1 of 22 (chance about 2) |
| 13 | 1014b1 | 1 | Q⁷⁶: 3 of 22 (chance about 2) |

- **The prediction for 13 failed:** its curve gives the first digit only. Only E7 reaches lane 49, and nothing reaches a
  third digit (counted). The weight rule names the only exponent that could work; it does not make a curve go deeper.
- **Why E7 (my reading):** its level holds ℓ once (42 = 6·7: multiplicative at 7, itself 7-ordinary), while 726 = 6·11² and
  1014 = 6·13² hold ℓ squared: those curves are twists by a character of conductor ℓ, which imitates Q^j on lane ℓ only.
- So the elliptic curves (and 5's Eisenstein shape) are the families' shadows on the first lane; E7 alone is a shadow on
  two. The families' own numbers carry every lane.

## Dig 1, a sixteenth time: deeper with FLINT, and three forms from 61 on (`partitions_flint_deeper`, `partitions_hida_digits_x12`, `partitions_three_forms`, `partitions_three_check`)

**FLINT** computes p(n) exactly, and on 8 threads about 5 times faster than on one (p(10¹⁴) in 7.7 s against 38.2 s;
measured). The literature check is in `reports/Partition tower literature check.md`; its column of digits confirmed by
partition numbers predates this section.

**Deeper checks (counted; `partitions_flint_deeper.txt`; every PREDICTION set down in the script before the run).**
- The laws: 10 of 10 members are multiples of ℓ for 37 at Q = 6599, 7541, 9547 and 2999 (the last up to
  p(9,522,596,052,022)) and for 43 at 4049 and 8941.
- The twins: H_(b+4) − T·H_(b+2) + N·H_b is divisible by ℓ, ℓ, ℓ², ℓ², ℓ³ at the levels (1,3,5) to (5,7,9), for each of 37,
  41 and 43: lane ℓ³ at levels 5, 7, 9, as predicted.
- The constants: the forms' digits (13 weights for 13 to 23, `partitions_hida_digits_x12.txt`; 11 for 29 and 31) against
  the ratio of exact values two levels apart, as deep as 10¹⁷:

  | ℓ | digits confirmed by partition numbers (before → now) | deepest value used |
  |---|---|---|
  | 13 | 8 → **13** (all the forms give: 6, 3, 1, 2, 5, 8, 0, 9, 5, 9, 9, 6, 8) | p(23,460,200,964,791,597), about 170 million digits |
  | 17 | 5 → **11** | p(2,888,835,259,597,565) |
  | 19 | 5 → **11** | p(8,761,038,221,303,554) |
  | 23 | 4 → **11** | p(21,001,515,080,686,141) |
  | 29 | 2 → **5** (3, 14, 17, 22, 28: the 17 is the digit in Boylan–Webb's 14706) | p(9,658,736,897,850,448) |
  | 31 | 2 → **5** (30, 6, 2, 28, 1) | p(17,997,671,134,953,422) |

  56 digits in all, each agreeing; at every pair the agreement grew by about one digit per level, as before. The 12th
  and 13th digits of 13 came from the longer forms run only, so they were fresh predictions.
- For 23 the levels come in equal twos (24·268181 − 1 = 23⁵ exactly, and so on), so each even pair repeats the odd one.

**Three forms from 61 on (counted).**
- The census (`partitions_three_forms.txt`): the tower's class of S_{ℓ−3}^new(Γ0(6)) has dimension 2 for 53 and 59, 3 for
  61, 67, 71, 73, 79, 83 and 4 for 89, each time one Galois orbit: exactly ⌊(ℓ + 12)/24⌋, the rank Boylan–Webb computed,
  as PREDICTED.
- 61, against the cubic polynomials of its class's numbers (`partitions_three_check.txt`, FLINT): along ℓ the four-term
  law (H₇ from H₅, H₃, H₁) holds for 10 of 10 values of n, and no two-term law or single constant fits; along Q the
  three-step law (c(M), c(Q²M), c(Q⁴M), c(Q⁶M)) holds on 48 of 48 chains (Q = 5, 7, 11, 13), and no law of degree 2 fits
  the 12 chains of any Q. The control 53 (two forms): 10 of 10 and 48 of 48.
- So the number of forms blended, and of constants, grows with the prime: none for 5, 7, 11 (the laws), one for 13 to
  31, two for 37 to 59, three for 61 to 83, four from 89.

## Dig 1, a seventeenth time: modular forms on FLINT, and the constants of 47 and 61 (`forms6_flint`, `partitions_flint_family`, `partitions_hida_more_check`)

**No more PARI (the user's instruction, 12 Nov 2026).** FLINT has no modular-forms routines, so the forms side was rebuilt
on it (`forms6_flint.py`):
- **The basis.** For even k ≥ 4, M_k(Γ0(6)) (dimension k + 1) has the basis G_k·t^(k−m), m = 0..k, from two eta-quotients:
  G_k = (η(z)η(6z)⁶/(η(2z)²η(3z)³))^k (all k zeros at the cusp ∞) and t = η(z)⁵η(3z)/(η(2z)η(6z)⁵) (a simple pole at ∞, a
  simple zero at 0); their orders at the four cusps come from the standard formula (derived). A form is fixed by its
  coefficients q⁰ … q^k (Sturm's bound for Γ0(6) is k).
- **Hecke operators** on q-expansions modulo ℓ^N, as matrices in those coordinates.
- **The ordinary part first.** The forms whose number at ℓ is a unit make a piece whose size does not grow with the weight
  (Hida): 11 for 13, 35 for 37 at every weight tried. Its projector is a polynomial in T_ℓ, lifted from mod ℓ by Newton's
  e ← 3e² − 2e³; inside it the tower's forms are picked out by their Hecke numbers on lane ℓ and lifted by Newton's method.
  (Picking them out of the whole space first failed: other forms share the numbers on lane ℓ by chance, about one in ℓ.)
- ~~For ℓ ≥ 7 no lower-level or Eisenstein form can share the tower's numbers (a level-lowering congruence would need
  a_p ≡ ±(p + 1)·p^((k−2)/2), impossible at p = 2, and at p = 3 unless ℓ = 5; derived), so the Hecke numbers alone suffice.~~
  **WRONG** (found the same day; the eighteenth time below). Exact equality is impossible (Deligne's bound), but agreement
  on lane ℓ is not: it happened for 17, 29, 37 and 499.

**Validation (counted).** 13: all 13 digits known from partition numbers, then 4 more: 6, 3, 1, 2, 5, 8, 0, 9, 5, 9, 9, 6,
8, 2, 8, 5, 3 (17 weights, 2 s). 37: both twins exactly as PARI had them for 9 digits, now 12: 34, 30, 35, 25, 2, 25, 28,
28, 1, 14, 0, 0 and 15, 4, 4, 3, 1, 10, 35, 5, 36, 21, 9, 32 (17 s). Every n-th difference divisible by exactly ℓⁿ.

**47 and 61 (counted; `partitions_flint_family_47.txt`, `_61.txt`).**
- **47** (the twins' field leaves 47 whole): two constants **A ± B·√5**, √5 not existing on lane 47:
  A = 38, 4, 33, 43, 39, 34, 9, 21, 12, 37, 23 and B = 12, 22, 14, 7, 15, 43, 4, 38, 21, 1, 21 (11 weights, to weight 504).
  Their sum and product on lane 47 are 29 and 19, the numbers the partition numbers gave (partitions_split.txt).
- **61** (61 breaks into one piece of degree 1 and one of degree 2 in the cubic field): one constant on lane 61,
  32, 40, 3, 28, 17, 11, 0, 28, 10, 17, and a pair **A ± B·√2** with A = 60, 24, 40, 28, 20, 18, 11, 21, 6, 16 and
  B = 49, 14, 33, 55, 8, 59, 49, 27, 27, 33 (10 weights, to weight 598). On lane 61 the three give the cubic 34, 15, 31, 1,
  the one the partition numbers gave.
- A trial on PARI before the instruction had agreed in its first 6 (47) and 5 (61) digits (its B was −B, the same pair).
- **On exact partition numbers (FLINT; `partitions_hida_more_check.txt`), the laws built from these constants hold as
  PREDICTED:** 47: lane 47 at levels (1,3,5), (2,4,6), lane 47² at (3,5,7), (4,6,8), lane **47³** at (5,7,9) (n = 0; for
  n = 1 at least as deep); 61: lane 61 at (1,3,5,7), (2,4,6,8), lane **61²** at (3,5,7,9).

## Dig 1, an eighteenth time: forms of lower level, the exception 607, and every constant to 20 or 30 digits (`forms6_flint`, `partitions_flint_family`, `partitions_flint_family_check`, `partitions_bw_check`, `partitions_flint_607`, `partitions_flint_exceptions`, `partitions_flint_exceptions_new`, `partitions_flint_607_span`)

**A correction.** The seventeenth time said no form of lower level can share the tower's numbers. That was wrong: my
derivation mixed up "equal" and "equal on lane ℓ". A form of lower level never has exactly the class's U₂, U₃ values
(Deligne's bound), but it can have them on lane ℓ (a level-raising congruence), about once in ℓ. Found (counted):
- 17 (weight 14) and 29 (weight 26): the class found on lane ℓ has 2 forms, the exact class 1.
- Along the families, the piece with the class's Hecke numbers on lane ℓ had 3 forms for 17, 4 for 29, 6 for 37 (classes of
  1, 1, 2).
- 499 (weight 496): 22 forms on lane ℓ, 21 exactly (below).

**The fix** (`forms6_flint.py`: `kernel_part`, `padic_kernel`, `class_padic`; the driver's "new part"). The tower's forms are
new at 2 and 3: U₂ and U₃ act on them as exactly (2|ℓ)·2^((k−2)/2) and (3|ℓ)·3^((k−2)/2), the signs carried along the
family. Inside the piece, the ℓ-adic kernel of U₂ − c₂ and U₃ − c₃ keeps exactly these forms. The precision lost is the
depth of the congruence; none was lost here. Checked: 17 and 29 agree with PARI's records, 13 and 37 are unchanged.

**The constants, from modular forms only (counted; FLINT; for each prime the `partitions_flint_family_ℓ*.txt` with the most
digits).** Every n-th difference along the weights was divisible by ℓⁿ, sometimes once more. Every digit shared with an
earlier record agrees.

| ℓ | constants | digits | earlier, all agreeing |
|---|---|---|---|
| 13 | one | 30 | 17 (FLINT), 13 from partition numbers |
| 17, 19, 23 | one | 30 | 13 (PARI), 11 from partition numbers |
| 29, 31 | one | 30 | 11 (PARI), 5 from partition numbers |
| 37, 41 | two on lane ℓ | 20 each | 9 each (PARI; 37 also 12 on FLINT) |
| 43 | two on lane 43 | 20 each | 5 each (PARI) |
| 47 | A ± B·√5 | 20 | 11 (FLINT) |
| 53 | two on lane 53 | 20 each | none (new) |
| 59 | A ± B·√2 | 20 | none (new) |
| 61 | one, and A ± B·√2 | 20 | 10 (FLINT) |

The digits, lowest first:
- 13: 6, 3, 1, 2, 5, 8, 0, 9, 5, 9, 9, 6, 8, 2, 8, 5, 3, 0, 8, 10, 9, 11, 0, 1, 0, 10, 8, 8, 7, 2
- 17: 6, 12, 12, 7, 7, 5, 16, 3, 10, 14, 0, 0, 15, 0, 5, 15, 16, 14, 1, 13, 1, 12, 4, 1, 0, 15, 1, 1, 8, 13
- 19: 10, 2, 10, 8, 10, 1, 16, 5, 7, 10, 6, 10, 8, 2, 18, 0, 8, 12, 12, 11, 14, 14, 0, 1, 3, 14, 10, 17, 10, 18
- 23: 5, 20, 8, 4, 19, 3, 6, 6, 20, 20, 3, 8, 14, 11, 21, 0, 10, 7, 6, 22, 15, 17, 1, 4, 18, 7, 14, 6, 16, 22
- 29: 3, 14, 17, 22, 28, 6, 14, 19, 1, 14, 0, 17, 2, 0, 26, 24, 3, 0, 21, 17, 16, 21, 12, 4, 11, 9, 22, 20, 8, 15
- 31: 30, 6, 2, 28, 1, 8, 3, 15, 15, 13, 16, 5, 3, 6, 27, 18, 2, 9, 28, 2, 17, 23, 10, 26, 8, 3, 1, 4, 0, 19
- 37: 34, 30, 35, 25, 2, 25, 28, 28, 1, 14, 0, 0, 35, 4, 16, 35, 18, 15, 1, 35 and
  15, 4, 4, 3, 1, 10, 35, 5, 36, 21, 9, 32, 7, 11, 17, 4, 20, 20, 22, 18
- 41: 31, 11, 10, 38, 14, 13, 39, 40, 7, 39, 3, 25, 31, 40, 38, 16, 34, 0, 5, 31 and
  23, 14, 7, 38, 12, 26, 7, 1, 6, 21, 27, 31, 12, 14, 29, 15, 37, 0, 14, 2
- 43: 31, 13, 8, 5, 40, 6, 9, 39, 26, 19, 18, 34, 27, 34, 9, 35, 15, 32, 24, 1 and
  6, 20, 10, 25, 3, 13, 24, 4, 39, 22, 32, 8, 38, 8, 11, 39, 19, 41, 17, 42
- 47: A = 38, 4, 33, 43, 39, 34, 9, 21, 12, 37, 23, 4, 0, 31, 12, 12, 23, 9, 31, 19;
  B = 12, 22, 14, 7, 15, 43, 4, 38, 21, 1, 21, 21, 28, 22, 5, 41, 41, 2, 2, 20
- 53: 43, 26, 7, 26, 45, 51, 7, 34, 42, 33, 29, 49, 21, 45, 0, 5, 14, 47, 18, 15 and
  7, 8, 0, 38, 7, 49, 30, 5, 25, 31, 3, 24, 17, 47, 42, 25, 5, 52, 24, 6
- 59: A = 3, 43, 25, 51, 13, 53, 41, 7, 14, 15, 24, 11, 23, 24, 48, 15, 16, 12, 8, 20;
  B = 25, 30, 4, 37, 24, 44, 45, 8, 4, 32, 51, 54, 21, 18, 36, 20, 43, 50, 31, 39
- 61: 32, 40, 3, 28, 17, 11, 0, 28, 10, 17, 38, 29, 38, 40, 48, 10, 1, 38, 23, 2; A = 60, 24, 40, 28, 20, 18, 11, 21, 6, 16,
  53, 60, 4, 18, 4, 15, 41, 38, 17, 7; B = 49, 14, 33, 55, 8, 59, 49, 27, 27, 33, 6, 15, 3, 18, 36, 33, 44, 42, 17, 41

**53 and 59 on exact partition numbers** (`partitions_flint_family_check_53_59.txt`, counted). The laws built from their
constants hold as PREDICTED: lane ℓ at levels (1,3,5) and (2,4,6), lane ℓ² at (3,5,7) and (4,6,8), lane **ℓ³** at (5,7,9),
for n = 0 and n = 1 (53, n = 1: one power deeper at (2,4,6)).

**Against what Boylan–Webb printed** (`partitions_bw_check.txt`; their numbers QUOTED in the literature notes). Each one is
as PREDICTED:
- 53, their relation p(53n+42) ≡ 22·p(53³n+117861) + 25·p(53⁵n+331071432) (mod 53): scaled, it is our polynomial 36, 3, 1.
- 53, their relation (mod 53²) on levels 3, 5, 7, with 2672 and 2304: scaled, it is our polynomial on lane 53², 89, 957, 1.
- Their periods c₃₇ = 36 and c₄₁ = 10: the orders of our polynomials' companion matrices on lane ℓ are 36 and 10.

Their numbers come from partition numbers, ours from modular forms alone, so this is a check from outside.

**The exception 607.** Boylan–Webb (QUOTED): the rank is ⌊(ℓ + 12)/24⌋ for every prime up to 1297 except 607, and
"r₆₀₇(m) = R₆₀₇ − 1", that is 24 instead of 25.
- PREDICTED first (my reading): one of 607's 25 forms is non-ordinary. **REFUTED** (`partitions_flint_607_599_601_607_613.txt`,
  counted): all 25 forms are there and none has T₆₀₇ ≡ 0; likewise for the controls 599, 601, 613.
- What is there instead (counted): T₆₀₇ on the class has a **repeated root** on lane 607. 599, 601, 613 have none.
- PREDICTED next (derived): on lane ℓ the levels act on the class as (3|ℓ)·T_ℓ does (each constant's first digit is
  (3|ℓ)·a_ℓ; the eleventh time). So the space they fill is spanned by v, Tv, T²v, …, and has at most as many directions as
  T_ℓ's minimal polynomial has degree. If other Hecke numbers tell the two forms with the same root apart, T₆₀₇ acts on
  both as one number, and that degree is 24.
- **FOUND** (`partitions_flint_607_span_599_601_607_613_499.txt`, counted): at 607 the minimal polynomial has degree 24, the
  repeated root has two independent eigenvectors, and T₅ tells the two forms apart. So the levels fill at most 24
  directions, exactly Boylan–Webb's rank. 599, 601, 613 and 499 (exact class): the degree equals the number of forms.
- In plain words: the staircase counts the forms, and the tower uses all of them unless two share their number at ℓ
  itself. Then the tower cannot tell them apart and loses one. The chance is about 1 in ℓ per prime (ℓ dividing the
  discriminant), about 0.7 primes expected from 37 to 1297 (Σ 1/ℓ; derived); one was found.
- The sweep on lane ℓ only (`partitions_flint_exceptions_37_700.txt`, counted) flagged 499 and 607. For 499 a form of lower
  level agreeing on lane 499 was PREDICTED, and that is what was found (`partitions_flint_exceptions_new_499_607_all.txt`):
  21 of the 22 forms are the class, and the 22nd is a copy of one of them on lane 499. On the 21: no repeated root, no zero
  root.
- **The full sweep, 37 to 1297, with the class taken exactly** (`partitions_flint_exceptions_new_37_1297.txt`; 200 primes,
  counted). Every class has exactly ⌊(ℓ + 12)/24⌋ forms; a form of lower level interferes on lane ℓ only at 499; no form
  is non-ordinary; and T_ℓ has a repeated root only at **607**. That is Boylan–Webb's list of exceptions up to 1297,
  exactly, as PREDICTED.

## Dig 1, a nineteenth time: the lucky primes, counted to a million (`partitions_flint_lucky`)

**The question.** A law at a prime Q (p(N) ≡ 0 on lane ℓ for a whole class of N built from Q) needs every form of ℓ's
class to take one particular number at Q: (3|Q)·a_Q ≡ ε·(−3|Q)·Q⁻² for a class law (ε = ±1), or ≡ 0 for Ono's kind (the
eleventh time). How often does that happen?

**PREDICTED before the run (derived).** By Chebotarev's density theorem, if each form's mod-ℓ Galois image is as large as
possible and the forms are independent, a form takes a given number c with chance (q + χ(c² − 4D))/(q² − 1). Here q = ℓ
(or ℓ² for a conjugate pair), D ≡ Q⁻³ is the determinant, and χ = ±1 says whether c² − 4D is a square. A law needs every
form at once, and there are three kinds (two classes and Ono's). So about 3 primes Q in ℓ for 13 to 31, 3 in ℓ² for 37 to
59 (1 in 456 for 37), and 3 in ℓ³ for 61. The expected numbers were summed prime by prime.

**Counted** (every prime Q from 5 to 10⁶, 78,495 of them; `partitions_flint_lucky_37_43_41_47_53_59_61_13_17_19_23_29_31.txt`):

| ℓ | forms | expected (class ε=+1, ε=−1, Ono's) | counted |
|---|---|---|---|
| 13 | 1 | 6036, 6036, 6074 | 5940, 6192, 6189 |
| 17 | 1 | 4616, 4616, 4633 | 4658, 4568, 4611 |
| 19 | 1 | 4131, 4131, 4143 | 4085, 4220, 4114 |
| 23 | 1 | 3413, 3413, 3420 | 3389, 3414, 3354 |
| 29 | 1 | 2707, 2707, 2710 | 2702, 2809, 2680 |
| 31 | 1 | 2532, 2532, 2535 | 2512, 2494, 2428 |
| 37 | 2 | 57.4, 57.4, 57.5 | 57, 63, 54 |
| 41 | 2 | 46.7, 46.7, 46.8 | 45, 59, 52 |
| 43 | 2 | 42.5 each | 35, 49, 37 |
| 47 | a conjugate pair | 35.5, 35.5, 35.6 | 32, 31, 35 |
| 53 | 2 | 28.0 each | 32, 34, 22 |
| 59 | a conjugate pair | 22.6 each | 20, 15, 27 |
| 61 | 1 and a pair | 0.3 each | 0, 0, 0 |

- Over the 36 counts with more than one expected, χ² = 37.7 against 36 (derived). That is the scatter chance gives.
- The finer prediction for one form: the share is 1/(ℓ − 1) where χ = +1 and 1/(ℓ + 1) where χ = −1. Over the 48 cells
  Σ(O − E)²/E = 55.0; a flat 1/ℓ gives 239.2 (derived). *Corrected 20 Nov 2026: this read 58.0, which the saved counts
  do not reproduce, and it went on to conclude that the Galois images are the full ones. A fit cannot show that. The
  reason is Ahlgren–Allen–Tang's Proposition 3.3 and Lemma 3.7 with a check at Q = 5; see STATEMENT.md, E2.*
- Single forms check out too ("some form has the number 0", counted against expected): 37: 4196 / 4188.5; 41: 3801 / 3784.3;
  43: 3594 / 3610.3; 53: 2928 / 2935.1; 61: 1294 / 1307.8.
- Every law known from partition numbers is among them (37: 2999, 6599, 7541, 9547; 43: 4049, 8941). Below 10,000 there are
  no others, exactly the lists of `partitions_split.txt`. All the lucky Q up to 10⁶ for 37 to 61 are listed in the output.

**In plain words.** Which primes Q are lucky for ℓ looks random, one by one. But their share follows a law: about 3 in ℓ^d,
where d is the number of forms, corrected by ±1 according to a square-or-not test on Q. Counted over 78,495 primes for each
of 13 primes ℓ, the law fits as well as chance allows.

## Dig 1, a twentieth time: the tower's loops (`partitions_loops`, `partitions_loops_survey`, `partitions_loops_stats`)

**The question (the user's).** Two levels up multiply the tower by its constants, so on lane ℓ^m the tower comes back to
itself after a fixed number of two-level steps: the order of the constants there. Boylan–Webb printed two such loops
(36 for 37, 10 for 41; `partitions_bw_check.txt` agreed). Do the loops show deeper patterns?

**Every prime from 5 to 61, at every depth** (`partitions_loops.txt`; counted from the constants of the eighteenth time).

| ℓ | 5 | 7 | 11 | 13 | 17 | 19 | 23 | 29 | 31 | 37 | 41 | 43 | 47 | 53 | 59 | 61 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| loop on lane ℓ (two-level steps) | 2 | 6 | 5 | 12 | 16 | 18 | 22 | 28 | 2 | 36 | 10 | 21 | 2208 | 26 | 870 | 3720 |

- **Each further digit makes the loop exactly ℓ times longer, for every one of these primes** (C^c = I + ℓ·E with E not 0
  on lane ℓ). So on lane ℓ^m the loop is c·ℓ^(m−1) two-level steps: Boylan–Webb's bound 2c·ℓ^(m−1) levels is reached, and
  no tower comes back early at a deeper digit (one such case was expected by chance; none was found).
- **The shortest loops close on exact partition numbers** (FLINT, n = 0..9, b = 1 and 2; H′ divided by the laws' power
  for 5, 7, 11). 5 and 31 have the constant −1 on lane ℓ: two levels flip the sign, four bring the tower back. 11 comes
  back after 10 levels, 7 after 12 (numbers up to p(7¹⁴·9 + δ₁₄)). All eight checks: exactly as PREDICTED.

**Is there a pattern across primes?** From 7 to 29 every loop is the longest possible (the constant is a primitive root),
six in a row, about a 1-in-300 chance. PREDICTED (derived): luck. The survey of every prime from 13 to 700
(`partitions_loops_survey_13_700.txt`: the constants on lane ℓ as the roots of (3|ℓ)·T_ℓ on the class; counted) and its
statistics (`partitions_loops_stats.txt`; expectations derived for random units):
- 119 constants on lane ℓ: 44 primitive roots against 45.7 expected. The run from 7 to 29 was luck.
- The loop as a share of the longest possible (1, 1/2, 1/3, 1/4, 1/5, 1/6, smaller): 44, 40, 7, 9, 2, 4, 13 counted against
  45.7, 34.0, 7.4, 8.9, 1.8, 5.7, 15.6 expected.
- 59 conjugate pairs: 21 with the longest possible loop, against 16.0 expected.
- Special values (±1, ±2, ±3, ±6 and their inverses; 0.92 expected each): at most 3 counted (6, at 13, 17 and 43; none
  above 43), within chance for 12 values tried.
- Towers with two constants on lane ℓ whose loops are equally long: 9 of 27, against 7.4 expected.

**In plain words.** One firm law: each digit deeper, the loop grows exactly ℓ-fold. Across primes the loops behave like
those of random numbers, as the lucky primes do: no prime is predictable on its own, but the statistics are exact.

## Dig 1, a twenty-first time: 5, 7 and 11 on FLINT, and new laws checked (`partitions_flint_slopes`, `partitions_flint_lawprimes_check`, `partitions_flint_lucky_check`)

**5, 7 and 11 on FLINT** (`partitions_flint_slopes_5.txt`, `_7.txt`, `_11.txt`; counted, from modular forms only). For these
primes the class is empty and the tower follows the smallest positive slope (STATEMENT.md, C). New in `forms6_flint.py`:
`berkowitz` (the characteristic polynomial without division, on any lane) and `slope_root` (the Newton polygon and Hensel's
lemma for the root of slope s). At each weight −2 + (ℓ − 1)·ℓ·(x + 1): the piece with the system's Hecke numbers (5:
a_Q ≡ Q²(1 + Q); 7 and 11: the weight-16 and weight-28 forms, each the only form of its class), lifted over the ℓ-adic
integers; inside it the forms new at 2 and 3 with the class's signs; T_ℓ on them; the one root of slope s.
- **At every weight exactly one root of slope s**, as PREDICTED: 30 weights for 5 (to weight 598), 20 for 7 (to 838), 12
  for 11 (to 1318). The other forms there have larger slopes: at weight 198 for 5, 7, 9 and five more above 9, as
  Calegari's printed list begins (2, 7, 9, …).
- **5:** 33 digits, all 33 settled by the partition numbers agree (the series reaches 30; PARI had 19).
- **7:** 20 digits: 5, 6, 1, 1, 2, 1, 3, 1, 0, 6, 5, 3, 1, 1, 5, 3, 3, 6, 4, 2. All 14 settled ones agree (PARI had 13 of 14).
- **11:** 12 digits: 9, 3, 2, 0, 1, 7, 9, 6, 8, 4, 5, 0. The 4 settled ones agree. PARI's three-weight series had printed
  9, 3, 2, 0, 5, 8, past its reach.

**The new digits of 7 and 11 on deeper partition numbers** (`partitions_flint_lawprimes_check.txt`; FLINT, n = 0, H′ divided
by the laws' power; the ratio two levels apart, r_b = H′_(b+2)/H′_b). PREDICTED: at least 17 digits for 7 and 6 for 11.
- 7: r shares 1, 2, 4, 5, 7, 8, 10, 11, 13, 14, 16, 17, 19, 20, 20, 20 digits with the forms' constant at the level pairs
  (1,3) to (16,18). **All 20 confirmed**, and the deepest pairs agree with each other on two digits more: 6, 2.
- 11: r shares 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7 digits at (1,3) to (13,15), up to p(11¹⁵ + …) ≈ p(4·10¹⁵). **7
  digits confirmed: 9, 3, 2, 0, 1, 7, 9.** The forms' 5th and 6th digits (1, 7) are right; PARI's (5, 8) were not.
- So the forms side now agrees with every digit the partition numbers settle, for every prime from 5 to 61.

**New lucky primes on exact partition numbers** (`partitions_flint_lucky_check.txt`; counted). For each of 37, 41, 43, 47,
53, 59 and each kind of law, the two smallest new Q from `partitions_flint_lucky_*.txt`, two members each:
**36 of 36 laws hold** (all members multiples of ℓ), as PREDICTED. The members reach p(16,702,856,902,164,337) (37, Ono's
kind, Q = 64171). Of the 18 controls one is 0 by chance (1/ℓ each). Among these laws: 53 at Q = 467 and 3797 (class
ε = −1), the smallest found for any prime from 37 on, and 59 at Q = 4091 (Ono's kind).

## Direction 1: is the recipe a general law? 3-colored partitions (`forms2_flint`, `colors3_census`, `colors3_tower`, `colors3_compare`, `colors3_family`, `colors3_digits_check`)

**The question (the user's choice).** p₃(n), the coefficients of ∏(1 − qⁿ)^(−3), are partitions in 3 colours. Belmont, Lee,
Musat and Trebat-Leder (Monatsh. Math. 173 (2014); QUOTED from arXiv 1510.01202) proved that their towers
p₃((ℓᵇn + 3)/24) settle like the partition tower, with rank at most R_ℓ(3) = ⌊(ℓ − 1)/6⌋ − ⌊(ℓ² − 1)/(8ℓ)⌋ for ℓ ≥ 8. As far
as their text was read, no constants, forms, Shimura lifts or families appear.

**PREDICTED before the runs (derived).** The recipe carries over with the Shimura weight of η^(−3): forms of weight ℓ − 5
and their families read at weight −4. And since η³ is a theta series of 8 (Jacobi), not of 24, on **level 2** instead of
level 6.

**The census (`colors3_census.txt`; exact, counted).** The newforms of weight ℓ − 5 on Γ0(2) with U₂ = (2|ℓ)·2^((ℓ−7)/2)
(the same sign rule at 2 as for ordinary partitions) number exactly R_ℓ(3) for **all 50 primes from 11 to 251**. The four
classes on level 6 match only 11 to 16 of 32. The engine for level 2 is new (`forms2_flint.py`: the basis η(2z)¹⁶/η(z)⁸,
2E₂(2z) − E₂(z), (η(z)/η(2z))²⁴; Hecke matrices; the class exactly, also on lane ℓ).

**The towers from exact values (`colors3_tower_N150000000_….txt/.json`; counted).** p₃(n) for n < 1.5·10⁸ on lane ℓ^K, as
1/∏(1 − qⁿ)³ (Jacobi's series inverted by FLINT). This ran in WSL with the same FLINT 0.9.0. The Windows build crashes on
products of more than about 2·10⁷ terms, and 3·10⁸ terms needed more memory than WSL has (counted).
- **11: p₃(11ᵇn + δ_b) is divisible by exactly 11ᵇ** at every level up to 7. **17: by 17^⌈b/2⌉** up to level 6. These
  are laws like Atkin's and Watson's for ordinary partitions (slopes 2 and 1, by the same reading), where R_ℓ(3) = 0. The
  first steps (mod 11, mod 17) are classical Ramanujan-type congruences for p₃; whether the powers are in print was not
  checked.
- 13, 19, 23, 29, 41 (one form): **one constant** for every level pair and n (45, 36, 37, 27, 24 checks).
- 31, 37 (two forms): no single constant fits (0 of 3 level pairs), but a three-term law does (below).

**The first digits against the forms (`colors3_compare.txt`; counted).** The constant is ≡ t(ℓ)·a_ℓ(f) on lane ℓ with
t(ℓ) = ±1 for every one-form prime, and the three-term law holds with the class's two numbers for 31 (6 of 6) and 37 (2 of
2) for one sign. The signs are **t(ℓ) = (−1|ℓ) at all 7 primes** where they can be read. That is the character of η³'s
theta series Σ(−1)ⁿ(2n + 1)q^((2n+1)²), as (3|ℓ) is η's (Calegari's χ₁₂). Every other character tried fits at most 5.

**The constants from the forms (`colors3_family_ℓ.txt`; counted, FLINT only).** The Hida family through each class form, on
level 2, at the weights ℓ − 5 + (ℓ − 1)x, read at x = −1 (weight −4). Every n-th difference was divisible by ℓⁿ.
- (−1|ℓ)·u(−4), lowest digit first: 13: 4, 4, 12, 7, 6, 1, 9, 11, 12, 8, 10, 2, 2, 12, 6; 19: 7, 11, 9, 18, 8, 10, 12, 14,
  2, 7, 14, 18, 0, 3, 15; 23: 18, 18, 2, 11, 6, 16, 11, 20, 19, 18, 19, 15, 3, 12, 8; 29: 5, 11, 0, 20, 20, 22, 28, 25, 21,
  17, 11, 8, 2; 41: 7, 32, 37, 12, 10, 18, 10, 8, 34, 1, 30, 31, 14.
- For 31, 37, 43, 47 the polynomial of the two constants, to 13 digits.
- **Against the exact towers (`colors3_digits_check.txt`):** the ratio two levels apart shares 1, 2, 3, 4, 5 digits with
  the forms' constant at the level pairs (1,3) to (5,7) for 13, and 1, 2, 3, 4 for 19 and 23: one more digit per level,
  exactly as for ordinary partitions. 29 and 41 reach 2. Every settled digit agrees.

**In plain words.** The recipe is not a fact about partitions alone. For 3 colours the same construction works with three
changes, each with a reason: level 2 for level 6, weight −4 for −2, and the twist of η³ for that of η. It predicts the
rank for 50 primes and every digit the exact numbers settle. 11 and 17 play the part that 5, 7 and 11 play for ordinary
partitions.

**5 colours, the census (`colors5_census.txt`; exact, counted).** PREDICTED (derived): η⁵ has no theta-series shortcut (5 is
prime to 3), so the forms return to level 6, at weight ℓ − 7 (≡ −6), in one Atkin–Lehner class with a fixed sign rule.
Found: U₂ = −(2|ℓ)·2^((k−2)/2), U₃ = −(3|ℓ)·3^((k−2)/2) gives exactly R_ℓ(5) forms for **all 33 primes from 11 to 157**.
The other three rules on level 6 match 13 to 17 of 33, and level 2 matches 16 or 17. Both signs are flipped from the
ordinary partitions' rule, as (2|5) = (3|5) = −1 would flip them (a reading, not tested beyond r = 1, 5).

**So far, for r = 1, 3, 5 colours:** the rank of the tower (proven bound) equals the number of forms of weight ℓ − r − 2 in
one class, on level 6 when 3 does not divide r and on level 2 when it does. For r = 1 and 3 the constants are the class's
families read at weight −(r + 1), twisted by the character of η^r's theta series. That is now checked against exact values
for both.

## Direction 1, continued: 5 colours, and the laws of 3 colours (`colors_tower`, `colors6_family`, `colors_check`, `colors5_check23`, `colors3_slopes`, `colors3_laws_check`)

**5 colours, the towers from exact values** (`colors_tower_r5_N120000000_….txt/.json`; p₅(n) for n < 1.2·10⁸ as
1/(J·E₁²), J Jacobi's cube, E₁ Euler's product; WSL; counted). The first values 1, 5, 20, 65, 190, 506, 1265, 2990 are p₅'s.

**5 colours, the forms** (`colors6_family_r5_ℓ.txt`; level 6, weight ℓ − 7, the class U₂ = −(2|ℓ)·…, U₃ = −(3|ℓ)·… of
`colors5_census.txt`; families read at weight −6; FLINT only). 15 digits for 13, 17, 19, 13 for 23, 29, 31, the two
constants' polynomial for 37 and 41. For 17 a form of lower level shared the numbers on lane 17 and was set aside exactly.

**5 colours, the comparison** (`colors_check_r5_N120000000_….txt`, `colors5_check23.txt`; counted).
- PREDICTED twist (3|ℓ), the character of η⁵ as of η: **it fits all 8 primes** where it can be read (13, 17, 19, 23, 29, 31,
  37, 41); every other character tried fits at most 5.
- One constant for every level pair and n at 13, 17, 19, 29, 31; the digits shared with (3|ℓ)·u(−6) grow with the level:
  13: 1, 2, 4, 4, 5 at (1,3) to (5,7); 17: 1, 2, 3, 4; 19: 1, 3, 3, 5; 29, 31: up to 2.
- 37 and 41 (two forms): the three-term law with the family's polynomial holds on lane ℓ for the predicted sign only.
- **11: a law** — p₅(11ᵇn + δ_b) is divisible by 11^(⌊b/2⌋+1) at levels 1 to 7 (R_11(5) = 0): Watson's pattern, slope 1.
- **23, the "superexceptional" pair (5, 23)** of Belmont–Lee–Musat–Trebat-Leder (QUOTED: rank 1 although a congruence holds):
  every value is divisible by 23 exactly once, at every level. Divided by 23, one constant from level 2 on, and it agrees
  with u(−6) (2 to 3 digits at (3,5)).

**3 colours, the laws at 11 and 17 from the forms** (`colors3_slopes_11.txt`, `_17.txt`; level 2, the weights
−4 + (ℓ − 1)·ℓ·(x + 1); FLINT only).
- **17: exactly one root of slope 1 at every weight** (7 weights, to 1900), as PREDICTED. (−1|17)·u(−4)/17 = 4, 12, 7, 2,
  10, 4, 8 (the series' reach). The exact towers (`colors3_laws_check.txt`, 60-bit lanes) confirm 4, 12, gaining about one
  digit per two levels, as for Watson's 7.
- **11: two roots of slope 2 at every weight** (8 weights, to 876), not one. Their constants agree on lane 11 (the
  polynomial is (x + 5)² there): a repeated root. So the tower should obey a three-term law built from the pair, and no
  single constant. **Found on exact values:** the law H′_(b+4) + c₁·H′_(b+2) + c₀·H′_b ≡ 0 (c₀, c₁ from the forms, to 11⁸)
  holds to 11² at the levels (1,3,5), 11⁴ at (2,4,6) and 11⁶ at (3,5,7), for every n; the ratio two levels apart takes
  different values for different n, as two constants must. This is new against ordinary partitions, where 5, 7, 11 each
  follow one form.

**Where the generality stands (r = 1, 3, 5 colours).** The recipe holds for all three: the rank of the tower equals the
number of forms of weight ℓ − r − 2 in one class (level 6 when 3 ∤ r, level 2 when 3 | r); the constants are those forms'
families read at weight −(r + 1), twisted by the character of η^r ((3|ℓ) for r = 1 and 5, (−1|ℓ) for r = 3). Where the class
is empty, families of positive slope give the laws, and sometimes, as at 11 for 3 colours, a pair of them does.

**607's two forms, one digit deeper** (`partitions_flint_607_pair.txt`; counted). On the two-dimensional piece where
T₆₀₇ (weight 604) has its repeated root, T₆₀₇ = 367·I + 607·C. C's polynomial x² + 70x + 484 is irreducible on lane 607,
and so is T₅'s there (x² + 561x + 66). So the two forms are a conjugate pair over the field with 607² elements: at a prime
of residue degree 2 above 607 in their Hecke field, each is the other's Frobenius conjugate. Their numbers at 607, two
conjugates in F_(607²), happen to fall into the small field F₆₀₇ and so coincide there (367). For an element of F_(607²)
that happens with chance 607/607² = 1/607, and about 0.7 such primes were expected from 37 to 1297. One digit deeper they
are again conjugates, and different. Reading (derived): on lane 607^m the tower then holds the 25th direction only with a
factor 607 (generated by v and 607·C·v), which fits Boylan–Webb's "rank R − 1 for all m". [13 Nov: that reading used the forms at weight 604. The tower's own constants were then computed to part on the second digit ("607's pair, one digit deeper in the tower itself"). So modulo 607² the levels need 25 generators. Whether that fits Boylan–Webb's "for all m" depends on what their rank counts; their second-digit computations stopped at ℓ = 89.]

## Where else: 24 colours, the reciprocal of Δ (`forms1_flint`, `forms1_family`, `colors_tower` r = 24, `colors24_check`, `colors24_check43`)

**Why this case.** Π(1 − qⁿ)^(−24) = q/Δ, so the tower p₂₄(ℓᵇn + 1) is c(ℓᵇn) for the coefficients c of 1/Δ: the operator
U_ℓ applied b times to a weakly holomorphic form of weight −12 and level 1. There is no half-integral weight and no
Shimura step. The classical route of PROOF.md then predicts the constants directly: the U_ℓ-values at weight −12 of the
ordinary families of level 1, one level per step.

**PREDICTED before the runs (derived).**
- The number of constants is dim M_(ℓ−13)(SL₂(ℤ)), cusp forms and the Eisenstein series together. This equals
  Belmont–Lee–Musat–Trebat-Leder's bound R_ℓ(24) = ⌊(ℓ − 1)/12⌋ for every ℓ ≥ 17 (a short derivation).
- One constant is exactly 1, from the Eisenstein family, whose U_ℓ is 1 at every weight. It is visible in p₂₄(1) = 24 at
  every level.
- The others are A_f(−12) for the ordinary cusp forms of weight ≡ −12 (mod ℓ − 1).

**The forms** (`forms1_family_w-12_ℓ.txt`; level 1, new engine `forms1_flint.py`: the basis Δ^(d−1)·E_w·j^i, checked
against τ(2) = −24, τ(3) = 252 and the weight-24 pair 540 ± 12√144169):
- The factor x − 1 appears exactly at every prime.
- 5, 7, 11, 17, 19, 23: only the Eisenstein family is ordinary there, so the only constant is 1.
- 13: 1, and Δ's own family read at −12: 8, 2, 2, 10, 6, 1, 2, 9, 2, 1, 1.
- 29: 1, and 23, 26, 17, 22, 15, 20, 4, 4, 15, 4, 9 (the family of ΔE₄).
- 31: 1, and 21, 6, 20, 7, 7, 21, 15, 28, 19, 1, 21 (ΔE₆).
- 37: 1 and a conjugate pair. 41: 1, 33, 34 on lane 41.
- 43: 1, 30, 30. The two cusp forms of weight 30 share their number at 43.

**Against exact values** (`colors24_check.txt`; p₂₄(n) for n < 10⁷ and 10⁸, WSL; counted). The law Σ [xⁱ]P·H_(b+i) ≡ 0
holds at every prime and every n, one digit deeper per level:
- 13, 17, 19, 23, 29, 31: depths 1, 2, 3, 4 from levels 1, 2, 3, 4;
- 5: two digits deeper per level;
- 37: depths 1, 2; 41 and 43: depth 1, the deepest the data reach.
- At 5, 7, 11, 17, 19, 23 that law is simply H_(b+1) ≡ H_b: **the 24-colour tower stops changing.** For example
  p₂₄(17ᵇ⁺¹n + 1) ≡ p₂₄(17ᵇn + 1) (mod 17ᵇ).

**A 607-type exception at 43, predicted and confirmed** (`colors24_check43.txt`; counted). The pair of weight 30 has Hecke
field ℚ(√51349), inert at 43, and its number at 43 happens to fall in F₄₃, as 607's pair did. So the shorter law with only
the constants 1 and 30, H_(b+2) − 31·H_(b+1) + 30·H_b ≡ 0, should hold on lane 43. **It holds for all 24 (b, n) checked.**
At 41, where the constants differ, none of the three shorter laws holds. So on lane 43 the 24-colour tower fills 2
directions, not 3: the forms found the exception before the partition numbers were looked at.

## Where else: spt, Andrews' smallest parts (`spt_tower`, `spt_check`)

**Why this case.** spt(n) counts the smallest parts over all partitions of n. 12·spt(n) + (24n − 1)·p(n) are the
coefficients of a mock modular form of weight 3/2. Its Shimura weight is 2, a classical weight, and there the recipe
predicts the simplest possible constants: ±1.

**PREDICTED before the run (derived; the header of `spt_tower.py`).**
- Every constant is ±1.
- On lane ℓ, H_(b+2) ≡ (3|ℓ)·H_b, where H_b(n) = spt(ℓᵇn + δ_b) and 24δ_b ≡ 1 (mod ℓᵇ).
- At 5, 7 and 13 there are laws.

**Measured** (`spt_tower_N16000000_….txt`, `spt_check.txt`; counted). spt(n) for n < 1.6·10⁷ comes from Andrews' identity,
checked against brute force for n ≤ 20. Each level uses n = 0…11. Depths are of H_(b+2) − (3|ℓ)·H_b, smallest over n:

| primes | depths at levels b = 1, 2, 3, … |
|---|---|
| 11, 17, 19, 23, 29, 31, 37, 41, 43 | 1, 1, 2, 2, as far as the data reach |
| 5, 7, 13 | the laws ℓ^⌊(b+1)/2⌋ hold at every level |

After dividing out those laws (H′):

| prime | depths of H′_(b+2) − (3\|ℓ)·H′_b |
|---|---|
| 5 | 1, 3, 4, 6, 7, then the end of the lane |
| 7 | 1, 3, 3, 5, 5 |
| 13 | 0, 1, 1, 2 |

**Already proven (quoted; read 12 Nov 2026).**
- Ahlgren, Bringmann and Lovejoy, "ℓ-adic properties of smallest parts functions", Adv. Math. 228 (2011).
  - Thm 1.1 (ii), for every prime ℓ ≥ 5 and m ≥ 1: spt((ℓ^(2m+1)n + 1)/24) ≡ (3/ℓ)·spt((ℓ^(2m−1)n + 1)/24) (mod ℓ^m).
  - Thm 1.2: M | T(ℓ^(2m)) − (3/ℓ)·M | T(ℓ^(2m−2)) ≡ 0 (mod ℓ^m).
  - Their Prop. 2.1 (i) turns Thm 1.2 into the even levels too: H_(2m+2) ≡ (3|ℓ)·H_(2m) (mod ℓ^m) for every n. This is a
    two-line step of ours, derived. So at every level: **H_(b+2) ≡ (3|ℓ)·H_b (mod ℓ^⌈b/2⌉)**.
- Garvan, "Congruences for Andrews' spt-function modulo powers of 5, 7 and 13", Trans. AMS 364 (2012), arXiv 1011.1955,
  for a ≥ 3:
  - (1.6) spt(5ᵃn + δ_a) + 5·spt(5^(a−2)n + δ_(a−2)) ≡ 0 (mod 5^(2a−3));
  - (1.7) the same with 7, modulo 7^⌊(3a−2)/2⌋;
  - (1.8) with 13 and the sign −, modulo 13^(a−1).
  - (1.9)–(1.11): the laws ⌊(a+1)/2⌋.
- The two papers index n exactly as we do. For the same n, levels b and b + 2 read spt at ℓᵇ·(24n + t), with t fixed by
  the parity of b (derived).

**Against the theorems (counted).** The measured smallest depth equals the proven modulus at every prime and level
computed:
- 11–43: ⌈b/2⌉ = 1, 1, 2, 2, and measured 1, 1, 2, 2.
- In Garvan's form, with the law's power added back:
  - 5: 3, 5, 7, 9, 11 for a = 3…7, and measured the same;
  - 7: 3, 5, 6, 8, 9, and measured the same;
  - 13: 2, 3, 4, 5, and measured the same.

So for some n ≤ 11 each congruence fails one power higher. Neither paper's modulus can be raised at these levels. We have
not seen this stated in the parts we read.

**Reading (derived).** All of it is one rule. Two levels multiply the tower by (3|ℓ)·r, where r is a root of
x² − (1 + ℓ)x + ℓ. That is the Hecke polynomial of the Eisenstein series of weight 2, the Shimura weight. Its roots are:
- 1: the ordinary family, whose U_ℓ is 1 at every weight. At 11–43 the tower follows it, which gives ABL's (3|ℓ).
- ℓ: slope 1. At 5, 7 and 13, Andrews' congruence spt(ℓn + δ) ≡ 0 removes the part on the root 1, and the tower follows
  ℓ. Garvan's −5, −7 and +13 are exactly (3|ℓ)·ℓ.

**What this means for the study.**
- The recipe predicted the right constant. In this case the answer is a theorem: spt is the Eisenstein case, and p(n),
  whose constants are values of cusp-form families, is the cusp-form case.
- Nothing about spt is new here except the count showing the moduli are sharp.
- ABL already noted the likeness to Folsom–Kent–Ono's periodicity for p(n). They wrote that the phenomenon "will hold for
  a much wider class of modular forms and mock modular forms".

## Where else: two more forms of level 1, j and 1/E₄ (`level1_tower`, `forms1_family` at w = 0 and −4, `level1_check`, `level1_j_cusp`, `level1_pole_check`)

**Why these.**
- j = E₄³/Δ has weight 0 and its only pole is at the cusp, like 1/Δ (24 colours).
- 1/E₄ has weight −4. It is holomorphic at the cusp but has a pole inside, at ρ = e^(2πi/3), the point of the curve
  y² = x³ + 1 (j = 0).
- Both towers are H_b(n) = c(ℓᵇn), one level per step.
- `level1_tower.py` builds E₄ from Euler's series alone: E₄ = η(τ)¹⁶/η(2τ)⁸ + 256·η(2τ)¹⁶/η(τ)⁸, checked against
  1 + 240σ₃(n) for n < 2000. The first values were checked against the known ones: j: 744, 196884, 21493760, …;
  1/E₄: 1, −240, 55440, −12793920, ….
- n < 1.6·10⁷, primes 5…43.

**The forms** (`forms1_family_w0_ℓ.txt`, `forms1_family_w-4_ℓ.txt`; counted).
- At both weights and every prime, one constant is exactly 1 (the Eisenstein family); the others are cusp forms'
  families.
- One coincidence stands out at weight −4 and ℓ = 37: all three constants are ≡ 1 on lane 37. 37 divides the Bernoulli
  number B₃₂, so a cusp form of weight 32 ≡ −4 (mod 36) is congruent to the Eisenstein series modulo 37.

### j

**PREDICTED before the run (derived).** The forms' law Σ [xⁱ]P·H_(b+i) ≡ 0 holds, deeper level by level. At 5, 7 and 11
it is H_(b+1) ≡ H_b.

**Measured** (`level1_check_j_….txt`; counted). It holds wherever the data reach:
- 5, 7, 11, 13, 17, 19, 23: depth b from level b; at 5, b + 1.
- 29, 31: depth 1 from level 1.
- 37, 41, 43: too few levels for the full law.
- In every case the least law the data alone allow on lane ℓ is the forms' P, except at 31, where 744 = 24·31 hides the
  constant 1 on lane 31.

**Sharper: for n ≥ 1, the cusp forms alone** (`level1_j_cusp_….txt`; counted).
- Why (derived): at weight 0 the Eisenstein family's member is the constant 1, since ζ*(1 − k) has its pole at k = 0.
  So that family carries only j's constant term 744.
- Prediction: for n ≥ 1, the law is C = P/(x − 1).
- Result: it holds at every prime 5…43, with depth b from level b (b + 1 at 5). 37, 41 and 43 are checkable this way
  (depth 1 from level 1).
- At 5, 7 and 11, C = 1. The law then says c(ℓᵇn) → 0: c(5ᵇn) ≡ 0 (mod 5^(b+1)), c(7ᵇn) ≡ 0 (mod 7ᵇ),
  c(11ᵇn) ≡ 0 (mod 11ᵇ). **These are Lehner's congruences (1949): "no ordinary cusp form of weight ℓ − 1".**
- At 13, 17, 19, 23 there is one constant. For every n ≥ 1 and b ≥ 1, c(ℓ^(b+1)n) ≡ A·c(ℓᵇn) (mod ℓᵇ), where A is the
  family of the one cusp form (Δ, ΔE₄, ΔE₆, ΔE₁₀) read at weight 0. A on lane ℓ is 8 at 13, 7 at 17, 6 at 19, 5 at 23.
  13-adic digits of A: 8, 10, 11, 9, 11, 0, 9, 7.
- The step from level 0 is not part of the law: c(13n)/c(n) on lane 13 varies with n.
- Atkin and O'Brien (Trans. AMS 126, 1967) studied c(n) modulo powers of 13. Their paper is unread: the AMS site refuses
  automated reads, and the 2025 proof of their Hecke-operator conjecture for p(n) by Garvan and Shomanov (arXiv
  2504.10824) does not restate the c(n) part. The 13 case may be classical. For ℓ ≥ 17, a law with these constants was not found in what we read.

### 1/E₄: the pole brings one constant of its own

**PREDICTED before the run (derived).**
- ℓ ≡ 2 mod 3: the forms' law alone. The curve at ρ is then supersingular, E₄ divides a power of the Hasse invariant on
  lane ℓ, and 1/E₄ is an ℓ-adic modular form.
- ℓ ≡ 1 mod 3: the curve is ordinary, 1/E₄ is not an ℓ-adic modular form, and the pole may add constants.

**What the first small run showed** (n < 3·10⁵, at 7, 13 and 17; a scratch run whose numbers the kept run below repeats).
- At 13 the forms' law x − 1 fails even on lane 13. The data alone ask for x² − 1.
- At 7 the law holds on lane 7 but never gets deeper.
- The remainders one level apart, G_(b+1)/G_b, give the extra constant. A search over ±u^m, m = −60…60, finds **u⁻⁶**:
  9 digits at 7, 3 at 13.
- Here u is the unit root of Frobenius of y² = x³ + 1 at ℓ, the root of x² − aₗx + ℓ (4ℓ = aₗ² + 3b²) that is an ℓ-adic
  unit. u⁶ is the same for all six twists of the curve.
- Before the full run this became the prediction for 19, 31, 37 and 43, and for 7 and 13 deeper: the law is
  P(x)·(x − u⁻⁶).

**Measured** (`level1_pole_check_invE4_….txt`, `level1_check_invE4_….txt`; n < 1.6·10⁷, counted).
- ℓ ≡ 2 mod 3, ρ supersingular: the forms' law alone holds, deeper level by level.
  - 5: depths 3, 4, …, 11 from levels 0, 1, ….
  - 11: 1…6. 17: 1…4. 23: 1…4. 29: 1, 2.
  - The least law the data allow on lane ℓ is P each time.
  - 41: too few levels.
- ℓ ≡ 1 mod 3, ρ ordinary:

| ℓ | the constant from the data | the law P(x)·(x − u⁻⁶), depth from levels 0, 1, … | the forms' law alone |
|---|---|---|---|
| 7 | u⁻⁶, 9 digits | 2, 4, 6, 8, then all of lane 7¹⁰ | depth 1 at every level |
| 13 | u⁻⁶, 5 digits | 1, 2, 3, 4, 5 | 0 |
| 19, not used to find it | u⁻⁶, all 3 digits the data give; the search picks m = −6 | 1, 2, 3, 4 | 1 at every level |
| 31, not used to find it | 8 = u⁻⁶ on lane 31 | 1, 2 | 0 |
| 43, not used to find it | 35 = u⁻⁶ on lane 43 | 1 (the only level span) | 0 |
| 37 | cannot be read | — | — |

- Found from the data alone, with no forms: the least law on lane ℓ is x² − 1 at 13, and 17 + 16x + 28x² + x³ at 31. The
  latter is exactly P(x)·(x − 8) on lane 31.
- At 7 and 19, u⁻⁶ ≡ 1 on lane ℓ, so the extra constant hides there on lane ℓ and shows only deeper.
- 37 cannot tell: all four constants are ≡ 1 on lane 37. The forms' three are ≡ 1 because of B₃₂, and u⁻⁶ ≡ 1 too.

**Reading (derived; a heuristic of ours, not checked against a paper).**
- Hecke's operator T_ℓ multiplies the residue of a weight-w form at ρ by π^(w−2) + π̄^(w−2). These come from the two
  endomorphisms of norm ℓ of the curve at ρ; the step is to write those two cosets of T_ℓ as multiplication by π and by π̄.
- So the pole behaves like a pair of U_ℓ-roots π^(w−2) and π̄^(w−2). Only the unit one, u^(w−2), can live in a tower of
  integers. For w = −4 the other has ℓ-adic valuation −6.
- That gives u⁻⁶, and it is canonical because 6 divides w − 2.
- At ℓ ≡ 2 mod 3 no element of ℤ[ζ₃] has norm ℓ, so T_ℓ moves the pole away and nothing stays.

**Against the literature (read 12 Nov 2026).** Recent papers treat meromorphic forms of POSITIVE weight with poles at CM
points, and their constants are of the same kind: powers of Frobenius of the pole's curve.
- Zhang, arXiv 2510.23200: E₄/j has the single constant (−3|p)·p.
- Bordignon, arXiv 2601.12157: a law whose constants come from Frobenius on Sym^(w−2) of the curve at the pole.
- Bringmann–Guerzhoy–Kane–Mertens–Rolen, arXiv 2606.14020: at supersingular primes the iterates approach a holomorphic
  cusp form. They report that "such congruences do not seem to occur for primes p≡1 (mod 4)" for Δ/E₆.

The negative weight of 1/E₄, its constant u⁻⁶, and the rule "the forms' constants and u^(w−2)" were not found in what we
read.

## How wide the recipe reaches (12 Nov 2026)

**The group all the cases belong to (quoted theory; the arrangement is ours).** Every tower studied here is
F | U_ℓᵇ for a q-series F that is, or becomes after one U_ℓ, an ℓ-adic modular form (Serre 1973, Katz 1973). On those,
U_ℓ is compact (Coleman). Iterating it keeps only the ordinary part, which has finite rank (Hida's projector
e = lim U_ℓ^(n!)), and decays the part of slope s like ℓ^(s·b). So a law with a few fixed constants must hold, and the
constants are values of Hida families at F's weight.

The known theorems are corners of this one mechanism:
- p(n): Folsom–Kent–Ono, Boylan–Webb;
- r colours: Belmont–Lee–Musat–Trebat-Leder;
- spt: Ahlgren–Bringmann–Lovejoy; at 5, 7, 13, Garvan;
- j at 5, 7, 11: Lehner;
- p(n) at 5, 7, 11: Watson and Atkin.

What this study adds is the concrete constants, and the rules for reading them:
- the Shimura step and its twist at half-integral weight;
- the Eisenstein member being the constant 1 at weight 0;
- an ordinary CM pole adding u^(w−2);
- an empty ordinary part handing over to the next slope.

**Layers, from narrow to wide.**
1. η-quotients: partition functions in any number of colours. Tested here.
2. Weakly holomorphic forms of any level and weight. j is tested here. Also distinct parts, overpartitions and the
   moonshine function T_2B (`group_members_….txt`; with the data-only detector).
3. Mock modular forms. spt is tested here and proven there. Mock theta functions and class numbers are untested.
4. Meromorphic forms with poles at CM points: 1/E₄ is tested here, and 1/E₆ is the test below. In positive weight there
   are Zhang, Bordignon, Bringmann–Guerzhoy–Kane–Mertens–Rolen, and Allen and coauthors (arXiv 2511.05718).
5. Non-congruence modular forms. Here there is no Hecke theory at all, yet Scholl (Invent. Math. 79, 1985) proved
   Atkin–Swinnerton-Dyer congruences. These are recurrences among a(npˢ) with constants from Frobenius on de Rham
   cohomology: the same shape as our tower laws. Kazalicki–Scholl (arXiv 1301.5876) extended them to weakly holomorphic
   forms.
6. Untested: automorphic forms of higher rank, where Hida theory also exists. Examples are Siegel forms, such as
   1/χ₁₀, which counts the black-hole states of string theory (1/Δ, our 24 colours, counts its simplest states), and
   Hilbert forms.

**Reading (derived).** In every layer that has been tested or proven, the constants are unit roots of Frobenius on
something geometric:
- Hida families for congruence forms;
- the curve at the pole for CM poles;
- Scholl's motive for non-congruence forms.

On this reading, a tower law with fixed constants exists exactly when the numbers come from geometry. The boundary test
below asks what happens when they do not.

**The boundary, tested** (`boundary_powers_N3000000_5_7_11_13.txt`, `group_members_N3000000_5_7_11_13.txt`; counted).
- The method: the same data-only search for a law on lane ℓ (degree ≤ 6, unknown constants solved over F_ℓ, from level 1
  and from level 2). No forms are computed. n < 3·10⁶, ℓ = 5, 7, 11, 13.
- Two functions with no modular origin, partitions into squares ∏(1 − q^(k²))⁻¹ and into cubes ∏(1 − q^(k³))⁻¹: **no law
  at any prime, from either level.**
  - That is up to degree 6 at 5, and at 7 from level 1.
  - Where the levels run short the search stops earlier: degree 5 at 7 from level 2, 5 and 4 at 11, 5 and 4 at 13.
- Three modular functions never looked at before. A unique law at 5, 7 and 11 for all three, from 14 to 82 equations:
  - **T_2B:** H_(b+1) ≡ H_b at 5 and 7; degree 2 at 11 and 13.
  - **overpartitions:** H_(b+2) ≡ H_b at 5 (from level 2), 7 and 11; degree 4 at 13, from only 9 equations.
  - **distinct parts:** the tower vanishes on lane 5 (from level 2), the mark of a law prime; H_(b+2) ≡ −H_b at 7;
    x(x² + 2) at 11 (from level 2); at 13 the data are too short for the degree needed.
- So the line between "a law" and "no law" falls where modularity ends, as the reading above says.

### 1/E₆: the pole rule on a second curve (`level1_tower` invE6, `forms1_family` at w = −6, `level1_pole_check2`, `level1_check`)

**Why.** 1/E₆ has weight −6 and its pole at i, the point of y² = x³ − x (complex multiplication by ℤ[i]). That curve is
ordinary at ℓ ≡ 1 mod 4 and supersingular at ℓ ≡ 3 mod 4. E₆ is built from Ramanujan's identity E₆ = E₂E₄ − 3θE₄; the
first values are checked by hand: 1, 504, 270648, 144912096.

**PREDICTED before the run** (from the rule u^(w−2) found on 1/E₄, written in the header of `level1_pole_check2.py`).
- ℓ ≡ 3 mod 4: the forms' P alone.
- ℓ ≡ 1 mod 4: P(x)·(x − u⁻⁸). On lane ℓ, u⁻⁸ ≡ 1 at 5 and 17 (so hidden there), 9 at 13, 7 at 29, 12 at 37, 18 at 41.

**Measured** (`level1_pole_check2_invE6_….txt`, `level1_check_invE6_….txt`; n < 1.6·10⁷, counted).
- ℓ ≡ 3 mod 4: P alone, deeper level by level: 7: 2…9. 11: 2…7. 19: 1…4. 23: 1…4. 31: 1, 2. 43: 1.
- ℓ ≡ 1 mod 4:

| ℓ | the constant from the data | the law P(x)·(x − u⁻⁸), depth from levels 0, 1, … | the forms' law alone |
|---|---|---|---|
| 5 | u⁻⁸ to all 10 digits; the next-best power ±u^m matches 3 | 4, 7, 10, then all of lane 5¹¹ from level 3 | 1 at every level |
| 13 | 9, 2, 7, 6 = u⁻⁸ to 4 digits | 1, 2, 3, 4, 5 | 0 |
| 17 | u⁻⁸ to 3 digits | 1, 2, 3, 4 | 1 at every level |
| 29 | 7, 7 = u⁻⁸ to 2 digits | 1, 2 | 0 |
| 37 | 12 = u⁻⁸ | 1 | 0 |
| 41 | 18 = u⁻⁸ | 1 | 0 |

- Found from the data alone, with no forms: on lane 13 the least law is x² + 3x + 9 = (x − 1)(x − 9). On lane 29 it is
  4 + 19x + 5x² + x³ = P(x)·(x − 7).
- In `level1_check_invE6_….txt` the label "≡ 1 / 2 mod 3" by each prime is the split for 1/E₄; it means nothing for 1/E₆.
- So the rule u^(w−2) holds on a second curve, at every prime that can test it. Every lane-ℓ prediction was written before
  the data were read.

## The second kind of exception, reduced to one form of weight 6 (`second_kind`)

**The question.** 607 is the "first kind": two class forms whose numbers at ℓ agree on lane ℓ. The "second kind" would
be one class form whose number at ℓ is itself divisible by ℓ, a non-ordinary form. Up to 1297 there is none. Can there
ever be one?

**The reduction (derived; mod-ℓ theory of theta cycles, Jochnowitz; Serre weights; Deligne–Serre lifting).**
- Let f be an eigenform of weight k = ℓ − 3 on Γ0(6), with a_ℓ(f) ≡ 0 (mod ℓ). Then f's theta cycle passes through
  filtration ℓ + 3 − k = 6, at g = θ⁴f.
  - Check: for a supersingular ρ̄ with inertia ω₂^(k−1), the twist by ω⁴ has Serre weight 1 + 0·ℓ + 5 = 6.
- g is again an eigenform, non-ordinary, with U₂(g) ≡ 2⁴·a₂(f) and U₃(g) ≡ 3⁴·a₃(f).
- With a₂(f) = −w₂·2^((k−2)/2) and 2^((ℓ−1)/2) ≡ (2|ℓ), this gives U₂(g) ≡ −w₂(f)·(2|ℓ)·4, and likewise
  U₃(g) ≡ −w₃(f)·(3|ℓ)·9.
- g lifts to an eigenform of weight 6 on Γ0(6). Those are (counted, `class_basis_exact`):
  - the Eisenstein series, which are ordinary, so excluded;
  - the level-3 newform η(z)⁶η(3z)⁶ = q − 6q² + 9q³ + …, stabilized at 2 (U₂ a root of x² + 6x + 32, U₃ = +9);
  - **one newform h₆ = q + 4q² − 9q³ + 16q⁴ − 66q⁵ − 36q⁶ + 176q⁷ + …**, with U₂ = +4 and U₃ = −9 (w₂ = −1, w₃ = +1).
    The other three sign classes are empty at weight 6.
- So, for the four classes at weight ℓ − 3, named by (w₂, w₃) relative to ((2|ℓ), (3|ℓ)):
  - **the partition class (−, −) needs U₂ ≡ 4 and U₃ ≡ 9.** No lift has these for any ℓ ≥ 5: h₆ has U₃ = −9, and the
    level-3 form would need ℓ | 72.
  - (+, ±) needs U₂ ≡ −4: no lift (the level-3 form would need ℓ | 24).
  - (−, +) needs U₂ ≡ 4 and U₃ ≡ −9: exactly h₆. So it holds a non-ordinary form exactly when ℓ | a_ℓ(h₆).

**PREDICTED, then:**
- every partition-class form is ordinary at every prime ℓ ≥ 7, so the second kind never happens for p(n);
- non-ordinary forms of weight ℓ − 3 on Γ0(6), in these classes, occur only in (−, +), only where h₆ is non-ordinary.

**Measured** (`second_kind_N10000000_L300_at853_857_859_863.txt`; counted).
- h₆ to 10⁷ terms: below 10⁷ it is non-ordinary at **exactly one prime, 859** (859 | a₈₅₉).
- All four classes for every prime 13 to 300: no form has T_ℓ-eigenvalue 0 on lane ℓ.
- **At ℓ = 859, the prediction set down before the run held exactly.** The class (−, +) has a form whose T₈₅₉-eigenvalue is
  0 on lane 859: a non-ordinary form, h₆'s partner. The partition class (36 forms), (+, −) and (+, +) have none.
- The controls 853, 857 and 863: none in any class.
- So, every prime tested from 13 to 300, plus 853–863, agrees with the reduction. A non-ordinary form appears exactly once,
  in the class and at the prime that h₆ named.

**What this means.**
- **The second kind of exception never happens for p(n), at any prime ℓ ≥ 7** (derived, with its mechanism confirmed by
  the test above). Every partition-class form is ordinary.
- In the proof this replaces "ordinary, checked for every ℓ ≤ 1297" with a short argument. It also explains why
  Boylan–Webb found no such case.
- The classes are cut by Atkin–Lehner signs at 2 and 3, and the one form of weight 6 decides the question for all ℓ
  (to 10⁷ it is non-ordinary only at 859).
- Remaining care: the theta-cycle facts are used for weight k = ℓ − 3 with 3 ≤ k ≤ ℓ − 1, so ℓ ≥ 7. For ℓ ≤ 11 the class
  is empty anyway.

**The same for 3 and 5 colours** (`second_kind_colors.txt`; exact linear algebra, counted; the reduction derived).
- For r colours the class has weight ℓ − r − 2, and the theta cycle leads through θ^(r+3) to weight r + 5. There it needs
  U_p ≡ ε·p^((r+3)/2), with ε = +1 for 1 and 3 colours and −1 for 5.
- Every cusp eigensystem of that weight, on the class's level, misses the required values except modulo a few small
  primes:
  - **1 colour** (weight 6, Γ0(6), needs (4, 9)): h₆ misses modulo 18, the level-3 pair modulo 72. Primes 2, 3.
  - **5 colours** (weight 10, Γ0(6), needs (−16, −81)):
    - the one newform h₁₀ has (−16, +81): misses modulo 162;
    - the level-2 newform (U₃ roots of x² + 156x + 3⁹): misses modulo 8;
    - the level-3 newforms (U₂ roots of x² − 18x + 512, and x² + 36x + 512 with U₃ = −81): miss modulo 6 and 192.
    - Primes 2, 3.
  - **3 colours** (weight 8, Γ0(2), needs U₂ ≡ 8): h₈ = η(z)⁸η(2z)⁸ has U₂ = −8, which misses modulo 16. Prime 2.
- So **the classes of 1, 3 and 5 colours are ordinary at every prime where they exist**: ℓ ≥ 7, 11 and 11, where the
  theta-cycle steps hold. A colour law (as at 11 and 17 for 3 colours) always comes from a family of positive slope at
  higher weight, never from a non-ordinary form of the class's own weight.

## Link 1 of the proof, checked on numbers (`cdhs_check`)

Chen–Du–Hou–Sun's Theorem 3.1, with its exact γ and λ, was checked at (ℓ, b) = (13, 1), (13, 2), (13, 3), (17, 1),
(17, 2), (19, 2), (23, 2), (29, 2) (`cdhs_check_….txt`; counted). Modulo ℓᵇ, φ = F_b/η(24z)^γ is fitted to the basis of
M_λ(SL₂(ℤ)) from its first d coefficients, and **every further coefficient up to 400 agrees: 2,942 of 2,942.**
- The largest case is ℓ = 13, b = 3: λ = 1008, d = 85, lane 13³.
- K_b = γ + 2λ − 1 ≡ −2 (mod (ℓ − 1)ℓ^(b−1)) in every case, as the formulas say: K_b = (ℓ − 1)ℓ^(b−1) − 2 for odd b and
  2(ℓ − 1)ℓ^(b−1) − 2 for even b.

**T_2B against the forms (12 Nov; counted on lane ℓ).** The data-only laws of the moonshine function T_2B
(`group_members_….txt`) match the recipe of j, now at level 2.
- At 5 and 7 the law is x − 1: only the Eisenstein constant, since the families through weights 4 and 6 on Γ0(2) have
  no cusp form.
- At 11 the law is x² + 3x + 7, the constants 1 and 7. The newform of weight 10 on Γ0(2) is
  q + 16q² − 156q³ + 256q⁴ + 870q⁵ − … − 56148q¹¹ + …, and **a₁₁ = −56148 ≡ 7 (mod 11)**: the first digit of its family.
- At 13 the law is x² + 4x + 8, exactly j's at 13: Δ's family, from the oldforms at level 2, and the Eisenstein 1.

## Overpartitions against the forms (`overpartitions_check`)

**Why.** Σ p̄(n)qⁿ = ∏(1 + qⁿ)/(1 − qⁿ) = 1/θ₄, with θ₄ = η(z)²/η(2z). This has weight −1/2, like p(n), but no
η-multiplier shift: it is θ-type, level 16. The data-only laws (`group_members_….txt`), in y = the two-level step, are:
- y − 1 at 5, 7 and 11;
- (y − 1)(y − 7) on lane 13.

**PREDICTED before the run at 17–31 (written 12 Nov, derived from 5–13).**
- The constants are the Eisenstein 1 and ψ(ℓ)·a_ℓ(f) for the newforms f of weight ℓ − 3 on Γ0(2) in one Atkin–Lehner
  class at 2.
- Evidence so far:
  - at 13 the weight-10 newform (w₂ = −1, a₁₃ = 178094 ≡ 7) gives the 7;
  - at 11 the weight-8 form h₈ (w₂ = +1) is absent;
  - at 5 and 7 the class is empty.
- Both "w₂ = −1 always" and "w₂ = (2|ℓ)" fit 5–13, and so do ψ = 1 and ψ = (−1|·). The primes 17 ((2|17) = +1), 19, 23
  and 31 (≡ 3 mod 4) separate them.

**Measured at 13–31** (`overpartitions_check_N16000000_13_17_19_23_29_31.txt`, and a data-only law search; counted on
lane ℓ, n < 1.6·10⁷).
- **13 and 17 (≡ 1 mod 4): the level-2 recipe holds on every check (16 and 12).**
  - 17 decides the sign: the class is **w₂ = (2|ℓ)**, not "w₂ = −1 always".
  - At 17 the law is (y − 1)(y − 15), with 15 ≡ a₁₇ of the w₂ = +1 newform of weight 14.
- **19 and 23 (≡ 3 mod 4): every level-2 hypothesis fails**, for every sign and twist.
  - The data alone give (y − 1)(y − 13) at 19 (7 equations, unique).
  - **13 = a₁₉(ΔE₄) mod 19 (1563257180 ≡ 13)**: the cusp form of LEVEL 1 and weight 16. The level-2 newform of weight 16
    (a₁₉ ≡ 3) is absent, as w₂ = (2|19) = −1 picks the empty class.
- 23, 29 and 31: the predicted laws have degree 3–4 in y. That needs 7–9 levels, far beyond n < 1.6·10⁷.

**The rule that fits 5–19 (derived; the level-1 part was added after seeing 19, so 19 does not test it independently).**
- The constants are:
  - the Eisenstein 1;
  - the cusp forms of level 1 and weight ℓ − 3;
  - the newforms of level 2 and weight ℓ − 3 with w₂ = (2|ℓ);
  - all with no twist.
- Every prime fits: 5 and 7 (nothing), 11 (h₈ has w₂ = +1 ≠ (2|11): excluded), 13 (level-2 form, 7), 17 (level-2 form,
  15), 19 (ΔE₄, 13).
- Reading: 1/θ₄(z) = 1/θ(z + ½) lives on Γ0(4), and Niwa's correspondence matches S_(k+1/2)(Γ0(4)) with S_2k(Γ0(2)):
  level-1 forms and level-2 newforms together. So overpartitions follow the recipe, with these forms in place of
  level 6.

## 5 colours at 11: the law's family (`colors5_slope11`; STATEMENT G3's open entry)

- The 5-colour class at 11 is empty at weights 4 and 14. At 24, 34, 44 it holds exactly one form, of slope 1, with T₅,
  T₇, T₁₃ ≡ 10, 5, 10 on lane 11. At 54 a slope-3 form joins.
- Along the weights 104 + 110x (x = 0…10) the piece with those numbers keeps **exactly one root of slope 1 at every
  weight**. Newton's series reaches weight −6 with differences divisible by 11 to 0, 1, 3, 5, …, 19.
- **The constant (3|11)·u(−6)/11, digits lowest first: 2, 8, 5, 2, 10, 5, 2, 0, 1, 1, 3** (11 digits).
- **Against the exact 5-colour numbers divided by the law 11^(⌊b/2⌋+1)** (`colors_tower_r5_….json`; counted): the ratios
  H′_(b+2)/H′_b agree to 1–2, 2–4, 4–5, 4 and 4 digits at b = 1…5. At b = 4 and 5 that is every digit the data hold.
  **PREDICTED, and confirmed:** the 5-colour law at 11 follows one slope-1 family, with this constant.

## Double poles: E₄/j (a proven check) and 1/E₄² (`level1_tower` E4overj/invE4sq, `forms1_family` at w = −8, `double_pole_check`)

**Why.** 1/E₄ and 1/E₆ have simple poles. E₄/j = Δ/E₄² (weight 4) and 1/E₄² (weight −8) have a double pole at ρ.
- Zhang (arXiv 2510.23200, Thm 1.2) proved, for E₄/j: a_(npˡ) ≡ (−3|p)·p·a_(npˡ⁻¹) (mod p^(3l)) for every p ≥ 5.
- So E₄/j checks our tools on a proven case, and 1/E₄² is new.
- First values, checked by hand: E₄/j = q − 504q² + 180252q³ − …; 1/E₄² = 1 − 480q + 168480q² − ….

**The residue reasoning extended to poles of order m (derived).**
- T_ℓ multiplies the coefficient of (τ − ρ)^(−j) by π^(j−1)·π̄^(w−1−j) + π̄^(j−1)·π^(w−1−j). For j = 1 this is the simple-pole
  rule π^(w−2) + π̄^(w−2).
- The stabilizer of ρ allows only the orders with −j ≡ w (mod 3). For both forms that leaves only j = 2: the 1/(τ − ρ) term
  must vanish.
- At ℓ ≡ 1 mod 3 the integral root of the pair is π·π̄^(w−3) = **ℓ·u^(w−4)**, of slope 1.
- For E₄/j this gives ℓ·u⁰ = ℓ = (−3|ℓ)·ℓ, Zhang's constant at the split primes.

**PREDICTED before the run, for 1/E₄².**
- ℓ ≡ 1 mod 3: the law P(x)·(x − ℓ·u⁻¹²), with P the weight −8 families.
- ℓ ≡ 2 mod 3: P with possibly one constant of slope ≥ 1, which the data must name. Zhang's −p at weight 4 suggests one.
- The weight −8 families on lane ℓ: 5, 7, 13, 17, 19, 23 have only the Eisenstein 1; 11 has (x − 1)², since τ(11) ≡ 1
  (mod 11); 29 and 31 have two constants, 37 to 43 three or four.

**Measured** (`double_pole_check.txt`; n < 1.6·10⁷, counted).
- **E₄/j reproduces Zhang's theorem at every prime 5–43.** The data alone give one constant of slope 1:
  - +ℓ at ℓ ≡ 1 mod 3, as the residue formula ℓ·u^(w−4) = ℓ·u⁰ says and as Zhang says;
  - −ℓ at ℓ ≡ 2 mod 3, as Zhang's (−3|ℓ)·ℓ says.
  - The law (x − 1)(x − (−3|ℓ)·ℓ) holds to depth 5–12, most or all of each lane, from level 1.
  - The tools pass a proven case.
- **1/E₄²: the simple prediction P(x)·(x − ℓu⁻¹²) does not explain the data.**
  - The remainders of P's law do carry a slope-1 constant.
  - But its unit part, read as a single ratio, is not u⁻¹². At 7 it is 1, 3, 5 against 1, 5, 0; at 13, 19 and 31 even the
    first digit differs.
- **Refined, and supported.**
  - A slope-1 holomorphic family (non-ordinary, so not in P) shares the remainders with the pole.
  - Fitting a two-term law x² − s·x + p to the remainders over every level and n gives v(s) = 1 and v(p) = 2: two slope-1
    roots, at 7, 13 and 19 (46, 25 and 17 equations, all consistent).
  - **At 7, ℓ·u⁻¹² is a root to all 6 digits the data hold.** At 13 and 19 it is a root to one digit beyond chance
    (residual valuation 3 of 5 and 3 of 4).
- Reading: the double-pole constant ℓ·u^(w−4) is right at weight 4 (Zhang) and fits at weight −8 at 7. A sharper test at
  13 and 19 needs more levels (WSL), or the slope-1 families at weight −8 computed from the forms.

**The two slope-1 constants separated, and why 13 cannot be read yet** (`double_pole_fit_…txt`; WSL towers to n < 10⁸;
counted).
- The remainders of P's law fit a two-term law x² − s·x + p with v(s) = 1 and v(p) = 2. With x = ℓ·y its two roots are
  the unit parts of the two slope-1 constants.
- **At 7, one root is u⁻¹² to every digit the data hold (1, 5, 0, 4), from 57 equations.** The other is 4, 0, 1, 3
  (its 4th digit was wrong: the forms give 4, 0, 1, 6; see "What 13 needs").
  Since the first digit ≡ 1 is automatic (Fermat), chance would match the rest about 1 time in 343.
- At 13 the roots are 1, 9, 9 and 8, 0, 0, against u⁻¹² = 1, 0, 6: no match beyond the automatic first digit. At 19 the two
  roots agree on lane 19 and cannot be separated.
- **The reason, from the forms** (Newton polygons of T_ℓ on M_k(SL₂(ℤ)) at weights close to −8 ℓ-adically):
  - **at 7**, every such weight (34, 76, 118, 160) has exactly ONE holomorphic form of slope 1. So the remainders hold the
    pole and one family, and a two-term fit is exactly right.
  - **at 13**, weight 148 already has THREE forms of slope 1. With the pole that makes four slope-1 parts, which a
    two-term fit cannot separate. (At weights 304 and beyond the polygons need more than ℓ¹⁴ to read.)
- **Status of the double-pole rule ℓ·u^(w−4):**
  - weight 4 (E₄/j): every prime, and Zhang's theorem;
  - weight −8 (1/E₄²): confirmed at 7, the one prime where the data can separate it;
  - 13 and 19 need the holomorphic slope-1 families computed from the forms (Coleman families at weight −8), or many more
    levels.

## A pole at a point with no symmetry: 1/(E₄³ − 2Δ) (`level1_tower` poleJ0, `forms1_family` at w = −12, `pole_j0_check`)

**Why.** Every pole tested so far sat at a point with complex multiplication (ρ for 1/E₄ and 1/E₄², i for 1/E₆). There the
curve at the pole has an endomorphism of degree ℓ, which maps the pole to itself. A point with no complex multiplication has
none. Does its curve still reach the tower?
- F = 1/(E₄³ − 2Δ) = 1/(Δ·(j − 2)), weight −12, holomorphic at the cusp. Its one pole is at the point where j = 2.
- First values, checked by hand: 1 − 718q + 336196q² − 129585280q³ + ….
- The curve with j = 2: C: y² = x³ + 10356x + 11916304. Good reduction at every ℓ from 5 to 43 (its only bad prime above 3
  is 863). Counted: C is supersingular at 29 and 31 (aₗ = 0), ordinary at the other ten primes.

**The reasoning (derived, not proved).**
- ℓ-adically, U_ℓ sends a pole at a curve E to a pole at E/H, H its canonical subgroup (the lift of Frobenius, φ).
- At a CM point with ℓ split, φ fixes the point, and the residue is multiplied by u^(w−2) at each step: the simple-pole rule.
- At j = 2, φ moves the pole: level b has it at φᵇ(C). In Serre–Tate's coordinate φ is q ↦ q^ℓ, so φᵇ(C) converges, about
  one ℓ-adic digit per level, to the canonical lift of C mod ℓ. That lift is a CM curve with the same point count as C mod ℓ.
- So the residue factor converges to u_C^(w−2) = **u_C⁻¹⁴**, u_C the unit root of x² − aₗx + ℓ, aₗ = ℓ + 1 − #C(𝔽ₗ).
  It is the same for every twist of C (aₗ → ±aₗ).

**PREDICTED before the run.**
- At the ten primes where C is ordinary: the remainders of the weight −12 families' law P carry one more constant, a unit,
  equal to u_C⁻¹⁴. Its agreement may grow with the level (the pole is still moving), rather than hold to every digit at once.
- At 29 and 31 (C supersingular): P alone, deeper level by level, as at the inert primes of 1/E₄ and 1/E₆.
- If this holds, the tower of F reads the point count of a curve with no special symmetry: aₗ(C) comes out of the
  coefficients of F.

**Measured** (`level1_tower_poleJ02_N16000000_….txt/.json`, `pole_j0_check_poleJ02_….txt`; n < 1.6·10⁷; counted).
- **The prediction holds at every prime the data can read.**
  - At the seven ordinary primes with enough levels (5, 7, 11, 13, 17, 19, 23), the data's extra constant equals u_C⁻¹⁴.
  - Its agreement grows by one ℓ-adic digit per level, as the moving pole predicts. At 5: 1, 2, 3, 4, 5, 6, 7 digits at
    levels 1→2 … 7→8. At 7: 1 to 6 digits.
  - Of all ±u_C^m with |m| ≤ 60, m = −14 fits best at every one of these primes.
  - The law P(x)·(x − u_C⁻¹⁴) holds to depth b + 1 from level b at every one of them (at 5: 2, 3, …, 9).
  - At 29 and 31, where C is supersingular, P's law alone gets one digit deeper per level (depth 1, 2, 3): no extra constant.
  - 37, 41 and 43 need more levels than n < 1.6·10⁷ gives (P has degree 3 there).
- **So the curve with j = 2, which has no complex multiplication, reaches the tower through its canonical lift.** The rule
  "a simple pole adds u^(w−2)" holds with u the unit root of the curve at the pole, whatever that curve is. The constant is
  reached one digit per level, not at once.
- **The blind reading: the tower finds the curve.** j₀ was not used for this step.
  - From the coefficients alone, at each prime: the traces whose u⁻¹⁴ equals the constant read (or aₗ ≡ 0 where the
    remainders die), then the residues j mod ℓ of curves with that trace.
  - Over the nine readable primes, 1728 residues mod 33426748355 remain. The integers |j| ≤ 10⁶ that fit every prime:
    **2** and −575358. About 0.1 chance fits were expected.
  - −575358 is ruled out by any of 37, 41, 43, 47, 53: its traces there (±6, ±4, ±1, ±8, ±4) differ from C's (±8, ±5, ±10,
    ±2, ±8). The tower would show this once those primes have enough levels. That is derived, not measured.
- Fixed after the first check run, before any j was read (recorded in the script header):
  - the constant is read only from unit ratios (the first run read the dying remainders at 31 as a constant 0);
  - the blind reading uses the deepest read that at least 3 values of n confirm (the first used a single n at 5, whose
    last 2 of 11 digits were not yet settled).

## What 13 needs: the slope-1 families at weight −8, computed directly (`slope1_families`, `double_pole_full_law`)

**Why.** 1/E₄²'s tower at 13 holds four constants of slope 1: the double pole's 13·u⁻¹² and three slope-1 holomorphic
families. A two-term fit cannot separate four. The way out is to compute the three families' constants at weight −8 from
the forms, so the whole law is fixed in advance and the tower only tests it.

**The tool** (`slope1_families.py`; Lauder's algorithm, LMS J. Comput. Math. 14 (2011), Algorithm 1).
- U_ℓ on overconvergent forms of level 1 at any weight, negative ones included, through Katz's expansions. No interpolation
  from far-off weights.
- Its characteristic series modulo ℓ^N. The slope-0 and slope-1 factors come out by Hensel's lemma, the cut of the series set
  by a bound.
- **Checked against the classical forms** at weight 148 for 13 (T₁₃ on M₁₄₈(SL₂(ℤ)), the whole polynomial, no cut). Same
  Newton polygon (0, 1, 1, 1, 3, 4). The slope-1 cubic agrees to all 7 digits claimed at ℓ¹², and all 13 at ℓ²⁰.
  - A first version cut the series too early. A check that cut both sides at the same place agreed to 8 digits but was
    wrong from the 6th. The untruncated check caught it.
  - Dwork's bound, which the method needs, holds on every matrix (smallest margin −0.71 ≥ −1).

**At weight −8** (`slope1_families_{7,13,19}_w-8_N….txt`; counted):
- **7:** slope 0: exactly 1 (the Eisenstein family). Slope 1: one family, unit part y = 4, 0, 1, 6, 1, 2, 0, 0, 2, 5, 6, 3, 0,
  4, 1, 3. Next slope 3.
- **13:** slope 0: exactly 1. Slope 1: **three families**. Their cubic in y is known to 17 digits: one root in ℤ₁₃ (11, 10,
  10, 7, 3, 12, …) and a conjugate pair. Next slopes 3, 4, 5, 6.
- **19:** slope 0: exactly 1. Slope 1: two families (y = 1, 13, 15, … and 7, 15, 5, …). Next slope 2.
- At 7 the earlier data fit (README "Double poles") gave the family's unit part as 4, 0, 1, **3**. The forms give 4, 0, 1,
  **6**. They differ at the 4th digit.

**PREDICTED before the run** (`double_pole_full_law.py`). The law
  P(x)·(x − ℓ·u⁻¹²)·S(x),  S(x) = ∏(x − ℓ·y) over the slope-1 families at weight −8 computed above,
has every constant fixed by the forms and the pole, with no tower data used. It holds on 1/E₄²'s tower at 7, 13 and 19, to a
depth that grows with the level: about 3 digits per level at 7 and 13 (next slope 3), 2 at 19 (next slope 2).
- Dropping either factor (the pole's, or the families') makes it fail.
- At 7 the law decides between 4, 0, 1, 6 (forms) and 4, 0, 1, 3 (fit).

**Measured** (`double_pole_full_law_invE4sq_N100000000_7_13_19.txt`; WSL towers to n < 10⁸; counted). The smallest depth
over n ≥ 1, by starting level:

| ℓ (lane) | P alone | pole only | families only | **pole and families (predicted)** |
|---|---|---|---|---|
| 7 (7¹⁰) | 1, 2, 3, … | 3, 4, 5, … | 2, 3, 4, … | **7, 10, 10, 10, 10, 10** (57 checks) |
| 13 (13⁸) | 1, 2, 3, … | 2, 3, 4, … | 4, 5, 7 | **8, 8** (12 checks) |
| 19 (19⁷) | 1, 2, 3, … | 2, 3, 4, … | 3, 4, 5 | **6, 7** (13 checks) |

- **The prediction holds at all three primes.** With every constant fixed in advance, the law reaches the whole lane: from
  level 1 at 13, from level 2 at 7 and 19. Each factor is needed: without the pole, or without the families, the depth
  stays several digits short.
- **What 13 needed:** four slope-1 constants. They are the pole's 13·u⁻¹² and the three slope-1 families' 13·y at weight
  −8, which the forms give (one in ℤ₁₃, a conjugate pair). With them, 13 agrees to every digit the data hold.
- **So the double-pole rule ℓ·u^(w−4) holds at 7, 13 and 19** at weight −8, and (Zhang) at weight 4 for every prime.
- **A correction.** At 7 the earlier two-term fit's family root, 4, 0, 1, 3, was wrong in its 4th digit. With it the law
  stops at 8, 9 digits where the forms' 4, 0, 1, 6 reaches 10. The fit's pole root (1, 5, 0, 4) was right. The slope-3
  part at low levels most likely pulled the fit's last digit.

## Curves that repeat themselves (`level1_tower` poleJ0/poleH, `slope1_families` at w = −12, `forms1_family` at w = −24, `repeat_check`)

**Why** (12 Nov).
- U_ℓ moves a pole by φ, the lift of Frobenius (a curve E goes to E/H, H its canonical subgroup).
- The curves φ sends back to themselves are the canonical lifts. Those are curves with complex multiplication (CM) in which
  ℓ splits.
- A CM curve whose prime above ℓ has order r in its class group comes back after r steps. The other CM curves of that
  discriminant go round the same cycle.
- A curve without CM never repeats. Its pole drifts toward the canonical lift of its reduction mod ℓ.

**Facts used** (counted from point counts and the CM j-values; derived where marked).
- **Where the pole of 1/(E₄³ − 2Δ) drifts to,** at the primes where that lift has class number 1. With s = v_ℓ(2 − j_lift),
  the drift's starting distance in digits:
  - 5: j = −32768 (CM by √−11), s = 1;
  - 7: j = 54000 (CM by √−3, conductor 2), s = 2;
  - 11: j = −3375 (CM by √−7), s = 1;
  - 19: j = 54000, s = 1;
  - 43: j = 8000 (√−2), s = 1.
  At 13, 17, 23, 37 and 41 the lift has class number 2 or more.
- **Discriminant −15** (class number 2; H(x) = x² + 191025x − 121287375, checked numerically):
  - 17, 23 (and 47, 53): the two CM curves **swap at every step**;
  - 19, 31: each curve maps to itself;
  - 7, 11, 13, 29, 37, 41, 43: supersingular; 5 ramified.

**PREDICTED before the runs.**
- **A. A pole on a repeating curve against a drifting pole.**
  - Put the pole on the lift itself: 1/(E₄³ − j_lift·Δ) at its prime, j_lift = −32768 at 5, 54000 at 7 and 19, −3375 at 11.
    Then the constant u⁻¹⁴ is exact from the first level.
  - The law P(x)·(x − u⁻¹⁴)·S(x) should then gain at least 2 digits per level, the next slope's rate. S is the slope-1
    families at weight −12, computed directly.
  - For j = 2 the same law should gain about one digit per level, from a start near s: the drift.
  - This also tests the earlier reading of j = 2's one-digit-per-level growth, which slope-1 families could have caused
    instead. Including S separates the two.
- **B. Two poles at once: 1/(Δ²·H₋₁₅(j)) = 1/(E₄⁶ + 191025E₄³Δ − 121287375Δ²)**, weight −24, so u^(w−2) = u⁻²⁶.
  - **7, 11, 13, 29, 37, 41, 43 and 5:** the forms' P alone.
  - **19 and 31 (each curve fixed):** one extra constant u⁻²⁶, with u from the trace ±4 at 19 and ±8 at 31. The law is
    P(x)·(x − u⁻²⁶), once even though there are two poles, since both curves have the same Frobenius up to sign.
  - **17 and 23 (the curves swap):** the factor after two steps is U⁻²⁶, with U the unit root of x² − t·x + ℓ² and t = ±14
    (17), ±34 (23). That is the trace of Frobenius over 𝔽_ℓ², the endomorphism of degree ℓ² that the 2-cycle builds.
    - So the law is P(x)·(x² − U⁻²⁶): a pair of constants +U⁻¹³ and −U⁻¹³. Neither one alone will do.
    - The remainders repeat every two levels: the one-level ratio is not constant, the two-level ratio is.
  - The cycle length r would then be read off the law's shape, x^r − c. It is the order of ℓ's prime in the class group.

**Measured** (`repeat_check_A.txt`, `repeat_check_B.txt`; towers `level1_tower_poleJ0-32768_…`, `…poleJ054000_…`,
`…poleJ0-3375_…`, `…poleHD15_…` to n < 1.6·10⁷; families `slope1_families_{5,7,11,19}_w-12_…`,
`slope1_families_{17,19,23,31}_w-24_…`, `forms1_family_w-24_ℓ.txt`; counted).

**A. The repeating curve against the drifting one.** Law P(x)·(x − u⁻¹⁴)·S(x), smallest depth by starting level:

| ℓ (lane) | slope-1 families | pole ON the lift (repeats) | pole at j = 2 (drifts), s |
|---|---|---|---|
| 5 (5¹¹) | none | **5, 7, 9, 11, 11, …** | 2, 3, 4, 5, 6, 7, 8, 9 (s = 1) |
| 7 (7¹⁰) | 1 | **6, 8, 10, 10, 10** | 4, 5, 6, 7, 8 (s = 2) |
| 11 (11⁸) | 1 | **5, 7, 8** | 3, 4, 5 (s = 1) |
| 19 (19⁷) | 2 | **6** | 4 (s = 1) |

- **The prediction holds at all four primes.**
  - A pole on the repeating curve gives its constant exactly: two digits per level (the next slope is 2), the whole lane
    within three or four levels.
  - The drifting pole gains one digit per level, and starts deeper where it starts closer (7, s = 2).
- **A correction to "A pole at a point with no symmetry."** There the law without S grew one digit per level at every prime,
  and I credited the drift. Without S, the repeating pole grows the same way at 7, 11 and 19 (b + 1), so a slope-1 family
  alone causes it there. The drift is real: with S included, the drifting pole falls one digit per level behind at all four
  primes. But only at 5, which has no slope-1 family at weight −12, was the earlier observation the drift alone.

**B. Two poles: 1/(Δ²·H₋₁₅(j)), weight −24.**
- **Supersingular primes** (7, 11, 13, 29, 37, 41, 43; 5 ramified): P alone, deeper by one digit a level, at every one.
- **19 and 31 (each curve fixed).** Both roots of H give the same trace (−4, −4 at 19; 8, 8 at 31), so there is one
  constant u⁻²⁶.
  - P alone: depth 0.
  - P·(x − u⁻²⁶): 2, 3.
  - With S at 19: 5.
- **17 and 23 (the curves swap).**
  - The curve over 𝔽_ℓ² has trace −14 (17) and 34 (23). Counted on the curve, and equal to the norm-form prediction ±14,
    ±34.
  - **The remainders repeat every two levels.** One level apart their ratio is scattered: 10 different values on lane 17 at
    each level. Two levels apart it is one value at every level and every n: 15 at 17, 16 at 23. That is U⁻²⁶ on the lane.
  - **Only the pair works.** P·(x − U⁻¹³)·S and P·(x + U⁻¹³)·S: depth 0. P·(x² − U⁻²⁶): 2, 3 at 17, and 2 at 23 (the
    levels run out; with S the law is too long for them).
- So the tower reads **how often the curve at the pole repeats**: every step at 19 and 31, every second step at 17 and 23.
  That cycle length is the order of ℓ's prime in the class group of discriminant −15, read off the coefficients.

**What this answers.** The curves that "end up repeating themselves" under the step U_ℓ makes are the CM curves where ℓ
splits; every other curve drifts toward one of them.
- The tower tells which kind sits at its pole: an exact constant from the start, or one digit per level.
- For a curve without CM it names the repeating curve it drifts to (here −32768, 54000, −3375).
- For CM curves it gives the cycle length (x − c for 1, x² − c for 2).

**C. A 3-cycle: 1/(Δ³·H₋₂₃(j)), weight −36** (H₋₂₃ = x³ + 3491750x² − 5151296875x + 12771880859375, checked numerically;
class number 3).
- Counted: at 13, 29, 31, 41 and 47, H is irreducible mod ℓ and ℓ's prime has order 3, so **the three CM curves cycle with
  period 3**. At 59 each is fixed. At 5, 7, 11, 17, 19, 37, 43 and 53 they are supersingular; 23 is ramified.
- The traces over 𝔽_ℓ³ (norm form x² + 23y² = 4ℓ³): ±74 at 13, ±282 at 29, ±344 at 31, ±426 at 41.
- **PREDICTED before the run:**
  - At the cycling primes the law is P(x)·(x³ − U⁻³⁸)·S(x), U the unit root of x² − t·x + ℓ³.
  - The remainders' ratio is constant only three levels apart (U⁻³⁸). One or two levels apart it scatters.
  - No single cube root of U⁻³⁸, as a factor x − c, works alone.
  - At the supersingular primes: P alone.
  - The towers to n < 1.6·10⁷ reach this at 13. 29 and 31 need more levels (WSL, n < 10⁸).

**Measured, first look at C** (`repeat_check_C_N16000000.txt`; towers to n < 1.6·10⁷; counted).
- **13 (the three curves cycle):**
  - P alone: depth 0.
  - The remainders' ratio scatters one level apart (5, 7, 3 values on lane 13 at levels 1, 2, 3) and two levels apart (7,
    then 2 values).
  - **Three levels apart it is one value, 9, which is U⁻³⁸ on the lane.**
  - P·(x³ − U⁻³⁸) holds to depth 2 on the one level span the data allow.
  - 9 is not a cube mod 13, so x³ − U⁻³⁸ has no root in ℤ₁₃: the three constants are not 13-adic numbers, and only the
    cubic as a whole can appear in the law.
- **59 (each curve fixed):** one ratio, 48 = U⁻³⁸, and P·(x − U⁻³⁸) holds (depth 2, one check).
- **41 (cycle):** a single cube root as a factor fails (depth 0). Three levels are beyond these data.
- **29 and 31:** too few levels here; the WSL towers to 10⁸ are running.
- **The supersingular primes** (5, 7, 11, 17, 19, 43; 23 ramified): P alone, one or two digits deeper per level.

**Measured at 10⁸** (`repeat_check_C.txt`; WSL towers `level1_tower_poleHD23_N100000000_13_29_31`; counted).
- **13:**
  - Three levels apart the ratio is **9 = U⁻³⁸ at both levels 1 and 2**.
  - One level apart it scatters (5, 7, 7 values at levels 1–3), and two levels apart (7, then 5 values).
  - P·(x³ − U⁻³⁸): depth 2, 3 (12 checks). With the slope-1 family, P·(x³ − U⁻³⁸)·S: depth 5 (1 check, all that level 7
    holds).
  - P alone and P·S: depth 0 everywhere.
- **29:** the ratio is not constant one level apart (9 values) or two levels apart (4 values). So the period is neither 1
  nor 2. Three levels need n beyond 29⁶ ≈ 6·10⁸.
- **31:** the ratio one level apart scatters (3 values); two and three levels apart are beyond 10⁸.
- **So the tower sees a 3-cycle where the class group says one**, at 13 directly, and at 29 as far as the data reach.
  Together with B: the cycle length of the curves at the poles, 1, 2 or 3, shows in the tower as the step at which the
  remainders' ratio settles.

**B at 10⁸** (`repeat_check_B_N100000000.txt`; WSL towers `level1_tower_poleHD15_N100000000_17_19`; counted).
- **17 (the curves swap):**
  - The remainders' ratio two levels apart is **15 = U⁻²⁶ at all three levels it can be read**, for every n. One level apart
    it scatters at every level (10, 10, 10, 4 values).
  - P·(x² − U⁻²⁶): 2, 3, 4 (26 checks).
  - **With the families, P·(x² − U⁻²⁶)·S: depth 6 of 7 (4 checks).**
  - Either half with S, P·(x ∓ U⁻¹³)·S: depth 0 at both levels (15 checks each).
- **19 (each curve fixed):** P·(x − u⁻²⁶)·S: **5, 7**, the whole lane from level 2 (13 checks).
- The prediction holds in full where the data reach: the pair, and only the pair, at 17; one constant at 19.

**Why FLINT or Sage cannot simply go deeper at 29 and 31.**
- FLINT already does every computation here. The wall is length, not speed: three levels at 29 need every coefficient up to
  n ≈ 29⁶ ≈ 6·10⁸.
- At that length one series is about 5 GB, and the computation holds several at once. That is more than WSL's 31 GB.
- Sage does this arithmetic with the same FLINT routines, so it would not help with this. (Its Lauder implementation could
  cross-check `slope1_families` some day.)
- The way around it is mathematical: choose a discriminant whose 3-cycle falls at a small prime.

**D. A 3-cycle at small primes: 1/(Δ³·H₋₃₁(j)), weight −36.**
- H₋₃₁ = x³ + 39491307x² − 58682638134x + 1566028350940383. Computed from the three CM values to 60 digits; rounding error
  below 10⁻⁴⁵. Class number 3.
- Counted:
  - the three curves cycle with period 3 at **5, 7, 19 and 41** (H irreducible mod ℓ, ℓ = 2x² ± xy + 4y²);
  - traces over 𝔽_ℓ³: ±2, ±16, ±156, ±278;
  - supersingular at 11, 13, 17, 23, 29, 37, 43; 31 ramified.
- **PREDICTED before the run:**
  - At 5, 7 and 19 the remainders' ratio settles only three levels apart, at U⁻³⁸ (U the unit root of x² − t·x + ℓ³), at
    every level the data reach. One and two levels apart it scatters.
  - P(x)·(x³ − U⁻³⁸)·S(x) holds and deepens level by level. No single cube root, as a factor x − r, does.
  - At 5 (ℓ ≡ 2 mod 3) the cubic has one root in ℤ₅ and a conjugate pair. Even the ℤ₅ root alone fails.
  - At the supersingular primes: P alone.

**Measured, D** (`repeat_check_C_D31_first.txt`; towers `level1_tower_poleHD31_N16000000_5_7_19`; weight −36 families
`slope1_families_{5,7,19}_w-36_…`; counted).
- **5 (the curves cycle):**
  - The remainders' ratio three levels apart is **4 = U⁻³⁸ at all six levels it can be read**, for every n. One level apart
    it scatters (3 or 4 values at every level), and two levels apart (3 or 4 values).
  - **P·(x³ − U⁻³⁸)·S: 9, 11, 11, 11, 11, the whole lane 5¹¹ from level 2 (42 checks).** Without S: 4, 5, …, 9, one digit
    per level (the slope-1 family).
  - The ℤ₅ cube root alone, P·(x − r)·S: depth 0 at all seven levels (64 checks). P alone and P·S: 0.
- **7:**
  - Three levels apart the ratio is **2 = U⁻³⁸ at all four levels**; one and two levels apart it scatters.
  - **P·(x³ − U⁻³⁸): 4, 6, 8, 10, two digits per level to the whole lane 7¹⁰** (35 checks). 7 has no slope-1 family at
    weight −36; the next slope is 2.
  - 2 is not a cube mod 7, so the cubic has no root in ℤ₇: the three constants live only in an extension.
- **19:** the ratio scatters one level apart (6, 5 values) and two levels apart (5 values), so the period is not 1 or 2.
  Three levels need n < 10⁸ (WSL, running).
- **So the 3-cycle is confirmed in full at 5 and 7.** The ratio settles exactly at the third level, at the predicted value,
  at every level. The cubic law reaches the whole lane, and no single constant can replace it.
- **Completed** (`repeat_check_C_D31.txt`; with the WSL tower `level1_tower_poleHD31_N100000000_19` and the controls):
  - **19:** three levels apart the ratio is **6 = U⁻³⁸**. One level apart it scatters (6, 7, 2 values), and two levels
    apart (9, 2). P·(x³ − U⁻³⁸): depth 2 (2 checks; the levels run out for the law with S).
  - **41:** a single cube root as a factor fails (depth 0, 16 checks). Three levels are beyond n < 1.6·10⁷.
  - **Controls:** P alone holds and deepens at every supersingular prime (11, 13, 17, 23, 29, 43) and at the ramified 31.
    The prediction holds everywhere the data reach.

**Where curves return to themselves, predicted in advance** (the user's question, 12 Nov; derived from the class group, counted
by the script in the session log).
- **For given curves, which primes.** The three curves of discriminant −31 (ℓ < 200):
  - return every step at 47, 67, 131, 149, 173 (ℓ = x² + xy + 8y²);
  - every third step at 5, 7, 19, 41, 59, 71, 97, 101, 103, 107, 109, 113, 157, 163, 191, 193;
  - never at the other 22 primes (supersingular).
  - In the long run (Chebotarev): half the primes never, a sixth every step, a third every third step.
- **For a given prime, which curves.** The curves that return every step at ℓ are the canonical lifts, one for each ordinary
  j mod ℓ. They are the CM curves whose discriminant D solves 4ℓ = t² − D·s² with ℓ ∤ t.
  - At 5 there are exactly four: j = 287496 (D = −16), −32768 (−11), 1728 (−4), −884736 (−19). Every other pole drifts to
    the one with the same j mod 5, as j = 2 drifted to −32768.
  - At 7 there are six: j = 0, −884736, 54000, −12288000, and the two curves of discriminant −24.
- A returning curve need not have a whole-number j: the j = 2 pole drifted at 13, 17 and 23 to lifts of class number 2, 3
  and 3, and the rule held there.
- This is classical (Deuring's lifting, Serre–Tate's canonical lifts). What is new here is that a tower's coefficients
  show it.

## What the tower view reveals: the distance from a curve to its returning twin (`level1_tower` poleJ0 at 5, `drift_distance_check`)

**Why** (12 Nov).
- A pole at a curve without CM drifts toward the curve's returning twin, its canonical lift.
- In Serre–Tate's coordinate t the step is t ↦ (1 + t)^ℓ − 1, which adds exactly one ℓ-adic digit of closeness per level.
  So a curve that starts s digits from its twin is s + b digits from it at level b (derived).
- At 5, j = 2 starts s = v₅(2 − (−32768)) = 1 digit away, and its law depth was b + 1. The pole placed on the twin (−32768)
  gave 2b + 3, the next slope's rate.

**PREDICTED before the runs.**
- For 1/(E₄³ − j₀Δ) at 5 with j₀ = −32768 + 5^s (s = 2, 3, 4, 6), all reducing to the same curve as j = 2: the law
  P(x)·(x − u⁻¹⁴) holds to depth min(2b + 3, b + s, 11) from level b.
- So a curve s digits from its twin behaves like the twin itself for about s − 3 levels, then falls behind one digit a level.
  - s = 2: 3, 4, 5, …;
  - s = 3: 4, 5, 6, …;
  - s = 4: 5, 6, 7, …;
  - s = 6: 5, 7, 9, 10, 11.
- If this holds, a tower's coefficients measure how close a curve is to the curve that repeats, a quantity no single
  coefficient shows.

**Measured** (`drift_distance_check.txt`; towers `level1_tower_poleJ0{-32743,-32643,-32143,-17143}_N16000000_5`; counted).

| j₀ | digits from the twin, s | depth from levels 1–8 |
|---|---|---|
| −32768 (the twin itself) | — | 5, 7, 9, 11, 11, 11, 11, 11 |
| 2 | 1 | 2, 3, 4, 5, 6, 7, 8, 9 |
| −32743 | 2 | 3, 4, 5, 6, 7, 8, 9, 10 |
| −32643 | 3 | 4, 5, 6, 7, 8, 9, 10, 11 |
| −32143 | 4 | 5, 6, 7, 8, 9, 10, 11, 11 |
| −17143 | 6 | 5, 7, 9, 10, 11, 11, 11, 11 |

- **Every level of every run is exactly as predicted**: 48 of 48 levels, 75 checks each.
- A curve s digits from its twin starts s digits deep. At s = 6 the tower is indistinguishable from the twin's for three
  levels, then falls behind one digit a level.
- **So a tower's coefficients measure the ℓ-adic distance between the curve at its pole and the curve that repeats**: the
  valuation of its Serre–Tate coordinate. No single coefficient shows this.

## Distinct parts: the constants from forms (`forms_eta2`, `distinct_parts_check`)

**Why.** Distinct parts, ∏(1 + qⁿ) = q^(−1/24)·η(2z)/η(z) (weight 0), had laws found from the data alone (README "How wide
the recipe reaches"):
- the tower vanishes on lane 5;
- H_(b+2) ≡ −H_b at 7;
- x(x² + 2) at 11.
Their constants had not come from forms. The level-1152 build that seemed necessary is not.

**The forms** (`forms_eta2.py`; derived, and checked).
- A form of weight k with the multiplier v^r of (η(2z)/η(z))^r is (η(2z)/η(z))^r·g, with g ∈ M_k(Γ0(2)) vanishing at the
  cusp 0. That is a space of dimension ⌊k/4⌋.
- T_ℓ carries class r to class rℓ mod 24: (T_ℓF)(n) = A(ℓn) + (2|ℓ)ℓ^(k−1)A(n/ℓ) on the forms scaled by z → 24z, which lie in
  M_k(Γ0(1152), (2|·)).
- So one tower level changes class, and two levels return. The two-level step is the composite C = T_ℓ∘T_ℓ.
- On the ordinary part, U_ℓ² is the unit root X of X² − (C − 2c)X + c² = 0, with c = (2|ℓ)ℓ^(k−1). Derived from
  U f = T f − c·V f between the two classes.
- Checked: at every ℓ from 5 to 23 the images of T_ℓ lie in the target space, on 8 coefficients beyond those used.

**On lane ℓ, at weight ℓ − 1** (counted). The composite's characteristic polynomial in y (two levels):
- 5: y. 7: y + 1. 11: y(y + 2).
- **These are exactly the three laws the data had found alone:** dying at 5, −1 at 7, −2 and a dying part at 11.
- 13: (y + 3)(y + 5)(y + 6).
- 17, 19: irreducible quartics.
- 23: (y − 1)·(a quartic).

**PREDICTED before the tower runs at 13 and beyond, and before the deep digits:**
- The towers obey, in two-level steps, χ(y), the unit part of that polynomial carried by Hida's families to weight 0.
- On lane ℓ that is the weight ℓ − 1 polynomial above.
- Every digit is the families' value at weight 0, from the weights ℓ − 1, 2(ℓ − 1), ….
- At 7 and 11 the law holds to a depth that grows with the level. At 13 it holds on lane 13 wherever the levels allow.

**Measured, first run** (`distinct_parts_check_N16000000_5_7_11_13_17_19_23.txt`; counted). In two-level steps; depth by
starting level:
- **5:** the composite is 0 on lane 5, so there is no ordinary part and the law is "the tower dies". Depth 0, 0, 1, 1, 2, 2,
  3, 3, 4, 5: one digit every two levels.
  - As I recall, that is the rate of the known congruences for distinct parts mod powers of 5 (Rødseth; Gordon–Hughes); not
    re-read here.
- **7:** on lane 7, H_(b+2) ≡ −H_b at levels 1–6 (57 checks).
  - With every digit of the weight-0 constant (Hida, 13 weights; differences to 7¹⁰): depth 1, 2, 3, 4, 5, 6, one digit a
    level.
  - So the constant from forms matches the tower to 6 digits, as deep as the tower's own convergence allows. The rate is set
    by the next family, which these forms do not give.
- **11:** on lane 11, χ(y) = y + 2 holds from level 3 (depth 0, 0, 1, 1). y·χ(y) holds from level 1, which is the data's law
  x(x² + 2) again.
  - The deep digits are not confirmed here: the dying part leaves depth 1 at levels 3 and 4.
- **13, 17, 19, 23:** the predicted laws have degree 3 to 5 in y and need 7 to 11 levels. The towers to n < 1.6·10⁷ have
  6 or fewer.
- **Every lane-ℓ polynomial is recorded above as a prediction.** At 13: χ(y) ≡ (y + 3)(y + 5)(y + 6); every digit to 13⁸
  in the output.

## What a tower's coefficients reveal (a summary, 12 Nov; each line points to its measured section)

A tower applies one step, U_ℓ, again and again. Its coefficients therefore record what that step does: its eigenvalues
and slopes on the forms the start touches, and its motion on the curves at the poles. Read that way, a list of integers
shows things no single coefficient does. Each item below was predicted before its run and then counted.

1. **The point count of a curve the form never names.** A pole at j₀ gives the unit root u of the curve with j = j₀ at
   every ordinary prime. Read blind, nine primes give back j₀ = 2 (A pole at a point with no symmetry).
2. **Whether that curve repeats.** On a CM curve (a canonical lift) the constant is exact: depth 2 digits a level, the
   whole lane in 3–4 levels. Elsewhere the pole drifts and the tower gains 1 digit a level (Curves that repeat themselves,
   A).
3. **How far the curve is from the one that repeats.** At 5 the depth is min(2b + 3, b + s, 11), with s the 5-adic distance
   to the twin; 48 of 48 levels as predicted (What the tower view reveals).
4. **How the CM curves move among themselves.** Their cycle length under Frobenius, 1, 2 or 3, is the level at which the
   ratios settle. It is the order of ℓ's prime in the class group (B, C, D; −15 at 17 and 23, −23 at 13, −31 at 5, 7, 19).
5. **Families no classical weight shows.** The U_ℓ-values of overconvergent families at negative weights, e.g. 13's four
   slope-1 constants for 1/E₄², confirmed to the whole lane (What 13 needs).
6. **The slopes, from the speed.** How fast the depth grows is the next slope of U_ℓ: 2 digits a level at 5 for weight
   −12, 3 for weight −8 at 7 and 13, one digit every two levels for distinct parts at 5.
7. **For p(n) itself** (Dig 1): the constants are a weight −2 family's values (Shimura's partner of weight −1/2); every form
   behind them is ordinary; 607 is a repeated root.

In one sentence: a tower's coefficients carry the spectrum of U_ℓ on the holomorphic side, and Frobenius's motion on the
moduli of elliptic curves at the poles. That motion has fixed points at the CM curves, cycles from the class groups, and a
pull of one digit per step toward the fixed points.

**Measured to 10⁸** (`distinct_parts_check_N100000000_7_11_13.txt`; WSL towers; counted).
- **7:** with every digit of the weight-0 constant, depth 1, 2, 3, 4, 5, 6, 7 at levels 1–7 (68 checks). The constant from
  forms matches the tower to 7 digits, still one digit a level. On lane 7 the law holds at every level.
- **11:** the deep ordinary law gives depth 0, 0, 1, 1, 2 at levels 1–5. That is one digit every two levels, as a dying part
  of slope 1 per two levels would allow. With the dying factor, y·χ(y): 1, 1, 2.
- **13 (out of sample):** the predicted law (y + 3)(y + 5)(y + 6) holds on lane 13 on the one check level 7 allows. With every
  digit it holds to depth 2. Chance alone would pass that about 1 time in 169. A real test needs n near 13⁸ ≈ 8·10⁸.

## More moonshine: T_2A and T_3C (`moonshine_check`, `forms_eta3`)

**Why.** T_2B's constants came from forms (README "T_2B against the forms"). Two more McKay–Thompson series fit the engines
we have:
- **T_2A** = (η(τ)/η(2τ))²⁴ + 24 + 4096·(η(2τ)/η(τ))²⁴ = q⁻¹ + 4372q + 96256q² + … (checked by hand), the Hauptmodul of
  Γ0(2)+;
- **T_3C** = j(3τ)^(1/3) = E₄(3τ)/η(3τ)⁸ = q⁻¹ + 248q² + 4124q⁵ + …, E₈'s series.

**The forms** (derived).
- **T_2A** is fixed by the Fricke involution W₂. At weight 0 the Eisenstein piece that gave T_2B its constant 1 (the limit of
  (E_k(z) − E_k(2z))/c_k) changes sign under W₂, so T_2A's tower cannot use it. Its constants are then:
  - the level-1 cusp forms (through the W₂-symmetric oldforms);
  - the level-2 newforms of one W₂ sign. Which sign, fixed or carried by (2|ℓ), the data decide.
- **T_3C** has the multiplier of η⁻⁸ = that of η¹⁶. Forms with it are η¹⁶·M_(k−8)(SL₂(ℤ)), and with the other class η⁸·M_(k−4).
  - None has a constant term, so there is no Eisenstein constant.
  - T_ℓ keeps the class when ℓ ≡ 1 mod 3 and swaps the two when ℓ ≡ 2 mod 3. That gives one-level laws for the first, and
    two-level laws (T_ℓ∘T_ℓ) for the second, as for distinct parts.
  - dim η¹⁶·M_(ℓ−9) = 0 for ℓ = 5, 7, 11.

**PREDICTED before the run.**
- **Lehner-type dying towers:** T_2A at 3, 5 and 7, and T_3C at 5, 7 and 11. No law there but "the tower dies" (no
  constant 1, unlike T_2B at 5 and 7).
- Above those primes, the lane-ℓ laws are the forms' ordinary eigenvalues as described: for T_2A one W₂ class, the same rule
  at every prime; for T_3C, T_ℓ or T_ℓ∘T_ℓ by ℓ mod 3.

**Measured** (`moonshine_check_N16000000_3_5_…_43.txt`; towers to n < 1.6·10⁷; first values checked: T_2A 1, 0, 4372, 96256,
1240002; E₄/η⁸ 1, 248, 4124, 34752, 213126; counted).

**T_2A:**
- **It dies at 3, 5 and 7**: every value from level 1 on is ≡ 0 on lane ℓ. That is the Lehner-type prediction: no Eisenstein
  constant, unlike T_2B's 1 at 5 and 7.
- From 11 on, the constants are the level-1 cusp forms and the level-2 newforms with **w₂ = (2|ℓ)**, at every prime the data
  test:

  | ℓ | law (forms) | data alone |
  |---|---|---|
  | 11 | x + 4 | x + 4 |
  | 13 | x + 5 (Δ's family) | x + 5 |
  | 17 | (x + 10)(x + 1) | same |
  | 19 | x² + 16x + 1 | same |
  | 23 | x² + 19x + 18 | same |
  | 29 | the cubic [12, 23, 20, 1] | same |
  | 31 | the cubic [28, 6, 20, 1] | same |

- "w₂ = +1 always" fails at 19; "w₂ = −1 always" fails at 17. Only (2|ℓ) fits every prime, the rule overpartitions also
  followed.
- At 37, 41 and 43 the laws (degree 4, 5, 5) are recorded as predictions; the levels run short.

**T_3C:**
- **It dies at 5, 7 and 11**, as predicted (from level 2 at 5 and 11, where level 1 still sits in the other class, η⁸·M₀).
- From 13 on, the forms η¹⁶·M_(ℓ−9) give every law the data found alone:
  - one-level steps for ℓ ≡ 1 mod 3: 13: x + 8; 19: x + 10; 31: x² + 4x + 23; 37: x³ + 11x + 34; 43: x³ + 12x² + 33x + 35;
  - two-level steps for ℓ ≡ 2 mod 3: 17: x² + 15; 23: x² + 22.
- 29 and 41: the two-level laws are recorded as predictions; the levels run short.
- In every space computed, the Hecke images lie in the target space (the engine's own check).

**Not checked against the literature:** whether these Lehner-type congruences for T_2A and for j^(1/3) are in print. Lehner's
for j are; congruences for McKay–Thompson series have been studied.

## A mock theta function: Ramanujan's f(q) (`mock_theta_tower`, `mock_theta_digits`)

**Why** (NEXT.md item 11). f(q) = Σ q^(n²)/(−q;q)_n² is not modular: q^(−1/24)·f(q) is the holomorphic part of a weight-1/2
harmonic Maass form. Its shadow is a weight-3/2 unary theta function g = Σ χ₁₂(n)·n·q^(n²/24). Do its towers have laws, and
what sets them?
- Computed by the Appell–Lerch identity f(q) = (1 + 4Σ_(n≥1) (−1)ⁿq^(n(3n+1)/2)/(1 + qⁿ))/(q;q)_∞. It agrees exactly with the
  definition for n < 400: 1, 1, −2, 3, −3, 3, −5, 7, −6, 6, −10, ….
- Tower α(ℓᵇn + δ_b), 24δ_b ≡ 1 (mod ℓᵇ), as for p(n).

**The data alone** (`mock_theta_tower_N16000000_5_…_43.txt`; counted). A law on lane ℓ at every prime where the levels
allow the degree:
- 5: x² + 1 (73 equations);
- 7: x² − 1 (56);
- 11: (x² + 1)² (19);
- 13: x⁴ − 1 (13);
- from 17 on: no law of the degrees the levels allow.
In two-level steps the constants are −1, +1, −1 (twice), ±1. That is (−3|ℓ) at all four primes, plus one more −1 at 13.

**A reading** (derived). On weight 3/2, T(ℓ²) multiplies a unary theta function Σψ(m)m·q^(m²) by ψ(ℓ)ℓ + (Nebentypus
character)·(−1|ℓ). With the character that makes g an eigenform, g's eigenvalue is χ₁₂(ℓ)·(ℓ + (−1|ℓ)) or χ₁₂(ℓ)(ℓ + 1). The
first is ≡ χ₁₂(ℓ)(−1|ℓ) = (−3|ℓ) (mod ℓ), which is what the data show; the second would be (3|ℓ), which fails at 7. So the mock
form's tower seems to be steered by its shadow, a "p-adic coupling", as Guerzhoy–Kent–Ono call it for other weights.

**PREDICTED before the deeper digits are read** (`mock_theta_digits`). At 5 and 7, where the law on lane ℓ has a single
two-level constant, that constant is exactly the shadow's eigenvalue: **−6 at 5 and −6 at 7**. The alternatives
(−3|ℓ)(1 + ℓ) (−6 and 8) and (−3|ℓ)(1 − ℓ) (4 and −6) are to be told apart by the digits.

**First digits read, at 5: the prediction above is WRONG** (`mock_theta_digits.txt`, which stopped after 5 on an empty level;
counted).
- At 5 the two-level ratio is 4, 4, 4, 4, … = **−1** in 5-adic digits, agreeing on 1, 2, 2, 3, 3, 4, 4 digits at levels
  1 → 3 … 7 → 9. It is not −6.
- **Corrected reading** (derived after seeing 5, before seeing 7). χ₁₂(ℓ)(ℓ + (−1|ℓ)) is the sum of the two roots of g's Hecke
  polynomial y² − χ₁₂(ℓ)(ℓ + (−1|ℓ))·y + (−1|ℓ)ℓ:
  - the unit root **(−3|ℓ)** exactly;
  - the root of slope 1, χ₁₂(ℓ)·ℓ.
  The tower's constant is the unit root, as for every family so far. The slope-1 root's piece fades one digit every two
  levels, which is the rate the digits at 5 come in.
- **PREDICTED for 7, not yet read:** the two-level constant is exactly +1 = (−3|7). The tower obeys the shadow's whole
  polynomial in two-level steps, y² + 6y − 7 = (y − 1)(y + 7), deeper than y − 1 alone. At 5, y² + 6y + 5 = (y + 1)(y + 5),
  is a check of this reading after the fact, not a prediction.

**Measured** (`mock_theta_digits.txt`; counted). Smallest depth over n ≥ 1, by starting level:

| law in two-level steps y | 5 | 7 |
|---|---|---|
| y − χ₁₂(ℓ)(ℓ + (−1\|ℓ)) (first prediction: −6) | 1, 1, 1, … | 1, 1, 1, … |
| y − (−3\|ℓ)(1 ± ℓ) (alternatives) | 1, 1, 1, … | 1, 1, 1, … |
| **y − (−3\|ℓ), the unit root** | 1, 2, 2, 3, 3, 4, 4 | 1, 1, 2, 2, 3, 3 |
| **the shadow's whole polynomial** | 2, 3, 3, 4, 4 | 2, 2, 3, 3 |

- **At 7, read after the corrected reading was written:** the ratio two levels apart settles on 1, 0, 0, … = **+1 = (−3|7)**,
  one digit more every two levels, like −1 at 5.
- The first prediction and the alternatives stop at depth 1 at both primes: refuted.
- The shadow's whole polynomial, (y + 1)(y + 5) at 5 and (y − 1)(y + 7) at 7, is one digit deeper than its unit root alone at
  both primes. The slope-1 root's piece is present, as the reading says.
- Beyond that the depth still grows one digit every two levels. So another piece of the same slope is in the tower, not yet
  identified. A holomorphic family through weight ℓ − 1 on level 6 would be the first place to look.
- **So:** the mock theta function's tower obeys a law, and its unit constant is exactly the unit root (−3|ℓ) of its shadow's
  Hecke polynomial. Its other root is present too.
  - That is the coupling of a mock form to its shadow, here at weight 1/2 with a unary theta shadow.
  - The literature (Guerzhoy–Kent–Ono's "p-adic coupling"; Garthwaite–Penniston on f(q)) is not checked here for this exact
    statement.
- 11 and 13 have one more constant on lane ℓ (a second −1 at 11, a −1 at 13), also not identified yet.

## Curves that never come back (`level1_tower` poleH at 5, weight −48, `never_back_check`)

**Why** (12 Nov). Under the step U_ℓ makes,
there are three kinds of curve that never return to their start, and each can be named in advance (derived):
1. **Curves without CM**: at every ordinary prime they approach their returning twin, one digit a step, forever (measured,
   "What the tower view reveals").
2. **CM curves where they are supersingular**: they are not in the ordinary dance at all, and the tower ignores them
   (measured: "P alone" at the inert primes, everywhere above).
3. **CM curves whose order's conductor contains ℓ** (Gross's quasi-canonical lifts): their Serre–Tate coordinate is an ℓᵏ-th
   root of unity, which the step raises to the ℓ-th power. So they **land exactly on the returning curve after k steps**
   and stay there. They never come back to their start, yet they stop drifting.

**The test, at 5.** Discriminant −275 = −11·5².
- H₋₂₇₅ was computed from the four CM values to 220 digits; the coefficients are integers to 10⁻¹⁸⁰.
- **Mod 5 it is (x − 2)⁴**: the four curves are the quasi-canonical lifts of level 1 above the returning curve j = −32768
  (CM by √−11), which reduces to 2.
- Three poles at weight −48, all reducing to j ≡ 2 mod 5:
  - (a) 1/(Δ⁴·H₋₂₇₅(j)): four curves that fall onto the twin after one step;
  - (b) 1/(Δ⁴·(j − 2)(j − 7)(j − 12)(j − 17)): four curves without CM, starting 1, 2, 1, 1 digits from the twin;
  - (c) 1/(E₄⁹·(E₄³ + 32768Δ)): a simple pole on the twin itself. Its other pole is at ρ, supersingular at 5, which the tower
    ignores.

**PREDICTED before the runs.** Law P(x)·(x − u⁻⁵⁰)·S(x), with u the unit root at 5 (trace −3) and P, S the weight −48 families.
- (c) is exact from level 1: depth grows at the next slope's rate.
- (a), never returning, **grows exactly like (c) from level 1**, since after one step its poles sit on the twin.
- (b) grows one digit a level (b + 1 + a constant), forever.

**Measured** (`never_back_check.txt`; towers `level1_tower_poleH{LIFT4,D275,NC4}_N16000000_5`; families
`forms1_family_w-48_5.txt` and `slope1_families_5_w-48_N20.txt`: slopes 0, 1, then 4, 5, 6; counted). Law
P(x)·(x − u⁻⁵⁰)·S(x), smallest depth over n ≥ 1, by starting level (lane 5¹¹, 64 checks each):

| pole at | 1 | 2 | 3 | 4 | 5 | 6 | 7 |
|---|---|---|---|---|---|---|---|
| (c) the returning curve j = −32768 | 10 | 11 | 11 | 11 | 11 | 11 | 11 |
| (a) the four curves of discriminant −275 | 6 | 10 | 11 | 11 | 11 | 11 | 11 |
| (b) four curves without CM | 4 | 5 | 6 | 7 | 8 | 9 | 10 |

- **(a) never returns to its start, yet from level 2 it is indistinguishable from the returning curve.** It runs one level
  behind it, which is the one step its poles take to land. It gains 4 digits in that one level, the next slope's rate.
  So it does not drift.
- **(b) gains one digit a level, forever**, the drift of curves without CM.
- Without S all three grow one digit a level (3…10, 2…9, 1…8): the slope-1 family hides the difference until it is in
  the law.
- **So the tower predicts and shows all three kinds of curve that never come back:**
  - those without CM, drifting forever;
  - the supersingular ones, absent;
  - the CM curves of conductor divisible by ℓ, which land on a returning curve after k steps and then repeat.
  Together with the returning curves (period 1, 2, 3 above), this is the whole of Frobenius's motion on the ordinary
  curves at ℓ, read from coefficients.

**The unidentified piece: a candidate from the forms** (12 Nov; counted, then predicted).
- In the class of level-6 newforms with U₂, U₃ signs ((2|ℓ), (3|ℓ)) at weight ℓ − 1 (p(n)'s class rule, one weight up),
  there is exactly one newform at every ℓ from 5 to 23. Each has **a_ℓ ≡ 1 (mod ℓ)**, the mark of an Eisenstein congruence.
  Counted with `class_basis_exact`, its common factor removed.
- If its family enters the tower with the sign (2|ℓ), it explains every lane-ℓ law. (2|ℓ) is the only simple character
  that fits all four:
  - at 5 and 7 its constant (2|ℓ)·1 equals the shadow's (−3|ℓ) on the lane, so it was hidden;
  - at 11 it doubles the shadow's −1, which is why the data needed (y + 1)²;
  - at 13 it is the extra −1.
- **PREDICTED before the family's digits are computed:** with u_f(0) = that family's U_ℓ-value at weight 0 (Hida, from the
  weights (ℓ − 1)(1 + x)), the law (y − (−3|ℓ))·(y − χ₁₂(ℓ)ℓ)·(y − (2|ℓ)·u_f(0)) holds deeper than the shadow's polynomial
  alone, wherever the levels reach.
- If it does not go deeper at 5 and 7, the family is absent there, or hidden below the depths reached. The piece of
  slope 1 seen at 5 and 7 then stays open.

**Measured, the candidate family** (`mock_theta_family_5_7_11_13.txt`; counted).
- The class family is one ordinary form at every weight (ℓ − 1)(1 + x), x = 0..8, at all four primes. **Its U_ℓ-value at
  weight 0 is 1 to every digit computed** (9 digits; Newton differences of valuations 0, 1, 3, 3, 4, 5, 7, 7, 8 at 5): a cusp
  family meeting the Eisenstein value exactly at weight 0.
  - A possible reason (a reading, not checked): an Euler factor at 2 or 3 vanishes there (a "trivial zero"), where cusp and
    Eisenstein families can meet.
- So the tower's unit constants are two pure signs: the shadow's (−3|ℓ) and the family's (2|ℓ). They are equal at 5, 7 and
  11 and differ at 13. That explains every lane law: hidden at 5 and 7, the double root (y + 1)² at 11, the extra −1 at 13.
- **The "deeper" half of the prediction failed at 5 and 7**: S(y)·(y − (2|ℓ)) goes no deeper than S(y) (2, 3, 3 at 5; 2, 2
  at 7). At 11 and 13 the full law needs more levels than n < 1.6·10⁷ gives.
- **A correction to "Measured" above.** "The shadow's whole polynomial is one digit deeper, so its slope-1 root is present"
  was over-read: any factor y − ℓ·(unit) gains one digit on any slope-1 piece.
  - A two-term fit of the remainders (after y − (−3|ℓ)) shows the slope-1 root that is really there. Its unit part is 2
    (mod 25) at 5 and 4 (mod 7) at 7; the shadow's χ₁₂(ℓ)ℓ would give −1. The fit also shows a root of slope 3 at 5 and 2
    at 7.
  - So the slope-1 piece is a family of its own, not on level 6 at weight ℓ − 1 (there is no such form at weight 4). It is
    not identified.
  - Whether the shadow's slope-1 root is in the tower at all is open.
- **The unidentified slope-1 constant, to the digits the data hold** (two-term fit of the remainders after y − (−3|ℓ), then
  Hensel on y = ℓ·z; counted):
  - at 5, ℓ·z with z ≡ 502 (mod 5⁵) (digits 2, 0, 0, 4, 0);
  - at 7, z ≡ −3 (mod 7²) (digits 4, 6).
  - A simple character value is ruled out: at 5, −3 would give 622 (mod 5⁵).
  - These digits are the target for any candidate family. Levels 12, 18 or 36, through the Shimura level 72 of Γ0(144), are
    the first places to look.

**Measured to 10⁸** (`mock_theta_tower_N100000000_5_7_11_13_17_19_23.txt/.json`).
- How it was computed: the Appell–Lerch numerator was found exactly with the venv's numpy (`mock_theta_appell.py`; 8164
  terms, largest coefficient 132) and passed to WSL's FLINT as a file. Counted.
- **The laws at 5, 7, 11, 13 hold with more equations**: x² + 1, x² − 1, (x² + 1)², x⁴ − 1, from 87, 67, 26, 23 equations.
- **The two signs, exactly, on lane ℓ:** (y − (−3|ℓ))(y − (2|ℓ)) holds at every level and n at 5, 7, 11, 13 (65, 45, 26, 23
  checks).
- **Depth with the exact signs:**
  - at 5, y + 1 gives 1, 2, 2, 3, 3, 4, 4, 5, 7 (levels 1–9);
  - at 7, y − 1 gives 1, 1, 2, 2, 3, 3, 4.
  - Both deepen one digit every two levels, the slope-1 piece's rate, as exact ±1 constants would. At 11 and 13 the levels
    allow only depth 1.
- **At 17, 19 and 23 the two-sign law fails on lane ℓ, even from level 2**, where every slope-1 piece is already ≡ 0. So the
  tower has more unit constants there.
  - They first appear where the other level-6 classes first have forms at weight ℓ − 1. At 17 those are a ≡ 5 and a ≡ 4
    (mod 17), beside the family's 1.
  - The data-only search finds no law of degree ≤ 4. A law of degree 3 or 4 in y needs 7 to 9 levels: n near 17⁸.

## The drift is the next Frobenius root: a test of the mechanism behind Zhang's Conjecture 1.5 (`drift_roots_check`)

**Why** (12 Nov).
- Zhang, "Elliptic curves and Fourier coefficients of meromorphic modular forms" (arXiv 2510.23200), conjectures (Conj.
  1.5) Atkin–Swinnerton-Dyer-type congruences for E₄/(j − j(C))^r. They are governed by the characteristic polynomial of
  Frobenius on Sym²C, modulo about p^(3l).
- Read through a fetched summary of the paper; the exact statement is to be checked against the text.

**The mechanism, derived.**
- At level b the pole sits at φᵇ(C), ε_b ≈ ℓ^(b+s) from the twin x̃ (Serre–Tate).
- Expand Φ(x_b) = Φ(x̃) + ε_b·Φ′(x̃) + ε_b²/2·Φ″(x̃) + …. The i-th term is a pole of order i + 1 at the CM point x̃.
- By the residue formula, U_ℓ multiplies it by π̄ⁱπ^(w−2−i) = **ℓⁱ·u^(w−2−2i)**.
- So a pole at a curve without CM carries every constant ℓⁱ·u^(w−2−2i), i = 0, 1, 2, …, not only u^(w−2). The drift
  measured before (depth b + s with one factor) is the i = 1 term.
- At weight 4 the first three are u², ℓ and ℓ²u⁻² = ū²: the eigenvalues of Frobenius on Sym²C, as Zhang's conjecture has
  them. The next, ℓ³u⁻⁴, has slope 3, which is why the depth there grows about 3 digits a level.

**PREDICTED before the run** (j₀ ≡ 2 mod 5, weight −12, the towers of "What the tower view reveals").
- With the i = 1 factor, P(x)·(x − u⁻¹⁴)·(x − 5u⁻¹⁶), the depth from level b is min(2b + 2s, 2b + 3, 11).
  - That is 4, 6, 8, 10, 11 for j = 2 (s = 1), against b + 1 with one factor.
  - It is the twin's 2b + 3 for every s ≥ 2.
- With i = 2 as well (x − 25u⁻¹⁸), j = 2 also reaches the twin's 2b + 3.
- The pole ON the twin is unchanged by the extra factors (it has no i ≥ 1 terms).
- If this holds, our towers show the Sym² structure of Zhang's conjecture as the Taylor expansion of a moving pole. That
  is a mechanism a proof could follow.

**Measured** (`drift_roots_check.txt`; counted). Depth by starting level, lane 5¹¹:

| pole | one root (i = 0) | two roots (i ≤ 1) | three roots (i ≤ 2) |
|---|---|---|---|
| the twin, −32768 | 5, 7, 9, 11, … | 6, 8, 10, 11, … | 9, 11, 11, … |
| j = 2 (s = 1) | 2, 3, 4, 5, 6, 7, 8, 9 | **5, 7, 9, 11, …** | 8, 10, 11, … |
| s = 2, 3, 4, 6 | b + s (until 2b + 3) | 6, 8, 10, 11, … | 9, 11, 11, … |

- **The drift is the next root.** With the i = 1 root 5·u⁻¹⁶ in the law, j = 2 gains two digits a level, not one. Every
  curve without CM then follows the twin's pattern: one digit behind it at s = 1, level with it from s = 2. The main claim
  holds.
- **Two parts of the prediction were off.**
  - j = 2 with two roots reached 2b + 3, one digit more than the 2b + 2 predicted.
  - The twin was predicted to be unchanged by the extra factors; it is not. Any factor x − c with v(c) = 1 gains one digit
    against the weight −12 families of slope 2 (their λ − c has valuation 1), and I did not account for that.
  - With the third root the twin jumps to 9 at level 1, three digits more. 25·u⁻¹⁸ is close to a slope-2 family's value
    there: not understood.
- **What this says about Zhang's Conjecture 1.5.**
  - A pole at a curve without CM carries the constants ℓⁱ·u^(w−2−2i), i = 0, 1, 2, …, from the Taylor expansion of a
    moving pole. At weight 4 the first three are Frobenius on Sym²C.
  - The data at weight −12 behave that way: each added root removes one order of drift.
  - This is a mechanism and evidence, not a proof. The conjecture's exact wording, read so far only through a summary, is
    still to be checked. The next test is Zhang's own form, E₄/(j − 2), at weight 4.

## Zhang's conjectures for E₄/(j − j(C))^r, tested and extended (`zhang_check`)

**The statements** (P. Zhang, "Elliptic curves and Fourier coefficients of meromorphic modular forms", arXiv 2510.23200v3).
Quoted from the paper's text, extracted from its PDF here. Hypotheses throughout: C/ℚ with j(C) ∉ {0, 1728}, p ≥ 5 good, and
v_p(j(C)) = 0 = v_p(j(C) − 1728).
- **Theorem 1.3 (proved):** a_p(E₄/(j − j(C))) ≡ a_p(C)² (mod p).
- **Conjecture 1.4:** if p is supersingular for C, then a_p(E₄/(j − j(C))) ≡ 0 (mod p²).
- **Conjecture 1.5:** for 1 ≤ r ≤ 3, with F^(r) = E₄/(j − j(C))^r and all n, l ≥ 1,
  a_(npˡ) + c₁·a_(npˡ⁻¹) + c₂·a_(npˡ⁻²) + c₃·a_(npˡ⁻³) ≡ 0 (mod p^(3l−2−r)),
  where X³ + c₁X² + c₂X + c₃ = (X − α²)(X − p)(X − β²), α, β the roots of X² − a_p(C)X + p (Frobenius on Sym²C). Here
  c₁ = −(a² − p), c₂ = p(a² − p), c₃ = −p³.

**Our reading** (derived; README "The drift is the next Frobenius root").
- A simple pole at a curve without CM carries the constants pⁱ·u^(w−2−2i). At weight 4 these are u², p, ū² (Zhang's Sym²),
  then p³u⁻⁴, p⁴u⁻⁶, ….
- **The roots from i = 3 on are p³ times the roots of a weight −2 pole** (u⁻⁴, p·u⁻⁶, …). By Bol's identity a pole of
  order ≥ 4 in weight 4 is θ³ of a weight −2 object, and θ³ multiplies each coefficient by n³, so its U_p-values by p³.
- So the cubic leaves out first a root of slope 3, which is why Zhang's modulus grows by 3 a level.

**PREDICTED before the run** (j(C) = 2, our curve; F^(r) for r = 1, 2, 3; towers to n < 1.6·10⁷; lane p^K < 2⁶²).
- **Zhang's law holds at every prime, level and n tested,** to depth at least 3l − 2 − r. Its depth grows about 3 a level
  at the ordinary primes.
- **Extension (new):** at the ordinary primes, adding the root p³u⁻⁴ makes the law grow about 4 digits a level; adding
  p⁴u⁻⁶ as well, about 5. This holds until the weight-4 cusp families take over, at a slope these runs will show.
- **Conjecture 1.4 at the supersingular primes of C (29 and 31):** v_p(a_p) ≥ 2. The reading above does not decide this;
  it is Zhang's.
- **Noted before the run (derived, a caution on the extension).** Coleman: a weight-4 overconvergent form of slope below 3
  is classical. Below slope 3 the only candidates are E₄'s ordinary root 1 (slope 0) and the p-new forms of level p (U_p = ±p,
  slope 1). Zhang's 3-a-level rate therefore says F^(r) has no part along either. At slope 3 more than our p³u⁻⁴ may sit:
  E₄'s critical root p³, and θ³ of weight −2 forms of slope 0. If one of them is present, L₄ = L₃·(X − p³u⁻⁴) gains
  nothing over L₃. The run tests L₃·(X − p³) as a probe for this.
- **A trial run (N = 2·10⁵, not kept; levels 1–7 at 5, 1–3 at 29) changed the probes.**
  - At 5, L₃ met Zhang's bound exactly, and L₄ fell back to 3 digits a level. L₄·(X − p³), which kills both slope-3 roots,
    grew 4 a level.
  - At the supersingular 29, the shorter law X² − p² reached 3l − r, 2 deeper than Zhang's cubic. That would make
    Conjecture 1.4 its first level (l = 1, r = 1: a(29n) ≡ 0 mod 29²).
- **PREDICTED for the full run** (N = 1.6·10⁷; j₀ = 2; p = 5…43 and the supersingular 307):
  - Zhang's law holds everywhere.
  - At every ordinary prime: L₄·(X − p³) grows 4 a level, and L₅·(X − p³) grows 5.
  - At 29, 31 and 307, at every level there: a(npˡ) ≡ p²·a(npˡ⁻²) (mod p^(3l−r)).

### Measured (`zhang_check_j02_N8000000_5_…_307.txt` and `.json`; n < 8·10⁶; COUNTED)

The first full run at n < 1.6·10⁷ crashed in the Windows FLINT build: the 63-bit lanes double the size of each product.
It left a header-only file, which was removed.

| test | result |
|---|---|
| Theorem 1.3 | holds at all 13 primes |
| Conjecture 1.4 | v_p(a(p)) = 2 at 29, 31 and 307; for every n ≤ 200 prime to p, v_p(a(np)) = 2 or more |
| Conjecture 1.5 | holds at every prime, level and n, for r = 1, 2, 3. Levels tested: 9 at 5, 8 at 7, 6 at 11 and 13, 5 at 17–23, 4 at 29–43, 2 at 307. The least depth over n is exactly 3l − 2 − r at every level from 2: **the bound is sharp** |
| supersingular, X² − p² | exactly 3l − r at 29, 31 and 307, every level, r = 1, 2, 3: **held as predicted** |
| ordinary, L₄·(X − p³) | 4 a level at 5, 7, 11, 13 and 17: **held**. No gain at 19 and 23: **failed**. Not testable at 37–43 (degree 5, only 4 levels) |
| ordinary, L₅·(X − p³) | 5 a level at 5 (levels 6–8) and one step at 7, then the lane's full width; too few levels elsewhere |
| our root alone, L₄ | no lasting gain anywhere: **wrong as first predicted**. The caution written before the run was right |

Depth of L₄·(X − p³) (counted): 4l − 6 − r at 5, 11, 13 and 17; 4l − 5 − r at 7.

**What it means:**
- Zhang's Conjecture 1.5 holds and is sharp in every case tested: the depth grows exactly 3 a level.
- **At supersingular primes Zhang's cubic is more than needed.** Its repeated root −p needs one factor, not two. The
  minimal polynomial X² − p² gives a(npˡ) ≡ p²·a(npˡ⁻²) (mod p^(3l−r)), 2 digits deeper. At l = 1, r = 1 this is
  Conjecture 1.4, and it holds for every n prime to p, not only n = 1. There, the two conjectures are one law.
- Below slope 3 the tower holds exactly Zhang's three roots. At slope 3 it holds at least our p³u⁻⁴ and E₄'s p³.

**Why 19 and 23 (derived after the run; one necessary condition counted):**
- Weight p − 3 is in weight −2's class mod p − 1. At 5, 7, 11, 13 and 17 there are no level-1 cusp forms of weight
  p − 3 (weights 2, 4, 8, 10, 14).
- At 19 and 23 there is one each, ΔE₄ (weight 16) and ΔE₈ (weight 20), and it is ordinary: a₁₉ ≡ 13 (mod 19) and
  a₂₃ ≡ 7 (mod 23).
- Hida's family through it has a weight −2 member with a unit U_p-value λ ≡ a_p (mod p). θ³ of that member is a
  weight-4 form of slope 3 with U_p-value p³λ: a third slope-3 root that the laws above do not kill.
- No family passes through weight 4 itself (S₄(1) = 0; at 19, weight 22 is in weight 4's class and a₁₉(ΔE₁₀) ≡ 0). That
  is why nothing of slope 0 spoils Zhang's law.
- This is the shape of the partition constants again (a family through weight ℓ − 3, read at weight −2), here at
  level 1.

**PREDICTED next:**
- At 19 and 23, the law L₄·(X − p³)·(X − p³λ) grows 4 a level, with λ from Lauder's algorithm at weight −2. This needs
  level 6 at 19 (WSL, n < 10⁸).
- At 37, 41 and 43, extra slope-3 roots appear wherever weight p − 3 (34, 38, 40) has ordinary forms.
- **PREDICTED before the Lauder runs** (`slope1_families.py 19 −2 14`, `23 −2 13`): at weight −2, the slope-0 factor is
  (z − 1)(z − λ). The 1 is the Eisenstein family; θ³ of it gives E₄'s critical root p³. λ ≡ 13 (mod 19) and
  λ ≡ 7 (mod 23), the a_p of ΔE₄ and ΔE₈.
- **Measured** (`slope1_families_19_w-2_N14.txt`, `slope1_families_23_w-2_N13.txt`): at weight −2 the slope-0 factor is
  (z − 1)(z − λ) at both primes. **Held as predicted:** λ = …, 9, 2, 17, 17, 10, 6, 17, 18, 6, 6, 6, 13 (base-19 digits,
  last digit first) ≡ 13 (mod 19), and λ ≡ 7 (mod 23).
  - Weight −2 also has 1 root of slope 1 at 19 and 2 at 23. θ³ of them gives slope-4 roots at weight 4.
  - So the full picture, derived: the weight-4 tower's roots are Zhang's three (u², p, ū²) plus p³ times the weight −2
    roots, which are the moving pole's u⁻⁴, p·u⁻⁶, … and the overconvergent spectrum 1, λ, ….
- **The six-root law cannot be tested yet.** Its first complete level at 19 is level 6. A 63-bit lane holds only 14
  digits of 19, less than Zhang's own depth there (15).
- **PREDICTED instead, from data already kept** (the n < 8·10⁶ towers):
  - At 19, the residual R(n) of L₄·(X − p³) at level 5 (its one complete level) is left with the single slope-3 piece
    θ³f. Here f is the weight −2 family member, whose coefficients are ≡ those of ΔE₄ (mod 19) for n prime to 19.
  - θ³ multiplies the n-th coefficient by n³. So R's leading digit runs like n³·a_n(ΔE₄): R(2)/R(1) ≡ 8·216 ≡ −1 and
    R(3)/R(1) ≡ 27·(−3348) ≡ 6 (mod 19), for r = 1, 2 and 3 alike.
  - A control: the E₄-critical piece would run like σ₃(n) instead (9 and 28 ≡ 9).
- **Measured** (`zhang_hida_check.txt`): **held, exactly.** At 19, level 5, R(2)/R(1) ≡ 18 ≡ −1 and R(3)/R(1) ≡ 6
  (mod 19), for r = 1, 2 and 3 alike. The control σ₃(n) would give 9 and 9.
  - The depths are 12, 11 and 10: Zhang's bound, so the leftover sits right at it.
  - So the piece that stops the extension at 19 is θ³ of the weight −2 member of the Hida family through ΔE₄. The chance
    that both ratios match by luck is about 1 in 361.
- **The picture of Zhang's towers, now (derived, and checked as above).** The weight-4 tower holds Zhang's three
  Frobenius roots u², p, ū², and then p³ times everything a weight −2 tower holds:
  - the moving pole's roots u⁻⁴, p·u⁻⁶, …;
  - the Eisenstein family's 1 (giving E₄'s critical p³);
  - the weight −2 members of Hida and Coleman families.
  This is why his modulus grows exactly 3 a level (Bol's identity: θ³ carries weight −2 to weight 4).
- **Next:** the same n-test at 23 needs n = 2 at level 5 (n < 1.3·10⁷: WSL, since Windows crashed there with 63-bit
  lanes); a second curve (j = 4); the six-root law needs a lane wider than 64 bits.

## The towers of 73 and 97, by the forms method (13 Nov 2026; `partitions_flint_family.py`, then checks)

The user asked for the laws and constants of 73 and 97, the two primes below 100 that are 1 more than a multiple of 24.
There 1/24's digits repeat every digit, and a one-level zoom fits the indexing (`side_prime_horns`, point 16). The
method is the one of Statements A and B: the class S_ℓ (level 6, weight ℓ − 3, signs w₂ = −(2|ℓ), w₃ = −(3|ℓ)), its
Hida families read at weight −2, P_ℓ(x) = ∏(x − (3|ℓ)·u_f(−2)), and the law Σ [xⁱ]P_ℓ·H_(b+2i) ≡ 0.

**PREDICTED before the runs:**
- dim S₇₃ = R₇₃ = ⌊85/24⌋ = 3 and dim S₉₇ = R₉₇ = ⌊109/24⌋ = 4. Every form is ordinary (Statement A2: counted to 1,297
  on lane ℓ).
- The constants' first digits are (3|ℓ)·a_ℓ(f) mod ℓ (A3), and (3|73) = (3|97) = +1.
- **On exact partition numbers (FLINT):** 73's law, built from the forms alone, holds on lane 73 at levels (1, 3, 5, 7),
  for n = 0 and n = 1.
- **The one-level zoom, a guess (not derived):** at 73 and 97 a law in single-level steps of degree R holds,
  Σ κᵢ·H_(b+i) ≡ 0 (mod ℓ) for i = 0 … R. At the controls 71 and 79 (R = 3) none of degree 3 holds. A degree-2R law,
  P(x²), always would.

**From the forms (`partitions_flint_family_73.txt`, `_97.txt`; 105 s and 630 s; COUNTED):**
- **73:** dim S₇₃ = 3, found exactly over the integers. **Held.**
  - P₇₃ on lane 73 is x³ + 65x² + 71x + 29: (3|Q)·a_Q at Q = 73, the U₇₃ polynomial at every weight 70 … 790.
  - On lane 73¹⁰ it is x³ + 1274766274151810209·x² + 1391144417610593507·x + 640148276439656497.
  - It does not factor over the 73-adic integers: the three constants are conjugate.
- **97:** dim S₉₇ = 4, found exactly over the integers. **Held.**
  - P₉₇ on lane 97 is x⁴ + 22x³ + 63x² + 3x + 46.
  - One constant is a 97-adic integer, digits lowest first 28, 26, 81, 49, 38, 24, 67, 78, 84, 46. The other three are
    the roots of the cubic [67174571304055357052, 48281256314335232937, 18624890627406768819, 1] on lane 97¹⁰.
- Every form is ordinary at both primes: the class's piece keeps its dimension at every weight.

**PREDICTED before reading the check (derived from these constants, not a guess): no law one level apart at 73 or 97.**
- A law in single-level steps of degree R needs a one-level operator B with B² = A, the two-level step whose
  eigenvalues are the constants α_f.
- If some α_f is not a square in its residue field, B would need both √α_f and −√α_f as eigenvalues. Then α_f would be
  a double eigenvalue of A, but A's eigenvalues are distinct (Statement D2).
- 97's integer constant is ≡ 28 = 4·7, and (7|97) = −1, so it is not a square.
- 73's constants lie in the cubic field of 73³ elements, where an element is a square exactly when its norm is a square
  mod 73. The norm is −29 ≡ 44 = 4·11, and (11|73) = −1, so they are not squares.
- So the least law in single steps should be P(x²), the two-level law, of degree 2R:
  - at 73: no fit for d = 1 … 5, and exactly one for d = 6;
  - at 97: none up to d = 4;
  - at the controls 71 and 79: none up to d = 3.

**On exact partition numbers (`partitions_73_97_check.py`, `.txt`; FLINT, 100 s; COUNTED). Every prediction held.**
- **73 along ℓ:** the 4-term law Σ [xⁱ]P₇₃·H_(1+2i)(n) ≡ 0 (mod 73), built from the forms alone, holds for 10 of 10
  values of n. That uses levels 1, 3, 5, 7, so partition numbers up to about 11 million digits. No law with fewer
  terms two levels apart fits (2 and 3 terms: dimension 0).
- **73 along Q:** the 3-step law from the forms holds on 12 of 12 chains at each of Q = 5, 7, 11, 13; no law of
  degree 2 fits any Q.
- **97 along Q:** the 4-step law from the forms holds on 12 of 12 chains at each of Q = 5, 7, 11, 13; no law of degree 3
  fits. (Along ℓ, 97's law needs level 9, partition numbers of about a billion digits, and was not run.)
- **One level at a time:**
  - At 73, the laws in single-level steps fitting n = 0 … 9 have dimension 0 for d = 1 … 5 and exactly 1 for d = 6.
    The one law is P₇₃(x²), the two-level law, since that law fits and the space is one-dimensional.
  - At 97, nothing up to d = 4. At the controls 71 and 79, nothing up to d = 3.

**So the towers of 73 and 97 are like every other prime's: laws two levels apart.** The one-level zoom that the
indexing allows at ℓ ≡ 1 (mod 24) is not used. The reason comes from the forms: a one-level law would need square roots
of the constants, and at both primes the constants are not squares.

**The laws, as numbers** (lowest coefficient first; base-ℓ digits, lowest first):
- **73:** H_(b+6) + c₂·H_(b+4) + c₁·H_(b+2) + c₀·H_b ≡ 0, with
  - c₀ = 29, 36, 72, 47, 23, 53, 44, 56, 63, 10 …
  - c₁ = 71, 44, 28, 21, 12, 1, 5, 0, 46, 23 …
  - c₂ = 65, 49, 54, 40, 63, 71, 45, 50, 47, 21 …
  - (73-adic, known to 10 digits; the three constants are its roots, conjugate in the unramified cubic extension.)
- **97:** H_(b+8) + c₃·H_(b+6) + c₂·H_(b+4) + c₁·H_(b+2) + c₀·H_b ≡ 0, with
  - c₀ = 46, 30, 23, 24, 51, 64, 90, 72, 74, 2 …
  - c₁ = 3, 9, 64, 22, 40, 73, 93, 94, 26, 91 …
  - c₂ = 63, 70, 67, 40, 26, 95, 93, 96, 26, 39 …
  - c₃ = 22, 17, 55, 22, 16, 46, 26, 57, 60, 74 …
  - It factors as (x − α)·(cubic), with α = 28, 26, 81, 49, 38, 24, 67, 78, 84, 46 … a 97-adic integer. Checked: the
    product equals the quartic to all 10 digits.

## Lucky primes for 73: Ramanujan-type laws modulo 73 (13 Nov 2026; `partitions_flint_lucky_dual`)

**The question.** Statement B gives c(Q²M) = (λ + εδ)·c(M) on lane ℓ. When the forms' number at Q is the same for every
form of the class, λ + εδ can be 0 and then p(N) ≡ 0 (mod ℓ) on whole progressions. At 73 the class is one conjugate
triple, so such a Q is rare (about 3 in 73³ ≈ 1.2 million primes). The scan of the nineteenth time stopped at 10⁶.

**A lighter test (derived before the run).** For every form f of the class and every n, a_n(f) = a_1(T_n·f). So if an
operator t of the Hecke algebra has a_1(t·f) ≡ 0 for every f of the class, then every coefficient of t·f is ≡ 0, and t
is 0 on the class modulo ℓ. With t = (3|Q)·T_Q − c: the class law needs only (3|Q)·a_Q(S_j) ≡ c·a_1(S_j) for the
class's basis S_j, so q-expansions to Q_max terms, not to Q_max·(largest pivot), which was 5 to 7 times longer.

**PREDICTED before the run:**
1. The lighter test gives exactly the lists of the nineteenth time: the same Q for 37, 41, 43, 47, 53, 59 and 61 up to
   10⁶, and the same counts for 13 to 31 (derived; any difference means the test is wrong).
2. 73, every prime Q up to 5·10⁷ (3,001,134 primes): one piece over the field of 73³ = 389,017 elements, so each kind
   of law (class ε = +1, class ε = −1, Ono's) has chance (q + χ(c² − 4D))/(q² − 1) per prime, about 7.7 laws of each
   kind and 23 in all (derived; the exact sums are printed with the counts). 61 to 5·10⁷ as a second prime with three
   forms (pieces 1 and a pair; about 13 of each kind).
3. Each class law at Q is a congruence p(73·Q³·n + β) ≡ 0 (mod 73) for (Q − 1)/2 values of β. Each Ono-type law is
   p(73·Q⁴·n + β) ≡ 0 (mod 73) for Q − 1 values of β (derived from 24N − 1 = Q²·M).
4. On exact partition numbers (FLINT): for the smallest class laws, the two smallest members are ≡ 0 (mod 73), and a
   member of the other class at the same Q is not (chance 72/73). Ono-type members at 73 need p(N) with N ≈ 3·Q³, out of
   reach unless Q is small.
5. Literature: the published lists of such congruences stop at modulus 31 (Weaver, Johansson) and 37 (Ahlgren–Allen–
   Tang); no list for 73 is expected.

**Counted.**
1. **The lighter test is right** (`partitions_flint_lucky_dual_37_43_41_47_53_59_61_13_17_19_23_29_31_max1000000.txt`):
   39 of 39 comparisons with the whole-matrix test of the nineteenth time are the same (every list and every count, 13
   primes, every Q up to 10⁶). It took 16 s instead of about 6 minutes.
2. **73 and 61, every prime Q up to 5·10⁷** (3,001,131 primes each; `partitions_flint_lucky_dual_73_61_max50000000.txt`,
   4 minutes):

   | ℓ | the class on lane ℓ | class ε = +1 | class ε = −1 | Ono's kind |
   |---|---|---|---|---|
   | 73 | one piece of degree 3 | 10 (7.7 expected) | 7 (7.7) | 8 (7.7) |
   | 61 | 1 and a pair | 12 (13.2) | 14 (13.2) | 9 (13.2) |

   - The six counts give χ² = 2.3 on 6 (derived): within chance. This is the first test of the "3 in ℓ³" law with
     laws actually found; to 10⁶, 61 had 0.3 expected of each kind and 0 counted.
   - The smallest laws at 73: class ε = −1 at Q = 11,367,817, class ε = +1 at 11,421,973, Ono's kind at 1,109,167.
3. **As congruences** (`partitions_lucky_progressions.txt`; derived by arithmetic from the lists), each law gives
   (Q − 1)/2 progressions (Q − 1 for Ono's kind). The smallest of the first laws:
   - p(107239327321929851108449·n + 9040523964878168) ≡ 0 (mod 73) for every n ≥ 0; the modulus is 73·11367817³.
   - p(108779296991224840285141·n + 18650553916953050) ≡ 0 (mod 73); the modulus is 73·11421973³.
   - p(110486857112773659249495433·n + 70558828191079138941) ≡ 0 (mod 73); the modulus is 73·1109167⁴ (Ono's kind).
   - p(27241927820165199472207·n + 1633490502450505) ≡ 0 (mod 61); the modulus is 61·7643683³.
4. **Literature (QUOTED; a web search and three papers read):**
   - Ryan, Scherr, Sirolli and Treneer (arXiv 1911.05799) list congruences of this Q³ kind only up to modulus 17
     (Q < 10⁵).
   - Ahlgren, Beckwith and Raum (arXiv 2006.07645) prove that p(ℓQ^m·n + β) ≡ 0 (mod ℓ) has infinitely many examples for
     m ≥ 3, for every ℓ ≥ 5. They know none for m = 1 or 2, and Ahlgren and Beckwith (arXiv 2508.19512) prove that the
     primes Q with an m = 1 congruence have density 0.
   - Our laws are m = 3 (class laws) and m = 4 (Ono's kind), where they occur at a positive density of primes Q. So
     they do not clash with that scarcity.
   - No explicit congruence modulo 61 or 73 was found in print: new as far as we read.

**PREDICTED before the second check (Ono's kind, derived from Statement B).** An Ono-type law at Q means λ ≡ 0 for every
form of the class. Then the step gives c(Q²M) ≡ ε·δ·c(M) = −(M|Q)·(−3|Q)·Q⁻²·c(M) for every M with ℓ | M, Q ∤ M, in both
classes (as checked for 2999 at 37 on 12 Nov). For the two smallest Ono-type Q at 73 (1,109,167 and 7,725,269) and at 61
(2,144,719 and 2,146,489): two M of each class, all 8 values as stated. One such value holds by chance with probability 1/ℓ.

**Counted on exact partition numbers (FLINT; every prediction held):**
- **The control rule on laws known before** (`partitions_flint_lucky_check_73.txt`): 14 of 14 controls at 37, 43 and 53
  take exactly 2ε·(−3|Q)·Q⁻²·p((M′ + 1)/24), the value Statement B gives. They are not merely "not 0".
- **73's class laws** (same file; 73 minutes in all, numbers of about 10⁸ digits):
  - Q = 11,421,973 (ε = +1): p(18,650,553,916,953,050) and p(28,174,241,023,482,267) are multiples of 73. The control
    p(9,126,866,810,423,833) ≡ 61, as predicted.
  - Q = 11,367,817 (ε = −1): p(9,040,523,964,878,168) and p(27,907,704,413,319,562) are multiples of 73. The control
    p(18,474,114,189,098,865) ≡ 54, as predicted.
  - The next laws (Q = 22,076,447 and 24,700,997) need N above 3·10¹⁶ and were not computed.
- **61's class laws:** at Q = 4,425,853, 7,643,683, 11,452,607 and 11,534,903, all 8 members are multiples of 61, and
  all 4 controls take exactly the predicted value (29, 30, 17, 23). The largest was p(28,068,847,387,161,157).
- **Ono's kind, through its step** (`partitions_flint_lucky_check_ono.txt`, 33 minutes): 16 of 16 values as predicted
  (73 at Q = 1,109,167 and 7,725,269; 61 at Q = 2,144,719 and 2,146,489; 2 M of each class each; up to
  p(17,244,970,086,815,564)).
- **Chance:** a wrong law would pass each value with probability 1/ℓ. Here 12 members, 6 class controls and 16 Ono-type
  values all passed.

## 97's tower is two towers: the simple part split off (13 Nov 2026; `partitions_97_split`)

**The idea.** P₉₇ = (x − α)·(cubic): the class of 97 has one form g whose constant α = 28, 26, 81, … is a 97-adic
integer, and three forms conjugate in the cubic field. Each step along Q acts on the tower form by form, with the forms'
own numbers. So a combination of the steps at one prime Q′ can cancel the three forms and keep only g:
G(M) = Σ_j κ_j·c(Q′^(2j)·M), j = 0 … 3, with κ from h(T), where h is the cubic factor of the polynomial of
(3|Q′)·a_Q′ on the class and T is the step operator (c(M) = p((M + 1)/24), ℓ | M, Q′ ∤ M; κ depends on (M|Q′)). Q′ = 5
here; the numbers at 5 come from the forms alone, and the four-step law at 5 held on 12 of 12 chains at 97.

**PREDICTED before the run (derived):**
1. From the forms: the class has pieces of degrees 1 and 3 on lane 97, and the piece of degree 1 has (3|97)·a₉₇(g) ≡ 28,
   the first digit of α.
2. **Along 97:** G(97²M) ≡ 28·G(M) (mod 97) for M = 97·(24n + 23), n = 0 … 19 (5 ∤ M). This uses levels 1 and 3 only,
   partition numbers below 10¹². By contrast, no single number links c(97²M) to c(M) for the whole tower.
3. **Along other primes:** G(Q²M) ≡ (λ_g(Q) + ε·δ)·G(M), with λ_g(Q) = (3|Q)·a_Q(g), ε = (M|Q), δ = −(−3|Q)·Q⁻², for
   Q = 7, 11, 13, 17, 19, 23 and 10 values of M each: a one-term law, where the whole tower needs four steps.
4. **Lucky primes of the simple part:** Q with λ_g(Q) ≡ −ε·δ (class laws) or λ_g(Q) ≡ 0 (Ono's kind). About 3 primes in
   97 should qualify, against about 1 in 30 million for the whole tower. For the smallest of each kind, G(Q²M) ≡ 0 on
   the class (two M each); a control M′ of the other class gives exactly −2ε·δ·G(M′); and G(Q³m) ≡ 0 for Ono's kind.
Each check of 2 to 4 would hold by chance with probability 1/97.

**Counted (`partitions_97_split.txt`; FLINT, 7 minutes; every prediction held):**
1. **The forms:** the class has pieces of degrees 1 and 3 on lane 97 (pivot coefficients 1, 5, 7, 11), and
   (3|97)·a₉₇(g) ≡ 28. At Q′ = 5: λ_g(5) = 39, the cubic factor is h = x³ + 13x² + 13x + 22, and h(39) = 81 ≠ 0. So
   κ = (49, 26, 44, 1) where (M|5) = +1, and (46, 93, 79, 1) where (M|5) = −1.
2. **Along 97: G(97²M) ≡ 28·G(M) for 20 of 20 values of M = 97·(24n + 23), and G(M) is not 0 for any of them.**
   - The partition numbers used are below 3·10¹¹: levels 1 and 3 only.
   - For the whole tower, c(97²M)/c(M) takes 17 different values over the same M.
   - So the 97-adic constant's first digit is now seen on partition numbers themselves. Yesterday's law for the whole
     tower would have needed level 9.
3. **Along other primes: the one-term law holds for 60 of 60** (Q = 7, 11, 13, 17, 19, 23, ten values of M each).
   - λ_g = 80, 87, 8, 69, 41, 50 there.
   - For the whole tower, c(Q²M)/c(M) takes 2 to 7 different values within one class of M.
4. **The simple part's lucky primes below 5000:**
   - Counted: 6, 4 and 4 of the three kinds, against 6.9 expected for each (derived).
     - Class ε = +1: Q = 53, 607, 1567, 2311, 2659, 3517.
     - Class ε = −1: Q = 1361, 1657, 2789, 2999.
     - Ono's kind: Q = 263, 1193, 1439, 1637.
   - Checked for the two smallest of each kind: every predicted 0 is 0 (12 of 12). All 6 controls take exactly the
     predicted value (−2εδ·G(M′) for the class kinds, ε·δ·G(M′) for Ono's kind).
   - The largest partition number used: p(3,324,019,351,171,224).

**As a congruence** (derived from 4; the case Q = 53, checked above for two values of N):
κ₀·p(N) + κ₁·p(25N − 1) + κ₂·p(625N − 26) + p(15625N − 651) ≡ 0 (mod 97)
- **Where it holds:** every N with 24N − 1 = 53²·M, 97 | M, (M|53) = +1 and 5 ∤ M.
- **The coefficients:** κ as in 1, chosen by (M|5).
- **Where the arguments come from:** they are (5^(2j)·(24N − 1) + 1)/24.

**In plain words.** 97's tower is two towers on top of each other.
- **The simple part:** one form with one constant.
- **The braided part:** three forms whose constants are conjugate in a cubic field.
- **Taking them apart:** one step at 5, whose effect on each form is known, cancels the braided part. What is left
  obeys one-term laws: two levels deeper along 97 multiply it by 28, and a step at any other prime multiplies it by one
  number.
- **Its lucky primes:** about 3 primes in 97 give zeros (53 is the first). The whole tower needs about 1 in 30 million.

**Next: the second digit, 26 (PREDICTED before the runs).** α = 28 + 26·97 + … (mod 97²: 2550).
- Testing it needs the step numbers at 5 to two digits too. They are the weight −2 member's, carried along the family
  like the constants (`partitions_flint_family_steps.py`, a copy of the family script that also keeps T₅).
- Predicted for that run (derived): on lane 97 the polynomial of λ(5) is the base's, [15, 0, 88, 71, 1]. Its root on
  lane 97 starts with 39, the rest is [22, 13, 13, 1], and the constants' digits are yesterday's (28, 26, 81, 49).
- Then, with κ to two digits, at levels 3 → 5 (M = 97³·(24n + 23), n = 0 … 9; numbers up to about 10¹⁵):
  G(97²M) ≡ 2550·G(M) (mod 97²) for all 10.
- At levels 1 → 3 only the first digit should hold. Level b carries about b digits (the Hida steps of 12 Nov), so
  most of these fail modulo 97².
- Also measured, not predicted: the multiplier modulo 97³ at levels 3 → 5 (the third digit 81 is on the edge of what
  level 3 carries).

**Counted (`partitions_flint_family_steps_97.txt`, 13 s; `partitions_97_split_digit2.txt`, 11 minutes; as PREDICTED):**
- **The forms, to four digits:** λ(5)'s polynomial on lane 97 is the base's, [15, 0, 88, 71, 1]. Its root is
  λ_g(5) = 39, 71, 61, 45 … (base-97 digits, lowest first), and the rest is [22, 13, 13, 1] on lane 97. The constants'
  digits are yesterday's: 28, 26, 81, 49.
- **Levels 3 → 5:** the multiplier G(97²M)/G(M), read off exact partition numbers (up to p(1,604,535,581,874,349)), has
  the digits 28, 26 for 10 of 10 values of n. So the 97-adic constant's second digit, 26, is now seen on partition
  numbers.
  - The third digit is 81 for only 1 of 10 (55, 92, 8, 61, 1, 55, 7, 7, 81, 91): level 3 does not reach it.
- **Levels 1 → 3:** the first digit 28 for 10 of 10, and the second digit 26 for 0 of 10, as predicted. Level b reaches
  one digit fewer than the level b + 2 above it, as the twins showed on lane ℓ² from level 3 (the thirteenth time).

## The twins taken apart: each constant alone, digit by digit (13 Nov 2026; `partitions_twins_split`)

**The idea.** At 37, 41, 43 and 53 the class is two forms f₁, f₂, each with its own ℓ-adic constant (the thirteenth time
saw only their sum and product on partition numbers; "one constant alone fails" on the whole tower). One step at 5
cancels either form: G_i(M) = c(25M) + (κ − λ_j(5))·c(M), with κ = (M|5)·(−3|5)·5⁻² and j the other twin. So G_i is
f_i's tower alone, and G_i(ℓ²M)/G_i(M) should be α_i.

**PREDICTED before the runs:**
1. Forms (`partitions_flint_family_steps.py`, 7 digits): the constants' digits agree with the table of the thirteenth
   time (37: 15, 4, 4, 3, 1, 10, 35 and 34, 30, 35, 25, 2, 25, 28; 41: 23, 14, 7, 38, … and 31, 11, 10, 38, …; 43: 31,
   13, 8, 5, 40 and 6, 20, 10, 25, 3). λ(5)'s polynomial on lane ℓ is the base's.
2. Partition numbers (levels b → b + 2, b = 1, 3, 5, 7 as far as p(N) with N ≤ 10¹⁶ allows):
   - Each twin's multiplier carries its own constant's digits, not the other's.
   - The digits that agree grow with the level, at least 1, 2, 3, 4 at b = 1, 3, 5, 7 (97 gave 1 and 2 at b = 1, 3).

**Counted (`partitions_flint_family_steps_37/41/43/53.txt`, seconds each; `partitions_twins_split.txt`, 20 minutes, up to
p(11,236,275,419,782,699); as PREDICTED):**
1. **The forms:** the constants' digits are those of the thirteenth time for 37, 41 and 43 (now to 8 digits). 53 is new:
   43, 26, 7, 26, 45, 51, 7, 34 … and 7, 8, 0, 38, 7, 49, 30, 5 ….
2. **Each twin alone:** the multiplier G_i(ℓ²M)/G_i(M), read off exact partition numbers.

   | levels | 37 | 41 | 43 | 53 |
   |---|---|---|---|---|
   | 1 → 3 | 1 | 1 | 1 | 1 |
   | 3 → 5 | 2 | 2 | 2 | 2 |
   | 5 → 7 | 3 | 3 | 3 | 3 |
   | 7 → 9 | 4 | 4 | (past 10¹⁶) | (past 10¹⁶) |

   - **The rule:** at levels b → b + 2 every twin agrees with its own constant in at least (b + 1)/2 digits, for every
     value of M (2 to 6 per level; at 37 and 41 one value per level is skipped, because G(M) ≡ 0 on lane ℓ). Sometimes there is one digit
     more, at particular M (the same n at every level).
   - **The other twin's constant:** 0 digits agree, every time (147 multipliers).
   - **The whole tower:** c(ℓ²M)/c(M) agrees with neither constant (0 digits, or 1 by chance).
3. So each twin's constant is now seen alone on partition numbers: 4 digits for 37 and 41, 3 for 43 and 53.
   - 37: 15, 4, 4, 3 and 34, 30, 35, 25.
   - 41: 23, 14, 7, 38 and 31, 11, 10, 38.
   - Before, only their sum and product were seen.

**In plain words.** Each of these towers is two simple towers added together, one per form. A single step at 5 removes
either one, and what is left grows by a fixed number, its own constant, every two levels. Two levels deeper show one
more digit of it.


## Laws at the other primes with three forms: 67, 71, 79, 83 (13 Nov 2026; `partitions_flint_lucky_dual`, checks)

**PREDICTED before the run (derived):**
- Every prime Q < 5·10⁷, with the lighter test. Each kind of law has about 3,001,131/ℓ³ expected: about 10.0 at 67,
  8.4 at 71, 6.1 at 79 and 5.2 at 83, about 89 in all. The exact sums are printed with the counts.
- The counts fall within chance of them.
- The smallest class law of each prime holds on exact partition numbers (two members ≡ 0, the control exactly
  2ε·(−3|Q)·Q⁻²·p((M′ + 1)/24)), where its members stay below N = 3·10¹⁶.

**Counted (`partitions_flint_lucky_dual_67_71_79_83_max50000000.txt`, 9 minutes):**

| ℓ | the class on lane ℓ | class ε = +1 | class ε = −1 | Ono's kind | expected, each kind |
|---|---|---|---|---|---|
| 67 | 1 and a pair | 11 | 10 | 7 | 10.0 |
| 71 | one piece of degree 3 | 8 | 12 | 8 | 8.4 |
| 79 | 1 and a pair | 10 | 5 | 3 | 6.1 |
| 83 | 1 and a pair | 6 | 6 | 6 | 5.2 |

- 92 laws against 89.3 expected. χ² = 7.2 on the 12 counts (derived): within chance.
- With 73 and 61, all six primes with three forms give 152 laws against 151.9 expected.
- The smallest of each prime:
  - 67: 665,279 (ε = −1);
  - 71: 351,929 (Ono's kind);
  - 79: 3,846,229 (Ono's kind);
  - 83: 1,260,437 (ε = +1).
- Every list is in the output.

**On exact partition numbers** (`partitions_flint_lucky_check_67_83.txt`, 28 minutes; as PREDICTED). The smallest class
law of each kind where N stays below 3·10¹⁶:

| ℓ | Q (kind) | both members | control: counted, predicted |
|---|---|---|---|
| 67 | 2,654,833 (ε = +1) | ≡ 0 | 61, 61 |
| 67 | 665,279 (ε = −1) | ≡ 0 | 60, 60 |
| 71 | 7,701,497 (ε = −1) | ≡ 0 | 51, 51 |
| 79 | 5,041,271 (ε = +1) | ≡ 0 | 1, 1 |
| 79 | 8,240,447 (ε = −1) | ≡ 0 | 18, 18 |
| 83 | 1,260,437 (ε = +1) | ≡ 0 | 16, 16 |

- The largest number used was p(25,442,829,664,257,819).
- 71's first ε = +1 law and 83's first ε = −1 law need N above 3·10¹⁶; not computed.
- **Simplest example found:** p(19728144816847346813·n + 6177904563614) ≡ 0 (mod 67) for every n (modulus 67·665279³;
  the first member, p(6,177,904,563,614), was computed above).

## The braided forms taken apart, with numbers from their own field (13 Nov 2026; `partitions_braid_split`)

**The idea.** A conjugate group (a pair at 47, 59 and 61; a triple at 73 and 97) cannot be split by steps with ℓ-adic
numbers: whatever cancels one member cancels its conjugates too. But their constants live in a bigger number field, the
unramified extension of degree 2 or 3, and steps with numbers from that field can split them.
- Write x for a root of the group's factor A of the constants' polynomial. So the constant of f₁ is x itself, and
  its conjugates are the other roots.
- λ₁(5) is the root of λ(5)'s polynomial paired with x (on lane ℓ through T₅ = R(T_ℓ) on the class), lifted by Newton's
  method.
- With h(y) = (λ(5)'s polynomial)/(y − λ₁), G₁(M) = Σ κ_j·c(5^(2j)·M) is f₁'s tower alone. Its κ_j are numbers
  a + b·x (+ c·x²) modulo ℓ^K.

**PREDICTED before the run (derived):**
- At levels b → b + 2 the multiplier G₁(ℓ²M)/G₁(M), a number in that field, agrees with x in at least (b + 1)/2
  digits, and with each conjugate of x in none.
- The other members need no separate check. Partition numbers are whole numbers, so their multipliers are the
  conjugates of f₁'s.
- Levels as far as p(N) with N ≤ 1.2·10¹⁶ allows. The single forms (61's and 97's) are checked the same way, as a
  control of the code against the twins' run.

**Counted (`partitions_flint_family_steps_47/59/61/73.txt`; `partitions_braid_split.txt`, 35 minutes, up to
p(10,721,336,237,467,474); every prediction held):**

| group | its field | levels 1 → 3 | 3 → 5 | 5 → 7 | 7 → 9 | digits shared with a conjugate |
|---|---|---|---|---|---|---|
| 47, a conjugate pair | degree 2 | 1 (6 of 6) | 2 (6 of 6) | 3 (6 of 6) | 4 (1 of 1) | 0, every time |
| 59, a conjugate pair | degree 2 | 1 | 2 | 3 | (past 10¹⁶) | 0 |
| 61, a conjugate pair | degree 2 | 1 | 2 | 3 | (past 10¹⁶) | 0 |
| 73, the conjugate triple | degree 3 | 1 | 2 | 3 (1 of 1) | (past 10¹⁶) | 0, against both conjugates |
| 97, the conjugate triple | degree 3 | 1 | 2 | (past 10¹⁶) | | 0, against both conjugates |
| 61, the single form | ℓ-adic | 1 | 2 | 3 (one value 4) | (past 10¹⁶) | none to compare |
| 97, the single form g | ℓ-adic | 1 | 2 | (past 10¹⁶) | | none to compare |

- The table gives the digits each group's multiplier shares with its own constant x, at least, for every value of M.
- "Digits" are digits in the group's own field, each one an element of the field with ℓ² or ℓ³ elements.
- **108 multipliers in all.** Every one agrees with its own constant in at least (b + 1)/2 digits at levels b → b + 2,
  the twins' rule.
- **The conjugates:** 79 multipliers belong to conjugate groups, and none shares even one digit with a conjugate.
- **The pairing is right everywhere.** λ₁(5) from R(x) is a root of λ(5)'s polynomial on lane ℓ (7 of 7 groups).
- **The code agrees with today's earlier runs.** The single form of 97 gives 28, 26 again, and 61's single form gives
  the twins' rule.
- **NEXT item 23, answered:** a conjugate group's members can be told apart on partition numbers, with numbers from
  their own field.

**In plain words.** Some forms come braided: their constants are not ordinary numbers, but numbers in a bigger number
system (as √2 is for the whole numbers). A filter made of ordinary numbers cannot pull them apart. A filter made of
numbers from that bigger system can, and then each braided constant shows up alone in the partition numbers, digit by
digit, like the others.
- **At 73:** all three of its forms, each in the cubic field.
- **At 97:** its triple.
- **At 47, 59 and 61:** their pairs.



## Pictures of the chords: the data behind the page (13 Nov 2026; `chords_page_data`, `partitions_flint_lucky_dual` 97)

**PREDICTED, written down right after the 97 scan was started and before any of its output existed (derived):**
1. **97's whole tower, every prime Q < 5·10⁷:** about 0.03 laws of each kind. A law needs all four forms at once,
   (1/97)·(1/97³) per prime, so most likely 0.
2. **The chord is silent only where every note is silent.** At 37, a prime Q is lucky for one twin alone when
   (3|Q)·a_Q(f_i) takes the law's number. The whole tower's laws must be exactly the Q where both twins are lucky
   with the same kind, the lists of the nineteenth time (57, 63, 54 up to 10⁶). Each twin alone should be lucky at
   about 3 primes in 37.
3. **On 37's cycle** (the angle of a number is its power of a primitive root, as in the page "Every prime at zero"):
   - a note's values turn by the same angle, its constant's, from level b to b + 2, for every n;
   - the whole tower's values do not.
4. **On the cycle of the field with 73³ elements:** each braided strand of 73 turns by its own angle. The three
   strands' angles are θ, 73θ and 73²θ (Frobenius).

**Counted.**
- Records: `partitions_flint_lucky_dual_97_max50000000.txt` (5 minutes) and `chords_page_data.txt` (13 minutes, up to
  p(11,236,275,419,782,699)).
- Predictions 1 to 4 held.

1. **97's whole tower:** 0 laws of each kind among 3,001,131 primes, against 0.03 expected each.
2. **37's silences** (primes below 10⁶; 6,365.9 expected for each twin, derived):
   - Twin A is silent at 6,363 primes (z = −0.04), twin B at 6,411 (z = +0.57).
   - The whole tower's laws are exactly the primes where both twins are silent with the same kind: 57, 63 and 54, the
     lists of the nineteenth time, all three the same.
3. **On 37's cycle** (primitive root 2), from level b to b + 2:
   - The whole tower's values turn by 15, 14 and 13 different angles at b = 1, 3, 5 (19 values of n), and by 3
     different angles at b = 7 (3 values).
   - Each twin turns by exactly one angle at every level: 13 for twin A (15 = 2¹³) and 8 for twin B (34 = 2⁸).
4. **73's braid:** strand 1 turns by θ = 95,533 of 389,016 for 13 of 13 values of n. The strands' turns are θ, 73θ and
   73²θ = 95,533, 360,637 and 262,429.
5. **97's single note** (not predicted beyond "about 3 in 97"):
   - Silent at 2,326 of 78,495 primes, against 2,427.8 expected (z = −2.07).
   - By kind: z = −0.85, −1.06 and −1.66. A little low; tested next.

The page: "Towers as Chords", https://claude.ai/artifact/JfyuxTXWjw6jZxAx9tHdi9 (private). It is built from
`chords_page_data.json` and from the 97 scan.

**The full forms (version 2 of the page; `chords_page_data2.py`, `.txt`, `.json`, 2 minutes; COUNTED).** The user asked
for the points completed round the full circle.
- **37, all 288 values of n from 0 to 359** (5 ∤ M), levels 1, 3, 5:
  - each twin turns by one angle, 13 and 8, every time;
  - the whole tower turns by all 36 angles there are.
- **73, 240 values of n:** strand 1 turns by one angle for 240 of 240.
- **Each note's star.** A note followed two levels at a time traces a star on its cycle (the field with ℓ^e elements for
  a braided group), with as many points as the multiplicative order of its constant.
  - The least common multiple of a tower's stars equals the loop measured on 12 Nov (`partitions_loops.txt`) for all 13
    primes from 13 to 61.
  - Examples: 31's note is a 2-point line (constant −1); 23's constant 5 is a primitive root (a 22-gon); 37's notes have
    36 and 9 points; 47's braided pair has 2,208 points (every place of its field's cycle).
  - New: 73's tower loops after 389,016 two-level steps on lane 73, every place of the cycle of its field.
  - New: 97's tower loops after 304,224 = lcm(32, 76,056): its single note has 32 points, its triple 76,056.
- **Corrected (13 Nov, after the user's questions about 31, 61 and 23): what a star shows.**
  - **The note's own:** the number of points, and their even spacing (a note of order K visits K evenly spaced places,
    whatever number walks the cycle).
  - **Our choice:** which point joins which. With walking number g, the step is θ_g = log_g(α), and another primitive
    root g′ = g^u turns the step into θ_g·u⁻¹ (mod K).
  - **Artifacts from that choice:**
    - 23's "plain 22-gon": its constant 5 is the smallest primitive root itself.
    - 61's 12-gon round the rim: 32 = 2⁵ and the walk is by 2.
    - The braids at 47, 59 and 61, drawn round the rim: the program walked the big field by the braid's own constant, so
      every step is one place.
  - **Not an artifact:** 31's line, since −1 is half a turn for every walk.
  - **The fix:** version 4 of the page draws the points themselves, completes the paths up to 6,000 steps, and has a
    button that joins the same points in other orders.

**PREDICTED before the test of 97's note (derived).** If the single note's Galois image modulo 97 is the full one, as
the twins' images at 37 are:
- the values of λ_g(Q) = (3|Q)·a_Q(g) for all primes Q < 10⁷ spread over the 97 residues as Chebotarev's weights say
  (χ² near 96 on 96 degrees of freedom);
- the silences up to 10⁷ come within about 2 standard deviations of expected, so the low count up to 10⁶ was chance.

A clearly smaller image (dihedral, or reducible) would show a χ² in the thousands.

**Counted (`note97_check.txt`, 5 minutes; as PREDICTED: the low count was chance):**

| note | primes | spread χ² (degrees of freedom) | silences: counted, expected (z) |
|---|---|---|---|
| 97's g | below 10⁶ | 107.6 (96) | 2,326, 2,427.8 (−2.07) |
| 97's g | below 10⁷ | 97.0 (96) | 20,827, 20,554.5 (+1.90) |
| 37's twin A | below 10⁷ | 48.6 (36) | 53,866, 53,896.7 (−0.13) |
| 37's twin B | below 10⁷ | 40.5 (36) | 53,765, 53,896.7 (−0.57) |

- The values spread over the residues as Chebotarev's weights say, so the notes' Galois images behave as full ones.
- The low count below 10⁶ turned into a slight excess by 10⁷. Both are within about 2 standard deviations.


## The speed of the towers: how many digits two levels add (13 Nov 2026; `tower_speed`, then the forms)

**The question.** At 37 to 97 every note's multiplier between levels b and b + 2 shares exactly (b + 1)/2 digits with its
constant: two levels add one digit. The fourth time measured 5, 3, 1 and 2 digits per two levels at 5, 7, 11 and 13
(NEXT item 5: at 5 the speed equals the gap between the law's slope and the next slope, at 7 and 11 half of it, and why
half was open). What sets the speed?

**PREDICTED before the run:**
- 13 adds 2 digits per two levels, as the fourth time measured.
- 17, 19, 23, 29 and 31 are measured, not predicted. If one rule held from 13 to 97, all of them would add 1, like 37 to
  97; 13's 2 says that is not the whole story.

**Counted (`tower_speed.txt`, 23 minutes, up to p(2,898,129,680,321,528)).** The digits the multiplier shares with α at
levels b → b + 2, smallest over n = 0 … 7:

| ℓ | b = 1 | 3 | 5 | 7 | 9 | 11 | the rule |
|---|---|---|---|---|---|---|---|
| 13 | 1 | 3 | 5 | 7 | 9 | 11 | b (one digit per level) |
| 17 | 1 | 3 | 5 | 7 | 9 | 11 | b |
| 19 | 1 | 3 | 5 | 7 | 9 | | b |
| 23 | 1 | 3 | 5 | 7 | 9 | | b |
| 29 | 1 | 2 | 3 | 4 | | | (b + 1)/2 (one digit per two levels) |
| 31 | 1 | 2 | 3 | 4 | | | (b + 1)/2 |
| 37 … 97, each note | 1 | 2 | 3 | 4 | | | (b + 1)/2 (twins and braid runs) |

- **13 is as predicted**, and 17, 19 and 23 are as fast. 29 and 31 are slow, although their class is also one form, so
  "one form or several" is not the reason.
- **Literature (QUOTED, through Garvan and Shomanov, arXiv 2504.10824, §1).**
  - Atkin and O'Brien's congruence for 13, their own theorem of 1967 (restated there as Theorem 1.1): P(13^(α+2)·N) ≡
    K_α·P(13^α·N) (mod 13^α). Level α carries α digits, our b. *Corrected 14 Nov (a relayed review):* this line first
    said "proved there". Garvan and Shomanov proved a different conjecture of Atkin and O'Brien, the one on Hecke
    operators T_(p²), which Atkin and O'Brien had shown for α = 1, 2.
  - For 5 ≤ ℓ ≤ 31, Folsom, Kent and Ono's Theorem 1.4 holds from level b_ℓ(α) = 2(⌊(ℓ − 1)/12⌋ + 2)α − 3. "Boylan and
    Webb showed that … we can take b_ℓ(α) = 2α − 1", which is our slow rule.
- **So:**
  - 13's speed is Atkin–O'Brien's (rediscovered).
  - 17, 19 and 23 are as fast as 13. That is not in what we read.
  - 29 and 31, and every note from 37 to 97, sit exactly on Boylan–Webb's bound.

**A reading (derived after the data, so not yet a test).** The speed follows the level-1 cusp forms of weight ℓ − 1, the
space where Folsom, Kent and Ono place the tower on lane ℓ.

| ℓ | dimension of S_(ℓ−1)(SL₂(ℤ)) | R_ℓ |
|---|---|---|
| 13, 17, 19, 23 | 1 | 1 |
| 29, 31 | 2 | 1 |
| 37 to 47 | 3 | 2 |

- Where the space has room beyond the class, the extra forms would be the part of the tower that fades, one ℓ per two
  levels. Where it has none (13 to 23), the tower is fast.

**PREDICTED before the test (derived from that reading).** On lane ℓ, the number of eigenforms of S_(ℓ−1)(SL₂(ℤ)) that
are not ordinary (a_ℓ ≡ 0 mod ℓ, counted with multiplicity: the generalized kernel of T_ℓ) is exactly
dim S_(ℓ−1) − R_ℓ, for every prime ℓ from 13 to 400.
- 0 for 13 to 23, so there all forms are ordinary.
- 1 for 29 and 31, 1 for 37 to 47, and so on.
- If this fails for some ℓ, the reading is wrong as stated.

**Counted (`level1_ordinary_check.txt`, 4 s): NOT as predicted; the reading is wrong as stated.**
- For every prime ℓ from 13 to 400 (73 primes), T_ℓ has full rank on S_(ℓ−1)(SL₂(ℤ)) on lane ℓ. Every eigenform there is
  ordinary, so there is no non-ordinary room.
- The prediction held only where it said 0 (13, 17, 19, 23).
- So the extra dimensions from 29 on are not non-ordinary forms, and what makes 13 to 23 fast is open again.
- The split itself stands: 13, 17, 19, 23 fast; 29, 31 and every note from 37 to 97 at Boylan–Webb's rate.

**PREDICTED before the next run (derived; for the page's split circles).** 97's braided triple is silent at Q when
(3|Q)·a_Q of its eigenform, a number of the field with 97³ elements, equals the law's number. By Chebotarev in that
field, about 0.26 such primes are expected below 10⁶, most likely none. So 97's whole tower (all four notes together)
is silent at none of them (`chords_page_data3.py`).

**Counted (`chords_page_data3.txt`, 23 s).**
- 97's braided triple is silent at 1 of the 78,495 primes below 10⁶: Q = 782,147, Ono's kind, against 0.26 expected
  (a chance of about 1 in 4 for one or more).
- 97's single note is not silent there, so the whole tower is silent at none, as predicted. This agrees with the scan
  to 5·10⁷.

## What fades: the speed's own note (13 Nov 2026; `tower_fading`)

**The question.** A tower settles on its constant because the rest of it fades. Write
E_b(n) = H(ℓ²M) − α·H(M) (M = ℓ^b·(24n + k)), the part of level b + 2 not explained by the constant. If what fades is
itself one note, with its own multiplier β per two levels, then E_(b+2)(n) = β·E_b(n) for every n, and v_ℓ(β) is the
speed.

**PREDICTED before the run (derived).**
- **One note:** for each of 13, 17, 19, 23, 29 and 31, the ratio E_(b+2)(n)/E_b(n) is the same number for every n,
  to as many digits as the data carry. If several fading notes of the same size mix, it varies with n instead.
- **Its size:** v_ℓ(β) = 2 at the fast primes (13, 17, 19, 23) and 1 at the slow ones (29, 31).
- **Data:** α to 30 digits (13 from `partitions_flint_family_13_N30.txt`); partition numbers modulo ℓ^30; n = 0 … 7;
  levels as in `tower_speed`.

**Counted (`tower_fading.txt`, 19 minutes, up to p(2,898,129,680,321,528)).**

| ℓ | v_ℓ(β) (the speed) | one fading note? | β's unit part, digits lowest first | digits shared by all n, by level |
|---|---|---|---|---|
| 13 | 2 | yes | 10, 11, 6, 2, 10, 10, … | 0, 3, 4, 6, 8 (levels 1 … 9) |
| 17 | 2 | yes | 3, 9, 1, 13, 4, 3, … | 1, 3, 5, 7 (levels 1 … 7) |
| 19 | 2 (from level 7; 1 to 3 below) | no | — | 0 at every level |
| 23 | 2 | no | — | 0 at every level |
| 29 | 1 | yes | 23, 3, 10, … | 0, 1, 2 (levels 1, 3, 5) |
| 31 | 1 | yes | 2, 28, 13, … | 0, 1, 2 |

- **The size is as PREDICTED at all six.** What fades shrinks by ℓ² every two levels at the fast primes and by ℓ at the
  slow ones, so the speed is the size of the fading part.
- **"One note" is as predicted at 13, 17, 29 and 31, NOT at 19 and 23.** There the ratio differs from n to n even where
  its size is settled, so the fading part is at least two notes of the same size.
- **Open:** what the fading notes are in terms of modular forms. A first guess is the next families after the class,
  slope 2 at 13 and 17 and slope 1 at 29 and 31, read at weight −2. β's digits above are the target to match.

## The braids' angles, and walking each note by its own constant (13 Nov 2026; `braid_angles`)

**The user's questions.** Do the braided strands' turns add up to something, or share a ratio? And can the walking
number of each star be derived from what the note itself does?

**PREDICTED before the run (derived from Frobenius and the norm).** For each braided group (pairs at 47, 59, 61;
triples at 73, 97), with the strands' turns θ_i measured on the big field's cycle:
- **Ratio:** each turn is exactly ℓ times the one before (mod a full turn).
- **Sum:** the turns add up to the angle of the product of the group's constants, an ordinary number mod ℓ, on ℓ's own
  cycle. So the sum is always a multiple of 360°/(ℓ − 1): 5° at 73, 3.75° at 97.
- **Odd or even:** the sum is an odd number of those steps exactly when that product is not a square mod ℓ. At 73 the
  product is 44, not a square (the fact behind "no one-level law" at 73), so odd.
- **Difference:** two neighbouring strands differ by a multiple of 360°/(ℓ² + ℓ + 1) for triples and 360°/(ℓ + 1) for
  pairs.
- **The twins:** at 37, 41, 43 and 53, one constant is a power of the other whenever its order divides the other's.
- **Walking by the note's own constant:** each star becomes a plain K-gon, one place per two levels, K the note's
  order. Only K is the note's own.

**Counted (`braid_angles.txt`; every prediction held, at all five braids).**

| braid | each turn ℓ × the one before | sum of the turns | product of the constants (mod ℓ) | a square? | odd/even steps |
|---|---|---|---|---|---|
| 47, pair | yes | 29 × 360°/46 | 19 | no | odd |
| 59, pair | yes | 26 × 360°/58 | 57 | yes | even |
| 61, pair | yes | 1 × 360°/60 | 18 | no | odd |
| 73, triple | yes | 5 × 5° = 25° | 44 | no | odd |
| 97, triple | yes | 60 × 3.75° = 225° | 64 = 8² | yes | even |

- The angles here use one walking number per field (the first generator in the field's listing other than the braid's
  own constants). The rules and the parities do not depend on that choice.
- Neighbouring strands differ by a multiple of the finer grid every time: 360°/(ℓ + 1) for pairs, 360°/(ℓ² + ℓ + 1)
  for triples.
- **97's triple has a square product (64),** so its constants are squares in their field. 97's tower has no one-level
  law because of its single note alone (28 is not a square), not the triple.
- **The twins as powers of each other:** 34 ≡ 15²⁰ (mod 37), 23 ≡ 31⁷ (mod 41), 6 ≡ 31⁷ (mod 43), 7 ≡ 43³ (mod 53).
- **The page, version 5:** every star can be walked by its own constant (its own clock, a plain polygon) or on the
  first note's clock (37's B jumps 20 ticks of A's clock). A dial per braid shows the strands' turns and their sum.

**Stars drawn as the tower goes (version 6 of the page; counted from `chords_page_data2.json`).** At the user's
suggestion, the stars can be drawn with no walking number at all.
- **How:** every remainder sits at its own place on a clock with ℓ hours, and each value is joined to the value the
  tower takes two levels down. The lines form the note's multiplication table, r → α·r.
- **The real steps at 37** (288 values of n, levels 1 → 3 → 5, from exact partition numbers):
  - note A: 566 of 566 steps lie on the ×15 table;
  - note B: 552 of 552 lie on the ×34 table;
  - the whole tower: 24 of its 558 steps lie on the ×15 table.
- **So** "joined as the tower goes" and "the multiplication table of the note's constant" are the same picture.
  Braided notes have no place on ℓ's clock (their values are numbers of a bigger field).

## The braids in size and direction (13 Nov 2026; `braid_polar_data`; the user's idea: layers as circles getting larger or smaller)

**The idea.** A number v of a braid's field has a size and a direction, as a complex number does.
- **Size:** its norm, the product of its strands' values, an ordinary number mod ℓ.
- **Direction:** v^(ℓ − 1), on the norm-one circle of (ℓ^e − 1)/(ℓ − 1) places.
- Two levels down (v ↦ α·v) the size is multiplied by the product of the braid's constants, and the direction turns by
  a fixed step. Drawn with circles for sizes and angles for directions, the braid needs no walking number for its
  sizes. Its directions are counted on its own clock where its constant's direction runs through every place.

**PREDICTED before the run (derived).** At 73, for every real tower value of strand 1, from exact partition numbers:
- the size at level 3 is 44 × the size at level 1;
- the direction turns by exactly one tick of the 5,403-place circle.

**Counted (`braid_polar_data.txt`, 20 s; as PREDICTED).**
- **73, the real tower** (strand 1, levels 1 → 3, 240 values of n, from exact partition numbers): the size is
  multiplied by 44 for 240 of 240, and the direction turns by exactly one tick of the 5,403-place circle for 240 of 240.
- **The others:**

  | braid | size multiplied by | direction turns by |
  |---|---|---|
  | 47 | 19 | 1 tick of 48 (own clock) |
  | 61 | 18 | 1 tick of 62 (own clock) |
  | 59 | 57 | 4 of 60 places (its direction has order 15, so it is counted by γ^58) |
  | 97's triple | 64 | 5,010 of 9,507 places (its direction has order 3,169) |

- **A braid is an ordinary note in its size and a pure turning in its direction.** The size note's constant is the
  product of the braid's constants: the same number the strands' angles add up to.
- **The page, version 7:**
  - "Before and now on 37's cycle" now places every number at its own hour of a 37-hour clock; the turning view is a
    button.
  - The animation follows the rule: each dot travels to where the note's constant sends it.
  - A new figure draws the braids in size and direction.

## Where the fading notes live (13 Nov 2026; `fading_forms`, a copy of `partitions_flint_slopes`)

**PREDICTED before the run (a guess to test, derived from the slopes picture of NEXT item 5).** The fading note is
the class's own steeper neighbour: a family congruent to the class form on lane ℓ, new at 2 and 3 with the class's
signs, of slope 2 at the fast primes and 1 at the slow ones.
- If so, at every weight −2 + ℓ(ℓ − 1)(x + 1), the class's piece has exactly one root of T_ℓ of that slope.
- Read at weight −2, (3|ℓ)·u(−2)/ℓ^S agrees with the measured unit part: 13 → 10, 11, 6, 2, 10, 10; 17 → 3, 9, 1, 13,
  4, 3; 29 → 23, 3, 10; 31 → 2, 28, 13.
- If there is no such root, or the digits disagree, the fading note lives elsewhere (another system or another sign).

**Counted (`fading_forms_13_s2.txt`): NOT as predicted.** At weight 154 the class's piece (new at 2 and 3, the class's
signs) has two forms, with slopes 0 and 9: there is no slope-2 root there. So the fading note is not the class's own
steeper neighbour.

**PREDICTED before the next run (widened).** The fading part is already divisible by 13, so it can belong to any system
on lane 13. The fading note is then a slope-2 family among all the forms new at 2 and 3 with the class's signs, whatever
their numbers on lane 13 (`fading_forms.py --any`). If there is exactly one slope-2 root at each weight, its
u(−2)/13² agrees with 10, 11, 6, 2, 10, 10.

**Counted (`fading_forms_13_s2_any.txt`, then `fading_forms_13_s2_any_N8.txt` with the digits printed).**
- At weight 154 the 7 forms new at 2 and 3 with the class's signs have slopes 6 (three), 4 (two), 2 (one) and 0 (the
  class). There is exactly one slope-2 root: (3|13)·u/13² = 10, 9, 9, 2, 9, 8, 12, 11, …
- At weight 310 there is again exactly one: 10, 7, 9 (3 digits known).
- At weight 466 it merges into a crowd of 19 roots of average slope 28/19, so the series stops at two weights. [Corrected 13 Nov: not a crowd but a precision artifact; see "What fades at 19 and 23".]
- Carried to weight −2 (two weights reach two digits): 10, 11. The fading note measured from partition numbers: 10, 11,
  6, 2, …
- **The first two digits agree.** A chance match would happen 1 time in 169. So 13's fading note looks like the slope-2
  form among those with the class's signs, a form not congruent to the class mod 13. Settling more digits needs more
  weights, beyond 466.
- **17** (`fading_forms_17_s2_any.txt`):
  - At weight 270 the 12 forms with the class's signs have slopes 13/5 (ten), 2 (one) and 0 (the class). The single
    slope-2 root: (3|17)·u/17² = 3, 13, 6, 12, 14, 5, …
  - The measured fading note: 3, 9, 1, 13, … The first digit agrees, which is all one weight carries (at 13 the single
    weight 154 also agreed in its first digit only).
  - At weight 542 the root merges into a crowd (9/7 × 21). [Corrected 13 Nov: a precision artifact; see "What fades at 19 and 23".]
- **Together:** the fading notes at 13 and 17 agree with the slope-2 form among those with the class's signs in 3
  digits out of 3 tested. A chance match would happen 1 time in 2,873.
- **29 and 31** (`fading_forms_29_s1_any.txt`, `fading_forms_31_s1_any.txt`; slope 1 wanted): at weights 810 and 928
  [corrected 13 Nov: a precision artifact, see "What fades at 19 and 23"] the forms with the class's signs give one long Newton segment besides the class, 33 roots of average slope 20/33
  at 29 and 38 of average slope 10/19 at 31. No single slope-1 root can be isolated there, so this is inconclusive, not
  a failure.
- **So far:** the fast primes' fading notes match a slope-2 form with the class's signs, one that is not congruent to
  the class. Settling more digits, and the slow primes, needs weights closer to −2 than these.

**The strands of each braid, separately (computed from `braid_polar_data.json` for the page "Notes on a Donut",
version 2).**
- **Shared size:** strand i of a braid has the same size steps as strand 1 (×N(α): the size is the product of all
  strands).
- **Its own direction:** its direction turns by step·ℓ^i (mod the circle's places).
- **The steps:** 47: 1 and 47 ≡ −1 (of 48); 59: 4 and 56 ≡ −4 (of 60); 61: 1 and 61 ≡ −1 (of 62); 73: 1, 73 and
  5,329 ≡ −74 (of 5,403); 97's triple: 5,010, 1,113 and 3,384 (of 9,507).
- **They always add up to whole turns** (derived: step·(1 + ℓ + … + ℓ^(e−1)) ≡ 0). A pair's two strands turn in
  exactly opposite directions.

## The fading form: which family, and more digits (13 Nov 2026; `fading_forms2`)

**PREDICTED before the run.**
1. **Which family (asked, not predicted):** the slope-2 form at 13 and 17 has its own Hecke numbers a_q on lane ℓ (q =
   5, 7, 11, 13). Is it the class form twisted by a power of q (a_q ≡ q^j·a_q(class) for one j, a θ-twist), or
   something unrelated?
   - The θ-twists that keep the weight class are j = 0 (the class itself, already excluded by its slope) and
     j = (ℓ − 1)/2, the quadratic twist.
   - Forms of that kind have slope (ℓ − 1)/2 = 6 at 13, not 2.
   - So no θ-twist is expected.
2. **More digits:** skip the weights where the slope-2 root is not alone, and carry the rest to weight −2 (Lagrange
   interpolation in x). 13's fading constant then agrees with 10, 11, 6, 2, … in as many digits as the weights carry.

**Counted (`fading_forms2_13_s2_any.txt`, `fading_forms2_17_s2_any.txt`).**
- **Not a θ-twist, as predicted.** The slope-2 form's Hecke numbers on lane ℓ:

  | ℓ | the slope-2 form | the class form |
  |---|---|---|
  | 13 | a₅ = 11, a₇ = 8, a₁₁ = 0 | 3, 5, 5 |
  | 17 | a₅ = 4, a₇ = 2, a₁₁ = 9, a₁₃ = 7 | 16, 12, 15, 12 |

  No power of q turns the class's numbers into these, so the fading note is a form of another family.
- **No more digits this way.** At the higher weights (13: 466 to 1,090; 17: 542 to 1,358), all forms with the class's
  signs fall into one long Newton segment (average slopes 1/2 to 28/19; [corrected 13 Nov: a precision artifact, see "What fades at 19 and 23"]), so the slope-2 root is never alone again;
  still 2 digits (13) and 1 (17).

**PREDICTED before the next run.** Taking the piece with the slope-2 form's own Hecke numbers (`--sys`) at every
weight, the slope-2 root stands alone at each one. Carried to weight −2, it agrees with the measured fading constants
(13: 10, 11, 6, 2, 10, 10; 17: 3, 9, 1, 13, 4, 3) in as many digits as the weights carry.

## 67, 71, 79 and 83 on the pages (13 Nov 2026; `partitions_flint_family_steps_67/71/79/83.txt`, `chords_page_data4`)

The user noticed that these four were missing from both pages: their laws had been counted, but their constants had
never been computed. They are now, from the forms (counted):

| ℓ | notes on lane ℓ | the whole tower's loop | braid: product of its constants | sum of the turns |
|---|---|---|---|---|
| 67 | constant 41 (66 points) and a braided pair (2,244 points) | 2,244 | 65, a square | 16 steps, even |
| 71 | a braided triple (357,910 points, the field's whole cycle) | 357,910 | 69, not a square | 51 steps, odd |
| 79 | constant 18 (13 points) and a braided pair (3,120 points) | 3,120 | 51, a square | 4 steps, even |
| 83 | constant 16 (41 points) and a braided pair (6,888 points, the whole cycle) | 6,888 | 79, not a square | 1 step, odd |

- Every braid's turns grow ×ℓ strand to strand, and the sum is the product's angle with the predicted parity (4 of 4).
- The sizes are multiplied by 65, 69, 51 and 79. The directions turn 46 of 68, 1 of 5,113 (own clock), 42 of 80, and
  1 of 84 (own clock).

## The fading form, found (13 Nov 2026; continued from "The fading form: which family")

**Counted (`fading_forms2_13_s2_sys.txt`, `fading_forms2_17_s2_sys.txt`, `fading_carry.txt`): as PREDICTED. The fading
note is found.**
- **Alone at every weight.** In the piece with the slope-2 form's own Hecke numbers, the slope-2 root stands alone at
  every weight: at 13, weights 154 to 1,090 (2 to 13 forms, the rest of slope 6 or about 26/9); at 17, weights 270 to
  1,358.
- **Carried to weight −2** through the weights known to 8 digits or more (Lagrange in x; `fading_carry.py`):

  | ℓ | (3\|ℓ)·u_g(−2)/ℓ² from the forms | the fading constant from partition numbers | digits that agree |
  |---|---|---|---|
  | 13 | 10, 11, 6, 2, 10, 10 | 10, 11, 6, 2, 10, 10 | 6 of 6 |
  | 17 | 3, 9, 1, 13, 4, 3 | 3, 9, 1, 13, 4, 3 | 6 of 6 |

  - Each added weight adds about two agreeing digits.
  - The script's own carry-back had used every weight, including those known to only 4, 2 or 0 digits, and so stopped
    at 2 digits.
  - A chance agreement would happen about 1 time in 10¹⁴.
- **So** what fades in the towers of 13 and 17 is the weight −2 member of a slope-2 family of level-6 forms, with the
  class's Atkin–Lehner signs, in another system on lane ℓ. Its multiplier per two levels is (3|ℓ)·u_g(−2), the same
  formula as the constants, now for the next family. Its size, ℓ², is the towers' speed.
- **17's run, finished with six weights** (270 to 1,630; `fading_forms2_17_s2_sys.txt`): its own carry-back through
  all six already gives 3, 9, 1, 13, 4, 3, which is 6 of 6 digits, the same as `fading_carry.txt`.

## What fades at 19 and 23: two quiet notes at once (13 Nov 2026; `fading_two`, `fading_forms3`)

**The question.** At 19 and 23 the fading part is not one note: its ratio from level to level changes with n
(`tower_fading`). If it is exactly two notes, with multipliers β₁ and β₂ per two levels, then every n obeys one
two-step rule: E_(b+4)(n) = s·E_(b+2)(n) − p·E_b(n), with s = β₁ + β₂ and p = β₁·β₂.

**PREDICTED before the runs (derived from what was found at 13 and 17).**
1. **From partition numbers alone (`fading_two`):** one s and one p fit every n, at 19 and at 23, to more digits the
   deeper the levels. Both notes have size ℓ², so v_ℓ(p) = 4 and v_ℓ(s) ≥ 2.
2. **From the forms alone (`fading_forms3 --any`):** at the family's first weight (19: 340; 23: 504), the forms new at
   2 and 3 with the class's signs have exactly two roots of slope 2 besides the class (13 and 17 had one).
3. **The two agree:** carried to weight −2, the two slope-2 forms give s = (3|ℓ)·(u₁ + u₂) and p = u₁·u₂ with the same
   digits as the partition numbers.
- If the two slope-2 forms are conjugate (a braided pair), the quadratic they share has no root on lane ℓ; the test of
  s and p is the same.
- If no single s and p fit (three notes or more), or the forms show one or three slope-2 roots, the picture of 13 and
  17 does not carry over as it stands.

**Counted (`fading_two.txt`; partition numbers up to p(2,898,129,680,321,528), 6 minutes): as PREDICTED (1).**
- **One rule fits every n.** At both primes and at every level triple, all 28 pairs of n (6 at 23's deepest triple)
  give the same s and p. None fails to fit.
- **The shared digits grow with depth:** at 19, s shares 2, 5, 6 digits and p 5, 6, 8 (levels 1–5, 3–7, 5–9); at
  23, s 2, 4, 6 and p 4, 6, 8. v(s) = 2 and v(p) = 4 at both, as predicted.
- **The rule then holds at every n** to every digit the data carry (for example ℓ¹³ of ℓ¹³ at levels 5–9).
- **The two notes are separate**, both on lane ℓ (the discriminant is a nonzero square there), not a braided pair:

  | ℓ | the two fading notes, β/ℓ² = (3\|ℓ)·u(−2)/ℓ², digits lowest first |
  |---|---|
  | 19 | 2, 1, 0, 17, …  and  13, 14, 2, 15, … |
  | 23 | 8, 20, 21, 8, …  and  17, 20, 19, 21, … |

**Counted (`fading_forms3_19_s2_any.txt`, `fading_forms3_23_s2_any.txt`, first weight).**
- **19, weight 340: as PREDICTED (2).** The 15 forms new at 2 and 3 with the class's signs have exactly two roots of
  slope 2 besides the class. Their first digits, (3|19)·u/19² ≡ 13 and 2, are the first digits of the two notes above.
  They lie in two different systems on lane 19 (a₅, a₇, a₁₁, a₁₃ ≡ 0, 7, 16, 8 and 16, 16, 18, 10), and neither is a
  θ-twist of the class.
- **23, weight 504: at first NOT as predicted.** The polygon showed one long segment, 21 roots of slope 38/21, and no
  slope-2 root.

**Correction (counted): the "crowds" were the working precision, not the forms.**
- In every "crowd" so far, the slopes add up to exactly the precision of that run:
  - 23: 38/21 × 21 = 38, at 23³⁸;
  - 13 and 17 above weight 154 and 270: for example 9/7 × 21 = 27 at 17²⁷, and 28/19 × 19 = 28 at 13²⁸;
  - 29 and 31: 20/33 × 33 = 20 and 10/19 × 38 = 20, at ℓ²⁰.
- So the last coefficient of T_ℓ's polynomial was 0 to every kept digit. The polygon's left end was then set by the
  precision, which swallowed the slope-2 (or slope-1) roots into one segment.
- The earlier readings "it merges into a crowd" (13, 17) and "one long Newton segment" (29, 31) were wrong: the roots
  were hidden, not merged. That is why `--sys` later found the slope-2 root alone at every weight (smaller pieces, smaller
  sums).
- **The fix** (`fading_forms3.py --split`): split the forms by their Hecke numbers on lane ℓ first, then read each
  piece's polygon. A piece whose last coefficient is 0 to the precision is marked "NOT determined".
- **23, weight 504, split** (counted; console run): the 22 forms fall into 13 systems, and every polygon but one is
  determined. Exactly **two** slope-2 roots, as PREDICTED (2), with first digits 8 (system 19, 20, 13, 5) and 17
  (system 15, 2, 22, 4): again the first digits of the two notes above.
- The slow primes 29 and 31 are therefore not stuck; the same split should isolate their slope-1 notes.

**PREDICTED before the next runs (the slow primes, now that the split removes the precision artifact).** At 29 and 31
the fading part is one note of size ℓ (`tower_fading`), so:
- **Split by systems** (`fading_forms3.py --any --split --slope=1`, weights 810, 1,622, 2,434 at 29 and 928, 1,858,
  2,788 at 31): among the forms new at 2 and 3 with the class's signs, exactly one root of slope 1 at each weight. Its
  first digit (3|ℓ)·u/ℓ is 23 at 29 and 2 at 31.
- **Carried to weight −2**, its digits agree with the fading notes measured from partition numbers: 29 → 23, 3, 10;
  31 → 2, 28, 13.
- If there are two slope-1 roots, the tower would have shown two notes (it showed one), unless one is silent in the
  tower.

**A fault of mine, fixed before any result was read.**
- The first `--split` runs skipped every piece whose polygon had its left end at the precision cap. From the second
  weight on, that included the pieces holding the slope-2 roots.
- Those runs (19, 23, 29, 31) were stopped and their partial outputs deleted.
- The rule now asks only that the slope-S segment's own two ends have valuations below the precision. Then the capped
  points lie to its left, on steeper segments, and cannot change it.

**Counted at 19 (`fading_forms3_19_s2_split.txt`; 5 weights, 340 to 1,708; 21 minutes): as PREDICTED (2 and 3).**
- **At every weight there are exactly two slope-2 roots with the class's signs**, one in each of two systems on lane 19
  (a₅, a₇, a₁₁, a₁₃ ≡ 16, 16, 18, 10 and 0, 7, 16, 8). Once the forms are split, each root is alone in its piece. The
  other forms in those pieces are much steeper (slope 12, or about 35/3).
- **Carried to weight −2** (Lagrange through the weights known to 8 digits or more):

  | the system on lane 19 | (3\|19)·u(−2)/19² from the forms | the fading note from partition numbers (`fading_two`) | digits that agree |
  |---|---|---|---|
  | 16, 16, 18, 10 | 2, 1, 0, 17, 13, 15 | 2, 1, 0, 17 | 4 of 4 |
  | 0, 7, 16, 8 | 13, 14, 2, 15, 13, 15 | 13, 14, 2, 15 | 4 of 4 |

  - Through 2 weights, 2 digits agree; through 3, all 4 (the partition numbers carry 4).
  - The 5th and 6th digits stay the same through 4 and 5 weights. That both notes end in 13, 15 there is not examined;
    by chance it happens 1 time in 361.
  - A chance agreement of both notes in all 8 digits would happen about 1 time in 19⁸ ≈ 1.7·10¹⁰.
- **So 19's two quiet notes are the weight −2 members of two slope-2 families of level-6 forms with the class's signs**,
  in two different systems on lane 19, neither a θ-twist of the class. That is the finding of 13 and 17, twice over.

**Counted at 23 (`fading_forms3_23_s2_split.txt`; 4 weights, 504 to 2,022; 31 minutes): as PREDICTED (2 and 3).**
- At every weight there are exactly two slope-2 roots with the class's signs, in two systems on lane 23 (a₅, a₇, a₁₁,
  a₁₃ ≡ 19, 20, 13, 5 and 15, 2, 22, 4). Once the forms are split, each is alone in its piece.
- **Carried to weight −2:**

  | the system on lane 23 | (3\|23)·u(−2)/23² from the forms | the fading note from partition numbers | digits that agree |
  |---|---|---|---|
  | 19, 20, 13, 5 | 8, 20, 21, 8, … | 8, 20, 21, 8 | 4 of 4 |
  | 15, 2, 22, 4 | 17, 20, 19, 21, … | 17, 20, 19, 21 | 4 of 4 |

  - Through 2 weights 2 digits agree, through 3 and 4 all four.
  - Digits 5 and 6 changed between 3 and 4 weights (to 19, 9 and 22, 16), so they wait for a fifth weight.
- **So the finding of 13, 17 and 19 holds at 23 too:** each quiet note is the weight −2 member of a slope-2 family with
  the class's signs. At 19 and 23 there are two such families, in two systems, and both sound in the tower.

**Counted at 29 and 31, two weights of three (`fading_forms3_29_s1_split.txt`, `fading_forms3_31_s1_split.txt`).**
- **Both runs were stopped** at their third weight (2,434 and 2,788) by Claude Code, because the PC ran critically low
  on memory. That was not a fault of the runs. They were not restarted; that waits for the user's word.
- **As PREDICTED at both weights done:** among the forms with the class's signs, exactly one root of slope 1, alone in
  one system on lane ℓ.
  - 29: weights 810 and 1,622; system a₅, a₇, a₁₁, a₁₃ ≡ 25, 20, 12, 5. The other form in its piece has slope 23, then 12.
  - 31: weights 928 and 1,858; system 15, 16, 21, 23. The other has slope 25, then 12.
  - At 31 the unsplit polygon showed "38 roots of slope 1". That was the precision artifact again; split, only one is
    of slope 1.
- **Carried to weight −2 through the two weights** (2·w(0) − w(1), from the printed digits; derived):

  | ℓ | (3\|ℓ)·u(−2)/ℓ from the forms | the fading note from partition numbers (`tower_fading`) | digits that agree |
  |---|---|---|---|
  | 29 | 23, 3, (4) | 23, 3 (the 10 was not settled; see the correction below) | 2 of 2 |
  | 31 | 2, 28, (23) | 2, 28 (the 13 was not settled) | 2 of 2 |

  - Two weights carry two digits, as at 13 to 23. The third digit needs the stopped third weight.
- **So the slow primes' quiet notes are found too, to two digits.** Each is the weight −2 member of a slope-1 family with
  the class's signs, in another system on lane ℓ. A chance agreement of both second digits would happen 1 time in 899.

**PREDICTED before the third weight was run again**.
- **The run:** `fading_forms3.py --from=2`, weight 2,434 at 29 and 2,788 at 31, one prime at a time to spare memory. The
  first two weights are taken from the kept outputs.
- **At the third weight:** exactly one slope-1 root with the class's signs, in the same system as before (29: 25, 20, 12,
  5; 31: 15, 16, 21, 23). Its first digit is again 23 at 29 and 2 at 31.
- **Carried through all three weights** (Lagrange at x = −1: 3·w(0) − 3·w(1) + w(2); `fading_slow_carry.py`): the
  digits agree with the fading notes from partition numbers in all three (29: 23, 3, 10; 31: 2, 28, 13). [Corrected 13 Nov:
  the third digits 10 and 13 were not settled on partition numbers; see "29's quiet note by the shortcut".]

## 607 on a page (13 Nov 2026; `chords_page_data607`)

The page "The Exception 607" draws 607's tower in the language of the other two pages. Its data
(`chords_page_data607.txt` and `.json`; 607's class at weight 604 modulo 607³, 30 seconds; counted unless marked):
- **The notes on lane 607**, from the constants' polynomial (3|607)·T₆₀₇ = −T₆₀₇:
  - a lone note, constant 406 (202 points);
  - the pair, constant 240 twice (303 points); T₆₀₇ ≡ 367 on both, and 240 = −367;
  - a braid of 22: its constants generate the field with 607²² elements; its order is (607²² − 1)/2.
- **The braid of 22 in size and direction:**
  - size ×53 every two levels (the product of its 22 constants; a square; 303 places);
  - its direction has order half of its circle's (607²² − 1)/606 places.
  - The tower's loop is the braid's order, about 8.5·10⁶⁰ two-level steps.
  - 607²² − 1 = 2⁶ · 3 · 19 · 23 · 89 · 101 · 7,459 · 153,913 · 444,841 · 5,167,686,547 · 1,976,083,289,819 ·
    4,319,117,608,586,023. Its logarithm was not computed (the largest prime factor is about 4.3·10¹⁵).
- **The pair in size and direction** (field with 607² elements; directions on the circle of 608 places, counted from
  form A's place at 5; each form followed by its own eigenvector):

  | where | form A (size, direction) | form B | |
  |---|---|---|---|
  | at 5 (roots of x² + 561x + 66) | 66, 1 | 66, 607 | mirrored |
  | at 607, first digit (both 367) | 542, 0 | 542, 0 | they meet |
  | one digit deeper (roots of x² + 70x + 484) | 484, 342 | 484, 266 | mirrored again |

- The deeper digit is that of the forms at weight 604. The tower's own constants differ in that digit; it would need
  weights near 368,000 and was not computed.
- **Directions recounted from a fixed reference** (`chords_page_data607b.txt` and `.json`, derived from the kept
  numbers).
  - The first count used form A's own place at 5 as its reference, so it drew the two forms at 5 side by side (+1 and −1).
    That came from the choice, not from the forms.
  - Counted from a fixed generator, 1 + 2·√3 of F₆₀₇(√3), the directions are: at 5, 409 and 199; at 607, 0 and 0; one
    digit deeper, 38 and 570.
  - Any reference keeps each pair mirrored and keeps the meeting at 0. The page uses these (version 2).

### 607's pair, one digit deeper in the tower itself (13 Nov 2026; `partitions_flint_607_deeper`; prompted by a review of the page)

**The point raised.** The page drew "one digit deeper they part again" from the forms at weight 604, while saying the
tower's own constants differ in that digit. So the parting was shown for the forms, not for the tower.

**The method (derived; README "a twelfth time").** The pair's U₆₀₇-numbers follow a Hida family. Its two members are
conjugate, and its local Hecke ring has rank 2 with residue field F_(607²), so it is W(F_(607²))[[T]]. Along weights
604 + 606·j, the number is a power series in j whose n-th differences are divisible by 607ⁿ. Hence, writing
u = 367 + 607·δ:
- u(−2) ≡ 2·u(604) − u(1210) (mod 607²), so δ(−2) = 2·δ(604) − δ(1210) in F_(607²);
- the tower's constants (3|607)·u(−2) part on the second digit exactly when δ(−2) is not in F₆₀₇.

**PREDICTED before the run (derived).**
1. On the pair at weights 1210 and 1816, T₅'s polynomial mod 607 is again x² + 561x + 66 (the residual system is the same
   along the family), and T₆₀₇ ≡ 367 on exactly a 2-dimensional piece of the class.
2. The second difference vanishes: δ(1816) − 2·δ(1210) + δ(604) ≡ 0 (mod 607), for each form.
3. Asked, not predicted: is δ(−2) in F₆₀₇? If not, the tower's own constants part on the second digit.

**Counted (`partitions_flint_607_deeper.txt` and `.json`; about 10 minutes).**
- **Prediction 1 holds.** At 604, 1210 and 1816, T₆₀₇ ≡ 367 on exactly a 2-dimensional piece of the class. T₅'s
  polynomial on it is x² + 561x + 66 each time.
  - The class has 25 forms at 604 and 50 at 1210, of which exactly 25 are ordinary (Hida's control).
  - At 604, C's polynomial is x² + 70x + 484 again: a third computation, at precision 607³.
- **Prediction 2 holds.** The second difference δ(1816) − 2·δ(1210) + δ(604) is 0 (0 predicted).
- **The answer to 3:** δ(−2) = 480 + 357·t for form A, not in F₆₀₇.
  - So the tower's own constants, (3|607)·u(−2) = 240 + 607·(−1 − δ(−2)), PART on the second digit: one digit deeper the tower hears two voices.
  - In size and direction (fixed model of `chords_page_data607b`), their second digits are: form A [13, 376], form B [13, 232].

**The review's other points (13 Nov 2026), and what changed on the page (version 3):**
- **"One digit deeper they part again" was shown for the forms at weight 604, not for the tower.**
  - The tower's own second digit is now computed above, and the page shows both.
  - Version 2 claimed it for the tower without showing it.
- **The chance estimate was one mechanism with one pair per prime.** It is now labelled a heuristic: a random polynomial
  mod ℓ has a repeated root with chance 1/ℓ, from any two notes landing together. Over the primes from 37 to 1,297
  that is 0.67 (derived).
- **Boylan and Webb's sentence is now checked in the paper itself** (the PDF's text, p. 2203): "For all primes ℓ ≤ 1297
  with the exception of ℓ = 607, our calculations show that r_ℓ(1) = R_ℓ, and hence, that r_ℓ(m) = R_ℓ for all m ≥ 1.
  For ℓ = 607, we find that r_607(m) = R_607 − 1 for all m ≥ 1."
  - Their §6.1 lists the bases they computed: m = 1 for the primes 13 ≤ ℓ < 1300, but m = 2 only up to ℓ = 89. So
    for 607 the second digit on was not computed there.
  - If "rank" counts generators, a pair whose constants part on the second digit would need a 25th generator mod 607²
    (derived). How their "for all m" fits that is left open here.
- **"Proven" for the first-digit rule (A3) and for E1 was too strong**: no argument was written down.
  - A3's argument is now written out, from the published statements our literature notes quote verbatim (PROOF.md, "The
    first digit, written out"). One point, the count of notes with multiplicity, still needs Yang 2014 itself.
  - E1's argument is not written out yet. Both are labelled "derived", with their checks (STATEMENT).
- **The sweep of 200 primes** is one run of our own code. The page now says so.

### The review of the page, assessed point by point (13 Nov 2026)

The user asked for each point to be tested, not accepted. Version 3 of the page had taken all five as corrections.
Tested:

| the point | legitimate? | what the data or the argument shows |
|---|---|---|
| "One digit deeper they part again" was shown with the forms at weight 604, not the tower | **yes, about the evidence** | The claim itself holds. The tower's own constants were computed: 2·u(604) − u(1210), checked by 1816 (second difference 0). They part on the second digit. x² + 70x + 484 was computed three times (607⁴ on 12 Nov, 607³ twice today), the same each time. |
| The 1-in-ℓ estimate counts only a mirrored pair, and one per prime | **half** | The wording described one mechanism; the number did not. 1/ℓ is the chance that a random polynomial of degree ≥ 2 mod ℓ has a repeated root, whatever the cause. Σ 1/ℓ over the 200 primes is 0.67 (derived). It is a heuristic and is labelled one. With 0.67 expected cases, no sweep of this size can test it sharply. Seeing exactly one has chance 0.34 under it. |
| Not confirmed that Boylan and Webb computed the rank at every prime up to 1,297 | **no** | Their paper says so (its text, extracted from the PDF): §6.1 lists bases for m = 1 at every prime 13 ≤ ℓ < 1300, and p. 2203 states r_ℓ(1) = R_ℓ for all ℓ ≤ 1297 except 607. New from that reading: for m = 2 they computed only up to 89. |
| "(proven)" for the first-digit rule, with no proof seen | **yes, about the record; no, about the claim** | No argument was written down in our files. It is now (PROOF.md, "The first digit, written out"): a few lines from AAT (2.4), (2.6) and (2.8) with AAD (1.4). A3 is labelled "proven" again *(since 20 Nov: "derived", as is E1; see STATEMENT.md)*. E1 had the same gap; it is now written out too (PROOF.md), its sign identity checked on 5,475 prime pairs. |
| The sweep of 200 primes was rechecked only at 13 of them | **not a flaw** | That is the reviewer's own coverage. Our sweep (the forms) and Boylan–Webb's computation (the tower's own levels) are two different methods. They agree at all 200 primes, so each checks the other. |

- **The steps of mine that over-accepted:**
  - lowering A3 and E1 to "derived" for want of a written argument (now restored, with the arguments written);
  - "one run of our own code" for the sweep, which ignored the independent agreement with Boylan–Webb;
  - logging all five points under "Corrected" in FINDINGS.
- **What stays corrected:** version 1's drawing of the pair (the reference); the estimate's wording; and the evidence
  for the second digit, which now comes from the tower itself.

- **Version 5 of the page (13 Nov 2026).**
  - The page now carries its own support, in the section "The numbers behind the pictures": 607's class, the pair
    weight by weight (604, 1210, 1816), the first-digit argument with its published inputs quoted, and all 200 primes of
    the sweep.
  - It is built from `chords_page_data607c.json` (`chords_page_data607c.py`, read from the kept outputs).
  - No local file is named on the page any more. The pages "Towers as Chords" (version 10) and "Notes on a Donut"
    (version 5) were cleared of file names too. The chords page now checks its braid-angle rules from its own data.

## Is one note a power of another, in a pattern across primes? (13 Nov 2026; `note_powers_check`)

**The question** (raised by another session about the donut views: "A pattern in those powers across primes would be worth
a look"). Where a tower has two or more ordinary notes, the "Together" view writes one constant as a power of the other.
- 37: 34 = 15²⁰.
- 41: 23 = 31⁷.
- 43: 6 = 31⁷.
- 53: 7 = 43³.

Is there a pattern in those exponents?

**PREDICTED before the run (derived).**
- **No pattern.** The exponent is a discrete logarithm. Two notes' constants are (3|ℓ)·a_ℓ of different forms, and
  nothing ties their logarithms together.
- **The algebraic tie is elsewhere.** Twin notes are the two roots of one quadratic, so their sum and product (T and N,
  the two-term law) carry the relation, not their logarithms.
- **Measured** on every prime from 37 to 500 whose class has two or more ordinary notes on lane ℓ:
  - the share of pairs where one constant is a power of the other matches what random pairs in a cyclic group would
    give;
  - the exponents, as fractions of the first constant's order, spread evenly over (0, 1).

**Counted (`note_powers_check_500.txt`; 3.5 minutes): as PREDICTED, no pattern.**
- **The pairs:** 24 primes from 37 to 500 have two or more distinct ordinary notes on lane ℓ, giving 70 ordered pairs.
- **Which are powers:** one constant is a power of the other in 43 pairs. Random pairs in the same groups would give
  43.4 (derived).
- **The exponents spread evenly**, as fractions of the first constant's order: by tenths 3, 7, 4, 8, 4, 3, 6, 2, 4, 2
  (χ² = 8.9 on 9 degrees of freedom).
- **The four known cases come out again:** 34 = 15²⁰ (37), 23 = 31⁷ (41), 6 = 31⁷ (43), 7 = 43³ (53).
- **So the exponents behave like discrete logarithms of unrelated numbers.** The tie between two notes is algebraic, in
  their sum and product, not in their powers.

## Notes in a line: patterns across rows, and one line through the levels (13 Nov 2026; `notes_line_data`)

**The user's questions.**
- "What happens if we plot all of the notes in a line, from start to back to start once, and line them all up starting
  at 1? ... are there patterns that are not right next to each other, like up one note and over three?"
- "Could we calculate single lines through all the different levels given we have the starting position?"

**The data.**
- **Every ordinary note on lane ℓ:** 13 to 97 and 607, 21 notes. Each walks 1, α, α², … back to 1, and becomes a row,
  sorted from the smallest prime to the largest.
- **The braids' sizes,** as extra rows.
- **The towers' own values on lane ℓ,** H_b(n) mod ℓ, from exact partition numbers: levels 1 to 9, n = 0 … 15, at the
  one-note primes 13 to 31; and 37's kept values at levels 1 to 9 (`chords_page_data.json`).

**PREDICTED before the run (derived).**
1. **Shifts between rows ("up one, over three") carry nothing beyond chance.**
   - Two rows of different primes walk mod different numbers, and nothing ties them.
   - Test: for every rows-apart Δr = 1 … 6 and steps-over Δc = −12 … 12, the circular correlation of the two rows'
     places on the clock, against the same statistic with each row's steps shuffled (2,000 shuffles).
   - Predicted: no offset beats the shuffles beyond what the number of offsets tested allows.
2. **The exact patterns are stretches within one prime.**
   - Where one note is a power of the other (b = a^e), row b at step j equals row a at step e·j.
   - So the twins match perfectly under a stretch: 37 ×20; 41 ×3 and ×7; 43 ×7; 53 ×3 and ×9.
3. **One note: one start gives the whole line.**
   - At 13, 17, 19, 23, 29 and 31, put each tower value at its logarithm.
   - Then for every n, H_(1+2j)(n) sits at log H₁(n) + j·log α (mod ℓ − 1), a straight line (a spiral on the rings).
   - The values that are 0 on lane ℓ stay 0 at every level.
4. **Two notes: one start is not enough, two are.**
   - At 37, no single line fits all n.
   - From levels 1 and 3, the two-step rule H_(b+4) ≡ 12·H_(b+2) − 29·H_b (mod 37) gives levels 5, 7 and 9 for every n.

**Counted (`notes_line_data.txt` and `.json`, 15 minutes, up to p(415,322,398,107,207); `notes_line_37.txt` and
`.json`): as PREDICTED in all four.**
1. **Shifts carry nothing beyond chance.**
   - Of the 150 shifts tested, the largest |z| is 3.17 (5 rows apart, +5 steps over).
   - In 19% of the 2,000 shuffles the largest |z| was higher than that (the 95th percentile is 3.59).
2. **The real patterns are stretches, all inside one prime.** There are exactly six, and all are power relations:
   - 37: 34 = 15²⁰;
   - 41: 31 = 23³ and 23 = 31⁷;
   - 43: 6 = 31⁷;
   - 53: 43 = 7⁹ and 7 = 43³.
3. **One note: one start gives the whole line.** At 13, 17, 19, 23, 29 and 31, every nonzero tower value at levels 1 to 9
   (n = 0 … 15) sits on the line from its own start, and none is off it.
   - The counts on the line: 70, 80, 80, 80, 80 and 75.
   - The zeros (10 at 13, 5 at 31) stay 0 at every level.
4. **Two notes: two starts.**
   - At 37 one multiplier fits all five levels for 1 of 19 values of n, the one where note 34's part is 0.
   - The two-step rule holds 41 of 41 times.
   - Each note's part follows its own constant: 60 of 60 (note 15) and 50 of 50 (note 34).
- **A misreading of mine, caught before the page.** `notes_line_data.py` read 37's kept values as numbers. They are
  logarithms to the root 2, as the chords page stores them, so its 37 line ("1 of 38") is wrong. `notes_line_37.py` reads
  them as logarithms and gives the counts in 4. The page uses `notes_line_37.json`.
- **The page "Notes in a Line"** shows all four, with its data inside it.

## The quiet notes without the loud note's constant (13 Nov 2026; `quiet_hankel_check`)

**The question**.
- On the first digit the quiet notes cannot be seen: they are multiples of ℓ.
- Deeper, they are the bend of the tower's path away from the loud note's straight line.
- That is the same information as E_b = H_(b+2) − α·H_b, which needs α to every digit used (30, from the forms).

**The trick (derived; Prony's method).**
- If H_(1+2j) = c·α^j + d·β^j + (smaller), then D_j = H_(j+2)·H_j − H_(j+1)² = c·d·(α − β)²·(αβ)^j + (smaller).
- So the loud note cancels without knowing α, and D shrinks by α·β every two levels.
- Larger determinants of the same kind peel off note after note.

**PREDICTED before the run.**
- **Test:** at 13 and 17, from exact partition numbers alone (levels 1 to 11, n = 0 … 7).
- **The size:** the ratio D_(j+1)(n)/D_j(n) has valuation 2.
- **Its unit part is α·β/ℓ²,** the same for every n, to more digits the deeper the levels. Known values to compare:
  - α to 30 digits from the forms;
  - β/ℓ² from `tower_fading` (13: 10, 11, 6, 2, 10, 10; 17: 3, 9, 1, 13, 4, 3).

**Counted (`quiet_hankel_check.txt`; partition numbers up to p(249,899,243,909,824), 6 minutes): as PREDICTED.**
- **The size:** the ratio D_(j+1)/D_j has valuation 2 for every n, at both primes.
- **The unit part, without α:** formed from partition numbers alone, it agrees with α·β/ℓ² from the known constants in
  every digit all 8 values of n share:

  | ℓ | digits shared by all n, levels 1–5, 3–7, 5–9 | they read | agree with α·β/ℓ² |
  |---|---|---|---|
  | 13 | 0, 3, 4 | 8, 9, 8, 2 | all |
  | 17 | 1, 3, 5 | 1, 6, 2, 7, 0 | all |

- **So the quiet note can be read without the loud note's constant.** The precision is the same as E_b gives at the same
  depth (13: 3 and 4 digits at levels 3–7 and 5–9 there too).
- **The gain** is that nothing from the forms is needed. Larger determinants of the same kind take a chord apart
  layer by layer, which opens the quiet notes of the towers from 37 up.

## The prism at 5: 37's tower as white light (13 Nov 2026; the page "The Prism at 5")

**The user's picture.** A prism splitting white light into colours: the tower arrives mixed, the step at 5 separates it
into notes, and each note then follows its own simple rule.

**What the page shows,** from the kept 37 data: the tower c(M), M = 37^b·(24n + 11), and the step at 5's two outputs
G = c(25M) + (κ − λ)·c(M), all mod 37 (`chords_page_data.py`, read through `notes_line_37.json`).
- **Derived and checked in the page's build** (counted; nothing new computed from partition numbers):
  - **The shares:** each note's share, c_A = G_A/(λ_A − λ_B) and c_B = G_B/(λ_B − λ_A), adds back up to the tower at 79 of
    79 points.
  - **The value heard at 5,** c(25M), rebuilt from one combination, gives the other exactly, 79 of 79.
  - **Each share follows its own rule:** note 15's share ×15 every two levels, 60 of 60; note 34's ×34, 50 of 50.
  - **The tower does not:** between levels 1 and 3 it has 15 different multipliers among 19 values of n (3 → 5: 15 of 19;
    5 → 7: 14 of 17).
  - **One colour only at n = 0:** there the tower is note 15's share alone (1, 15, 3, 8, 9) and note 34's share is 0
    throughout.
- **Why it works (derived):** the step at 5 and the step down the tower commute. So a share that is pure for one is
  pure for the other, the same reason a single colour stays one colour after the prism.

### 607 through the prism (13 Nov 2026; `prism607_data.txt` and `.json`)

**Why not from partition numbers.** At 37 one combination of c(M) and c(25M) separates two notes. To isolate one of 607's
25 forms, the step at 5 would need values up to 25²⁴ times further along. And 607's tower is within reach of exact
partition numbers only to level 5. So the 607 picture is drawn from the forms' exact numbers.

**Counted (weight 604, mod 607, 29 seconds).**
- **The step at 5 tells all 25 forms apart.** T₅ has 25 distinct roots on the class, and the three pieces' T₅
  polynomials are pairwise coprime:
  - the lone note: a₅ ≡ 142;
  - the pair: the conjugate roots of x² + 561x + 66;
  - the braid of 22: the roots of one irreducible polynomial of degree 22.
- **The tower's own step does not.** T₆₀₇ ≡ 367·I on the pair on lane 607, so down the tower, on the first digit, both of
  the pair's forms are multiplied by 240.
- **So the prism at 5 sees 25 colours and the tower hears 24** (T₆₀₇'s minimal polynomial has degree 24). One digit
  deeper the 25th returns: the tower's own constants for the pair differ mod 607² (`partitions_flint_607_deeper`).
- This is Boylan–Webb's exception, said in one line: two forms that differ at 5 but agree at 607 on the first digit.
- On the page "The Prism at 5", section "A colour the tower cannot see" (version 3). Version 2 added the complete
  picture for 37: every n at once, and the grids of every n and level.

## 29's quiet note by the shortcut, beside the forms run (13 Nov 2026; `quiet_hankel_29`)

**A correction first (counted, from `tower_fading.txt`).** At 29 and 31 the fading notes measured from partition numbers
are settled to two digits only:
- 29 → 23, 3; 31 → 2, 28 (shared by all 8 values of n at levels 5 → 7).
- The third digits quoted before, 10 and 13, came from the first value of n alone. `tower_fading.py` printed one digit
  more than all n shared.
- So in "Counted at 29 and 31" the forms agree with every settled digit (2 of 2, not 2 of 3).
- The prediction "in all three (23, 3, 10; 2, 28, 13)" for the third weight has no settled third digit to meet. A
  third digit needs deeper tower values too.

**PREDICTED before the run.**
- **The run:** at 29, H_b(n) for levels 1 to 9 (n = 0 … 7) and level 11 (n = 0 … 2, partition numbers up to 3.4·10¹⁶);
  D_b = H_(b+4)·H_b − H_(b+2)².
- **The shortcut:** D_(b+2)/D_b has valuation 1 (α a unit, β of slope 1). Its unit part, divided by α, gives β/29 with
  digits 23, 3, ….
- **A third digit:** the deepest ratio (levels 5 to 11, n = 0 … 2) settles it, shared by the three values of n.
- **The cross-check:** the leftover E_(b+2)/E_b at levels 7 → 9 gives the same third digit.
- **The test of the forms:** the forms' third weight (still running) then agrees with it.

**Counted (`quiet_hankel_29.txt`; partition numbers up to p(34,059,756,429,262,106), 30 minutes).**
- **The shortcut, partly as predicted.** D_(b+2)/D_b has valuation 1 at every depth. Divided by α, its digits shared by
  all n are:
  - 23 at levels 3 … 9;
  - 23, 3 at levels 5 … 11 (three values of n).

  It settles two digits, not three: the shortcut uses one more level for the same digits, because it needs H_(b+4).
- **The cross-check settles the third digit.** The leftover E_(b+2)/E_b at levels 7 → 9 gives 23, 3, 4, shared by all
  three values of n.
- **So the tower's own quiet note at 29 is 23, 3, 4, …** The third digit agrees with the forms. The small weights gave
  23, 3, 4, 28, 20 (found while this run was going); the far weights gave 23, 3, (4) from two weights.
- Three routes, partly independent, agree on three digits.

## Toward 607's quiet notes: do nearby weights carry the quiet note's numbers? (13 Nov 2026; `fading_forms3 --small`)

**The problem.** 607's quiet notes are out of reach of both methods used at 13 to 31.
- **Partition numbers:** the tower reaches only level 5, under 25 loud notes known to one or two digits.
- **The forms:** a slope-1 family is guaranteed smooth only across weights ℓ(ℓ − 1) apart, and from −2 the first such
  weight is 367,840.

**The test (derived idea; to be checked where the answer is known).** Take weights only ℓ − 1 apart, −2 + (ℓ − 1)·x, a
range where the family is not guaranteed to hold together. Do the quiet note's system on lane ℓ and its first digit
appear there all the same? If so, 607's quiet notes come within reach at weights 1,210 and 1,816. If not, the route is
U₆₀₇ computed directly on overconvergent forms of weight −2.

**PREDICTED before the run (a guess to test, not derived).**
- **29:** split by systems, at the weights 54, 82, … 222 (x = 2 … 8), some weights show a slope-1 root in the quiet
  note's system (a₅, a₇, a₁₁, a₁₃ ≡ 25, 20, 12, 5), with first digit (3|29)·u/29 ≡ 23.
- **31:** at 58, 88, … 238, the same for the system 15, 16, 21, 23 with first digit 2.
- If neither shows, the nearby weights do not carry the quiet note, and the shortcut is closed.

**Counted (`fading_forms3_29_s1_split_small.txt`, `_31_s1_`, `_13_s2_`, `_17_s2_split_small.txt`; seconds each).**
- **The system: as PREDICTED.** At every nearby weight the quiet note's own system holds exactly one root of the
  quiet note's slope:
  - 29: 25, 20, 12, 5 at weights 54 … 222;
  - 31: 15, 16, 21, 23 at 58 … 238;
  - 13: 11, 8, 0 at 34 … 94;
  - 17: 4, 2, 9, 7 at 46 … 142.
- **The first digit: NOT as predicted.** It does not stay at the tower's value. It moves smoothly with the weight; at
  29 it is 6, 12, 18, 24, 1, 7, 13, six more each time.
- **Carried back to weight −2, it reproduces the quiet notes,** about one digit per weight:

  | ℓ | slope | (3\|ℓ)·u(−2)/ℓ^s from the nearby weights | known before | agree |
  |---|---|---|---|---|
  | 29 | 1 | 23, 3, 4, 28, 20 (7 weights) | 23, 3 settled on partition numbers | 2 of 2 |
  | 31 | 1 | 2, 28, 23, 8, 3 (7 weights) | 2, 28 settled; the far-away weights also gave 23 third | 2 of 2 |
  | 13 | 2 | 10, 11, 6, 2 (6 weights) | 10, 11, 6, 2, 10, 10 | 4 of 4 |
  | 17 | 2 | 3, 9, 1, 13 (6 weights) | 3, 9, 1, 13, 4, 3 | 4 of 4 |

- **Why it can work (derived):** the characteristic series of U_ℓ is analytic in the weight everywhere. While the
  quiet note's root stays alone on its slope, it follows the weight analytically too. The "guaranteed" disc is a
  worst-case bound.
- **So the quiet notes can be read from weights only ℓ − 1 apart.** At 29 this also gives a third digit, 4: a
  prediction for the shortcut run on partition numbers, which was running when this was found.

**PREDICTED before the 607 run (a guess to test).**
- **The run:** `fading_forms3.py 607 3 4 --slope=1 --any --split --small --from=1 --extra=8 --minp=4`, weights 1,210,
  1,816 and 2,422, modulo 607^12.
- **The count:** among the forms with the class's signs, exactly 25 roots of slope 1 at each weight, as many as the class
  has forms. This follows the pattern of 29 and 31, where weight 2ℓ − 4 has one slope-1 form for the one class form.
- **Rational systems:** each sits in a piece of its own, and its carried first digits are 607's quiet notes. These
  cannot be checked on partition numbers.
- **Conjugate groups:** some slope-1 roots may come in conjugate groups like the pair and the braid. Those need the
  bigger fields, and are only counted here.

## 607's quiet braids (13 Nov 2026; `fading_forms3 --small` at 1,210, then `fading_forms4`)

**Counted at weight 1,210** (`fading_forms3_607_s1_split_small.txt`, modulo 607^12). The 50 forms with the class's
signs split by their systems on lane 607:
- **slope 0:** the lone note's family (1 form), the pair's (2) and the braid's (22);
- **slope 1, all 8 roots:** a conjugate group of 8 forms;
- **undecided:** a conjugate group of 17, whose polygon is "12/17 × 17". Its left end is not determined, because the slopes
  add up to more than the 12 digits kept; all 17 of slope 1 would add up to 17.
- So 607's quiet layer has no single note: its slope-1 forms come in groups, like braids.

**Counted at all three weights** (`fading_forms3_607_s1_split_small.txt`, modulo 607^12, 109 minutes): the forms with the
class's signs, by groups (one system on lane 607 each):

| weight | forms | slope 0 | slope 1, settled | not settled at 12 digits | settled, steeper |
|---|---|---|---|---|---|
| 1,210 | 50 | 25 (1 + 2 + 22) | the 8-group | the 17-group | — |
| 1,816 | 76 | 25 | the 8-group | the 17-group; groups of 10 and 13 | one form and a pair, slope 2 |
| 2,422 | 101 | 25 | the 8-group | two 17-groups; 10; 13; a 6-group at slope 2 | one form slope 2; pairs at slopes 2 and 3 |

- **The 8-group is of slope 1 at every weight**, a quiet braid of 8 that holds through the whole range.
- **The 17-group is there at every weight**, but its slopes add up to more than the 12 digits kept.
  - So "25 slope-1 roots at each weight" is not yet settled: it holds for 8, and 17 more if that group is of slope 1.
  - The run at 607^24 settles it.

**The tool for groups** (`fading_forms4.py`). It follows the whole group: the factor of T_ℓ's polynomial whose roots are
the group's u/ℓ^S, carried from weight to weight and taken to −2.
- **29** (a group of one): 23, 3, 4, 28, as `fading_forms3`.
- **47's loud braid** (a conjugate pair, slope 0): the constants' polynomial agrees with the record
  (`partitions_flint_family_steps_47.txt`, to 47⁷) in all 4 digits printed of both coefficients: 19, 29, 19, 0 and 18,
  37, 27, 6.

**PREDICTED before the 607 run** (`fading_forms4.py 607 3 4 --slope=1`, weights 1,210, 1,816, 2,422, modulo 607^24).
- **At each weight:** the 17-group is all of slope 1. Together with the 8-group that makes 25 slope-1 roots, as many as
  the class has forms.
- **Carried to weight −2:** each group's polynomial stays irreducible on lane 607, a quiet braid of 8 and one of 17. Its
  coefficients are known to about 2 digits from 3 weights.
- These are 607's quiet notes as the forms give them; partition numbers cannot reach them.

**Counted so far (`fading_forms4_607_s1.txt`, modulo 607^24; weights 1,210 and 1,816 done, 2,422 running): as PREDICTED.**
- **At both weights the 25 slope-1 forms are the 8-group and the 17-group,** every root of slope 1, and the two groups'
  polynomials cannot be split on lane 607 (irreducible of degrees 8 and 17).
  - The 8-group's polynomial is known to 607^16, the 17-group's to 607^7.
- **At 1,816 the steeper groups settle too:** one form, a pair and a 10-group, all of slope 2. The 13-group's slopes still
  add up to more than the 24 digits kept.
- **So 607's quiet layer is a braid of 8 and a braid of 17,** as many quiet forms as loud ones, none alone. The values at
  weight −2 follow when 2,422 is done.

**PREDICTED before the 2,422 result** (14 Nov 2026, written while the run was at 2,422; derived): how the run's own
output tells whether its last step is right.
- **The pattern behind it.** At 29 and 31 the quiet note's first digit is a constant times x at all 7 nearby weights
  (counted from `fading_forms3_29_s1_split_small.txt` and `_31_`: 6·x at 29, 29·x at 31). So it is 0 at x = 0, the
  weight ℓ − 3, where the class has no quiet form. At 607 that weight is 604, which has only the 25 loud forms.
  - If 607's strands do the same, a coefficient of a braid's polynomial that multiplies j strands together is a
    constant times x^j on the first digit.
- **The tool's last step** is a curve through x = 1, 2, 3, read at x = −1.
  - It reads x^j right for j = 0, 1, 2.
  - From j = 3 on it gives c_j = 6 − 8·2^j + 3·3^j in place of (−1)^j: 23 for x³, 121 for x⁴, and so on. No c_j is 0
    mod 607 up to j = 17.
  - With only two weights, c_j would be 3 − 2·2^j.
- **So, predicted:**
  1. The polynomials printed for weight −2 most likely do not factor as 8 and 17 on lane 607. Both would by chance
     about 1 time in 136 (a random polynomial of degree d is irreducible about 1 time in d).
  2. Multiplying each printed first digit by (−1)^j/c_j undoes the last step. The result factors exactly as 8 and 17:
     it is then the 1,210 polynomials, up to the signs of their roots.
- **If instead** the printed polynomials factor as 8 and 17 and the corrected ones do not, the straight-line picture is
  wrong for braids, and the printed first digits stand.
- **Either way,** the deeper digits need each weight's polynomial, which this run does not keep.

**Counted** (`fading_forms4_607_s1.txt`, all three weights, 2 hours 41 minutes; the check in
`fading_forms4_607_correct.py` and `.txt`). **As PREDICTED, both parts.**
- **At 2,422** the two quiet groups are again all of slope 1, and neither can be split on lane 607: the 8-group (G known
  to 607^16) and the 17-group (to 607^7). So at all three weights the quiet layer is these two braids, 25 forms
  (counted).
  - A second group of 17 at 2,422 is steeper. Its slopes add up to the 24 digits kept, so it is not settled.
- **1. The run's own carry to −2 splits both braids.** On lane 607 its polynomials factor as 3 + 5 and 1 + 1 + 6 + 9
  (counted).
- **2. Undone on the first digit, both are whole again.** Multiplied by (−1)^j/c_j, the polynomials are irreducible, of
  degrees 8 and 17 (counted). Both would be so by chance about 1 time in 136.
- **So the first digits of 607's quiet braids at the tower's weight** (derived): the polynomials of (3|607)·u/607,
  coefficients lowest first, are
  - for the braid of 8: 458, 176, 331, 412, 360, 410, 588, 195, 1. Its strands multiply to 458 and add up to 412.
  - for the braid of 17: 315, 423, 171, 125, 236, 389, 297, 460, 584, 262, 424, 472, 300, 473, 153, 541, 121, 1. Its
    strands multiply to 292 and add up to 486.
- **Why it works (derived).** Each strand's first digit is a constant times x, so on the first digit the value at −2
  follows from any one weight.
  - At 607 the two sign changes, x = −1 and (3|607) = −1, cancel. These are therefore the polynomials of u/607 at
    weight 1,210.
  - A run that prints each weight's polynomial would check that coefficient by coefficient.
- **The lesson for groups:** the plain curve through three weights follows x and x², not x³ or higher.
  - So it gets the first digit right for groups of one or two strands: 29's one strand, and 47's loud pair, whose digit
    does not move at all.
  - It gets it wrong for 607's braids. Carrying a braid needs each coefficient divided by x^j first (the scaled fit).
- **Still open:** the deeper digits. They need each weight's polynomial, which this run did not keep, and the scaled fit
  (3 digits from 3 weights at 29 and 31).

**The braids split into their strands, on the first digit** (`fading_forms4_607_split.py` and `.txt`, 14 Nov 2026;
counted). In each braid's own field, F₆₀₇[t]/(G), the strands are t, t^607, t^(607²), ….
- **The braid of 8** (G = the corrected polynomial above): all 8 strands are roots and distinct, and the 607th power steps
  through them, returning after 8 steps. They multiply to 458 and add up to 412. Strand 2 = t^607 = 445 + 159t + 605t²
  + 585t³ + 381t⁴ + 90t⁵ + 31t⁶ + 467t⁷.
- **The braid of 17:** the same checks hold for all 17 strands; they multiply to 292 and add up to 486.
- **As the run printed them** (no correction), the polynomials split into braids of 3 and 5, and into two lone notes
  (571 and 484) and braids of 6 and 9. That would change the quiet layer's shape at the tower's weight, against all three
  computed weights.

**PREDICTED before the next run** (`fading_forms5.py`: the same per-weight computation, each weight's polynomial kept in
a .jsonl file, both carries printed with their checks; weights 1,210, 1,816 and 2,422 as three processes, joined with
`--fit`; the script was checked by a review the user relayed, whose three points were fixed first).
- **The direct test at 607:** each coefficient divided by x^j has the same first digit at all three weights, for 8 of 8
  and 17 of 17 coefficients.
- **The first digits** of the scaled carry are exactly the corrected ones above (458, 176, … and 315, 423, …).
- **Digits:** the carries from every two of the three weights agree on about 2 digits, and the third digit is the full
  carry's own (as at 29, where it was right). Digits 2 and 3 also need each strand to be 0 at weight 604 beyond the first
  digit; at 29 and 31 the scaled fit held through all 5 known digits, and at 607 weight 604 has only the 25 loud forms.
- **The split** of the answer matches the one above on the first digit.

**Counted, weight 1,210 first** (`fading_forms5_607_s1_x1-1.txt` and `.jsonl`; 27 minutes): its own polynomials, read
on the first digit, equal the corrected polynomials above in **9 of 9** (braid of 8) and **18 of 18** (braid of 17)
coefficients.
- This is the coefficient-by-coefficient check named above: by the straight-line picture, the corrected polynomials
  at −2 are 1,210's own (the two signs cancel).
- So the first-digit correction is confirmed by direct computation, and the split built on it holds on the first digit.
  The old run's printed polynomials were wrong there.
- Each match, on its own, has odds of 1 in 607 if it were chance.
- The full per-weight test (each coefficient divided by x^j, with the same first digit at all three weights) comes with
  1,816 and 2,422.

**Counted, weight 1,816** (`fading_forms5_607_s1_x2-2.txt` and `.jsonl`; 53 minutes). The groups are again of slope 1
and whole on lane 607.
- **The direct test, at two weights: as PREDICTED.** Each coefficient divided by x^j has the same first digit at 1,210
  and 1,816, for 8 of 8 coefficients (braid of 8) and 17 of 17 (braid of 17).
- **A first look at the second digit** (both weights joined in the console only, so that no file blocks the final join;
  two weights, so no cross-check yet):
  - the scaled carry's second digits are 397, 107, 14, 20, 569, 139, 494, 544 (braid of 8) and 126, 0, 299, 26, 155,
    174, 294, 274, 45, 334, 402, 436, 431, 592, 105, 48, 388 (braid of 17);
  - the plain carry from the same two weights breaks both braids apart (1 + 2 + 5 and 1 + 1 + 3 + 12).
- 2,422 then gives the third digit, and the cross-check on the second.

**Counted, weight 2,422 and the join** (`fading_forms5_607_s1_x3-3.txt` and `.jsonl`, 85 minutes; the join in
`fading_forms5_607_s1_fit.txt`). **As PREDICTED in all four parts.**
- **The direct test, at all three weights:** each coefficient divided by x^j has the same first digit at 1,210, 1,816
  and 2,422, for 8 of 8 and 17 of 17 coefficients.
- **The first digits** of the scaled carry are the corrected ones.
- **Digits:** the carries from every two of the three weights agree with the full carry on the first 2 digits of every
  coefficient (confirmed), and the third digit is the full carry's own (likely). The plain carry's two-weight carries
  agree on 0 digits.
- **The split** on the first digit is the one above (8: product 458, sum 412; 17: product 292, sum 486).
- **607's quiet braids at the tower's weight,** the polynomials of (3|607)·u/607, coefficients lowest first, each as
  digits lowest first (the first two confirmed, the third likely):
  - **braid of 8:** [458, 397, 109], [176, 107, 591], [331, 14, 41], [412, 20, 251], [360, 569, 555], [410, 139, 116],
    [588, 494, 265], [195, 544, 29], [1];
  - **braid of 17:** [315, 126, 29], [423, 0, 103], [171, 299, 601], [125, 26, 411], [236, 155, 321], [389, 174, 185],
    [297, 294, 467], [460, 274, 194], [584, 45, 202], [262, 334, 548], [424, 402, 382], [472, 436, 27], [300, 431, 391],
    [473, 592, 92], [153, 105, 260], [541, 48, 545], [121, 388, 494], [1].
- These cannot be checked on partition numbers (the tower reaches only level 5 at 607).

**Each strand refined to three digits** (`fading_forms5_607_strands.py` and `.txt`; derived from the join's polynomials,
the third digit likely like theirs). In each braid's ring (Z/607³)[t]/(G), strand 1 is t and each next strand is the
root lying over the 607th power of the last, found by Newton's method (Hensel: the strands are distinct on the first
digit).
- **Checked (counted, modulo 607³):** every strand is a root; all 8 (all 17) are distinct; each strand's first digits are
  the first-digit split's.
- **Size and sum, digits lowest first** (both match the polynomials exactly):
  - braid of 8: the strands multiply to 458, 397, 109 and add up to 412, 62, 577;
  - braid of 17: they multiply to 292, 480, 577 and add up to 486, 218, 112.
- Example: the 8-braid's strand 2 = t^607 refined: [445, 415, 427] + [159, 195, 546]·t + [605, 106, 231]·t² + … (each
  coefficient's digits lowest first). All strands are in the .txt.

## What sets the settling rate: the class's growth along the weights (14 Nov 2026; a reading, then a test)

**The question** (open: NEXT item 5; the survey's "settling rate"). Why do 13, 17, 19 and 23 settle one digit per level,
and 29 onward one digit per two levels? The quiet notes set the rate: slope 2 gives the fast rate, slope 1 the slow one.
What sets their slope?

**A reading (derived from counts already kept, so not yet a test).** Along the row of weights ℓ − 3, 2ℓ − 4, 3ℓ − 5, …,
the class grows:
- **13 and 17 (fast):** 1 form at ℓ − 3, still 1 at 2ℓ − 4, then 2 at 3ℓ − 5. The new form there has slope 2
  (`fading_forms3_13_s2_split_small.txt`, `_17_`: weights 10, 22, 34 and 14, 30, 46).
- **29 and 31 (slow):** 1, then 2 already at 2ℓ − 4. The new form has slope 1 (weights 26, 54 and 28, 58).
- **607 (slow):** 25 at 604, then 50 at 1,210. All 25 new forms have slope 1, and they are the quiet braids.
- **So:** the quiet notes' slope is the number of steps the class stays at its size before it first grows.
  - If the class grows at the first step (2ℓ − 4), the quiet notes have slope 1, and the tower is slow.
  - If it waits one step, they have slope 2, and the tower is fast.
  - The same staircase shows in the kept slopes: at 29 each step adds one form of the next slope (0; 1, 0; 2, 1, 0;
    3, 2, 1, 0).
- **Why the split falls between 23 and 29 (derived):** the class gains about (ℓ − 1)/24 forms per step. That is below 1 for
  ℓ ≤ 23, so the first step can add nothing, and above 1 from 29 on, so it always adds. At 607, 606/24 = 25.25, and the
  first step added exactly 25 forms.

**PREDICTED before the test** (weights ℓ − 3, 2ℓ − 4, 3ℓ − 5; the class's systems and slopes; console only):
- **19 and 23 (fast, two quiet notes each):** 1 form at ℓ − 3 and still 1 at 2ℓ − 4 (34 and 42). Then 3 at 3ℓ − 5 (52
  and 64): two new forms, both of slope 2.
- **37, 41, 43 (slow):** the class (2 forms) grows already at 2ℓ − 4 (70, 78, 82), and the new forms there have slope 1.
- If any of these fails, the reading is wrong as stated.

**Counted** (`fading_forms5.py ℓ 2 6 --slope=9 --from=0`, console only, seconds each). **As PREDICTED at all five.**

| ℓ | rate | forms at ℓ − 3 | at 2ℓ − 4 | at 3ℓ − 5 | slopes of the new forms |
|---|---|---|---|---|---|
| 19 | fast | 1 | 1 | 3 | 2, 2 (weight 52) |
| 23 | fast | 1 | 1 | 3 | 2, 2 (weight 64) |
| 37 | slow | 2 | 3 | 5 | 1 (weight 70); then 2, 2 |
| 41 | slow | 2 | 3 | 5 | 1 (weight 78); then 2, 2 |
| 43 | slow | 2 | 3 | 6 | 1 (weight 82); then 2, 2, 2 |

- **The count of new forms is the count of quiet notes** measured from partition numbers: one at 13, 17, 29 and 31, two at 19
  and 23.
- **The rule, now tested** (a reading, counted at 13, 17, 19, 23, 29, 31, 37, 41, 43 and 607):
  - A tower settles one digit per level exactly when its class gains no form at the first step of weight, 2ℓ − 4. Its
    quiet notes then come one step later, with slope 2.
  - When the class gains forms at once, those forms are the quiet notes, with slope 1. The tower then settles one digit per
    two levels, Boylan and Webb's bound, met exactly.
- **Why 13 to 23 alone are fast (derived):** the class gains about (ℓ − 1)/24 forms per step. That is less than one for
  ℓ ≤ 23, so the first step can add nothing, and more than one from 29 on. At 607 it is 606/24 = 25.25, and the first
  step added 25 forms: its two quiet braids.
- **What it rests on:** counts and slopes of classical forms (counted) at the first three weights of the row.
  - That the slope equals the number of steps waited is, as far as I recall, the shape of Bergdall and Pollack's "ghost"
    description of slopes. That link is not checked against their papers here.
  - 5, 7 and 11, with empty classes, are outside this rule.

**A review of this (relayed by the user, 14 Nov 2026), assessed point by point.**
- **"The counts match an independent recount with PARI"** (theirs, not ours: newforms on Γ₀(6) with the class's signs).
  - It agrees with ours at every prime we both counted (13 to 43).
  - It adds 47 (2, 3, 6), 53 (2, 4, 6), 59 (2, 4, 7) and 61 (3, 5, 8).
  - Legitimate, and checked below with our own code.
- **"(ℓ − 1)/24 per step is true only on average: the gains alternate, because the class's signs flip each step."**
  - Legitimate. Our own code builds the signs that way (`sign(p, k)`: (p|ℓ) raised to a power that grows by one each
    step), and our counts alternate too (43: +1 then +3; 19: +0 then +2).
  - So the rule is stated by the first step's count itself, which is what was tested. The average only explains why the
    line falls near 24.
  - A proof that the first step always gains from 29 on needs the exact dimension formula for the class, not the average.
- **"Bergdall and Pollack's ghost conjecture, proved in a local form by Liu, Truong, Xiao and Zhao under a genericity
  condition; your forms are ordinary, the locally reducible setting."**
  - The link is the right one to test (quoted by the review; the papers not read here).
  - One point is not right as stated. The forms whose slopes the rule is about, the new forms at 2ℓ − 4, are not
    ordinary; they lie in other systems on lane ℓ than the class's ordinary forms (at 29, 25, 20, 12, 5, not the class's
    5, 25, 4, 28).
  - Each such system has its own residual representation, so the genericity condition has to be checked system by
    system, not for the class as a whole.

**PREDICTED before our own count of 47, 53, 59 and 61** (from the review's counts and the rule):
- **47:** 2 forms at 44, 3 at 90 (one new, of slope 1).
- **53 and 59:** 2 at ℓ − 3, 4 at 2ℓ − 4 (two new, both of slope 1).
- **61:** 3 at 58, 5 at 118 (two new, both of slope 1).
- All four are slow on partition numbers (`tower_speed`: every note from 37 to 97 at one digit per two levels).

**Counted** (`fading_forms5.py ℓ 2 6 --slope=9 --from=0`, console only, seconds each). **As PREDICTED at all four, and
the same counts as the review's recount with PARI.**
- **47:** 2, 3, 6; the first step adds one form, of slope 1.
- **53:** 2, 4, 6; it adds two, both of slope 1, in two separate systems.
- **59:** 2, 4, 7; two, both of slope 1, separate.
- **61:** 3, 5, 8; two, both of slope 1, braided together (one system of degree 2), as at 607.
- The second step adds forms of slope 2 at all four (the staircase).
- **So the rule holds at 14 primes, from 13 to 61, and at 607.** It is stated by the first step's count:
  - a gain of 0 means fast (13, 17, 19, 23);
  - a gain of g ≥ 1 means slow, with g quiet forms of slope 1 (1 at 29 to 47; 2 at 53, 59, 61; 25 at 607).
- **Not yet shown:** that the gain is at least 1 for every ℓ ≥ 29. That needs the class's exact dimension formula. The
  average (ℓ − 1)/24 makes it plain for large ℓ, but the step-to-step alternation needs the exact form.

**The rate from the quiet notes, back to 5, 7 and 11**.
- **The idea, as set down on 13 Nov:** the rate is set by the next family, which the forms alone can predict. On 13 Nov
  it was made exact: the rate in digits per two levels is the quiet note's slope minus the loud note's ("What fades: the
  speed's own note").
- **12 Nov:** that fitted 5 (slopes 2 and 7, rate 5). At 7 and 11 the rates 3 and 1 were half the gap to the next slopes
  then listed (7 and 4). The half was never explained (NEXT item 5).
- **The counts above show a staircase:** a form that first appears at step s of the row ℓ − 3, 2ℓ − 4, … has slope s
  (29: 0; 1; 2; 3. 13: 0; —; 2. 47, 53, 59, 61: new forms of slope 1, then 2).
  - If that holds where the class is empty, the slopes of 5, 7 and 11 are read from the steps at which forms first
    appear.
  - The law's family is the first to appear, and the rate is the gap to the second.

**PREDICTED before the count** (all forms with the class's signs; weights −2 + (ℓ − 1)·(x + 1), x = 0 … 8; console only):
- **The staircase holds:** every form, in its system, first appears at step s with slope s.
- **The law's families come first,** at step 2 at 5, 1 at 7 and 2 at 11: the slopes 2, 1, 2 of Watson's and Atkin's laws.
- **The next forms come at step 7 at 5, 4 at 7 and 3 at 11.** Their slopes, 7, 4 and 3, minus the law's give the measured
  rates 5, 3 and 1.
- If the next forms come elsewhere, the half at 7 and 11 stays unexplained.

**Counted** (`fading_forms5.py ℓ 8 6 --slope=9 --from=0`, console only; the script got a guard for empty spaces first,
which changes nothing where forms exist). **Partly as predicted.**

| ℓ | weights (x = 0 … 8) | the forms with the class's signs, by weight |
|---|---|---|
| 5 | 2, 6, …, 34 | none up to 14; from 18 one form, slope 2 (at 30 and 34 slope 3) |
| 7 | 4, 10, …, 52 | none at 4, 10; from 16 one form, slope 1; from 40 a second, in another system, slope 4 (at 52 the first shows slope 2) |
| 11 | 8, 18, …, 88 | none at 8, 18; from 28 one form, slope 2; from 48 a second, in another system, slope 3; from 68 a third, in the first one's system, slope 4 |

- **The rates: as PREDICTED at 7 and 11.**
  - The first new form outside the law's system has slope 4 at 7 and 3 at 11.
  - Minus the law's slopes 1 and 2, that gives 3 and 1 digits per two levels: the measured rates.
  - The "half" of 12 Nov came from taking the next slope inside the law's own system. At 11 that is the slope-4 form,
    seen here from weight 68.
  - As at 13 to 61, the quiet note lives in another system than the loud one.
- **The steps: NOT as predicted.** At 5, 7 and 11 forms first appear later than their slope says: at steps 4, 2 and 2,
  with slopes 2, 1 and 2, and the second at steps 6 and 4 with slopes 4 and 3.
  - So "slope = the step it first appears at" is not a rule. It held only at the first steps from 13 to 61.
  - The fast/slow rule above rests on its own counts, not on this.
- **At 5**, the second form is not reached by weight 34. The earlier slope 7 (12 Nov, gap 5, the measured rate) stands.
- **So the rate is derived from the quiet notes at every prime tested** (5, 7, 11, 13 to 61, 607): the quiet note's slope
  minus the loud note's, with the quiet note found in its own system.
  - At 7 and 11 only the rate is matched, not yet the note's digits on partition numbers.
  - Slopes vary at a few weights (the law's form at 5 and 7 near 30 and 52). Those weights lie closer, 7-adically or
    5-adically, to the empty start of the row. A reading.

**The quiet notes of 7 and 11, digit by digit**.

**PREDICTED before the runs:**
- **From partition numbers** (`tower_fading_law.py`: the tower divided by Watson's 7^(⌊b/2⌋+1) or Atkin's 11^b; α from
  `partitions_flint_slopes_ℓ.txt`, 20 digits at 7 and 12 at 11):
  - the fading part E_b = H(ℓ²M) − α·H(M) shrinks by 7³ every two levels at 7 and by 11 at 11, at every level, for every
    n;
  - its ratio E_(b+2)/E_b is one number for all n, to more digits the deeper the levels (one quiet note dominating).
- **From the forms** (`fading_forms5.py --plain`):
  - the slope-4 form at 7 (its own system, from weight 40) and the slope-3 form at 11 (from weight 48) are carried through
    12 small weights to weight −2, with the plain curve (these forms are not 0 at the row's start, so no division by x);
  - (3|ℓ)·u(−2)/ℓ^s of each is the unit part of the ratio above, in every digit the partition numbers settle.
- If the digits disagree, the quiet note at 7 or 11 is another form, or a mix.

**Counted** (`tower_fading_law.py` and `.txt`, 88 seconds, up to p(87,745,226,699,158); `fading_forms5.py` at 7, slopes 4
and 5, and at 11, slope 3, kept as `.txt`/`.jsonl`; the comparison in `fading_quiet_fit_7_11.py` and `.txt`).
**As PREDICTED at both.**
- **On partition numbers:**
  - at 7 the fading part shrinks by exactly 7³ every two levels, at every level and every n, and is one note: all 8 values
    of n share 7 digits of it, 6, 2, 3, 1, 4, 1, 3;
  - at 11 it shrinks by exactly 11 every two levels from level 3 on; all n share 4, 5, 0.
- **On the forms:** the quiet form at 7 (system 5, 6, 0) and at 11 (system 0, 4, 1) carried to weight −2 by the plain curve.

  | ℓ | (3\|ℓ)·u(−2)/ℓ^s from the forms | confirmed by the sub-carries | from partition numbers | agree |
  |---|---|---|---|---|
  | 7 (s = 4) | 6, 2, 3, then 1 | 2 (3rd likely) | 6, 2, 3, 1, 4, 1, 3 | 4 of the forms' 4 |
  | 11 (s = 3) | 4, 5, 0, 3, 1 | 4 (5th likely) | 4, 5, 0 | 3 of 3 |

- **Two things the fit needed:**
  - Weights whose steps differ by a multiple of ℓ (x = 6 and 13 at 7) make the curve divide by ℓ. So at most ℓ weights in
    a row can be used, and at 7 that limits the forms to about 3 digits.
  - At some weights the quiet form is steeper: at 7, slope 5 at x = 9 … 12 (and absent at slope 4 at 16, 17); at 11, not
    of slope 3 at x = 12, 13, 14.
    - Those are its zeros on the first digit. They fall at x ≡ 2, 3, 4, 5 (mod 7) and x ≡ 1, 2, 3 (mod 11): as many as its
      slope, at the steps just before it first appears.
    - At 7 it was read there at slope 5 and put on the same scale.
- **So the quiet notes of 7 and 11 are found, digit by digit.** Each is the next layer down from the law's family, in a
  system of its own, as at 13 to 31. With them, the settling rate is derived at every prime tested from 5 up: the quiet
  note's slope minus the loud note's.

**Verified independently** (14 Nov 2026; another session, given only the definitions and not our numbers; FLINT 0.9.0 for
exact partition numbers, PARI 2.17.2 for the forms; its report relayed by the user). Assessed point by point:
- **The forms: the same.**
  - At every weight exactly one ℓ-adic system had the stated residues.
  - Slopes: 4, 4, 4, 5, 5, 5, 5 at 7 (x = 6 … 12) and 3 at every x = 4 … 11 at 11.
  - Carried values: (3|7)·P(−1) = 6, 2, 3, 1, 0, 4, … (2 digits stable leaving one point out); P(−1) = 4, 5, 0, 3, 1, 2, …
    at 11 (4 stable). These are our values and our stabilities.
- **The partition numbers: the same where both reach.**
  - At 7: valuation 3 at every level; α from H₁₅/H₁₃ to 19 digits, equal to ours; shared digits 6, 2, 3, 1, 4, 1, 3 at b = 5.
  - At 11 their α, read off the partition numbers alone, has 6 digits, so E_(b+2)/E_b settles only one digit (4).
- **A weakness of our run, which the check exposed:** our 11 measurement took α from the forms (12 digits). So its 4, 5, 0
  was not wholly independent of the forms side.
- **Their remedy, a method without α:**
  - Fit H_(b+4) = σ·H_(b+2) − π·H_b from every pair of n, and take β as the root of x² − σx + π divisible by ℓ.
  - At 7: valuation 3, all 28 pairs agree on 16 digits, 6, 2, 3, 1, 4, 1, 3, 1, 5, 6, 1, 2, 5, 3, 1, 5 (b = 11).
  - At 11: valuation 1; 4, 5 from all 8 n at b = 7, and 4, 5, 0, 3 from n = 0 … 2 at b = 9.
- **The comparison:**
  - At 11, the forms (4, 5, 0, 3, 1) match the α-free partition numbers in all 4 digits both carry.
  - At 7 they match in 4 digits; the fifth (forms 0, partition numbers 4) lies past the 2 to 3 digits the forms support with
    7 weights.
- **Also from them:** at 7, H_b(n) is divisible by 7 for n = 2, 4, 5, 6 at every level, beyond Watson's power (counted;
  not examined further).

## When a form first appears, and its slope (14 Nov 2026; a reading, then a test)

**What 7 and 11 showed.** A form's slope is one higher at the weights whose step falls on certain residues mod ℓ:
- 7's quiet form (slope 4) at x ≡ 2, 3, 4, 5 (mod 7);
- 11's quiet form (slope 3) at x ≡ 1, 2, 3 (mod 11).

There are as many such residues as its slope, and they are the steps just before it first appears (x = 6 at 7, x = 4 at 11).

**The reading (derived; as far as I recall, the shape of Bergdall and Pollack's "ghost zeros"; not checked against their
papers):**
- A form that first appears at step a with slope s has s marked steps, a − s, …, a − 1.
- At any later weight whose step is congruent to a marked step mod ℓ, its slope is s + 1. Elsewhere it is s.
- So when a form appears and how steep it is are one fact: it waits one step per unit of slope.
- Two older counts already fit:
  - 5's law form (a = 4, s = 2) has slope 3 at x = 7 and 8, which are ≡ 2, 3 (mod 5);
  - 7's law form (a = 2, s = 1) has slope 2 at x = 8, which is ≡ 1 (mod 7).
- For a slope of ℓ or more the marked steps would cover every residue, so the reading is stated for s < ℓ, and for forms
  alone in their system.

**PREDICTED before the test** (each form named by its system on lane ℓ; console only):

| form | a, s | slope s + 1 predicted at | slope s predicted at |
|---|---|---|---|
| 5's law form | 4, 2 | x = 12, 13 (weights 50, 54) | x = 9, 10, 11 |
| 7's law form | 2, 1 | x = 15 (94) | x = 16, 17 |
| 7's quiet form (5, 6, 0) | 6, 4 | x = 16, 17 (100, 106) | x = 13, 14, 15 |
| 11's law form | 2, 2 | x = 11, 12 (118, 128) | x = 9, 10, 13, 14 |
| 11's quiet form (0, 4, 1) | 4, 3 | x = 12, 13, 14 and 23, 24, 25 | x = 9, 10, 11, 15, 22, 26 |
| 13's quiet form (11, 8, 0) | 2, 2 | x = 13, 14 (166, 178) | x = 12, 15 |
| 17's quiet form (4, 2, 9, 7) | 2, 2 | x = 17, 18 (286, 302) | x = 16, 19 |
| 29's quiet form (25, 20, 12, 5) | 1, 1 | x = 29 (838) | x = 28, 30 |

**Counted** (`fading_forms5.py … --labels`, console only; seconds each, 29 in 4 minutes). **As PREDICTED at all 37
weights.** Each form is followed in its own system; where a steeper form has joined the system, the form's slope is the
lower one.

| form (its system) | slopes counted, by step x |
|---|---|
| 5's law form (2, 2, 1) | 2, 2, 2 at x = 9, 10, 11; **3, 3** at x = 12, 13 |
| 7's law form (4, 5, 1) | 1, 1 at x = 13, 14; **2** at x = 15; 1, 1 at x = 16, 17 |
| 7's quiet form (5, 6, 0) | 4, 4, 4 at x = 13, 14, 15; **5, 5** at x = 16, 17 |
| 11's law form (7, 3, 7) | 2, 2 at x = 9, 10; **3, 3** at x = 11, 12; 2, 2 at x = 13, 14 |
| 11's quiet form (0, 4, 1) | 3, 3, 3 at x = 9, 10, 11; **4, 4, 4** at x = 12, 13, 14; 3 at x = 15; 3 at x = 22; **4, 4, 4** at x = 23, 24, 25; 3 at x = 26 |
| 13's quiet form (11, 8, 0) | 2 at x = 12; **3, 3** at x = 13, 14; 2 at x = 15 |
| 17's quiet form (4, 2, 9, 7) | 2 at x = 16; **3, 3** at x = 17, 18; 2 at x = 19 |
| 29's quiet form (25, 20, 12, 5) | 1 at x = 28; **2** at x = 29; 1 at x = 30 |

- **So the link is exact where tested:**
  - A form whose slope is s waits s marked steps, the s steps just before it first appears.
  - From then on its slope is s, except at weights whose step is congruent to a marked step mod ℓ, where it is s + 1.
  - When it appears and how steep it is are one fact.
- **The ordinary (slope-0) forms never step up** (for example 13's (3, 5, 5) and 17's (16, 12, 15, 12) at every weight
  above). They have no marked steps.
- **Where the marked steps start differs:**
  - at step 0 (the weight ℓ − 3) for the quiet forms of 13 to 61 and for 11's law form;
  - at step 1 or 2 for 5's and 7's law forms and for 7's and 11's quiet forms.
  - At 13 to 61 the start at step 0 is what makes "slope = the step it first appears at" (and the fast/slow rule) hold.
- **In ghost terms (a reading, as I recall Bergdall and Pollack's picture):** the marked steps would be the weights where
  the form's system already has room at level 6ℓ but none yet at level 6. Testing that directly needs the forms of level
  6ℓ, split by system.
- **Not covered:** forms that share their system from the start, and slopes of ℓ or more. In those cases the steps above
  move in other ways, for example the second form in 11's law system (4, then 3 at x = 11, then 5 from x = 12).

### Why the waiting starts where it does: the forms of level 6ℓ (14 Nov 2026)

**The question.** The loud notes do not start the waiting (above). The ghost picture, as I recall Bergdall and Pollack's
"ghost series", puts the waiting in each form's own system:
- the weight k is a waiting step of the first form of a system exactly when the system has no form of level 6 there, but
  at least 2 forms of level 6ℓ that are new at ℓ;
- in their terms, d_k < 1 < d_k^Iw − d_k, where d_k counts level 6 and d_k^Iw counts level 6ℓ, the old forms giving
  2·d_k of it.

**The tool** (`level6l_counts.py`, new; level 6 is unchanged in `forms6_flint.py`):
- **The space:** forms of level 6ℓ over the ℓ-adic integers, in three steps.
  1. A pool of integral forms: products X(z)·Y(ℓz) of level-6 forms, then the same times one or two extra factors. The
     factors are weight-2 Eisenstein series E(z) − d·E(dz), products a(d1·z)·a(d2·z), and brackets
     a(d1·z)·θa(d2·z) − θa(d1·z)·a(d2·z), where a = Σ q^(x²+xy+y²) is the weight-1 theta series and θ = q·d/dq.
  2. From the pool, as many forms as the space's dimension that are independent over Q. The script stops if the pool
     falls short.
  3. Saturation: a combination that is 0 mod ℓ is divided by ℓ, until the forms are independent mod ℓ. Then they are a
     basis over Z_ℓ.
- **Changed while building, before any counted run** (each found by a trial at ℓ = 7, weights 4 and 10):
  - The plan was products alone. They cannot be a basis over Z_ℓ for k ≥ ℓ − 1 (derived): mod ℓ they fall into the
    forms of level 6 and weight kℓ, a space of dimension kℓ + 1, smaller than dim M_k(Γ0(6ℓ)) = kℓ + k − ℓ + 3. Hence
    the saturation. At weight 10 it took one round of 11 divisions.
  - At weight 4 the Eisenstein products span only 25 of 28 dimensions over Q. A cusp form meets them only through its
    L-value at the centre, which is 0 for the forms of sign −1 (Rankin–Selberg). Hence the theta factors, which give
    all 28.
- **The counts:** the class's signs at 2 and 3 (exact, to ℓ^N), split by the systems on lane ℓ (T_q mod ℓ, q = 5, 7, 11,
  13 other than ℓ). Each system is counted at level 6 and at level 6ℓ; new at ℓ = level 6ℓ − 2·level 6.
- **Two checks:**
  1. the rank of the products equals the dimension formula for M_k(Γ0(6ℓ));
  2. summed over all four sign classes, the counts at level 6ℓ equal dim S_k(6ℓ) − 2·dim S_k(3ℓ) − 2·dim S_k(2ℓ) +
     4·dim S_k(ℓ), the forms new at 2 and 3.

**PREDICTED before the run** (the forms of "When a form first appears, and its slope"; their level-6 count is 0 at every
step listed):

| form (system) | waiting steps: at least 2 new at ℓ | earlier steps: at most 1 new at ℓ |
|---|---|---|
| 5's law form (2, 2, 1) | x = 2, 3 (weights 10, 14) | x = 1 (weight 6) |
| 7's law form (4, 5, 1) | x = 1 (10) | x = 0 (4) |
| 7's quiet form (5, 6, 0) | x = 2, 3, 4, 5 (16, 22, 28, 34) | x = 0, 1 (4, 10) |
| 11's law form (7, 3, 7) | x = 0, 1 (8, 18) | — |
| 11's quiet form (0, 4, 1) | x = 1, 2, 3 (18, 28, 38) | x = 0 (8) |
| 13's quiet form (11, 8, 0) | x = 0, 1 (10, 22) | — |

If a waiting step has fewer than 2, or an earlier step has 2 or more, the ghost reading is wrong as stated.

**Also PREDICTED, for every system of the class, not only these six** (written after the trial at ℓ = 7, x = 0 and 1,
before the counted runs):
- Call z a zero of a system if, at step z before the system's first form appears, the system has at least 2 forms new
  at ℓ.
- The system's smallest slope at a later step x is then Σ over its zeros z of (1 + v_ℓ(x − z)).
- This is the marked-steps rule again, read from level 6ℓ alone, with no slope put in by hand.
- Runs: ℓ = 5 at x = 1…9 (weight 2 left out), ℓ = 7 at x = 0…9, ℓ = 11 at x = 0…5, ℓ = 13 at x = 0…3.

**Counted** (`level6l_counts.py`; outputs `level6l_counts_5_x1-9`, `_7_x0-9`, `_11_x0-5` and `_13_x0-3`, each `.txt`
and `.jsonl`; 17, 167, 179 and 75 seconds). **As PREDICTED at all 19 steps of the table, and for every system: 30 of 30
slopes.**

- **The checks held at all 29 weights:**
  - the pool spans the formula's dimension over Q;
  - the saturation ends in a basis over Z_ℓ (up to 9 rounds; the products' span had index up to 7^267);
  - U_2, U_3 and T_q are right on all B + 1 coefficients, and they commute;
  - the four sign classes add up to the formula, at level 6ℓ and at level 6 (58 of 58);
  - every system of level 6 is found at least twice at level 6ℓ.
- **The counts** (forms at level 6 / at level 6ℓ, by step x; **bold** = the waiting steps; COUNTED):

| form (system) | x = 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 |
|---|---|---|---|---|---|---|---|---|---|---|
| 5's law form (2, 2, 1) | (weight 2) | 0/1 | **0/2** | **0/3** | 1/4 | 1/5 | 1/6 | 1/7 | 1/8 | 1/9 |
| 7's law form (4, 5, 1) | 0/1 | **0/2** | 1/3 | 1/4 | 1/5 | 1/6 | 1/7 | 1/8 | 1/9 | 1/10 |
| 7's quiet form (5, 6, 0) | 0/0 | 0/1 | **0/2** | **0/3** | **0/4** | **0/5** | 1/6 | 1/7 | 1/8 | 1/9 |
| 11's law form (7, 3, 7) | **0/2** | **0/4** | 1/6 | 1/8 | 1/10 | 1/12 | | | | |
| 11's quiet form (0, 4, 1) | 0/1 | **0/2** | **0/3** | **0/4** | 1/5 | 1/6 | | | | |
| 13's quiet form (11, 8, 0) | **0/2** | **0/4** | 1/6 | 1/8 | | | | | | |
| 13's ordinary form (3, 5, 5) | 1/2 | 1/4 | 1/6 | 1/8 | | | | | | |

- **The slopes from level 6ℓ alone (DERIVED), each equal to the slope counted at level 6:**
  - 5: 2, 2, 2, 3, 3, 2 at x = 4…9; the 3s fall at x ≡ 2, 3 (mod 5).
  - 7's law form: 1 at x = 2…9, except 2 at x = 8 ≡ 1 (mod 7).
  - 7's quiet form: 4, 4, 4, 5 at x = 6…9; the 5 falls at x = 9 ≡ 2 (mod 7).
  - 11's forms: 2 at x = 2…5 and 3 at x = 4, 5. 13's forms: 2 at x = 2, 3, and 0 at x = 0…3.
- **What the counts show (COUNTED):**
  - Level 6ℓ gains a system's forms on a straight line, one or two per step. That is one per step for 5's form, both of
    7's, 11's quiet form and the three waiting systems below; two per step for 11's law form and 13's (11, 8, 0) and
    (3, 5, 5).
  - The line does not bend when the form reaches level 6. It goes on, now counting the form twice (as old forms) along
    with the new ones.
  - The waiting starts at the first step where the line reaches 2:
    - one per step from 0 at step 0 gives step 2 (5's form, 7's quiet form);
    - one per step from 1 gives step 1 (7's law form, 11's quiet form);
    - two per step from 2 gives step 0 (11's law form, 13's quiet form).
- **In plain words:**
  - Before a form reaches level 6, its system is already one level down, at level 6ℓ. There it gains one or two forms
    per step.
  - The waiting starts as soon as level 6ℓ holds two of them.
  - Each step of waiting adds one to the slope the form has when it arrives.
- **Systems still waiting at the end of the runs (COUNTED), and what the reading PREDICTS for them** (written before
  checking at level 6; it assumes the straight lines go on):
  - 11's (2, 6, 9) and (6, 4, 1): 0/1, 0/2, 0/3, 0/4, 0/5 at x = 1…5, waiting since x = 2. Whichever step a they reach
    level 6 at, their first slope is predicted to be a − 2 (for a ≤ 12).
  - 13's (2, 12, 5): 0/2, 0/4, 0/6, 0/8, waiting since x = 0. Predicted first slope: a (for a ≤ 12).
  - 13's (7, 6, 3): 0/1, 0/2, 0/3 at x = 1…3, waiting since x = 2. Predicted first slope: a − 2.
- **Checked at level 6: all four as PREDICTED** (COUNTED, seconds each):
  - Runs: `fading_forms5.py 11 12 8 --slope=1 --from=6 --labels --keep`, and the same for 13 from x = 4. Outputs:
    `fading_forms5_11_s1_x6-12.txt` and `fading_forms5_13_s1_x4-12.txt`.
  - 11's (2, 6, 9) and (6, 4, 1) reach level 6 at x = 10 (weight 108), both with slope 8 = 10 − 2.
  - 13's (2, 12, 5) reaches it at x = 4 (weight 58), with slope 4.
  - 13's (7, 6, 3) reaches it at x = 12 (weight 154), with slope 10 = 12 − 2.
  - So level 6ℓ gives a form's slope before the form exists at level 6. Its slope is the step at which it arrives,
    minus the step at which level 6ℓ first held two forms of its system.
- **Still open:** why level 6ℓ grows by one or two per step, and at which step a form reaches level 6 (the end of the
  waiting). Both are dimension counts of the two levels, which the ghost picture takes as given. Also untested here: the
  second forms of a system (for example 13's (11, 8, 0) gains a slope-6 form at x = 8). In the ghost picture they have
  zeros of their own.

### A relayed review of the level-6ℓ test (14 Nov 2026), assessed point by point

- **"The rule is Bergdall and Pollack's ghost series; the mechanism is not new."** Legitimate. Our text already called it
  their picture "as I recall". The review cites Bergdall and Pollack, "Slopes of modular forms and the ghost conjecture"
  (IMRN 2017) and its part II (QUOTED from the review, not read here). From now on the counting rule (zeros where level 6
  has nothing and level 6ℓ has at least 2 new; each adds 1 + v_ℓ(k − k′)) is theirs, and our part is the test.
- **The review's own check at 7.** The ghost slopes from our level-42 counts are 4, 4, 4, 5, 5, 5, 5 at x = 6 … 12, the
  same as its PARI run. Consistent with ours: 4, 4, 4, 5 at x = 6 … 9 above, and slope 5 at x = 9 … 12 in "The quiet
  notes of 7 and 11, digit by digit".
- **"The one-or-two growth is a dimension-formula fact. Liu, Truong, Xiao and Zhao note that d^Iw increases by 2 per
  step."** Partly. It is a dimension fact, but which systems gain 1 and which gain 2 is not said by "+2". The section below
  derives it: 2 for a system whose seed has weight 4 … ℓ, 1 for a seed of weight ℓ + 1. If that holds, the "+2" is the
  generic case, and our "+1" systems fall outside it.
- **"Your forms are ordinary, so locally reducible; whether the conditions of the proofs hold is the question."** The
  question is right, the wording not: the quiet forms have slopes 1 to 10. What matters is the residual representation at
  ℓ. A single seed (one Serre weight) means it is reducible there. The seeds below say which case each system is.
- **"What to claim: the link to the settling rate (r = s_Q − s_L), the ⌊ℓ/24⌋ first-step count, the computations."**
  Agreed. The first-step gains in our record (0, 0, 1, 1, 1 at 19, 23, 37, 41, 43, and 25 at 607) all equal ⌊ℓ/24⌋.
  Our record states the rule as the first step's own count. The closed form ⌊ℓ/24⌋ is the review's and is not checked
  here at every prime.
- **"Next: a plain ghost-series computation compared with the 30 slopes."** Our "same reading for every system" is
  exactly the ghost series' first coefficient, so the 30 of 30 is that test. The full series also predicts the second and
  later forms of each system, a stronger test. It is next after the seeds.

## Why level 6ℓ gains one or two forms per step (14 Nov 2026)

**First, a correction to what the class counts.** The trial of `level6l_counts2.py` (ℓ = 7, x = 0 … 2) counts two
things at each level:
- the class itself, forms new at 2 and 3 with the class's signs exactly;
- its whole generalized eigenspace mod ℓ, which adds every form whose U_2 and U_3 are congruent to the class's mod ℓ.
For 7's (5, 6, 0) they differ: the class has 0, 1, 2 forms at level 42, the whole eigenspace 0, 4, 8. The extra forms
are old at 2 or 3, congruent to the class mod 7. The ghost test above used the class itself, which is the clean count.

**The reading (DERIVED from standard facts, recalled, not checked against papers):**
1. **Forms new at 2 and 3 live on a Shimura curve.** By Jacquet–Langlands they are the forms on the Shimura curve of
   discriminant 6. That curve has no old forms at 2 and 3. So the class is a whole generalized eigenspace there, with
   nothing congruent mixed in.
2. **Level ℓ on top is a coefficient system (Shapiro's lemma).** Forms of level Γ0(ℓ) on that curve, of weight k, mod ℓ,
   are forms of level 1 at ℓ with coefficients Ind_B^G Sym^(k−2). Here G = GL2(F_ℓ) and B is its upper-triangular group.
3. **Counting splits over composition factors.** Counted by system, the forms add up over the composition factors, the
   Serre weights σ = Sym^r ⊗ det^m (0 ≤ r ≤ ℓ − 1). Each factor contributes μ(σ): the number of class forms at level 6 of
   weight r + 2 in the system twisted by ω^(−m) (T_q → q^(−m)·T_q, U_2 → 2^(−m)·U_2, U_3 → 3^(−m)·U_3). I call the forms
   with μ(σ) > 0 the system's **seeds**.
4. **Every step adds the same pieces.** Ind_B^G Sym^(k−2) is built from k − 1 principal series Ind(a^(k−2−i)·d^i),
   i = 0 … k − 2. A step (k → k + ℓ − 1) adds ℓ − 1 more, one for every character: the same set at every step. So the
   count at level 6ℓ grows by the same amount at every step, which is the straight line.
5. **Why 2 or 1.** In those ℓ − 1 principal series, each Serre weight with 1 ≤ r ≤ ℓ − 2 appears twice. Each with r = 0
   or r = ℓ − 1 appears once. Weight 2 (r = 0) has no forms at level 6. So the gain per step is
   2 × (the system's seeds of weight 4 … ℓ) + 1 × (its seeds of weight ℓ + 1).
6. **Why level 6 waits.** Level 6 is the same sum over the factors of Sym^(k−2) alone. Those come in less than one
   principal series per step, with shifting characters, so a system's forms reach level 6 only now and then.

**PREDICTED, before `serre_weights.py` runs.** From the class forms of weights 4 … ℓ + 1 at level 6 alone, it computes:
- every system's class count at level 6ℓ at every weight counted (`level6l_counts` and `level6l_counts2` at 5, 7, 11, 13,
  17, 19, 23);
- every system's class count at level 6 at those weights;
- the same pair for the whole generalized eigenspace, from all forms of weights 4 … ℓ + 1 at level 6 (Γ0(6) in place of
  the Shimura curve, so old forms included).
A single mismatch means the reading is wrong or incomplete.

**A correction, found by a trial before the full run** (`serre_weights.py 7` on `level6l_counts_7_x0-9.jsonl`):
- 7's law system (4, 5, 1) came out exactly at all 10 weights and both levels. Its one seed has weight 8 = ℓ + 1, so it
  gains 1 per step.
- 7's quiet system (5, 6, 0) was predicted to gain 2, and gains 1. It is of Eisenstein type mod 7: T_q ≡ q⁴ + q⁵ for
  q = 5, 11, 13, and the class's U_2 ≡ 2 = 2⁴ and U_3 ≡ 4 = 3⁴ (mod 7). That is the Hecke action on the constant factor
  det⁴.
- The fix (DERIVED). On the compact Shimura curve the forms counted are minus the Euler characteristic h⁰ − h¹ + h²,
  which adds up exactly over the factors. Only the one-dimensional factors det^m′ have h⁰ and h²: each is 1 there and
  h¹ = 0. So every factor det^m′ carrying the class's numbers takes one form away. With step 5 it gives: gain 2 − 1 = 1.
- Checked by hand on (5, 6, 0): it gives level 6 counts 0, 1, 1, 1, 1 at x = 5 … 9, exactly as counted. The level-6
  form arrives at x = 6 because det⁴ is missing from Sym^38 but present at x = 5 and 7 … 9.
- Everything else in the full run is a fresh test of the reading with this correction: the other systems and primes. The
  whole-eigenspace counts are predicted only for systems not of Eisenstein type, since Γ0(6) also has Eisenstein series.

**Counted** (`level6l_counts2.py` at 5, 7, 11, 13, 17, 19 and 23; `serre_weights.py`, outputs `serre_weights_ℓ.txt`;
seconds each). **As PREDICTED, every count at every weight:**
- the class: **132 of 132** (system, weight) pairs, each at both levels at once;
- the whole generalized eigenspace: **88 of 88**;
- 39 systems in all, at 29 + 10 weights.
- The checks inside `serre_weights.py` all held: the composition factors of every principal series and of Sym^n, by
  Brauer characters on every semisimple class; and Ind_B^G Sym^(k−2) = Sym^(k−2) ⊗ Ind(1).
- `level6l_counts2.py` repeats `level6l_counts.py`'s class counts at 5, 7, 11 and 13 exactly.

**Each system's seed, and its gain per step** (COUNTED; the seed is a class form at level 6 of the weight shown, twisted
by ω^m):

| ℓ | gain 1 | gain 2 |
|---|---|---|
| 5 | (2, 2, 1): weight 4, Eisenstein type | — |
| 7 | (4, 5, 1): weight 8; (5, 6, 0): weight 6, Eisenstein type | — |
| 11 | (0, 4, 1) and (6, 4, 1): weight 12; (2, 6, 9): weight 10, Eisenstein type | (7, 3, 7): weight 8 |
| 13 | (7, 6, 3): weight 12, Eisenstein type | (11, 8, 0): 8; (2, 12, 5): 12; (3, 5, 5): 10 |
| 17 | (4, 0, 15, 2): weight 18; (5, 13, 6, 3): weight 16, Eisenstein type | five systems, weights 6 … 16 |
| 19 | two of weight 20; (1, 18, 18, 12): weight 18, Eisenstein type | six systems, weights 4 … 16 |
| 23 | two quadratic pairs of weight 24 (1 each); two of weight 22, Eisenstein type | eight systems, weights 8 … 20 |

- **So the answer (COUNTED at all 39 systems, DERIVED as above):** a system gains 2 forms per step at level 6ℓ when its
  one seed has weight 4 … ℓ. It gains 1 in two cases:
  - its seed has weight ℓ + 1 (the Steinberg weight, which a step brings in once rather than twice);
  - or it is of Eisenstein type mod ℓ, with seed weight ℓ − 1, and a constant factor det^m′ that carries its numbers
    takes one form back each step (2 − 1).
- Every system here has exactly one seed, of multiplicity one. That is why the gains are only 1 or 2.
- **Level 6 too:** the same seeds, summed over Sym^(k−2) alone, give every level-6 count, including the step at which
  each system's first form arrives. For example (5, 6, 0) at 7 arrives at x = 6 because the factor det⁴, which takes a
  form away, is missing from Sym^38. So the waiting's start and its end both come from the seeds.
- **The Eisenstein-type systems** include 5's law form (2, 2, 1), 7's quiet form (5, 6, 0) and 13's (7, 6, 3). Their
  seeds are forms of weight ℓ − 1 congruent mod ℓ to Eisenstein series. For example the class form of weight 6 at
  level 6, mod 7, has T_q ≡ 1 + q⁵ = σ₅(q) for q = 5, 11, 13 (the label's primes), as E₆ does.
- **One run stopped and was redone.** `level6l_counts2.py 17` stopped at weight 46 on its own check: a combination was
  0 mod 17 on the first 829 coefficients but not at one later one, which Sturm's bound rules out for true forms.
  - The same computation in a scratch copy (same code, same inputs, 4 threads) passed.
  - Redone single-threaded (`level6l_counts2_17_x2-2`), it passed every check.
  - The likeliest cause is a one-off fault in a threaded routine. The checks are there to catch exactly that.

### The full ghost series from the seeds, against every slope (14 Nov 2026; the review's next run, made stronger)

**The computation (`ghost_series.py`, new):**
- From the seeds alone, it computes each system's class counts d_z (level 6) and D_z (level 6ℓ) at every step z, as in
  `serre_weights.py`.
- Bergdall and Pollack's ghost series (as the review states it): the coefficient g_i has a zero at step z of multiplicity
  min(i − d_z, D_z − d_z − i) whenever d_z < i < D_z − d_z.
- At step x, v_ℓ(g_i) = Σ over the zeros of their multiplicity × (1 + v_ℓ(x − z)). The first d_x slopes of the Newton
  polygon of the points (i, v_ℓ(g_i)) are the predicted slopes of the system's d_x forms at level 6.
- **Counted beside it:** the class at level 6, split by system, and the slopes of T_ℓ from its characteristic
  polynomial.

**PREDICTED before the run:**
- Every slope of every class form at level 6, at every step x = 0 … X: X = 25 at 5, 20 at 7, 14 at 11, 12 at 13, 9 at
  17, 8 at 19 and 7 at 23.
- This includes the second and later forms of a system, which the first-coefficient test above did not reach. Examples:
  11's (7, 3, 7) gains a second form from x = 6, and 13's (11, 8, 0) one at x = 8.
- The review notes that the ghost conjecture's original statement fails in some locally reducible cases. The
  Eisenstein-type systems are reducible mod ℓ, so a mismatch there would be informative rather than a fault of the tool.

**Counted** (`ghost_series.py`, outputs `ghost_series_ℓ_x0-X.txt` and `.jsonl`; 3 to 10 seconds each). **As PREDICTED:
225 of 225 (system, step) pairs, every slope of every form**, at 5, 7, 11, 13, 17, 19 and 23, with no inexact cases.
Steps 0 … 25, 20, 14, 12, 9, 8 and 7, up to weight 174.
- **Second and third forms included** (ghost slopes from the seeds = slopes counted at level 6):

| system | slopes, both ways |
|---|---|
| 7's (4, 5, 1) | 1, 7 (7 steps); 2, 6; 1, 7, 13; 1, 7, 14 (2 steps) |
| 7's (5, 6, 0) | 4, 10; 4, 11; 5, 9; 5, 10; 5, 11 |
| 11's (7, 3, 7) | 2, 4 (5 steps); 3, 3; 3, 5; 2, 5; 2, 5, 12 |
| 13's (11, 8, 0) | 2, 6 (5 steps) |
| 13's (3, 5, 5) | 0, 9 (3 steps) |
| 13's (2, 12, 5) | 4, 4 (7 steps) |

- **The case the marked-steps reading did not cover now fits.** 11's law system's second form went 4, then 3 at
  x = 11, then 5 ("When a form first appears", "Not covered"). The ghost series gives exactly the pairs (2, 4), (3, 3),
  (3, 5), (2, 5) at x = 10 … 13, and (2, 5, 12) when a third form arrives at x = 14.
- **The Eisenstein-type systems fit as well**, with the corrected counts (5's law form, 7's (5, 6, 0), 11's (2, 6, 9),
  13's (7, 6, 3), and those of 17, 19, 23).
- **So every slope of the class at level 6, at every weight tested, follows from the seeds alone.** The seeds are a few
  forms of weight at most ℓ + 1 at level 6. From them come the counts at both levels (Serre weights) and then every slope
  (Bergdall and Pollack's ghost series).
- **What is ours and what is not.** The ghost series is Bergdall and Pollack's. The review says Liu, Truong, Xiao and
  Zhao proved it under a genericity condition (QUOTED, not read here). Whether these systems meet that condition is not
  settled here; they match anyway, including the Eisenstein-type ones. Ours: the seeds of the partition class, the
  Eisenstein-type correction, the tool for level 6ℓ, and the test. (Settled in the next section.)

## Do the families meet the conditions of the proofs? (14 Nov 2026)

**What the papers require** (QUOTED, as read through arXiv's and ar5iv's HTML pages; theorem numbers as given there):
- **Liu, Truong, Xiao and Zhao**, "Slopes of modular forms and geometry of eigencurves" (arXiv:2302.07697v4; Annals).
  Theorem 1.3, the ghost conjecture, assumes:
  - p ≥ 11;
  - r̄ absolutely irreducible;
  - r̄ restricted to Gal(Q̄_p/Q_p) reducible (split or not) and "very generic", Serre weight Sym^a ⊗ det^b with
    2 ≤ a ≤ p − 5.
  - The spaces are S_k(Γ0(Np), ψ) localized at r̄, with no conditions at the primes dividing N stated.
  - Remark 1.4(2): "We expect that Theorem 1.3 also holds for a = 1, p − 4, p − 3, as well as a = 0 and r̄_p très
    ramifié" (the last with a modified formulation). For a = p − 3 they are "not sure" one of their main results holds.
  - Remark 1.4(3): for reducible r̄, "the only difference is some additional dimension computation".
- **Liyan Wang**, "The Local Ghost Theorem in the Très Ramifié Case" (arXiv:2609.16714, Sep 2026). It proves the local
  ghost theorem for p ≥ 11 and the Steinberg weight Sym^(p−1) (très ramifié). No global statement was found in what was
  read.
- **Bergdall and Pollack**, "Slopes of modular forms and reducible Galois representations: an oversight in the ghost
  conjecture" (arXiv:2110.07973; Proc. AMS 2022).
  - The original conjecture excluded ρ̄ ≅ ω^j ⊕ ω^(j+1), the "twists of E2", "for convenience".
  - The salvaged Conjecture 3.1 covers every **regular** ρ̄ and "includes ρ̄ ≅ ω^j ⊕ ω^(j+1)" (Remark 3.2).
  - Regular means (ρ̄|_D)^ss is not a twist of ind(ω₂^s). Equivalently it is neither irreducible nor a twist of an
    unramified representation with Frobenius trace zero (Definition 2.1, Proposition 2.2). The irregular ones always
    fail (Proposition 3.8).
  - Their dimensions come from weights 2 … p + 1 and twists, "plus Eisenstein multiplicities" when ρ̄ is Eisenstein
    (Section 3.1). That is our seeds and our det correction, in their setting.

**Measured** (`seed_types.py`, outputs `seed_types_ℓ.txt`, at 5, 7, 11, 13, 17, 19, 23; COUNTED):
- **Every seed is ordinary:** a_ℓ ≢ 0 mod ℓ in 39 of 39. For seeds of weight 4 … ℓ + 1 that means ρ̄ is reducible at ℓ
  (a non-ordinary form of such a weight has ρ̄ irreducible there). So the review's "ordinary, so locally reducible" holds
  for the seeds, though not for the quiet forms at high weight.
- **Truly Eisenstein, on 23 primes up to 97:** 8 families, all of the form ω^(j+1) ⊕ ω^j. They are exactly the ones whose
  det factor takes a form back.
- **11's (7, 3, 7) is not Eisenstein.** The 3-prime test had matched by chance.
- **Every family is regular in Bergdall and Pollack's sense** (DERIVED): it is reducible at ℓ, and its a is even, never
  ℓ − 2, so its two characters never differ by an unramified one.

**The families against the theorems** (DERIVED from the hypotheses above):

| kind | families | LTXZ Theorem 1.3 | other |
|---|---|---|---|
| very generic (2 ≤ a ≤ ℓ − 5), irreducible, ℓ ≥ 11 | 21: 1 at 11, 2 at 13, 4 at 17, 6 at 19, 8 at 23 | **covered** | — |
| Steinberg weight (a = 0, très ramifié): seed weight ℓ + 1 | 8: one at 7, two each at 11, 19, 23, one at 17 | no | Wang's local theorem at ℓ ≥ 11 |
| a = ℓ − 3, irreducible | 2: 13's (2, 12, 5), 17's (2, 13, 13, 9) | no: "expected", but for a = p − 3 the authors are also "not sure" one of their main results holds | — |
| a = ℓ − 3, Eisenstein ω^(j+1) ⊕ ω^j | 8: one at each of 5 … 19, two at 23 | not in the stated theorem; reducible r̄ needs "some additional dimension computation" (Remark 1.4(3)), which the seeds' correction supplies | Bergdall–Pollack Conjecture 3.1 includes them |

- **The partition tower's own forms:**
  - Covered by the theorem: 11's law form (7, 3, 7) (a = 6); 13's quiet form (11, 8, 0) (a = 6); 13's ordinary form
    (3, 5, 5) (a = 8); 17's quiet form (4, 2, 9, 7) (a = 6).
  - Not covered:
    - 5's law form (2, 2, 1): Eisenstein, and ℓ < 11;
    - 7's law form (4, 5, 1): Steinberg, and ℓ < 11;
    - 7's quiet form (5, 6, 0): Eisenstein, and ℓ < 11;
    - 11's quiet form (0, 4, 1): Steinberg, local theorem only.
- **All 39 fit the ghost series** (225 of 225 slopes above), including the 18 outside the proven range. That is evidence
  for Liu–Truong–Xiao–Zhao's expected cases (a = ℓ − 3, a = 0 très ramifié) and for Bergdall and Pollack's inclusion of
  ω^j ⊕ ω^(j+1). For those, our corrected dimensions are the "additional dimension computation".
- **Caveats:**
  - The theorem is stated for the whole r̄-part of S_k(Γ0(Np)). Our class, with fixed signs at 2 and 3, is a summand of
    it. The local theorem applies to "primitive" projective modules of type σ. A single seed of multiplicity one fits
    "primitive". That the class's module (from the Shimura curve of discriminant 6) is projective is the standard input,
    assumed here, not shown. (Closed by subtraction in the next section, without this assumption.)
  - "Covered" means the hypotheses as quoted hold for the family, not that we checked the proofs.

### Three more relayed reviews of the families' standing, assessed point by point (14 Nov 2026)

**What the reviews said, and the verdict on each:**
- **"Your seeds are Serre weights; the 1-versus-2 rule is the known Serre-weight bookkeeping."** Legitimate. Bergdall
  and Pollack's part II gives the per-step growth as 2m₁, m₁ + m₂, 2m₁ + 2m₃ or m₁ + m₂ + 2m₃, by whether ρ̄ is
  split at p and whether its weight is ≡ 2 mod p − 1, with the m_i from weights ≤ p + 1 (Corollary 6.17). This is
  QUOTED from a summary of their §6 read through ar5iv, not verbatim. The rule is theirs; ours are the computation and
  the check.
- **"The Eisenstein-type families are exactly the case Bergdall and Pollack excluded and later corrected."**
  Legitimate, and checked in the papers:
  - part II §7.1 (verbatim): they "cannot efficiently determine the slopes of the abstract ghost series for ρ̄ = 1 ⊕ ω
    (or its cyclotomic twists)", "omitted making a conjecture in this case", and "know of no reason" for it to fail;
  - the correction paper's Conjecture 3.1 includes them (Remark 3.2) and was checked for 1 ⊕ ω at levels 1 and 11.
  - Their part II's Lemma 6.1 (summary only) treats H²_c of the characters a^s·d^r, nonzero only when s ≡ r. Those are
    the characters whose factor det^s carries our correction.
- **"Three primes is evidence, not proof, that 7's system is Eisenstein."** Right for three primes. Resolved by a proof:
  `eisenstein_proof.py` (outputs `eisenstein_proof_ℓ.txt`), below.
- **"Ordinary: I was wrong; what matters is ρ̄."** Agreed. Settled above: all 39 seeds are ordinary.
- **"⌊ℓ/24⌋: counted directly at the 22 primes from 13 to 73, with 67, 71 and 73 predicted first; it fits 607."** The
  review's own count. It agrees with every first-step count in our record (19, 23, 37, 41, 43, 607). Not recounted here
  at the other primes.
- **"Check Theorem 1.3's exact range of a; 'expect a = ℓ − 3' needs a qualifier; the Eisenstein row undersells
  Remark 1.4(3)."**
  - The range was read verbatim: §1.2 of arXiv:2302.07697v4, "very generic (i.e. 2 ≤ a ≤ p − 5)". A review then
    confirmed it.
  - Both wordings are now fixed in the table above.
- **"The projectivity gap: recover the class by subtraction."** Legitimate. Argument and checks below.

**Proved: the 8 Eisenstein-type families are Eisenstein** (`eisenstein_proof.py`; COUNTED coefficients, DERIVED
conclusion; *in the surveys' labels: the coefficient check is computed, and the conclusion from Sturm's bound is derived*):
- **The method.** Let f be the seed (weight ℓ − 1, normalized). Let G be the weight-(ℓ + 1) Eisenstein series of level
  6 with U_2 = U_3 = 1. Then θ^(ℓ−1)·f ≡ θ^(ℓ−2)·G mod ℓ in every coefficient up to Sturm's bound ℓ² + ℓ − 2.
  - The comparison covers 29, 55, 131, 181, 305 and 379 coefficients at 5, 7, 11, 13, 17 and 19, and 551 for each of
    the two at 23.
  - Mod ℓ both sides are modular forms of weight ℓ² + ℓ − 2 on Γ0(6) (θ raises the weight by ℓ + 1), so they are equal.
  - Hence a_q(f) ≡ 1 + q^(ℓ−2) for every prime q ∤ 6ℓ, and ρ̄ ≅ ω^m ⊕ ω^(m−1) exactly.
- **Results.** 8 of 8 proven: ω³ ⊕ ω² at 5; ω⁵ ⊕ ω⁴ at 7; ω⁹ ⊕ ω⁸ at 11; ω¹¹ ⊕ ω¹⁰ at 13; ω¹⁵ ⊕ ω¹⁴ at 17;
  ω¹⁷ ⊕ ω¹⁶ at 19; ω¹⁰ ⊕ ω⁹ and ω²¹ ⊕ ω²⁰ at 23.
- **Control.** Every other seed fails the same comparison, at n = 5 or at a_2 and a_3.

**The direct-summand data** (from `level6l_counts2`'s two counts; COUNTED):
- **29 of the 39 families:** the class is its whole generalized eigenspace mod ℓ at both levels and every weight, so
  Hecke operators at 2 and 3 cut it out.
- **The other 10:** the 8 Eisenstein ones, 17's (16, 12, 15, 12) and 19's (12, 14, 15, 13). Each shares its numbers
  mod ℓ with old forms; the last two have twice the class's count.
- **On the Shimura curve of discriminant 6** every form is new at 2 and 3. There the Atkin–Lehner idempotents
  (1 ± w_p)/2 (ℓ odd) cut the class out for every family.

**Closing the gap by subtraction (the review's argument, checked here):**
1. **Within one ρ̄, the signs at 2 and 3 are fixed (DERIVED).**
   - At 2, a form new at 2 has Frobenius eigenvalues {α, 2α} with α fixed by its U_2 sign.
   - Another sign would need {α, 2α} = {β, 2β} with β ≠ α, so 3α ≡ 0: impossible for ℓ ≠ 3.
   - At 3 the same with {α, 3α}: 8α ≡ 0, impossible for ℓ ≠ 2.
2. **Checked** (`sign_check.py`, which runs `level6l_counts2.py`'s own functions; outputs `sign_check_ℓ.txt`): at all 7
   primes, 2 or 3 weights each, both levels, no ρ̄ occurs in two sign classes.
   - Two 4-prime labels at 19 (weight 34) did occur in both (+, +) and (−, −). Split by more primes (a scratch run up to
     q = 29), each pair is ρ̄ and its twist by the quadratic character mod 19.
   - That twist flips both signs, since (2|19) = (3|19) = −1. It keeps the 4-prime label, since a₁₃ ≡ 0 and only q = 13
     has (q|19) = −1 among them. It shows at q = 29: a₂₉ ≡ 5 against 14, and 7 against 12 (each is −a₂₉). So they are
     two different ρ̄.
3. **So within one ρ̄,** the forms at level 6ℓ are the class plus forms old at 2 or 3, and the U_ℓ-slopes subtract:
   slopes(class) = slopes S(6ℓ)_ρ̄ − 2·slopes S(3ℓ)_ρ̄ − 2·slopes S(2ℓ)_ρ̄ + 4·slopes S(ℓ)_ρ̄, as lists with
   multiplicity.
4. **Theorem 1.3 applies to each term**, with tame level 6, 3, 2 or 1 and the same ρ̄. Each gives the same ghost series
   G_σ (it depends only on σ), stretched by that level's multiplicity.
   - So the class has the slopes of G_σ with multiplicity m(6) − 2m(3) − 2m(2) + 4m(1).
   - That multiplicity is 1, which our counts show: one seed of multiplicity one, gaining 2 per step.
- **Hence (DERIVED):** for the 21 covered families, the class's slopes follow from Theorem 1.3 itself, with no
  projectivity assumption.
- **Not done:** the review's numerical form of this (the ghost series of the whole family at each lower level). Step 4
  makes it automatic once steps 1 and 2 hold.
- **What is left:** the theorem's own conditions, which hold as quoted, and its slope-zero part, which the paper treats
  separately (its Theorem 8.7, not read).

**Two lines for the write-up, from a further review (both legitimate, both already so):**
- **The Sturm comparison used the level-6 series, adjusted at 2 and 3.** G = (1 − 2^ℓ·B₂)(1 − 3^ℓ·B₃)E_(ℓ+1),
  normalized, with U_2 = U_3 = 1. It was compared at every n up to the bound except the multiples of ℓ, multiples of 2
  and 3 included.
  - The level-1 E_(ℓ+1) would fail at n = 2. Its θ^(ℓ−2) has a₂ ≡ (1 + 2^ℓ)/2 ≡ 3/2, while the seeds have U_2 ≡ 1/2.
- **Like with like.** Theorem 1.3 describes all of S_k(Γ0(Mℓ))_ρ̄: both ℓ-stabilizations (slopes s and k − 1 − s) and
  the forms new at ℓ. The subtraction is of those full lists (their characteristic series).
  - Our ghost series was built from the level-6ℓ counts. Its first d_x slopes, the smaller one of each pair, are what
    was compared with the level-6 slopes.
- The review saw arXiv v3 as the latest version of Liu–Truong–Xiao–Zhao. Ours read the HTML of v4
  (arxiv.org/html/2302.07697v4), where §1.2 gives the range quoted.

## Why exactly 13, 17, 19 and 23 settle fast (14 Nov 2026)

**The question.**
- At 13, 17, 19 and 23 the tower settles one digit per level; at 29 and every prime after it, one per two levels
  ("What sets the settling rate").
- The rule found there: for ℓ ≥ 13, fast exactly when the class gains no form at the first weight step, ℓ − 3 → 2ℓ − 4.
- A review counted that first-step gain as ⌊ℓ/24⌋ at the 22 primes 13 … 73. Why?

**The derivation (DERIVED, from the Serre-weight bookkeeping of "Why level 6ℓ gains one or two forms per step").**

Write c(w; ε₂, ε₃) for the number of newforms of level 6 and weight w with U_2 = ε₂·2^((w−2)/2) and
U_3 = ε₃·3^((w−2)/2) (characteristic 0).

1. **The class at level 6 is a sum over the composition factors of Sym^(k−2).** The class count at weight k equals the
   sum, over the factors Sym^r ⊗ det^m of Sym^(k−2)(F_ℓ²), of c(r + 2; signs twisted by m). Here U_2 ≡ 2^(−m)·(the
   class's U_2) mod ℓ, and likewise at 3. This holds whenever no factor det^m′ carries the class's numbers. That is the
   Euler-characteristic count on the Shimura curve, tested above at 132 of 132.
2. **The first two weights.**
   - At x = 0, Sym^(ℓ−5) is irreducible: the class itself.
   - At x = 1, Glover's recursion gives Sym^(2ℓ−6) = {Sym^(ℓ−5), Sym⁴ ⊗ det^(ℓ−5), Sym^(ℓ−7) ⊗ det} (dimensions
     (ℓ − 4) + 5 + (ℓ − 6) = 2ℓ − 5).
   - So for ℓ ≥ 11 (no factor of type det^m′ there), the twists work out to:

     **gain(ℓ) = c(ℓ − 5; (2|ℓ), (3|ℓ)) + c(6; +, +).**
   - c(6; +, +) = 0: the only newform of weight 6 at level 6 has signs (+, −) ("Sym⁴ ⊗ det⁵" in `serre_weights_7.txt`).
     So gain(ℓ) = c(ℓ − 5; (2|ℓ), (3|ℓ)), the newforms of level 6 and weight ℓ − 5 with one fixed pair of signs.
3. **Why 1/24.** dim S_w^new(Γ0(6)) grows like w/6, shared by four sign classes, so each grows like w/24.
   - The exact count is w/24 plus a part periodic in w. This is the shape of the trace formula: the Atkin–Lehner traces
     on level 6 come from elements of trace 0, ±2 (det 2), 0, ±3 (det 3) and 0 (det 6), whose eigenvalues are √Q times
     roots of unity of order 4, 8 or 12.
   - So c(w + 24; ε) = c(w; ε) + 1 for every even w ≥ 4, if the slope is exactly 1/24. Then gain(ℓ + 24) = gain(ℓ) + 1,
     since (2|ℓ) and (3|ℓ) depend only on ℓ mod 24. And gain(ℓ) = ⌊ℓ/24⌋ follows for every prime once it holds at one
     prime in each of the 8 classes mod 24.
4. **Which families are slow (DERIVED, from the ghost zeros).**
   - Take a family whose one seed is Sym^r ⊗ det^m with 1 ≤ r ≤ ℓ − 2 and m ≤ ℓ − 5. Its count at level 6ℓ starts at 2
     at step 0, since the seed appears twice in the principal series of Ind Sym^(ℓ−5). So it waits from step 0, and its
     slope is the step it arrives at.
   - That step is the first Glover layer holding the seed: layer m, or layer ℓ − 5 − m, which holds it one step later.
     So **slope = min(m, ℓ − 4 − m)**.
   - Slope 1 means m = 1 or m = ℓ − 5: exactly the two new factors at x = 1. So a family of slope 1 exists iff the class
     gains at the first step.
   - For ℓ ≥ 13 the other seeds (m > ℓ − 5, the Steinberg weight, the Eisenstein type) wait longer.
5. **So, for ℓ ≥ 13:** the tower is slow iff ⌊ℓ/24⌋ ≥ 1, iff ℓ > 24. Fast exactly at 13, 17, 19, 23.

**PREDICTED before `first_step.py` runs** (exact counts c(w; ε) for even w ≤ 610 and all four sign pairs, by ranks
modulo two large primes):
- (a) c(w + 24; ε) = c(w; ε) + 1 for every w and ε.
- (b) gain(ℓ) = c(2ℓ − 4; +, +) − c(ℓ − 3; (2|ℓ), (3|ℓ)) equals c(ℓ − 5; (2|ℓ), (3|ℓ)) at every prime 11 ≤ ℓ ≤ 307.
- (c) Both equal ⌊ℓ/24⌋, also at 607 (25, the record's count) and at every prime up to 613 through c(ℓ − 5).
- (d) c(6; +, +) = 0.
- (e) The slope rule min(m, ℓ − 4 − m) fits every generic family of `seed_types` (twists m from `serre_weights`)
  against the first slopes counted in `ghost_series`.

**Counted** (`first_step.py`, output `first_step_w610.txt` and `.jsonl`, 524 seconds; `slope_rule.py`, output
`slope_rule_200.txt`, seconds). **As PREDICTED, all five:**
- **Checks:**
  - c(w; ε) was counted for every even w = 4 … 610 and all four sign pairs. The two primes agree at every weight.
  - The four classes add up to the newform dimension dim S_w(6) − 2 dim S_w(3) − 2 dim S_w(2) + 4 dim S_w(1) at every
    weight.
- **(a)** c(w + 24; ε) − c(w; ε) = 1 at every weight and every sign pair.
- **(b)** The reduction holds at every prime where the direct count reaches (11 … 307): 59 of 59.
- **(c)** gain(ℓ) = ⌊ℓ/24⌋ at all 108 primes from 11 to 613, including 25 at 607 (the record's own count, from 607's
  first step). The first prime of each class mod 24 fits (73, 29, 31, 11, 13, 17, 19, 23), so with (a) it holds at
  every prime of every class.
- **(d)** c(6; +, +) = 0. The weight-6 newform has signs (+, −).
- **(e)**
  - From the representation theory alone, at every prime 13 … 200 and for every possible single seed: the slope is
    min(m, ℓ − 4 − m) for the generic seeds with m ≤ ℓ − 5.
  - Slope 1 occurs exactly for Sym^(ℓ−7) ⊗ det and Sym⁴ ⊗ det^(ℓ−5).
  - Every other seed (Steinberg, Eisenstein type, late twists) has slope at least 2.
  - The families actually found at 13, 17, 19, 23: predicted arrival step and slope against those counted, 23 of 23 (the
    others arrive beyond the counted steps).
- **The loud notes:** c(ℓ − 3; (2|ℓ), (3|ℓ)) ≥ 1 at every prime 13 … 613. So from 13 on the loud note is an ordinary
  form at weight ℓ − 3 (slope 0); there are 1, 1, 1, 1, 1, 1, 2, … of them at 13, 17, 19, 23, 29, 31, 37, ….

**So, the answer (DERIVED, every step checked as above):**
1. At every prime ℓ ≥ 13 the class's first step gains exactly ⌊ℓ/24⌋ forms: the newforms of level 6 and weight ℓ − 5
   with the signs ((2|ℓ), (3|ℓ)).
2. Those forms' families, and only they, have slope 1 at the tower's weight. Every other family below the loud one
   waits at least two steps.
3. The loud note has slope 0, so the tower settles at rate 1 (one digit per two levels) when a slope-1 family exists,
   and at rate ≥ 2 otherwise.
4. Hence the tower settles fast exactly when ⌊ℓ/24⌋ = 0: **at 13, 17, 19 and 23, and at no prime after them.**
- **In plain words:**
  - Going up one weight step, the class picks up the level-6 newforms of weight ℓ − 5 with one fixed pair of signs.
  - Such forms come at a rate of one per 24 units of weight, so there are ⌊ℓ/24⌋ of them.
  - Below 24 there are none, so nothing with slope 1 joins the loud note and the tower settles at full speed.
  - From 29 on there is always at least one, and it slows the tower to one digit per two levels.
- **What it rests on:**
  - the Serre-weight count (derived; checked at 59 primes);
  - the trace formula's shape for (a) beyond weight 610 (standard);
  - the ghost slopes. For the slope-1 families these are proved: their seed Sym^(ℓ−7) ⊗ det has a = ℓ − 7, which is very
    generic for ℓ ≥ 11, so Liu–Truong–Xiao–Zhao cover them if irreducible. Above 23, that makes the slow side proved.
  - The fast side at 13 … 23 needs every other family's slope ≥ 2. Those slopes were counted directly at those four
    primes (225 of 225).
  - Not shown here: that every family of the class carries a share of the tower (the record's working assumption,
    checked at 13 … 31 in "What sets the settling rate").

### A further review of this proof, assessed, and its gap closed (14 Nov 2026)

**The review's points, and the verdict on each:**
- **Confirmed by the review independently.** Its own PARI counts at the 22 primes 13 … 73 equal ⌊ℓ/24⌋. The periodic
  shape (linear plus a part depending on the weight mod 24) is what it derived earlier.
- **"Two assumptions carry the weight: every family carries a share of the tower; and r = s_Q − s_L rests on A1 and A4,
  which are conjectured."** Legitimate.
  - STATEMENT.md marks A1 (the tower's two-level operator has characteristic polynomial P_ℓ) conjectured, and A4 (what
    fades is the next family) counted, at 13 and 17.
  - The result is restated below with these conditions.
- **"The slope-1 side needs ρ̄ irreducible; check it."** Legitimate. It also needs ρ̄ reducible at ℓ, that is the seeds
  ordinary. A non-ordinary seed would make the family irregular in Bergdall and Pollack's sense, where the ghost series
  always fails. Both are now checked (below).
- **"Weight −2 is not classical; quote the line."** The line, verbatim from arXiv:2302.07697v4 (23 Feb 2026; the review
  had seen v3): "Then for every w⋆∈𝔪_{ℂ_p}, the Newton polygon NP(C_r̄(w⋆,−)) is the same as the Newton polygon
  NP(G_σ(w⋆,−)), stretched in both x- and y-directions by m(r̄) times, except possibly for the their slope zero parts."
  - So it holds on the whole weight disc. Weight −2 is w = exp(−4ℓ) − 1, inside it.
  - The slope-zero part, which the slope-1 claim does not touch, is covered by their Theorem 8.7 (Remark 1.4(1), not
    read).
- **"Drop the slope rule from the proof."** Agreed.
  - The slow side needs only that the first-step families have slope 1 at −2, now shown directly (below).
  - The fast side rests on the measured rates.
  - The rule min(m, ℓ − 4 − m) stays as an observation.

**The checks** (COUNTED; `first_step_check.py`, output `first_step_check_29-613.txt`, 432 s; `first_step_ghost.py`,
output `first_step_ghost_613.txt`):
- **Every first-step system is ordinary and irreducible, at every prime 29 … 613.**
  - 103 of 103 primes have no exceptions at all.
  - Ordinary: T_ℓ (≡ U_ℓ mod ℓ) is invertible on every system, read from the coefficients a_(ℓp).
  - Irreducible: either some eigenvalue lies outside F_ℓ, or each of the ℓ − 1 possible pairs ω^a ⊕ ω^b is missed at
    some prime q ≤ 43, which is a proof.
  - *Why only these pairs (the conductor reason; added 14 Nov, later, at a relayed review's request):* if ρ̄ were
    reducible, its semisimplification would be χ₁ ⊕ χ₂ with cond(χ₁)·cond(χ₂) dividing the level's part away from ℓ,
    6 (semisimplifying never raises the conductor). The only nontrivial character of conductor dividing 6 is the
    quadratic one of conductor 3. It cannot occur once (det ρ̄ = ω^(w−1), the forms having trivial character) nor twice
    (that needs conductor 9). Characters unramified away from ℓ with values in F̄_ℓ are powers of ω. So χ₁ = ω^a and
    χ₂ = ω^b with a + b ≡ w − 1 (mod ℓ − 1), and a_q ≡ q^a + q^b (mod ℓ) at every q ∤ 6ℓ. One q where this fails rules
    out the pair.
  - At 29 … 47, where the first step adds a single form, that one form passes.
  - At 607 the first-step forms make two systems, of 8 and 17: 607's two quiet braids, found this morning by another
    route.
- **Their slope at the tower's weight −2 is exactly 1, at 103 of 103 primes.**
  - The ghost series of their Serre weight Sym^(ℓ−7) ⊗ det has v(g₁(w₋₂)) = 1, from one zero at step 0.
  - It has v(g₂(w₋₂)) = ℓ − 5 ≥ 24, since the family's second form arrives only at step ℓ − 4. So the first segment of the
    Newton polygon has slope 1.
- **A remark, elementary:** ⌊(ℓ − 1)/12⌋ = ⌊(ℓ + 12)/24⌋ + ⌊ℓ/24⌋ for every prime ℓ ≥ 5 (check the 8 residues mod 24).
  - That is, dim S_(ℓ−1)(SL₂(ℤ)) = R_ℓ + ⌊ℓ/24⌋.
  - *Reworded (14 Nov, later; a relayed review):* this is Hermite's identity ⌊x⌋ + ⌊x + ½⌋ = ⌊2x⌋ at x = ℓ/24 (for a
    prime ℓ ≥ 5, ⌊ℓ/12⌋ = ⌊(ℓ − 1)/12⌋). The numbers agree; they do not say that Folsom, Kent and Ono's space splits
    into the R_ℓ loud notes and the ⌊ℓ/24⌋ first-step families under the Hecke operators. That would need a
    Hecke-equivariant correspondence, not shown here. An earlier wording ("room for exactly") said more than this.
  - The earlier reading in "The speed of the towers" counted the same numbers but gave them the wrong mechanism.

**The result, restated with its conditions** (the review's labels):

| claim | status |
|---|---|
| The first step gains ⌊ℓ/24⌋ forms: the level-6 newforms of weight ℓ − 5 with signs ((2\|ℓ), (3\|ℓ)) | derived (Serre-weight count, periodic trace formula); counted at 108 primes 11 … 613, and by the review at 22 |
| Those families are ordinary and irreducible, of Serre weight Sym^(ℓ−7) ⊗ det (a = ℓ − 7, very generic) | counted, with a proof for each system, at every prime 29 … 613 |
| Their slope at the tower's weight −2 is exactly 1 | proved for 29 ≤ ℓ ≤ 613: Theorem 1.3 on the whole weight disc, and the ghost series |
| Every ℓ ≥ 29 is slow (one digit per two levels at most) | proved for 29 ≤ ℓ ≤ 613, given A1, A4 and that these families carry a share of the tower; measured on partition numbers at 29, 31 and 37 … 97 |
| 13, 17, 19 and 23 are fast | measured on partition numbers (13: Atkin–O'Brien (1967)) |
| A family's slope is min(m, ℓ − 4 − m) | an observation, checked for every seed at 13 … 200 |

- **Beyond 613**, "slow" also needs one first-step family that is ordinary and irreducible. That held at every prime
  checked, but is not proved in general.
- **The share condition is the one this record cannot prove.** It holds wherever the rate was measured (29 … 97).
  - *Corrected (14 Nov, later):* an earlier line here said 607's 25 slope-1 forms "were seen in its tower's quiet
    braids". That was wrong. They were found among the forms (at weights 1,210 to 2,422), not in 607's partition numbers,
    which reach only level 5. Whether 607's tower carries a share of them is not measured.

## The share condition, tested on partition numbers (14 Nov 2026)

**The condition.** The slow side needs the tower to carry a share of at least one first-step (slope-1) family. Folsom,
Kent and Ono and Boylan and Webb bound the space the tower lives in; they do not say which families it uses. At 607
Boylan and Webb found the tower missing one of its loud notes (stable rank R − 1), so a family can be left out, and the
condition has to be checked.

**Where a share can be seen (DERIVED).**
- At level 1, modulo ℓ, the tower is loud notes only (Yang: its space is the class at weight ℓ − 3). So every other
  family enters level 1 with a share divisible by ℓ.
- Two levels multiply a family's share by its U_ℓ number, whose valuation is the slope. At level 3 a family of slope s
  has a share of valuation at least 1 + s.
- Slope 1 belongs to the first-step families alone. The loud notes' own next families have slope ℓ − 4 (the ghost series
  of the loud seed Sym^(ℓ−5), every prime 13 … 613: 9 at 13, 25 at 29, 603 at 607). The first-step systems' next family
  has slope ℓ − 6. Every other family has slope at least 2.
- So at level 3, modulo ℓ³, the tower is the loud notes (units), plus ℓ² times the first-step families' shares, and
  nothing else. A share at depth one shows as exactly ℓ²; a missing share shows as nothing.
- Boylan and Webb's "sharp in all computed cases with ℓ ≥ 29" means b_ℓ(2) = 3 at 31 … 89. Any family with a share at
  depth one gives that, of any slope (13 has it too, through its slope-2 family). So their sharpness supports the share
  condition only at 29, the one slow prime where they computed modulo ℓ³.

**The test (no forms used).** Rows are squarefree t in the tower's class mod 24, columns are primes Q ≥ 5 (Q ≠ ℓ):
M_b[t, Q] = c(Q²·ℓᵇt) + (ℓᵇt|Q)·(−3|Q)·Q⁻²·c(ℓᵇt), with c(M) = p((M + 1)/24). This is the T(Q²) rule of
`partitions_split` and `partitions_multipliers`. Each family i adds λ_i(Q)·c_i(ℓᵇt), so M_b = (shares by t) ×
(eigenvalues by Q), and its invariant factors over ℤ/ℓ⁵ sort the families by the depth of their share: a unit for each
loud note, ℓ² for each first-step family that carries a share.

**PREDICTED before the run** (grid S × S with S = R + K + 2, R = ⌊(ℓ + 12)/24⌋ loud notes, K = ⌊ℓ/24⌋ first-step
forms; for the largest primes S = R + 3):
1. Level 1: exactly R unit factors, the rest divisible by ℓ. At 607, R − 1 = 24 (Boylan and Webb's missing note).
2. Level 2 (the even tower): exactly R unit factors. Its space modulo ℓ also holds the first-step systems, but they
   reach it only through U(ℓ) from level 1, so they arrive divisible by ℓ.
3. Level 3: R unit factors, then exactly K factors ℓ² at every ℓ ≥ 29 (every first-step family carries a share at depth
   one: the full share condition), none at 13, 17, 19 and 23 (K = 0), and the remaining factors divisible by ℓ³. On the
   smaller grid S = R + 3: three factors ℓ².
- Fewer than K factors ℓ² would mean some first-step family's share is missing or deeper; the prime is slow as soon as
  one is present.
- At 29 … 97 the result must agree with the measured rates (slow, the share at depth one).

**PREDICTED before the deeper run (derived from the depths at levels 1 to 3; `share_deeper.py`).** If the first-step
shares are slope-1 families, their depth is ⌊(b + 1)/2⌋ at level b: 1, 1, 2, then 2 at level 4 and 3 at level 5. A
slope-2 family has depth b. So on the same grids, modulo ℓ⁶:
- level 4: R unit factors, K factors ℓ², the rest divisible by ℓ⁴;
- level 5: R unit factors, K factors ℓ³, the rest divisible by ℓ⁵.

**Amended before 607's result came in (14 Nov, evening).** Prediction 1's "at 607, R − 1 = 24" took Boylan and Webb's
rank for ours. Theirs counts the span of the levels, one operator X; ours counts Hecke directions (the T(Q²)). At 607 the
two forms of the pair share their X-number modulo 607 (367) but not their T₅ numbers, which are conjugates in F_(607²)
(README "607's two forms, one digit deeper"). So this test should see all 25 loud notes at 607: **25 unit factors at
levels 1, 2 and 3**, not 24. The script's own check still uses 24 there; its verdict at 607 is read with this amendment.

**Counted: which families give the factors** (`share_identify.py`, output `share_identify_29-200.txt`, 502 s). The
Smith form's directions at level 3, read modulo ℓ on the grid's primes Q:
- at **every prime 29 … 199 (37 of 37)**, the unit directions are exactly the loud notes' Hecke numbers
  Q ↦ (3|Q)·a_Q(f), f of weight ℓ − 3;
- and the directions of valuation ≤ 2 are exactly those plus the first-step forms' numbers Q ↦ (3|Q)·Q·a_Q(g), g of
  weight ℓ − 5 (θ-twisted).
- So each factor ℓ² at level 3 is a first-step family, and each unit is a loud note. A chance match of two subspaces this
  size is out of the question (an (R + K)-dimensional subspace of F_ℓ^(R+K+2) landing on a given one).

**What Boylan and Webb's "R − 1 at 607" counts (14 Nov, evening; their paper read again, the authors' reprint; symbols
restored from the extracted text).**
- QUOTED, (1.13): r_ℓ(m) := rank over ℤ/ℓ^mℤ of Ω_odd(m), the span of the odd levels. §6.2 says what independence means:
  "a_1f_1 + ⋯ + a_nf_n ≡ 0 (mod ℓ^m) implies … a_i ≡ 0 (mod ℓ^m); if m ≥ 2, we further require that not all a_i ≡ 0
  (mod ℓ)". So only relations with a unit coefficient count.
- QUOTED, how "for all m" is reached: "We calculate r_1297(1) = 54. Hence, for all m ≥ 1, we have r_1297(m) = 54." That
  step works when r_ℓ(1) equals the bound R_ℓ. At 607, where r_607(1) = R − 1, it does not apply, and their modulo-ℓ²
  computations stopped at ℓ = 89. So "r_607(m) = R − 1 for all m" rests on m = 1.
- **The first digit, explained (DERIVED from our records):**
  - Two of 607's loud notes have the same two-level multiplier to one digit, 367 ("607's two forms").
  - They don't cancel. Stepping up the levels moves them together, so modulo 607 the levels see them as one direction:
    24 instead of 25. That is Boylan and Webb's R − 1.
  - Nothing is missing from the tower: the share test sees all 25 loud notes at full size, because T₅ tells the two
    apart where the levels can't.
- **The second digit (DERIVED; not computed at 607 directly, which would take about 50 levels):**
  - The two multipliers part at the second digit ("607's pair, one digit deeper in the tower itself").
  - Both shares are units (the share test).
  - Modulo 607², no relation with a unit coefficient joins 25 consecutive odd levels: it would need a double root at
    367 modulo 607 on top of the other 23 roots, 25 in all, from a polynomial of degree 24.
  - So by their own definition r_607(2) = 25. Their "R − 1 for all m" holds at m = 1 only.

**Counted: one and two levels deeper** (`share_deeper.py`, output `share_deeper_13-113.txt`, 7,323 s): **as predicted at
25 of 25 primes, 13 … 113.**
- Level 4: R units, K factors ℓ², and the rest ℓ⁴ or deeper.
- Level 5: R units, K factors ℓ³, and the rest ℓ⁵ or deeper.
- So the first-step shares go 1, 1, 2, 2, 3 digits deep at levels 1 … 5: one digit per two levels, the slope-1 rate, seen
  on partition numbers. At 13, 17, 19 and 23 every share below the loud note goes 1, 2, 3, 4, 5: one digit per level.

**Counted: the share test** (`share_test.py`; outputs `share_test_13-250.txt` 521 s, `share_test_251-400.txt` 13,673 s,
`share_test_only_409_499_607_613.txt` 3,578 s, `share_test_only_401_…_491.txt` 12,082 s, with their .jsonl; 12 or 4
workers at idle priority).
- **Every prime 13 … 397 on the full grid (73 of 73) as predicted:**
  - level 1: R units, the rest ℓ (or ℓ² for the last one at 13, 17, 23, 29, 37);
  - level 3: R units, then exactly K factors ℓ², then ℓ³ or ℓ⁴ (or ℓ⁵, which here means ≡ 0);
  - at 13, 17, 19 and 23 no factor ℓ²: the shares below the loud note sit 1, 2, 3 digits deep at levels 1, 2, 3.
- **Not predicted, counted at all 73:** level 2 has exactly K factors ℓ, then ℓ² or deeper. So the first-step shares sit
  1, 1, 2 digits deep at levels 1, 2, 3: the step from an odd level to the next even one leaves them as deep as they were,
  and the step from even to odd adds the digit.
- **The smaller grid (S = R + 3), at 401 … 491, 409, 499, 607, 613:** R units at levels 1, 2, 3, and all three free
  slots filled, by ℓ at level 2 and by ℓ² at level 3. So at least three first-step families carry a share at each.
- **607:** 25 units at levels 1, 2 and 3, as the amended prediction said. The script's own check, written with 24, flags
  607; with 25 units only three slots are free, and all three hold first-step shares.
- **What this settles:**
  - **The share condition holds at every prime measured.** On the full grid (29 … 397) every one of the ⌊ℓ/24⌋
    first-step families carries a share, exactly one digit deep at level 1.
  - **A4 ("what fades is the next family")** is no longer counted at two primes only. The part of the tower two digits
    deep at level 3 is exactly the first-step families, named by their Hecke numbers at 29 … 199.
  - **Slow, measured on partition numbers:** at every prime from 29 on that was tested, level 3 carries a part two
    digits deep, where a fast tower would carry none. At 29 … 113 it then fades one digit per two levels through level 5.
  - **What still rests on A1:** that this goes on at every level, i.e. that the first-step shares' multiplier per two
    levels is the families' U_ℓ number at weight −2, of valuation exactly 1 (Theorem 1.3 and the ghost series).

## Why depth one? Every family starts as a unit (14 Nov 2026, late)

**The reduction (DERIVED, with the steps measured above).** One level up is U(ℓ). From an odd level to the next even
one a first-step share keeps its depth; from even to odd it gains one digit. Depth one at level 1 therefore means: the
share of each first-step family in level 0, which is 1/η(24z) itself, is an ℓ-adic unit. That it is at least divisible by
ℓ at level 1 is proved (Yang: the odd tower modulo ℓ holds only loud notes). That it is not divisible by ℓ² is the
content.

**Chance alone (DERIVED from the first-step orbit sizes in `first_step_check_29-613.txt`).** If the share were a random
number modulo ℓ, an orbit of d conjugate forms would miss depth one with chance ℓ^(−d). Over 29 … 613 that predicts 0.65
misses (0.56 over 29 … 397), so seeing none has a chance of about 52% (57%). Testing more primes cannot separate a law from
luck: the expected count grows like log log ℓ. Testing more families can.

**Found by looking (COUNTED; a scratch comparison, not predicted before):**
- **Level 2 shows the class's slope spectrum.** For each system of the class (seeds Sym^r ⊗ det^m on the row, with
  c(r + 2; ε) newforms each, from `first_step_w610.jsonl`; slopes at −2 by slope_rule's formula, or its ghost count where
  the formula does not apply), the slopes sorted are, at every full-grid prime 13 … 199, exactly the depths at level 2:
  41 of 41. That is, every family seen starts as a unit.
  - 37 first looked like an exception on its 5 × 5 grid ([0, 0, 1, 2, 3] against the slopes 0, 0, 1, 2, 2). On 9 × 9 and
    12 × 12 grids, also with shifted rows and columns, both slope-2 families sit at depth 2: the small grid was degenerate.
- **The two half-steps split the slope.** At 37 on the 12 × 12 grid, levels 1 … 4 read [0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4],
  [0, 0, 1, 2, 2, 3, 4, 4, 5, …], [0, 0, 2, 3, 3, 5, 6, …], [0, 0, 2, 4, 4, 6, …]. A family of slope s starting as a unit
  sits at depth k·s at level 2k and k·s + ⌈s/2⌉ at level 2k + 1: the step from even to odd costs ⌈s/2⌉ digits, from odd to
  even ⌊s/2⌋. For the first-step families (s = 1) that is 1, 1, 2, 2, 3, the measured depths.
- **So depth one is a unit start plus a first step of ⌈1/2⌉ = 1.** The question becomes why 1/η(24z) starts every family as
  a unit.

**PREDICTED before the next two runs** (data or grids not yet looked at):
1. At every full-grid prime 211 … 397 (share_test's saved level-2 depths), level 2 equals the slope spectrum as above.
2. On 12 × 12 grids at 13 … 61 (and 16 × 16 at 53 … 61), levels 1 … 4 are exactly the spectrum with unit starts and the
   half-steps ⌈s/2⌉, ⌊s/2⌋ (depths capped at 6, the precision). A miss, if confirmed on a shifted grid, is a family that
   starts deeper, the first one found.
1b. At the same primes, the saved levels 1 and 3 follow the half-step rule from the same spectrum (depth ⌈s/2⌉ at level 1,
    s + ⌈s/2⌉ at level 3), capped at 5, share_test's precision.

**Counted: prediction 1 and 1b** (`share_spectrum.py --saved`, output `share_spectrum_saved_211-397.txt`, 7 s): **as
predicted at 32 of 32 primes 211 … 397, levels 1, 2 and 3 all exact.** Every family the grids see starts as a unit, and
each single step costs ⌈s/2⌉ (even to odd) or ⌊s/2⌋ (odd to even) digits.
- The same comparison at 13 … 199, run before (unsaved, not a test): 38 of 41 exact. The three misses, at 23, 37 and
  103, are all in the last slot of a grid, one digit too deep; at 37 bigger and shifted grids remove it (above). 23 and 103
  are rerun on 14 × 14 shifted grids.

**Counted: prediction 2, big grids** (`share_spectrum.py --grid`; outputs `share_spectrum_grid12_13-61.txt` 62 s,
`share_spectrum_grid16_53-61.txt` 104 s, `share_spectrum_grid14s3_only_23_103.txt` 207 s,
`share_spectrum_grid16s5_only_59.txt` 101 s).
- **Exact at all four levels:** 23, 29, 31, 37, 41, 43, 47, 53 and 61 on 12 × 12 grids; 53 and 61 on 16 × 16; 103 on a
  shifted 14 × 14 grid (so its saved-grid miss was the grid).
- **13, 17 and 19: level 1 only.** These primes have few systems (4, 7 and 9), so the extra slots hold each system's
  second family. That family has the same Hecke numbers modulo ℓ as the first, so the grid sees it only through
  higher-order differences, deeper than its own depth (as with 607's pair). The prediction did not model this. Levels 2,
  3 and 4 are exact.
- **23 on the shifted grid: one slot at level 1 one digit deeper.** It is a small prime with second families in play; not
  resolved.
- **59: a family one digit deeper, the same on a shifted grid.** Named by `share_identify_seed.py` (outputs
  `share_identify_seed_59_grid16_s0_level2_prec11.txt`, `…_level4_prec11.txt`, modulo 59¹¹):
  - the system of the seed Sym⁴⁶ ⊗ det⁴ with T₅ ≡ 18 (one of the two rational newforms of weight 48) sits at depth 5 at
    level 2 and 9 at level 4, against 4 and 8 for a unit start;
  - so its slope is 4, as the ghost series says, and **its share starts one digit deep**;
  - each of the other systems at 59 (both loud notes, both first-step families, the slope-2, slope-3 and slope-5
    systems, the other two slope-4 systems) sits exactly where a unit start puts it, at levels 2 and 4.
- **So unit starts are the rule, not a law.** 59 is the first exception found, for a slope-4 family. The first-step
  families' depth one at every prime to 613 is the generic case: chance predicts about 0.65 first-step exceptions up to
  613, and none was seen.
- **Slow is unaffected.** It needs a nonzero share, not depth one. A first-step family starting a digit deeper would only
  move the start of the slow rate one level up.

**Counted: the rest of the range** (`share_test_only_503_…_557.txt` 14,520 s, `share_test_only_563_569_571.txt` 8,923 s,
`share_test_only_577_587_593_599_601.txt` 6,132 s; smaller grid): as predicted at 15 of 15. **With these, every prime
29 … 613 (103 of 103) carries a first-step share at depth one: level 3 has a factor ℓ² at each, so each is slow,
measured on partition numbers. 13, 17, 19 and 23 have none.**
- Update to the label table in "A further review of this proof": the row "Every ℓ ≥ 29 is slow" now reads: measured at
  every prime 29 … 613 at level 3 (and at 29 … 113 through level 5), and at every level given A1. The share condition is
  no longer assumed for these primes. Beyond 613 it is a measured rule, not a proof.
