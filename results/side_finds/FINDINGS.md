# Findings register (results/side_finds, 11–13 Nov 2026)

One line per finding, with its status and where its record is. Each section of README.md names its scripts and kept
outputs (.py beside .txt/.json). README.md's opening "New? (strict)" was written on 11 Nov, before the digs below.
STATEMENT.md states the main claims precisely (A–H), and PROOF.md has the proof route. What is still open (finishing the
proof, and the places to dig further) is in NEXT.md.

**Status words.**
- **new**: not found in the sources we read. Our search was limited, so "new as far as we read".
- **rediscovered**: found from the data, then found in print.
- **known, sharpened**: in print; we add a measurement.
- **reading**: our derivation, not a proof.
- **proved here**: proved by us; for an expert to check.
- **corrected**: said wrong first, then fixed in the record.
- **open**: unexplained.

## New on 14 Nov 2026, in plain words

23. **The tower carries its slow quiet notes, at every prime to 613 (a test on partition numbers alone).**
   - **The question:** the proof that every prime from 29 on is slow needed the tower to actually contain its slowest
     quiet notes, the first-step families. Nothing published says it does, and the slow rate could not be proved
     without it.
   - **What we did:** a quiet note the tower contains shows up at level 3 exactly two digits deep; nothing else does. A
     table of partition numbers, with rows for square classes and columns for primes Q, sorts the notes by that depth. No
     forms are used.
   - **What we found:**
     - At every prime 29 … 613 (103 of 103) the slow quiet notes are there, and from 29 to 397 every one of them is.
     - At 13, 17, 19 and 23 there are none, as fast primes need.
     - Read against the forms (29 … 199), the notes two digits deep are exactly the first-step forms, and the full-size
       ones exactly the loud notes.
     - At 29 … 113 the slow notes then fade one digit per two levels through level 5, seen directly.
   - **Why "depth one" (a second finding):**
     - Every family the tower carries starts as a unit. Each step up one level then costs a family of slope s ⌈s/2⌉ or
       ⌊s/2⌋ digits.
     - So the depths at level 2 are exactly the slopes of the class's families. This was predicted and then counted at
       211 … 397 (32 of 32), and counted at 13 … 199 and on big grids at 13 … 61.
     - Depth one for the slow notes is that rule with s = 1.
   - **The first exception:** at 59 one slope-4 family (seed Sym⁴⁶ ⊗ det⁴, T₅ ≡ 18) starts one digit deep. It was
     confirmed on a shifted grid and modulo 59¹¹. So unit starts are the rule, not a law. For the slow notes, chance
     would predict about 0.65 exceptions up to 613; none was seen. Slowness needs only a nonzero share, so it is
     unaffected.
   - **On the side:** Boylan and Webb's "one note missing" at 607 is two loud notes with the same multiplier to one
     digit. The levels cannot tell them apart modulo 607, but the tower holds all 25 at full size.
   - **Status:** counted. It removes the share condition and A4 for 29 … 613. What remains for "slow at every level" is
     A1: that each quiet note's multiplier is its family's number at weight −2.
   - **Records:** README "The share condition, tested on partition numbers" and "Why depth one? Every family starts as a
     unit"; `share_test.py`, `share_identify.py`, `share_deeper.py`, `share_spectrum.py`, `share_identify_seed.py`.

22. **Why exactly 13, 17, 19 and 23 settle fast (a derivation, then a test at every prime to 613).**
   - **The question:** why do the towers of 13, 17, 19 and 23 settle one digit per level, and every prime from 29 on
     only one digit per two levels? It was open in our record from the start (NEXT 5).
   - **What we found:**
     - When the weight goes up one step, the class picks up exactly the newforms of level 6 and weight ℓ − 5 with one
       fixed pair of signs at 2 and 3. This comes from how Sym^(2ℓ−6) breaks into pieces.
     - Those forms arrive at a rate of one per 24 units of weight, so there are ⌊ℓ/24⌋ of them.
     - They, and only they, give families of slope 1 at the tower's weight. A family's slope is min(m, ℓ − 4 − m) for a
       seed twisted by ω^m.
     - A slope-1 family next to the loud note (slope 0) is what makes a tower slow.
     - So the tower is fast exactly when ⌊ℓ/24⌋ = 0: at 13, 17, 19 and 23 and nowhere after.
   - **The test:**
     - the first-step gain equals ⌊ℓ/24⌋ at all 108 primes from 11 to 613 (607: 25);
     - the reduction matches the direct count at 59 of 59 primes;
     - the slope rule holds for every possible seed at every prime 13 … 200, and for 23 of 23 families counted;
     - the newform counts step by exactly 1 per 24 in weight, so the formula carries to every prime.
   - **Status:**
     - The count ⌊ℓ/24⌋ was the review's. The reduction to newforms of weight ℓ − 5, its proof route and the slope rule
       are ours.
     - It rests on the Serre-weight bookkeeping (standard, checked), the trace formula's periodic shape (standard), and
       the ghost slopes: proved for the slope-1 families (very generic), counted directly for the rest at 13 … 23.
   - **Since then (a further review, assessed in README):**
     - Every first-step family at every prime 29 … 613 is ordinary and irreducible, with a proof for each. So the proved
       theorem gives them slope exactly 1 at the tower's weight; it holds on the whole weight disc, as quoted.
     - At 607 those families are the two quiet braids, 8 and 17, found among the forms. 607's partition numbers reach
       only level 5, so whether its tower carries them is not measured.
     - The result with its conditions: **every ℓ from 29 to 613 is slow, given A1, A4 and that the slope-1 families
       carry a share of the tower; 13, 17, 19 and 23 are fast, measured.** The share condition is the one part we cannot
       prove; it held wherever the rate was measured.
     - The identity dim S_(ℓ−1)(SL₂(ℤ)) = R_ℓ + ⌊ℓ/24⌋ is Hermite's identity; the numbers agree, but they do not show a
       Hecke-structural split (reworded at a relayed review).
   - **Records:** README "Why exactly 13, 17, 19 and 23 settle fast".

21. **Which families the proved theorem covers (a reading of the papers, then measurements).**
   - **The question:** Liu, Truong, Xiao and Zhao proved the ghost conjecture under conditions. Do our families meet them?
   - **The conditions** (read in the papers, QUOTED in README): ℓ ≥ 11; the mod-ℓ representation irreducible, but
     reducible locally at ℓ; and its Serre weight "very generic", 2 ≤ a ≤ ℓ − 5.
   - **Measured:**
     - every family's seed is ordinary (39 of 39), so every family is reducible locally at ℓ, as required;
     - 8 families are truly Eisenstein mod ℓ (tested on 23 primes), all of the form ω^(j+1) ⊕ ω^j.
   - **So:**
     - 21 of the 39 families are covered: very generic, irreducible, ℓ ≥ 11.
     - 18 are not: 8 Eisenstein, 8 of the Steinberg weight (seed of weight ℓ + 1, covered only by a 2026 local
       theorem), and 2 with a = ℓ − 3.
     - Among the partition tower's own forms, 11's law form and the quiet forms of 13 and 17 are covered. 5's and 7's
       forms and 11's quiet form are not.
   - **They all fit the ghost series anyway** (item 20). The uncovered ones are evidence for the cases the authors
     "expect" (for a = ℓ − 3 they are also "not sure") and for Bergdall and Pollack's corrected conjecture, which
     includes the Eisenstein-type ones.
   - **Since then (three relayed reviews, assessed in README):**
     - The 8 Eisenstein families are now **proved** Eisenstein, by Sturm's bound (θ^(ℓ−1)·seed ≡ θ^(ℓ−2)·E_(ℓ+1) mod ℓ
       in every coefficient up to ℓ² + ℓ − 2).
     - The "projectivity" caveat is gone. Within one ρ̄ the signs at 2 and 3 are fixed: no ρ̄ was found in two sign
       classes at any of the 7 primes, and the two apparent cases at 19 were a quadratic twist. So the class's slopes
       are the whole family's slopes minus the old forms', and Theorem 1.3 covers each of those.
     - For the 21 covered families the class's slopes follow from the proved theorem itself.
   - **Status:** the hypotheses and the subtraction argument are the papers' and a review's; the measurements, the
     Eisenstein proof and the sorting are ours.
   - **Records:** README "Do the families meet the conditions of the proofs?".

20. **Every slope of the class comes from a few small forms (a test of the ghost series).**
   - **What we did:** from each system's seed alone (item 19), we computed its counts at both levels at every weight.
     From those we built Bergdall and Pollack's full ghost series and read off its slopes.
   - **The test:** every slope of every class form at level 6, second and third forms of a system included. Predicted,
     then counted at 7 primes (5 to 23), up to weight 174: 225 of 225 (system, weight) pairs agree.
   - **It also settles the one case the waiting rule did not cover:** the second form of 11's law system, whose slope
     went 4, 3, 5. The ghost series gives it exactly.
   - **Status: known (the ghost conjecture), tested here.** The review says it is proved under a genericity condition.
     Our Eisenstein-type systems may lie outside that, and they fit anyway.
   - **What it means for the towers:** the quiet notes' slopes, and so the settling rates (item 16), are fixed by a few
     forms of weight at most ℓ + 1 at level 6.
   - **Records:** README "The full ghost series from the seeds, against every slope".

19. **Why level 6ℓ gains one or two forms per step (a derivation, then a test).**
   - **The question:** at level 6ℓ each system gains forms on a straight line, one or two per step (item 18). Why
     straight, and why 1 or 2?
   - **What we found:** every system has one **seed**, a single small form at level 6 of weight between 4 and ℓ + 1,
     twisted (its Serre weight).
     - Each weight step adds to level 6ℓ the same complete set of ℓ − 1 building blocks (principal series of GL2(F_ℓ)),
       hence the straight line.
     - A seed of weight 4 … ℓ sits in two of those blocks: 2 forms per step.
     - A seed of weight ℓ + 1 sits in one: 1 per step.
     - A system of Eisenstein type (its numbers ≡ q^a + q^(a+1) mod ℓ) with a seed of weight ℓ − 1 also gains 1: one
       block's constant piece takes a form back each step.
     - The same seeds, summed over Sym^(k−2) instead, give level 6 itself, including the step at which each form
       arrives. So the waiting's start and its end both come from the seeds.
   - **The test:** predicted, then counted, for 39 systems at 5, 7, 11, 13, 17, 19 and 23, every count at every weight,
     both levels: 132 of 132. Counted again including every congruent form (old at 2 or 3): 88 of 88.
   - **Status: known in principle.** It is the standard machinery (Shapiro's lemma, Serre weights, Jacquet–Langlands to
     the Shimura curve of discriminant 6). It is the kind of dimension count the ghost conjecture takes as input.
     Bergdall and Pollack compute their dimensions the same way: from weights 2 … p + 1 and twists, "plus Eisenstein
     multiplicities" (QUOTED, item 21). Ours are the seeds of the partition class, the Eisenstein-type correction as it
     shows there, and the test.
   - **One correction was found on the way, not assumed.** 7's quiet form, and 5's law form, are of Eisenstein type
     mod ℓ. For them the plain count is too big by one per step until the constant factors are taken out (an Euler
     characteristic argument). That was derived after the first trial at 7 and then held at every other prime.
   - **Records:** README "Why level 6ℓ gains one or two forms per step".

18. **Why a form waits where it does: the forms of level 6ℓ (a reading, then a test).**
   - **The question:** a form of slope s waits s steps before it first appears (item 17). What starts the waiting? The
     loud notes do not.
   - **What we found:**
     - The form's system is already there one level down, among the forms of level 6ℓ, before it reaches level 6.
     - There it gains one or two forms per step, on a straight line.
     - The waiting starts at the first step where level 6ℓ holds two of them. Each step of waiting adds one to the slope.
   - **Status: rediscovered.** This is Bergdall and Pollack's ghost series (first coefficient). A relayed review
     confirmed it and cites "Slopes of modular forms and the ghost conjecture" (IMRN 2017) and part II (QUOTED, not read
     here). The rule is theirs; ours are the test on the partition class and the link to the settling rate (item 16).
   - **The test:**
     - A new tool builds all the forms of level 6ℓ over the ℓ-adic integers and counts them by system
       (`level6l_counts.py`).
     - Predicted, then counted at 29 weights for 5, 7, 11 and 13. All 19 predicted steps came out as predicted.
     - The slopes of every system at every later weight, read from level 6ℓ alone, match the counts: 30 of 30, including
       the steps up at the congruent weights.
     - Four forms still waiting at the end of those runs had their slopes predicted before they reached level 6:
       8, 8, 4 and 10. All four are right.
   - **Checks:** the dimension formulas at every weight (58 of 58 sums), and the operators right on every coefficient.
   - **Open:** why level 6ℓ grows by one or two per step, and when a form reaches level 6 (the end of its waiting).
   - **Records:** README "Why the waiting starts where it does: the forms of level 6ℓ".

17. **When a form first appears tells its slope (a reading, then a test).**
   - **The rule:** a form of slope s waits s steps along the weights before it first appears. From then on its slope is s,
     except at weights whose step is congruent to one of those s steps mod ℓ: there it is s + 1.
   - **The test:** predicted, then counted at 37 weights, for 8 forms at 5, 7, 11, 13, 17 and 29. All 37 came out as
     predicted. Ordinary forms (slope 0) never step up.
   - **What it explains:** from 13 on, the waiting starts at the very first weight. So a quiet form's slope is the step at
     which it first appears, which is the fast/slow rule (item 16).
   - **Open:** why the waiting starts where it does. The likely answer is the forms of level 6ℓ (the "ghost" picture of
     Bergdall and Pollack). **Answered (item 18):** the waiting starts where level 6ℓ first holds two forms of the
     system.
   - **Records:** README "When a form first appears, and its slope".

16. **What makes a tower settle fast (a reading, then a test).**
   - **The question:** 13, 17, 19 and 23 settle one digit per level; 29 and every prime after it, one digit per two levels.
     Why? (Open in our record; for 13 the fast rate is Atkin and O'Brien's theorem, 1967.)
   - **What we found:** follow the class's forms as the weight steps up from ℓ − 3. If the class gains no form at the first
     step, its next forms come one step later, fade by ℓ² every two levels, and the tower is fast. If it gains forms at
     once, those fade by only ℓ, and the tower is slow.
   - **Why the line falls between 23 and 29:** the class gains about (ℓ − 1)/24 forms per step on average: fewer than one
     below 25, more than one above. The gains alternate from step to step (the class's signs flip each step), so the rule
     itself is the first step's own count, not the average (a relayed review's correction, accepted).
   - **The test:** predicted, then counted at 19, 23, 37, 41, 43, 47, 53, 59 and 61. It also fits 13, 17, 29, 31 and 607,
     where the first step adds 25 forms, all of slope 1: 607's two quiet braids. At 13 to 31 the number of new forms is the
     number of quiet notes measured from partition numbers.
   - **Checked independently:** a review recounted the forms with PARI, from 13 to 61, and its counts equal ours.
   - **The rate itself comes from the quiet notes everywhere we looked:** it is the quiet note's slope minus the loud note's,
     the quiet note sitting in a system of its own. This also settles an old puzzle at 7 and 11. Their rates, 3 and 1,
     looked like half a gap only because the next slope was taken in the law's own system. Their quiet forms, in other
     systems, have slopes 4 and 3, which gives exactly 3 and 1.
   - **And those quiet notes are confirmed digit by digit:**
     - Carried to the tower's weight, 7's quiet form gives 6, 2, 3, 1 and 11's gives 4, 5, 0, 3, 1.
     - The fading parts measured on exact partition numbers are 6, 2, 3, 1, 4, 1, 3 and 4, 5, 0.
     - They agree in every digit both carry.
     - **Verified independently** by another session, from the definitions alone, with its own code (PARI for the forms):
       - the same slopes and carried values;
       - and, with a method that needs no loud constant, the partition numbers agree in 4 digits at both primes;
       - that method gives 7's quiet note to 16 digits from partition numbers alone.
   - **Open:** a proof that the first step gains a form at every ℓ ≥ 29 (the exact dimension formula). Also open is
     whether the ghost theorem of Liu, Truong, Xiao and Zhao covers these systems, which would derive the slopes.
     **Answered since (items 21, 22):** the first step gains exactly ⌊ℓ/24⌋ forms (≥ 1 from 29 on), and the theorem
     covers the slope-1 families that make the towers slow.
   - **Records:** README "What sets the settling rate".

## New on 13 Nov 2026, in plain words

1. **Divisibility laws for 73 (and 61).**
   - **What they are:** Ramanujan found that p(5n + 4) is always a multiple of 5. Laws like it exist for every prime,
     but beyond 31 they are rare and enormous, and none had been listed for 73.
   - **What we found:** 25 below Q = 50 million. For example, p(107239327321929851108449·n + 9040523964878168) is a
     multiple of 73 for every n.
   - **The check:** partition numbers about 100 million digits long are multiples of 73. The control numbers next to
     them land on exactly the value the theory gives.
   - **61:** the same, with 35 laws.
   - **Records:** README "Lucky primes for 73"; STATEMENT E3, E4.
2. **A tower is a chord, and it can be taken apart.**
   - **The tower:** zooming into p(n) by a prime ℓ again and again gives a tower, and its law mixes several numbers,
     the constants.
   - **The chord:** the tower is a sum of simple towers, one per modular form of ℓ's class, like notes in a chord.
   - **The filter:** one extra step at the prime 5 works like a filter. It silences chosen forms, and what is left is
     one simple tower, multiplied by its own constant every two levels.
   - **Records:** README "97's tower is two towers", "The twins taken apart"; STATEMENT B2.
3. **So each constant can be read off partition numbers, digit by digit.**
   - **97:** 28, then 26. The whole tower would have needed numbers a billion digits long.
   - **The twins:** four digits for each twin at 37 and 41, three at 43 and 53.
   - Two levels deeper show one more digit, and a twin never shows the other twin's digits (147 of 147).
4. **Each simple part has its own lucky primes, far more of them.** 97's simple part has laws at about 3 primes in 97
   (the first is 53). The whole tower has them at about 1 prime in 30 million.
15. **The quiet notes can be read from small weights (13 Nov).**
   - **The old way:** a quiet note's family was followed only through weights ℓ(ℓ − 1) apart, where it is guaranteed to
     be smooth. At 607 the first such weight is 367,840.
   - **Small weights work:** weights only ℓ − 1 apart work too. At 13, 17, 29 and 31 the quiet note's own family shows up
     at every small weight, and carried back to the tower's weight it gives the right digits, about one per weight.
   - **At 29, a third digit:** 23, 3, 4, 28, 20.
   - **Confirmed on partition numbers.** Up to p(3.4·10¹⁶), the tower's own quiet note at 29 is 23, 3, 4, settled by all
     n. The shortcut without α settles 23, 3, one digit fewer at the same depth.
   - **This brings 607's quiet notes within reach,** at weights 1,210 to 2,422.
   - **607's quiet layer:** 25 quiet forms, as many as loud ones, in two braids of 8 and 17, with no lone quiet note
     (counted at weights 1,210, 1,816 and 2,422).
   - **607's quiet braids at the tower's weight, to three digits (14 Nov):**
     - Each strand's first digit grows in step with the weight (counted at all three weights, 25 of 25 coefficients).
     - So a braid is carried to the tower's weight with each coefficient divided by its power of the step.
     - The first fit (fading_forms4) had not done that, and its answer broke the braids apart.
     - Both braids' polynomials are known to 3 digits (2 confirmed, the 3rd likely), and every strand is refined to 3 digits.
     - The strands multiply to 458, 397, 109 (braid of 8) and 292, 480, 577 (braid of 17), digits lowest first.
     - Records: README "607's quiet braids"; fading_forms5_607_s1_fit.txt, fading_forms5_607_strands.txt.
   - **Records:** README "Toward 607's quiet notes".
14. **The step at 5 is a prism, and 607 is a colour the tower cannot see (13 Nov).**
   - **At 37, a prism.** Hearing each tower value together with the one 25 times further along splits the tower into its
     two notes. Each share is then multiplied by its own constant (60 of 60 and 50 of 50), and the shares add back up to
     the tower (79 of 79).
   - **At 607, a colour the tower cannot see.** The same step at 5 tells all 25 of 607's forms apart (counted on the
     forms). Down the tower, on the first digit, two of them are multiplied by the same number (240) and sound as one.
     That is the exception: 24 voices, not 25. One digit deeper they part.
   - **Records:** README "The prism at 5" and "607 through the prism"; the page "The Prism at 5".
13. **Notes in a line (13 Nov).**
   - **Across the rows, no pattern beyond chance.** Laid out as rows from 1 back to 1, the notes were tested on 150 shifts
     ("up one, over three") against 2,000 shuffles. The only real patterns are stretches inside one prime, where one
     note is a power of the other.
   - **Down the tower, one note is one straight line.** Placed by its logarithm, a one-note tower's value at every level
     follows from its start alone. Every nonzero value at 13 to 31 sits on its line (465 of 465), and the zeros stay
     zero.
   - **Two notes need two starts:** at 37, the two-step rule held 41 of 41 times.
   - **Records:** README "Notes in a line"; the page "Notes in a Line".
12. **The quiet notes at 19 and 23 come in twos, and all four are found.**
   - At 19 and 23, the leftover that fades is not one note but exactly two, each shrinking by ℓ² every two levels. One
     two-step rule fits every n, to 6 to 8 digits (counted).
   - At 19 and 23 each of the two is a steeper (slope-2) family of modular forms, from two different families on
     lane ℓ. Computed from the forms alone, their constants match the ones measured from partition numbers in every
     digit there is (4 of 4 each): 19 → 2, 1, 0, 17 and 13, 14, 2, 15; 23 → 8, 20, 21, 8 and 17, 20, 19, 21.
   - Along the way, the "crowds" that seemed to hide these forms at 13, 17, 29 and 31 turned out to be too few digits,
     not crowds (see Corrected).
   - **Records:** README "What fades at 19 and 23: two quiet notes at once".
11. **A braid in size and direction.**
    - **The picture:** like a complex number, each of a braid's numbers has a size (an ordinary number mod ℓ, the
      product of its strands' values) and a direction.
    - **Two levels down:** the size is multiplied by the product of the braid's constants (44 at 73) and the direction
      turns by one tick. That held for all 240 real values at 73.
    - **So:** a braid is an ordinary note in its size and a pure turning in its direction.
    - **Records:** README "The braids in size and direction".
10. **What a braid's angles obey.**
    - **The ratio:** each strand turns exactly ℓ times as far as the one before.
    - **The sum:** all the strands' turns add up to the angle of the product of the braid's constants, an ordinary
      number mod ℓ, so the sum always lands on a coarse grid (5° steps at 73).
    - **Odd or even:** the sum is an odd number of steps exactly when that product is not a square. At 73 it is odd
      (44 is not a square), the same fact that rules out a one-level law there.
    - **The twins** are powers of each other: 34 = 15²⁰ (mod 37).
    - **Records:** README "The braids' angles".
9. **Each note's full form is a star.**
   - **The star:** follow one note down two levels at a time and it turns by the same angle again and again, tracing a
     star.
   - **Its points:** how many steps the note needs to come back to where it started.
   - **Corrected:** how the points are joined, and so the star's look, depends on the number we walk the circle by,
     not on the note.
   - **The whole tower:** it comes back when all its notes do. That matches every tower loop measured before (13 of 13).
   - **New loops:** 389,016 steps at 73 and 304,224 at 97.
   - **Records:** README "Pictures of the chords"; the page "Towers as Chords".
8. **How fast a tower reaches its constants.**
   - **Fast:** at 13, 17, 19 and 23 each level adds one digit of the constant. 13's speed is Atkin and O'Brien's
     congruence (their theorem, 1967; an earlier line here said "proved in 2025", corrected 14 Nov).
   - **Slow:** at 29, 31 and every note from 37 to 97 it takes two levels per digit, the rate Boylan and Webb proved
     possible.
   - **The cause:** a first explanation (extra non-ordinary forms in weight ℓ − 1) was tested and is wrong; why is
     open.
   - **What fades is a quieter second note:** it shrinks by ℓ² every two levels at the fast primes and by ℓ at the slow
     ones. At 13 and 17 (and 29, 31) it is one note with a constant of its own, now read off partition numbers
     (13: 169 × (10, 11, 6, 2, …)). At 19 and 23 it is itself a little chord.
   - **What the fading note is (found):** at 13 and 17 it is the weight −2 value of a steeper (slope-2) family of
     modular forms with the class's signs, from another family on lane ℓ. Its constant, computed from the forms
     alone, matches the one measured from partition numbers in all 6 digits at both primes (13: 10, 11, 6, 2, 10, 10;
     17: 3, 9, 1, 13, 4, 3). At 29 and 31 (slope 1) it is found too, to the two digits two weights give (29: 23, 3;
     31: 2, 28); the third weight was stopped by low memory.
   - **Records:** README "The speed of the towers"; NEXT item 5.
5. **Two tools.**
   - **(13 Nov) The quiet notes without the loud note's constant.**
     - Take three levels and form H(b+4)·H(b) − H(b+2)². The loud note cancels exactly, and what is left shrinks by α·β
       every two levels.
     - At 13 and 17 this reads α·β from partition numbers alone, in every digit all n share (4 at 13, 5 at 17).
     - Bigger versions of the same combination take chords apart (README "The quiet notes without the loud note's
       constant").
   - A lighter test for lucky primes: the same answers as before, 39 of 39, about 20 times faster.
   - The exact control rule: at a lucky prime, the numbers next to the zeros take one predicted value (36 of 36).
6. **More divisibility laws: 67, 71, 79 and 83.**
   - **Found:** 92 laws, just as many as predicted (89). Each checked law held on exact partition numbers (6 laws,
     12 of 12 members, every control exact).
   - **The simplest:** p(19728144816847346813·n + 6177904563614) is always a multiple of 67.
   - **Records:** README "Laws at the other primes with three forms".
7. **Even braided notes can be heard alone.**
   - **Braided forms:** some forms come in braided groups (a pair, or a triple at 73 and 97). Their constants are not
     ordinary numbers but numbers of a bigger number system, as √2 is for whole numbers. A filter made of ordinary
     numbers silences a whole braid or nothing of it.
   - **The way through:** a filter made of numbers from that bigger system pulls the strands apart. Then each braided
     constant shows alone in the partition numbers, digit by digit.
   - **The check:** 108 of 108 cases, and a strand never shows a sister strand's digits (79 of 79).
   - **So:** every form of every class from 37 to 97 can now be seen on its own.
   - **Records:** README "The braided forms taken apart"; STATEMENT B3.


## The general law (STATEMENT A–H; README "How wide the recipe reaches")

| finding | status | record |
|---|---|---|
| A tower of any ℓ-adic modular form obeys a short law. The constants are the U_ℓ-values, at the form's weight, of the ordinary families it uses: Eisenstein family 1, cusp families A_f(w), twist ψ(ℓ) after the Shimura step | general shape known (Hida, FKO, Boylan–Webb); the explicit recipe new | STATEMENT A, G, H1; PROOF.md |
| Where no ordinary family is used, the next slope takes over: the laws of Ramanujan, Watson, Atkin, Lehner and Garvan are this case | reading | Dig 1, 14th and 21st; spt; j |
| At weight 0 the Eisenstein member is the constant 1, so for n ≥ 1 the cusp forms alone set j's constants | reading, checked | j |
| **A simple pole at a CM point whose curve is ordinary at ℓ adds the constant u^(w−2)** (u the unit root of that curve's Frobenius). 1/E₄: u⁻⁶. 1/E₆: u⁻⁸, predicted before its run | new; the derivation is a reading | 1/E₄, 1/E₆; STATEMENT H3 |
| **A pole at a point with no complex multiplication (j = 2) adds u_C⁻¹⁴, u_C the unit root of the curve at the pole** (y² = x³ + 10356x + 11916304). Predicted before the run. Met one ℓ-adic digit per level, since the pole moves to the canonical lift. Holds at 5, 7, 11, 13, 17, 19, 23; P alone at 29 and 31 (supersingular). **Blind: the tower alone finds j = 2** (the only fits with |j| ≤ 10⁶ are 2 and −575358, a chance fit that 37–53 rule out) | new; the derivation is a reading | A pole at a point with no symmetry |
| **Curves that repeat: a pole ON a CM curve (the canonical lift) gives its constant exactly, 2 digits per level; the pole at j = 2 drifts toward that lift and gains 1 digit per level** (5: lift j = −32768; 7, 19: 54000; 11: −3375). **Two CM poles of discriminant −15 swap at 17 and 23: the tower repeats every two levels (two-level ratio = U⁻²⁶, only the pair law x² − U⁻²⁶ holds); at 19 and 31 each is fixed (u⁻²⁶)**. **The three CM curves of discriminant −23 cycle with period 3 at 13: the remainders' ratio settles only three levels apart (U⁻³⁸, at two levels), and only the cubic law x³ − U⁻³⁸ fits (it has no root in ℤ₁₃)**. **At small primes the 3-cycle is confirmed in full (discriminant −31): at 5 the ratio settles three levels apart at all six levels and P·(x³ − U⁻³⁸)·S reaches the whole lane 5¹¹; at 7 the cubic law reaches 7¹⁰, two digits a level**. The cycle length, the order of ℓ's prime in the class group, is read from the coefficients. All predicted before the runs | new; the reasoning is a reading | Curves that repeat themselves |
| **The tower measures how far a curve is from its returning twin.** At 5, poles at j₀ ≡ 2 (mod 5) that start s = v₅(j₀ + 32768) digits from the twin j = −32768 give depth min(2b + 3, b + s, 11) at level b. Predicted before the runs; 48 of 48 levels exact for s = 1, 2, 3, 4, 6 and the twin itself | new; the Serre–Tate reasoning is a reading | What the tower view reveals |
| **Curves that never come back, predicted and seen.** At 5 (weight −48): poles at the four curves of discriminant −275 (quasi-canonical lifts: CM, conductor divisible by 5) never return to their start but land on the returning curve −32768 after one step, and from level 2 the tower is indistinguishable from a pole on that curve (depth 6, 10, 11, …). Four curves without CM gain one digit a level forever (4, 5, …, 10). Predicted before the runs | new; the reasoning (Serre–Tate, Gross) is a reading | Curves that never come back |
| The boundary: the data-only search finds laws in every modular function tried, none in partitions into squares or cubes | new (as a test) | How wide …; STATEMENT H4 |
| Non-congruence forms obey laws of the same shape with Frobenius constants (Scholl 1985) | known, placed here | How wide … |

## The partition function p(n)

| finding | status | record |
|---|---|---|
| Ramanujan's three laws mod 125 found blind | rediscovered (his manuscript, eq. 21.3) | Dig 1 |
| Mod 125: p(N) is a multiple exactly when (24N − 1)/25 is a square on lane 5 (and the analogue on lane 7) | reading, counted | Dig 1 |
| p(2401n + B) ≡ 0 (mod 7³) for B = 929, 1615, 1958; p(78125n + B) ≡ 0 (mod 5⁷) for B = 46,224, 77,474; p(117649n + B) ≡ 0 (mod 7⁴) for B = 45,519, 79,133, 95,940 | new as stated; they follow from Watson's and Garvan's theory | Dig 1 and Dig 1, deeper |
| Ramanujan's wrong 7³ guess holds exactly on the square classes; OEIS A340957 and A340757 have the numbers but not the law | new as stated | Dig 1, deeper |
| Ramanujan's tally for 7 (n ≤ 200) has slips: 51, 34, 23, 27, 23, 23, 19 | counted correction | Dig 1, again |
| The extra multiples of 5, 7, 13 come from Atkin's 1968 law. A formula μ ≡ Q⁻²[(3\|Q)(1 + Q) − (−3\|Q)(m\|Q)], 0 wrong of 187,989 | law rediscovered; the formula not seen written | Dig 1, again |
| One more factor where μ = 0, e.g. p(1375n + 149), p(1375n + 424) ≡ 0 (mod 121) | rediscovered (Lovejoy–Ono; Jameson Cor. 1.3) | Dig 1, third |
| Two levels apart, the laws divided out, the constant settles: 5 → 4, 2, 3, 1, 0, … (not a simple fraction), 7 → 5, 6, 1, 1, 2, …, 11 → 9, 3, … | new as stated | Dig 1, third |
| p(25N − 1) ≡ 25·N·(N − 3)·p(N) (mod 125) for every N, at every depth | new as stated; follows from Ramanujan's eq. 21.4 | Dig 1, fourth |
| A lane-49 form of the step rule for 7 | new as stated | Dig 1, fourth |
| Curves behind the towers: 7 → 42a, 11 → 726d (Jameson's), **13 → 1014f (Q¹⁰) and its twist 1014b (Q⁴)** | 7, 11 rediscovered; 13 new | Dig 1, fifth and sixth |
| The 13 statement proved with Sturm bounds, Jameson's method, at every level | proved here | Dig 1, seventh and eighth |
| No elliptic curve can do this for 17, 19, 23, 29, 31: at some prime Q the required number of points breaks Hasse's bound | proved here | Dig 1, eighth |
| Abelian varieties take over: at 31, two of dimension 4 at level 5766 = 2·3·31² (fields Q(ζ₁₅)⁺, Q(ζ₁₆)⁺), predicted first; at 37, dimension 3 at level 8214 with Q(ζ₁₈)⁺, predicted first | new | Dig 1, ninth and tenth |
| The first digits of the constants for 13 to 31 are in print: Ono 2000 (13 to 23), Boylan–Webb 2013 (29), Narayanan's 2019 Harvard thesis (13 to 31; found 20 Nov 2026, after an outside review) | known | reports/Partition tower literature check.md |
| **The constants: α = (3\|ℓ)·u_f(−2)**, the families through the level-6 forms of weight ℓ − 3, read at weight −2. Known to 30 digits (13–31) and 20 digits (37–61) | new | Dig 1, 12th to 18th; STATEMENT A |
| 56 digits of those constants confirmed on exact partition numbers, including one with 404,903 digits | counted | Dig 1, 12th and 16th |
| 47 and 61 (three forms from 61 on) predicted, then confirmed on lanes 47³ and 61² | new | Dig 1, 16th and 17th |
| **73 and 97** (the primes ≡ 1 mod 24 below 100): classes of 3 and 4 forms, found exactly; the laws from the families to 10 digits (97 has one 97-adic constant, 28, 26, 81, …, and a cubic). 73's law holds on exact partition numbers along ℓ (10 of 10, levels 1-7, numbers up to ~11 million digits) and both primes' laws along Q = 5, 7, 11, 13 (12 of 12 chains each). **No one-level law**, although the indexing allows one there: the constants are non-squares, so no one-level operator exists (derived before the check; then counted: at 73 the least single-step law is P(x²)) | new; counted | The towers of 73 and 97 |
| **Laws modulo 73 and 61** (Ramanujan-type congruences, Q³ and Q⁴ kinds): 25 and 35 primes Q below 5·10⁷, e.g. p(73·11367817³·n + 9040523964878168) ≡ 0 (mod 73); the count follows the 3-in-ℓ³ law; 6 class laws checked on exact partition numbers (12 of 12 members, controls exact), Ono-type ones through their step (16 of 16) | new as far as we read (published lists stop at 31 and 37) | Lucky primes for 73 (`partitions_flint_lucky_dual`, checks) |
| **Laws at 67, 71, 79, 83**: 92 primes Q below 5·10⁷ (89.3 expected); with 73 and 61, 152 against 151.9; e.g. p(19728144816847346813·n + 6177904563614) ≡ 0 (mod 67); 6 laws checked on exact partition numbers (12 of 12 members, controls exact) | new as far as we read | Laws at the other primes with three forms (`partitions_flint_lucky_dual`, `partitions_flint_lucky_check_67_83`) |
| **The braided forms taken apart**: steps with numbers from the group's own field (degree 2 or 3) cancel all members but one; each member's multiplier carries (b + 1)/2 digits of its own constant and none of its conjugates' (108 of 108; 79 of 79), at 47, 59, 61, 73, 97 | new as far as we read; derived from Statement B | The braided forms taken apart (`partitions_braid_split`) |
| **97 split in two**: one step at 5 cancels the three conjugate forms; what is left is one form's tower, which obeys one-term laws (two levels along 97 multiply it by 28: 20 of 20 at levels 1 and 3; each other prime one number: 60 of 60) and has lucky primes at about 3 in 97 (53, 263, 607, …; 18 of 18 zeros and controls); its second digit 26 read off partition numbers at levels 3 → 5 (10 of 10) | new as far as we read; derived from Statement B | 97's tower is two towers (`partitions_97_split`) |
| **What fades is the next family** (13, 17): the tower's fading part, multiplied by ℓ²·u every two levels, is the weight −2 member of a slope-2 family of level-6 forms with the class's signs (another system on lane ℓ); its constant from the forms equals the one read off partition numbers in 6 of 6 digits at both primes; its size ℓ² is the towers' speed | new as far as we read | What fades; Where the fading notes live; The fading form (`tower_fading`, `fading_forms2`, `fading_carry`); STATEMENT A4 |
| **Each twin alone** (37, 41, 43, 53): one step at 5 removes either form; each twin's multiplier between levels b and b + 2 carries (b + 1)/2 digits of its own constant and none of the other's (147 of 147), 4 digits at 37 and 41 | new as far as we read; derived from Statement B | The twins taken apart (`partitions_twins_split`) |
| Steps along another prime Q: the T(Q²) multipliers are the weight −2 member's. New laws: 2999 for 37; 4049 and 8941 for 43 | new | Dig 1, 13th; STATEMENT B |
| 5, 7, 11: the laws are the slopes 2, 1, 2 of their families; constants to 33, 20 and 12 digits | new as stated | Dig 1, 14th and 21st; STATEMENT C |
| Boylan–Webb's printed relations for 53, and their periods 36 and 10 for 37 and 41, come out of our constants | known, matched | Dig 1, 18th |
| **607**, the one prime ≤ 1297 with a rank one short: two conjugate forms whose numbers at 607 fall into F₆₀₇, so T₆₀₇'s minimal polynomial has degree 24 | explanation new (the rank deficit is Boylan–Webb's) | Dig 1, 18th; STATEMENT D |
| Lucky primes (where such collisions happen) follow a Chebotarev density; counted to a million | new | Dig 1, 19th; STATEMENT E |
| **The second kind (a non-ordinary class form) can never happen**: theta cycles turn it into a weight-6 form on Γ0(6) with U₂ ≡ 4, U₃ ≡ 9, and none exists. The one weight-6 newform h₆ predicts non-ordinary forms only in the class (−(2\|ℓ), +(3\|ℓ)), only where h₆ is non-ordinary; below 10⁷ that is only at 859. **Confirmed: at 859 exactly that class has one; nothing at 13–300 or 853, 857, 863** | new; the reduction is derived, the test counted | The second kind of exception |
| Yang's published Theorem 1 (Nagoya 2014) works out only n ≡ 1, 11 (mod 24); that already covers the tower primes ℓ ≡ 1, 11 (mod 24) | counted (read), derived | NEXT.md, link 4 |
| The same ordinarity for 3 and 5 colours (weights 8 on Γ0(2) and 10 on Γ0(6) decide it; only primes 2, 3, 5 could fail) | new; derived, inputs exact | Second kind … (colours); STATEMENT G4 |
| Chen–Du–Hou–Sun's representatives with their exact weights confirmed on numbers (2,942 of 2,942 coefficients, up to ℓ = 13, b = 3) | counted | Link 1 … checked on numbers |
| A complete draft note: Theorems A and B with every ingredient cited | drafted | NOTE.md |
| The proof's eigenvalue step rests on Ahlgren–Andersen–Dicks (Adv. Math. 2024), valid at every prime. Yang's omitted cases are not needed | quoted, applied | NOTE.md §4.3 |
| The loops: return depth always 1; the orders c_m = c₁·ℓ^(m−1); otherwise random | new | Dig 1, 20th; STATEMENT F |
| PARI 2.15.4 gives a wrong p(56,628,877,314); the tower law caught it, and sympy has the right value | new (bug) | Dig 1, 13th (line ~1009) |

## Other numbers of colours, and other functions

| finding | status | record |
|---|---|---|
| 3 and 5 colours: the same recipe, at level 2 with twist (−1\|ℓ) when 3 \| r, else level 6 with (3\|ℓ); censuses equal the rank bound; digits confirmed | levels and twists follow from Yang's theorems; the constants are new | Direction 1; STATEMENT G |
| 3 colours: at 11 a pair of slope-2 forms whose constants agree on lane 11, giving a three-term law checked to 11⁶; at 17 one slope-1 form, law 17^⌈b/2⌉. 5 colours: the law 11^(⌊b/2⌋+1) at 11; (5, 23) superexceptional | new as stated (superexceptional is Boylan's term, 2004) | Direction 1, continued |
| 24 colours (1/Δ): Eisenstein 1 plus level-1 cusp families. The tower stops changing at 5, 7, 11, 17, 19, 23. **A 607-type exception at 43** (field Q(√51349)), predicted from the forms, then confirmed | new | Where else: 24 colours |
| spt: constant (3\|ℓ) (ABL 2011) and (3\|ℓ)·ℓ at 5, 7, 13 (Garvan 2012). **Our measured depths equal their moduli exactly** | known, sharpened | Where else: spt |
| j: for n ≥ 1 the cusp forms alone. Lehner's congruences are the case "no cusp form"; constants 8, 7, 6, 5 on lane ℓ at 13, 17, 19, 23 | 5–11 rediscovered (Lehner); 13 possibly Atkin–O'Brien (unread); 17+ new | j |
| j at 31: 744 = 24·31 hides the Eisenstein constant on lane 31 | counted | j |
| 1/E₄ at 37: all four constants ≡ 1 on lane 37, because 37 divides B₃₂ (and u⁻⁶ ≡ 1) | counted, explained | 1/E₄ |
| The laws of distinct parts, overpartitions and T_2B, found from the data alone (no forms) | new as stated | How wide … (boundary) |
| T_2B's constants are the forms' (at 11: a₁₁ of the weight-10 newform on Γ0(2)); overpartitions' constants are the Eisenstein 1, level-1 cusp forms, and level-2 newforms with w₂ = (2\|ℓ) (fits 5–19) | new; counted on lane ℓ | T_2B …; Overpartitions … |
| **Distinct parts' constants are the forms'**: the two-level step is T_ℓ∘T_ℓ on (η(2z)/η(z))·{forms of Γ0(2) vanishing at the cusp 0}, the Hecke operators carrying the multiplier class r to rℓ mod 24 and back. It reproduces the three laws the data found alone (dying at 5, −1 at 7, −2 with a dying part at 11). At 7 the weight-0 constant matches the tower to 7 digits (n < 10⁸). At 13–23 the laws are predicted (13: (y + 3)(y + 5)(y + 6) on lane 13, holding on the one check 10⁸ allows) | new; derived and partly confirmed | Distinct parts: the constants from forms |
| **More moonshine: T_2A and T_3C = j(3τ)^(1/3).** T_2A dies at 3, 5, 7 (no Eisenstein constant: its Fricke symmetry shuts it out); from 11 its constants are the level-1 cusp forms and the level-2 newforms with w₂ = (2|ℓ), matching the data-only law at every prime 11–31. T_3C dies at 5, 7, 11; from 13 its constants are the forms η¹⁶·M_(ℓ−9), one-level for ℓ ≡ 1 mod 3 and two-level for ℓ ≡ 2 mod 3, matching every data-only law at 13–43 | new as far as we read (literature not checked) | More moonshine |
| **A mock theta function's tower (Ramanujan's f(q)) has laws; its unit constant is exactly the unit root (−3|ℓ) of its shadow's Hecke polynomial** y² − χ₁₂(ℓ)(ℓ + (−1|ℓ))y + (−1|ℓ)ℓ, (its slope-1 root's presence was over-read: the slope-1 piece really there has another unit part, 2 at 5 and 4 at 7, unidentified): −1 at 5 and +1 at 7 to every digit the levels settle. A second unit constant, (2|ℓ), comes from the level-6 newform family of weight ℓ − 1 in the class ((2|ℓ), (3|ℓ)), whose value at weight 0 is exactly 1 to 9 digits; it explains the double root at 11 and the extra −1 at 13. The first prediction (the whole eigenvalue, −6) was wrong and is recorded; the corrected one held at 7, read afterwards. A further piece of the same slope is not identified | new as far as we read (literature not checked) | A mock theta function |
| 5 colours at 11: the law follows one slope-1 family; constant 2, 8, 5, 2, 10, 5, 2, 0, 1, 1, 3, …, matching the exact numbers | new; counted | 5 colours at 11 |
| Double poles: E₄/j reproduces Zhang's theorem ((−3\|ℓ)·ℓ) at every prime 5–43. **For 1/E₄² the whole law P(x)·(x − ℓ·u⁻¹²)·S(x) holds to the whole lane at 7, 13 and 19, with no constant read from the tower.** S is the slope-1 families at weight −8, computed directly by Lauder's algorithm and checked against the classical forms at weight 148. At 13 that is four slope-1 constants. | E₄/j known (Zhang), our check counted; 1/E₄² new | Double poles; What 13 needs |
| Zhang's conjectures for E₄/(j − j(C))^r (arXiv 2510.23200) at j(C) = 2, 13 primes from 5 to 307: Theorem 1.3 and Conjectures 1.4 and 1.5 hold in every case, and 1.5's bound 3l − 2 − r is met exactly (sharp). **At supersingular primes the minimal polynomial X² − p² reaches 3l − r (2 deeper) and contains Conjecture 1.4 as its first level.** At ordinary primes the next roots are our p³u⁻⁴ and E₄'s p³ (4 a level at 5–17), plus, at 19 and 23, θ³ of a weight −2 Hida form (λ ≡ 13 and 7 from Lauder; at 19 the leftover runs exactly like n³·a_n(ΔE₄), as predicted). So the tower is Zhang's three roots plus p³ times the weight −2 world, which is why his bound grows 3 a level. | conjectures Zhang's; sharpness, the supersingular law and the slope-3 roots not in his paper (literature not checked); counted | Zhang's conjectures |

## Open

- The leftover extra multiples: 20.64% of p(N) are multiples of 5 and 14.46% of 7 where no law reaches (24N − 1
  square-free), and 8.94% of 13 in its class. It is not class numbers. Delaunay's Tate–Shafarevich heuristic gives
  0.2067 and 0.1454, a resemblance only. (README "The leftover, dug into".)
- At 1/E₄ and ℓ = 37 the pole's constant cannot be read: everything is ≡ 1 there.

## Corrected (said wrong first; the record keeps both)

- "At 13 and 17 the slope-2 root merges into a crowd; at 29 and 31 one long Newton segment" (Where the fading notes
  live): WRONG reading. Each such segment's slopes add up to exactly the working precision, so the roots were hidden by
  too few digits, not merged. Splitting the forms by their Hecke numbers on lane ℓ shows them (README "What fades at
  19 and 23").
- The page "The Exception 607":
  - Version 1 drew the pair at 5 side by side only because of the reference chosen.
  - Its chance estimate was explained by one mechanism, though the number counts every way; the wording is fixed, the
    number unchanged.
  - Its "they part one digit deeper" was shown with the forms at weight 604. Computed for the tower itself, it holds
    (version 4).
  - README "The review of the page, assessed point by point" has the rest.

- "Adding our next root p³u⁻⁴ makes Zhang's law grow 4 a level" (a prediction): wrong alone. E₄'s critical root p³
  must go with it, and at 19 and 23 a third root as well (README "Zhang's conjectures").
- "No form of lower level can share the tower's numbers": WRONG. It happens at 17, 29, 37 and 499, and is separated
  exactly (Dig 1, 17th and 18th).
- "At 7 the slope-1 family of 1/E₄² has unit part 4, 0, 1, 3" (a two-term fit, Double poles): its 4th digit was wrong. The forms give 4, 0, 1, 6, and only that value makes the whole law reach the full lane (What 13 needs).
- "j = 2's constant is reached one digit per level because the pole drifts" (A pole at a point with no symmetry): true at 5 only; at 7, 11 and 19 a slope-1 family alone gives that growth. The drift itself is confirmed once the families are in the law (Curves that repeat themselves).
- "The mock theta tower's two-level constant is the shadow's whole T(ℓ²)-eigenvalue χ₁₂(ℓ)(ℓ + (−1|ℓ))" (−6 at 5): WRONG at 5; it is the unit root (−3|ℓ) (A mock theta function).
- "The shadow's slope-1 root χ₁₂(ℓ)ℓ is present in the mock theta tower" (one digit deeper): over-read; the slope-1 root really there has unit part 2 at 5 and 4 at 7 (A mock theta function, "Measured, the candidate family").
- "607 has a non-ordinary form": REFUTED (Dig 1, 18th).
- "25 levels deeper" became two levels (Dig 1, third).
- Berndt–Ono's 0.174 for 11 is laws plus chance (Dig 1, again).
