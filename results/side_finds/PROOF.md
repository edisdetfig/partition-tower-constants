# Toward a proof: the chain, its sources, and its gaps (12 Nov 2026)

**In plain words.** A proof of the formula for the partition tower's constants would chain together four published
results. Three of them are proven theorems. The fourth, which says *where* the relevant forms live (level 6, with a twist),
is stated in Calegari's lecture notes without proof, though its classical version is a theorem of Yang. **Update, same day:**
that open step can be bypassed by a classical argument (the section "A classical route around step 3" below). At every
finite level the tower is congruent to classical forms, where Yang's theorem applies, and Hida's classical theory carries
the numbers to weight −2. What remains is for an expert to check that argument and two published details. Two further
points, that the tower uses *all* the forms and that each family has a well-defined value at weight −2, are known
computationally for every prime up to 1297 but not proven in general. So the formula is "what the published theory
predicts", not yet a theorem. The gaps are specific, and each has a plausible route. Everything below is sourced in
`research_notes/Proof ingredients/ingredients.md` (quotes, theorem numbers, links); "[open]" marks a gap.

## The statement (r = 1, ordinary partitions)

For a prime ℓ ≥ 13, let S_ℓ be the newforms of weight ℓ − 3 on Γ0(6) with w₂ = −(2|ℓ), w₃ = −(3|ℓ) (STATEMENT.md). Assume
every f in S_ℓ is ℓ-ordinary.
- Checked: all ℓ ≤ 1297.
- *Update 12 Nov:* true for every ℓ ≥ 7, by a theta-cycle argument through the one weight-6 newform h₆ of level 6 (README
  "The second kind of exception"; derived, and its prediction at 859 confirmed). So the assumption is now an argument to
  write out, not a computation.

Then the operator that moves the tower's ℓ-adic limit two levels
up is diagonalizable over the ℓ-adic closure, with eigenvalues (3|ℓ)·u_f(−2) for f in a subset of S_ℓ. It is all of S_ℓ
exactly when the limit has full rank R_ℓ = ⌊(ℓ + 12)/24⌋: every ℓ ≤ 1297 except 607, by Boylan–Webb's computation.

## The chain

**Step 1. The tower is U(ℓ²) applied to η(24z)⁻¹, and its limit is ordinary. Proven.**
- Folsom–Kent–Ono, Lemma 2.1 and (2.1)–(2.5), with Boylan–Webb's X(ℓ), Y(ℓ): P_ℓ(b + 2; 24z) = P_ℓ(b; 24z)|U(ℓ²). The
  even levels start from η(24z)⁻¹ (weight −1/2, Γ0(576), χ₁₂), the odd ones from η(24z)⁻¹|U(ℓ).
- The appendix to FKO (Theorem 7.1, by Ramsey; proof due to Calegari): the iterates converge ℓ-adically into the ordinary
  part of overconvergent forms of that weight, a finite-rank Hecke-stable module Ω, via the projector lim U(ℓ²)^(s!).
- So the constants are U(ℓ²)-eigenvalues on the part of Ω the tower spans.

**Step 2. Ramsey's overconvergent Shimura lift matches those eigenvalues with U_ℓ at weight −2. Proven.**
- Ramsey, IMRN 2009: a rigid-analytic map from the half-integral weight cuspidal eigencurve (tame level 4N) to the
  integral weight one (tame level 2N); weight κ ↦ κ², so −1/2 ↦ −2; T_ℓ ↦ T_ℓ², U_p ↦ U_p², and g(κ, α) = (κ², α).
- So the U(ℓ²)-eigenvalue of a half-integral eigen-system equals the U_ℓ-eigenvalue of its lift, with no twist.

**Step 3. The lifts are level-6 forms twisted by χ₁₂, in the class of S_ℓ. [open in the overconvergent setting]**
- Calegari's Lemma 4.6.4 (AWS 2013 notes): the image of S_(−1/2)(1, r) lies in S†_(−2)(Γ0(6), r) ⊗ χ₁₂, with U₂ = 2⁻²,
  U₃ = 3⁻². No proof; he notes "fairly soft proofs" via the trace formula.
- The classical version is a theorem: Yang (Nagoya Math. J. 215 (2014), Theorem 1). For (r, 6) = 1, 0 < r < 24,
  η(24z)^r·M_s ≅ S^new_(r+2s−1)(6, −(8/r), −(12/r)) ⊗ (12/·) as Hecke modules, with multiplicity one.
- With r ≡ −ℓ (mod 24) this is exactly S_ℓ ⊗ χ₁₂ (derived). The twist gives the factor χ₁₂(ℓ) = (3|ℓ) on U_ℓ.
- **Route to close it:** take the classical points of Ramsey's half-integral eigencurve at the weights of the tower's
  component; there, Yang's Theorem 1 places their lifts in S^new(6, …) ⊗ χ₁₂. Classical points are dense (Ramsey's and
  Coleman's control theorems), and "U₂ = ±2^(k/2−1), U₃ = ±3^(k/2−1), twisted by χ₁₂" is a closed condition on the
  eigencurve. So it should hold on the whole component, including weight −1/2 → −2. (A route, not a proof: Yang's theorem
  is for holomorphic forms of positive weight; the density argument must be run on the right component; and Yang's proof
  in print omits some residue classes.)

**Step 4. At weight −2 the ordinary U_ℓ-eigenvalues are the values u_f(−2) of the Hida families through S_ℓ.
[assembled, not in one place]**
- Hida (Astérisque 298, Theorem 3.1 and the argument before it): the ordinary part of Katz's ℓ-adic forms is Λ-free.
  Its specialization at every weight κ in the component (z ↦ z⁻² included) has the same rank as at classical weights.
- Coleman (Invent. Math. 127, Theorem B3.2): the characteristic series of U is analytic in the weight, with Iwasawa-function
  coefficients. So the slope-0 eigenvalues at κ = −2 are the values there of the slope-0 eigenvalues as functions of
  κ: the u_f(κ). This is what we computed by Newton's series through k = ℓ − 3 + (ℓ − 1)x, x = 0, 1, 2, … to x = −1.
- Needed besides (derived from our checks, not from print):
  - every f ordinary: checked to 1297, and argued for every ℓ ≥ 7 on 12 Nov (second_kind.py);
  - each family a single analytic branch through x = −1 (every n-th difference divisible by ℓⁿ in every run).

**Step 5. Which f occur. [open in general; computed to 1297]**
- FKO's Theorem 7.1 gives "⊆" only. The tower spans everything iff its rank is R_ℓ.
- Boylan–Webb computed rank R_ℓ for every ℓ ≤ 1297 except 607 (rank R − 1). Our reading (STATEMENT.md, D1): the levels
  fill at most deg(minimal polynomial of T_ℓ on S_ℓ, on lane ℓ) directions, and at 607 that degree is 24.
- So "all of S_ℓ" holds for the computed range; for 607 one eigen-system's place is taken by the repeated root.

**Further gaps** (from the ingredient notes):
- Ramsey's lift matches eigen-systems, not multiplicities. Yang's multiplicity one is classical; a Λ-adic or half-integral
  Hida theory (Hida 1995, not read) would carry it over.
- Cuspidality: η(24z)⁻¹ has poles. FKO's appendix shows U(ℓ²) strictly lowers them, so the iterates become ℓ-adic forms;
  Calegari asserts S†. Needs a line of proof.

## A classical route around step 3 (derived, 12 Nov 2026; for an expert to check)

**The idea.** Calegari's lemma is needed only if one works with the tower's limit, an overconvergent form of weight −1/2.
But at every finite level the tower is congruent to a classical form of large positive weight, where Yang's theorem
applies directly. Those weights converge ℓ-adically to −(r + 1), and Hida's theory at classical weights carries the
eigenvalues along. No overconvergent forms, no Ramsey lift and no Calegari lemma are needed. For r colours
(r = 1 is ordinary partitions), with P_r(b; z) = Σ p_r((ℓᵇn + r)/24)·qⁿ:

1. **Classical representatives** (Chen–Du–Hou–Sun, arXiv 1206.6642, Theorem 3.1; Ahlgren–Boylan for r = 1):
   P_r(b; 24z) ≡ η(24z)^γ_b · φ_b(24z) (mod ℓᵇ), with φ_b ∈ M_λ_b(SL₂(ℤ)), γ_b ≡ −r·ℓ^(−b) (mod 24) and 0 < γ_b < 24.
2. **The tower step** is exact on q-series: P_r(b + 2; z) = P_r(b; z) | U(ℓ²).
3. **U and T agree to high order.** On η(24z)^γ·M_λ, of weight λ′ + 1/2, T(ℓ²) = U(ℓ²) + (terms divisible by
   ℓ^(λ′−1)), by the formula for half-integral Hecke operators. The weights λ_b grow like ℓᵇ, so modulo any fixed ℓ^m
   U(ℓ²) ≡ T(ℓ²) there.
4. **The tower is a Hecke orbit.** Hence, modulo ℓ^m, P_r(b + 2j) ≡ g_b | T(ℓ²)^j for every j ≥ 0, with
   g_b = η(24z)^γ_b·φ_b(24z) (induction on j, using steps 2 and 3).
5. **Yang's theorems name the eigenvalues.** On η(24z)^γ·M_λ, T(ℓ²) has the eigenvalues ψ(ℓ)·a_ℓ(f), with multiplicity
   one, where f runs over the newforms of weight K_b = γ + 2λ − 1 with Atkin–Lehner signs set by γ. Theorem 1 covers
   (γ, 6) = 1: level 6, ψ = (12/·), signs −(8/γ), −(12/γ). Theorem 2 covers 3 | γ: level 2, ψ = (−4/·).
6. **Only the ordinary part survives.** Iterating T(ℓ²) kills, modulo ℓ^m, every f with a_ℓ(f) divisible by ℓ. For
   ordinary f, a_ℓ(f) agrees with its family's U_ℓ-value A_f(K_b) to ℓ^(K_b − 1).
7. **The weights converge to −(r + 1).** From Chen–Du–Hou–Sun's λ: K_b = r·(ℓ^(b−1)(ℓ − 1) − 1) − 1 for odd b (similarly
   for even b), so K_b ≡ −(r + 1) (mod (ℓ − 1)ℓ^(b−1)). Hida families are continuous in the weight, so
   A_f(K_b) ≡ A_f(−(r + 1)) to about ℓᵇ.
8. **The families are the class's families.** Hida's control theorem at classical weights matches the ordinary f at
   weight K_b with the families through weight ℓ − r − 2: the same component, the same count. Carried along the family,
   the signs of step 5 are exactly the class's (derived):
   - r = 1: w₂ = −(2|ℓ), w₃ = −(3|ℓ);
   - r = 5: both flipped;
   - r = 3: level 2, U₂ sign (2|ℓ).

**Conclusion.** Modulo every ℓ^m, the characteristic polynomial of the two-level operator on the tower's limit divides
∏_f (x − ψ(ℓ)·A_f(−(r + 1))), and is equal to it exactly when the tower has full rank. That is Statements A1 and G2,
with step 3 replaced by classical theorems. It also explains why the census equals the rank bound: by step 5 the count
is the number of ordinary forms.

**What an expert would still check.**
- Yang's Theorems 1–2 in full. The arXiv text works out only n ≡ 1, 11 (mod 24) and omits Theorem 2's proof; the Nagoya
  version should be read.
- The exact weights λ_b in Chen–Du–Hou–Sun (even b too).
- The limit with congruences between eigenforms (non-semisimple modulo ℓ^m): work with characteristic polynomials and
  the ordinary projector, which is a polynomial in T(ℓ²).
- Hida theory at tame level 6 (or 2) with the U₂, U₃ conditions. This is standard: they are Hecke conditions.
- Full rank (step 5 of the chain) stays a computation, and 607 shows it can fail.

## The first digit, written out (STATEMENT A3; 13 Nov 2026, after a review asked for the proof)

**Claim.** Let ℓ ≥ 13 be prime. On lane ℓ, the tower's two-level step acts as (3|ℓ)·T(ℓ) does on the class. So each
note's constant has first digit (3|ℓ)·a_ℓ(F), where F is the class's newform (level 6, weight ℓ − 3) the note belongs to.

**Published inputs** (quoted verbatim in `research_notes/Partition tower literature check/constants_and_families.md`,
Q5 and Q6):
1. **AAT (2.6):** "f_ℓ ∈ S_{(ℓ−2)/2}(1, ν_η^{−ℓ}, Z) with f_ℓ ≡ Σ p((ℓn+1)/24) q^{n/24} (mod ℓ)".
2. **AAT (2.4):** for f = Σ a(n) q^{n/24} in S_k(1, ν_η^r),
   "f|T(Q²) = Σ (a(Q²n) + Q^{k−3/2} (−1/Q)^{k−1/2} (12n/Q) a(n) + Q^{2k−2} a(n/Q²)) q^{n/24}".
3. **AAD Theorem 1.1, (1.4):** "S_t(T_{p²}F) = (12/p) T_p S_t(F) for each prime p ≥ 5", for the Shimura lifts from weight
   λ + 1/2 and level N (here N = 1).
   - **Yang 2014, Theorem 1, as quoted by ABR:** each eigenform has a newform in S^new_{2k−1}(6) "with the same Hecke
     eigenvalues as Sh_t(f) ⊗ (12/•) at all primes Q ≥ 5".
4. **AAT (2.8):** "each Shimura lift of f_ℓ is in the space S^new_{ℓ−3}(6, −(8/−ℓ), −(12/−ℓ)) ⊗ (12/•)".
   - (8/−ℓ) = (2/ℓ) and (12/−ℓ) = (3/ℓ), and a_p = −p^{k/2−1}·ε_p for p ∥ 6.
   - So these are the newforms with U₂ = (2|ℓ)·2^{(ℓ−5)/2} and U₃ = (3|ℓ)·3^{(ℓ−5)/2}: our class.

**The argument (derived).**
1. **The step is U(ℓ²), exactly.**
   - The tower at level b is P_ℓ(b) = Σ p(N)·q^{(24N−1)/ℓ^b}, over N ≡ δ_b (mod ℓ^b).
   - U(ℓ²) keeps the exponents divisible by ℓ² and divides them by ℓ², so P_ℓ(b + 2) = P_ℓ(b)|U(ℓ²).
   - P_ℓ(1) is AAT's Σ p((ℓn+1)/24)·q^{n/24}, in the variable q^{1/24}.
2. **On lane ℓ, U(ℓ²) is T(ℓ²).** Let S = S_k(1, ν_η^{−ℓ}, Z) with k = (ℓ − 2)/2.
   - In (2.4) with Q = ℓ, the second term carries ℓ^{k−3/2} = ℓ^{(ℓ−5)/2} and the third ℓ^{2k−2} = ℓ^{ℓ−4}.
   - Both are multiples of ℓ for ℓ ≥ 7, and the coefficients are integers. So T(ℓ²)g ≡ U(ℓ²)g (mod ℓ) for every
     g ∈ S, and T(ℓ²) maps S to S.
3. **By induction, P_ℓ(1 + 2j) ≡ f_ℓ|T(ℓ²)^j (mod ℓ).**
   - j = 0 is (2.6).
   - The step: P_ℓ(3 + 2j) = P_ℓ(1 + 2j)|U(ℓ²) ≡ (f_ℓ|T(ℓ²)^j)|U(ℓ²) ≡ f_ℓ|T(ℓ²)^{j+1}. This uses 2 for g = f_ℓ|T(ℓ²)^j,
     and the fact that reducing mod ℓ commutes with U(ℓ²).
4. **So the notes on lane ℓ are T(ℓ²)'s eigenvalues mod ℓ**, on the part of S ⊗ F_ℓ that f_ℓ generates. The odd levels
   are f_ℓ, T f_ℓ, T² f_ℓ, ….
5. **Those eigenvalues are (12/ℓ)·a_ℓ(F).**
   - The Hecke operators are normal for the Petersson product, so S ⊗ C has a basis of eigenforms.
   - By input 3 at p = ℓ (allowed: ℓ ≥ 5 and ℓ ∤ 6) and input 4, each eigenvalue of T(ℓ²) on S is (12/ℓ)·a_ℓ(F) for a
     newform F of the class.
   - T(ℓ²)'s characteristic polynomial on the lattice S has integer coefficients and these roots. Reduced mod ℓ, it
     gives the eigenvalues on lane ℓ.
   - Since (12/ℓ) = (3|ℓ), this is A3.

**One point to read in the source.** For each note on its own, the eigenvalue form of Yang's theorem quoted above is
enough. The stronger use, "the constants' polynomial on lane ℓ is T(ℓ)'s characteristic polynomial on the class, with
multiplicity" (used to count notes, for example 607's 25), needs the correspondence to be one-to-one.
- The dimensions agree: both are ⌊(ℓ + 12)/24⌋ (FKO and Boylan–Webb for S; counted by us for every class from 37 to
  1,297).
- Yang's theorem should still be read in Yang 2014 itself (Nagoya Math. J. 215) for that statement. Our notes quote it
  only through ABR.

**Checks (counted):**
- the first digits read off partition numbers at 13, 17, 19, 23, 29 and 31 (six of six);
- the twins' (T, N) at 37, 41, 43 and 47;
- 73's real steps (240 of 240) and the braids' data;
- at 607, the repeated root, which gives Boylan–Webb's rank 24.

## E1, the laws at a lucky prime, written out (13 Nov 2026)

**Claim.** Let ℓ ≥ 13 and Q ≥ 5 be distinct primes. Suppose every newform f of the class has (3|Q)·a_Q(f) ≡ c (mod ℓ)
for one number c.
- If c ≡ ε·(−3|Q)·Q⁻² (ε = ±1), then p(N) ≡ 0 (mod ℓ) for every N with 24N − 1 = Q²·M, ℓ | M and (M|Q) = ε.
- If c ≡ 0, then p(N) ≡ 0 (mod ℓ) whenever 24N − 1 = ℓ·Q³·m with Q ∤ m.

**Published inputs:**
- AAT (2.6) and (2.4), as in "The first digit, written out";
- AAT Lemma 5.1: if g|T(Q) ≡ λ_Q·g (mod ℓ) for every g in the class, then f_ℓ|T(Q²) ≡ (12/Q)·λ_Q·f_ℓ (mod ℓ).

**The argument (derived).**
1. **The coefficients.** Write f_ℓ = Σ a(n)·q^{n/24}. By (2.6), a(n) ≡ p(N) (mod ℓ) when 24N − 1 = ℓ·n.
2. **The relation.** Lemma 5.1 with λ_Q = a_Q(f), read coefficient by coefficient through (2.4), with k = (ℓ − 2)/2 and
   c = (12/Q)·a_Q = (3|Q)·a_Q:
   a(Q²n) + Q^{k−3/2}·(−1/Q)^{k−1/2}·(12n/Q)·a(n) + Q^{2k−2}·a(n/Q²) ≡ c·a(n) (mod ℓ).
3. **Ono's kind (c ≡ 0).** Take n = Q·m with Q ∤ m. Then (12n/Q) = 0 and a(n/Q²) = 0, so a(Q³m) ≡ c·a(Qm) ≡ 0. With
   24N − 1 = ℓ·Q³·m, that is p(N) ≡ 0.
4. **The class laws.** Take n with Q ∤ n. Then a(n/Q²) = 0, so a(Q²n) ≡ (c − μ(n))·a(n), where
   μ(n) = Q^{(ℓ−5)/2}·(−1/Q)^{(ℓ−3)/2}·(3|Q)·(n|Q), using (12n/Q) = (3|Q)·(n|Q).
   - Mod ℓ, Q^{(ℓ−5)/2} ≡ (Q|ℓ)·Q⁻² (Euler's criterion).
   - If (ℓn|Q) = ε, that is (n|Q) = ε·(ℓ|Q), then μ(n) ≡ (Q|ℓ)·(ℓ|Q)·(−1|Q)^{(ℓ−3)/2}·(3|Q)·ε·Q⁻².
   - Quadratic reciprocity gives (Q|ℓ)·(ℓ|Q) = (−1)^{(Q−1)(ℓ−1)/4}, and (−1|Q)^{(ℓ−3)/2} = (−1)^{(Q−1)(ℓ−3)/4}. The
     product is (−1)^{(Q−1)(ℓ−2)/2} = (−1|Q), because ℓ − 2 is odd. So μ(n) ≡ ε·(−3|Q)·Q⁻².
   - Hence c ≡ ε·(−3|Q)·Q⁻² gives a(Q²n) ≡ 0 for every n with Q ∤ n and (ℓn|Q) = ε. With M = ℓn and
     24N − 1 = Q²·M, that is the claim.
5. **Checked (computed):** the sign identity in 4, and the whole middle coefficient against ε·(−3|Q)·Q⁻², for all 5,475
   pairs of primes with 13 ≤ ℓ < 400 and 5 ≤ Q < 400, Q ≠ ℓ. No mismatch.

## Other numbers of colours (Statement G)

- Yang's two theorems already contain the classical structure we found numerically:
  - Theorem 1, (r, 6) = 1: level 6, twist (12/·), signs −(8/r), −(12/r).
  - Theorem 2, 3 | r: level 2, twist (−4/·), sign −(8/(r/3)).
- With the tower's exponent r ≡ −(colours)·ℓ⁻¹ (mod 24), this gives:
  - 3 colours: level 2, twist (−1|ℓ), U₂ sign (2|ℓ);
  - 5 colours: level 6, twist (3|ℓ), both signs flipped against ordinary partitions.
  These are exactly the outcomes of our censuses (50 and 33 primes) and twist fits (7 and 8 primes) (derived).
- The same chain (Steps 1–5) then predicts Statement G, with the same gaps.

## What our computations add to the chain

- The constants themselves, to 20–30 digits for 13 to 61 (STATEMENT.md A; README, the eighteenth time), checked against
  exact partition numbers and Boylan–Webb's printed relations.
- Ordinarity of every class form for ℓ ≤ 1297; the single-branch behaviour (all differences divisible).
- The minimal-polynomial reading of 607.
- The r-colour analogues: censuses, twists, digits, and the slope families behind the laws (STATEMENT.md G).

## What it would take to finish

1. Step 3: check the classical route above (or, the original way, write the density argument on Ramsey's eigencurve).
   The classical route needs only Yang's Theorems 1–2 in full, Chen–Du–Hou–Sun's representatives and Hida's theory at
   classical weights.
2. A precise reference for ordinary ⇒ overconvergent at tame level 6 and weight −2, or an argument through Coleman's
   characteristic series. This is bookkeeping.
3. Multiplicity one along the families (half-integral Hida theory). This needs a source.
4. Accept Step 5 as a computation for ℓ ≤ 1297, or prove full rank in general. The latter is open: 607 shows full rank
   can fail.

**Update, 12 Nov (later; the full table is in NEXT.md).** Along the classical route, items 2 and 3 are not needed, and
item 4 becomes part of the statement (our reading). What remains:
- **Chen–Du–Hou–Sun's Thm 3.1 is read and confirmed:** F ≡ η(24z)^γ·φ(24z) (mod ℓᵇ), with 0 < γ < 24 and φ a level-one form
  of explicit weight.
- **U ≡ T** follows from Yang's formula for T(p²).
- **The gap is Yang's Thm 1 with its whole proof.** arXiv v1 works out the trace comparison only for n ≡ 1 and 11
  (mod 24), and the journal text is not compared.
- **Status at the end of 12 Nov: a complete draft with no known gap (NOTE.md).** Every step is a published theorem or a
  short argument of ours (Theorem A; the assembly; §8's routine points). What is left is a specialist's reading.
- **Later the same day: the gap is closed.** Ahlgren–Andersen–Dicks (Adv. Math. 447, 2024, Thm 1.1) give a lift into
  the level-6 new space with our signs, Hecke-equivariant at every prime p ≥ 5. With T(ℓ²) semisimple, that gives the
  tower's minimal polynomial dividing the class polynomial for every ℓ (NOTE.md §4.3). Yang's omitted cases are not
  needed.
- **To do ourselves:** precise citations for Hida's classical control and continuity; the theorem stated with the Λ-adic
  characteristic polynomial ("divides, with equality at full rank"); a complete note; a specialist's check.
