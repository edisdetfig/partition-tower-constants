# The ℓ-adic constants of the partition tower — draft note for a specialist (12 Nov 2026)

This note states two results about the partition function modulo powers of a prime ℓ and sets out their proofs from
published theorems. **Theorem A** (every form behind the tower is ordinary) has a complete argument here.
**Theorem B** (the tower's constants are values of Hida families at weight −2) has a proof from published theorems:
- Chen–Du–Hou–Sun for the classical representatives;
- Ahlgren–Andersen–Dicks (Adv. Math. 447, 2024) for the eigenvalues, valid at every prime;
- Hida for the families.

An earlier draft leaned on Yang's Theorem 1, whose written proof omits six of eight residue classes. That is no longer
needed (§4.3). Everything is labelled as one of:
- **[cited]**: a published statement, quoted;
- **[ours]**: an argument given here, to be checked;
- **[computed]**: a numerical check, with kept outputs in this folder.

## 1. Setting and statements

Let ℓ ≥ 13 be prime. For b ≥ 1 put

  F_b(z) = Σ p((ℓᵇn + 1)/24)·qⁿ, the sum over n ≥ 0 with ℓᵇn ≡ −1 (mod 24).

As q-series, F_(b+2) = F_b | U(ℓ²), where U(ℓ²): Σ a(n)qⁿ ↦ Σ a(ℓ²n)qⁿ. So the odd levels form one tower and the even
levels another (the "partition tower"). Folsom–Kent–Ono [FKO] proved that, modulo any ℓᵐ, the levels eventually lie in a
module of rank at most R_ℓ = ⌊(ℓ + 12)/24⌋. So each tower obeys a linear recurrence two levels apart with at most R_ℓ
constants.

Let k = ℓ − 3. Let S_ℓ be the set of newforms of weight k on Γ0(6) with Atkin–Lehner eigenvalues w₂ = −(2|ℓ),
w₃ = −(3|ℓ). Its size is ⌊(ℓ + 12)/24⌋ = R_ℓ [cited: Yang 2011, Remark; Guo–Ono]. For an eigenform f with coefficients
in a number field and a prime λ | ℓ, f is **ordinary at λ** if a_ℓ(f) is a λ-adic unit. For ordinary f, let u_f(κ)
denote the U_ℓ-eigenvalue, at ℓ-adic weight κ, of the Hida family through f's ordinary ℓ-stabilization.

**Theorem A [ours].** For every prime ℓ ≥ 7, every f ∈ S_ℓ is ordinary at every prime above ℓ.

**Theorem B [ours, from (I1), (I2′), (I3)–(I5)].** Let Q_ℓ(x) = ∏_{f ∈ S_ℓ} (x − (3|ℓ)·u_f(−2)). For every
m ≥ 1 there are b₀ and j₀ such that, for each parity class of b ≥ b₀ and all j ≥ j₀,

  Σ_i c_i·F_(b+2j+2i) ≡ 0 (mod ℓᵐ), where Q_ℓ(x) = Σ_i c_i xⁱ.

That is, the two-level step on the tower satisfies Q_ℓ modulo ℓᵐ. Its minimal polynomial on the tower is Q_ℓ itself
exactly when the tower has full rank R_ℓ. Boylan–Webb [BW] computed full rank for every ℓ ≤ 1297 except 607, where the
rank is R_ℓ − 1. §5 gives our account of 607.

## 2. Ingredients [cited]

**(I1) Chen–Du–Hou–Sun [CDHS, Thm 3.1].** For a prime m ≥ 5 and positive integers k, r there is
φ ∈ M_λ(SL₂(ℤ)) with F_(m,k,r)(z) ≡ η(24z)^γ·φ(24z) (mod mᵏ), where:
- F_(m,k,r) = Σ_{mᵏn ≡ −r (24)} p_r((mᵏn + r)/24) qⁿ;
- γ = (24β − r)/mᵏ, with β ≡ r/24 (mod mᵏ) and 0 ≤ β < mᵏ, so 0 < γ < 24 for odd r < mᵏ;
- λ = ((mᵏ − mᵏ⁻¹)/2)·r − (γ + r)/2 for odd k, and λ = (mᵏ − mᵏ⁻¹)·r − (γ + r)/2 for even k.

For r = 1 and the prime ℓ: γ_b ≡ −ℓ⁻ᵇ (mod 24), which is prime to 6.

[computed: checked with these exact γ and λ at (ℓ, b) = (13, 1–3), (17, 1–2), (19, 2), (23, 2), (29, 2). φ fitted to
M_λ(SL₂(ℤ)) from its first dim M_λ coefficients agrees with every further coefficient up to 400: 2,942 of 2,942
(`cdhs_check_….txt`).]

**(I2) Yang [Y, Thm 1].** For (r, 6) = 1, 0 < r < 24 and even s ≥ 0, the space S_(r,s) = {η(24τ)^r·g(24τ) : g ∈ M_s(1)}
is isomorphic, as a module for the Hecke operators T(n²) with (n, 6) = 1, to S^new_(r+2s−1)(6, −(8/r), −(12/r)) ⊗ (12/·).
Multiplicity one holds. Here T(n²) corresponds to T_n twisted by (12/n). The proof compares traces:
"tr T_(n²) = (12/n)·tr T_n … The verification is done case by case according to the residue of n modulo 24"
(Nagoya Math. J. 215, p. 58), and only n ≡ 1, 11 (mod 24) are written out. See §4.3.

**(I2′) Ahlgren–Andersen–Dicks [AAD, Thm 1.1]** (Adv. Math. 447, 2024; arXiv 2307.07438).
- Take λ, N ≥ 1 with (N, 6) = 1, r with (r, 6) = 1, ψ mod N, and F = Σ_{n ≡ r (24)} a(n)q^(n/24) ∈ S_(λ+1/2)(N, ψν^r),
  where ν is the η-multiplier. For λ = 1, F must also be orthogonal to theta series.
- For squarefree t > 0, the lift S_t(F) = Σ b(n)qⁿ, with Σ b(n)n^(−s) = L(s − λ + 1, ψ(·/t))·Σ (12/n)a(tn²)n^(−s),
  lies in S^(new 2,3)_(2λ)(6N, ψ², ε₂, ε₃), with ε₂ = −ψ(2)(8/r) and ε₃ = −ψ(3)(12/r).
- And "S_t(T_{p²}F) = (12/p)T_p S_t(F) for each prime p ≥ 5": no restriction on p mod 24.
- With N = 1 and ψ = 1 this covers our spaces: F(z) = η(z)^γ·φ(z), weight γ/2 + λ_b, has multiplier ν^γ on SL₂(ℤ), and
  2λ = γ + 2λ_b − 1 = K_b.
- Boylan and Swati (arXiv 2505.01018, 2025) also give "a constructive proof of Yang's result" (their abstract).

**(I3) Hecke operators of half-integral weight.** On S_(r,s), of weight κ + ½ = r/2 + s, for p ∤ 6:

  T(p²): a(n) ↦ a(p²n) + (12/p)((−1)^κ n/p)·p^(κ−1)·a(n) + p^(2κ−1)·a(n/p²)   [Y, §2].

So T(ℓ²) ≡ U(ℓ²) (mod ℓ^(κ−1)) on S_(r,s).

**(I4) Atkin–Lehner.** A newform of weight k on Γ0(N) with p ∥ N has a_p = −w_p·p^((k−2)/2).

**(I5) Hida theory [H86a, H86b; W88].**

Citations as given by Ghate–Kumar (arXiv, "Control theorems for ordinary 2-adic families of modular forms", read 12 Nov):
- the ℤ_p-rank of the ordinary part is constant along the weights, for p ≥ 5: "proved in [Hid86a, Thm. 3.1 and
  Cor. 3.2]";
- the ordinary Hecke algebra "is independent of the weight, for all k ≥ 2" [Hid88b, Thm. 3.2];
- the control isomorphism h⁰(N) ⊗_Λ Λ_P/P ≅ h⁰_k is their Theorem 8.1, stated for p = 2 after Hida's for p ≥ 5.

What is used:
- For tame level 6 the ordinary parts of S_K(Γ0(6)) ⊗ ℤ_ℓ, for K ≥ 2 in a fixed class mod ℓ − 1, are the
  specializations of one finite flat Hecke algebra over Λ = ℤ_ℓ[[T]].
- Hence their rank is constant in K.
- Their U_ℓ-eigenvalues A(K) are values of elements of that algebra, so they are continuous in K: K ≡ K′
  (mod (ℓ − 1)ℓ^(t−1)) gives A(K) ≡ A(K′) to an ℓ-adic precision growing with t.
- The value at κ = −2 is u(−2).

**(I6) Mod-ℓ forms [Katz 77; Jochnowitz 82; Edixhoven 92, §3].** The drop statement is recalled, e.g., in de Shalit–Goren
(arXiv 1412.5494): "ω(θ(f)) < ω(f) + p + 1 if and only if ω(f) ≡ 0 mod p", and after a drop the filtration is
ω(f) + p + 1 − a(p − 1) with a > 0.

Katz's theorem is for a fine level. Γ₁(6) is one, and division by the Hasse invariant respects the Γ₀(6) structure. So it
applies to our forms on Γ₀(6).

For a nonzero mod-ℓ form h of level N prime to ℓ, let w(h) be its filtration and θ = q·d/dq. Then:
- w(θh) ≤ w(h) + ℓ + 1, with equality if and only if w(h) ≢ 0 (mod ℓ);
- w(h) is ≡ the weight of h (mod ℓ − 1);
- θ^(ℓ−1)h = h − V(Uh), where U, V are the coefficient operators n ↦ ℓn;
- U(θh) = ℓ·θ(Uh) ≡ 0;
- for p ∤ Nℓ, T_p(θh) = p·θ(T_p h).

**(I7) Lifting [Deligne–Serre, Lemme 6.11].** A system of Hecke eigenvalues mod λ occurring on M_K(Γ0(N); F̄_ℓ) occurs,
up to reduction, on M_K(Γ0(N); Q̄_ℓ). For K ≥ 2 and ℓ ∤ 6N, M_K(Γ0(N); F_ℓ) is the reduction of M_K(Γ0(N); ℤ_(ℓ)).

## 3. Proof of Theorem A [ours]

Let ℓ ≥ 7, f ∈ S_ℓ, λ | ℓ, and suppose a_ℓ(f) ≡ 0 (mod λ). Write f̄ for its reduction, a mod-λ eigenform of weight
k = ℓ − 3 and level 6.

**Step 1 (the filtration of f̄).** w(f̄) = ℓ − 3: it is ≤ ℓ − 3, ≡ ℓ − 3 (mod ℓ − 1), and positive.

**Step 2 (f̄ lies on a theta cycle).** T_ℓ f̄ = U f̄ + ℓ^(k−1)V f̄ ≡ U f̄, so U f̄ = a_ℓ(f)·f̄ = 0, and θ^(ℓ−1) f̄ = f̄ by (I6).
In particular θ^j f̄ ≠ 0 for all j.

**Step 3 (the first four steps).** Write w_j = w(θ^j f̄). By (I6):
- w₀ = ℓ − 3, w₁ = 2ℓ − 2, w₂ = 3ℓ − 1, w₃ = 4ℓ, since ℓ − 3, 2ℓ − 2 and 3ℓ − 1 are ≢ 0 (mod ℓ).
- As w₃ ≡ 0, w₄ = 5ℓ + 1 − s(ℓ − 1) for some s ≥ 1.

**Step 4 (s ≥ 5, from closing the cycle).** Suppose s ≤ 4.
- Then w₄ ≡ 1 + s (mod ℓ), and w_(4+j) ≡ 1 + s + j (mod ℓ). For 0 ≤ j ≤ ℓ − 6 this residue lies in [2, ℓ − 1], so it is
  never 0, and θ adds exactly ℓ + 1 at each of these steps.
- So w_(ℓ−1) = w₄ + (ℓ − 5)(ℓ + 1) = ℓ² + ℓ − 4 − s(ℓ − 1).
- But w_(ℓ−1) = w₀ = ℓ − 3 by Step 2, which forces s(ℓ − 1) = ℓ² − 1, i.e. s = ℓ + 1 > 4. A contradiction.

**Step 5 (s ≤ 5, from positivity).**
- g = θ⁴ f̄ is nonzero with constant term 0, so it is not constant and w₄ ≥ 1.
- Also w₄ ≡ w₀ + 4·2 ≡ 6 (mod ℓ − 1), since θ raises the weight by ℓ + 1 ≡ 2.
- For ℓ ≥ 7 the least positive integer ≡ 6 (mod ℓ − 1) is 6, so w₄ ≥ 6, i.e. s ≤ 5.

Hence s = 5 and **w(g) = 6**.

**Step 6 (the eigenvalues of g).** By (I6), g is a mod-λ eigenform of weight 6 and level 6 with:
- T_p g = p⁴a_p(f)·g for p ∤ 6ℓ;
- U g = 0;
- U₂ g = 2⁴a₂(f)·g and U₃ g = 3⁴a₃(f)·g, since U_p θ = p θ U_p.

By (I4), a₂(f) = −w₂·2^((ℓ−5)/2) = (2|ℓ)·2^((ℓ−5)/2). By Euler's criterion 2^((ℓ−1)/2) ≡ (2|ℓ) (mod ℓ), so

  U₂(g) ≡ (2|ℓ)·2^((ℓ+3)/2) ≡ (2|ℓ)²·4 = 4,  and likewise  U₃(g) ≡ 9  (mod λ).

**Step 7 (no form of weight 6 has these eigenvalues).** By (I7) the system of eigenvalues of g occurs on
M₆(Γ0(6); Q̄_ℓ). That space has dimension 7 [computed, `second_kind.py`, with the dimension formula]:
- **four Eisenstein series**, with T_p-eigenvalue 1 + p⁵. Their T_ℓ-eigenvalue is ≡ 1, not 0. Excluded.
- **two oldforms from h₃ = η(z)⁶η(3z)⁶ = q − 6q² + 9q³ + …**, the newform of weight 6 on Γ0(3). They have U₃ = +9, and
  U₂-eigenvalues the roots of x² + 6x + 32. A root ≡ 4 would need 4² + 6·4 + 32 = 72 ≡ 0 (mod λ), impossible for ℓ ≥ 5.
  Excluded.
- **one newform h₆ = q + 4q² − 9q³ + 16q⁴ − 66q⁵ − 36q⁶ + 176q⁷ + …** (U₂ = +4, U₃ = −9; the sign classes
  (U₂, U₃) = (±4, ±9) other than (+4, −9) are empty at weight 6). U₃ = −9 ≡ 9 would need ℓ | 18. Excluded.

There are no other cusp forms: S₆(SL₂(ℤ)) = S₆(Γ0(2)) = 0. So no eigensystem matches, which contradicts Step 6.
Hence a_ℓ(f) is a λ-unit. ∎

**Remarks.**
- **The other sign classes.** The same steps apply to every class of newforms of weight ℓ − 3 on Γ0(6). With
  (w₂, w₃) written relative to ((2|ℓ), (3|ℓ)), g needs U₂ ≡ −4·w₂(2|ℓ) and U₃ ≡ −9·w₃(3|ℓ):
  - (+, ±) needs U₂ ≡ −4: no lift (h₃ would need ℓ | 24), so always ordinary;
  - (−, +) needs U₂ ≡ 4 and U₃ ≡ −9: exactly h₆.
- **The converse for (−, +), and the test [computed].** If ℓ ≥ 7 divides a_ℓ(h₆), the cycle of h̄₆ read backwards gives
  a non-ordinary form of weight ℓ − 3 in the class (−, +). h₆ is non-ordinary below 10⁷ only at ℓ = 859. At 859 the
  class (−, +) has a form with T₈₅₉-eigenvalue 0 on lane 859, and the other three classes, the partition class among
  them, have none. No class has one at any prime from 13 to 300, or at 853, 857, 863.
  (`second_kind_N10000000_L300_at853_857_859_863.txt`.)
- **Other numbers of colours.** For r colours the class has weight ℓ − r − 2, and the same cycle leads to weight r + 5,
  through θ^(r+3), on Γ0(6) or on Γ0(2) when 3 | r. There it needs U_p ≡ ε·p^((r+3)/2) (ε = +1 for r = 1, 3; −1 for
  r = 5). Every cusp eigensystem of weight 6 and 10 on Γ0(6), and of weight 8 on Γ0(2), misses these values except
  modulo 2, 3 or 5 [computed, exact: `second_kind_colors.txt`]. So **Theorem A holds for 1, 3 and 5 colours**, at every
  prime ℓ ≥ 7, 11, 11 respectively.

## 4. Proof of Theorem B [ours, from (I1)–(I5)]

**4.1 Each level is classical.**
- By (I1) with r = 1 and k = b, F_b ≡ G_b (mod ℓᵇ), with G_b = η(24z)^(γ_b)·φ_b(24z) ∈ S_(γ_b, λ_b) (if λ_b is odd,
  M_(λ_b) = 0 and F_b ≡ 0).
- Its weight under Yang's map is K_b = γ_b + 2λ_b − 1. From (I1)'s formulas for λ:
  - K_b = (ℓ − 1)ℓ^(b−1) − 2 for odd b;
  - K_b = 2(ℓ − 1)ℓ^(b−1) − 2 for even b.

  So K_b ≡ −2 (mod (ℓ − 1)ℓ^(b−1)) exactly, K₁ = ℓ − 3 is the class's own weight, and K_b → −2 ℓ-adically [ours, from the
  formulas]. Example, ℓ = 13: b = 1 gives γ = 11, λ = 0, K = 10; b = 2 gives γ = 23, λ = 144, K = 310 = 2·156 − 2.
- The signs are −(8/γ_b) = −(2|ℓ)ᵇ and −(12/γ_b) = −(3|ℓ)ᵇ. That is the class S_ℓ's for odd b, and (−1, −1) for even b.

**4.2 The tower is a Hecke orbit.**
- F_(b+2j) = F_b | U(ℓ²)^j exactly. Since U(ℓ²) preserves congruences mod ℓᵇ, and by (I3) U(ℓ²) ≡ T(ℓ²) on S_(γ_b,λ_b)
  modulo ℓ^(κ_b−1) with κ_b − 1 ≥ b, which is stable under T(ℓ²) [Y; Garvan]:
  F_(b+2j) ≡ G_b | T(ℓ²)^j (mod ℓᵇ) for all j ≥ 0.

**4.3 The eigenvalues (from (I2′), for every ℓ ≥ 5).**
- Let F ∈ S_(γ_b,λ_b) be a Hecke eigenform, with T(p²)F = λ_F(p)·F for all p ≥ 5. Some a(n) ≠ 0 with n ≡ γ_b (mod 24),
  so n is prime to 6. Write n = tm², t squarefree. The Hecke relations make a(tm²) a multiple of a(t), so a(t) ≠ 0 and
  S_t(F) ≠ 0 [standard].
- By (I2′), (12/p)·T_p S_t(F) = λ_F(p)·S_t(F) for every p ≥ 5.
- By multiplicity one for newforms at level 6, S_t(F) is a multiple of one newform G of weight K_b, in the class
  (−(8/γ_b), −(12/γ_b)). So **λ_F(ℓ) = (3|ℓ)·a_ℓ(G)**, since (12/ℓ) = (3|ℓ).
- T(ℓ²) is semisimple on cusp forms (Hecke operators are normal for the Petersson product). So its minimal polynomial
  μ_b on S_(γ_b,λ_b) has simple roots, each of the form (3|ℓ)·a_ℓ(G).
- Hence μ_b divides P_(K_b)(x) = ∏_G (x − (3|ℓ)·a_ℓ(G)), the characteristic polynomial of (3|ℓ)·T_ℓ on the class at weight
  K_b. Both are monic in ℤ[x] and μ_b is squarefree.
- So **P_(K_b)(T(ℓ²))·G_b = 0 exactly**, and the tower obeys P_(K_b) modulo ℓᵇ:
  Σ_i [xⁱ]P_(K_b)·F_(b+2j+2i) ≡ 0 for all j ≥ 0.
- Yang's isomorphism (I2) is not needed. It would add equality of characteristic polynomials (multiplicity one on the
  half-integral side), which the theorem does not use.

**4.4 Only the ordinary part survives.**
- Split P_(K_b) = P^ord·P^non over ℤ_ℓ (Hensel): the roots that are ℓ-adic units, and those with positive valuation.
  The two factors are coprime modulo ℓ.
- On the P^non part, T(ℓ²)^(deg P^non) ≡ 0 (mod ℓ). So for each m there is j₀ with the P^non component of
  G_b | T(ℓ²)^j ≡ 0 (mod ℓᵐ) for j ≥ j₀.
- So for j ≥ j₀ the tower obeys P^ord modulo ℓ^(min(b, m)).
- The roots of P^ord are (3|ℓ)·a_ℓ(f) for the ordinary newforms f of weight K_b in the class.
- For such f, a_ℓ(f) = α_f + ℓ^(K_b−1)/α_f with α_f the unit root, so a_ℓ(f) ≡ α_f (mod ℓ^(K_b−1)), and K_b − 1 ≥ m once b
  is large.

**4.5 The values move to weight −2.**
- By (I5), the ordinary newforms of weight K_b in the class are the members, at K_b, of the ordinary families through the
  class at weight ℓ − 3. The Atkin–Lehner signs move along a family by (2|ℓ)ˣ and (3|ℓ)ˣ, since 2^((ℓ−1)/2) ≡ (2|ℓ); this
  is how both parities of b land in the families through S_ℓ.
- By Theorem A all of S_ℓ is ordinary, so these are exactly the families through S_ℓ, R_ℓ of them.
- And α_f(K_b) = A_f(K_b) → A_f(−2) = u_f(−2), since K_b → −2.
- So modulo ℓᵐ, for b ≥ b₀ and j ≥ j₀, the sequence F_(b+2j) satisfies the recurrence of ∏_f (x − (3|ℓ)·u_f(−2)) = Q_ℓ(x). ∎

**4.6 What is not used.** Overconvergent forms, Ramsey's lift, Calegari's Lemma 4.6.4, half-integral multiplicity one
(in families or at a single weight), cuspidality of the limit, and Yang's isomorphism. Only divisibility of the minimal
polynomial is needed (§4.3). Each G_b is a cusp form because γ_b > 0.

## 5. Full rank, and 607 [ours, computed]

The minimal polynomial of the step on the tower divides Q_ℓ. It equals Q_ℓ when the tower spans all R_ℓ directions.
- On lane ℓ the levels span at most as many directions as the degree of T_ℓ's minimal polynomial on S_ℓ mod ℓ.
- At ℓ = 607 two Frobenius-conjugate forms of S_607 have the same a_607 on lane 607 (in F₆₀₇), so that degree is 24, not
  25. This matches Boylan–Webb's rank R − 1.
- Theorem A shows the only possible obstruction is of this first kind (coinciding numbers), never a non-ordinary form.
- The chance of such a coincidence is about 1/ℓ per prime, so about 0.7 were expected between 37 and 1297; one was seen.

## 6. Numerical confirmation [computed; files in this folder]

- The constants u_f(−2), from the families, against exact partition numbers: up to 30 digits for 13–31 and 20 digits for
  37–61; 56 digits confirmed on partition numbers, including one of 404,903 digits. (README, Dig 1, 12th–18th.)
- The steps along other primes Q, the laws at 5, 7, 11 as slopes 2, 1, 2, and 3 and 5 colours. (STATEMENT B, C, G.)
- Theorem A's mechanism at 859, and the absence of non-ordinary forms in every class to 300. (§3 Remarks.)

## 7. What remains

1. ~~Yang's trace comparison for n ≡ 5, 7, 13, 17, 19, 23 (mod 24).~~ Not needed: (I2′) gives the eigenvalues at every
   prime (§4.3).
2. ~~Hida's restriction to Atkin–Lehner classes; the continuity modulus.~~ Written out in §8 (R2, R3).
3. ~~(I6) on Γ₀(6) through Γ₁(6).~~ Written out in §8 (R1).
4. **A specialist's reading of the whole. This is the only item left.**

## 8. Routine points, written out [ours]

**(R1) Katz's theorem on Γ₀(6).**
- Γ₁(6) is a fine level (6 ≥ 4), so (I6) holds for mod-ℓ forms on Γ₁(6) with ℓ ≥ 5.
- A form on Γ₀(6) is a form on Γ₁(6) fixed by the diamond operators ⟨d⟩, d ∈ (ℤ/6)^×. If f = A^j·g with A the Hasse
  invariant (level 1) and g on Γ₁(6), then g = f/A^j is again fixed by every ⟨d⟩. So it lies on Γ₀(6), and the
  filtrations of f on Γ₀(6) and on Γ₁(6) agree.
- θ, U, V and T_p (p ∤ 6ℓ) commute with the ⟨d⟩. So every statement of (I6) used in §3 holds on Γ₀(6).
- Taking ⟨d⟩-invariants is exact here, since |(ℤ/6)^×| = 2 is prime to ℓ.

**(R2) Hida's theory on one Atkin–Lehner class.**
- The ordinary Λ-adic Hecke algebra of tame level 6 [H86a, H86b; W88] contains U₂ and U₃, which commute with everything
  else.
- The part new at 2 and 3 is a direct factor: Λ-adic newforms of tame level 6, whose specializations at classical weights
  k ≥ 2 are the ordinary ℓ-stabilized newforms of level 6.
- On a Λ-adic newform new at p ∈ {2, 3}, the U_p-eigenvalue specializes at weight k to −w_p(k)·p^((k−2)/2), and it is
  ℓ-adically continuous in k. As p^((ℓ−1)/2) ≡ (p|ℓ) (mod ℓ), the sign moves by (p|ℓ) per step of ℓ − 1 in k.
- So the class at weight K, with signs set by K, is the set of specializations at K of a fixed set of Λ-adic newforms:
  the families through the class at weight ℓ − 3. Restricting Hida's control isomorphism to that set is restricting it to
  a sum of components, which commutes with specialization.

**(R3) The continuity modulus, with no assumption on ramification.**
- Q_ℓ is not built from individual families' values but from the Λ-adic characteristic polynomial of U_ℓ on those
  components: χ(T; x) = Σ_i c_i(T)·xⁱ with c_i ∈ Λ = ℤ_ℓ[[T]], where weight k corresponds to T = (1 + ℓ)^(k−k₀) − 1.
- If K ≡ K′ (mod (ℓ − 1)ℓ^(t−1)), then (1 + ℓ)^K ≡ (1 + ℓ)^K′ (mod ℓ^t), so the specializations of every c_i agree modulo
  ℓ^t.
- Hence the characteristic polynomials at weights K_b and −2 agree modulo ℓ^b: K_b ≡ −2 (mod (ℓ − 1)ℓ^(b−1)), §4.1.
- No étaleness of the Hecke algebra over Λ is needed, so congruent families (607) and ramified points are covered.

**(R4) Semisimplicity.** For p ≥ 5 the operator T(p²) is normal for the Petersson product on S_(λ+1/2)(1, ν^r), since
these are Hecke operators of index prime to the level of the multiplier. So it is diagonalizable over ℂ, and its minimal
polynomial over ℚ is squarefree. That is what §4.3 uses.

## References

- [AAD] S. Ahlgren, N. Andersen, R. Dicks, The Shimura lift and congruences for modular forms with the eta multiplier,
  Adv. Math. 447 (2024); arXiv:2307.07438.
- [BS] M. Boylan, Swati, Explicit images for the Shimura correspondence, arXiv:2505.01018 (2025).
- [BW] M. Boylan, J. Webb, The partition function modulo prime powers, Trans. AMS 365 (2013).
- [CDHS] W. Y. C. Chen, D. K. Du, Q.-H. Hou, L. H. Sun, Congruences of multipartition functions modulo powers of primes,
  arXiv:1206.6642.
- [DS] P. Deligne, J.-P. Serre, Formes modulaires de poids 1, Ann. Sci. ENS 7 (1974).
- [E92] B. Edixhoven, The weight in Serre's conjectures on modular forms, Invent. Math. 109 (1992).
- [FKO] A. Folsom, Z. Kent, K. Ono, ℓ-adic properties of the partition function, Adv. Math. 229 (2012), with an appendix
  by N. Ramsey.
- [dSG] E. de Shalit, E. Goren, A theta operator on Picard modular forms modulo an inert prime, arXiv:1412.5494 (recalls
  the elliptic drop criterion).
- [GK] E. Ghate, N. Kumar, Control theorems for ordinary 2-adic families of modular forms (TIFR preprint), for the Hida
  citations.
- [H86a] H. Hida, Iwasawa modules attached to congruences of cusp forms, Ann. Sci. ENS (4) 19(2) (1986) 231–273.
- [H88b] H. Hida, On p-adic Hecke algebras for GL₂ over totally real fields, Ann. of Math. 128 (1988) 295–384.
- [H86b] H. Hida, Galois representations into GL₂(ℤ_p[[X]]) attached to ordinary cusp forms, Invent. Math. 85 (1986).
- [J82] N. Jochnowitz, A study of the local components of the Hecke algebra mod l, Trans. AMS 270 (1982).
- [K77] N. Katz, A result on modular forms in characteristic p, LNM 601 (1977).
- [W88] A. Wiles, On ordinary λ-adic representations associated to modular forms, Invent. Math. 94 (1988).
- [Y] Y. Yang, Modular forms of half-integral weights on SL(2, ℤ), Nagoya Math. J. 215 (2014) 1–66; arXiv:1110.1810.
- [Y11] Y. Yang, Congruences of the partition function, IMRN 2011.
