# Spectral expansion units: what the literature says

14 Nov 2026. Question: is the component of P(2) = (1/η)|U(ℓ²) along each family exactly of the family's slope?
Equivalently, are 1/η's own spectral coefficients ℓ-adic units?
Labels: PROVED (theorem in the source), OBSERVED (printed numbers), GUESS (the source's conjecture),
DERIVED (my step from printed numbers, not in the source), MEASURED (computed here).

**In one line:** I found no source that states or explains the pattern. One printed number fits it (Calegari,
ℓ = 5). For j at p = 2 the components sit at the slope + 8 instead. Watson's theorem appears to break the pattern
at ℓ = 7.

## 1. Loeffler, IMRN 2007 (math/0701168)
- PROVED (Thm 6): "U is self-adjoint with respect to this form". The form is the residue pairing
  ⟨α,β⟩ = Res_∞ w*(α)dβ. Hence πᵢ(h) = (⟨h,φᵢ⟩/⟨φᵢ,φᵢ⟩)φᵢ when the eigenspace is a line and ⟨φᵢ,φᵢ⟩ ≠ 0.
- PROVED (Cor. 10, p = 2, r ∈ (5/12, 7/12)): "the finite slope eigenfunctions form an orthonormal basis of the
  space S0(r)".
- OBSERVED (§5: p = 5, h = 1/j, φⱼ = q + … from §8): v₅(cⱼ) = 0, 4, 4, 8, 9, 10, 12, 14, 18, 22 for j = 1…10.
  His only comment: "These numbers appear to be tending 5-adically to zero extremely rapidly".
- DERIVED: the a₅ entries in his §8 give the slopes 1, 4, 5, 8, 9, 10, 13, 14, 19, 20. So
  v(cⱼ) − slope = −1, 0, −1, 0, 0, 0, −1, 0, −1, +2: near the slope, not equal, and not units.
- His only conjecture (Conj. 7) is about convergence. Nothing on units, nonvanishing or the size of coefficients.

## 2. Coefficients through a pairing
- Calegari, AWS 2013 notes, §3.11: F = Σ⟨φᵢ,F⟩/⟨φᵢ,φᵢ⟩·φᵢ. Exercise 3.11.1: "Show that, suitably normalized,
  the invariant pairing ⟨φ, φ⟩ for a finite slope eigenform coincides with the p-adic L-function" Lₚ(1, ad⁰).
  This is a p-adic congruence number. NOT PROVED: §6.6.1 lists it as a project ("Prove Exercise 3.11.1").
- GUESS 5.1.10 (p = 2, N = 1; "based on an embarrassingly small amount of data (for n = 1 to 4)"):
  - It states |⟨φₙ,φₙ⟩| = |n| when φₙ ≡ gⁿ mod 2, and ‖φₙ‖₀²/‖φₙ‖² = |2⁹n²λₙ|.
  - The second formula is my reading of a badly extracted line. It matches his checks for φ₈ (2⁴⁶ = 2⁹·8²·2³¹)
    and φ₆₄.
  - He adds: "the final identity is similar to one guessed by Loeffler when N = 1 and p = 5".
- DERIVED: with a₁ = 1, this gives v⟨φₙ,φₙ⟩ = −(9 + v(n) + slopeₙ). The slope enters each coefficient through the
  denominator. So "valuation = slope + constant" is the same as "⟨F,φₙ⟩ has one valuation for all n". This
  restates the pattern; it does not explain it.
- Level 6, weight −2. None of these three computes expansion coefficients:
  - Bruinier–Ono (1104.1182) build "a theta lift from spaces of weight −2 harmonic weak Maass forms to spaces of
    weight −1/2 vector-valued" forms. Their source form on Γ₀(6) is F = q⁻¹ − 10 − 29q − ….
  - Yang (Nagoya 2014) proves Hecke isomorphisms between η(24τ)ʳf(24τ) and level-6 newforms with fixed
    Atkin–Lehner signs.
  - Calegari, Lemma 4.6.4: the overconvergent Shimura image of S₋₁/₂ lies in S₋₂(Γ₀(6), r)⊗χ₁₂.
- Cohomology: Candelori–Castella (1605.05765) show that a weakly holomorphic form "can be expressed as a linear
  combination fβ and fβ′" (the two p-stabilizations). Boxer–Pilloni (Épiga 2022) prove p-adic Serre duality; not
  read in detail.
- NOT FOUND: any such formula in half-integral weight, or for η⁻¹.

## 3. Cyclic vectors and nonzero projections
- NOT FOUND: a U_p-orbit spanning the slope ≤ α space, or a given form meeting every eigenspace.
- FKO appendix (Ramsey, after Calegari), Thm 7.1, remark (3): "If the rank of Ω is 1, then since Ω is
  Hecke-stable, it is spanned by a single nonzero Hecke eigenform."
- Whether even one dominant component can vanish was left open (GKO 2010, Thm 1.2: "For all but at most one
  choice of α"). It was later settled case by case:
  - Ahlgren–Samart (1602.01036), one CM weight-2 form, p ≡ 3 mod 4: "vp(C(p^{2m+1})) = m". This is your shape:
    a unit component times powers of the eigenvalue.
  - Guerzhoy (2506.07107): "γ ≠ 0 for all but possibly finitely many primes p ∤ N such that b(p) = 0".
  - Tajima (2307.01450), CM form, p inert: "If dimC Sk(Γ0(N)) = 1, then αg is a p-adic unit."

## 4. Partition literature
- Calegari, Ex. 4.6.1 (proved there): "there exists a constant c depending only on p such that"
  Uⁿη⁻¹ ≡ Uⁿ(η⁻¹)_ord mod p^{n−c}. No slope lies in (0, 1), so the non-ordinary parts fade by at least p per U(p²).
- Calegari, §4.7, OBSERVED at p = 5:
  - The U(25) slopes are 2, 7, 9, 15, 19, 22, 27, 29, 36, 39, and "α1 = 4 · 5^2 + 2 · 5^3 + 3 · 5^4 + 5^5 + …".
  - So the first component of P(2) has valuation 2, equal to its slope. This is the only printed instance of your
    pattern that I found. The normalization of φ₁ is not stated.
- Calegari, Ex. 4.7.2, asks to "Explain these congruences in terms of the eigenvalue of Up or Up2 of smallest
  slope". The congruences are Lehner's for j and Watson's 7^{[(d+2)/2]}.
- FKO 2012: rank bounds; "for 13 ≤ ℓ ≤ 23, the dimension is 1". On the rate, FKO say of Boylan–Webb: "Their
  work and numerics suggest that one can generically take bℓ(m) := 2m − 1". Garvan–Shomanov (2504.10824)
  restate Atkin–O'Brien: "an integral constant Kα not divisible by 13", i.e. a unit ordinary eigenvalue.
- NOT FOUND (read): no unit statement in FKO and its appendix, Belmont et al., Garvan–Shomanov, Bruinier–Ono or GKO.
- NOT ACCESSED: Boylan–Webb (Trans. AMS 2013; not on arXiv) and Webb's thesis (PDF refused, 403). Its abstract
  links partition-value ratios to central L-value ratios for ℓ ∈ {5, 7, 11}. Also Ahlgren–Allen–Tang (abstract
  only), Ahlgren–Boylan, Garvan's papers, Yang 2011, Gouvêa–Mazur 1995.
- **DERIVED and MEASURED, ℓ = 7. Check this against your data.**
  - Watson's theorem gives P(2) ≡ 0 mod 7².
  - MEASURED (`watson_check.py`, sympy): for b = 2, 4, 6, 8 the minimum of v₇(p(7^b n + δ_b)) is exactly 2, 3, 4, 5
    (n = 0…5; n = 0, 1 for b ≥ 7). For ℓ = 5 (b ≤ 6) and ℓ = 11 (b ≤ 4) it is exactly b.
  - So the even tower gains one factor of 7 per U(49): the leading slope is 1.
  - Since λ₁^{−m}P(2+2m) → c₁φ₁ and P(2+2m) ≡ 0 mod 7^{m+2}, it follows that c₁φ₁ ≡ 0 mod 7².
  - If there is a single slope-1 family and φ₁ is primitive, then v(c₁) ≥ 2 = slope + 1, and the pattern fails at
    ℓ = 7. The data at ℓ = 5 and 11 fit the pattern.

## 5. j and 1/Δ
- Calegari §4.2, p = 2 (φᵢ with leading coefficient 1): for U₂j − 744, "α1 = 2^11 + 2^12 + . . . and
  α2 = 2^16 + 2^17 + . . .". The slopes are 3 and 7, hence "v(c(2^m)) = 3m + 8". So the components sit at
  slope + 8 and slope + 9, and j's coefficients are NOT units.
- DERIVED from Lehner's exponents by the same limit argument: j's leading component sits at least 3 above the slope
  at p = 3, and at least 1 above it at p = 5.
- 1/j (p = 5): see §1. 1/Δ: NOT FOUND (Belmont et al. give rank bounds only, and need ℓ ≥ 29).

## Sources
- Loeffler: https://arxiv.org/abs/math/0701168 (§5 Thm 6 and table; §6 Cor. 10; §8)
- Calegari: https://swc-math.github.io/aws/2013/2013CalegariLectureNotes.pdf (§3.11, §4.2, §4.6–4.7, §5.1, §6.6)
- FKO with Ramsey's appendix: https://uva.theopenscholar.com/files/ken-ono/files/133_8.pdf
- Garvan–Shomanov https://arxiv.org/abs/2504.10824 · Belmont et al. https://arxiv.org/abs/1510.01202
- GKO https://uva.theopenscholar.com/files/ken-ono/files/125_8.pdf · Ahlgren–Samart https://arxiv.org/abs/1602.01036
- Guerzhoy https://arxiv.org/abs/2506.07107 · Tajima https://arxiv.org/abs/2307.01450
- Candelori–Castella https://arxiv.org/abs/1605.05765 · Bruinier–Ono https://arxiv.org/abs/1104.1182
- Yang https://arxiv.org/abs/1110.1810 · Boxer–Pilloni https://epiga.episciences.org/6112
- Webb's thesis https://scholarcommons.sc.edu/etd/1618
