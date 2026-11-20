# ℓ-adic limiting constants of the partition tower and their identification with weight −2 families: literature check (as of 12 Nov 2026)

How this was read. The PDFs of FKO, Ono 2000, Boylan–Webb, Garvan–Shomanov, AAT, AAD, ABR, the BFKO survey, Guerzhoy 2006 and Calegari's 2013 notes were downloaded and turned into text locally. Quotes of formulas below are copied from that text, with the maths rewritten in plain notation (exponents, fractions and Legendre or Kronecker symbols put back in). Where a glyph was lost in extraction, this is said.

## Q0. Overall verdict: are C1–C4 already in the literature?

### Takeaway
The four claims, in short:
- **C1:** the digits of α_ℓ, the constant with H_{b+2} ≡ α_ℓ·H_b, H_b(n) = p(ℓ^b n + δ_ℓ(b))/ℓ^v.
- **C2:** α_ℓ = (3|ℓ)·u(−2)/ℓ^s, with u the U_ℓ-eigenvalue of the Hida (s = 0) or Coleman (s > 0) family through the level-6 newform of weight ℓ−3.
- **C3:** α_ℓ ≡ (3|ℓ)·a_ℓ(f) (mod ℓ).
- **C4:** the T(Q²) multipliers deep in the tower equal (3|Q)·a_Q(F_{−2}) ℓ-adically.

What the literature has:
- **C1 (digits):** partly explicit.
  - ℓ = 5: a 2013 source gives the first 13 base-5 digits of α5. They match ours.
  - ℓ = 13, 17, 19, 23, 29: the first digit is published for each. The first two digits are published for ℓ = 13 and ℓ = 29.
  - ℓ = 7, 11, 31: no digits were found. A paper on 7 and 11 (Jameson) could not be read.
- **C2 (identification with weight −2):** implicit (b). Calegari's 2013 Arizona Winter School notes set out exactly this frame:
  - the overconvergent Shimura correspondence takes U_{p²}-eigenforms of weight −1/2 to weight −2 eigenforms on Γ0(6) ⊗ χ12, with U2 = 2^−2 and U3 = 3^−2;
  - the eigencurve links these to classical forms of weight ≡ −2 mod (p−1);
  - for p = 5 he computes the eigenvalue λ1.

  Nobody states the formula α_ℓ = (3|ℓ)·u(−2)/ℓ^s, with u read from the Hida or Coleman family through the level-6 newform of weight ℓ−3, or checks it numerically. Not found.
- **C3:** (b), implicit only. It was found stated nowhere, but it follows in a few lines from published results (AAT 2022, Lemma 5.1; AAD 2024, Theorem 1.1, both applied with Q = p = ℓ) once one notes T(ℓ²) ≡ U(ℓ²) (mod ℓ). The ℓ = 13, 17, 19, 23 first digits are published (Ono 2000) without any link to a_ℓ(f).
- **C4:** (b). The depth-1 (mod ℓ) version for Q ≠ ℓ is explicit as AAT Lemma 5.1: f|T(Q²) ≡ (12/Q)·λ_Q f, with λ_Q the T(Q)-eigenvalue of the level-6, weight ℓ−3 newforms. The ℓ-adic version at depth b, compared with a weight −2 family member, is implicit in the Ramsey/Calegari frame but not stated or computed (c).

### Cited Findings
- **C1, ℓ = 5.** Calegari: "Here one may numerically compute that λ1 = 4·5² + 2·5³ + 3·5⁴ + 5⁵ + 2·5⁷ + 3·5⁹ + 5¹¹ + 5¹² + 5¹⁴ + ⋯". Here λ1 is the slope-2 eigenvalue of U_{25} on S†_{−1/2}(1, r) (§4.7, p. 47) — [Calegari, *Congruences between modular forms*, AWS 2013 notes](https://swc-math.github.io/aws/2013/2013CalegariLectureNotes.pdf); same notes at [math.uchicago.edu/~fcale/papers/AWS.pdf](http://www.math.uchicago.edu/~fcale/papers/AWS.pdf)
- **C1, ℓ = 13.** "It turns out that A₁₃(1, 3, 1) = 6 … p(13³n + 1007) ≡ 6 p(13n + 6) (mod 13)" and "A₁₃(2, 4, 2) = 45, which in turn implies that p(13⁴n + 27371) ≡ 45 p(13²n + 162) (mod 13²)" (pp. 2–3) — [Folsom–Kent–Ono preprint](https://uva.theopenscholar.com/files/ken-ono/files/133_8.pdf)
- **C1, first digits for 13, 17, 19, 23.** Ono's §4 lists the "Ramanujan cycles" mod m for 5 ≤ m ≤ 23 (details in Q8) — [Ono, Distribution of the partition function modulo m, arXiv math/0008140](https://arxiv.org/abs/math/0008140)
- **C1, ℓ = 29.** Boylan–Webb give the ℓ = 29 constants mod 29, 29² and 29³ (details in Q2 and Q8) — [Boylan–Webb, TAMS 365 (2013)](https://people.math.sc.edu/boylan/reprints/partition.pdf)
- **C2, frame.**
  - Calegari, §4.6.3: "By the overconvergent Shimura correspondence [Ram09], there exists a corresponding overconvergent eigenform of weight −2…"
  - Calegari, Lemma 4.6.4: "The image of the Shimura correspondence from S_{−1/2}(1, r) lies in S†_{−2}(Γ0(6), r) ⊗ χ12…"

  Source: [Calegari AWS 2013 notes](https://swc-math.github.io/aws/2013/2013CalegariLectureNotes.pdf)
- **C3 and C4, ingredients.** AAT, Lemma 5.1: if g|T(Q) ≡ λ_Q g (mod ℓ) for all g ∈ S^new_{2k−1}(6, −(8/r), −(12/r)), "Then f|T(Q²) ≡ (12/Q) λ_Q f (mod ℓ)" — [Ahlgren–Allen–Tang, arXiv 2112.09481](https://arxiv.org/abs/2112.09481)
- **C3 and C4, Hecke equivariance at every p ≥ 5.** AAD, Theorem 1.1: "S_t(T_{p²}F) = (12/p) T_p S_t(F) for each prime p ≥ 5" — [Ahlgren–Andersen–Dicks, arXiv 2307.07438](https://arxiv.org/abs/2307.07438)

### Inferences
- The (3|ℓ) factor in C2 and C3 is the χ12 = (12/·) twist that appears in all of these sources. For odd ℓ, (12/ℓ) = (3/ℓ).
- C3 and C4 at b = 1 follow from AAT Lemma 5.1 or AAD (1.4) once one notes T(ℓ²) ≡ U(ℓ²) (mod ℓ) in weight (ℓ−2)/2. What is new in our work is the explicit eigenvalue equations, the ℓ-adic depth, and the numerical checks.

### Gaps
- Not readable:
  - Atkin–O'Brien 1967 (AMS PDF: HTTP 403);
  - Atkin 1968 (*Multiplicative congruence properties…*, Proc. LMS);
  - Jameson's 7 and 11 paper (connection reset three times);
  - Calegari's 2011 preprint "A remark on a theorem of Folsom, Kent, and Ono", which Boylan–Webb cite as [10] and which was not found online.
- Not reachable with the tools: MathOverflow (blocked to the fetch tool), Google Scholar, zbMATH (403), Semantic Scholar API (429). OpenAlex was used for "cited by" instead (see Q5).

## Q1. Does FKO (with Ramsey's appendix) identify the limit or eigenvalue, or only prove that it exists? What about the 13 example? Hida, weight −2, Shimura, level 6?

### Takeaway
FKO prove existence only:
- an integer A_ℓ(b1, b2, m) for 5 ≤ ℓ ≤ 31 (Theorem 1.1);
- the Hecke-eigenform property mod ℓ^m (Theorem 1.3);
- rank ≤ ⌊(ℓ−1)/12⌋ (Theorem 1.2).

Ramsey's appendix, with the proof due to Calegari, places this inside overconvergent half-integral weight forms: the limit is the ordinary (slope-0) projection, and rank 1 gives a single Hecke eigenform. The eigenvalue is never identified or computed, beyond the ℓ = 13 numbers 6 (mod 13) and 45 (mod 13²).

None of the following is mentioned: weight −2, Γ0(6), weight ℓ−3, Hida families for p(n). "Hida" appears only about Guerzhoy's work on level-1 functions and in a list of names behind p-adic modular forms. "Shimura's correspondence" appears only as a historical tool. For ℓ = 5, 7, 11, FKO's constant is 0. The normalised positive-slope constants lie outside their frame.

### Cited Findings
- Citation: Folsom, Kent, Ono, "ℓ-adic properties of the partition function", with an appendix by Nick Ramsey, Adv. Math. 229 (2012), no. 3, 1586–1609, DOI 10.1016/j.aim.2011.11.013. Read in the 21-page preprint "Celebrating the life of A. O. L. Atkin" — [preprint PDF](https://uva.theopenscholar.com/files/ken-ono/files/133_8.pdf); [DOI via OpenAlex](https://api.openalex.org/works?search=adic%20properties%20of%20the%20partition%20function&per-page=10&select=id,display_name,publication_year,cited_by_count,doi)
- Abstract: "Modulo any power of a prime ℓ ≥ 5, these sequences of generating functions ℓ-adically converge to linear combinations of at most ⌊(ℓ−1)/12⌋ − ⌊(ℓ²−1)/(24ℓ)⌋ many special q-series." Also: "Answering questions of Mazur, in the Appendix we give a new general theorem which fits these results within the framework of overconvergent half-integral weight p-adic modular forms. This result, which is based on recent works by N. Ramsey, is due to Frank Calegari." — [FKO](https://uva.theopenscholar.com/files/ken-ono/files/133_8.pdf)
- Theorem 1.1: "Suppose that 5 ≤ ℓ ≤ 31 is prime, and that m ≥ 1. If b1 ≡ b2 (mod 2) are integers for which b2 > b1 ≥ b_ℓ(m), then there is an integer A_ℓ(b1, b2, m) such that for every non-negative integer n we have p((ℓ^{b2} n + 1)/24) ≡ A_ℓ(b1, b2, m) · p((ℓ^{b1} n + 1)/24) (mod ℓ^m). If ℓ ∈ {5, 7, 11}, then A_ℓ(b1, b2, m) = 0." Here b_ℓ(m) := 2(⌊(ℓ−1)/12⌋ + 2)m − 3 — [FKO](https://uva.theopenscholar.com/files/ken-ono/files/133_8.pdf)
- Remark after Theorem 1.1: "Boylan and Webb [12] have recently lengthened the range on b in Theorem 1.1. Their work and numerics suggest that one can generically take b_ℓ(m) := 2m − 1." — [FKO](https://uva.theopenscholar.com/files/ken-ono/files/133_8.pdf)
- The ℓ = 13 example (pp. 2–3):
  - "It turns out that A₁₃(1, 3, 1) = 6, and so we have that p(13³n + 1007) ≡ 6 p(13n + 6) (mod 13)."
  - "b1 := 2 and b2 := 4 satisfy the conclusion of Theorem 1.1 with A₁₃(2, 4, 2) = 45, which in turn implies that p(13⁴n + 27371) ≡ 45 p(13²n + 162) (mod 13²)."
  - They check n = 0, 1, 2: both sides are ≡ 99, 89, 20 (mod 13²).

  Source: [FKO](https://uva.theopenscholar.com/files/ken-ono/files/133_8.pdf)
- §6 Examples (pp. 14–15):
  - "P₁₃(2k + 1; z) ≡ 11 · 6^k η(z)^11 (mod 13), P₁₃(2k + 2; z) ≡ 10 · 6^k η(z)^23 (mod 13)", citing "§ 4 of [28]" (Ono 2000);
  - mod 169: "154 P₁₃(4; 24z) ≡ 36q^23 + 17q^47 + 38q^71 + 155q^95 + … (mod 13²)" and "Indeed, it turns out that P₁₃(4; 24z) ≡ 45 P₁₃(2; 24z) (mod 13²)";
  - "Theorem 1.3 is illustrated by the fact that this is a Hecke eigenform modulo 169."

  Note: 154 ≡ 45⁻¹ (mod 169). Source: [FKO](https://uva.theopenscholar.com/files/ken-ono/files/133_8.pdf)
- Theorem 1.3: "If 5 ≤ ℓ ≤ 31 and m ≥ 1, then for b ≥ b_ℓ(m) we have that P_ℓ(b; 24z) (mod ℓ^m) is an eigenform of all of the weight k_ℓ(m) − 1/2 Hecke operators on Γ0(576)."
- Corollary 1.4: "…for every prime c ≥ 5 there is an integer λ_ℓ(m, c) such that for all n coprime to c we have p((ℓ^b n c³ + 1)/24) ≡ λ_ℓ(m, c) p((ℓ^b n c + 1)/24) (mod ℓ^m)."
- No λ_ℓ(m, c) values are given.

  Source: [FKO](https://uva.theopenscholar.com/files/ken-ono/files/133_8.pdf)
- Atkin's conjecture as reported by FKO (p. 5):
  - "Contrary to conventional thinking, he suspected that the P_ℓ(b; 24z) (mod ℓ^m), where the b, m → +∞, converge to Hecke eigenforms for ℓ = 13 and 17."
  - "For 13 (see Theorem 5 of [10]) they prove this observation modulo 13 and 13², and for 17 Atkin claims (see § 6.3 of [7]) to have a proof modulo 17, 17², and 17³."
  - Remark: "Atkin [7] found such congruences modulo 13², 17³, 19², 23⁶ [sic, as printed], 29, and 31."

  Source: [FKO](https://uva.theopenscholar.com/files/ken-ono/files/133_8.pdf)
- Only Hida mention in the main text: "Remark. Guerzhoy [19, 20] has confirmed this speculation for level 1 modular functions using the theory of integer weight p-adic modular forms as developed by Hida, and refined by Wan." — [FKO](https://uva.theopenscholar.com/files/ken-ono/files/133_8.pdf)
- Only Shimura mention: "…modular equations, Hecke operators, Shimura's correspondence, and the Deligne-Serre theory of ℓ-adic Galois representations" (intro, as history) — [FKO](https://uva.theopenscholar.com/files/ken-ono/files/133_8.pdf)
- Appendix, Theorem 7.1 (Ramsey; proof outlined by Calegari): for fixed ℓ, odd k, N and χ, "There exists a module Ω ⊆ Z_ℓ[χ][[q]] of finite rank over Z_ℓ[χ] that is preserved by all half-integral weight Hecke operators", and Σ a(ℓ^{2s} n) q^n ≡ ω_s (mod ℓ^m) with ω_s ∈ Ω for large m, s. Remarks:
  - (1) "We note that, given ℓ, k, N, and χ, the module Ω can be explicitly computed."
  - (2) "If this module vanishes … This is exactly the situation of Watson's and Atkin's generalization of Ramanujan's congruences."
  - (3) "If the rank of Ω is 1, then since Ω is Hecke-stable, it is spanned by a single nonzero Hecke eigenform. This theorem then implies that the images of F under sufficiently many iterates of U(ℓ²) are congruent modulo ℓ^m to multiples of this single form."

  Source: [FKO appendix](https://uva.theopenscholar.com/files/ken-ono/files/133_8.pdf)
- Appendix proof:
  - "e : F ↦ lim_{s→∞} U(ℓ²)^{s!} F … is the projector onto the finite-dimensional ordinary subspace";
  - "Taking σ = 0, we see that the space of ordinary overconvergent forms is finite-dimensional";
  - "U(ℓ²) is invertible on the space of ordinary forms."
  - "Hida" appears only in "the work of Serre, Hida, Katz, Coleman, and Mazur".
  - Acknowledgement: "Calegari outlined the proof of Theorem 7.1".

  Source: [FKO appendix](https://uva.theopenscholar.com/files/ken-ono/files/133_8.pdf)
- The Bruinier–Folsom–Kent–Ono survey repeats the 6 and 45 example and the Calegari/Ramsey remark. Separately, it says the Bruinier–Ono formula uses "the Γ0(6) weight −2 meromorphic modular form F(z)" and a theta lift that "maps spaces of weight −2 harmonic weak Maass forms to spaces of weight −1/2 vector valued harmonic weak Maass forms". That is an archimedean weight −2 ↔ −1/2, level-6 link, not an ℓ-adic one — [BFKO, "Recent work on the partition function"](https://uva.theopenscholar.com/files/ken-ono/files/143_8.pdf)

### Inferences
- FKO's A_ℓ(b, b+2, m) is the U(ℓ²)-eigenvalue mod ℓ^m on the rank-1 ordinary space. Two things show this: P_ℓ(b+2) = P_ℓ(b)|U(ℓ²) in the 24z normalisation, and appendix remark (3). So our α_ℓ for 13 ≤ ℓ ≤ 31 is "the U(ℓ²)-eigenvalue of the ordinary weight −1/2 eigenform" in FKO and Ramsey's terms. They never name it as a value of any other object.
- For ℓ = 5, 7, 11, FKO and Ramsey see only the vanishing of the ordinary part (A = 0). Our normalised α5, α7, α11 are positive-slope eigenvalues and outside FKO.

### Gaps
- The published Adv. Math. version was not compared line by line with the preprint read here. Theorem numbering matches the Boylan–Webb and Garvan–Shomanov citations.
- Where I looked: the full FKO preprint (all 21 pages, searched for Hida, Coleman, Shimura, weight, Γ0(6), slope, eigenvalue); the BFKO survey.

## Q2. Do Boylan–Webb (TAMS 2013) or Belmont–Lee–Musat–Trebat-Leder (Monatsh. 2014) compute such constants, eigenvalues, ranks or stable forms? Any tables?

### Takeaway
- **Boylan–Webb** prove b_ℓ(m) ≤ 2m − 1 (when d_ℓ = 0) and compute the stable bases:
  - m = 1 for 13 ≤ ℓ < 1300; m = 2 for ℓ ≤ 89; m = 3 for ℓ ≤ 29; m = 4 for ℓ = 13.
  - They print a few explicit constants, the useful one being ℓ = 29 mod 29, 29² and 29³. It matches our α29 digits (3, 14) and seems to give a third digit.
  - They also print ranks, the period c_ℓ, and examples of rank 2, 4 and 54.
  - There is no table of constants for 13 ≤ ℓ ≤ 31, and no eigenvalue identification, overconvergent forms (beyond citing Calegari and Ramsey), Hida or Shimura.
- **Belmont et al.** extend the frame to p_r and spt. They give constants only for p_2 (ℓ = 13, constant 10) and spt, quote FKO's 45, and have no tables or families.

### Cited Findings
- Boylan–Webb, "The partition function modulo prime powers", Trans. Amer. Math. Soc. 365 (2013), no. 4, 2169–2206. Abstract: "we prove, with a mild condition on ℓ, that b_ℓ(m) ≤ 2m − 1. Our bound is sharp in all computed cases with ℓ ≥ 29 … This work sheds further light on a question of Mazur posed to Folsom, Kent, and Ono." — [Boylan–Webb reprint](https://people.math.sc.edu/boylan/reprints/partition.pdf)
- Rank bound R_ℓ = ⌊(ℓ+12)/24⌋: "For primes 13 ≤ ℓ ≤ 31, (1.3) gives R_ℓ = 1." Also "Calculations reveal, for all primes 5 ≤ ℓ ≤ 1300 … that d_ℓ = 0." — [Boylan–Webb](https://people.math.sc.edu/boylan/reprints/partition.pdf)
- §6.1: "we calculated bases for the spaces Ω^odd_ℓ(m) and Ω^even_ℓ(m) in the following cases: m = 1, primes 13 ≤ ℓ < 1300; m = 2, primes 13 ≤ ℓ ≤ 89; m = 3, primes 13 ≤ ℓ ≤ 29; m = 4, ℓ = 13." — [Boylan–Webb](https://people.math.sc.edu/boylan/reprints/partition.pdf)
- §6.1, Example 1 (ℓ = 29), with "b₂₉(1) = 1, b₂₉(2) = 3, and b₂₉(3) = 5":
  - "p(29n + 23) ≡ 10 p(29³n + 806) (mod 29)"
  - "p(29³n + 806) ≡ 329 p(29⁵n + 19308) (mod 29²)"
  - "p(29⁵n + 19308) ≡ 14706 p(29⁷n + 13656152078) (mod 29³)"

  Source: [Boylan–Webb](https://people.math.sc.edu/boylan/reprints/partition.pdf)
- Other examples:
  - ℓ = 53: rank 2, "p(53n + 42) ≡ 22 p(53³n + 117861) + 25 p(53⁵n + 331071432) (mod 53)" and a mod 53² relation;
  - ℓ = 89: rank 4; ℓ = 1297: rank 54; ℓ = 607 is the only ℓ ≤ 1297 with r_ℓ(1) < R_ℓ;
  - Corollary 1.3 (periodicity): c₄₁ = 10, c₃₇ = 36, c₁₃₇ = 177423288, where "c_ℓ is the order of a matrix in GL_{r_ℓ(1)}(Z/ℓZ)".

  Source: [Boylan–Webb](https://people.math.sc.edu/boylan/reprints/partition.pdf)
- Remark 2: "Theorem 7.1 of the appendix to [12] describes work of Calegari [10] on how the stability and finiteness results … can be generalized using aspects of the theory of half-integral weight overconvergent p-adic modular forms developed by Ramsey [20], [21]. However, bounds on the stability threshold … require explicit analysis…". Reference [10] is "F. Calegari, A remark on a theorem of Folsom, Kent, and Ono, 2011, preprint." The words overconvergent (as a term), Hida and Shimura do not otherwise occur — [Boylan–Webb](https://people.math.sc.edu/boylan/reprints/partition.pdf)
- Belmont, Lee, Musat, Trebat-Leder, "ℓ-adic properties of partition functions", Monatsh. Math. 173 (2014), no. 1, 1–34 (arXiv posted Oct 2015). Abstract: "We generalize their work to analyze powers p_r of the partition function as well as Andrews's spt-function." — [arXiv 1510.01202](https://arxiv.org/abs/1510.01202)
- Belmont et al. contents (from a reading of the ar5iv HTML, summarised rather than quoted): Corollary 1.2 gives constants C for rank ≤ 1. Example 1 (ℓ = 13, r = 2, m = 1) relates p_2 at levels 4 and 2 with constant 10. Examples 2 and 3 concern spt. There are no tables and no Hecke eigenvalue values, and no Hida, Coleman, overconvergent forms, weight −2, Shimura lift or Γ0(6). They cite FKO's 45 example — [ar5iv 1510.01202](https://ar5iv.labs.arxiv.org/html/1510.01202)

### Inferences
- The ℓ = 29 constants, derived here by hand arithmetic from the printed numbers:
  - 10·3 ≡ 1 (mod 29) and 329·409 = 134561 = 160·841 + 1, so 329·409 ≡ 1 (mod 29²). The first two displays read H_b ≡ c·H_{b+2}, so α29 ≡ 3 (mod 29) and ≡ 409 = 3 + 14·29 (mod 29²). This matches C1's α29 digits (3, 14).
  - The third printed constant 14706 = 3 + 14·29 + 17·29² is ≡ α29 itself mod 29² (14706 ≡ 409), not its inverse (329). So the third display has its sides in the other order from the first two, probably a typo or a change of direction.
  - Read as p(29⁷n + …) ≡ 14706 p(29⁵n + 19308) (mod 29³), it gives a third base-29 digit of α29 equal to 17. This is unverified. Our list has only two digits for α29.
- The two exact relations behind this: H_{b+2}(n) = p(ℓ^{b+2}n + δ) is the U(ℓ²)-image of the level-b series, and Boylan–Webb's c_ℓ is the multiplicative order of the U(ℓ²)-eigenvalue mod ℓ when the rank is 1. They print c_ℓ only for primes of rank 2 or more (37, 41, 137).

### Gaps
- Boylan–Webb computed mod 13⁴ (m = 4) and mod ℓ³ for ℓ ≤ 29, but print constants only for ℓ = 29. The ℓ = 13 constants to 13⁴ and the ℓ = 17, 19, 23 constants to ℓ³ are not printed.
- Where I looked: the full Boylan–Webb reprint (38 pages, text-searched); the ar5iv HTML of Belmont et al. (both parts).

## Q3. Do Atkin and O'Brien (TAMS 126, 1967) give constants relating levels two apart mod powers of 13?

### Takeaway
Yes for existence. Atkin–O'Brien already proved in 1967, for every power 13^α, that there is a constant K_α (13 ∤ K_α) with P(13^{α+2}N) ≡ K_α P(13^α N) (mod 13^α). That is the 13-case of the "two apart" constant decades before FKO. They also proved Hecke-eigenform behaviour mod 13 and 13².

Whether they printed numerical values of K_α (6, 45, …) could not be checked: the paper was not accessible. Their Hecke-multiplier conjecture for all α was proved in 2025 (Garvan–Shomanov, Adv. Math. 2026), again with no values.

### Cited Findings
- Garvan–Shomanov restate the Atkin–O'Brien theorem: "Theorem 1.1. For all α ≥ 1 there exists an integral constant K_α not divisible by 13 such that for all N, P(13^{α+2}N) ≡ K_α P(13^α N) (mod 13^α)." Here P(N) = p(n) if N = 24n − 1, and 0 otherwise — [Garvan–Shomanov, arXiv 2504.10824](https://arxiv.org/abs/2504.10824)
- Atkin–O'Brien conjecture (Conjecture 1.2 there): for α ≥ 1 and primes p ≠ 13, p ≥ 5, "there exists a constant k = k(p, α) such that for all N, P(p²·13^α N) − {k − (−3·13^α N / p) p^{−2}} P(13^α N) + p^{−3} P(13^α N/p²) ≡ 0 (mod 13^α)". Also: "Note: Atkin and O'Brien proved the conjecture for α = 1, 2." — [Garvan–Shomanov](https://arxiv.org/abs/2504.10824)
- Garvan–Shomanov abstract: "In 1967, Atkin and O'Brien conjectured congruences for the partition function involving Hecke operators modulo powers of 13. In this paper, we provide a simple proof of this conjecture." The proof uses Atkin–O'Brien's modular-equation bases and Newman's results. No numerical constants, families or weight −2 — [arXiv 2504.10824](https://arxiv.org/abs/2504.10824)
- Published version: "A simple proof of the Atkin-O'Brien partition Hecke congruence conjecture for powers of 13", Adv. Math. (2026), DOI 10.1016/j.aim.2026.110840 — [OpenAlex citing-works list for FKO](https://api.openalex.org/works?filter=cites:W2049683464&per-page=50&select=id,display_name,publication_year,doi&sort=publication_year:desc)
- FKO: "For 13 (see Theorem 5 of [10]) they prove this observation [Hecke eigenforms] modulo 13 and 13²" — [FKO](https://uva.theopenscholar.com/files/ken-ono/files/133_8.pdf)
- Atkin–O'Brien, "Some properties of p(n) and c(n) modulo powers of 13", Trans. AMS 126 (1967), 442–459. Landing page: [AMS](https://www.ams.org/tran/1967-126-03/S0002-9947-1967-0214540-7/). The PDF returned HTTP 403.
- For j at p = 13, Calegari's notes: "Atkin and O'Brien conjectured something much more precise, namely, they conjecture that Σ c(13^n k) q^k mod 13^n is already a Hecke eigenform. That is, the convergence of U^n_13 j to the ordinary projection is linear." (§2.3) — [Calegari AWS notes](https://swc-math.github.io/aws/2013/2013CalegariLectureNotes.pdf)

### Inferences
- K_α ≡ α13 (mod 13^α) is implied: in our normalisation v = 0 for ℓ = 13, and K_α is the level-α two-step constant. So Atkin–O'Brien own the existence half of C1 for ℓ = 13.

### Gaps
- Not seen: Atkin–O'Brien's own text and any table of K_α or k(p, α); also Atkin's "Congruences for modular forms" (1968, *Computers in Mathematical Research*), which may hold 13-adic tables. MathSciNet review MR0214540 was not accessible.

## Q4. Calegari: does he (blog or elsewhere) identify the limit as a U_ℓ-eigenvalue of an overconvergent eigenform, a Hida family member, or a form of a specific weight, or compute it?

### Takeaway
No blog post on FKO or the partition function was found. The blog's own full-text search finds no post containing "Folsom" or "half-integral", and "partition" hits only two unrelated posts.

His identification is in two places:
1. **FKO appendix (Theorem 7.1).** He outlined it; Ramsey wrote it up. The limit is the ordinary projection under compact U(ℓ²) on overconvergent half-integral weight forms.
2. **His 2013 Arizona Winter School lecture notes, §4.5–4.8.** He:
   - writes the tower as U^n η^{−1} approaching its ordinary projection mod p^{n−c};
   - passes, via Ramsey's overconvergent Shimura correspondence, to an overconvergent eigenform of **weight −2** on **Γ0(6) ⊗ χ12** with U2 = 2^−2 and U3 = 3^−2, and uses the eigencurve to reach classical forms of weight ≡ −2 mod (p−1);
   - for p = 5, computes the U_{25} slopes (2, 7, 9, …) and the slope-2 eigenvalue λ1 5-adically to 13 digits after the factor 5²;
   - sets an exercise on slopes of weight −2, level-6 eigenvalues.

This is the closest prior statement of C2. He does not write α_ℓ = (3|ℓ)u(−2)/ℓ^s, does not use the Hida or Coleman family through the weight ℓ−3 newform, and gives no numbers for ℓ ≥ 7.

### Cited Findings
- Blog searches on galoisrepresentations.org (the blog moved there from galoisrepresentations.wordpress.com):
  - "partition" → only "The stable cohomology of SL(F_p)" (19 Jun 2019) and "Scholze on Torsion, Part III" (22 Jun 2013) — [search "partition"](http://galoisrepresentations.org/?s=partition)
  - "Folsom" → "Nothing Found" — [search "Folsom"](http://galoisrepresentations.org/?s=Folsom)
  - "half-integral" → "Nothing Found" — [search "half-integral"](http://galoisrepresentations.org/?s=half-integral)
  - "Ramanujan" → posts on the "Ramanujan Machine" and similar; none of the excerpts concern p(n) — [search "Ramanujan"](http://galoisrepresentations.org/?s=Ramanujan)
- §4.5.1 Theorem (Ramsey): "Let k be a half-integer. There exists a space of overconvergent forms M†_k(Γ, r) of weight k, and U_{p²} acts compactly on this space." — [Calegari AWS 2013 notes](https://swc-math.github.io/aws/2013/2013CalegariLectureNotes.pdf)
- §4.6.1 Example (η^{−1} "the inverse of Dedekind's eta function" = Σ p((n+1)/24) q^n). Paraphrase; one symbol was lost in extraction: with [g η^{−1}] the projection of η^{−1} to the ordinary subspace, "there exists a constant c depending only on p such that U^n η^{−1} ≡ U^n [g η^{−1}] mod p^{n−c}". Then: "Let e denote Hida's idempotent operator, and let F = U_{p²} η^{−1}(24τ). Then, formally, there is an equality F = e_p F + H + H_{≥1}, where H is a finite sum of generalized eigenforms of slope strictly between 0 and 1 … It suffices to show that H is zero." — [Calegari AWS notes, p. 45–46](https://swc-math.github.io/aws/2013/2013CalegariLectureNotes.pdf)
- §4.6.3, "The Shimura correspondence": "Suppose that H ≠ 0. It follows that there exists, in S†_{−1/2}, an eigenform of slope between 0 and 1. By the overconvergent Shimura correspondence [Ram09], there exists a corresponding overconvergent eigenform of weight −2, also with slope between 0 and 1. A priori, one might expect the level to be Γ0(288). However, after twisting, this form lies in Γ0(6), a fact that requires proof but follows from the underlying symmetry of η." — [Calegari AWS notes, p. 46](https://swc-math.github.io/aws/2013/2013CalegariLectureNotes.pdf)
- §4.6.4 Lemma: "The image of the Shimura correspondence from S_{−1/2}(1, r) lies in S†_{−2}(Γ0(6), r) ⊗ χ12, where χ12 is the quadratic character of conductor 12. Moreover, the eigenvalues of the corresponding form of level Γ0(6) of U2 and U3 are given by 2^−2 and 3^−2 respectively." The exponent's minus sign was lost in extraction; −2 is the reading consistent with weight −2. Then: "From the existence of the eigencurve [CM98] (or by Coleman's results), it follows that there exists a classical form of weight k ≡ −2 mod p − 1 and slope between 0 and 1…". The argument ends in a contradiction via Buzzard–Gee and the unique S_6(Γ0(6)) form "h = q + 4q² − 9q³ + 16q⁴ − 66q⁵ − 36q⁶ + ⋯ has U2 = 4 but U3 = −9". So H = 0: there are no slopes in (0, 1) — [Calegari AWS notes, p. 46–47](https://swc-math.github.io/aws/2013/2013CalegariLectureNotes.pdf)
- §4.7, "An example: congruences for the partition function modulo powers of 5":
  - "A routine computation shows that the corresponding eigenvalues λ1, λ2, λ3, … of U_{p²} on S†_{−1/2}(1, r) have slope 2, 7, 9, 15, 19, 22, 27, 29, 36, 39, …"
  - "The fact that φ1 has slope 2 corresponds to the congruence for the partition function modulo powers of 5. In particular, it follows that: Σ_{n≥0} p((25^m n + 1)/24) q^n ≡ [const]·λ1^m · φ1 + O(5^{7m}), where φ1 is an eigenform for the Hecke operators T_{ℓ²} with ℓ ≠ 5 as well as the operator U_25." (The constant's symbol was lost.)
  - "Here one may numerically compute that λ1 = 4·5² + 2·5³ + 3·5⁴ + 5⁵ + 2·5⁷ + 3·5⁹ + 5¹¹ + 5¹² + 5¹⁴ + ⋯"

  Source: [Calegari AWS notes, p. 47](https://swc-math.github.io/aws/2013/2013CalegariLectureNotes.pdf)
- §4.7.1 Exercise (★): "The eigenvalues of weight −2 of level N = 6 with p = 5 which are new at 2 and 3 with w2 = w3 = +1 have slope: v(λ_n) = v(5^{2n−1} · (3n)!(3n)!/((3n+1)!(3n−1)!) · (6n+2)!(6n−2)!/((2n)!(2n)!))". The formula is transcribed from extracted text — [Calegari AWS notes, p. 47](https://swc-math.github.io/aws/2013/2013CalegariLectureNotes.pdf)
- §4.7.2 Exercise: "(2) If 24n ≡ 1 mod 5^c 7^d 11^e, then p(n) ≡ 0 mod 5^c 7^{[(d+2)/2]} 11^e. Explain these congruences in terms of the eigenvalue of U_p or U_{p²} of smallest slope. Compute the slope of the next smallest eigenvalue in each case to give convergence results as above for p = 5." — [Calegari AWS notes](https://swc-math.github.io/aws/2013/2013CalegariLectureNotes.pdf)
- §4.8, on Watson's method: "it does not seem obvious how one can use this approach to understand the second eigenvalue (and eigenvector) of U_25." — [Calegari AWS notes, p. 49](https://swc-math.github.io/aws/2013/2013CalegariLectureNotes.pdf)
- The notes cite [Ram09] = N. Ramsey, "The overconvergent Shimura lifting", IMRN 2009, no. 2, 193–220, and [Ram08] = "The half-integral weight eigencurve", Algebra Number Theory 2 (2008), 755–808. They do not cite FKO: "Folsom" does not occur — [Calegari AWS notes](https://swc-math.github.io/aws/2013/2013CalegariLectureNotes.pdf); [Ramsey, arXiv 0906.3255](https://arxiv.org/abs/0906.3255)

### Inferences
- **α5 matches.** Direct reading: λ1/5² = 4 + 2·5 + 3·5² + 1·5³ + 0·5⁴ + 2·5⁵ + 0·5⁶ + 3·5⁷ + 0·5⁸ + 1·5⁹ + 1·5¹⁰ + 0·5¹¹ + 1·5¹² + ⋯. Its digits are 4, 2, 3, 1, 0, 2, 0, 3, 0, 1, 1, 0, 1, identical to the first 13 digits of our α5. Our digits 14–19 (4, 0, 1, 3, 2, 3) are not in the notes. "Slope 2" matches our s = 2 for ℓ = 5.
- **Sign rules agree.** Calegari's χ12 twist and U2 = 2^−2, U3 = 3^−2 at weight −2 are consistent with C3's rule (a2 > 0 iff (2|ℓ) = 1, a3 > 0 iff (3|ℓ) = 1) at weight ℓ−3. In a family the ℓ-adic continuation of a_p = ±p^{k/2−1} changes sign by p^{(ℓ−1)/2} ≡ (p|ℓ) (mod ℓ) between k = ℓ−3 and k = −2 (Euler's criterion). So the sign (p|ℓ) at weight ℓ−3 becomes + at weight −2. This is reasoning, not a computation.
- **Status of C2.** Our C2 is the explicit, quantitative form of Calegari's §4.6–4.7 frame: one equation for the eigenvalue at weight −2, computed through the family through the classical weight ℓ−3 newform, for every 5 ≤ ℓ ≤ 31. The frame is published; the equation, the interpolation and the digit checks were not found.

### Gaps
- The compiled blog PDF ([math.uchicago.edu/~fcale/papers/blog.pdf](http://www.math.uchicago.edu/~fcale/papers/blog.pdf)) exceeded the 10 MB fetch limit and was not searched. The site search above covers the live blog.
- Calegari's 2011 preprint "A remark on a theorem of Folsom, Kent, and Ono" (cited by Boylan–Webb) was not found online.
- MathOverflow could not be searched (blocked to the fetch tool).

## Q5. Is there any paper or thesis (2012–2026) that relates the FKO limit, or p(n) mod prime powers, to Hida or Coleman families, the Shimura lift from weight −1/2 to weight −2, or level-6 newforms of weight ℓ−3?

### Takeaway
- **Level-6 newforms of weight ℓ−3:** firmly in the literature, but mod ℓ only and for congruences of Atkin type, not for the FKO constants. Sources: Yang 2014 (Shimura lifts of level-1 η-multiplier forms are twisted level-6 newforms); AAT 2022 (the lift of P_ℓ(1) mod ℓ lies in S^new_{ℓ−3}(6, −(8/−ℓ), −(12/−ℓ)) ⊗ (12/·)); AAD 2024 (level 6N with Atkin–Lehner signs, Hecke-equivariant at all p ≥ 5); ABR 2023.
- **Weight −2 with a p-adic family:** only in Calegari's 2013 notes, built on Ramsey's eigencurve and overconvergent Shimura lift (Q4).
- **The j-function analog:** Guerzhoy 2006 and 2010 use Hida's control theorem for the U_p-limit of j at weight 0, level 1.
- Not found: any paper or thesis that names a Hida or Coleman family for p(n), evaluates one at weight −2, or ties FKO's A_ℓ to a level-6 eigenvalue.

### Cited Findings
- AAT (2.6) and (2.8), from Ahlgren, Allen, Tang, "Congruences like Atkin's for the partition function", Trans. AMS Ser. B 9 (2022) 1044–1064:
  - "f_ℓ ∈ S_{(ℓ−2)/2}(1, ν_η^{−ℓ}, Z) with f_ℓ ≡ Σ p((ℓn+1)/24) q^{n/24} (mod ℓ)" (2.6);
  - "each Shimura lift of f_ℓ is in the space S^new_{ℓ−3}(6, −(8/−ℓ), −(12/−ℓ)) ⊗ (12/•)" (2.8).

  Source: [arXiv 2112.09481](https://arxiv.org/abs/2112.09481); [journal PDF](https://www.ams.org/journals/btran/2022-09-33/S2330-0000-2022-00128-X/S2330-0000-2022-00128-X.pdf)
- AAT notation: "We use the notation S^new_k(6, ε2, ε3) to denote the new subspace of modular forms of integral weight k on Γ0(6) with eigenvalues ε2 and ε3 under the Atkin-Lehner involutions W2 and W3." Theorem 1.4 assumes "there is no congruence modulo any prime above ℓ between distinct newforms in S^new_{ℓ−3}(6, −(8/−ℓ), −(12/−ℓ))" — [AAT](https://arxiv.org/abs/2112.09481)
- AAT text search: "Folsom" 0 hits; "Hida" only as a bibliography entry ([Hid00]); no "weight −2", no U(ℓ²) eigenvalue, no ℓ-adic limits. Theorem 2.1(3) and Lemma 3.1 discuss whether ι_ℓ(a_ℓ) is an ℓ-adic unit (the shape of ρ_f at ℓ) — [AAT](https://arxiv.org/abs/2112.09481)
- AAD (Adv. Math. 447 (2024) 109655) abstract: "We prove that the lift of a cusp form of weight λ + 1/2 and level N has weight 2λ and level 6N, and is new at the primes 2 and 3 with specified Atkin-Lehner eigenvalues." Theorem 1.1 adds "S_t(T_{p²}F) = (12/p) T_p S_t(F) for each prime p ≥ 5" (1.4). No FKO or U(ℓ) in the text — [arXiv 2307.07438](https://arxiv.org/abs/2307.07438)
- ABR (Amer. J. Math. 145 (2023) 1509–1548): "By [Yan14, Thm 1], there is a newform F ∈ S^new_{2k−1}(6) with the same Hecke eigenvalues as Sh_t(f) ⊗ (12/•) at all primes Q ≥ 5". [Yan14] = Yifan Yang, "Modular forms of half-integral weights on SL(2, Z)", Nagoya Math. J. 215 (2014), 1–66 — [arXiv 2006.07645](https://arxiv.org/abs/2006.07645)
- Ramsey's eigencurve papers:
  - "The half-integral weight eigencurve" (ANT 2008) builds an eigencurve for finite-slope overconvergent half-integral weight eigenforms;
  - "The overconvergent Shimura lifting" (IMRN 2009) builds a rigid-analytic map from the half-integral weight cuspidal eigencurve to the integral weight one, interpolating the classical Shimura lift.

  Neither mentions the partition function in its abstract — [arXiv 0906.3249](https://arxiv.org/abs/0906.3249); [arXiv 0906.3255](https://arxiv.org/abs/0906.3255)
- Guerzhoy, "On a conjecture of Atkin for the primes 13, 17, 19, and 23", Math. Res. Lett. 13 (2006):
  - "We show how a combination of Serre's theory of p-adic modular forms and Hida's Control Theorem explains the phenomenon."
  - Theorem 1: for m ≥ M "there exists a simultaneous Hecke eigenform f ∈ S_k(Γ0(p); Z_p) … such that f ≡ (F − c(F; 0))|U^m mod p^α".
  - Proof: "Hida's Control Theorem … implies that rank S^ord_k(Γ0(p); Z_p) does not depend on k ≥ 2 if k is constant modulo p − 1 … Thus rank … = 1".

  Source: [Guerzhoy MRL 2006 PDF](http://www.math.hawaii.edu/~pavel/publications/atkin.pdf)
- Guerzhoy 2010 (Crelle 639, 121–130): per the publisher abstract seen in search, it proves for primes p ≤ 23 that coefficients of U_m(j − 744) behave like a Hecke eigenform modulo growing powers of p, with special attention to p = 5, 7, 11 — [De Gruyter](https://www.degruyterbrill.com/document/doi/10.1515/crelle.2010.013/html). The page itself refused the fetch (HTTP 405); the abstract was seen only in the search snippet.
- Newer work on U_p-congruences: Bringmann, Guerzhoy, Kane, Mertens, Rolen, "On U_p-congruences for meromorphic modular forms with supersingularity" (arXiv June 2026). Integral weight, uses overconvergent forms; no partition function or half-integral weight in the abstract — [arXiv 2606.14020](https://arxiv.org/abs/2606.14020)
- OpenAlex "cited by" for FKO (W2049683464): 34 citing works. Titles include Garvan–Shomanov (Adv. Math. 2026), "ℓ-Adic properties and congruences of ℓ-regular partition functions" (Res. Number Theory 2024), "Newman's conjecture for the partition function modulo integers with at least two distinct prime divisors" (Adv. Math. 2025), Belmont et al., and "Ramanujan's Unpublished Manuscript on the Partition and Tau Functions" (2012). None has Hida, Coleman, eigencurve, Shimura or weight −2 in its title. The list is incomplete: it omits Boylan–Webb, which cites FKO — [OpenAlex citing works](https://api.openalex.org/works?filter=cites:W2049683464&per-page=50&select=id,display_name,publication_year,doi&sort=publication_year:desc)
- Other recent level-6 work: Ono, "The partition function and elliptic curves" (arXiv 2508.09608, v4 July 2026) recasts the Bruinier–Ono formula as traces at CM points on X0(6), reduced to the supersingular locus. It touches 5, 7, 11 and Ramanujan's congruences, but no ℓ-adic limits or families — [arXiv 2508.09608](https://arxiv.org/abs/2508.09608)
- Bringmann–Craig–Ono, "Ramanujan's partition generating functions modulo ℓ" (Ramanujan J. 68 (2025), art. 85): P_ℓ mod ℓ ≡ c_ℓ·T_ℓ(q)/(q^ℓ; q^ℓ)_∞ with Hecke traces of weight ℓ−1 forms. For ℓ = 13: "c_13 ≡ 6 mod 13". This is a normalising constant, not the two-step constant. No U(ℓ²) eigenvalues, Shimura lift or families — [PMC12586220](https://pmc.ncbi.nlm.nih.gov/articles/PMC12586220/)

### Inferences
- Atkin–Lehner signs. As Kronecker symbols, (8/−ℓ) = (2/ℓ) and (12/−ℓ) = (3/ℓ). So AAT's space is S^new_{ℓ−3}(6, ε2 = −(2/ℓ), ε3 = −(3/ℓ)). With a_p = −p^{k/2−1}ε_p for p ∥ 6, a2 > 0 iff (2|ℓ) = 1 and a3 > 0 iff (3|ℓ) = 1. That is exactly how C2 and C3 pick f. The literature therefore fixes "the" form f in C3 (AAT 2022, building on Yang 2014).
- Precedents and novelty for C2:
  - Guerzhoy is the integral-weight precedent (j at weight 0, level 1, Hida control);
  - Calegari 2013 is the half-integral weight precedent (weight −1/2 → −2, level 6, via Ramsey).

  What is new in ours: specialising the family through the weight ℓ−3 newform at k = −2 and matching it to the measured α_ℓ.

### Gaps
- No full-text search engine for arXiv or theses was available. Searches run (all found nothing that combines the partition tower with Hida or Coleman families or weight −2 eigenvalues):
  - "partition function ℓ-adic limit U_ℓ eigenvalue Hida family weight −2";
  - "partition function p-adic family weight -2 Γ0(6) Hida OR Coleman OR eigencurve";
  - "partition function Hida family / Coleman family … weight −1/2";
  - "half-integral weight eigencurve partition function Shimura lift overconvergent";
  - "overconvergent half-integral weight slopes U_{p^2} partition function 1/eta 5-adic eigenvalue weight -2 level 6";
  - "Ramsey overconvergent half-integral weight partition function … Atkin congruence Hecke operators";
  - blog site searches (Q4).
- No theses were found. Theses of Webb, Belmont et al. or Ramsey's students were not checked individually.

## Q6. Is C3 (first digit α_ℓ ≡ (3|ℓ)·a_ℓ(f) mod ℓ, f the level-6 newform of weight ℓ−3) stated anywhere?

### Takeaway
Not found as a statement. Every ingredient is published, and C3 follows in a few lines:
1. P_ℓ(3) = P_ℓ(1)|U(ℓ²).
2. T(ℓ²) ≡ U(ℓ²) (mod ℓ) in weight (ℓ−2)/2.
3. AAT Lemma 5.1 or AAD (1.4), Hecke equivariance with the (12/ℓ) twist, applied at Q = p = ℓ.
4. AAT (2.8) puts the lift in S^new_{ℓ−3}(6, …).

The mod-ℓ first digits themselves (6, 6, 10, 5 for 13, 17, 19, 23) are published in Ono 2000 without any link to a_ℓ(f). The j-function analog of C3 is explicit (Serre; Guerzhoy).

### Cited Findings
- AAT Lemma 5.1, verbatim apart from symbol layout: "Let ℓ ≥ 5 be prime, and suppose that f ∈ S_k(1, ν_η^r, Z) where (r, 24) = 1 and k ∈ ½Z∖Z. Suppose that Q ≥ 5 is a prime and that λ_Q is an integer with g|T(Q) ≡ λ_Q g (mod ℓ) for all g ∈ S^new_{2k−1}(6, −(8/r), −(12/r), Z). Then f|T(Q²) ≡ (12/Q) λ_Q f (mod ℓ)." AAT apply it with Q ≡ ±1, −2 (mod ℓ) and Q ≠ ℓ — [AAT](https://arxiv.org/abs/2112.09481)
- AAT's Hecke operator (2.4): "f|T(Q²) = Σ (a(Q²n) + Q^{k−3/2} (−1/Q)^{k−1/2} (12n/Q) a(n) + Q^{2k−2} a(n/Q²)) q^{n/24}". For Q = ℓ and k = (ℓ−2)/2, the last two terms vanish mod ℓ — [AAT](https://arxiv.org/abs/2112.09481)
- AAD Theorem 1.1 (1.4), Hecke equivariance "for each prime p ≥ 5" (includes p = ℓ when ℓ ∤ N) — [AAD](https://arxiv.org/abs/2307.07438)
- Ono 2000, §4 "Examples", where F(m, k; z) is the generating function of p((m^k n + 1)/24):
  - "Using Sturm's theorem again, one readily verifies that η^11(24z)|U(13) ≡ 8 η^23(24z) (mod 13), η^23(24z)|U(13) ≡ 4 η^11(24z) (mod 13)".
  - "F(13, 2k + 1; z) ≡ 11·6^k η^11(24z) (mod 13), F(13, 2k + 2; z) ≡ 10·6^k η^23(24z) (mod 13). These two congruences appear in Ramanujan's unpublished manuscript on τ(n) and p(n)".

  Source: [Ono, Annals 151 (2000), arXiv math/0008140](https://arxiv.org/abs/math/0008140)
- j-analog (Guerzhoy 2006, citing Serre): "the congruences (5) c(p) F_{p−1} ≡ (j − 744)|U mod p with c(p) ≢ 0 mod p are proved by Serre in [13, 6.16]". It lists b(p) = a_p(F_{p−1}): "b(13) = −577738, b(17) = 1646527986, b(19) = 1487499860, b(23) = −73845437470344", to check p-ordinarity — [Guerzhoy 2006](http://www.math.hawaii.edu/~pavel/publications/atkin.pdf)

### Inferences
- Deriving C3:
  - Apply Lemma 5.1 or (1.4) with Q = ℓ to f_ℓ ≡ P_ℓ(1). Then f_ℓ|U(ℓ²) ≡ f_ℓ|T(ℓ²) ≡ (12/ℓ)·a_ℓ(F)·f_ℓ (mod ℓ), with F the untwisted newform in S^new_{ℓ−3}(6, −(2/ℓ), −(3/ℓ)), assuming that space is one-dimensional or F is unique mod ℓ.
  - (12/ℓ) = (3/ℓ), so α_ℓ ≡ (3|ℓ)·a_ℓ(F) (mod ℓ), which is C3.
  - It therefore counts as "implicit in published results", not "stated". I found no paper that writes it or checks it for 13 ≤ ℓ ≤ 31. Ono 2000 computes the first digits by Sturm-bound checks.
- The j-analog (U-eigenvalue mod p = a_p of the weight p−1 cusp form) is explicit and classical. The half-integral weight version with the χ12 twist and level 6 is not.

### Gaps
- I did not read Yang 2014 (Nagoya) or Garvan's and Yang's papers on Hecke-invariant subspaces (cited in FKO §5); one of them might state T(ℓ²)-eigenvalues on η^r·E-type forms mod ℓ. Not checked: Atkin 1968, where FKO say Atkin found "accidental" eigenvalues.

## Q7. Is C4 (deep in the tower, T(Q²) multipliers λ(Q) ≡ (3|Q)·a_Q(F_{−2}) to about b digits) stated anywhere?

### Takeaway
- Existence of the ℓ-adic Hecke multipliers is published: FKO Theorem 1.3 and Corollary 1.4 (λ_ℓ(m, c), 5 ≤ ℓ ≤ 31); Atkin–O'Brien (13, α = 1, 2); Garvan–Shomanov (13, all α).
- Their identification mod ℓ with (12/Q)·a_Q of the level-6 weight ℓ−3 newform is AAT Lemma 5.1.
- Their identification at depth ℓ^b with a_Q of the weight −2 family member is implicit in the Hecke-equivariant overconvergent Shimura lift (Ramsey 2009) as used by Calegari, but no statement or computation was found. Calegari only notes, for p = 5, that φ1 "is an eigenform for the Hecke operators T_{ℓ²} with ℓ ≠ 5".

### Cited Findings
- FKO Corollary 1.4 (existence of λ_ℓ(m, c) mod ℓ^m, no values); remark "Atkin [7] found such congruences modulo 13², 17³, 19², 23⁶ [sic], 29, and 31" — [FKO](https://uva.theopenscholar.com/files/ken-ono/files/133_8.pdf)
- Garvan–Shomanov Corollary 1.3: for α ≥ 1 and primes p ≥ 5, p ≠ 13, "there is a constant k = k(p, α) such that for all n coprime to p, p((13^α p³ n + 1)/24) ≡ k p((13^α p n + 1)/24) (mod 13^α)". No values — [arXiv 2504.10824](https://arxiv.org/abs/2504.10824)
- AAT Lemma 5.1 (mod ℓ, quoted in Q6) — [AAT](https://arxiv.org/abs/2112.09481)
- Calegari §4.7: "φ1 is an eigenform for the Hecke operators T_{ℓ²} with ℓ ≠ 5 as well as the operator U_25" — [Calegari AWS notes](https://swc-math.github.io/aws/2013/2013CalegariLectureNotes.pdf)
- Ramsey's overconvergent Shimura lifting interpolates the classical Shimura lift and is Hecke-equivariant, mapping half-integral weight eigensystems to integral weight ones — [arXiv 0906.3255](https://arxiv.org/abs/0906.3255)

### Inferences
- C4's mod-ℓ case is AAT Lemma 5.1. The deep version follows from Ramsey's lifting applied to the ordinary (or slope-s) eigenform, as in Calegari §4.6, together with C2's choice of family. No source states the result, computes λ(Q) ℓ-adically, or compares λ(Q) with a_Q of a weight −2 family member.

### Gaps
- Not seen: any table of λ_ℓ(m, Q) values (Atkin 1968 may have some mod 13², 17³ and so on, inaccessible); and the Weaver (2001) and Johansson (2012) data on "accidental eigenvalues" (≡ 0 mod ℓ cases only).

## Q8. C1 digit check: which sources give digits of α_ℓ, and for which primes?

### Takeaway
| ℓ | Digits in print (lowest first) | Source | Our list |
|---|---|---|---|
| 5 | 4, 2, 3, 1, 0, 2, 0, 3, 0, 1, 1, 0, 1 (13 digits, as λ1/5²) | Calegari 2013 | first 13 of 19 match; digits 14–19 not found |
| 7 | none | — | — |
| 11 | none | — | — |
| 13 | 6 | Ono 2000; credited to Ramanujan's unpublished manuscript | match |
| 13 | 45 = 6 + 3·13 (mod 13²) | FKO | match |
| 17 | 6 | Ono 2000 | match |
| 19 | 10 | Ono 2000 | match |
| 23 | 5 | Ono 2000 | match |
| 29 | 3 (mod 29), 3 + 14·29 (mod 29²), from the printed inverses 10 and 329 | Boylan–Webb | match |
| 29 | possibly 17 as a third digit, from 14706 if its display is read with sides swapped | Boylan–Webb | not in our list |
| 31 | none | — | — |

Atkin–O'Brien 1967 proved that the 13-constants K_α exist mod 13^α for all α. Whether they print values was not checked.

### Cited Findings
- ℓ = 5: λ1 = 4·5² + 2·5³ + 3·5⁴ + 5⁵ + 2·5⁷ + 3·5⁹ + 5¹¹ + 5¹² + 5¹⁴ + ⋯ — [Calegari AWS notes §4.7](https://swc-math.github.io/aws/2013/2013CalegariLectureNotes.pdf)
- ℓ = 13:
  - "F(13, 2k + 1; z) ≡ 11·6^k η^11(24z) (mod 13)", credited to "Ramanujan's unpublished manuscript on τ(n) and p(n)" — [Ono 2000](https://arxiv.org/abs/math/0008140)
  - A₁₃(1, 3, 1) = 6 and A₁₃(2, 4, 2) = 45 — [FKO](https://uva.theopenscholar.com/files/ken-ono/files/133_8.pdf)
- ℓ = 17: "η^7(24z)E4(24z)|U(17) ≡ 7 η^23(24z)E4(24z) (mod 17), η^23(24z)E4(24z)|U(17) ≡ 13 η^7(24z)E4(24z) (mod 17)". The corollary for even levels carries "15·6^k" (mod 17) — [Ono 2000](https://arxiv.org/abs/math/0008140)
- ℓ = 19: "F(19, 2k + 1; z) ≡ 5·10^k η^5(24z)E6(24z) (mod 19), F(19, 2k + 2; z) ≡ 11·10^k η^23(24z)E6(24z) (mod 19)" — [Ono 2000](https://arxiv.org/abs/math/0008140)
- ℓ = 23: "F(23, 2k + 1; z) ≡ 5^k η(24z)E4(24z)E6(24z) (mod 23), F(23, 2k + 2; z) ≡ 5^{k+1} η^23(24z)E4(24z)E6(24z) (mod 23)". Ono's §4 covers "5 ≤ m ≤ 23" only — [Ono 2000](https://arxiv.org/abs/math/0008140)
- ℓ = 29: the three congruences in Q2 — [Boylan–Webb §6.1](https://people.math.sc.edu/boylan/reprints/partition.pdf)

### Inferences
- The two-step factor equals the product of the two one-step U(ℓ) factors: 8·4 = 32 ≡ 6 (mod 13) and 7·13 = 91 ≡ 6 (mod 17). This is derived arithmetic, consistent with the 6^k in Ono's corollaries.
- The ℓ = 5 match to 13 digits is the strongest prior overlap with our results. Calegari's λ1 is defined as the U_{25}-eigenvalue on S†_{−1/2}(1, r), which is the same object as our normalised two-step constant times 5².

### Gaps
- ℓ = 7 and 11: Jameson's "A refinement of Ramanujan's congruences modulo powers of 7 and 11" ([PDF](https://www.math.utk.edu/~jameson/research/03.pdf)) could not be retrieved (connection reset three times). Per the search snippet, it gives Lovejoy–Ono-type refinements using "the elliptic curves E_7 and E_11" to describe the relevant Hecke eigenvalues mod 7 and 11. These may be the classical first digits of α7 and α11, but this is unverified.
- Lovejoy–Ono (Crelle 542, 2002) on powers of 5 was not read.
- ℓ = 31: nothing found. Ono 2000 stops at 23, and Boylan–Webb print only ℓ = 29 among 13–31.
- More 13-adic digits may be in Atkin–O'Brien 1967 or Atkin 1968; both inaccessible.
