# The Ramanujan–Watson–Atkin laws (ℓ = 5, 7, 11) as U_ℓ slopes in ℓ-adic families: literature check

Notes written 12 Nov 2026. How the search was done: web search with mathematical terms only, full texts read where they could be opened, two citation lists, and OEIS. No mathematical computation was run. One small text-extraction script was used only to compare two copies of the same PDF. Each number below is marked as **quoted** (copied from a source) or **derived** (my own short arithmetic on quoted material).

Grading key, as the brief asked: **(a)** explicitly stated in the literature; **(b)** partly or implicitly known (I say what and where); **(c)** not found after the searches listed.

---

## Q1. Has anyone described the congruences for 5, 7, 11 as slopes of U_ℓ on overconvergent forms, or as a property of a Coleman family (claim C5)? What did Loeffler 2007 compute?

### Takeaway
**For ℓ = 5 this is stated explicitly.** Calegari's Arizona Winter School 2013 lecture notes (§4.6–4.7) say three things:
- the powers-of-5 congruences come from the **slope-2** eigenvalue of U_25 on overconvergent forms of weight −1/2;
- the overconvergent Shimura lift moves the problem to **weight −2, level Γ0(6), twisted by χ12**;
- the eigencurve then supplies classical forms of weight k ≡ −2 (mod p−1).

**For ℓ = 7 and 11 the same reading is only set as an exercise** (Ex. 4.7.2), so it is implicit, grade (b). Two further sources state the principle in general form:
- Folsom–Kent–Ono 2012, in an appendix by Ramsey whose theorem is credited to Calegari: the Watson/Atkin situation is the case where the module of ordinary (slope-0) overconvergent forms vanishes.
- Raum (arXiv 2105.13170): "U_p-congruences" are linked to "slopes of p-adic modular forms".

**Loeffler 2007** works only in weight 0, tame level 1 (the j-function side). It never mentions partitions, Ramanujan, Watson, or Atkin's matrices.

Grade for each part of C5:
- **The laws are slopes, ℓ = 5:** (a), Calegari §4.7.
- **The laws are slopes, ℓ = 7, 11:** (b). Calegari Ex. 4.7.2 asks the reader to "explain these congruences in terms of the eigenvalue of U_p or U_{p²} of smallest slope". FKO and Ramsey give the general principle.
- **The tower is the weight −2 Shimura partner on level 6 (with the χ12 twist):** (a), Calegari Lemma 4.6.4, for every p > 3. A different weight −2, level-6 partner of p(n) (Bruinier–Ono, theta-lift / Maass sense) is (a) in Ono 2026.
- **Slope 2 for ℓ = 5:** (a). **Slope 1 for ℓ = 7 and slope 2 for ℓ = 11:** (c), not written anywhere I found. They follow directly from Watson's and Atkin's exponents (derived, see Inferences).
- **The formula α_ℓ = (3|ℓ)·u(−2)/ℓ^slope:** (b). For ℓ = 5 Calegari's λ1/25 equals α5 to 13 digits (see Q3), and the χ12 twist is in his Lemma 4.6.4. The formula itself is not written.
- **The specific family members** (ℓ = 7 at weights 40, 82, …, 292; ℓ = 5 at weights 18, 38, …, 258; ℓ = 11 at weights 108, 218, 328), **and Newton interpolation back to weight −2:** (c). Calegari uses classical members of weight ≡ −2 (mod p−1) only inside a proof by contradiction. He gives λ1 as "numerically" computed and does not say how.

### Cited Findings

**Calegari, "Congruences between modular forms", AWS 2013 lecture notes.** Two copies were read in full: [AWS 2013 copy](https://swc-math.github.io/aws/2013/2013CalegariLectureNotes.pdf) and [author's later copy](http://www.math.uchicago.edu/~fcale/papers/AWS.pdf). The later copy has the same §4.7 text; only its acknowledgements gain one name.
- §1.1, p. 3, quoted: "We shall argue that the systematic use of overconvergent modular forms is a direct descendant of the modular equation method." A table on p. 3 pairs "The compactness of the U operator" with "The modular equations method".
- §4.6, Example 4.6.1, p. 45, quoted: "Let η̃⁻¹ denote the projection of η⁻¹ to the ordinary subspace. Then there exists a constant c depending only on p such that Uⁿη⁻¹ ≡ Uⁿη̃⁻¹ mod p^{n−c}." Also: "Note that this has no content if p = 2 or p = 3, so assume that p > 3."
- §4.6.3, p. 46, quoted: "By the overconvergent Shimura correspondence [Ram09], there exists a corresponding overconvergent eigenform of weight −2, also with slope between 0 and 1. A priori, one might expect the level to be Γ0(288). However, after twisting, this form lies in Γ0(6), a fact that requires proof but follows from the underlying symmetry of η."
- Lemma 4.6.4, p. 46, quoted: "The image of the Shimura correspondence from S_{−1/2}(1, r) lies in S†_{−2}(Γ0(6), r) ⊗ χ12, where χ12 is the quadratic character of conductor 12. Moreover, the eigenvalues of the corresponding form of level Γ0(6) of U2 and U3 are given by 2^{−2} and 3^{−2} respectively."
- p. 46, quoted: "From the existence of the eigencurve [CM98] (or by Coleman's results), it follows that there exists a classical form of weight k ≡ −2 mod p − 1 and slope between 0 and 1, as well as the indicated eigenvalues for U2 and U3." This is used, together with Buzzard–Gee and the one-dimensional space S6(Γ0(6)), to rule out slopes in (0, 1).
- §4.7, p. 47, quoted: "A routine computation shows that the corresponding eigenvalues λ1 ≥ λ2 ≥ λ3 ≥ … of U_{p²} on S†_{−1/2}(1, r) have slope 2, 7, 9, 15, 19, 22, 27, 29, 36, 39, …"
- §4.7, quoted: "The fact that φ1 has slope 2 corresponds to the congruence for the partition function modulo powers of 5." Also: "Σ p((25^m n + 1)/24) qⁿ ≡ λ1^m · φ1 + O(5^{7m}), where φ1 is an eigenform for the Hecke operators T_{ℓ²} with ℓ ≠ 5 as well as the operator U25. Here one may numerically compute that λ1 = 4·5² + 2·5³ + 3·5⁴ + 5⁵ + 2·5⁷ + 3·5⁹ + 5¹¹ + 5¹² + 5¹⁴ + …"
- Exercise 4.7.1 (starred as tricky), quoted: "The eigenvalues of weight −2 of level N = 6 with p = 5 which are new at 2 and 3 with w2 = w3 = +1 have slope: v(λn) = v(5^{2n−1} · (3n)!(3n)!/((3n+1)!(3n−1)!) · (6n+2)!(6n−2)!/((2n)!(2n)!))".
- Exercise 4.7.2, quoted: "(2) If 24n ≡ 1 mod 5^c 7^d 11^e, then p(n) ≡ 0 mod 5^c 7^{[(d+2)/2]} 11^e. Explain these congruences in terms of the eigenvalue of U_p or U_{p²} of smallest slope. Compute the slope of the next smallest eigenvalue in each case to give convergence results as above for p = 5."
- Acknowledgements, p. 4, quoted: "thanks also to David Loeffler for some conversations and for making available some of his previous computations."

**Folsom, Kent, Ono, with an appendix by N. Ramsey, "ℓ-adic properties of the partition function", Adv. Math. 229 (2012), DOI 10.1016/j.aim.2011.11.013.** [Full text](https://uva.theopenscholar.com/files/ken-ono/files/133_8.pdf) read.
- Abstract, quoted: "Ramanujan's famous partition congruences modulo powers of 5, 7, and 11 imply that certain sequences of partition generating functions tend ℓ-adically to 0." Also: "For ℓ ∈ {5, 7, 11} we have ⌊(ℓ−1)/12⌋ − ⌊(ℓ²−1)/(24ℓ)⌋ = 0, thereby giving a conceptual explanation of Ramanujan's congruences." And: "Answering questions of Mazur, in the Appendix we give a new general theorem which fits these results within the framework of overconvergent half-integral weight p-adic modular forms. This result, which is based on recent works by N. Ramsey, is due to Frank Calegari."
- Remark (1) after Thm 1.2, quoted: "We do not see how to obtain the full strength of the Ramanujan congruences using only the ideas in the proof of Theorem 1.2. The extra information provided by the modular equations employed by Atkin, Ramanujan, and Watson seem to be necessary for this task."
- Appendix, Thm 7.1, Remark (2), quoted: "If this module vanishes, then this result implies that the coefficients a(ℓ^{2s}n) tend to zero as s → ∞ uniformly in n. This is exactly the situation of Watson's and Atkin's generalization of Ramanujan's congruences." Also, §7.1: "Taking σ = 0, we see that the space of ordinary overconvergent forms is finite-dimensional."
- No slopes beyond 0 are named for 5, 7 or 11, and no limiting constant is given for them. For ℓ = 13 the paper gives the unnormalised analogue (slope 0), quoted: "A13(1, 3, 1) = 6", and "P13(4; 24z) ≡ 45 P13(2; 24z) (mod 13²)".

**Raum, "Relations among Ramanujan-type congruences II" ([arXiv 2105.13170](https://arxiv.org/pdf/2105.13170)).** Full text read.
- p. 2, quoted: "the special case of U_p-congruences, which enjoy a much more arithmetic geometric theory via their connection to slopes of p-adic modular forms."
- p. 3, quoted: "The development set off by Ramanujan bifurcated into a geometric theory of U_p-congruences and a theory of much more combinatorial flavor around Ramanujan-type congruences." Also: "Ramanujan's original congruences lie at the intersection of U_p-congruences and general Ramanujan-type congruences, since M in his cases is a power of ℓ."
- No slopes are computed. On p. 15 he cites Yang (Nagoya Math. J. 215, 2014) for narrowing the level of the eta-character Shimura lift.

**Ono, "The partition function and elliptic curves" ([arXiv 2508.09608v4](https://www.arxiv.org/pdf/2508.09608), 6 Jul 2026).** Full text read.
- Quoted: "we revisit Ramanujan's partition congruences, which follow from a global property of the Uℓ-operator, the Watson–Atkin contraction."
- §5.2 states the inclusions "U5(M⁰5) ⊆ 5M5", "U7²(M⁰7) ⊆ 7M⁰7" and "U11(M⁰11) ⊆ 11M11".
- p. 18, quoted: "This is also the language used in modern treatments of ordinary and finite-slope Up-phenomena". Its footnote: "Calegari's exposition of overconvergent modular forms and Up-operators [8] discusses this."
- It uses the Bruinier–Ono weight −2 meromorphic form on Γ0(6), quoted: "F(τ) := ½ · (E2(τ) − 2E2(2τ) − 3E2(3τ) + 6E2(6τ)) / (η(τ)²η(2τ)²η(3τ)²η(6τ)²) = q⁻¹ − 10 − 29q − ⋯".
- It does not name slopes as valuations of eigenvalues, Coleman families, or any limiting constant.

**Loeffler, "Spectral expansions of overconvergent modular functions", IMRN 2007, no. 16, rnm050 ([arXiv math/0701168](https://arxiv.org/abs/math/0701168)).** The PDF could not be decoded, so I read the [ar5iv HTML version](https://ar5iv.labs.arxiv.org/html/math/0701168) through a summarising fetch. The quotes below are fragments from that summary.
- Setting: N = 1, k = 0, small primes p with X0(p) of genus 0. The main theorem is an instance of a Gouvêa–Mazur conjecture for p = 2, r ∈ (5/12, 7/12).
- p = 5: truncations give eigenvalues of slopes 1, 4, 5. The conjectured valuations i + 2ν5((3i)!/i!) are said to agree with Calegari's conjecture. Constants are given for the spectral expansion of a function h (j = 1..10), and q-expansions of the first 20 eigenfunctions (to O(5¹⁰), per the summary).
- p = 7 and 13: "the pattern is much less clear".
- "Atkin" appears only as "Atkin-Lehner operator U". The U-matrix recurrence is said to be "apparently initially due to Kolberg". **The partition function, Ramanujan, Watson and Atkin's matrices do not appear.**
- Weights other than 0 and Coleman families are mentioned only in passing.

**Bringmann–Guerzhoy–Kane–Mertens–Rolen, "On U_p-congruences for meromorphic modular forms via supersingularity" ([arXiv 2606.14020](https://arxiv.org/html/2606.14020v1), June 2026).** It treats meromorphic forms as overconvergent forms and uses slopes, but makes no mention of p(n) or of the Ramanujan congruences.

### Inferences
- **The (3|ℓ) factor in C5 is the χ12 twist (derived).** χ12(ℓ) = (12|ℓ) = (4|ℓ)(3|ℓ) = (3|ℓ) for odd ℓ. So the factor (3|ℓ) in α_ℓ = (3|ℓ)·u(−2)/ℓ^slope is exactly the "⊗ χ12" twist of Calegari's Lemma 4.6.4: the U_ℓ eigenvalue of a twisted level-6 weight −2 form is χ12(ℓ)·u(−2).
- **The C5 picture for ℓ = 5 is Calegari's §4.6–4.7.** The weight −2 Shimura partner of η⁻¹, at level 6 with the χ12 twist, has slope 2. The ℓ-adic constant per two tower levels is λ1/25 (Q3 shows the digits match). The project's route differs: classical family members at weights ≡ −2 (mod (ℓ−1)ℓ), then Newton interpolation back to weight −2. No source describes that route.
- **Slopes 1 (ℓ = 7) and 2 (ℓ = 11) follow from the classical exponents (derived).** Two tower levels are one application of U_{ℓ²} in weight −1/2, which corresponds to one application of U_ℓ in weight −2.
  - Watson's 7^{⌊b/2⌋+1} gains 7¹ per two levels, so slope 1.
  - Atkin's 11^b gains 11² per two levels, so slope 2.
  - This matches Ono 2026's "U7²(M⁰7) ⊆ 7M⁰7" (one factor of 7 per two steps) and "U11(M⁰11) ⊆ 11M11" (one factor of 11 per step).
- **A consistency check of Exercise 4.7.1 (derived, by hand, from the exercise's formula).** Evaluating the 5-adic valuation for n = 1, 2, 3, 4, 5 gives 2, 7, 9, 15, 19. These equal the first five slopes listed for U_25 on S†_{−1/2}(1, r) in §4.7. So Calegari's half-integral-weight slope list and his weight −2 level-6 slope formula agree, as the Shimura correspondence predicts.

### Gaps
- No source states slope 1 for ℓ = 7 or slope 2 for ℓ = 11 in the overconvergent sense, or a second slope for them. Calegari leaves both to Ex. 4.7.2. Grade: (c) for the explicit statements.
- No source uses classical family members at weights ≡ −2 (mod (ℓ−1)ℓ), or Newton interpolation in the weight, to get α_ℓ. Grade: (c).
- I did not open Ramsey's "The overconvergent Shimura lifting" ([arXiv 0906.3255](https://arxiv.org/pdf/0906.3255), IMRN 2009). The search summary describes a rigid-analytic map from the half-integral weight eigencurve to the integral weight one, with no partition example mentioned. Whether it contains a p(n) example is unchecked.
- **Searched:**
  - Web queries: "Loeffler Spectral expansions…"; "Ramanujan partition congruences slope U_p overconvergent" (extended mode); "overconvergent partition function Ramanujan congruences slope weight −2 level 6 Shimura eigencurve" (extended); "mathoverflow Ramanujan partition congruences p-adic modular forms U_p slope"; "mathoverflow … U_25 overconvergent Calegari"; "Persiflage Calegari blog partition function…"; "site:galoisrepresentations.wordpress.com partition function congruences"; "Arizona Winter School 2013 Calegari project partition function slopes… U_49"; "Ramsey overconvergent Shimura lifting".
  - No MathOverflow thread and no Persiflage post on this topic surfaced.
  - Citation lists: [Semantic Scholar citations of Loeffler 2007](https://api.semanticscholar.org/graph/v1/paper/arXiv:math/0701168/citations?fields=title,year,externalIds&limit=100) gave 23 citing works; none has partitions or Ramanujan congruences in its title. [Semantic Scholar citations of FKO 2012](https://api.semanticscholar.org/graph/v1/paper/DOI:10.1016/j.aim.2011.11.013/citations?fields=title,year,externalIds&limit=200) gave 47 citing works; none mentions slopes, overconvergent forms, eigencurves or Coleman families in its title. Relevant-looking but unread: "A simple proof of the Atkin-O'Brien partition Hecke congruence conjecture for powers of 13" (arXiv 2504.10824; Adv. Math. 2026); "On a Problem of Gordon and Hughes: The 11-adic Properties of Meromorphic Modular Functions on X0(22)" (arXiv 2611.09681); "Explicit constructions for Ramanujan-type congruences" (arXiv 2606.19008).
  - Google Scholar was not used.

---

## Q2. Do Smithline, Buzzard–Calegari, or "Atkin's U_p matrices" for the partition function give the normalized dominant eigenvalue (α5, α7, α11 or their first digits)?

### Takeaway
**Smithline and Buzzard–Calegari: no.** Both work in weight 0, tame level 1, on the j-function side, with no partitions.

**The classical partition-side matrices (Watson, Atkin) give only the first ℓ-adic digit:**
- **ℓ = 7:** first digit 5, from Watson 1938, p. 124 (as restated by Jameson).
- **ℓ = 11:** first digit 9, from Atkin's U11 computations (as presented by Jameson).
- **ℓ = 5:** first digit 4; the full value is Calegari's λ1.

Calegari notes that the Watson method does not readily reach beyond the first eigenvalue.

### Cited Findings

**Smithline, "Bounding slopes of p-adic modular forms" ([arXiv 0705.3614](https://ar5iv.labs.arxiv.org/html/0705.3614), dated 21 June 2001).**
- Treats p ∈ {2, 3, 5, 7, 13}, N = 1, k = 0. It builds a matrix for U from a rational generating function, with sharp results proved only for p = 3.
- "Atkin" appears only in "Atkin's U operator is an endomorphism of…".
- No mention of p(n), Ramanujan, Watson, Lehner, Kolberg or eta-quotients.
- His 2000 thesis and "Compact operators with rational generation" (CRM Proc. 2004) are known here only through secondary descriptions: [Kilford's U5 paper](https://arxiv.org/pdf/math/0606363) and Loeffler. They describe weight-0, level-1 slope bounds (for U5, "a lower bound of x² for the sum of the first x slopes", per a search summary).

**Buzzard–Calegari, "Slopes of overconvergent 2-adic modular forms", Compos. Math. 141 (2005) 591–604.**
- As restated in [Calegari's notes, §5.1.2](https://swc-math.github.io/aws/2013/2013CalegariLectureNotes.pdf), quoted: "the slope of the nth eigenvalue λn of the eigenfunction φn is v(λn) = 1 + 2v((3n)!/(n)!)". This is for N = 1, p = 2, weight 0.
- Calegari §4.2 relates this to Lehner's j-congruences, quoted: "the first few slopes of the spectrum of U2 to be 0, 3, 7, and 13". The setting is the j-function, not partitions.

**Watson's 7-adic matrix: first digit of α7.** [Jameson, Thm 2.1](https://web.math.utk.edu/~mjameso2/research/03.pdf), citing "page 124 of [11]" (Watson 1938, J. reine angew. Math. 179):
- For the generating function of p(7^j n + β7(j)), quoted: "x_{j,i} ≡ 5^{⌊j/2⌋}7^{⌊j/2⌋+1} (mod 7^{⌊j/2⌋+2}) if i = 1, 0 (mod 7^{⌊j/2⌋+2}) otherwise".
- Corollary 2.2, quoted: "p(7^j n + β7(j)) ≡ 5^{⌊j/2⌋}7^{⌊j/2⌋+1} a7(24n + 17) (mod 7^{⌊j/2⌋+2}) if j is odd", with b7(24n + 23) if j is even.

**Atkin's 11-adic matrix: first digit of α11.** [Jameson, Cor. 2.3 and proof](https://web.math.utk.edu/~mjameso2/research/03.pdf):
- Quoted: "p(11^j n + β11(j)) ≡ 4^{3(j−1)/2} 11^j a11(24n + 13) (mod 11^{j+1}) if j is odd, 4^{(3j−4)/2} 11^j b11(24n + 23) (mod 11^{j+1}) if j is even".
- From Atkin's basis on Γ0(11), quoted: "g2(z) | U11 ≡ 4·11 g2(z) (mod 11²), (φ(z)g2(z)) | U11 ≡ 5·11 g2(z) (mod 11²)".

**Calegari §4.8 (Watson approach for 5), quoted.** "it does not seem obvious how one can use this approach to understand the second eigenvalue (and eigenvector) of U25. What one needs to show is that the action on U25 on some natural space is divisible by 5⁷ — yet this is only possible if one can somehow project away from the eigenform φ1 of slope 2."

**Ono 2026 ([arXiv 2508.09608v4](https://www.arxiv.org/pdf/2508.09608)), quoted.** "Atkin's explicit U11-matrix calculation proves U11(M⁰11) ⊆ 11M11". This gives divisibility only, with no eigenvalue constant.

### Inferences
- **First digit of α7 is 5 (derived from Watson's formula).** Going from level j to j+2, with the power of 7 divided out, multiplies the leading constant by 5^{⌊(j+2)/2⌋}/5^{⌊j/2⌋} = 5 (mod 7). This matches the project's first digit of α7, which is 5.
- **First digit of α11 is 9 (derived from Atkin as presented by Jameson).** For either parity, j to j+2 multiplies by 4³ = 64 ≡ 9 (mod 11). Equivalently, the two-step factor is 4·5 = 20 ≡ 9 (mod 11). This matches the first digit of α11, which is 9.
- **First digit of α5 is 4**, the leading coefficient of Calegari's λ1/25. This matches.

### Gaps
- I did not read Smithline's thesis or his CRM paper in full text, nor Buzzard–Calegari 2005 in full text. Secondary descriptions show no partition content.
- I did not read Watson 1938, Atkin 1967, Hirschhorn–Hunt 1981 (powers of 5), or the Paule–Radu / Smoot "localization" papers. Any of them could, in principle, contain a second ℓ-adic digit of the dominant eigenvalue. I found no citation suggesting so.
- **Searched:** "Smithline rational generation compact operators slopes U_p"; "Loeffler spectral expansions"; Kilford U5 ([arXiv math/0606363](https://arxiv.org/pdf/math/0606363)) and Kilford–McMurdy U7 (arXiv 1110.6801) seen only in result lists. Both concern weight k with characters of conductor 25 or 49, not partitions.

---

## Q3. Are the normalized constants α5, α7, α11 computed anywhere beyond the first digit?

### Takeaway
**α5: yes, 13 digits are published, and they match the project's value digit for digit.** Calegari's λ1 in AWS 2013 §4.7 satisfies λ1/25 = α5 to 13 five-adic digits. The project's digits 14–19 (4, 0, 1, 3, 2, 3) are not found anywhere: Calegari's expansion stops at 5¹⁴ with "+ …".

**α7 and α11: only the first digits (5 and 9) are implicit** in Watson and Atkin (see Q2). No further digits were found. OEIS has no hit for either digit string.

### Cited Findings
- **Calegari's λ1, quoted.** "λ1 = 4·5² + 2·5³ + 3·5⁴ + 5⁵ + 2·5⁷ + 3·5⁹ + 5¹¹ + 5¹² + 5¹⁴ + …". It is the U_25 eigenvalue of slope 2 on S†_{−1/2}(1, r), with "Σ p((25^m n + 1)/24) qⁿ ≡ λ1^m · φ1 + O(5^{7m})". Source: [AWS 2013 notes §4.7, p. 47](https://swc-math.github.io/aws/2013/2013CalegariLectureNotes.pdf). The same text appears in the [later copy](http://www.math.uchicago.edu/~fcale/papers/AWS.pdf).
- **OEIS searches: no matches.**
  - α5 digits "4,2,3,1,0,2,0,3,0,1,1,0,1,4": [OEIS search](https://oeis.org/search?q=4%2C2%2C3%2C1%2C0%2C2%2C0%2C3%2C0%2C1%2C1%2C0%2C1%2C4&fmt=short) returned "Sorry, but the terms do not match anything in the table."
  - α7 digits "5,6,1,1,2,1,3,1,0,6,5,3,1": [OEIS search](https://oeis.org/search?q=5%2C6%2C1%2C1%2C2%2C1%2C3%2C1%2C0%2C6%2C5%2C3%2C1&fmt=short) returned no matches.
- **Precedent for ℓ = 13 (slope 0, so no normalization).** [FKO 2012](https://uva.theopenscholar.com/files/ken-ono/files/133_8.pdf) give "A13(1, 3, 1) = 6" and "P13(4; 24z) ≡ 45 P13(2; 24z) (mod 13²)".

### Inferences
- **Digit-by-digit match for α5 (derived).** Dividing Calegari's λ1 by 5² gives the five-adic digits, lowest first: 4, 2, 3, 1, 0, 2, 0, 3, 0, 1, 1, 0, 1 (positions 0–12).
  - The project's α5 begins 4, 2, 3, 1, 0, 2, 0, 3, 0, 1, 1, 0, 1, 4, 0, 1, 3, 2, 3.
  - The first 13 digits agree exactly.
  - The normalizations agree as well: Calegari's U_25 acts once per two tower levels, and his normalizing power is 5² = 5^slope.
- **The sign convention matches C5 (derived).** Calegari's λ1 is the half-integral-weight U_25 eigenvalue. By his Lemma 4.6.4 it equals χ12(5) = (3|5) = −1 times the U5 eigenvalue u(−2) of the untwisted weight −2 level-6 form. So α5 = λ1/25 = (3|5)·u(−2)/5², exactly the C5 formula.
- **The ℓ = 13 precedent (derived).** 45 = 6 + 3·13, so FKO's two-level constant for 13 is known to two 13-adic digits: 6, 3. It is the slope-0 analogue of α_ℓ.

### Gaps
- No source found for α5 digits beyond position 12, or for any α7 or α11 digit beyond the first.
- Calegari's notes do not say how λ1 was computed ("one may numerically compute"). Since the acknowledgements credit Loeffler with computations, Loeffler's unpublished data might hold more digits. Not verified.
- I did not search OEIS for α11 (9, 3, 2, 0): four digits are too few for a meaningful search.
- I did not search for the 5-adic or 7-adic integers written in base 10.

---

## Q4. Is the dependence of slopes and first digits on the weight modulo (ℓ−1)ℓ (Gouvêa–Mazur-type periodicity) documented for these Hecke systems?

### Takeaway
**Not found for these specific systems** (level 6 with Atkin–Lehner signs fixed by the η-multiplier, ℓ = 5, 7, 11). The project's observed patterns are grade (c):
- for ℓ = 5, slope 2 exactly at weights ≡ 18, 2, 6 (mod 20) and slope 3 at ≡ 10, 14 (mod 20);
- for ℓ = 11, the first digit has period 110 in the weight and is symmetric about weight 68;
- for ℓ = 7, (3|7)·a7/7 mod 7 rises by 1 per weight step of 6, and the other members have slopes 7, 14, 19, 24, 31.

Only general theory exists: Hida's theorem for slope 0, the Gouvêa–Mazur conjecture and its refutation and repair, and the recent ghost-conjecture results.

### Cited Findings
- **Hida's theorem, slope-0 analogue of the first-digit periodicity.** [Calegari, Thm 2.2.2 (3)](https://swc-math.github.io/aws/2013/2013CalegariLectureNotes.pdf), quoted: "The minimal polynomial of Up on epMk(Γ, R/pⁿ, 0) only depends on k mod p^{n−1}(p − 1)." Also Thm 2.2.2 (1): the ordinary dimension "depends only on k mod p − 1".
- **Calegari §6.3, quoted.** "φ_{2ⁿ}, which has slope 2^{n+2} − 1, lives in a Coleman family — and if the radius of the family with constant slope is very large (exponential rather than linear as predicted by the Gouvêa–Mazur conjecture) then it will pass through the evil Eisenstein series". This is for N = 1, p = 2.
- **Liu–Truong–Xiao–Zhao, "Slopes of modular forms and geometry of eigencurves"** ([arXiv 2302.07697](https://arxiv.org/abs/2302.07697), v4 Feb 2026). The abstract lists "(d) Gouvea--Mazur conjecture" among corollaries "in the reducible and strongly generic case". According to a search-result summary, not checked in full text, a variant states that for weights above a threshold with v_p(k1 − k2) ≥ m, the slope sequences agree on the localized spaces.
- **Status of the original conjecture,** per a search-result summary of [arXiv 1708.03663](https://arxiv.org/pdf/1708.03663) and related sources, not read in full. The original Gouvêa–Mazur conjecture (constant slope-α multiplicities for weights congruent modulo (p−1)p^{⌈α⌉}) was refuted by Buzzard–Calegari; a modified version is still expected.
- **Loeffler 2007** ([ar5iv](https://ar5iv.labs.arxiv.org/html/math/0701168)) says, quoted, "the complete list of slopes is known only for p=2, tame level 1 and weight 0", and that conjectures for "a general weight, prime and level … appear to be rather inaccessible at present".

### Inferences
- **The ℓ = 5 and ℓ = 11 slope patterns are finer than the original conjecture's guarantee (derived).** The original Gouvêa–Mazur radius for slope 2 would need weights congruent mod (p−1)p², which is mod 100 for p = 5 and mod 1210 for p = 11. The observed period is (ℓ−1)ℓ: 20 and 110.
- **First-digit periodicity matches the ordinary case by analogy only (derived).** Mod-p periodicity of the unit part with period (p−1)p is what Hida gives for slope 0 at n = 2. I found no statement for positive slope.
- **Linear growth of a_p/p for ℓ = 7 is not found.** That (3|7)·a7/7 mod 7 rises linearly with the weight in the slope-1 case was not found in any source read. Literature on reductions of slope-1 crystalline representations might be relevant, but I did not check it (see Gaps).

### Gaps
- Not checked: Gouvêa–Mazur 1992 (Math. Comp.), Wan 1998 (Invent. Math.), Bergdall–Pollack (ghost conjecture), Buzzard's slope conjectures note, Clay's thesis.
- Not checked: papers on reductions of crystalline representations of slope 1 (for example Bhattacharya–Ghate–Rozensztajn). I name these from memory as possibly relevant, unverified.
- The LTXZ hypotheses (p ≥ 11?, reducible generic residual representation) were not confirmed from the full text. Whether they cover the ℓ = 11 system here is unknown.
- **Searched:** "Gouvêa–Mazur conjecture statement slope α weights k ≡ k' mod (p−1)p^n…"; the LTXZ abstract page.

---

## Q5. Is the E = 42a1 mod-49 link (Jameson) explained anywhere in terms of families or weight −2?

### Takeaway
**No, not found (c).** Jameson's paper proves the statement modulo 7^{⌊j/2⌋+2} on the partition values. That is mod 49 for j = 1, 2, and amounts to determining the Hecke eigenvalues mod 7. It uses:
- E7: y² + xy + y = x³ + x² − 4x + 5, which she states has conductor 42;
- a weight-16 Shimura lift on Γ0(288);
- the Θ-operator (a cyclotomic twist), noting that ΘH7 is congruent mod 7 to a form in S16(Γ0(6)).

The paper never mentions weight −2, families, slopes, U_7, or 7³ = 343. The two ingredients of an explanation exist separately: Calegari's weight −2, Γ0(6) ⊗ χ12 eigencurve framework (Q1), and Jameson's level-6 weight-16 residual form. No source joins them.

The same holds for ℓ = 11 (Jameson's E11, conductor 726, at weight 28 on Γ0(6)) and, in Eisenstein form, for ℓ = 5 (Lovejoy–Ono).

### Cited Findings
All from Jameson, "A refinement of Ramanujan's congruences modulo powers of 7 and 11" ([undated PDF on the author's page](https://web.math.utk.edu/~mjameso2/research/03.pdf), Emory address; it thanks a referee; no journal or arXiv number is shown).
- **Abstract, quoted.** "We require the arithmetic properties of two special elliptic curves."
- **The curves, quoted.** "E7 : y² + xy + y = x³ + x² − 4x + 5 and E11 : y² + xy + y = x³ − 14x + 20". Later: "the elliptic curve E7 (which has conductor 42)" and "E11 (which has conductor 726)". No Cremona label is given.
- **Theorem 3.3, quoted.** "If ℓ ≥ 5 is prime, then λ7,a(ℓ) ≡ λ7,b(ℓ) ≡ ℓ(3|ℓ)(1 + ℓ − #E7(F_ℓ)) (mod 7) and λ11,a(ℓ) ≡ λ11,b(ℓ) ≡ ℓ³(3|ℓ)(1 + ℓ − #E11(F_ℓ)) (mod 11)". These λ are the T(ℓ²) eigenvalues of F7 = η¹⁷(24z), G7 = η²³(24z), F11 = η¹³(24z)E8(24z) and G11 = η²³(24z)E8(24z), the mod-ℓ reductions of the normalized tower (Cor. 2.2, 2.3).
- **Proof of Thm 3.3, ℓ = 7, quoted.** The Shimura lifts are "F7(z) = q + 114810q⁵ + ⋯ ∈ S16(Γ0(288))" and "G7(z) = q + 23245050q⁵ + ⋯ ∈ S22(Γ0(288))". Also: "ΘH7 is congruent (mod 7) to a modular form in S16(Γ0(42)). In fact, more is true; it is congruent to a modular form in S16(Γ0(6))". The check uses Sturm's bound and p(n) values "up to 26,161,203".
- **Proof of Thm 3.3, ℓ = 11, quoted.** "Θ³H11 … is congruent to a modular form in S38(Γ0(726)). In fact, more is true; it is congruent to a modular form in S28(Γ0(6))".
- **Theorem 1.1, quoted.** The refined congruences hold "(mod 7^{⌊j/2⌋+2})" and "(mod 11^{j+1})". No higher modulus is claimed.
- **Lovejoy–Ono for 5, as restated by Jameson pp. 1–2.** The coefficient in their mod-5^{j+1} congruence is "(15|ℓ)(1 + ℓ − ℓ²((−24n − 19)|ℓ))" for odd j. Jameson: "their result is equivalent to determining the eigenvalues of the half integer weight Hecke operator T(ℓ²) (mod 5)".
- **No later source joins Jameson to families or weight −2.**
  - [Ono 2026](https://www.arxiv.org/pdf/2508.09608) recasts the congruences through CM points and supersingular points on X0(6). It does not cite Jameson or any conductor-42 curve.
  - [Singh, "The smallest-parts function and elliptic curves" (arXiv 2611.04387, 3 Nov 2026)](https://arxiv.org/abs/2611.04387) follows Ono for spt(n). Its abstract does not mention slopes, families or weight −2.

### Inferences
- **The project's mod-49 multipliers agree with Jameson's theorem (derived).** For a prime Q ≠ 7, Q⁶ ≡ 1 (mod 7), so Q¹⁹ ≡ Q (mod 7). Hence (3|Q)·Q¹⁹·a_Q(E) ≡ (3|Q)·Q·a_Q(E7) (mod 7), which is Jameson's eigenvalue. Two readings:
  - If "lane 49" means mod 49 on the raw p-values (mod 7 on the normalized tower), the observation is exactly Jameson's Theorem 1.1 and 3.3. "Fails mod 343" is consistent with her claiming nothing beyond 7^{⌊j/2⌋+2}.
  - If it means mod 49 on the eigenvalues, the Q¹⁹ refinement is new relative to everything found.
- **The weights line up (derived).** Jameson's weight-16 level-6 form satisfies 16 ≡ −2 (mod 6 = ℓ−1). The project's family weights 40 + 42x satisfy ≡ −2 (mod 42). The weight-28 form for 11 satisfies 28 ≡ −2 (mod 10). So Jameson's forms are the residual members of the families in C5. She identifies them as Θ-twists of E7 and E11 (cyclotomic twists ω¹ and ω³), but never as weight −2 objects.
- **The mod-5 system fits as well (derived).** By Euler's criterion (ℓ|5) ≡ ℓ² (mod 5), and by reciprocity (15|ℓ) = (3|ℓ)(ℓ|5). So Lovejoy–Ono's eigenvalue (15|ℓ)(1 + ℓ) ≡ (3|ℓ)·ℓ²(1 + ℓ) (mod 5). This matches C5's Eisenstein-congruent system a_Q ≡ Q²(1 + Q) (mod 5) with the χ12 = (3|·) twist. So all three residual Hecke systems of C5 are in the literature (Lovejoy–Ono; Jameson). What is missing is their ℓ-adic lift, the family through weight −2, and the constants.

### Gaps
- The label "42a1" was not confirmed. Jameson gives only the equation and conductor 42, and LMFDB was excluded.
- No source found that lifts Jameson's mod-7 eigenvalue statement to mod 49 (for example with Q¹⁹), or that explains E7 or E11 through a Coleman family or weight −2. Grade: (c).
- **Searched:** "partition function congruences powers of 7 elliptic curve conductor 42"; "refinement of Ramanujan's congruences modulo powers of 7 and 11 elliptic curve"; "y^2 + xy + y = x^3 + x^2 - 4x + 5 partition Jameson" (no hits); "Jameson Lovejoy–Ono subprogressions … E_7 E_11" (extended). Read in full: Jameson's paper (a second author URL; the first URL, math.utk.edu/~jameson/research/03.pdf, refused connections), and Ono 2026. Abstract only: Singh 2026.
- Not read: Lovejoy–Ono 2002 in full text ([PDF](https://lovejoy.perso.math.cnrs.fr/5powers.pdf), seen in results only); Yang's "Congruences of the partition function" (IMRN) and "Modular forms of half-integral weight on SL(2,Z)" (Nagoya Math. J. 215, 2014). The latter may document the level-6 Shimura target with Atkin–Lehner signs.
