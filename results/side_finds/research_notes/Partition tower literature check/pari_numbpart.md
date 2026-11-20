# PARI/GP numbpart for large n: known bugs, the 2.15.4 code, cypari2/WSL reports, and how reliable sympy and FLINT/Arb are

Search done 12 Nov 2026. Primary sources read: the PARI CHANGES files (2.2 to 2.19.0, from the 2.19.0 release tarball), the whole PARI bug tracker list (2,302 reports, #1 to #2729), the PARI source of `numbpart` in 2.15.4, 2.17.4 and 2.19.0 (release tarballs), the pari-dev and pari-users messages that web search found, the GitHub trackers of cypari2, sympy, Sage and WSL, the sympy 1.12 and 1.14 source, and Johansson's paper, page, blog and the Arb/FLINT docs. Nothing was computed (no p(n) evaluated). Numbers are labelled: *quoted* (from a source), *counted* (I counted them in a source), *derived* (my own arithmetic from a formula; rough).

## 1. Is there a known bug in PARI's numbpart giving wrong values for large n (any version, in particular 2.15.4)?

### Takeaway
No. The whole PARI change history and the whole bug tracker hold only one wrong-value report for numbpart: bug #557 from 2007 (off by one for n ≥ 52602 in 2.3.1/2.4.1), fixed in 2.4.2 and in the 2.3 branch. There are also two 2003 faults (a crash near 10^15 and numbpart(0) = 0). Nothing has been reported for versions 2.5 to 2.19. The 2.15.5, 2.16.x, 2.17.x and 2.18.x changelogs have no numbpart item at all. 2.19.0 (28 Sep 2026) has only "speed up numbpart", under *Changed*, not *Fixed*. If the PARI value is the wrong one, it is an unreported bug.

### Cited Findings
- The function was added in 2.2.5. CHANGES-2.2, "Done for version 2.2.5 (released 11/02/2003)": `RS  9- routine numbpart()` (RS = Ralf Stephan) — [PARI 2.19.0 tarball, file CHANGES-2.2](https://pari.math.u-bordeaux.fr/pub/pari/unix/pari-2.19.0.tar.gz)
- CHANGES-2.2, "Done for version 2.2.7 (released 18/12/2003)", Fixed: `11- numbpart(10^15+2) --> SEGV  [ check for x < 10^15 incorrect ]` and `15- numbpart(0) --> 0 [ should be 1 ]` — [PARI 2.19.0 tarball, CHANGES-2.2](https://pari.math.u-bordeaux.fr/pub/pari/unix/pari-2.19.0.tar.gz)
- pari-dev, 10 Jul 2003, reply by Karim Belabas to a report that `numbpart(10^15+2)` crashed a 2.2.7 session: "Fixed in CVS. [ in floating point arithmetic, x < 10^15 is not the same as log_10(x) < 15 ... ]" — [pari-dev-0307/msg00022](https://pari.math.u-bordeaux.fr/archives/pari-dev-0307/msg00022.html)
- CHANGES-2.4, "Done for version 2.4.2 (released 22/12/2007)", Fixed: `12- numbpart inaccurate (and slow). E.g numbpart(52602) off by 1 [#557].` — [PARI 2.19.0 tarball, CHANGES-2.4](https://pari.math.u-bordeaux.fr/pub/pari/unix/pari-2.19.0.tar.gz); same text online at [fossies CHANGES-2.4](https://fossies.org/linux/pari/CHANGES-2.4)
- Bug #557, "numbpart(n) returns incorrect results for n>=52602", reported 28 Mar 2007 by Max Alekseyev: "The problem is reproduced in both the stable version 2.3.1 and the development version 2.4.1", on "amd64 running linux (x86-64/GMP-4.2.1 kernel) 64-bit version". Bill Allombert could reproduce it only in 64-bit mode (the 32-bit build was right) — [PARI bug #557](https://pari.math.u-bordeaux.fr/cgi-bin/bugreport.cgi?bug=557)
- PARI bug tracker list: 2,302 reports (*counted*: active and archived, numbers up to #2729). Only one title names numbpart: #557. Two others concern partitions: #671 (wishlist: "explicit partitions, number of restricted partitions") and #787 ("partitions()") — [PARI bug list, archive=both](https://pari.math.u-bordeaux.fr/cgi-bin/pkgreport.cgi?archive=both;package=pari)
- No numbpart item in the changelogs after 2.4.2 until 2.19.0:
  - 2.15.5 changelog (covers 2.15.0 to 2.15.5): no match for numbpart, partition or Rademacher — [pari-2.15.5.changelog](https://pari.math.u-bordeaux.fr/pub/pari/OLD/2.15/pari-2.15.5.changelog)
  - 2.16.2.beta (covers 2.16.0 to 2.16.2): none — [pari-2.16.2.beta.changelog](https://pari.math.u-bordeaux.fr/pub/pari/testing/OLD/2.16/pari-2.16.2.beta.changelog)
  - 2.17.4 (covers 2.17.0 to 2.17.4): none — [pari-2.17.4.changelog](https://pari.math.u-bordeaux.fr/pub/pari/unix/pari-2.17.4.changelog)
  - 2.18.1.alpha (covers 2.18.0 and 2.18.1): none — [pari-2.18.1.alpha.changelog](https://pari.math.u-bordeaux.fr/pub/pari/testing/pari-2.18.1.alpha.changelog)
  - 2.19.0 (released 28/09/2026), section *Changed*: `     3- speed up numbpart` (no author initials). Its *Fixed* section has no partition item — [pari-2.19.0.changelog](https://pari.math.u-bordeaux.fr/pub/pari/unix/pari-2.19.0.changelog)
  - A grep of every CHANGES* file in the 2.19.0 tarball for "numbpart|partition|rademacher" finds only the items above, plus `BA 138- partitions(1,[0,5],[3,4]) -> incorrect result` (2.10.0, a different function) — [PARI 2.19.0 tarball](https://pari.math.u-bordeaux.fr/pub/pari/unix/pari-2.19.0.tar.gz)
- PARI's own regression test (`src/test/in/partition`, 2.15.4) checks `numbpart(n)` against the pentagonal-number recurrence for every n below 10^4 (`test(10^4)`). It then evaluates `numbpart(0)`, `numbpart(52602)`, `numbpart(147007)` and `numbpart(10^15+2)`. No larger n is tested — [PARI 2.15.4 tarball](https://pari.math.u-bordeaux.fr/pub/pari/OLD/2.15/pari-2.15.4.tar.gz)
- The latest pari-users thread on numbpart (Apr 2026, "Re: numbpart(n, {a = k})") reports no bug. Bill Allombert: "numbpart use Rademacher formula, which is much faster than #partitions(n)" and "pari can compute numbpart(1000000) in 5 ms." — [pari-users-2604/msg00009](https://pari.math.u-bordeaux.fr/archives/pari-users-2604/msg00009.html)
- The author of the reference implementation used PARI's numbpart as a trusted alternative. In Sage ticket #29600, Fredrik Johansson (28 Apr 2020): "I propose removing the 'bober' implementation and supporting algorithm='pari' instead. Pari's numbpart() was improved in the recent past and it's now twice as fast as Bober's code on my computer." — [sagemath/sage#29600](https://github.com/sagemath/sage/issues/29600)

### Inferences
- If the PARI value p(56628877314) is the wrong one, nobody has reported this bug. It is not in the tracker, the CHANGES files or the mailing lists found.
- PARI's own tests check numbpart only below 10^4, plus three spot values up to 147,007. Correctness at n ≈ 10^10 to 10^11 rests on the algorithm's design, not on tests.

### Gaps
- I found no full-text search of the pari-users/pari-dev archives, and both web searches (standard and extended) found no thread about wrong numbpart values after 2007. A report in a thread that search engines did not index cannot be fully ruled out.
- The PARI gitweb returned HTTP 500 (for CHANGES and for source files), so the commit history of `src/modules/part.c` between releases was not read. The release tarballs were compared instead (see section 2).

## 2. Cause, affected n and fixing version — and does the 2.15.4 code explain our observation?

### Takeaway
The one known wrong-value bug (#557, 2007) had a numerical cause: the large sum was evaluated starting from the largest terms. It was fixed in May 2007 by reversing the order of summation and using the correct bit accuracy. Nothing is known for 2.15.4, so no version fixes our case. The 2.15.4 algorithm takes some shortcuts:
- a fixed number of terms, ⌊0.24√n + 5⌋, not a proven error bound;
- 32 guard bits;
- no check that the final sum lies close to an integer.

2.17.4 has the same code (one line differs). 2.19.0 rewrote it for speed but kept the same term count and guard bits.

Rough estimates (*derived*) put both the truncation error and the rounding error at n ≈ 5.7×10^10 far below 0.5. Also, the observed residues imply that the two values differ by at least 357, not by ±1. All documented failures of this formula (PARI #557, sympy, Maple) were off by one. So a marginal precision or truncation failure in numbpart is an unlikely explanation.

### Cited Findings
- Cause and fix of #557, from Karim Belabas's closing message (May 2007): "It is now gone in both versions. The problem was due to the unstable evaluation of a large sum [ starting from the *largest* terms ]. Reversing the order of summation fixes the problem. I also optimized the computation to use the correct bit accuracy throughout [ was often too large, even though the end result was wrong for large values ]." Earlier, the reporter had written (1 May 2007): "The bug seems to have gone in the development version 2.4.2 (CHANGES-1.1737). But it still exists in the stable version 2.3.2." — [PARI bug #557](https://pari.math.u-bordeaux.fr/cgi-bin/bugreport.cgi?bug=557)
- Header of `src/modules/part.c` in 2.15.4: "Original code contributed by: Ralf Stephan … Updated by Bill Allombert (2014) to use Selberg formula for L following http://dx.doi.org/10.1112/S1461157012001088". It implements the script `part(n) = round(sum(q=1,5 + 0.24*sqrt(n),L(n,q)*Psi(n,q)))` "only faster" — [PARI 2.15.4 tarball](https://pari.math.u-bordeaux.fr/pub/pari/OLD/2.15/pari-2.15.4.tar.gz)
- `numbpart()` in 2.15.4, `src/modules/part.c` lines 107–140 (quoted) — [PARI 2.15.4 tarball](https://pari.math.u-bordeaux.fr/pub/pari/OLD/2.15/pari-2.15.4.tar.gz):
  - input limit `if (cmpii(n, uu32toi(0x38d7e, 0xa4c68000)) >= 0) pari_err_OVERFLOW("numbpart [n < 10^15]");`
  - working precision `bitprec = (long)(rtodbl(est)/M_LN2) + 32;` where `estim(n)` is the low-precision log of exp(π√(2n/3))/(4n√3), i.e. about log2 p(n) + 32 bits;
  - number of terms `for (q = (ulong)(sqrt(gtodouble(n))*0.24 + 5); q >= 3; q--, set_avma(av))`, summing the smallest terms first, then q = 2 and q = 1;
  - a term is skipped when `if (abscmprr(t, mpexp(divru(est,q))) < 0) continue;` (with `est` negated, the threshold is p(n)^(−1/q); this drops A_q(n) values that are numerically zero);
  - each term `t = mulrr(t, psi(gprec_w(C, nbits2prec(bitprec / q + 32)), q, prec));`, and `L(n,k,…)` (Selberg's sum for A_k(n)) works at `nbits2prec(bitprec / k + k)`;
  - the result `return gerepileuptoint (ltop, roundr(mulrr(D,sum)));` is plain rounding, with no check of how far the sum lies from an integer.
- In `part.c` (2.15.4), doubles appear only for the term count (`gtodouble`, `sqrt`) and the bit estimate (`rtodbl`). The file has no threads and no `long double`. All other work is PARI t_REAL/t_INT arithmetic — [PARI 2.15.4 tarball](https://pari.math.u-bordeaux.fr/pub/pari/OLD/2.15/pari-2.15.4.tar.gz)
- The 2.15.4 manual entry: "$n$ must be of type integer and $n<10^{15}$ … The algorithm uses the Hardy-Ramanujan-Rademacher formula." It makes no accuracy statement — [PARI 2.15.4 tarball, src/functions/combinatorics/numbpart](https://pari.math.u-bordeaux.fr/pub/pari/OLD/2.15/pari-2.15.4.tar.gz)
- 2.15.4 and 2.17.4 differ in `part.c` by one line only: `sum = cgetr (prec); affsr(0, sum);` became `sum = utor(0, prec);` (*counted*, from diffing the two release tarballs) — [2.15.4](https://pari.math.u-bordeaux.fr/pub/pari/OLD/2.15/pari-2.15.4.tar.gz), [2.17.4](https://pari.math.u-bordeaux.fr/pub/pari/unix/pari-2.17.4.tar.gz)
- 2.19.0 rewrote `part.c` for speed — [PARI 2.19.0 tarball](https://pari.math.u-bordeaux.fr/pub/pari/unix/pari-2.19.0.tar.gz):
  - The new `SelbergL` builds the cosines from roots of unity (`rootsof1u_cx(12*q, prec)`), halves the range for odd q, and returns `NULL` when no l solves the congruence. This replaces the numerical skip.
  - The term count is unchanged: `for (q = (ulong)(sqrt(gtodouble(n))*0.24 + 5); q >= 3; q--)`. So are the guard bits: `prec = (long)log2pn(n) + 32; /* ~ log2 p(n) + 32 */`.
  - New limit and comment: "Impose n < 2^54 since the computation is in O~(n) anyway. This also ensures that 6(0.24 sqrt(n) + 5) fits in an ulong".
  - The manual now says: "$n<2^{54}$ … runs in time $\tilde{O}(n)$, which is roughly quadratic in the output size. Note that a more complicated implementation of the formula would run in time $\tilde{O}(n^{1/2})$." The same text is in the [online manual](https://pari.math.u-bordeaux.fr/dochtml/html/Combinatorics.html).
- Rademacher's proven remainder bound (Johansson eq. (1.8)), as coded in sympy 1.14, is `c1/sqrt(N) + c2*sqrt(N/(n - 1))*math.sinh(c3*sqrt(n)/N)` with `c1 = 44*math.pi**2/(225*math.sqrt(3))`, `c2 = math.pi*math.sqrt(2)/75`, `c3 = math.pi*math.sqrt(2/3)` — [sympy 1.14 partitions_.py](https://github.com/sympy/sympy/blob/sympy-1.14.0/sympy/ntheory/partitions_.py)
- Johansson (§3.7): Rademacher's bound "overshoots the actual magnitude of the terms with large k". In §1 he notes that the bound implies O(n^{1/2}) terms suffice. Per the fetch tool's reading of Table 2, his term counts correspond to about 0.26√n at n = 10^12, falling to about 0.18√n at 10^19 (the tool computed that ratio; I did not re-check it) — [Johansson 2012, ar5iv](https://ar5iv.labs.arxiv.org/html/1205.5991)
- The 0.24√n rule is not a proven bound. sympy ≤ 1.12 used the same rule with one term fewer (`M = max(6, int(0.24*n**0.5 + 4))`, then `for q in range(1, M)`). It returned values off by one from n = 1001. A commenter in that issue noted: "pari/gp seems to be using the same series but with a couple of more terms" — [sympy#24267](https://github.com/sympy/sympy/issues/24267), [sympy 1.12 source](https://github.com/sympy/sympy/blob/sympy-1.12/sympy/ntheory/partitions_.py)
- Every documented failure of this formula was off by a small amount: PARI #557 "off by 1"; sympy #24267 "off by one, sometimes up, sometimes down"; sympy #26750 "one off". Johansson also wrote: "past versions of Maple gave the wrong value of p(n)" for some n, presumably from floating-point roundoff and/or too few series terms — [bug #557](https://pari.math.u-bordeaux.fr/cgi-bin/bugreport.cgi?bug=557), [sympy#24267](https://github.com/sympy/sympy/issues/24267), [sympy#26750](https://github.com/sympy/sympy/issues/26750), [Johansson blog 2014](https://fredrikj.net/blog/2014/03/new-partition-function-record/)

### Inferences
- **Number of terms for our n** (*derived*). For n = 56,628,877,314, √n ≈ 237,968, and PARI 2.15.4 sums q = 1 … 57,117.
  - At N = 57,117, Rademacher's bound is about 1.3, so it does not certify PARI's truncation at this n. The bound drops below 0.5 only at about 63,100 terms (≈ 0.265√n). That is the count sympy 1.14 uses.
  - For n = 132,114,388,022 the bound at PARI's N (87,239, *derived*) is also above 0.5 (≈ 1.06, *derived*). Yet PARI and sympy agreed there. This fits the bound being very loose.
- **Truncation, estimated directly** (*derived*, rough). Using only |A_k(n)| ≤ k (the bound Johansson's §3.1 uses), an integral estimate of the terms PARI leaves out at n = 56,628,877,314 gives about 0.15 to 0.17. That is already below 0.5. With the usual cancellation inside A_k(n), the true tail should be orders of magnitude smaller.
- **Rounding** (*derived*, rough). The largest loss is in exp(C) for the q = 1 term, with C = π√(2n/3) ≈ 6.1×10^5 ≈ 2^19.2. With 32 guard bits this leaves an absolute error of order 2^(32−19.2) smaller than 1, i.e. about 2^−12.8 ≈ 10^−4, far below 0.5. The terms for q ≥ 3 have smaller absolute errors still.
- **So, as written, numbpart 2.15.4 should round to the right integer at this n** (by these rough estimates). Where the shortcuts do fail (small n, too few terms), the result is off by ±1.
- **The size of the disagreement.** The given residues differ: 63332 − 62975 = 357 (mod 43³ = 79,507), and 36 − 23 = 13 = 357 mod 43 (*derived*). So the two values differ by at least 357 in absolute value; they cannot differ by ±1. This is not the off-by-one pattern of every documented HRR failure. It points to a gross fault on one side, such as:
  - a wrong or missing big term;
  - a wrong limb or memory corruption in the 880,000-bit numbers;
  - a fault in the arithmetic library of that particular build;
  - or a slip in how one of the two values was obtained or reduced.
- **The same wrong residue came out in two separate sessions.** That argues for a deterministic cause (software, or data handling), not a random hardware fault.
- **No version fixes it.** No PARI version records a fix that could apply. 2.19.0 changes the inner routines (cosines from roots of unity, no numerical skip), so running another build (2.17.x or 2.19.0, a different GMP) is an independent check of the build and of the inner routines. It is not a check of the 0.24√n / 32-bit design, which is unchanged.
- **Sizes involved** (*derived*): p(56628877314) has about 880,600 bits, about 265,000 decimal digits. PARI works at about log2 p(n) + 32 bits for q = 1 and 2, and at about log2 p(n)/q + 32 bits for term q.

### Gaps
- I found nothing (bug, CHANGES line or message) that names a fault able to give an error of this size in PARI 2.15.4. The sources cannot tell which of the two values is wrong.
- Lehmer-type bounds on |A_k(n)| (sharper than |A_k(n)| ≤ k) were not sourced here, so they are not used.
- I found no documented faults in PARI's `mpexp`/`mpcos`/`mulrr` or in GMP at about 880,000-bit precision. None were found, but I did not search GMP's history in depth.

## 3. Are there reports specific to the cypari2 wheels or to WSL builds?

### Takeaway
None for PARI's numbpart or for cypari2. cypari2's tracker has no issue that mentions numbpart, and its integer-conversion issues are about compatibility, not wrong values. The one WSL precedent is Sage ticket #29600 (2020): wrong partition counts under WSL from Sage's old "bober" code, which relied on x87 `long double`. It was linked to WSL issue #830: WSL1 starts processes with the x87 control word set to 53-bit precision, where Linux expects 64-bit. PARI 2.15.4's source contains no `long double` at all, so this mechanism cannot reach PARI's own code.

### Cited Findings
- cypari2 2.2.0 wheels (manylinux2014 x86_64, i686, macOS arm64, …) were uploaded to PyPI on 2024-08-13 (*counted*: 13 files for 2.2.0). The next releases with wheels are 2.2.3 and 2.2.4 (Oct 2025) — [PyPI JSON for cypari2](https://pypi.org/pypi/cypari2/json)
- cypari2 2.2.0 wheel workflow: the sdist job sets `PARI_VERSION: pari-2.15.4`. For the wheels, "In the CIBW_BEFORE_ALL phase, we install libraries using the Sage distribution", with `TARGETS_PRE: pari` (Sage `develop` checked out, `make … pari`, `cibuildwheel==2.17.0`) — [cypari2 2.2.0 dist.yml](https://github.com/sagemath/cypari2/blob/2.2.0/.github/workflows/dist.yml)
- GitHub search `repo:sagemath/cypari2 numbpart` gives 0 issues or PRs. Searches for conversion problems found only compatibility issues: #114 "Incompatible with python 3.11", #137 "Incompatibility with Python 3.12", #169 "PyInt_AS_LONG is no longer defined". None reports a wrong value — [cypari2 issues](https://github.com/sagemath/cypari2/issues)
- Sage #29600, "number_of_partitions returns wrong results under the Windows Subsystem for Linux" (opened 2020-04-27, closed 2020-08-23, milestone sage-9.2) — [sagemath/sage#29600](https://github.com/sagemath/sage/issues/29600):
  - The report: "installing Sage under WSL 2 (Ubuntu 18.04) produces 8 doctests failures". The 'bober' algorithm differed from FLINT at n such as 241…271, 9894, 9908, 16564, 16573, …, 16969.
  - mkoeppe: "This might be related to: https://github.com/microsoft/WSL/issues/830 (long double floating point calculations give inexact results)" and "The bug is clearly in the ad-hoc float precision code in that file! It has ifdefs for the behavior of `long double`."
  - Another participant: "The WSL ticket you point to seems to claim that floating point issues were fixed in WSL2, which is what you specifically mention as failing the tests."
  - Resolution: the commit "#29600 : remove partitions_c.cc C code" and the redirection of `number_of_partitions` to FLINT.
- WSL issue #830, "long double floating point calculations give inexact results" (opened 2016-08-09, closed 2020-07-24) — [microsoft/WSL#830](https://github.com/microsoft/WSL/issues/830):
  - A Microsoft developer: "it is a difference between the default floating point control register … between Linux and NT. glibc assumes 0x37f but NT uses 0x27f."
  - Later comments note that WSL2 runs a real Linux kernel, while WSL1 users still complained in 2024 and 2025.
- Johansson (§3.4) on hardware floating point: "accuracy may be lost, for example, if the compiler generates long-double instructions which trigger double rounding" — [Johansson 2012, ar5iv](https://ar5iv.labs.arxiv.org/html/1205.5991)
- grep of the whole PARI 2.15.4 `src/` tree for "long double": 0 files (*counted*) — [PARI 2.15.4 tarball](https://pari.math.u-bordeaux.fr/pub/pari/OLD/2.15/pari-2.15.4.tar.gz)

### Inferences
- **Under WSL2**, which runs a real Linux kernel, the x87 issue should not apply. **Under WSL1** it would still not touch PARI's numbpart, which uses only integer-mantissa t_REAL arithmetic and plain doubles, never `long double`. It could touch other code that relies on `long double`.
- **A coincidence worth noting:** n = 16564 is one of the agreement checks between PARI and sympy, and it is also in Sage's WSL failure list for the old long-double code. PARI and sympy agreed there, so this gives no sign that the environment shows the x87 problem.
- **The cypari2 role.** numbpart runs entirely inside libpari; cypari2 only passes n in and the result out. If the residues were taken in Python after converting the ~880,000-bit result to a Python int, that conversion is one more step where a fault could enter. No cypari2 report describes such a fault.

### Gaps
- I did not determine the exact GMP version and configure flags that the 2.2.0 wheels bundle (PARI is built by Sage's build system at wheel-build time).
- The sources do not say whether the environment used was WSL1 or WSL2.

## 4. Is sympy's partition() (Hardy–Ramanujan–Rademacher with mpmath, versions 1.12–1.14) reliable for n around 10^10–10^11?

### Takeaway
sympy 1.12.x and earlier are not reliable: they used the same 0.24√n rule as PARI, with one term fewer, and gave off-by-one results from n = 1001. This was fixed for 1.13 (PR #24293, merged 2023-10-28). Versions 1.13 and 1.14 take the number of terms from Rademacher's proven bound and work with about 10% + 100 extra bits. That is a generous margin, but it is not interval arithmetic, so it is not proven. The docstring line "tested through 10^10" was shown false for the old code and has no published test record for the new one. Hard limit: n ≤ about 1.498×10^11 (both of our n are below it). I found no report of wrong values from 1.13 or 1.14.

### Cited Findings
- sympy 1.12 code: `M = max(6, int(0.24*n**0.5 + 4))` and `for q in range(1, M):`, docstring "The correctness of this implementation has been tested through $10^{10}$." — [sympy 1.12 partitions_.py](https://github.com/sympy/sympy/blob/sympy-1.12/sympy/ntheory/partitions_.py)
- Issue #24267 (opened 2022-11-16) — [sympy#24267](https://github.com/sympy/sympy/issues/24267):
  - "`npartitions(n)` gives incorrect results for n=1001, 1055, 1077, 1110, 1111, ... . The results are off by one, sometimes up, sometimes down." It was checked "against Sage, Pari, OEIS, and WolframAlpha".
  - A commenter, on raising the constant: "replacing 0.24 with 0.26 … seems to work correctly up to at least n=100000. But it would of course be preferable to have a method that provably works for all n."
- PR #24293, "Fix wrong results from npartitions" — [sympy PR #24293](https://github.com/sympy/sympy/pull/24293):
  - "Takes more terms in the series approximation and increases the arithmetic precision." Release note: "fixes the numerical evaluation of `npartitions`", filed for the 1.13 release notes.
  - Merged 2023-10-28 (GitHub API: `merged: True`, `merged_at: 2023-10-28T19:01:48Z`).
  - A comment (2022-11-23): "Maybe remove the docstring line about testing to $10^{10}$? (Apparently, not even the current version was tested properly.)"
- Issue #26750 (2024-06-28): in sympy 1.12.1, `npartitions(11160)` was "one off" compared with PARI's `numbpart(11160)`. Replies: "Looks like it's fixed in 1.13" (shown with 1.13.0rc3) and "This was fixed in gh-24293" — [sympy#26750](https://github.com/sympy/sympy/issues/26750)
- sympy 1.14 `_partition(n)` (quoted from source) — [sympy 1.14 partitions_.py](https://github.com/sympy/sympy/blob/sympy-1.14.0/sympy/ntheory/partitions_.py):
  - for n ≤ 200,000 near the cache it uses the exact pentagonal recurrence (`_partition_rec`); otherwise it uses HRR;
  - `prec = p = int(pbits*1.1 + 100)`;
  - the number of terms is found by bisection on Rademacher's bound `_M(n, N) < 0.5` ("using Rademacher's bound M(n, N) for the remainder after a partial sum of N terms (https://arxiv.org/pdf/1205.5991.pdf, (1.8))");
  - `if M > 10**5: raise ValueError("Input too big")  # i.e. n > 149832547102`;
  - the summation loop is `for q in range(1, M):`, and after each term `p = bitcount(abs(to_int(d))) + 50`;
  - A_k(n) comes from Whiteman's multiplicative formulas (`_a`, reference https://msp.org/pjm/1956/6-1/pjm-v6-n1-p18-p.pdf), using modular square roots and a factor table of size `maxn = 10**5`;
  - the result is `int(to_int(mpf_add(s, fhalf, prec)))`;
  - the deprecated `npartitions` wrapper still says "tested through $10^{10}$".
- `sympy.functions.combinatorial.numbers.partition.eval` returns `S(_partition(as_int(n)))` (numbers.py line 1648, 1.14) — [sympy 1.14 numbers.py](https://github.com/sympy/sympy/blob/sympy-1.14.0/sympy/functions/combinatorial/numbers.py)

### Inferences
- **One term short** (my reading of the code). The bisection leaves M as the smallest N whose bound is < 0.5, but `range(1, M)` sums only q = 1 … M−1. So the truncation is not strictly covered by the bound. Because the bound is very loose, this does not matter in practice.
- **Our n** (*derived*): for n = 56,628,877,314, sympy 1.14 sums about 63,000 terms (PARI: 57,117). Its spare precision is about 10% of ~880,600 bits, roughly 88,000 bits, against PARI's 32 bits. Its numerical margin at this n is much larger than PARI's.
- **Two independent algorithms.** sympy computes A_k(n) by Whiteman's formulas; PARI uses Selberg's cosine sum. A shared A_k bug cannot explain an agreement, and the two do not share code.
- **What a sympy fault would look like.** The difference is at least 357, so if sympy were the wrong side, its error would also have to be gross (e.g. a wrong A_k(n) for a small k, multiplying a huge term), not a precision effect.
- **Versions.** sympy 1.12 is unreliable at every size above about 10^3. Versions 1.13 and 1.14 are bound-based, but they are unproven and untested at 10^10 to 10^11.

### Gaps
- I found no published test record for sympy ≥ 1.13 at n ≥ 10^10, and no issue reports either way at these sizes.
- I did not find a stated accuracy guarantee for mpmath's `mpf_cos`, `mpf_sin` and `mpf_cosh_sinh` at these precisions.

## 5. The reference implementation (Johansson, FLINT/Arb) and its verification record

### Takeaway
Fredrik Johansson's implementation is the reference:
- FLINT, 2012 (used for p(10^19)). This version was "heuristically sound" but not proven.
- Arb, 2014 (used for p(10^20)). This version uses ball arithmetic, a proven truncation bound and a check that the final ball holds exactly one integer.
- In FLINT 3, `partitions_fmpz_fmpz` is documented as exact.

It is the accepted arbiter for a disputed p(n). Sage's `number_of_partitions` calls FLINT.

### Cited Findings
- Johansson, "Efficient implementation of the Hardy–Ramanujan–Rademacher formula", LMS J. Comput. Math. 15 (2012) 341–359, DOI 10.1112/S1461157012001088; arXiv 1205.5991 (v1 27 May 2012, v2 6 Jul 2012). The abstract reports a speed-up of more than a factor 500 over previously published software, and the computation of p(10^19) (paraphrased) — [arXiv:1205.5991](https://arxiv.org/abs/1205.5991)
- What the paper verified and what it left unproven — [Johansson 2012, ar5iv](https://ar5iv.labs.arxiv.org/html/1205.5991):
  - §3.7: "All n ≤ 10^6 have been checked explicitly"; a large number of isolated n ≫ 10^6 were "compared against known congruences and values computed with Sage and Mathematica".
  - §3.3: "our implementation currently relies on routines that, although heuristically sound, have not yet been proved correct".
  - §3.6: "We do not attempt to prove that the internal precision management of Algorithm 3 is correct."
  - §4.1: a multithreaded test "failed when computing 10^19 for reasons yet to be determined".
  - §3.1: N is chosen so that |R(n,N)| < 0.25, and n > 2000 is assumed for the error analysis.
- Johansson's page — [fredrikj.net/math/hrr.html](https://fredrikj.net/math/hrr.html):
  - "A newer implementation is available as part of the Arb library. This implementation is provably correct (by use of ball arithmetic) and slightly faster for huge $n$."
  - "Version 5.11 and later of Sage now uses this implementation by default."
  - Errata: per David Broadhurst, "the formula for the number of partitions into odd (or distinct) parts given at the end of the paper is incorrect"; two typos reported by a reader.
- Blog post, 2 Mar 2014, on p(10^20) (11,140,086,260 digits, *quoted*) — [fredrikj.net blog 2014](https://fredrikj.net/blog/2014/03/new-partition-function-record/):
  - the code "uses ball (mid-rad interval) arithmetic to bound the contribution of all rounding errors";
  - "The error that results from truncating the infinite series after a finite number of terms is also accounted for";
  - the final interval "contains exactly one integer, which then provably is the correct value", though "there could be bugs in the code (or even hardware errors)";
  - the distance to the nearest integer was "smaller than 10−10";
  - p(10^10) through p(10^16) were re-run on a 32-bit system and "all came out with the correct value";
  - a read-back mismatch was traced to an overflow in `mpz_out_raw`;
  - the run took 110 hours and 130 GB (Xeon E5-2650, 256 GB RAM), using MPIR 2.6.0 and MPFR 3.1.1.
- Arb documentation — [arblib.org partitions](https://arblib.org/partitions.html):
  - `partitions_fmpz_fmpz` "computes a numerical ball containing $p(n)$ and verifies that the ball contains a unique integer".
  - `partitions_hrr_sum_arb`: "If use_doubles is nonzero, doubles and the system's standard library math functions are used to evaluate the smallest terms … unless the system's standard library is broken, use of doubles can be considered safe." and "Setting use_doubles to zero gives a fully guaranteed bound."
  - Module introduction: it "uses ball arithmetic throughout to guarantee a correct error bound"; the older FLINT version "relies on some numerical subroutines that have not been proved correct."
- FLINT 3 documentation — [flintlib.org partitions](https://flintlib.org/doc/partitions/index.html):
  - "Set p to p(n), which is 0 for n < 0. The argument use_doubles is ignored."
  - "The result is exact: the formula is evaluated in ball arithmetic with a proved bound for the truncation error," and "the function throws if the final ball does not determine a unique integer (which is not expected to happen)".
  - `partitions_rademacher_bound`: "The constants are rounded upwards and the bound is evaluated with directed rounding at about 30 bits."
  - Timing example: "p(10^9), a number of 117 thousand bits, takes a few hundredths of a second."
  - Caveat: the page says "Parts of this text have been generated by an LLM."
- Sage redirected `number_of_partitions` to FLINT in 2020 (Sage 9.2), when the old 'bober' code was removed — [sagemath/sage#29600](https://github.com/sagemath/sage/issues/29600)

### Inferences
- **The documented arbiter.** For p(56628877314), FLINT ≥ 3 / Arb `partitions_fmpz_fmpz` (or Sage's `number_of_partitions`, which calls FLINT) is the documented arbiter. Unlike PARI and sympy, it checks itself: it proves that the result ball holds a unique integer.
- **Run time** (*derived*, rough). FLINT's complexity is near Õ(n^{1/2}), against PARI's Õ(n). From the FLINT doc's p(10^9) timing, n ≈ 5.7×10^10 should take far less time in FLINT than PARI's ~450 s.

### Gaps
- I found no independent published audit of Arb's partitions code beyond the author's own statements and tests.
- I did not check which Python bindings expose `partitions_fmpz_fmpz` (e.g. python-flint), or their exact function names.
