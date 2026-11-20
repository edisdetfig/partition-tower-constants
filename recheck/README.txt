Independent recheck of the partition-tower record (review of 13-14 Nov 2026)
===========================================================================

These are the reviewer's own programs and their outputs. None of the project's
scripts was run. Requirements: Python 3, python-flint 0.9.0, numpy, sympy.
Exact partition numbers come from FLINT (fmpz.partitions_p), which evaluates the
Hardy-Ramanujan-Rademacher formula with a proved error bound.

scripts/
  consts.py      constants of 5..31 from exact p(N): pairs of levels two apart, n = 0   -> logs/consts.log
  twins.py       the laws at 37..73 with the project's P_l, lanes l, l^2, l^3, plus controls -> logs/twins.log
  three.py       the four-term law on lane l at 67, 71, 79, 83                           -> logs/three.log
  third.py       third-digit control at 37, levels (5,7,9)
  lucky.py       zeros and controls at lucky primes (37, 43; 2999; the mod-67 example)   -> logs/lucky.log
  ctrl67.py      controls for the mod-67 example (caution: bound the loop before running)
  big73.py       p(9040523964878168) mod 73 (about 29 minutes, 1 GB)                     -> logs/big73.log
  vr2.py         class S_l and T_l modulo l from eta quotients on Gamma0(6): python3 vr2.py <l>
  census.sh      vr2.py for every prime 13..199                                          -> logs/census_now.log
  cmpP.py        compares the project's P_l mod l with vr2's result at 13 primes
  tq.py          T_Q on the class modulo l: python3 tq.py <l> <Q,Q,...> [full]            -> logs/tq37/43/67/73.log
  loops.py       orders of the constants modulo l and l^2 (5..61)
  colours.py     constants of 3- and 5-coloured partitions from exact values             -> logs/colours.log
  colours511.py  5 colours at 11 (the last case of colours.py, run alone)
  level1.py      1/E4 at 7 and 13 (pole constant u^-6) and j at 13                       -> logs/level1.log
  earlier/       checks made earlier in the review (page values, 24 colours, the split at 37)

Notes
  * In lucky.py the sign "eps" is for M = (24N-1)/Q^2. The record lists signs for n = M/l.
    They differ by (l|Q), which is -1 at Q = 6599 (l = 37) and Q = 4049 (l = 43).
  * logs/tq37.log and logs/tq43.log were written before a fix to tq.py for primes Q smaller
    than the dimension; their rows for Q = 5, 7, 11, 13 are superseded (rerun tq.py to get
    the corrected polynomials). Their rows for Q > 1000 are unaffected.
  * colours.log ends with "Killed": the last case ran out of memory alongside other jobs
    and was rerun alone as colours511.py.
  * The polynomials P_l inside twins.py, three.py and loops.py are the project's numbers,
    copied from its output files; they are inputs to the check, not results of it.

scripts/week2/  (checks of 14 and 20 Nov 2026; see section 8 of the survey)
  example13.py   the worked example at 13: level 3 = 6 x level 1 (mod 13), level 5 = 45 x level 3 (mod 169)
  dims.py        PARI count of the class at weights l-3 and 2l-4: the first-step gain floor(l/24), 13..73
  partA.py       exact p(N) for the towers of 7 and 11; anaA.py reads the quiet ratio E_(b+2)/E_b from them
  twoterm.py     the same quiet constant from a two-term fit, with no separately estimated alpha
  partB.py, sysB.gp, lagB.py   the 7-adic and 11-adic Hecke systems on Gamma0(6) (PARI) carried to weight -2
  share.py       the share test from its definition, at 13, 23, 29, 31, 37, 53
  braid8.py      607: the first-digit polynomial of the group of 8 is irreducible; its roots cycle
  Requirements beyond the first week: cypari2 (PARI 2.17).

Paths: the scripts were written in the reviewer's sandbox. Absolute sandbox paths were replaced with relative ones on
20 Nov 2026 (cmpP.py, census.sh, twins.py, consts.py, and one traceback line in logs/lucky.log); nothing else in them
was changed. consts.py and twins.py now write their .json beside wherever they are run. The files in logs/ are the
original outputs and are part of the record.
