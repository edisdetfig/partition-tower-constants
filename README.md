# Partition tower constants

A computational study of the partition function p(n) modulo powers of a prime ℓ, and an independent recheck of it.

For each prime ℓ ≥ 5, the numbers p(ℓ^b·n + δ_b), with 24·δ_b ≡ 1 (mod ℓ^b), form a "tower" with one level for each b.
Going up two levels multiplies the deep levels by fixed ℓ-adic constants. For ℓ ≥ 13 this project computes those
constants from newforms of weight ℓ − 3 on Γ0(6), followed in ℓ-adic families to weight −2, and tests them against
exact partition numbers. The primes 5, 7 and 11 are treated separately: no such newform exists there, the tower is first
divided by the powers of ℓ that Watson's and Atkin's congruences give, and its one constant comes from a family of
positive slope (statement C).

## Status, read this first

- **AI-assisted.** The computations and notes were made in sessions with Claude (Anthropic), 11 to 20 Nov 2026,
  directed by the repository owner, who is not a mathematician.
- **Not peer reviewed.** As far as we know, no specialist has read any of it. The central statement (A1) is a
  conjecture with an unrefereed draft proof.
- **Independently rechecked in part.** A separate Claude session recomputed a sample of every kind of result with its
  own programs (`recheck/`) and wrote the surveys. It found no numerical disagreement, and several errors of wording
  and attribution, which are corrected and logged.

Every statement carries one of four labels: **published**, **derived** (a short argument from published results, not
refereed), **computed**, or **conjectured**. `STATEMENT.md` and `NOTE.md` use these four words. Two working files keep
their own words, each mapped to these four in its opening key: the notebook (counted, measured, derived) and
`FINDINGS.md` (reading, derived here). Nothing the project argued itself is called "proved".

## Where to start

| If you want | Read |
|---|---|
| One calculation, start to finish | [`recheck/scripts/week2/example13.py`](recheck/scripts/week2/example13.py) (ten lines; needs python-flint) |
| A plain-language guide (start here if you are new to the subject) | [Partition towers in plain words (PDF)](surveys/Partition%20towers%20in%20plain%20words.pdf) |
| The technical survey: what is known, new and unproved, with sources | [The partition tower modulo prime powers: a computational survey (PDF)](surveys/The%20partition%20tower%20modulo%20prime%20powers%20a%20computational%20survey.pdf) |
| The exact statements, A to H | [`STATEMENT.md`](results/side_finds/STATEMENT.md) |
| The draft proof of the main statement | [`NOTE.md`](results/side_finds/NOTE.md) |
| The short arguments for A3 and E1 | [`PROOF.md`](results/side_finds/PROOF.md) |
| Every finding with its status, and corrections | [`FINDINGS.md`](results/side_finds/FINDINGS.md) |
| What is open | [`NEXT.md`](results/side_finds/NEXT.md) |
| The day-by-day working notebook | [`results/side_finds/README.md`](results/side_finds/README.md) |
| What is in print | [literature check](results/side_finds/reports/Partition%20tower%20literature%20check.md) |

## Layout

- `results/side_finds/`: the project's record as it was made. Scripts (`.py`) sit beside their kept outputs
  (`.txt`, `.json`, `.jsonl`). Usage lines in the scripts assume this path.
- `recheck/`: the reviewer's own programs and logs. None of the project's scripts was used for these.
- `surveys/`: the two write-ups, one for specialists and one in plain words.

## Running things

Python 3.11 or later with the packages in `requirements.txt`. Exact partition numbers come from FLINT
(Hardy–Ramanujan–Rademacher with a proved error bound). A few scripts use PARI through `cypari2`, which has no Windows
wheel; those were run under WSL.

## Known limits

- A1, B and C are conjectures supported by computation.
- The settling-rate result (every prime from 29 to 613 settles one digit per two levels) assumes A1 for the quiet
  families.
- Some first digits were already in print (Ono 2000; Boylan–Webb 2013; Narayanan 2019). The literature report lists
  them.
- Full rank of the tower (that it uses every form of the class) is computed, not proved. See `NOTE.md`, section 5.
- The folder holds more than the surveys cover: digs on mock theta functions, moonshine series, poles of modular
  functions and others. Those were not part of the independent recheck, except where the survey says so (section 7).
- The notebook keeps a few references to material that is not in this repository (sums of three cubes, and other
  studies in the same workspace).

## Licence

Code: MIT (`LICENSE`). Text and data: CC BY 4.0 (`LICENSE-docs.md`).

Corrections are welcome: please open an issue.
