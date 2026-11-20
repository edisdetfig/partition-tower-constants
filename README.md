# Partition tower constants

A computational study of the partition function p(n) modulo powers of a prime ℓ, and an independent recheck of it.

For each prime ℓ ≥ 5, the numbers p(ℓ^b·n + δ_b), with 24·δ_b ≡ 1 (mod ℓ^b), form a "tower" with one level for each b.
Going up two levels multiplies the deep levels by fixed ℓ-adic constants. This project computes those constants from
newforms of weight ℓ − 3 on Γ0(6), followed in ℓ-adic families to weight −2, and tests them against exact partition
numbers.

## Status, read this first

- **AI-assisted.** The computations and notes were made in sessions with Claude (Anthropic), 11 to 20 Nov 2026,
  directed by the repository owner, who is not a mathematician.
- **Not peer reviewed.** As far as we know, no specialist has read any of it. The central statement (A1) is a
  conjecture with an unrefereed draft proof.
- **Independently rechecked in part.** A separate Claude session recomputed a sample of every kind of result with its
  own programs (`recheck/`) and wrote the surveys. It found no numerical disagreement, and several errors of wording
  and attribution, which are corrected and logged.

Every statement carries one of four labels: **published**, **derived** (a short argument from published results, not
refereed), **computed**, or **conjectured**.

## Where to start

| If you want | Read |
|---|---|
| One calculation, start to finish | `recheck/scripts/week2/example13.py` (ten lines; needs python-flint) |
| The overview, with what is known, new and unproved | the surveys in `surveys/` |
| The exact statements, A to H | `results/side_finds/STATEMENT.md` |
| The draft proof of the main statement | `results/side_finds/NOTE.md` |
| The short arguments for A3 and E1 | `results/side_finds/PROOF.md` |
| Every finding with its status, and corrections | `results/side_finds/FINDINGS.md` |
| What is open | `results/side_finds/NEXT.md` |
| The day-by-day working notebook | `results/side_finds/README.md` |
| What is in print | `results/side_finds/reports/Partition tower literature check.md` |

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
- The notebook keeps a few references to material that is not in this repository (sums of three cubes, and other
  studies in the same workspace).
- The folder holds more than the surveys cover: digs on mock theta functions, moonshine series, poles of
  modular functions and others. Those were not part of the independent recheck, except where the survey
  says so (section 7).

## Licence

Code: MIT (`LICENSE`). Text and data: CC BY 4.0 (`LICENSE-docs.md`).

Corrections are welcome: please open an issue.
