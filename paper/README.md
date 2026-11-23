# Standalone draft note (`paper/`)

This directory holds a self-contained LaTeX extraction of the core mathematical note (`results/side_finds/NOTE.md`) for specialist readers who want a compact, referee-ready document without navigating the full 700-file working archive in `results/side_finds/`.

- [`partition_tower_constants.tex`](partition_tower_constants.tex): standalone 6--8 page `amsart` draft containing:
  - **Theorem 1** (ordinarity of every $f \in S_\ell$ for $\ell \ge 7$ via theta cycles),
  - **Theorem 2** (the classical-weight proof of the $\Lambda$-adic recurrence, with the saturated-module Hida argument avoiding level-raising issues and explicit quantifier uniformity in $b_0, j_0$),
  - A compact table of computed constants ($5 \le \ell \le 31$), the Boylan--Webb $\ell = 29$ erratum, the $\ell = 607$ rank-deficit explanation, and sample explicit congruences modulo $37, 43, 67, 73$.

## Building

```bash
pdflatex partition_tower_constants.tex
pdflatex partition_tower_constants.tex
```
