# Standalone draft note (`paper/`)

This directory holds a self-contained LaTeX extraction of the core mathematical draft (`results/side_finds/NOTE.md`) for specialist readers who want a compact document without navigating the full 700-file working archive in `results/side_finds/`.

- [`partition_tower_constants.tex`](partition_tower_constants.tex): standalone 7-page `amsart` draft (unrefereed, AI-assisted) containing:
  - **Theorem 1** (*derived*): ordinarity of every $f \in S_\ell$ for $\ell \ge 7$ via theta cycles,
  - **Theorem 2** (*conjectured*, with a proposed classical-weight proof): the $\Lambda$-adic recurrence modulo $\ell^m$, using a saturated-kernel Hida construction and explicit quantifier uniformity in $b_0, j_0$,
  - A compact table of computed constants ($5 \le \ell \le 31$), the apparent Boylan--Webb $\ell = 29$ display misprint, the $\ell = 607$ rank-deficit explanation, and sample explicit congruences modulo $37, 43, 67, 73$.

## Building

```bash
pdflatex partition_tower_constants.tex
pdflatex partition_tower_constants.tex
```
