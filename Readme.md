# The Non‑Reversibility Selector (NRS)
**Subtitle:** *Gibbs Asymmetry as a Principle of Epistemic Fitness*

This repository contains the canonical LaTeX source for **NRS**, the Mathematica support code for computation/experiments, and Lean/Coq skeletons for formalization. It integrates with the broader idiom‑crossing runtime series:
**QFI**, **Idiom Projections**, **Modal Shell**, **Ergodic Metric Drift**, **Parallel Murphy Budgeting**, **Measure‑Theoretic Flow Centrality**, and **Large A: Graded Effect Calculus & Safety Contracts**.

---

## Source hierarchy (authoritative order)

1. **`latex/main.latex` — canonical manuscript body.**  
   - If restructuring is needed, this file is replaced as a *monolithic unit*.
   - Uses `bib/main.bib` for citations and the `alpha` style.
2. **Mathematica support** (preferred for modeling/symbolics/simulations):
   - `mathematica/NRS_Core.wl` — EPR, DB checks, stationary distributions, QFI metric, centrality.
   - `mathematica/NRS_Murphy.wl` — greedy/DP schedulers, submodularity certificate.
3. **Proof skeletons** (optional, for hardening key results):
   - **Lean 4:** `lean/NRS/Kernel.lean`, `lean/NRS/Selector.lean`
   - **Coq:** `coq/NRS.v`

> **Policy:** Always model or compute in **Mathematica** first (derivations of asymmetries in Gibbs free energy, Murphy‑budget optimization, flow‑curvature, morphospace metrics). Lean/Coq serve as formal landing pads for core theorems once definitions stabilize.

---

## Build & use

### LaTeX (Overleaf or local)
- Root file: `latex/main.latex`  
- Bib file: `bib/main.bib`  
- Recommended engine: `pdfLaTeX`  
- Sequence (local):  
  ```bash
  pdflatex main.latex && bibtex main && pdflatex main.latex && pdflatex main.latex