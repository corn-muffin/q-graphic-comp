# Deletion-Contraction for q-Arrangements on Simplicial Complex Facets

Computational exploration of a generalization of Nian's edge deletion-contraction
theorem (for the q-deformed arrangements S^q_Δ) to higher-dimensional facets.
Built with an LLM (Claude Code) — see caveat below.

## Contents

- `scripts/qarr.sage` — core library: builds S^q_Δ and computes χ(S^q_Δ,t) via a fast
  matroid lattice-of-flats method.
- `scripts/isolating_utils.sage` — privacy/isolation/connectivity checks.
- `scripts/*.sage` — individual experiments (counterexample search, bridge families,
  random sweeps at n=5,6, modularity check, etc.); filenames are descriptive.
- `scripts/*.log` — captured output from the longer runs.
- `data/*.json` — raw computed polynomials.

## Main result (Conjecture, unproven)

For a maximal face F of Δ with an **isolating** edge e={i,j} ⊂ F — private to F
(no other face contains both i,j) *and* such that i,j land in different components
of the 1-skeleton once F's own edges are removed —

    χ(S^q_Δ,t) = χ(S^q_(Δ\F),t) − (q−1)^(|F|−1) χ(S^q_(Δ/e),t).

Reduces exactly to Nian's Theorem 3.4 when |F|=2. Verified computationally on
thousands of connected complexes (n=5 exhaustively, n=6 by random sampling),
zero counterexamples to sufficiency found; known **not** necessary, and known to
fail for complexes with no isolating edge (e.g. a "prism" example: two triangles
plus three connecting edges, sides not filled in).

## Requirements

SageMath (10.x). Run any script with `sage scriptname.sage` from inside `scripts/`.

## Caveat

This was produced with heavy LLM assistance (code and computation). Verify
independently before relying on it.
