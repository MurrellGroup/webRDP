# Changelog

## 0.5.1 — 2026-08-11

- Made every reusable analysis panel own an explicit bounded scroll body on
  desktop, with a persistent scrollbar, keyboard focus, and scroll chaining
  back to the center workspace at its edges.
- Gave the center workspace a definite zero-basis flex height so its overflow
  cannot be absorbed by an intrinsic grid row on short viewports.
- Removed the tutorial dialog's clipping override and made its navigation rail
  independently safe on short desktop and narrow mobile screens.
- Expanded rendered CSS regressions to cover all primary scroll owners and the
  responsive handoff from nested desktop scrolling to document scrolling.

## 0.5.0 — 2026-08-11

- Rebuilt viewport sizing and overflow ownership so the center workspace,
  settings, and inspector scroll reliably on short, desktop, and mobile screens.
- Promoted local phylogenies to a dedicated event-centered workspace with
  flank/tract trees, adaptive nearest-sequence cohorts, role highlighting,
  topology summaries, and per-tree/batch Newick export.
- Added explicit method result tabs with named supporting methods, raw and
  adjusted p-values, statistic/calibration details, correction scope, and an
  honest per-method fidelity warning.
- Added eight deterministic synthetic example datasets spanning clean triplets,
  virus families, circular and segmented genomes, nested mosaics, an 80 kb
  bacterial core panel, and a 520-genome performance stress test.
- Made sequence-changing exports accepted-and-fresh by default, preserved full
  scan multiplicity during event recalculation, invalidated breakpoint bounds
  after manual edits, retained alternative parents, and introduced an immutable
  project-level audit ledger in `.rdpweb` schema 0.5.
- Avoided disabled BootScan bootstrap and 3SEQ exact-DP work, replaced the fixed
  500-candidate ceiling with an adaptive per-recombinant bound, and surfaced
  HMM informative-site breakpoint intervals.

## 0.4.0 — 2026-08-10

- Added exact bounded 3SEQ HGRW first-passage probabilities with exhaustive
  small-path and published 30/30 regression cases.
- Added deterministic native p-distance bootstrap support and reproducible
  seed/replicate controls.
- Added four-gamete, PHI-style proximity, rate-density and parent-conflict
  challenge diagnostics plus null/rate-shift fixtures.
- Added manual event creation, duplication, deletion, reassignment, grouping,
  draggable breakpoints, exact-hypothesis recalculation, undo/redo and audit
  history in `.rdpweb` schema 0.4.
- Added local neighbor-joining trees with Newick export, circular breakpoint
  matrices and seeded hotspot permutation tests.
- Added GFF3, GenBank FEATURES and BED annotation import, a feature track,
  GFF3 export and CDS-phase-aware masking.
- Preserved static GitHub Pages operation at root and repository subpaths.
