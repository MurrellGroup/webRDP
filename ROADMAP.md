# RDP Web replacement roadmap

RDP Web will not call itself a validated RDP5 replacement merely because the
interface and method names exist. The replacement claim is gated by fixtures,
simulation results, workflow parity, and browser-scale performance.

## Track A — primary detection methods

- RDP: reproduce reference selection, informative-window enumeration, and
  local/global p-values on documented triplets.
- GENECONV: add mismatch/G-scale scoring and clean-room KA calibration; add an
  optional native permutation implementation without bundling GENECONV.
- BOOTSCAN/RECSCAN: seeded native triplet bootstrapping is operational; retain
  multi-taxon NJ/substitution-model detection and RDP5 cutoff parity work.
- MAXCHI and CHIMAERA: finish polymorphic-window semantics, peak pairing, and
  corrected-p parity.
- SISCAN: implement nearest/random outgroup selection, 15 site categories,
  horizontal/vertical randomization, reproducible seeds, and fast/full modes.
- 3SEQ: exact HGRW first-passage dynamic programming is operational within a
  browser work budget; finish large-case lookup/algorithm comparisons and
  RDP5 corrected-p validation without copying restricted code or tables.

## Track B — false-positive control and breakpoint inference

- Extend the operational candidate-seeded two-state, windowless Viterbi kernel
  to BURT's documented 2–20-state step-up fitting, repeated random starts, EM
  training, and validated confidence intervals.
- Bounded four-gamete, seeded proximity-permutation, PHI-style contrast,
  rate-variation, parent-conflict, gap and misalignment diagnostics are exposed
  as review evidence; add exact PHI calibration and tree-conditioned homoplasy
  tests.
- The dual-origin circular-event path now retains alternative parents and
  overlapping/nested candidates in automatic ancestry groups; expand this to
  explicit reassortment, unknown-parent, and recombinant-parent models.
- Parametric and column-preserving null simulations with reproducible seeds.

## Track C — interactive analysis parity

- Manual events, draggable/numeric breakpoints, scan-scope exact-hypothesis
  recalculation, and seven bounded method-specific review profiles are
  operational; add independent method-specific candidate scan canvases.
- Grouping/ungrouping, recombinant/parent reassignment, unresolved-sequence
  rescans, undo/redo, per-event history, an immutable project ledger,
  duplication, deletion, IndexedDB autosave, and partial-scan checkpoints are
  operational. The ordered Global reconstruction queue, stale downstream
  propagation, mosaic map, possible-overprint links, and recombinant-parent
  dependencies are operational. A tunable heuristic now resolves fresh events
  in order, holds ambiguous dependent branches, and schedules targeted/adaptive
  rescans; calibrate it against authorized RDP5 outputs and add full signal
  erasure/splitting plus richer ancestral-event propagation semantics.
- Dedicated flank/tract/right-flank NJ comparison, nearest-context cohorts,
  role highlighting, linked leaf marking, connected SVG geometry, Newick
  export, circular breakpoint pair matrices and seeded hotspot permutations
  are operational; ML, RF/SH/AU, clade editing and association tests remain.
- Arbitrary 2–6-parent alignment affinity highlighting, shared/novel/missing
  states, informative-site filtering, coordinate navigation and bounded
  large-alignment behavior are operational; add amino-acid/codon overlays and
  scalable overview minimaps.
- GenBank/GFF3/BED annotations, feature-coordinate mapping, GFF3 export and
  CDS-phase-aware masking are operational; add ORF calling and complex
  circular/multi-record feature mappings.
- PDB/SCHEMA, LDHat interoperation, and ancestral-sequence workflows as
  clearly isolated companion modules where licensing or browser constraints
  make a monolith undesirable.

## Track D — validation and performance

- Positive, independent-mutation negative, gap-block, circular, rate-shift and
  large-path fixtures are automated; expand low/high-diversity and
  recombinant-of-recombinant corpora.
- Numerical comparison ledger against RDP5 outputs supplied by authorized
  users; no proprietary code or restricted binary is redistributed.
- Sensitivity/specificity and breakpoint-error curves across published
  simulation designs.
- Chrome, Firefox, Safari, and Edge benchmarks at 100 × 10 kb, 4,000 × 10 kb,
  and long-genome workloads with time and memory ceilings.
- The `.rdpweb` 0.5 migration path and same-engine Node batch runner are
  operational; add a long-term migration corpus and standalone reproducible
  analysis manifests.

## Release gates

- **0.x scientific alpha:** useful interactive hypotheses, explicit calibration
  labels, no parity claim.
- **1.0 candidate:** all primary methods implemented; every result field has a
  fixture or simulation test; secondary gaps are explicit.
- **Validated replacement:** RDP5 task/workflow parity, numerical comparison,
  false-positive challenge suite, cross-browser scale targets, and independent
  scientific review are all complete.
