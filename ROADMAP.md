# RDP Web replacement roadmap

RDP Web will not call itself a validated RDP5 replacement merely because the
interface and method names exist. The replacement claim is gated by fixtures,
simulation results, workflow parity, and browser-scale performance.

## Track A — primary detection methods

- RDP: reproduce reference selection, informative-window enumeration, and
  local/global p-values on documented triplets. Multi-excursion source
  enumeration is operational; complete the desktop tie/missing/circular corpus.
- GENECONV: source finite-G integer mismatch scoring and the
  `CalcKMaxP`/`GCCalcPValP` KA path are operational, including G=0; finish
  indel-run modes, overlap suppression, pair scans and native permutations.
- BOOTSCAN/RECSCAN: seeded native triplet bootstrapping is operational; retain
  multi-taxon NJ/substitution-model detection and RDP5 cutoff parity work.
- MAXCHI and CHIMAERA: compressed informative-site windows, independent peak
  pairs and the source half-window × three-pair correction are operational;
  finish `GrowMChiWin` peak growth, smoothing/ban-window edge behavior and
  desktop table parity.
- SISCAN: source nearest-tree/direct, most-divergent, manual and randomized
  fourth-sequence paths; 15 site categories and topology sums;
  horizontal/vertical randomization; reproducible seeds; topology-run
  enumeration and shrinkage; and bounded-memory fast/full controls are
  operational. Finish desktop tie/missing-data and plot-by-plot corpus parity.
- 3SEQ: exact HGRW first-passage dynamic programming is operational within a
  browser work budget; finish large-case lookup/algorithm comparisons and
  RDP5 corrected-p validation.

## Track B — false-positive control and breakpoint inference

- Source-compatible fixed-three-state BURT fitting, repeated random starts,
  historical backtrace behavior, and posterior confidence intervals are
  operational; validate exact circular/tie behavior against a desktop corpus
  and finish the optional manual 2–20-state mode's selection fixtures.
- Direct `PHITest2`/`PHI` multistate incompatibility scoring, analytic moments
  and normal-tail p-values now accompany the bounded four-gamete permutation,
  rate-variation, parent-conflict, gap and misalignment diagnostics. Add a
  desktop PHI golden corpus, exact VB thinning above 6,000 informative sites,
  and tree-conditioned homoplasy tests.
- The dual-origin circular-event path now retains alternative parents and
  overlapping/nested candidates in automatic ancestry groups; expand this to
  explicit reassortment, unknown-parent, and recombinant-parent models.
- Parametric and column-preserving null simulations with reproducible seeds.

## Track C — interactive analysis parity

- Manual events, draggable/numeric breakpoints, scan-scope exact-hypothesis
  recalculation, seven bounded method-specific review profiles, and a
  source/manual-guided six-stage Review studio are operational. The primary
  surface now includes characterization-order navigation, best-unresolved and
  method-count filtering, source-weighted role-polarity audition, group-level decisions, and
  tract/background tree verification. All seven families now own independent
  discovery intervals and co-location gates; add multi-signal raw queues and
  richer source plot controls for every non-RDP family.
- The desktop-default `MakeConsensusC` role path now covers its 18 standalone
  statistics, final-trim penalty, and six joint rules, including `FindSets`-driven parsimony
  fallbacks, O:E/O:EDist, SSDist/OUIndex, Conflict, SetDistT/P, final-trim
  penalties, and packed-WASM dMax with the source 60% ambiguity rule. Add the
  optional logistic/neural selectors and calibrate all strategies against an
  authorized desktop corpus.
- Grouping/ungrouping, recombinant/parent reassignment, unresolved-sequence
  rescans, undo/redo, per-event history, an immutable project ledger,
  duplication, deletion, IndexedDB autosave, and partial-scan checkpoints are
  operational. The ordered Global reconstruction queue, stale downstream
  propagation, mosaic map, possible-overprint links, and recombinant-parent
  dependencies are operational. A tunable heuristic now resolves fresh events
  in order, holds ambiguous dependent branches, and schedules targeted/adaptive
  rescans. Recursive signal erasure, tract-component extraction, crossing-signal
  splitting, and structural breakpoint uncertainty are operational; calibrate
  iteration order and per-method uncertainty windows against authorized RDP5
  outputs and add richer ancestral-event propagation semantics.
- Dedicated tract/combined-background and optional flank/tract/right-flank NJ comparison, six-tree seeded bootstrap/collapse clustering evidence, nearest-context cohorts,
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
- Numerical comparison ledger against RDP5 outputs supplied by the original
  authors and collaborators; the permitted desktop source is the reference
  specification for browser ports.
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
