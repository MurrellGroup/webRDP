# RDP Web replacement roadmap

RDP Web will not call itself a validated RDP5 replacement merely because the
interface and method names exist. The replacement claim is gated by fixtures,
simulation results, workflow parity, and browser-scale performance.

## Track A — primary detection methods

- Manual §4.1.6 cyclical detection is operational and default: best-signal
  selection, co-recombinant erase/extract, recursively addressed components,
  unaffected-signal pooling, affected-origin redo scans, and current-component
  refresh before application. Finish desktop golden fixtures for exact
  redo-list compaction, minimum-size component dropping, and tie ordering.
- The parity scan now constructs one `a < b < c` record for every unordered
  concrete triplet. RDP, six-track GENECONV and combined MAXCHI/CHIMAERA each
  run once and resolve internal tracks/polarities; two-bit extraction skips
  invariant columns for that triplet only. Keep the exact packed/byte and
  kernel-call regressions.
- RDP: reproduce reference selection, informative-window enumeration, and
  local/global p-values on documented triplets. Multi-excursion source
  enumeration is operational; complete the desktop tie/missing/circular corpus.
- GENECONV: the default `AlistGC2`/`GCXoverD` six-track batch, fragment queue,
  `CalcKMaxP`/`GCCalcPValP`, overlap deletion, independent discovery and exact
  packed/byte paths are operational. Finish desktop golden fixtures and the
  non-default indel, overlap-count and permutation modes.
- BOOTSCAN/RECSCAN: the author-source default distance batch is operational and
  independently discovers events. It shares one `SEQBOOT2` table and every
  requested pair/window `FastBootDistIP` row across the unordered triplet set,
  applies `GetPltVal2` topology support and triplet-local
  `MakeScoresBS`/`ProbCalc` evidence, and supports exact recalculation. Finish
  optional `BSTypeFlag` UPGMA/NJ relationship modes, desktop warning/tie edges,
  circular boundary goldens and the broad executable corpus.
- MAXCHI and CHIMAERA: source compressed tracks, window rules, 11-position
  smoothing, basin destruction, GrowMChiWin expansion and multi-peak queues are
  operational; finish desktop lookup-table/rounding and edge-case corpus.
- SISCAN: source nearest-tree/direct, most-divergent, manual and randomized
  fourth-sequence paths; 15 site categories and topology sums;
  horizontal/vertical randomization; reproducible seeds; topology-run
  enumeration and shrinkage; once-per-unordered-triplet independent discovery;
  topology-transition role inference; exact telescoping permutation bands;
  deterministic band-moment reuse; and bounded-memory fast/full controls are
  operational. The 220-triplet source workload is 2.40× faster than 0.9.8 with
  direct-enumeration equality coverage. Finish desktop tie/missing-data and
  plot-by-plot corpus parity, then pursue larger-cohort worker parallelism/SIMD.
- 3SEQ: the fused all-three-role `FindSubSeqTS` triplet walk, strict excursion
  selection, `Seq3PVals`-equivalent exact probabilities, `GetTSPVal` large-walk
  path, source `CheckwrapC` origin/linear-complement handling,
  `CheckSplit3Seq`/`SubPVal` continuous-piece recalibration, independent
  discovery and edited-event recalculation are operational. Finish desktop
  coordinate/tie output fixtures, precomputed-table rounding comparison, and
  the broad executable golden corpus.

## Track B — false-positive control and breakpoint inference

- Source-compatible fixed-three-state BURT fitting, repeated random starts,
  historical backtrace/sentinel behavior, exact circular half-copy padding and
  cropping, `.995`/`.999` confidence scans, `MatchBPtoCI` VNP matching,
  non-reassortment `PolishBP` adoption/reversion rules, missing-data snapping,
  and an interactive posterior/switch workbench are operational with
  source-derived circular, tie, same-switch, missing-data and information-guard
  fixtures. Remaining gate: compare these deterministic outputs with a broad
  desktop-generated golden corpus and finish optional reassortment segment-
  boundary handling plus the manual 2–20-state mode's selection corpus.
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
  tract/background tree verification. Source-ready families own their interval
  provenance; finish and re-enable the remaining method batches one at a time.
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
- The C(5,3) worker regression guards complete concrete-triplet coverage and
  proves one RDP, one SiScan, one fused three-role 3Seq and one combined χ
  source call per triple whenever those families are enabled. BootScan records
  its batch calls, concrete triplets, unique reused pairs, windows, replicates
  and workspace bytes.
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
