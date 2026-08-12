# Changelog

## 0.9.5 — 2026-08-12

- Re-enabled BootScan/RecScan production discovery with a direct port of the
  author-supplied RDP5 default distance-mode path: `BSXoverR2`, `SEQBOOT2`,
  `FastBootDistIP`, `SingleToInt`, `GetPltVal2`, the source cutoff/overlap run
  logic, `BSSubSeq`, `MakeScoresBS`, and `ProbCalc`.
- Removed the retired triplet-at-a-time/reference-style behavior. Full mode
  now constructs every unordered concrete triplet exactly once, generates one
  deterministic MSVCRT bootstrap table, computes each requested sequence pair
  once per window, and reuses that row across every containing triplet. All
  three pair relationships are evaluated together; no sequence outside the
  triplet can act as a proxy.
- Added an optimized packed WASM batch with 16-bit quantized JC distances,
  sparse pair compaction for preview/query scans, a triangular distance lookup,
  mismatch-only bootstrap accumulation for complete high-identity windows, and
  an exact missing-data fallback. Whole-alignment baseline ordering uses the
  source `Distance` identity fractions, so JC saturation cannot scramble a
  divergent triplet's role assignment. The reproducible 24 × 2,000 all-triplet
  benchmark reports 2,024 triplets, 276 pairs, 102 windows × 100 replicates in
  307.01 ms, a 22× pair-row reuse factor and 0.71 MiB workspace.
- BootScan can independently create hypotheses, confirm co-located hypotheses
  from other detectors, and recalculate edited events. Projects retain source
  topology/baseline, bootstrap support and replicate count, run windows,
  triplet-local tract/background counts, window/step, raw probability,
  truncation count, batch calls, triplets, reused pairs and workspace bytes.
- Added source controls for window, step, topology cutoff and raw-run retention;
  raised the default retained-run buffer to 20,000 while keeping overflow
  explicit. The method result ledger names the source routine/calibration.
- Added byte/packed, deterministic-seed, unrelated-decoy, multi-triplet,
  divergent-baseline,
  independent worker discovery, recalculation and project-round-trip tests,
  plus `npm run bench:bootscan`.
- Optional desktop `BSTypeFlag` UPGMA/NJ relationship transformations and a
  broad RDP5 executable golden corpus remain explicit parity gates; this release
  does not claim those modes.

## 0.9.4 — 2026-08-12

- Re-enabled GENECONV with a direct WebAssembly port of the author-supplied
  RDP5 batch path. Each unordered concrete triplet is compressed locally,
  invariant and incomplete columns are removed, all-different columns are
  retained, and all six source tracks are evaluated together: the three
  pair-identity tracks and their three complementary outer tracks.
- Ported the source default finite-G fragment rules, mismatch penalties,
  `CalcKMaxP` lambda/K calibration, `GCCalcPValP` probability conversion,
  global cross-track p-value ordering, and source overlap deletion. Every
  retained row carries its actual track, target/minor/major slots, compressed
  ranks, nucleotide interval, score, match/mismatch counts, penalty, raw p,
  and routine provenance.
- Added independent GENECONV discovery, co-located method confirmation,
  exact-hypothesis recalculation, project/CSV round-trip support, per-triplet
  retention controls, truncation accounting, and unordered-triplet kernel-call
  provenance. No designated reference sequence or alignment-wide proxy is
  involved.
- Replaced the source routine's quadratic positive-run extension with an exact
  right-to-left monotone excursion index. The 100 × 10,000 benchmark processes
  2,632 complete six-track GENECONV triplets in about 0.36 s (about 74 million
  full triplet-sites/s) while preserving the native latest-maximum tie rule.
  Scratch memory for that index was reduced from roughly 64 to 16 bytes per
  alignment column.
- Added byte-versus-packed oracle tests, decoy independence, all-different-site
  retention, exact six-track role mapping, finite-G bridging, worker-level
  independent discovery, recalculation, diagnostics, and performance gates.
  A broad desktop-generated golden corpus and non-default indel/permutation
  modes remain validation/feature gates.

## 0.9.3 — 2026-08-12

- Made automated analysis follow the RDP5 manual §4.1.6 cycle by default:
  screen the intact alignment, select the strongest supported signal,
  identify its recombinant/co-recombinant lineage, erase the tract from each
  remainder, append gap-padded tract components, and continue until no
  supported signal remains. These internal cycle decisions remain
  `unreviewed` in the analyst-facing ledger; they never masquerade as manual
  acceptance.
- Added a source-style redo queue. Signals whose concrete triplets do not
  contain a sequence changed by the latest split stay in the detection pool;
  only triplets containing an affected origin are rescanned. Before an
  unaffected pooled signal is applied, its triplets are refreshed against the
  current component alignment so recombinant identification and
  co-recombinant grouping cannot remain stale.
- Added persisted cycle provenance: characterization order, component
  lineage, pass count, initial versus redo comparisons, stop reason, and a
  tunable 1–1,000 pathological-loop safety cap. Both erased remainders and
  extracted tracts remain eligible as separate evolutionary histories.

- Replaced the remaining generic MAXCHI/CHIMAERA discovery shortcut with a
  WebAssembly port of the author-supplied RDP5 scan path. MAXCHI now constructs
  all three pair-equality tracks from concrete triplets, including
  all-different sites; CHIMAERA constructs each of the three possible
  recombinant-oriented binary parent-match tracks and excludes all-different
  sites like the desktop `FSSRDP` path.
- Ported source half-window selection, adjacent 2×H χ² windows, the critical
  count-difference gate, missing-span/end-window bans, circular 11-position
  smoothing, strict peak selection, source-style peak-basin destruction, and
  `GrowMChiWinP`/`GrowMChiWinP2` symmetric peak growth with its failure ceiling.
- Added deterministic paired-peak queues so multiple disjoint MAXCHI and
  CHIMAERA tracts from one concrete sequence triplet become separate candidate
  hypotheses. A tunable 1–256 retention ceiling and omitted-signal provenance
  keep pathological scans bounded without silently hiding truncation.
- Cached the source profiles once per unordered triplet and, in exhaustive
  exploratory mode, emits all eligible recombinant polarities from that one
  scan. The same single-pass optimization now applies to the order-invariant
  RDP detector, eliminating three redundant detector scans per unordered
  concrete triplet while preserving the C(N,3) multiplicity family.
- Replaced the outer target-oriented scheduler with an explicit `a < b < c`
  source batch in role-agnostic mode. Kernel-call provenance and the C(5,3)
  regression now prove that RDP and the combined χ detector each run exactly
  once per unordered triplet—not once per presumed recombinant.
- Added production two-bit triplet extraction for RDP, MAXCHI and CHIMAERA.
  Sixteen alignment columns are decoded per word, only the current triplet's
  informative sites are materialized, and packed-versus-byte oracle tests are
  exact across non-word-aligned lengths and missing data.
- Added explicit decoy-sequence regressions proving that mutations outside a
  concrete triplet cannot leak alignment-global variable columns into its RDP,
  MAXCHI, or CHIMAERA compressed streams.
- Fixed the page-level default that had incorrectly selected query/reference
  mode. New analyses and every bundled example now start with role-agnostic
  all-vs-all enumeration; the targeted mode is explicitly labelled non-parity.
- Removed the former non-source GENECONV, BootScan, SiScan-prelocator and 3Seq
  discovery paths from production. GENECONV, BootScan and 3Seq are disabled
  until their full author-source batch routines are ported; source SiScan now
  confirms every retained source-detector candidate without a generic locator.
- Final evidence and edited-event recalculation now use the retained source
  peak basin rather than the legacy generic χ² locator. `ChiPVal2P` calibration
  receives the exact compressed-site count and selected source half-window.
- Persisted track, target orientation, compressed-site count, half-window,
  both grown boundary statistics and ranks, growth widths, direction, and the
  complete source routine chain in project JSON and CSV. The method panel and
  event dossier expose that ledger directly.
- Added exact two-tract, all-different-site, deterministic queue, worker
  provenance, recalculation, project-round-trip, and UI regressions. These are
  routine-level source fixtures; a broad desktop-executable golden corpus,
  especially missing-data/circular/tie edge cases, remains a validation gate.

## 0.9.2 — 2026-08-11

- Changed the scientific default from an eight-parent shortlist to exhaustive
  enumeration of every unordered combination of three concrete alignment
  sequences. Query/reference mode likewise tests every allowed query with
  every concrete reference pair. The worker records `all-concrete-triplets`
  and `concreteTripletInputs: true` in project provenance.
- Kept distance-pruned screening only as an explicitly named approximate
  preview. Its toggle and settings warn that unlisted triples are not tested;
  it never constructs a consensus/rest-of-alignment proxy, because even the
  preview passes three real sequence indexes to every method kernel.
- Added a concrete-triplet ledger to every method result, made the CLI
  exhaustive by default, and added a C(5,3)=10 worker regression so future
  optimizations cannot silently restore pair/consensus-proxy screening.
- Replaced the simplified BURT post-fit path with direct translations of the
  supplied RDP5 `BenHMM`, `MatchBPtoCI`, and non-reassortment `PolishBP`
  behavior. Switches are now enumerated at the same shifted informative-site
  indexes, and midpoint assignments use VB6 round-half-to-even semantics.
- Reproduced the desktop circular BURT working sequence, including its rotated
  half-copy, zero sentinel, central copy, trailing half-copy, and post-fit crop;
  circular origin-spanning tracts now have deterministic source-mode fixtures.
- Implemented the source `.995`/`.999` any-state confidence scans, wrapped
  confidence intervals, strict first-on-tie variable-site matching, negative
  outside-CI signals, half-tract adoption rule, same-switch conflict handling,
  three-inside/three-outside reversion guard, and missing-data-edge snapping.
- Added an always-available interactive BURT posterior plot to each selected
  event. It overlays the recombinant tract, every HMM state posterior,
  candidate breakpoints, source switches, and 95/99 intervals, and lets an
  analyst apply any switch as the left or right breakpoint.
- Persisted and exported the complete breakpoint audit: candidate and polished
  coordinates, selected switches, winning restart, forward and Viterbi
  likelihoods, source routine chain, state/category mapping, circular-padding
  dimensions, and every `PolishBP` decision/reversion reason.
- Added source-derived regressions for the circular sentinel layout, VB6 ties,
  `MatchBPtoCI` priority/ties/wrapped intervals, origin-spanning polishing,
  same-switch resolution, missing-data edges, insufficient-information
  reversion, deterministic fitting, project round trips, and UI wiring.

## 0.9.1 — 2026-08-11

- Replaced SiScan's final category-Z surrogate with the author-supplied RDP5
  Sister-Scanning workflow: `GetSSOL` nearest/tree/direct outgroup selection,
  the source 15 quartet-pattern categories and topology sums,
  `DoPerms3P`/`MakeZValue2` vertical randomization, `FindMaxZ` topology runs,
  `ShrinkRegionC` boundary refinement, and the desktop `NormalZ × L/region`
  probability convention.
- Added tunable nearest, most-divergent, horizontal-randomization, and explicit
  analyst-selected fourth-sequence modes; triplet/quartet/all-position modes;
  gap stripping/fifth-state modes; separate discovery/final permutation counts;
  and full settings/project sanitization.
- Preserved every significant source SiScan topology run for an ordered
  triplet as a distinct event rather than retaining only the global optimum.
  The fast WebAssembly category run is now only a candidate locator and can
  never by itself count as SiScan support.
- Added deterministic source provenance to every SiScan interval: fourth
  sequence, selection path, pattern/sum family, topology change, scan/final
  permutations, source routine chain, calibrated Z, raw probability, and a
  bounded saved topology-window trace used by the interactive method plot.
- Made the source workflow practical for long genomes: large `MakeVRand`
  tables are regenerated from the identical MSVC stream instead of
  materialized, fixed-outgroup window counts roll incrementally, and exact
  vertical-permutation prefix ranges are cached across triplets. The 80 kb
  source benchmark recovers a 15 kb tract in about 270 ms in the development
  container while avoiding an approximately 80 MB random table.
- Added direct quartet-table, outgroup, `NormalZ`, cached/streamed RNG
  equivalence, exact synthetic-breakpoint, worker provenance, and two-signal
  event-queue regressions, plus a long-genome performance gate.
- Fixed method toggling so the minimum-support threshold is reduced with the
  enabled-method count instead of leaving a hidden impossible threshold.
- Replaced the PHI-labelled proximity surrogate with direct translations of
  RDP5's `PHITest2`, `PHI`, `pair_score`, `GetFandG`, and
  `AnalyticMeanVariance`: parsimony-informative sites now feed the source
  multistate reticulation graph score and analytic normal tail. Large browser
  jobs use a deterministic position-balanced site ceiling and record both the
  retained and all-site counts instead of silently claiming an all-site test.

## 0.9.0 — 2026-08-11

- Added the RDP5 §4.1.6 working-alignment cycle: accepted co-recombinant tracts
  are erased, copied into gap-padded component sequences, recursively split by
  lineage for nested events, and rescanned. Signals crossing deleted sequence
  are split into continuous pieces, and breakpoints near structural gaps are
  persisted and displayed as uncertain.
- Replaced the clustering topology proxy with six seeded Jukes–Cantor
  neighbor-joining bootstrap trees per event, source-default 50% branch
  collapse, source-style tree-position scoring, deterministic large-cohort
  pruning, and tunable exact-site/balanced-block bootstrap controls.
- Extended the direct RDP5 detector to retain multiple distinct raw excursions
  per triplet, with a tunable 1–256 signal ceiling and explicit truncation
  provenance, while preserving the compatibility export for the strongest
  signal.
- Enforced distinct reference-group parents in grouped query/reference scans
  and moved rescanned event clustering onto the private disassembled component
  alignment before mapping all roles back to original user sequences.
- Ported the active RDP5 variable-nucleotide-position triplet detector from the
  author-supplied source, including exact two-equal/one-different categories,
  circular `2h+1` windows, pair ranking, medium-pair dominance,
  `FindNextP`/`DefineEventP2` delineation, role polarity, and source probability
  scaling.
- Added source-mode BURT breakpoint refinement with the desktop fixed-three-
  state model, 21 random starts, MSVC-compatible random stream, 100-iteration
  ceiling, historical backtrace/sentinel behavior, and 95/99% confidence
  intervals. The manual's tunable 2–20-state step-up mode remains available.
- Ported the RDP5 `GetSupers` event-similarity merge and the principal
  `MakeSDMP2`/`FillRmat`/`CalCR` co-recombinant logic: every retained event now
  tests every other sequence under all three presumed-recombinant orientations,
  applies the source six-cell distance correlations and SDM plausibility
  filter, combines phylogenetic, distance, and detectable-signal evidence, and
  retains sequences passing the configurable 2-of-3 rule even when they had no
  independently retained signal.
- Added an exact packed streaming maximum-distance kernel so source-normalized
  event clustering remains memory-bounded on browser-scale alignments.
- Stored method-specific full-alignment intervals for all seven methods and
  exposed their ownership in plots, the inspector, project JSON, and CSV.
- Removed the shared non-RDP candidate seed: every primary family now owns its
  locator, directional families scan both parent orientations, and a method can
  confirm an event only when its interval is co-located. Multiplicity now counts
  unique unordered exploratory triplets and caches order-invariant RDP results.
- Ported finite-G GENECONV fragment scoring plus `CalcKMaxP`/`GCCalcPValP`, and
  moved MaxChi/Chimaera discovery into compressed informative-site coordinates
  with the source probability multiplier.
- Added a persisted RDP5 recombinant-identification ledger using direct
  `MakePhPrScore`, `MakeTrpGroups`/`MakeTrpScore`, `MakeINList`/
  `MakeOUCheck`, `SimpleDist`, `MakeSSDistB`, `GetBadDists`, `FindSets`,
  `MakeEList`/`MakeListCorr`, `MakeLDist`/`MakeRCompat`, packed-WASM
  `CalcMaxD`/`CMaxD2P3`, and default `MakeConsensusC` behavior. All 18
  standalone source statistics, the final-trim penalty, six joint rules, VB6
  banker’s-rounding gates, four-stage parsimony cascade,
  historical-event set closure, PS gates, and final-trim penalties are shown
  per polarity with the source 60% ambiguity rule; the auto-resolver has
  tunable role-confidence gates and weights.
- Added the desktop `GetWinPPfromDists`/TBreak tie path and a production-cohort
  VisRD performance gate; the packed kernel evaluates a 30-taxon, 10 kb role
  cohort in the release benchmark instead of leaving role latency unmeasured.
- Added source-parity regressions for detector boundaries and probability,
  BURT fitting/backtrace behavior, co-recombinant orientation screening,
  unsignalled descendants, clustering, circular events, project round trips,
  source recombinant-role identification, unique-triplet comparison scope,
  independent method discovery, and the large-alignment worker path.

## 0.8.0 — 2026-08-11

- Audited the complete 52-page RDP5 manual plus the author-supplied Visual
  Basic/native-source archive with the RDP authors' permission, and added a
  source/manual-to-web workflow crosswalk with explicit parity boundaries.
- Replaced the primary event stack with an integrated Review studio organized
  around the desktop program's actual refinement sequence: signal,
  breakpoints, roles, trees, co-recombinant grouping, decision, and rescan.
- Added a bounded ordered queue with wrapping previous/next navigation,
  best-unresolved selection, skip-accepted/rejected switches, minimum-method,
  warning-only, and stale-only filters.
- Added a selected-event dossier and review checklist for method confirmation,
  breakpoint precision/freshness, role polarity, false-positive challenges,
  and common-ancestor grouping; all-method and best-evidence modes now share
  the same review surface.
- Added a transparent three-polarity role-assignment workbench that auditions
  each triplet member as the recombinant, exposes tract/background identities,
  switch signal and informative-site count, and clearly labels the score as a
  fast diagnostic rather than RDP5 classifier parity.
- Added distinct per-sequence and whole-event-group decisions plus one-action
  shared breakpoint transfer. Group edits are audited, undoable, and mark
  affected downstream hypotheses stale.
- Changed the local-tree default to recombinant tract versus combined
  background, retained optional left/tract/right comparison, and clarified
  that the connected rectangular NJ drawings are unrooted and
  midpoint-oriented only for readability.
- Added queue, grouping, checklist, polarity, workflow-wiring, and tree-guidance
  regression coverage.

## 0.7.0 — 2026-08-11

- Added an explainable ordered auto-resolver to Global reconstruction with
  conservative, balanced, and aggressive profiles plus a live dry-run of
  accept, reject, analyst-review, locked-decision, and next-rescan outcomes.
- Combined method concordance, adjusted significance, informative-site depth,
  breakpoint precision, parent-conflict, rate-density, diffuse-incompatibility,
  and warning evidence into a tunable score guarded by independent hard gates.
  Stale or uncalibrated evidence can never be auto-accepted.
- Added precise advanced controls for every evidence gate and score weight,
  whether reviewed decisions may be revisited, and all rescan thresholds and
  dependency-risk contributions.
- Added dependency-aware rescan barriers for overlapping tracts in one
  recombinant, recombinant-parent use, reciprocal nested-parent use, and
  co-recombinant groups. Unresolved events transitively block only their
  causally linked downstream branch rather than freezing independent events.
- Added off, impacted-target, adaptive, and full unresolved-target rescan
  strategies. Each rescan rebuilds the remaining plan, respects a tunable round
  cap, preserves resolved/provenance records, excludes accepted mosaic
  sequences from the parent pool, suppresses only close same-parent duplicate
  detections, and is reversible to the pre-run queue state through Undo.
- Indexed event dependencies and genomic-interval buckets so a 5,000-event
  dry-run remains interactive; the development-container stress check completed
  sparse and concentrated queues in roughly 64 ms and 25 ms, respectively.

## 0.6.1 — 2026-08-11

- Closed the main-page paint leak by making every panel a clipping and stacking
  boundary while leaving its content height unconstrained. The recombination
  map now has an explicit, paint-contained scroll viewport that expands into
  the panel body in full-screen mode.
- Replaced the incorrect event-by-event “breakpoint matrix” with a dense,
  interactive genome-position breakpoint-pair density matrix matching the
  semantics of RDP4 Figure 2c.
- Added a split-triangle recombination-region matrix: observed counts of events
  separating pairs of genomic windows above the diagonal, and signed residuals
  from a circular random-tract placement null below it.
- Added a split-triangle local distance-profile discordance matrix as a fast,
  explicitly labeled browser proxy for the SH/RF-style compatibility views in
  RDP4 Figure 2e.
- Reworked matrix rendering around canvas, adaptive 48/64/96-bin resolution,
  keyboard/pointer inspection, exact coordinate readouts, event selection,
  compact legends, and perceptually ordered sequential/diverging palettes.
- Increased the sequence-distance view from 24 to as many as 64 sequences when
  the worker has already calculated the matrix, while keeping each cell compact
  instead of stretching a sparse grid across the panel.

## 0.6.0 — 2026-08-11

- Replaced competing nested panel caps with a single desktop workspace scroll
  owner; panels now grow naturally, while intentional data grids retain local
  scrolling. Every titled panel can expand to an accessible full-screen view.
- Added a multi-parent alignment highlighter with arbitrary 2–6 parent sets,
  unique/shared/novel/missing affinity states, informative-site filtering,
  density controls, sticky labels, and bounded large-alignment behavior.
- Added a Global reconstruction workspace implementing an ordered RDP-style
  review queue, downstream stale propagation, rescan prompts, collection-wide
  mosaic map, event grouping, possible overprinting, recombinant-parent
  dependencies, and direct alignment/tree verification links.
- Replaced tree line fragments with a connected cumulative-distance SVG layout,
  explicit joints and zero-length nodes; added linked leaf marking across all
  regional trees and layout regression tests.
- Eliminated ambient locale/time-zone text from server rendering, fixed the
  deterministic tutorial timestamp, and added byte-identical double-SSR and
  hydration-risk regressions for React error #418.
- Expanded the five-step guide, scientific ledger, feature audit, and tutorial
  references to describe ordered review, rescan, nesting, and parent-proxy
  limitations explicitly.

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
