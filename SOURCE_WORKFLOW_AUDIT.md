# RDP5 manual/source workflow audit

## Scope and source-use rule

This audit records product behavior and active algorithm paths from the
author-supplied RDP5 manual and source archive. The original authors have
permitted direct use for this MIT project. The web implementation therefore
ports selected routines where compatibility matters while replacing the VB6 /
Win32 state architecture with typed project state, workers, and WebAssembly.
The legacy projects are read-only behavioral specifications: they are never
compiled, linked, or included in the web build.

Materials reviewed for this checkpoint:

- `RDP5Manual.pdf`, 52 pages, especially sections 4.1, 5.1–5.5, 8, 9, and
  10.3–10.4.
- `VB Source.zip`, including the 62,622-line `MainForm22.frm`,
  `DendrogramForm6.frm`, `OptionsForm2.frm`, `ManualSeqSelect.frm`, and the
  associated VB modules.
- `dnaDLLSource.zip` and `dna5DLLSource.zip`, traced into active exported
  operations for source-compatible compute ports.

The legacy archives themselves are not redistributed. Routine-level mappings,
known divergences, and regression fixtures are recorded in `SOURCE_PARITY.md`.

## The central workflow recovered from RDP5

The desktop product is not fundamentally a report viewer. It is an ordered,
editable hypothesis-reconstruction system:

1. Load and classify an alignment, including linear/circular structure and
   query/reference roles.
2. Detect preliminary signals with enabled exploratory methods.
3. Reconcile repeated triplet detections into putative unique events.
4. Characterize the strongest/earliest event first.
5. Confirm method support and challenge artefacts.
6. Refine both breakpoints.
7. test which triplet member is recombinant and which sequences are the best
   available major/minor parent proxies.
8. Compare local trees for the recombinant tract and the complementary
   background; track sequences across trees.
9. Group descendants of the same ancestral recombinant and handle overprinted
   events.
10. Accept or reject one occurrence or the full event group.
11. Rescan after a changed upstream hypothesis, because the original signal
    removal/grouping decision may have changed later detections.
12. Continue through the reconstruction queue, then export a project or a
    recombination-aware alignment.

The manual repeatedly describes automated output as a fallible preliminary
account. Its key operational warning is that an incorrect early event can
propagate mistakes into all subsequent characterization. That dependency—not
the visual styling of the old Windows forms—is the organizing principle of the
new Review studio.

## Observable behavior crosswalk

| RDP5 evidence | Observable behavior | RDP Web 0.9.5 implementation | Fidelity boundary |
| --- | --- | --- | --- |
| Manual §4.1 and §10.3 | Detection and unique-event inference are separate phases, and detection screens real sequence triplets. | Full mode emits one `a < b < c` record per concrete triple. RDP, six-track GENECONV and combined MAXCHI/CHIMAERA each receive it once. BootScan consumes the same set in one shared-pair `BSXoverR2`-shape batch, evaluates all three relationships together, and never substitutes a consensus/reference proxy. Packed extraction and triplet-local scoring ignore only the current triplet's invariant/incomplete sites where the source does. | BootScan optional UPGMA/NJ relationship modes and 3Seq remain pending; source SiScan presently confirms source-detector candidates rather than owning the initial batch. |
| Manual §5.1 and §10.4; `GoToNextEventMnu_Click` | Events are reviewed in stored characterization order; next/previous wrap and can skip accepted/rejected events. | Ordered queue, previous/next wrap, best-unresolved jump, skip accepted, skip rejected, warning-only, stale-only, and minimum-method filters. | RDP5's exact internal `BestEvent` ordering has not been numerically replicated. |
| `GoToBestMnu`, `Detect1Mnu`–`Detect7Mnu` | Jump to strongest unaccepted evidence and filter by number of detecting methods. | Best adjusted-p unresolved jump and 0–7 minimum-support filter. | Adjusted p is used as the transparent tie-breaker; exact desktop ordering remains a fixture target. |
| Manual §5.1–5.3; `ShowBestAllMnu` / `ShowAllAllMnu` | Schematic context, event information, and method plots work together; best evidence and all evidence can be selected. | Integrated dossier, recombination map, all-method tabs, and best-evidence mode in one continuous review surface. | The original method-specific graphic primitives are not pixel-cloned. |
| Manual §5.2 | Event information includes sequence roles, breakpoints, confidence, warnings, and recombinant-identification confidence. | Selected dossier plus checklist for method confirmation, breakpoint precision/freshness, role warnings, artefact warnings, and co-recombinant grouping. | Checklist rules are documented heuristics rather than the hidden RDP5 score. |
| Manual §8/§9; `PHITest2` and PHIPACK-derived helpers | PHI distinguishes spatially clustered incompatibilities from diffuse homoplasy/error and is reported as challenge evidence. | Direct multistate `pair_score` graph, PHI neighbourhood statistic, analytic mean/variance and normal lower-tail calculation appear beside a separate four-gamete permutation diagnostic. Site work ceilings and total informative-site counts are explicit. | Exact VB random thinning above 6,000 sites and desktop golden-corpus calibration remain. |
| Manual §5.3–5.4; `BenHMM`, `MatchBPtoCI`, `PolishBP` | Plot and sequence displays support breakpoint placement, source breakpoint refinement, and parent-match inspection. | Draggable/numeric breakpoints, confidence reset, stale-evidence marking, arbitrary 2–6-parent alignment affinity highlighting, and an interactive BURT posterior/state/switch plot. Source mode implements circular padding/crop, `.995`/`.999` confidence scans, VNP-space switch matching, and principal non-reassortment polish/reversion behavior. | Reassortment segment-boundary branches, exact RDP5 method-plot statistics, and desktop BURT golden outputs remain. |
| Manual §4.1.5 and §10.4; source role-classifier routines | Each member of a detection triplet must be considered as the recombinant; a weighted profile/tree consensus identifies polarity and roles remain editable. | Three-way workbench persists all 18 standalone source statistics, the final-trim penalty, and six joint rules: profile/tree movement, TrpScore/OuCheck, O:E/O:EDist, SSDist/OUIndex, four-stage ParsimonyO/I with `FindSets` event-history closure, distinct-category Conflict, SetDistT/P, packed-WASM dMax, VB6 banker’s-rounding gates, and the desktop default `MakeConsensusC` weighting/60% ambiguity rule. Auto-resolve holds challenged/ambiguous roles. | Optional logistic/neural selector modes and desktop golden-corpus calibration remain. |
| Manual §5.5; `DendrogramForm6.frm` | Paired trees show recombinant region versus non-recombinant region; marking follows a sequence across trees. | Default tract-versus-combined-background trees, optional left/tract/right mode, linked leaf marking, role colors, nearest-sequence cohorts, and Newick export. | Browser trees are NJ/p-distance aids. ML/Bayesian trees and SH/AU tests remain open. |
| Manual §5.5 | Drawn trees are unrooted even when visually oriented; absence of a topology change is not proof against recombination. | Tree caveat explicitly distinguishes midpoint-oriented drawing from an inferred root and warns against interpreting left-to-right ancestry. | No root inference is claimed. |
| Node menus in `DendrogramForm6.frm` | Mark/unmark clades, find parent candidates, and accept/reject events above a node. | Linked leaf selection and event-group decisions are operational; event-centered role reassignment is directly accessible in the inspector. | Editable internal-node/clade operations and parent search above a node remain open. |
| Manual §4.1.3–4.1.4 and §10.4 | Descendants of one ancestral recombinant are grouped; overprinting can obscure group membership. | Every retained event tests every other sequence under all three presumed-recombinant orientations. Six seeded JC/NJ bootstrap trees, source-default 50% branch collapse, source-style tree movement, distance correlation/SDM, and detectable-signal evidence feed the tunable 2-of-3 rule; qualifying unsignalled descendants enter `GetSupers`-style ancestral groups. Group decisions and dependency reconstruction are operational. | Exact Clearcut consensus and the complete `CheckBSTree` weighting tree need desktop-corpus parity; ancestral-sequence reconstruction remains open. |
| `AcceptMnu_Click` / `AcceptSMnu_Click` and reject counterparts | A decision can apply to the selected sequence or every sequence carrying the event. | Separate selected-hypothesis and common-ancestor-group accept/reject controls; each is one undoable, audited action. | Group membership is analyst/heuristic supplied rather than guaranteed to match RDP5 inference. |
| `FindBestRecSignal`, `BuildFirstXOList`, redo lists, `CheckDrop`/`DropSeqs`; manual §4.1.6/§10.4 | Automated detection repeatedly selects the best signal, characterizes it, erases/extracts co-recombinant tracts and re-screens; later manual corrections invalidate downstream assumptions. | Default Run performs the full sequential cycle. Unaffected signals remain pooled; changed-origin triplets alone enter the redo scan; a pooled signal is refreshed on the current component alignment before application. Scientific edits still mark later queue items stale and manual/heuristic rescans rebuild the lineage. Crossing signals split and gap-near breakpoints carry structural uncertainty. | Exact desktop redo-list compaction/minimum-size dropping, method-specific VNP uncertainty thresholds and desktop iteration-order golden fixtures remain. |
| `RCheckMnu`, tree recheck commands | Re-run methods or recheck a plot with a different candidate role. | Exact selected-hypothesis recalculation, role audition, method-specific views, and full/targeted unresolved rescan. | Recalculation preserves scan multiplicity conservatively but does not claim RDP5 p-value identity. |
| Matrix menu and manual §5.5 | Compatibility, recombination, breakpoint, LD, MaxChi, and LARD matrices are interactive context. | Dense genome-position breakpoint-pair, region-separation, local-discordance, and p-distance views with coordinate inspection. | Exact RF/SH/LD/MaxChi/LARD matrix implementations remain separate parity tasks. |
| Manual §9 | Save project, remove recombinants, remove/mask tracts, split mosaics, and partition at breakpoints. | Restorable project, CSV, input FASTA, remove, mask, codon-aware mask, split, and breakpoint-partition outputs. | Native `.rdp` project import/export is not implemented. |

## Source architecture observations

The main form is a large stateful coordinator. Menu handlers change shared
event arrays, rebuild event order, refresh schematic/tree/matrix displays, and
then navigate to a new event. The DLLs expose low-level scan, distance, tree,
matrix, grouping, HMM, and display-support operations. This reinforces two web
architecture choices:

- heavy calculation belongs in the worker/WebAssembly layer, while review state
  and provenance remain in an explicit project model;
- a scientific edit is a transaction affecting the reconstruction graph, not a
  local cosmetic change to one card.

RDP Web preserves those useful semantics without preserving the original
global-variable architecture. Group decisions and shared-breakpoint edits are
single undo frames, every action receives an audit entry, and stale downstream
evidence is visible rather than silently recomputed.

## Implemented through the 0.9.5 source-parity checkpoint

- Review studio replaces the disconnected main-page result stack.
- Six-stage selected-event rail: signal, breakpoints, roles, trees, grouping,
  decision.
- Ordered, filterable, wrapping reconstruction queue with best-unresolved jump.
- Compact event dossier and explicit five-part review checklist.
- Best-evidence/all-method confirmation modes.
- Three-polarity source role ledger with direct default-`MakeConsensusC`
  components, historical-event set closure, packed quartet dMax,
  transparent component weights, and
  ambiguity-aware auto-resolution.
- Individual versus whole-group decisions and group-wide breakpoint transfer.
- Tract-versus-combined-background tree comparison as the default, with an
  optional three-region view and explicit unrooted-tree guidance.
- Regression tests for queue navigation, filters, best-event selection,
  grouping, checklist state, role scoring, workflow wiring, and tree wording.
- Direct RDP5 VNP detection and source probability scaling.
- Complete role-agnostic concrete-triplet enumeration by default, one
  `a < b < c` scheduler record and one source-kernel call per family, with all
  internal polarities and no consensus/rest-of-alignment proxy. Packed and byte
  extraction are exact; targeted/approximate modes are opt-in and marked.
- Direct MAXCHI/CHIMAERA compressed tracks, window rules, smoothing, peak-basin
  destruction, GrowMChiWin expansion and multi-peak queues.
- Direct six-track GENECONV compression, finite-G fragment extension,
  lambda/K calibration, global p-value queue, source overlap suppression,
  independent discovery and exact recalculation.
- Direct distance-mode BootScan/RecScan batch with MSVCRT-seeded `SEQBOOT2`,
  pair/window `FastBootDistIP` reuse across the full concrete-triplet set,
  `GetPltVal2` topology support, source cutoff/overlap runs, triplet-local
  `BSSubSeq`/`MakeScoresBS`/`ProbCalc`, independent discovery and exact
  hypothesis recalculation.
- Fixed-three-state source BURT with exact circular working-set construction,
  shifted switch/CI scanning, `MatchBPtoCI`, principal non-reassortment
  `PolishBP`, and an interactive posterior workbench, plus the manual 2–20-state
  mode.
- `GetSupers` ancestral-event merging and all-sequence/all-three-orientation
  co-recombinant screening, including unsignalled descendants.
- Six-tree seeded JC/NJ bootstrap evidence with low-support branch collapse.
- Recursive §4.1.6 signal erasure, gap-padded tract components, component-aware
  rescans, crossing-signal splitting, and structural breakpoint uncertainty.
- Multi-excursion RDP raw-signal retention with a tunable cap and truncation
  audit provenance.
- Source Sister-Scanning with tree/direct/manual/randomized fourth-sequence
  selection, the 15-pattern/sum table, horizontal/vertical randomization,
  topology-run enumeration, region shrinkage, calibrated whole-region Z
  probability, all significant runs, and bounded-memory long-genome execution.

## Remaining workflow gaps

The following are not disguised as complete:

- the optional BootScan UPGMA/NJ relationship modes; the complete disabled
  3Seq source discovery batch; GENECONV's non-default indel/permutation modes;
  independent source SiScan discovery orchestration; and desktop numerical
  parity for every method;
- the optional RDP5 logistic/neural recombinant selectors and desktop
  golden-corpus calibration of the default decision tree;
- editable tree nodes/clades, RF distances, ML/Bayesian inference, SH/AU tests,
  and ancestral-sequence reconstruction;
- exact Clearcut/`CheckBSTree` parity and a desktop event corpus for
  overlapping/overprinted co-recombinant histories;
- native RDP4/RDP5 project conversion;
- per-method deleted-boundary uncertainty parity and desktop iteration-order
  fixtures for the component disassembly engine;
- an authorized RDP5 comparison corpus and cross-browser scientific E2E suite.

Until those gates are met, RDP Web remains a scientific alpha and every inferred
event remains a hypothesis requiring independent review.
