# RDP Web

RDP Web is a local-first, browser-based recombination analysis workbench. It
combines optimized WebAssembly kernels with an interactive workflow for
alignment/project import, method-specific screening, event verification,
breakpoint editing, matrices/topology checks, and recombination-aware exports.
Version 0.6 added a collection-level event reconstruction workspace, arbitrary
multi-parent alignment highlighting, linked and continuously connected local
trees, full-screen analysis panels, deterministic hydration-safe rendering,
and a single reliable workspace scrolling model. Version 0.6.1 also made
panel paint containment structural and adds dense RDP4-style breakpoint-pair,
recombination-region, and local-discordance genome matrices with interactive
coordinate readouts and accessible scientific palettes. Version 0.7 added a
tunable ordered auto-resolver with conservative/balanced/aggressive profiles,
explicit evidence gates, dependency-aware targeted rescans, transitive safety
holds, live dry-runs, rescan-round caps, and one-step workflow undo. Version
0.8 rebuilds the primary analysis surface from a direct audit of the RDP5 user
manual and author-supplied desktop source: a six-stage Review studio, ordered
and filterable reconstruction queue, best-unresolved navigation, integrated
event dossier/checklist, best-versus-all method evidence, three-polarity role
audition, tract-versus-background trees, and separate per-sequence and
whole-group decisions. Version 0.9 advances numerical source compatibility: the
active RDP5 VNP-window detector now retains multiple raw excursions per
triplet, fixed-three-state BURT fitting/backtrace provides posterior breakpoint
intervals, and accepted events drive the manual's erase/extract/rescan
component workflow. Signals crossing deleted tracts are split into continuous
pieces and gap-adjacent breakpoints are marked uncertain. Weighted
ancestral-event merging and the all-sequence/all-three-orientation
co-recombinant screen are ports from the author-supplied desktop source. Its
phylogenetic set now uses six seeded JC/NJ bootstrap trees with low-support
branches collapsed; distance-correlation and detectable-signal sets complete
the configurable 2-of-3 rule, including descendants without their own raw
detection. Version 0.9.1 replaces the remaining SiScan confirmation surrogate
with the supplied RDP5 Sister-Scanning path: source 15-category and sum scores,
tree/direct/manual fourth-sequence selection, horizontal or vertical
randomization, topology-run enumeration, region shrinkage, and whole-region Z
calibration. Every locally significant run becomes its own hypothesis and the
fast WASM category pass is only a locator. Long-genome scans reproduce the
desktop MSVC random stream without materializing its potentially multi-gigabyte
table and reuse exact permutation-prefix ranges across triplets. See
`SOURCE_WORKFLOW_AUDIT.md` for the source-to-web crosswalk. The
same checkpoint replaces the PHI-style warning surrogate with direct ports of
RDP5's multistate `pair_score`, `PHI`, and analytic mean/variance path; any
browser work-ceiling subset is explicitly reported with its all-site count.
Version 0.9.2 closes the previously simplified half of source BURT: exact
circular working-sequence padding/cropping, shifted switch enumeration,
`.995`/`.999` confidence scans, VNP-space `MatchBPtoCI`, and the principal
non-reassortment `PolishBP` adoption, same-switch, missing-data and
information-reversion rules are now active. Every selected event exposes an
interactive posterior/state/switch plot and a persisted source-decision
ledger.
Version 0.9.3 replaces the former MAXCHI/CHIMAERA peak-pair approximation with
the author-supplied RDP5 control flow in WebAssembly: three MAXCHI pair-equality
tracks, three recombinant-oriented CHIMAERA tracks, source compressed-site and
half-window rules, missing/end-window bans, 11-position smoothing, peak-basin
destruction, `GrowMChiWin` expansion, and bounded multi-tract queues. Every
retained call records both grown boundary statistics, ranks, growth widths and
its source routine chain. Exhaustive exploratory scans evaluate each unordered
concrete triplet once for RDP, MAXCHI and CHIMAERA and emit all eligible
polarities from that cached pass, avoiding three redundant detector scans.
The direct source fixtures cover multiple disjoint tracts and the important
MAXCHI-versus-CHIMAERA treatment of all-different sites; executable golden-
corpus validation remains required before a numerical-parity claim. The
production scheduler now mirrors the desktop `AList` batch shape: it walks
only `a < b < c`, invokes RDP once and the combined MAXCHI/CHIMAERA kernel once,
and resolves all internal pair tracks and recombinant orientations inside that
concrete triplet. Sixteen-site two-bit extraction skips invariant columns
without rereading three full sequence strings; byte-oracle regressions require
identical event rows, including missing data. Every detector receives three
explicit sequence records, never a pair plus an alignment-consensus/rest-of-
alignment proxy. The main page no longer overrides this role-agnostic default
with query/reference mode; designated references are unnecessary. Distance-
pruned and query/reference screening survive only as visibly non-parity opt-ins.
The default Run command now also executes the manual §4.1.6 detection cycle:
it selects the strongest remaining signal, splits every inferred
co-recombinant into erased remainder and gap-padded tract components, then
redoes only concrete triplets containing an origin changed by that split.
Unaffected signals stay in the pool, but each is refreshed against the current
component alignment before it is characterized and applied. The cycle stops
only when no supported signal remains (or a disclosed, tunable safety cap is
reached), and its internally applied events remain visibly unreviewed for the
analyst.
Version 0.9.4 re-enables GENECONV with the author-source six-track batch rather
than the retired simplified locator. One packed concrete triplet pass removes
only locally invariant/incomplete columns, retains all-different columns,
builds three pair-identity and three complementary outer tracks, applies the
source integer-G mismatch penalty and `CalcKMaxP`/`GCCalcPValP` calibration,
then drains one globally p-ordered overlap-suppressed fragment queue. The
native quadratic fragment extension is represented exactly by a linear-time
monotone excursion index. GENECONV can now independently create hypotheses,
confirm co-located hypotheses from other methods, and recalculate edited events;
the project, CSV and interface retain its full track/role/score/probability
ledger.
Version 0.9.5 re-enables the default RDP5 distance-mode BootScan/RecScan path.
The worker enumerates each unordered concrete triplet once, generates one
seeded `SEQBOOT2` table, computes every requested sequence-pair distance once
per window, and shares that pair matrix across all triplets. The packed WASM
kernel follows `FastBootDistIP` Jukes–Cantor distances and `GetPltVal2` tie
ordering, while the whole-alignment baseline follows the source `Distance`
identity matrix without JC-saturation ties. It retains supported topology runs and applies triplet-local
`BSSubSeq`/`MakeScoresBS`/`ProbCalc` scoring. It independently creates
hypotheses, confirms co-located signals, recalculates edited events, and
persists topology/support/window/probability provenance. Sparse preview/query
batches compact the requested pair set; complete high-identity windows visit
bootstrap weights only at mismatches. The optional desktop UPGMA/NJ
relationship transformations and broad executable golden corpus remain parity
gates and are disclosed in the interface.

Scientific integrity takes precedence over a seven-method checkbox count.
3Seq is disabled in production until its complete author-source batch path
replaces the retired simplified locator. BootScan/RecScan's author-source
distance mode is active; its optional tree relationship modes remain pending. SiScan is
run through its source confirmation path on candidates from the source
triplet detectors; its former oriented prelocator is retired. Direct
ports of `MakePhPrScore`, `MakeTrpGroups`/`MakeTrpScore`,
`MakeINList`/`MakeOUCheck`, `SimpleDist`, `MakeSSDistB`, `GetBadDists`,
`MakeEList`/`MakeListCorr`, `MakeLDist`/`MakeRCompat`, `CalcMaxD`/
`CMaxD2P3`, and the default `MakeConsensusC` decision tree now persist an
18-statistic recombinant-identification ledger, final-trim penalty, and six
joint rules for all three polarities;
the role workbench exposes every value/weight and auto-resolution holds source-
challenged or sub-60%-confidence assignments. The truth-annotated synthetic
dataset library, method-by-method result interpretation, safe-by-default
fresh/accepted exports, scan-scope event recalculation, and project audit ledger
remain available. Challenge diagnostics, breakpoint matrices, annotation
workflows, and the WebAssembly worker architecture remain available.

The project is MIT licensed. Source-compatible RDP5 routines are included with
the RDP authors' permission and collaboration; the author-supplied RDP5
VB/native source is the sole detector specification. See `SCIENTIFIC_BASIS.md` for the
implementation/validation ledger and `FEATURE_AUDIT.md` for the product-level
gap ledger. Details are in `THIRD_PARTY_NOTICES.md`.

## Local development

```sh
npm install
npm run wasm:build
npm run dev
```

The checked-in `public/wasm/rdp.wasm` is built from `assembly/index.ts`. The
module worker in `public/rdp-worker.js` owns the WebAssembly instance and sends
only compact progress/results messages to the interface. `public/rdp-statistics.js`
contains dependency-free probability calibrations.

## GitHub Pages

The repository exports to static files and is safe at both a domain root and a
project subpath. To verify a project-site build locally:

```sh
RDP_BASE_PATH=/your-repository npm run github:build
```

The artifact is written to `out/`. The included
`.github/workflows/github-pages.yml` rebuilds the WebAssembly binary, runs the
scientific regression suite, creates the static export, and deploys it when the
repository's Pages source is set to **GitHub Actions**. The workflow derives the
project path from `GITHUB_REPOSITORY`; repositories named `*.github.io` build at
the domain root.

Worker and WASM URLs are resolved relative to the document/worker locations,
so a repository name or later custom-domain change does not require source
edits.

## Tests and benchmarks

```sh
npm run test:wasm
npm run test:ui
npm run bench:wasm -- 100 10000
npm run lint
```

## Batch runner

The same WebAssembly worker can be run without a browser and writes a
restorable project:

```sh
npm run cli -- alignment.fasta results.rdpweb \
  --bootstrap 100 --seed 1511506142
```

Use `npm run cli -- --help` for method, circular-genome, query/reference, and
explicit approximate-parent-shortlist options. Full concrete-triplet screening
is the default. This is a local batch interface, not a separate algorithm
implementation.

The tests cover scalar/packed distance equivalence, the RDP5 multi-signal VNP event locator,
the C(5,3) full-concrete-triplet invariant, order-invariant cached triplets,
source-detector discovery with exact current-triplet site compression, co-location
gating, all six GENECONV tracks, finite-G fragment scoring/KA calibration,
global overlap-suppressed source queueing, informative-site
MAXCHI/CHIMAERA windows and source peak correction,
source SiScan quartet categories, outgroup selection, cached/streamed MSVC
randomization, topology-run calibration and multi-event expansion,
source PHI incompatibility graphs, analytic moments and lower-tail calibration,
exact HGRW tails, seeded bootstraps, sampled large-dataset parent selection,
source-mode BURT circular padding, VB6 ties, CI matching, missing-edge and
reversion behavior plus manual step-up BURT localization, recursive
erase/extract signal disassembly, affected-triplet redo scheduling,
current-component refresh, deleted-tract splitting and uncertainty,
RDP5 weighted event merging, bootstrapped-JC/NJ tree evidence, and
three-orientation co-recombinant clustering, direct source PhPr/SubDist/
TrpScore/OuCheck/O:E/dMax/parsimony/conflict/set-distance role identification,
historical-event set closure, origin-spanning circular events, rate-variation and
gap-block challenges, exact manual-event recalculation, NJ/hotspot routines, a
513-sequence large-path integration fixture, annotation-aware masking, and
lossless `.rdpweb` project round-trips, parent-affinity classification,
ordered/nested event reconstruction, connected tree geometry, global scroll
ownership, full-screen panel controls, deterministic server rendering, and
auto-resolution scoring/rescan barriers. The source-guided regressions also
cover filtered/wrapping event navigation, best-event selection, group
membership, checklist freshness, source-weighted role polarity, and tract/background
tree workflow wiring.

## Current status

This is a functional scientific alpha, not yet a validated drop-in replacement
for RDP5. The application itself surfaces that distinction prominently. See
`ROADMAP.md` for the replacement gates rather than treating the current method
names as a parity claim.
