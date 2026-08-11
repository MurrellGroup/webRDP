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
ledger. The
former cross-method CUSUM discovery shortcut has also been removed: every
enabled family now contributes its own interval and a method can confirm an
event only when that interval is co-located with the event. GENECONV uses the
source finite-G mismatch penalty and `CalcKMaxP`/`GCCalcPValP` probability
path (including the G=0 special case); MAXCHI and CHIMAERA scan compressed
informative-site coordinates and use the source peak multiplier. Unique
unordered triplets define both the default exhaustive scan and the multiplicity
family: every method receives three explicit sequence records, never a pair
plus an alignment-consensus/rest-of-alignment proxy. Distance-pruned parent
shortlisting survives only as a visibly approximate opt-in. Order-invariant RDP
triplets are cached to avoid redundant three-polarity source scans. Direct
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
remain available. The exact
bounded 3SEQ calibration, seeded p-distance bootstrapping, challenge
diagnostics, breakpoint matrices, annotation workflows, and WebAssembly worker
architecture remain available.

The project is MIT licensed. Source-compatible RDP5 routines are included with
the RDP authors' permission and collaboration; OpenRDP remains an additional
algorithmic reference. See `SCIENTIFIC_BASIS.md` for the
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
independent non-RDP discovery, co-location
gating, finite-G GENECONV fragment scoring/KA calibration, informative-site
MAXCHI/CHIMAERA windows and source peak correction,
source SiScan quartet categories, outgroup selection, cached/streamed MSVC
randomization, topology-run calibration and multi-event expansion,
source PHI incompatibility graphs, analytic moments and lower-tail calibration,
exact HGRW tails, seeded bootstraps, sampled large-dataset parent selection,
source-mode BURT circular padding, VB6 ties, CI matching, missing-edge and
reversion behavior plus manual step-up BURT localization, recursive
erase/extract signal disassembly, deleted-tract splitting and uncertainty,
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
