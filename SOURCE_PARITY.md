# RDP5 source-parity ledger

This ledger maps author-supplied RDP5 routines to browser ports and tests. A
row marked “ported” means the active control flow and important historical
quirks are represented; it does not replace golden-corpus comparison against
the desktop executable.

| Desktop source path | Browser implementation | State | Regression evidence / remaining edge cases |
| --- | --- | --- | --- |
| `FindSubSeqPB3` / `FindSubSeqPB6` | Full worker triplet enumerator plus `scan_rdp5_triplet` category compression in `assembly/index.ts` | Ported | Every unordered set of three concrete alignment records is scanned by default, every polarity is considered, and the kernel receives three explicit indexes rather than an alignment proxy. C(5,3)=10 orchestration, exact two-equal/one-different VNP categories, and all-three-category regressions pass; desktop missing-data/tie corpus remains. |
| `XOHomologyP2` | `scan_rdp5_triplet` circular 2h+1 sliding pair counts | Ported | Deterministic synthetic mosaic test; origin/tie fixtures remain |
| `FindNextP` + `DefineEventP2` | Medium-pair dominance and common/different event delineation in `scan_rdp5_triplet` / `scan_rdp5_triplet_all` | Ported | Multiple distinct excursions are retained per triplet with strongest-first bounded retention and explicit truncation provenance; polarity/window/tract regressions pass; desktop boundary-by-boundary multiplicity fixtures remain |
| `ProbCalcP` / `ProbCalcP2` | `rdp5SourceProbability` in `public/rdp-statistics.js` | Ported formula | Binomial tail, source `y-1` probability-site multiplier, 169-site scaling and lower clamp; factorial-table rounding fixtures remain |
| `BenHMM` + `DoHMMCyclesSerial` + `GetLaticePathP` + `ForwardCP` + `ReverseCP` | Source mode in `public/rdp-burt.js` | Ported | Fixed 3 states; inclusive 21 starts; shared historical convergence score; 100-iteration ceiling; 0.01 counts; Single transition initialization; rejection-sampled emission anchors; MSVC 15-bit RNG; predecessor-offset backtrace; inclusive zero sentinel; exact circular rotated-half/central/trailing-half working sequence and crop; deterministic linear/circular regressions pass. Desktop golden-corpus comparison remains. |
| `MatchBPtoCI` / non-reassortment `PolishBP` | Source switch/CI/matching/polishing path plus interactive posterior workbench in `public/rdp-burt.js` and `app/page.tsx` | Ported | Shifted switch enumeration; VB6 midpoint ties; `.995`/`.999` any-state scans; wrapped intervals; strict first-on-tie VNP matching; negative outside-CI signal; half-tract adoption; same-switch resolution; three-variable-sites-per-side reversion; missing-data-edge and complete-site snapping; full decision ledger, project/CSV persistence, circular/tie/missing/reversion fixtures. Reassortment segment-boundary branches and executable golden outputs remain. |
| `GetSupers` | Source similarity and weighted-average agglomeration in `public/rdp-clustering.js` | Ported | Daughter windows, overlap deficit, 0.1 default and WPGMA-like merge regression |
| `MakeNJTreesP` / `CollapseNodesXP3` / `CheckBSTree`; Manual §4.1.4 | `public/rdp-bootstrap-tree.js` plus tree-movement scoring in `public/rdp-clustering.js` | Substantially ported | Six seeded JC/NJ tree paths, bootstrap split support, source-default 50% collapse, all-three orientation movement, cohort caps, and exact-site versus balanced-block provenance are operational. Exact Clearcut consensus serialization, every `CheckBSTree` weight branch, and desktop corpus comparison remain. |
| `MakeSDMP2` / `FillRmat` / `CalCR`; Manual §4.1.4 | All-sequence, three-orientation inference in `public/rdp-clustering.js` | Substantially ported | Every other sequence is tested under every presumed-recombinant orientation; 60-VNP flanks, source breakpoint/tract regions, outside-flank averaging, six-cell category vectors, five relabeling/inversion checks, df=4 correlation test, `MakeProperRCorr` SDM plausibility filtering, detectable-signal overlap, retained unsignalled descendants, seeded bootstrapped-tree evidence, and tunable 2-of-3 rule are operational. |
| `MakePhPrScore`, `MakeTrpGroups`/`MakeTrpScore`, `MakeINList`/`MakeOUCheck`, `SimpleDist`, `MakeSSDistB`, `GetBadDists`, `FindSets`, `MakeEList`/`MakeListCorr`, `MakeLDist`/`MakeRCompat`, `CalcMaxD`/`CMaxD2P3`, `MakeConsensusC`; Manual §4.1.5 | `public/rdp-recombinant-identification.js`, packed dMax kernel in `assembly/index.ts`, persisted event ledger, role workbench and auto-resolve gate | Default decision tree substantially ported | Direct source correlation/sentinel and PS-gate behavior; VB6 banker’s quantization; NO/PI/NI mapping; tree-category movement; historical-event set closure; the four-stage uncollapsed/collapsed/ordinary/set parsimony cascade; O:E/O:EDist; SSDist/OUIndex; distinct-category Conflict; SetDistT/P; packed quartet dMax; all six joint rules; per-polarity final-trim penalties; and the 60% ambiguity rule are operational. The optional logistic/neural selectors and desktop golden-corpus calibration remain. |
| Manual §4.1.6 erase/extract/rescan | `public/rdp-disassembly.js` and component-aware worker scans | Substantially ported | Accepted co-recombinant tracts are erased, gap-padded tract components are appended, nested lineages recursively target prior components, clustering uses the working component alignment, crossing signals split, and gap-near breakpoints are marked uncertain in VNP space. Exact per-method uncertainty windows and desktop iteration-order fixtures remain. |
| `GetFragsP3` / `GetMaxFragScoreP` / `CalcKMaxP` / `GCCalcPValP` | GENECONV fragment locator and probability in `assembly/index.ts` / `public/rdp-statistics.js` | Substantially ported | Independent G=0 exact-run and finite-G integer mismatch scoring, interval recovery, tunable G-scale and KA root/K probability regressions pass. Indel-run modes, overlapping-fragment deletion, pair-scan/global-polymorphism mode and permutations remain. |
| `FastRecCheckMC` / `MakeWindowSizeP` / `ChiPVal2P` | MAXCHI discovery/statistic path | Substantially ported | Independent compressed triplet-polymorphic-site χ² peak pairs and informative-half-window × three-pair probability multiplier pass monomorphic-padding and breakpoint fixtures. Exact `GrowMChiWin`, smoothing/ban windows, lookup-table rounding and multi-peak queue remain. |
| `FastRecCheckChim` | CHIMAERA discovery/statistic path | Substantially ported | Independent recombinant/parent-match binary compression, variable-site χ² peak pairs and source correction are operational. Exact peak growth, missing-data ban windows, lookup-table rounding and multi-peak queue remain. |
| BootScan/RecScan native/VB orchestration | Independent distance-window locator + seeded triplet p-distance bootstrap | Partial | Directional parent-switch discovery, deterministic regional resampling and topology support are operational; multi-taxon NJ/UPGMA, substitution models, per-window source orchestration and source binomial/χ² modes remain. |
| `SetUpSiScan`; `GetSSOL`; `Get3Score`/`GetPScores2`; `DoPerms3P`; `MakeZValue2`; `DoSums`; `FindMaxZ`; `ShrinkRegionC`; `NormalZ` | `public/rdp-siscan.js`, source-confirmation path and multi-run event expansion in `public/rdp-worker.js` | Substantially ported | Exact 15-pattern table; recommended triplet-variable, quartet-variable and all-position modes; gap strip/fifth-state handling; tree-path then direct-distance nearest outlier, most-divergent, analyst-selected and horizontal-randomized fourth-sequence modes; MSVC 15-bit vertical random stream including discarded cells; population-variance Z scores; source topology sums; all disjoint topology runs; boundary shrinkage; and whole-region `NormalZ × L/region` probabilities are operational. Large tables use an equivalent streamed RNG and cached permutation prefixes. Desktop missing/tie fixtures and GUI/manual plot-by-plot comparison remain. |
| RDP5 3Seq orchestration | `method_stats` HGRW locator + `threeSeqExactP` | Partial | Independent maximum-descent interval and exact bounded first-passage DP; desktop large-table/correction comparison remains. |
| `PHITest2`; `PHI`; `pair_score`; `GetFandG`; `AnalyticMeanVariance` | `public/rdp-phi.js` and alignment challenge studio | Substantially ported | Parsimony-informative-site selection, multistate bipartite reticulation score, source PHI neighbourhood mean, analytic moments, VB banker-rounded window conversion and normal lower tail have direct formula regressions. Browser work ceilings use a disclosed position-balanced subset; exact VB random thinning above 6,000 sites and desktop golden fixtures remain. |

The worker no longer uses a common non-RDP discovery seed. Full exploratory
mode enumerates each unordered concrete sequence triplet, and target-oriented
passes evaluate all three possible recombinant assignments; no detector uses
an alignment consensus or the unselected sequences as a third-member proxy.
Directional method
families are scanned in both parent orientations, every retained interval
records its originating family, and evidence is forced to `P=1` when that
family has no co-located signal. Multiplicity counts unique unordered triplets
in exploratory mode; the order-invariant RDP result set is cached so the three
target passes do not rerun the native source scan. RDP and SiScan now both
retain several disjoint source runs per ordered triplet; the other non-RDP
families still expose only their strongest raw run.

## Compatibility policy

- Source mode is deterministic for a saved random seed and records the mode in
  every project.
- Historical quirks are retained when they affect output and are called out in
  code comments/tests rather than silently “fixed.”
- Source-derived and fallback calculations are labeled separately in the UI,
  project JSON, and CSV.
- A method is not called complete until its independent locator, statistic,
  calibration, circular/missing-data behavior, and desktop golden fixtures all
  pass.
