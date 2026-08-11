# RDP Web scientific basis and validation ledger

RDP Web is an MIT-licensed browser implementation of the RDP recombination
workflow. Published methods and permitted source-compatible RDP5 ports coexist;
this ledger distinguishes what is operational from what still needs desktop
corpus validation or implementation before the project can claim full parity.

## Ground rules

- Published papers, the RDP5 manual, and author-supplied RDP5 desktop source are
  authoritative references. Selected active routines are ported with the
  authors' permission; `SOURCE_PARITY.md` records routine-level provenance.
- Detector implementations are derived only from the author-supplied RDP5
  VB/native source and documented source behavior.
- No 3Seq or GENECONV executable is bundled; both have distribution terms that
  are unsuitable for an MIT-only project.
- A matching name in the interface means that an independent statistic informed
  by that published method is available. It does **not** yet mean numerical
  equivalence to the original implementation.
- Until the validation suite is complete, every inferred event is a hypothesis
  that must be independently checked before publication-grade use.

## Current implementation ledger

| RDP5 surface | State in RDP Web | Validation needed |
| --- | --- | --- |
| Aligned FASTA, CLUSTAL, PHYLIP, NEXUS import | Operational, local-only | Corpus of difficult interleaved/quoted files |
| Fully exploratory triplet scan | Role-agnostic default; the scheduler constructs exactly one `a < b < c` record per unordered combination, and the source kernels resolve their internal pair tracks/polarities. Packed extraction retains only that triplet's informative columns. Kernel-call provenance and C(5,3)=10 regressions guard the invariant. | RDP5 desktop fixture comparison and simulation ROC curves |
| Query vs reference scan | Operational roles, editable reference groups, every allowed concrete query × reference pair by default, and opt-in reference-as-recombinant testing | RDP5 naming/group parity fixtures |
| RDP evidence | Source-compatible `FindSubSeqPB3` → `XOHomologyP2` → `FindNextP` → `DefineEventP2` VNP-window locator, RDP5 common/different tract rule, polarity assignment, multiple distinct excursions per triplet with a tunable retention ceiling, and `ProbCalcP/P2`-equivalent binomial-tail scaling | Golden desktop corpus across ties, missing data, circular origins, raw-signal multiplicity, and correction modes |
| GENECONV evidence | Disabled in production; the simplified locator was retired. | Complete `AlistGC2` → `GCXoverDP2` → fragment/probability batch port, then desktop corpus |
| BootScan/RecScan evidence | Disabled in production; the triplet p-distance stand-in was retired. | Complete source multi-taxon window/tree/bootstrap orchestration and desktop corpus |
| MaxChi evidence | One source triplet pass builds all three pair-equality tracks; all-different triplet sites are retained. Source half-window selection, missing/end bans, adjacent χ² profiles, circular 11-position smoothing, peak-basin destruction, GrowMChiWin expansion and bounded multi-peak pairing are active. | Lookup-table rounding and broad missing/circular/tie desktop corpus |
| Chimaera evidence | The same source triplet pass builds all three recombinant-oriented binary tracks, excluding all-different sites. Source window, smoothing, basin, growth and multi-peak behavior is active. | Lookup-table rounding and broad missing/circular/tie desktop corpus |
| SiScan evidence | Source `GetSSOL` → 15-category/sum → `DoPerms3P`/`MakeZValue2` → `FindMaxZ` → `ShrinkRegionC` confirmation runs on source-detector candidates; the generic oriented prelocator is retired. Fourth-sequence, position/gap, permutation, multi-run and whole-region controls are operational. | Independent source discovery orchestration and desktop missing/tie/GUI corpus |
| 3Seq evidence | Disabled in production; the bounded approximation was retired. | Complete author-source discovery/probability-table path and desktop corrected-p corpus |
| Breakpoint polishing | Source-compatible fixed-three-state `BenHMM`/`DoHMMCyclesSerial` fitting with 21 random starts, MSVC RNG, source backtrace/sentinel behavior, 95/99% posterior intervals; plus manual-spec 2–20-state BIC/AIC step-up mode | Golden desktop corpus for tie/backtrace edge cases and circular events |
| Event clustering | Source `GetSupers` daughter-distance/overlap score with WPGMA-like average linkage, plus §4.1.4 characterization of every sequence under all three presumed-recombinant orientations. Source-style 60-VNP regions, outside-flank averaging, six-cell category vectors, `CalCR` relabeling/inversion checks, `MakeProperRCorr` SDM filtering, detectable-signal overlap, and the configurable 2-of-3 rule retain co-recombinant descendants even without an independent raw signal. Six seeded JC/NJ bootstrap trees are constructed per event, low-support branches are collapsed at the source 50% default, and tree-movement evidence is recorded with cohort/block provenance. | Exact `Clearcut` consensus/tree-distance equivalence, full `CheckBSTree` weighting, and desktop event-corpus parity; long regions use tunable balanced site blocks unless the block count covers every site |
| Circular genomes | Opt-in dual-origin triplet scan, circular-breakpoint deduplication, native wrapping coordinates, plots, editing, and masking | Broader circular simulation corpus and RDP5 numerical comparison |
| False-positive flags | Direct RDP5/PHIPACK `PHITest2` multistate incompatibility graph, PHI statistic, analytic mean/variance and normal tail; bounded four-gamete proximity permutation; gap/boundary, local rate-density, and parent-conflict diagnostics exposed as review evidence | Desktop PHI golden-corpus calibration, exact VB random thinning above the work ceiling, tree-conditioned homoplasy and broader challenge calibration |
| Alignment verification | Base view plus arbitrary 2–6-parent affinity highlighting; unique/shared/novel/missing states; parent-informative-site filtering; event-parent defaults; bounded memory guard | Very-large alignment virtualization, codon/amino-acid overlays, and user studies |
| Event verification and editing | Source-guided six-stage Review studio; draggable/numeric breakpoints; seven selected-triplet method profiles; best/all evidence modes; manual create/duplicate/delete; persisted three-polarity default-RDP5 ledger with 18 standalone source statistics, the final-trim penalty, and six joint rules, including historical-event set closure, source quantization, and packed-WASM quartet dMax; individual/group decisions; shared breakpoints; ancestry groups; ordered queue; role-aware tunable auto-resolution; dependency holds and targeted/adaptive rescans; recursive erase/extract reconstruction; structural-gap uncertainty; stale propagation; global mosaic; dependency links; exact-hypothesis recalculation; undo/redo; partial recovery; IndexedDB autosave; audit trails | Optional logistic/neural role selectors, desktop-corpus calibration and auto-decision calibration, exact per-method deleted-boundary thresholds, historical direction inference, unknown-parent/reassortment/ancestral propagation parity, and user studies |
| Pairwise matrix and local topology contrast | Dense sequence p-distance matrix for up to 64 already-calculated sequences, closest-pair contrast, dedicated tract/combined-background or left/tract/right native NJ comparison, explicit unrooted interpretation, adaptive nearest context, linked leaf marking, continuous branch geometry with explicit zero-length nodes, cohort display, and Newick export | ML/Bayesian trees, clade operations, and exact SH/AU/RF matrices |
| Genome-position pattern matrices | Symmetric breakpoint-pair counts with RDP4 Figure 2c semantics; split region matrix with observed tract-separation counts and analytical circular random-placement residuals; split sampled local p-distance-profile RMS/correlation-loss matrix | RDP4/RDP5 permutation calibration, SH/RF numerical parity, confidence overlays, and feature-association matrices |
| Breakpoint density | Descriptive density and seeded uniform-null hotspot permutation, integrated with the dense breakpoint-pair matrix | RDP5 hot/cold-spot model parity and association covariates |
| Recombination-free outputs | Remove, mask, CDS-phase-aware mask, split and partition, including wrapping tracts | More frame/feature edge cases and every RDP5 export variant |
| Project/results output | Restorable `.rdpweb` schema 0.5, legacy schema import, immutable project audit ledger, annotations, edit history, diagnostics, circular/HMM provenance, CSV, and a same-engine Node batch runner | Long-term migration corpus, analysis manifests and RDP5 project conversion |
| Annotation | GFF3, GenBank FEATURES and BED import; mapped feature track and GFF3 export | ORF calling, multi-record mapping and richer feature analyses |
| PDB/SCHEMA, LDHat, ancestral inference | Not implemented | Separate workstreams |

## Ordered event reconstruction semantics

The 0.6 workflow follows the public RDP5 manual’s distinction between raw
signals and unique events. It treats the array order as characterization/review
order, not as literal historical time. Earlier edits to recombinant identity,
parent proxies, breakpoints, circular status, or event grouping mark the edited
hypothesis and later retained characterizations stale. The safe export path then
excludes them until recalculation or an unresolved-sequence rescan.

The reconstruction model exposes three relationships as hypotheses:

1. overlapping tracts in the same recombinant are labeled as possible
   overprinting, never automatically asserted as historical direction;
2. a sequence that is itself recombinant and is used as a parent elsewhere is
   labeled as a recombinant-parent dependency; and
3. a shared analyst group labels co-recombinant descendants.

“Major parent” and “minor parent” remain sampled proxy relatives. They must not
be interpreted as proof that the exact historical donor was sampled. Linked
marking in the regional trees and the parent-affinity alignment are verification
aids for these assignments, not independent statistical tests.

The 0.7 auto-resolver is an explicitly labeled decision-support heuristic, not
a reverse-engineered RDP5 classifier. Automatic acceptance requires fresh,
calibrated evidence to pass configurable method-count, adjusted-P,
informative-site, breakpoint-uncertainty, parent-conflict, rate-density,
diffuse-incompatibility, and warning gates. A weighted score only acts inside
those gates. Ambiguous or stale hypotheses remain for the analyst; when such a
hypothesis has a strong dependency footprint, its causally linked downstream
branch is also held while independent branches may continue.

After a changed decision, the runner estimates downstream impact only from
explicit relationships: tract overlap in the same recombinant, use of a
resolved recombinant as a later parent proxy, reciprocal nested-parent use,
and shared co-recombinant groups. Crossing the configured risk threshold
creates a barrier. The remaining queue is not processed until an impacted-only,
adaptive, or full unresolved-target rescan completes and the plan is rebuilt.
During that rescan, every accepted tract is erased from the selected
co-recombinant remainder rows and copied into a gap-padded internal component,
including recursive lineage targeting for nested events. Extracted components
remain eligible as targets and parent proxies while accepted mosaic remainders
are withheld. New signals crossing structural gaps are divided into continuous
pieces and breakpoints within the RDP VNP-window distance of deleted sequence
are marked uncertain. Internal component indexes are never exposed as user
sequence identities; their lineage is retained in the event dossier and
project. This implements the manual's §4.1.6 state transition, while exact
method-specific uncertainty thresholds and historical direction still require
desktop-corpus validation.

The Review studio follows the manual and desktop source's user-facing
transaction loop rather than imitating the old form layout. Queue navigation
can skip accepted/rejected items and filter by method count; selected evidence,
breakpoint editing, role challenges, tract/background trees, grouping, and
decisions are kept in one continuous context. Rescanned events use the direct
desktop-default `MakeConsensusC` ensemble; imported/manual events without a
saved source ledger retain the clearly labeled bounded identity-switch
diagnostic until recalculated. Group-wide decisions and shared breakpoint edits are atomic
undo/audit actions and invalidate affected downstream evidence.

## Engine design

The WebAssembly engine uses three layers:

1. Small/medium alignments use an exact packed canonical-site distance matrix.
   Sixteen bases are compared per 32-bit word, with a separate validity mask.
2. Large alignments avoid materializing the O(N²L) display matrix. A sampled
   distance pass supports the UI and optional explicitly approximate preview,
   while the default scientific scan still enumerates every concrete triplet.
   The visible matrix is calculated exactly for the first 24 sequences. The
   clustering normalizer uses a separate streaming packed kernel to obtain the
   exact global maximum distance without materializing an N×N matrix.
3. The enabled RDP method has an independent O(L) source-compatible VNP scan
   that retains multiple distinct excursions per triplet up to a visible,
   tunable ceiling and reports any truncation. Its order-invariant result set
   is cached per unordered triplet to avoid repeating the native scan during
   the three target passes.
   Every other enabled family now contributes its own full-alignment interval:
   GENECONV fragments, topology windows, compressed-variable-site MAXCHI,
   binary CHIMAERA, a fast oriented SISCAN locator, and maximum-HGRW descent.
   Directional families are run in both parent orientations. A method is
   excluded from an event's evidence family unless its interval is co-located
   with that event. A WebAssembly bitmask skips disabled kernels; seeded column
   bootstraps run natively only when BootScan is enabled. A preliminary SiScan
   interval is then confirmed by the supplied-source fourth-sequence,
   15-category, vertical-permutation and topology-run workflow; its several
   significant runs are queued independently. Fixed-outgroup windows use
   rolling category counts and cached exact permutation-prefix ranges, while
   very large historical random tables are regenerated from the identical
   MSVC stream rather than materialized.
4. Exact 3SEQ first-passage probabilities are computed in the worker within a
   per-event and per-job work budget; larger cases retain a labeled
   conservative bound so one alignment cannot monopolize the browser.
5. Breakpoints can be refined with source-compatible fixed-three-state BURT or
   the manual's optional 2–20-state step-up model over parent-discriminating
   sites. Source mode reproduces circular informative-site half-copy padding,
   its sentinel/crop, `.995`/`.999` intervals, VNP CI matching and the principal
   non-reassortment polishing rules. Circular detection also scans a second
   coordinate origin, which is not double-counted in the multiplicity family.

Exploratory multiplicity counts unique unordered triplets rather than the
three target-oriented passes. Query/reference mode counts each query-specific
reference pair. This matches the hypothesis family more closely and prevents
the implementation detail of polarity evaluation or circular rotation from
inflating adjusted p-values.

The RDP5-parity default pays the true O(N³L) triplet cost. An explicit
"approximate parent shortlist" preview can reduce this toward O(NK²L), where K
defaults to eight, but its UI and saved provenance state that unlisted triples
were not tested. That preview never constructs a consensus/rest-of-alignment
proxy: every executed method call still receives three real sequence indexes.
The worker boundary keeps rendering and cancellation responsive; typed arrays
avoid per-site JavaScript objects.

GitHub Pages cannot guarantee cross-origin-isolation headers, so shared-memory
threads are not a baseline assumption. A later pool can use transferable
partitions, but only after browser memory/merge-order benchmarks show a net win.

## Primary RDP lineage

1. Martin D, Rybicki E. 2000. *RDP: detection of recombination amongst aligned
   sequences.* Bioinformatics 16:562–563.
   https://doi.org/10.1093/bioinformatics/16.6.562
2. Martin DP, Williamson C, Posada D. 2005. *RDP2: recombination detection and
   analysis from sequence alignments.* Bioinformatics 21:260–262.
   https://pubmed.ncbi.nlm.nih.gov/15377507/
3. Martin DP et al. 2010. *RDP3: a flexible and fast computer program for
   analyzing recombination.* Bioinformatics 26:2462–2463.
   https://pubmed.ncbi.nlm.nih.gov/20798170/
4. Martin DP et al. 2015. *RDP4: Detection and analysis of recombination
   patterns in virus genomes.* Virus Evolution 1:vev003.
   https://doi.org/10.1093/ve/vev003
5. Martin DP et al. 2021. *RDP5: a computer program for analyzing recombination
   in, and removing signals of recombination from, nucleotide sequence
   datasets.* Virus Evolution 7:veaa087.
   https://doi.org/10.1093/ve/veaa087

## Practical and tutorial sources

- Martin DP, Lemey P, Posada D. 2011. *Analysing recombination in nucleotide
  sequences.* Molecular Ecology Resources 11:943–955.
  https://doi.org/10.1111/j.1755-0998.2011.03026.x
- Martin DP, Murrell B, Khoosal A, Muhire B. 2017. *Detecting and Analyzing
  Genetic Recombination Using RDP4.* Methods in Molecular Biology 1525:433–460.
  https://doi.org/10.1007/978-1-4939-6622-6_17
- Martin DP. *RDP5 Instruction Manual.* Detailed settings, displays, method
  descriptions, and step-by-step workflow.
  https://web.cbio.uct.ac.za/~darren/RDP5Manual.pdf
- Sattar MN et al. 2025. *Recombination Analysis of Geminiviruses Using RDP.*
  Methods in Molecular Biology. Practical dataset preparation, event
  characterization, and recombination-free output workflow.
  https://pubmed.ncbi.nlm.nih.gov/40064777/
- Lemey P, Posada D. 2009. *Introduction to recombination detection.* In *The
  Phylogenetic Handbook*, 2nd ed., pp. 493–518.
  https://doi.org/10.1017/CBO9780511819049.017
- Salminen M, Martin DP. 2009. *Detecting and characterizing individual
  recombination events.* In *The Phylogenetic Handbook*, 2nd ed., pp. 519–546.

## Primary method sources

- GENECONV: Padidam M, Sawyer S, Fauquet CM. 1999. Virology 265:218–225.
  https://doi.org/10.1006/viro.1999.0058
- BootScan/RecScan: Martin DP et al. 2005. AIDS Research and Human
  Retroviruses 21:98–102. https://doi.org/10.1089/aid.2005.21.98
- MaxChi: Maynard Smith J. 1992. Journal of Molecular Evolution 34:126–129.
  https://doi.org/10.1007/BF00182389
- Chimaera evaluation: Posada D, Crandall KA. 2001. PNAS 98:13757–13762.
  https://doi.org/10.1073/pnas.241370698
- SiScan: Gibbs MJ, Armstrong JS, Gibbs AJ. 2000. Bioinformatics 16:573–582.
  https://doi.org/10.1093/bioinformatics/16.7.573
- 3Seq: Boni MF, Posada D, Feldman MW. 2007. Genetics 176:1035–1047; Lam HM,
  Ratmann O, Boni MF. 2018. Molecular Biology and Evolution 35:247–251.
  https://doi.org/10.1093/molbev/msx263
- LARD: Holmes EC, Worobey M, Rambaut A. 1999. Molecular Biology and
  Evolution 16:405–409.
- TOPAL/DSS: McGuire G, Wright F. 2000. Bioinformatics 16:130–134.
- PhylPro: Weiller GF. 1998. Molecular Biology and Evolution 15:326–335.
- VisRD: Lemey P et al. 2009. BMC Bioinformatics 10:126.
- PHI: Bruen TC, Philippe H, Bryant D. 2006. Genetics 172:2665–2681.

## Acceptance criteria for a replacement claim

The label “RDP5 replacement” should be used only after all of the following are
true:

1. Every primary method has fixture-level numerical tests and simulation-based
   sensitivity/specificity benchmarks.
2. False-positive filters reproduce the documented interpretation on
   misalignment, rate-variation, and homoplasy challenge sets.
3. Circular genomes, reassortment, recombinant ancestors, grouped events, and
   unknown parents are covered by regression fixtures.
4. The manual review workflow has parity for breakpoint, recombinant, parent,
   grouping, and event-state corrections.
5. Recombination-aware phylogenetics, all matrix/hotspot/association analyses,
   annotation, SCHEMA, and every documented export are implemented or clearly
   scoped into interoperable companion tools.
6. Performance is reported for the RDP4/RDP5 reference workloads and does not
   freeze the interface or exhaust browser memory.
