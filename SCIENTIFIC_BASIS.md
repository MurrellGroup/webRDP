# RDP Web scientific basis and validation ledger

RDP Web is a clean-room, MIT-licensed browser implementation of an RDP-style
recombination workflow. This document distinguishes what is operational from
what still requires numerical validation or implementation before the project
can claim parity with RDP5.

## Ground rules

- Published method papers and the public RDP5 manual are behavioral sources.
- No OpenRDP (GPL-3.0) or Windows RDP source code is copied or translated.
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
| Fully exploratory triplet scan | Operational | RDP5/OpenRDP fixture comparison and simulation ROC curves |
| Query vs reference scan | Operational roles, editable reference groups, group-diversified parent shortlist, and opt-in reference-as-recombinant testing | RDP5 naming/group parity fixtures |
| RDP evidence | Exact binomial tail on informative-site identity shifts with within-triplet window correction | RDP5 window/reference-mode corrected-p parity |
| GENECONV evidence | G-scale 0 concordant-fragment scan with a conservative run bound | Mismatch scoring, KA parameters and global permutation parity |
| BootScan/RecScan evidence | Seeded native p-distance triplet bootstrap at tract/flank windows plus an all-window topology sign statistic | Full multi-taxon NJ/substitution-model and cutoff parity |
| MaxChi evidence | Pairwise variable/non-variable χ² at shared candidate boundaries | Independent candidate scan, polymorphic-window and corrected-p parity |
| Chimaera evidence | Compressed binary-triplet χ² at shared candidate boundaries | Independent candidate scan, RDP5 window and event-region parity |
| SiScan evidence | Oriented site-category Z surrogate | Fourth-sequence selection and vertical/horizontal permutations |
| 3Seq evidence | Maximum HGRW descent with exact first-passage dynamic programming inside a bounded work budget and a conservative large-case fallback | Published large-table/algorithm scale comparisons and RDP5 corrected-p parity |
| Breakpoint polishing | Candidate-seeded, windowless two-state HMM/Viterbi path with local χ² fallback | BURT 2–20-state step-up fitting, random-start training, EM, and confidence-interval parity |
| Circular genomes | Opt-in dual-origin triplet scan, circular-breakpoint deduplication, native wrapping coordinates, plots, editing, and masking | Broader circular simulation corpus and RDP5 numerical comparison |
| False-positive flags | Gap/boundary checks, bounded four-gamete incompatibility, seeded four-gamete proximity permutation, PHI-style contrast, local rate-density shift, and parent-conflict diagnostics exposed as review evidence | Exact PHI p-values, tree-conditioned homoplasy and broader challenge calibration |
| Event verification and editing | Draggable/numeric breakpoints, seven selected-triplet method profiles, manual create/duplicate/delete, recombinant/parent reassignment, automatic ancestry groups, unresolved-target rescans, scan-scope recalculation, undo/redo, partial recovery, IndexedDB autosave, and event/project audit trails with deletion tombstones | User studies and explicit unknown-parent/reassortment/ancestral propagation parity |
| Pairwise matrix and local topology contrast | Pairwise matrix, closest-pair contrast, dedicated left/tract/right native NJ comparison, adaptive nearest context, role highlighting, cohort display, and Newick export | ML trees and SH/AU/RF matrices |
| Breakpoint density | Descriptive density, circular-aware pair matrix and seeded uniform-null hotspot permutation | RDP5 hot/cold-spot model parity and association covariates |
| Recombination-free outputs | Remove, mask, CDS-phase-aware mask, split and partition, including wrapping tracts | More frame/feature edge cases and every RDP5 export variant |
| Project/results output | Restorable `.rdpweb` schema 0.5, legacy schema import, immutable project audit ledger, annotations, edit history, diagnostics, circular/HMM provenance, CSV, and a same-engine Node batch runner | Long-term migration corpus, analysis manifests and RDP5 project conversion |
| Annotation | GFF3, GenBank FEATURES and BED import; mapped feature track and GFF3 export | ORF calling, multi-record mapping and richer feature analyses |
| PDB/SCHEMA, LDHat, ancestral inference | Not implemented | Separate workstreams |

## Engine design

The WebAssembly engine uses three layers:

1. Small/medium alignments use an exact packed canonical-site distance matrix.
   Sixteen bases are compared per 32-bit word, with a separate validity mask.
2. Large alignments avoid the O(N²L) matrix. A 64/128/256-site deterministic
   screen plus a stratified panel identifies plausible parents, while the
   visible matrix is calculated exactly for the first 24 sequences.
3. Each retained recombinant/parent triplet is scanned in O(L) with prefix sums
   and a detrended cumulative-sum excursion. The candidate tract is evaluated
   by method-specific GENECONV-run, topology-window, MAXCHI, CHIMAERA,
   SISCAN-category, and maximum-HGRW-descent kernels. Candidate discovery is
   still shared and is not yet seven independent scans. A WebAssembly bitmask
   skips disabled kernels; seeded column bootstraps run natively only when
   BootScan is enabled.
4. Exact 3SEQ first-passage probabilities are computed in the worker within a
   per-event and per-job work budget; larger cases retain a labeled
   conservative bound so one alignment cannot monopolize the browser.
5. Breakpoints can be refined with an O(L) two-state Viterbi pass over only
   parent-discriminating sites. Circular mode adds one half-genome rotation,
   maps candidates back to native coordinates, and corrects probabilities for
   both tested origins; the linear fast path pays none of that scan cost.

Default parent pruning changes the triplet cost from O(N³L) toward O(NK²L),
where K defaults to eight. For large N the parent screen is O(N²S), where S is
at most 256 and falls to 64 above 2,000 sequences. Exhaustive mode is available
for definitive analyses. The worker boundary keeps rendering and cancellation
responsive; typed arrays avoid per-site JavaScript objects.

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
