# Feature audit — 0.5.1 checkpoint

This ledger prevents interface completeness from being confused with validated
RDP5 replacement status. “Resolved” means the product defect is implemented and
covered by a local test or build check. It does not imply numerical parity with
the proprietary Windows package.

| Audit item | 0.5.1 state | Evidence / remaining gate |
| --- | --- | --- |
| Desktop page or a tall result panel cannot scroll | **Resolved** | Viewport ownership uses `100dvh`, a zero-basis center scroller, and a bounded scroll body in every reusable analysis panel. Sidebar, inspector, matrices, trees, modals, and long sublists all have explicit scroll owners; short/mobile layouts return scrolling to `body`. Rendered CSS regressions cover every primary owner. |
| Trees are hidden or too weak | **Resolved for exploratory NJ review** | Dedicated Local Trees tab; selected event, parent roles, left/tract/right regions, nearest-context cohort, topology summaries, highlighted leaves, cohort list, and individual/batch Newick export. ML, RF, SH and AU tests remain open. |
| Method results are impossible to distinguish | **Resolved** | Named support chips in the event table plus seven method tabs showing decision, raw p, adjusted p, experiment scope, statistic, calibration, limitation, primary paper, and a method-specific selected-triplet profile. |
| Only one small example | **Resolved** | Eight deterministic, truth-annotated synthetic datasets: clean triplet, tutorial control, three-lineage virus family, circular ssDNA family, four-segment family, nested/overlapping mosaics, 80 kb bacterial core genome, and 520-genome stress panel. Empirical sequences remain intentionally unbundled pending dataset-specific licensing/provenance. |
| All methods score one shared candidate screen | **Open — critical parity gate** | The UI now says this explicitly. Candidate retention preserves nested size classes, scales adaptively, and keeps alternatives, but independent RDP/GENECONV/BootScan/MaxChi/Chimaera/SiScan/3SEQ discovery passes remain required. |
| Exact RDP5 numerical parity | **Open — critical parity gate** | Exact small-case probability tests and simulation fixtures exist; GENECONV mismatch/permutation, multi-taxon BootScan, complete MaxChi/Chimaera, SiScan permutations, BURT, and large-case 3SEQ parity remain. See `SCIENTIFIC_BASIS.md`. |
| Edited-event multiple testing silently becomes one test | **Resolved** | Recalculation retains the original scan triplet count. Holm edits use a labeled conservative full-family Bonferroni correction until a full rescan can re-rank Holm hypotheses. Regression test uses 56 tests. |
| Breakpoint confidence remains stale after edits | **Resolved** | HMM informative-site bounds are retained for linear scan events. Any assignment/breakpoint edit collapses old bounds to the analyst’s exact coordinates and marks evidence stale; recalculation refreshes the state. |
| Disabled methods still consume full compute | **Resolved at statistic-kernel level** | A WebAssembly bitmask skips disabled GENECONV, BootScan, MaxChi, Chimaera, SiScan, 3SEQ, and polishing loops. BootScan resampling and exact 3SEQ DP are fully bypassed when disabled. Candidate discovery remains shared. |
| Fixed top-500 candidate loss | **Resolved as bounded adaptive retention** | Retention scales to 5,000, caps each recombinant at 12, preserves distinct tract-size classes, and still bounds browser work. Independent method candidate queues remain part of the critical parity gate. |
| Alternative parents discarded | **Resolved** | Alternative sampled parent indexes survive deduplication, project round-trip, inspector display, and one-click reassignment. Explicit unsampled-parent models remain open. |
| Unsafe export silently includes unreviewed/stale events | **Resolved** | Sequence-changing export defaults to accepted + fresh. “All fresh” and an explicitly unsafe “all retained” override are visible choices. Full project/CSV always preserve all hypotheses. |
| Undo/delete destroy provenance | **Resolved** | Per-event histories are supplemented by an append-only project ledger. Undo/redo append entries; deletion saves a serialized event tombstone. `.rdpweb` schema is 0.5. |
| Initial 7/7 demo evidence is hand-authored | **Resolved** | Initial event is explicitly known synthetic truth with no method p-values. It exists to expose review/tree UI and must be scanned or recalculated. |
| No iterative analysis / grouped ancestry | **Partially resolved** | Accepted fresh recombinants can be excluded as targets in an unresolved-sequence rescan; overlapping/nested retained signals receive automatic ancestry groups and joint-review warnings. Full disassembly/edit propagation and recombinant-parent ancestry remain open. |
| No formal false-positive statistic | **Partially resolved** | Added a seeded 199-replicate four-gamete proximity permutation p-value alongside bounded incompatibility, rate, ambiguity, gap, and parent-conflict diagnostics. Exact PHI and tree-conditioned homoplasy remain open. |
| Query/reference mode cannot group references or test them | **Resolved** | Editable/inferred reference groups diversify the parent shortlist; Advanced controls can test reference records as recombinants. RDP5 naming-rule parity fixtures remain. |
| No method-specific manual views | **Resolved for bounded selected-triplet profiles** | Each primary method has a distinct tract-aware profile. Independent scan canvases and interactive candidate generation per method remain tied to the critical parity gate. |
| No autosave or stopped-scan recovery | **Resolved** | IndexedDB autosave with restore/start-fresh prompt; worker checkpoints up to 100 partial candidates; Stop recovers them as stale hypotheses requiring recalculation. |
| No command-line batch path | **Resolved for same-engine local batch** | `npm run cli -- ...` runs the checked-in worker/Wasm and writes `.rdpweb` 0.5. Standalone packaging and analysis manifests remain. |
| Limited import/project compatibility | **Open** | FASTA, CLUSTAL, PHYLIP, NEXUS, GFF3, GenBank FEATURES, BED and `.rdpweb` are supported. Native RDP4/RDP5 project conversion and additional downstream formats remain. |
| Validation too narrow | **Partially resolved** | Scientific, circular, large-path, diagnostic, exact-probability, method-mask, example-generation, project, scrolling, static-render, Sites artifact, GitHub Pages subpath, and performance gates pass. RDP5 comparison corpus, sensitivity/specificity curves, and Chrome/Firefox/Safari/Edge E2E remain release blockers. |

## Current release judgment

Version 0.5.1 is a substantially clearer and more recoverable **scientific
alpha**, not a validated drop-in RDP5 replacement. The next highest-value work
is independent candidate discovery plus numerical comparison fixtures for each
primary method; no amount of additional UI should supersede those gates.
