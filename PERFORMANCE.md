# Performance notes

Performance work is split into reproducible kernel benchmarks and end-to-end
browser workloads. The former is available now; the latter belongs in the
cross-browser validation suite.

## Current kernel benchmark

Run:

```sh
node bench/benchmark.mjs 100 10000
```

The August 2026 development-container run for 100 sequences × 10,000 sites
(one million aligned nucleotides), eight parent candidates per recombinant,
2,800 O(L) triplet comparisons, method-family statistics on the 500 strongest
raw signals, and windowless HMM refinement on those retained signals. The same
command also runs an independent 4 × 80 kb source SiScan fixture with 100
window and 1,000 final-region permutations. It reported:

| Stage | Time |
| --- | ---: |
| Scalar distance reference | 136.45 ms |
| Packed production distance | 8.93 ms |
| Packed distance speed-up | 15.28× |
| Triplet candidate scans | 189.38 ms |
| Seven-family evidence + 100×3 bootstrap kernels | 733.27 ms |
| Legacy HMM compatibility kernel | 57.10 ms |
| Packed 30-taxon VisRD/dMax role statistic | 81.50 ms |
| Source PHI, 100 taxa × 256 retained / 9,928 informative sites | 41.07 ms |
| Production-kernel total | 1,111.24 ms |
| Triplet scan throughput | 147.9 million site-comparisons/s |
| VisRD throughput | 1.209 billion site-quartets/s |
| Source SiScan, 4 × 80 kb, 100/1,000 permutations | 287.30 ms |

`npm run bench:gate` enforces deliberately hardware-tolerant CI ceilings for
the same workload (150 ms packed distance, 500 ms for the production 30-taxon
VisRD/dMax cohort, 500 ms for bounded source PHI, 2 s production total, and at least 30 million triplet
site-comparisons/s), plus a 2 s ceiling and exact 30,000–45,000 tract-recovery
gate for the 80 kb source SiScan workload, so large regressions cannot ship
silently while ordinary runner variance remains harmless.

This measures WebAssembly kernels in Node, not an RDP5-equivalent seven-program
analysis. It now includes the production VisRD/dMax path used by source-parity
recombinant-role consensus, but must not be compared directly with RDP4/RDP5
wall-clock figures.
Browser, device, rendering, and file-parsing costs are not included. The
synthetic 12 × 2,400 end-to-end module-worker fixture completes in roughly
0.9–1.2 s with source SiScan, bootstraps, challenge diagnostics and exact
bounded probability work, and localizes both inserted mosaic tracts. A 600-site circular fixture
recovers a known 520→100 origin-spanning event. A separate 513-sequence fixture
forces the sampled/stratified large-data path and verifies a compact 24 × 24
display matrix without allocating an N² matrix.

## Scaling choices

- Canonical bases are represented both as scan bytes and as 2-bit lanes; the
  packed exact distance kernel uses the WebAssembly `i32.popcnt` instruction
  and is bit-for-bit tested against the scalar kernel.
- Full N² distances are skipped when N > 512 or the projected word comparisons
  exceed 50 million.
- Large-data parent panels combine nearest sampled references with stratified
  references so dense clone groups do not consume every candidate slot.
- Parent pruning defaults to K=8, changing the dominant triplet cost from
  O(N³L) toward O(NK²L).
- Prefix sums make each retained parent-pair scan O(L).
- The entire computation runs in a dedicated worker, and cancellation
  terminates that worker immediately.
- Sequence bytes are packed into a single typed array; the kernel allocates no
  per-site objects and reuses two prefix buffers for all triplets.
- HMM refinement reuses those same O(L) scratch buffers and stores one compact
  predecessor word per informative site.
- Circular mode alone allocates one rotated byte copy and runs the second
  origin; ordinary linear analyses retain the original memory and scan cost.
- Exhaustive mode remains available when K-pruning is scientifically
  inappropriate.
- Seeded bootstrap resampling is compiled into WebAssembly and capped at 1,000
  replicates; the default 100-replicate, three-region pass adds bounded work.
- Source SiScan keeps the desktop random stream but does not allocate the
  historical `(permutations+1) × (alignment length+1)` byte table once it would
  exceed 8 MB. It regenerates the exact MSVC values in the worker, rolls
  fixed-outgroup category counts between windows, and caches permutation-class
  prefix ranges shared by every triplet with the same window/replicate setup.
- The source PHI incompatibility graph is O(S²N), so the worker retains up to
  384 informative sites at ≤96 taxa, 256 at 97–256 taxa, and 160 above 256
  taxa. Selection is deterministic and position-balanced; projects record
  retained and total informative-site counts and never label a bounded result
  as an all-site PHI test.
- A method bitmask now prevents disabled GENECONV, BootScan, MaxChi, Chimaera,
  SiScan, 3SEQ, and local-polishing loops from running. In particular, disabled
  BootScan performs no resampling and disabled 3SEQ performs no exact DP.
- Exact 3SEQ dynamic programming has a four-million-operation event guard and
  a 20-million-operation job budget, with tuple caching and a labeled fallback.
- Pre-scan dataset summaries, uncomputed matrix fallbacks and local NJ trees
  use explicit stratified work bounds; the overview and alignment viewer cap
  rendered rows while preserving every event-bearing sequence and searchable
  role editing. This avoids main-thread O(N²L) work before the worker starts.
- Candidate retention scales from 500 to a bounded 5,000 with a 12-event
  per-recombinant guard, and discovery checkpoints at most 100 partial
  candidates so stopping a long job can recover useful hypotheses without an
  unbounded message or render cost.
- Parent-affinity indexing is a tight O(LP) string scan for P ≤ 6 and is
  disabled above two million sites; coordinate-window highlighting remains
  available without allocating an all-site index.
- Global reconstruction derives parent dependencies through sequence/group
  indexes instead of comparing every event pair. The queue pages 250 events,
  mosaic rendering caps at 300 recombinant rows, dependency rendering caps at
  500, and derived relationships are bounded at 20,000 while the project keeps
  every event.
- Auto-resolution scores events in O(E) and uses recombinant, parent, group,
  and genomic-interval indexes to visit plausible downstream dependencies
  rather than comparing every pair. It calculates only the next operational
  rescan barrier because the queue will be rebuilt after that rescan. Deferred
  dry-run rendering keeps sliders responsive. A 5,000-event development check
  took about 64 ms for a sparse 520-sequence queue and 25 ms for a deliberately
  concentrated same-recombinant queue.
- The Review studio filters and navigates the queue in O(E). Its three-polarity
  role challenge samples at most 4,096 positions per tract/background segment
  and advances directly by the sampling stride, so opening the inspector does
  not add a full-length pass on long bacterial alignments.
- Dense genome-position matrices render through one canvas per view rather than
  thousands of interactive DOM cells. Breakpoint pairs are O(E); the region
  matrix is O(ER²) with an adaptive cap of R=48 above 2,000 visible events and
  R=64 above 500; local discordance samples at most 24 sequences and 72 sites
  per window. A synthetic 5,000-event/520-sequence development-container check
  completed all three default bounded matrix calculations in roughly 80 ms.

## Next optimization gates

1. Establish browser benchmarks for the RDP5 reference workloads (100 × 10 kb,
   4,000 × 10 kb, long bacterial genomes) with memory ceilings and truth sets.
2. Benchmark transferable worker partitions on GitHub Pages; do not assume
   shared memory or cross-origin isolation.
3. Add SIMD kernels for triplet state extraction after
   scalar/SIMD equivalence tests exist.
4. Cache probability and bootstrap lookup tables by parameter tuple, following
   the performance strategy described for RDP5 while preserving numerical
   traceability.
5. Profile parsing and visualization separately so UI density never obscures
   engine throughput.
