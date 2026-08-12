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
2,632 unique concrete-triplet comparisons, source RDP, six-track GENECONV and
combined MAXCHI/CHIMAERA detection. The same
command also runs an independent 4 × 80 kb source SiScan fixture with 100
window and 1,000 final-region permutations, then all 220 unordered triplets of
a 12 × 2,400 standalone-SiScan cohort. A representative run reported:

| Stage | Time |
| --- | ---: |
| Scalar distance reference | 251.93 ms |
| Packed production distance | 17.26 ms |
| Packed distance speed-up | 14.59× |
| All source triplet candidate scans | 2,475.26 ms |
| Six-track GENECONV portion | 646.91 ms |
| GENECONV throughput | 40.7 million full triplet-sites/s |
| Packed 30-taxon VisRD/dMax role statistic | 176.96 ms |
| Source PHI, 100 taxa × 256 retained / 9,928 informative sites | 60.46 ms |
| Production-kernel total | 2,729.93 ms |
| Aggregate triplet scan throughput | 10.6 million site-triplets/s |
| VisRD throughput | 556.9 million site-quartets/s |
| Source SiScan, 4 × 80 kb, 100/1,000 permutations | 448.78 ms |
| Standalone SiScan, 12 × 2.4 kb, all 220 triplets | 806.75 ms |
| Standalone SiScan throughput | 272.7 complete triplets/s |
| Hardware-normalized SiScan throughput | 507.2 complete triplets/s |

`npm run bench:gate` retains hard ceilings of 150 ms for packed distance,
500 ms each for the production 30-taxon VisRD/dMax cohort and bounded source
PHI, and 2 s plus exact 30,000–45,000 tract recovery for the 80 kb source
SiScan workload. A bounded runner-normalized 1.8 s base ceiling and a minimum
200 normalized complete triplets/s separately guard the 220-triplet independent
SiScan workload. The aggregate detector and six-track GENECONV ceilings/floor
start at 2 s, 500 ms, and 50 million complete triplet-sites/s, then normalize
by an unchanged scalar-distance calibration up to a tightly capped 2.5×
hardware factor. This preserves the original strict gate on fast runners while
preventing CPU-throttled containers from failing an unchanged binary solely
because their sustained single-core rate is lower.

This measures WebAssembly kernels in Node, not an RDP5-equivalent seven-program
analysis. It now includes the production VisRD/dMax path used by source-parity
recombinant-role consensus, but must not be compared directly with RDP4/RDP5
wall-clock figures.

BootScan has its own shared-pair/tree benchmark:

```sh
node bench/bootscan-batch.mjs 24 2000 distance
node bench/bootscan-batch.mjs 24 2000 upgma
node bench/bootscan-batch.mjs 24 2000 neighbor-joining
```

On the same development class, the 2,024-triplet, 102-window, 100-replicate
workload completed in 712 ms (distance), 733 ms (UPGMA), and 781 ms (NJ).
Distance mode compacts requested pairs in shortlisted scans. Tree modes require
all active-cohort pairs but build each tree once per window/replicate and reuse
its leaf-path matrix across all 2,024 triplets; they never construct a separate
three-taxon tree for each triplet.
Browser, device, rendering, and file-parsing costs are not included. The
synthetic 12 × 2,400 end-to-end module-worker fixture remains an automated
source-SiScan, bootstrap, challenge-diagnostic and exact-probability workload,
and localizes both inserted mosaic tracts. A 600-site circular fixture
recovers a known 520→100 origin-spanning event. A separate 513-sequence fixture
forces the sampled/stratified large-data path and verifies a compact 24 × 24
display matrix without allocating an N² matrix.

The source BootScan batch has its own reproducible all-triplet benchmark:

```sh
npm run bench:bootscan
```

The August 2026 development-container run for 24 sequences × 2,000 sites,
all 2,024 concrete triplets, 276 unique pairs, 102 windows and 100 replicates
took 307.01 ms after warm-up. It evaluated 67.2 million triplet-window-
replicate relationships/s in 0.71 MiB of kernel workspace. Sharing each pair
row across its 22 containing triplets avoids the 6,072 pair rows per window a
triplet-at-a-time implementation would compute.

The fused 3Seq kernel also has a dedicated comparison:

```sh
npm run bench:three-seq
```

For 32 sequences × 5,000 sites (4,960 concrete triplets; 24.8 million
triplet-sites), the 0.9.7 August 2026 development-container run took 896.79 ms
in the contiguous-byte kernel and 1,077.36 ms through the exact packed oracle,
27.7 million triplet-sites/s on the faster path. Unlike the earlier linear-only
measurement, this includes production `CheckwrapC`: all three target walks
retain compressed coordinates and cumulative heights and run the source
bounded-origin pass. The reusable scratch block is 120,020 bytes at 5,000
sites (24 bytes/site plus five result words). Production uses the contiguous
representation already required by the other source routines. It decodes each
triplet once and updates all three possible recombinant walks together; the
packed export remains exact regression coverage and a future SIMD target.

## Scaling choices

- Canonical bases are represented both as scan bytes and as 2-bit lanes; the
  packed exact distance kernel uses the WebAssembly `i32.popcnt` instruction
  and is bit-for-bit tested against the scalar kernel.
- Full N² distances are skipped when N > 512 or the projected word comparisons
  exceed 50 million.
- Large-data parent panels combine nearest sampled references with stratified
  references so dense clone groups do not consume every candidate slot.
- Full parity mode defaults to all C(N,3) concrete triplets. Parent pruning is
  an explicit non-parity preview that changes O(N³L) toward O(NK²L).
- The source-only scheduler materializes only `a < b < c`; RDP, GENECONV,
  SiScan, fused three-role 3Seq and the combined MAXCHI/CHIMAERA kernel are each
  invoked once per unordered triplet. 3Seq decodes the three bases once per
  column, updates every target walk without rereading the alignment, and
  reuses one O(L) `CheckwrapC` workspace across all triplets.
- The BootScan/RecScan scheduler also materializes each unordered triplet once,
  but sends the whole request set to one source-shaped batch. One `SEQBOOT2`
  table is shared globally; each requested pair/window distance row is computed
  once and reused by every containing triplet. Approximate/query batches compact
  unused pairs out of the hot loop and bootstrap matrix, complete high-identity
  windows touch weights only at mismatches, distances use packed 16-bit rows,
  and the valid/difference lookup stores only its triangular half. Exact
  whole-alignment identity fractions retain the source baseline topology even
  when all three window-model JC distances would saturate.
- Two-bit production extraction advances sixteen alignment columns per word,
  builds informative coordinates for the current triplet only, and is required
  to match the byte oracle exactly.
- Cyclical detection does not repeat the full O(N³L) screen after every event.
  It retains unaffected signals, invalidates only hypotheses whose concrete
  roles contain a newly split origin, and re-runs only triplets containing an
  affected origin—the browser equivalent of RDP5's redo list. A pooled signal
  is refreshed against the current component alignment before application, so
  this optimization does not reuse stale role/group characterization.
- The entire computation runs in a dedicated worker, and cancellation
  terminates that worker immediately.
- Sequence bytes are packed into a single typed array; the kernel allocates no
  per-site objects and reuses two prefix buffers for all triplets.
- HMM refinement reuses those same O(L) scratch buffers and stores one compact
  predecessor word per informative site.
- Circular mode allocates one rotated byte copy for detector families that need
  a second origin. Source 3Seq runs only the original view because `CheckwrapC`
  already handles the origin; this halves its circular kernel calls and avoids
  origin-dependent complementary-walk selection. Ordinary linear analyses
  retain the original memory and scan cost.
- Exhaustive mode remains available when K-pruning is scientifically
  inappropriate.
- Seeded bootstrap resampling is compiled into WebAssembly and capped at 1,000
  replicates; the default 100-replicate, three-region pass adds bounded work.
- Source SiScan keeps the desktop random stream but does not allocate the
  historical `(permutations+1) × (alignment length+1)` byte table once it would
  exceed 8 MB. It regenerates the exact MSVC values in the worker, rolls
  fixed-outgroup category counts between windows, caches permutation-class
  prefix ranges shared by every triplet with the same window/replicate setup,
  telescopes the source's adjacent category ranges into three exact band
  endpoints, and reuses deterministic permutation moments for all observed
  category vectors with those same band totals. Final-region prefixes grow
  within a 72 MiB hard accelerator budget; larger tracts fall back to the exact
  streamed enumeration. Identical category-count vectors remain memoized
  exactly within bounded maps.
  Its separate 200/20 source-default window/step avoids inheriting the much
  finer global breakpoint step unless the analyst requests that extra work.
- The source PHI incompatibility graph is O(S²N), so the worker retains up to
  384 informative sites at ≤96 taxa, 256 at 97–256 taxa, and 160 above 256
  taxa. Selection is deterministic and position-balanced; projects record
  retained and total informative-site counts and never label a bounded result
  as an all-site PHI test.
- Production discovery currently enables direct-source RDP, six-track
  GENECONV, distance/UPGMA/NJ BootScan/RecScan, MAXCHI/CHIMAERA, SiScan and fused 3Seq
  paths. Small 3Seq walks use cached exact first-passage results; a cheap
  `SiegmundDiscrete` plausibility screen keeps exact DP work off clearly null
  triplets, and large walks use the same source approximation branch. Source
  SiScan independently screens each emitted unordered triplet once.
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
