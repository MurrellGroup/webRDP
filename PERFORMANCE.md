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
raw signals, and windowless HMM refinement on those retained signals reported:

| Stage | Time |
| --- | ---: |
| Scalar distance reference | 141.72 ms |
| Packed production distance | 8.80 ms |
| Packed distance speed-up | 16.11× |
| Triplet candidate scans | 191.47 ms |
| Seven-family evidence + 100×3 bootstrap kernels | 180.36 ms |
| Two-state HMM polishing | 40.86 ms |
| Production-kernel total | 421.49 ms |
| Triplet scan throughput | 146.2 million site-comparisons/s |

`npm run bench:gate` enforces deliberately hardware-tolerant CI ceilings for
the same workload (150 ms packed distance, 1.5 s production total, and at
least 30 million triplet site-comparisons/s) so large regressions cannot ship
silently while ordinary runner variance remains harmless.

This measures WebAssembly kernels in Node, not an RDP5-equivalent seven-program
analysis, and must not be compared directly with RDP4/RDP5 wall-clock figures.
Browser, device, rendering, and file-parsing costs are not included. The
synthetic 12 × 2,400 end-to-end module-worker fixture completes in roughly
40–70 ms with bootstraps, challenge diagnostics and exact bounded probability
work, and localizes both inserted mosaic tracts. A 600-site circular fixture
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
