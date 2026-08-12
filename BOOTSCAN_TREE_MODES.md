# BootScan/RecScan relationship-mode implementation

This note records the 0.10.0 implementation boundary for the three automatic
relationship choices described in the RDP manual. It is intentionally specific
about the unit of work because a tree constructed separately for each tested
triplet is not equivalent to RDP's scan.

## Shared batch order

For each ordinary or circular-coordinate view, the worker:

1. Enumerates each requested unordered concrete triplet once.
2. Builds one deterministic `SEQBOOT2` weight table for the requested window
   size, replicate count, and saved random seed.
3. Calculates every required pair once at each window and replicate with the
   `FastBootDistIP` Jukes–Cantor transform and source ×3200 quantization.
4. In distance mode, sends that shared pair matrix directly to `GetPltVal2`-
   compatible triplet interpretation.
5. In UPGMA/NJ mode, constructs one tree from the complete active-cohort pair
   matrix, replaces pair distances with leaf-to-leaf topology paths, and only
   then interprets every requested triplet.
6. Retains cutoff/overlap topology runs and applies triplet-local
   `BSSubSeq`/`MakeScoresBS`/`ProbCalc` evidence.

The same mode-aware batch is called by first-pass discovery, affected-origin
redo scans after signal disassembly, candidate characterization, and edited-
event recalculation. No sequence is designated as a reference in exploratory
mode, and no alignment consensus substitutes for a third triplet member.

## Mode semantics

| Mode | Pair set | Relationship stored per replicate | Expected scaling |
| --- | --- | --- | --- |
| Distance | Union of pairs in requested triplets | Quantized pair distance | Pair/window/bootstrap work plus triplet interpretation |
| UPGMA | Every pair in the active cohort | Unweighted leaf path in a size-weighted UPGMA tree | Full pair work plus one O(N³) tree per window/replicate |
| Neighbor joining | Every pair in the active cohort | Unweighted leaf path in a classic NJ tree | Full pair work plus one O(N³) tree per window/replicate |

Tree scratch memory is serially reused. The transformed pair rows, not separate
three-taxon trees, are shared by all triplets. An approximate parent shortlist
can reduce the number of interpreted triplets but deliberately cannot reduce
the cohort matrix required by a tree mode.

## Performance path

Packed complete windows validate and compare sixteen sites per 32-bit word.
The bootstrap table is visited only for mismatch lanes. If any lane is missing,
the pair is restarted through the scalar complete/missing-data path. Distance
mode retains sparse pair compaction. UPGMA/NJ reuse one mutable tree workspace
and retain only the final pair-path rows needed by `GetPltVal2`.

Run the reproducible gate with:

```sh
npm run bench:bootscan:gate
```

## Regression boundary

Automated coverage includes:

- legacy distance-ABI equality;
- byte/packed equality in all modes;
- missing data across word and circular-origin boundaries;
- a known four-taxon UPGMA topology;
- randomized NJ equality against the independent split-path implementation;
- a cohort-context fixture in which a fourth taxon changes the stored
  relationship of the focal triplet;
- full-pair coverage from a single shortlisted tree-mode triplet;
- independent worker discovery in all modes;
- edited-event recalculation and project/CSV provenance.

The remaining scientific gate is an authorized desktop corpus for exact native
UPGMA/NJ tie, midpoint-root/tree-position, rounding, warning, and boundary
behavior. Until that comparison passes, 0.10.0 claims the three source-shaped
workflows and explicit numerical provenance, not bit-for-bit executable parity.
