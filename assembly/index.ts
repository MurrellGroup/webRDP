// WebAssembly kernels for RDP Web.
//
// Source-compatibility routines are ports of the RDP5 implementation supplied
// by the original RDP authors and are used here with their permission.  The
// small, stable ABI keeps the analysis worker independent of the legacy VB6 /
// Win32 runtime while preserving the numerical and event-delineation rules.

@inline
function valid(base: u8): bool {
  return base < 4;
}

@inline
function seqBase(seqPtr: i32, nSites: i32, sequence: i32, site: i32): u8 {
  return load<u8>(usize(seqPtr + sequence * nSites + site));
}

export function distance_matrix(
  seqPtr: i32,
  nSeq: i32,
  nSites: i32,
  outPtr: i32,
): void {
  for (let a: i32 = 0; a < nSeq; a += 1) {
    store<f32>(usize(outPtr + (a * nSeq + a) * 4), 0.0);
    for (let b: i32 = a + 1; b < nSeq; b += 1) {
      let validSites: i32 = 0;
      let differences: i32 = 0;
      for (let site: i32 = 0; site < nSites; site += 1) {
        const left = seqBase(seqPtr, nSites, a, site);
        const right = seqBase(seqPtr, nSites, b, site);
        if (!valid(left) || !valid(right)) continue;
        validSites += 1;
        if (left != right) differences += 1;
      }
      const distance: f32 = validSites > 0
        ? f32(differences) / f32(validSites)
        : 1.0;
      store<f32>(usize(outPtr + (a * nSeq + b) * 4), distance);
      store<f32>(usize(outPtr + (b * nSeq + a) * 4), distance);
    }
  }
}

// Canonical bases are packed 16 sites per u32 (two bits per base). The
// validity mask uses the low bit of each two-bit lane. This exact p-distance
// kernel moves 16x fewer sequence bytes than distance_matrix and is the
// default path used by the browser worker.
@inline
function popCount32(value: u32): i32 {
  return i32(popcnt<u32>(value));
}

export function distance_matrix_packed(
  packedPtr: i32,
  validityPtr: i32,
  nSeq: i32,
  wordsPerSequence: i32,
  outPtr: i32,
): void {
  for (let a: i32 = 0; a < nSeq; a += 1) {
    store<f32>(usize(outPtr + (a * nSeq + a) * 4), 0.0);
    for (let b: i32 = a + 1; b < nSeq; b += 1) {
      let validSites: i32 = 0;
      let differences: i32 = 0;
      const aOffset = a * wordsPerSequence;
      const bOffset = b * wordsPerSequence;
      for (let word: i32 = 0; word < wordsPerSequence; word += 1) {
        const left = load<u32>(usize(packedPtr + (aOffset + word) * 4));
        const right = load<u32>(usize(packedPtr + (bOffset + word) * 4));
        const bothValid = load<u32>(usize(validityPtr + (aOffset + word) * 4))
          & load<u32>(usize(validityPtr + (bOffset + word) * 4));
        const xor = left ^ right;
        const differingLanes = (xor | (xor >> 1)) & 0x55555555 & bothValid;
        validSites += popCount32(bothValid);
        differences += popCount32(differingLanes);
      }
      const distance: f32 = validSites > 0
        ? f32(differences) / f32(validSites)
        : 1.0;
      store<f32>(usize(outPtr + (a * nSeq + b) * 4), distance);
      store<f32>(usize(outPtr + (b * nSeq + a) * 4), distance);
    }
  }
}

// Exact maximum full-alignment p-distance without materialising an N x N
// matrix. RDP5 GetSupers normalises daughter divergence by the most divergent
// sequence pair in the input alignment; this streaming packed implementation
// preserves that rule for large datasets at O(N^2 L / 16) work and O(1)
// additional memory.
export function maximum_packed_distance(
  packedPtr: i32,
  validityPtr: i32,
  nSeq: i32,
  wordsPerSequence: i32,
): f32 {
  let maximum: f32 = 0.0;
  for (let a: i32 = 0; a < nSeq; a += 1) {
    for (let b: i32 = a + 1; b < nSeq; b += 1) {
      let validSites: i32 = 0;
      let differences: i32 = 0;
      const aOffset = a * wordsPerSequence;
      const bOffset = b * wordsPerSequence;
      for (let word: i32 = 0; word < wordsPerSequence; word += 1) {
        const left = load<u32>(usize(packedPtr + (aOffset + word) * 4));
        const right = load<u32>(usize(packedPtr + (bOffset + word) * 4));
        const bothValid = load<u32>(usize(validityPtr + (aOffset + word) * 4))
          & load<u32>(usize(validityPtr + (bOffset + word) * 4));
        const xor = left ^ right;
        const differingLanes = (xor | (xor >> 1)) & 0x55555555 & bothValid;
        validSites += popCount32(bothValid);
        differences += popCount32(differingLanes);
      }
      const distance: f32 = validSites > 0 ? f32(differences) / f32(validSites) : 1.0;
      if (distance > maximum) maximum = distance;
    }
  }
  return maximum;
}

// Packed implementation of the VisRD dMax statistic used by RDP5's
// recombinant-identification consensus.  The author-supplied C++ routine
// (CMaxD2P/CMaxD2P3) enumerates quartets containing at least one member of the
// event triplet, constructs the three quartet-map coordinates on the tract and
// its complement, and averages their L1 displacement for each triplet member.
//
// VScoreMat awards two units to a 2/2 split and one unit to the matching split
// of a three-state pattern (the C++ table stores 1 and 0.5).  Working on the
// low bit of each packed two-bit nucleotide lane lets one machine word score
// sixteen sites without changing that table's algebra.
@inline
function packedEquality(left: u32, right: u32, validLanes: u32): u32 {
  const xor = left ^ right;
  return ~(xor | (xor >> 1)) & 0x55555555 & validLanes;
}

@inline
function quartetSplit0(
  ab: u32,
  ac: u32,
  ad: u32,
  bc: u32,
  bd: u32,
  cd: u32,
  validLanes: u32,
): i32 {
  const full = ab & cd & ~ac & validLanes;
  const halfAB = ab & ~ac & ~ad & ~cd & validLanes;
  const halfCD = cd & ~ac & ~bc & ~ab & validLanes;
  return popCount32(full) * 2 + popCount32(halfAB) + popCount32(halfCD);
}

@inline
function quartetSplit1(
  ab: u32,
  ac: u32,
  ad: u32,
  bc: u32,
  bd: u32,
  cd: u32,
  validLanes: u32,
): i32 {
  const full = ac & bd & ~ab & validLanes;
  const halfAC = ac & ~ab & ~ad & ~bd & validLanes;
  const halfBD = bd & ~ab & ~bc & ~ac & validLanes;
  return popCount32(full) * 2 + popCount32(halfAC) + popCount32(halfBD);
}

@inline
function quartetSplit2(
  ab: u32,
  ac: u32,
  ad: u32,
  bc: u32,
  bd: u32,
  cd: u32,
  validLanes: u32,
): i32 {
  const full = ad & bc & ~ab & validLanes;
  const halfAD = ad & ~ab & ~ac & ~bc & validLanes;
  const halfBC = bc & ~ab & ~bd & ~ad & validLanes;
  return popCount32(full) * 2 + popCount32(halfAD) + popCount32(halfBC);
}

export function dmax_visrd_packed(
  packedPtr: i32,
  validityPtr: i32,
  wordsPerSequence: i32,
  cohortPtr: i32,
  cohortCount: i32,
  candidate0: i32,
  candidate1: i32,
  candidate2: i32,
  tractMaskPtr: i32,
  backgroundMaskPtr: i32,
  outPtr: i32,
): void {
  let sum0: f64 = 0.0;
  let sum1: f64 = 0.0;
  let sum2: f64 = 0.0;
  let count0: i32 = 0;
  let count1: i32 = 0;
  let count2: i32 = 0;
  if (cohortCount >= 4) {
    for (let ai: i32 = 0; ai < cohortCount - 3; ai += 1) {
      const a = load<i32>(usize(cohortPtr + ai * 4));
      const aOffset = a * wordsPerSequence;
      for (let bi: i32 = ai + 1; bi < cohortCount - 2; bi += 1) {
        const b = load<i32>(usize(cohortPtr + bi * 4));
        const bOffset = b * wordsPerSequence;
        for (let ci: i32 = bi + 1; ci < cohortCount - 1; ci += 1) {
          const c = load<i32>(usize(cohortPtr + ci * 4));
          const cOffset = c * wordsPerSequence;
          for (let di: i32 = ci + 1; di < cohortCount; di += 1) {
            const d = load<i32>(usize(cohortPtr + di * 4));
            const has0 = a == candidate0 || b == candidate0 || c == candidate0 || d == candidate0;
            const has1 = a == candidate1 || b == candidate1 || c == candidate1 || d == candidate1;
            const has2 = a == candidate2 || b == candidate2 || c == candidate2 || d == candidate2;
            if (!has0 && !has1 && !has2) continue;
            const dOffset = d * wordsPerSequence;
            let outside0: i32 = 0;
            let outside1: i32 = 0;
            let outside2: i32 = 0;
            let inside0: i32 = 0;
            let inside1: i32 = 0;
            let inside2: i32 = 0;
            for (let word: i32 = 0; word < wordsPerSequence; word += 1) {
              const av = load<u32>(usize(packedPtr + (aOffset + word) * 4));
              const bv = load<u32>(usize(packedPtr + (bOffset + word) * 4));
              const cv = load<u32>(usize(packedPtr + (cOffset + word) * 4));
              const dv = load<u32>(usize(packedPtr + (dOffset + word) * 4));
              const allValid = load<u32>(usize(validityPtr + (aOffset + word) * 4))
                & load<u32>(usize(validityPtr + (bOffset + word) * 4))
                & load<u32>(usize(validityPtr + (cOffset + word) * 4))
                & load<u32>(usize(validityPtr + (dOffset + word) * 4));
              if (allValid == 0) continue;
              const ab = packedEquality(av, bv, allValid);
              const ac = packedEquality(av, cv, allValid);
              const ad = packedEquality(av, dv, allValid);
              const bc = packedEquality(bv, cv, allValid);
              const bd = packedEquality(bv, dv, allValid);
              const cd = packedEquality(cv, dv, allValid);
              const tract = allValid & load<u32>(usize(tractMaskPtr + word * 4));
              const background = allValid & load<u32>(usize(backgroundMaskPtr + word * 4));
              if (tract != 0) {
                inside0 += quartetSplit0(ab, ac, ad, bc, bd, cd, tract);
                inside1 += quartetSplit1(ab, ac, ad, bc, bd, cd, tract);
                inside2 += quartetSplit2(ab, ac, ad, bc, bd, cd, tract);
              }
              if (background != 0) {
                outside0 += quartetSplit0(ab, ac, ad, bc, bd, cd, background);
                outside1 += quartetSplit1(ab, ac, ad, bc, bd, cd, background);
                outside2 += quartetSplit2(ab, ac, ad, bc, bd, cd, background);
              }
            }
            const insideTotal = inside0 + inside1 + inside2;
            const outsideTotal = outside0 + outside1 + outside2;
            if (insideTotal <= 0 || outsideTotal <= 0) continue;
            const distance = Math.abs(f64(inside0) / f64(insideTotal) - f64(outside0) / f64(outsideTotal))
              + Math.abs(f64(inside1) / f64(insideTotal) - f64(outside1) / f64(outsideTotal))
              + Math.abs(f64(inside2) / f64(insideTotal) - f64(outside2) / f64(outsideTotal));
            if (has0) { sum0 += distance; count0 += 1; }
            if (has1) { sum1 += distance; count1 += 1; }
            if (has2) { sum2 += distance; count2 += 1; }
          }
        }
      }
    }
  }
  store<f64>(usize(outPtr), count0 > 0 ? sum0 / f64(count0) : 0.0);
  store<f64>(usize(outPtr + 8), count1 > 0 ? sum1 / f64(count1) : 0.0);
  store<f64>(usize(outPtr + 16), count2 > 0 ? sum2 / f64(count2) : 0.0);
  store<i32>(usize(outPtr + 24), count0);
  store<i32>(usize(outPtr + 28), count1);
  store<i32>(usize(outPtr + 32), count2);
}

// For very large alignments, materializing and calculating an N×N full-length
// distance matrix is the wrong asymptotic tradeoff. This kernel evaluates a
// deterministic, evenly-spaced site sketch against a supplied reference pool
// and retains only the nearest K candidates in insertion-sorted output memory.
// Definitive scans can still request the exact full matrix/exhaustive path.
export function nearest_candidates_sampled(
  seqPtr: i32,
  nSites: i32,
  target: i32,
  poolPtr: i32,
  poolCount: i32,
  requestedSamples: i32,
  limit: i32,
  outIndexesPtr: i32,
  outDistancesPtr: i32,
): i32 {
  if (limit <= 0 || poolCount <= 0 || nSites <= 0) return 0;
  const sampleCount = requestedSamples < nSites ? requestedSamples : nSites;
  let retained: i32 = 0;
  for (let poolPosition: i32 = 0; poolPosition < poolCount; poolPosition += 1) {
    const candidate = load<i32>(usize(poolPtr + poolPosition * 4));
    if (candidate == target) continue;
    let validSites: i32 = 0;
    let differences: i32 = 0;
    for (let sample: i32 = 0; sample < sampleCount; sample += 1) {
      const site = sampleCount <= 1
        ? 0
        : i32((i64(sample) * i64(nSites - 1)) / i64(sampleCount - 1));
      const left = seqBase(seqPtr, nSites, target, site);
      const right = seqBase(seqPtr, nSites, candidate, site);
      if (!valid(left) || !valid(right)) continue;
      validSites += 1;
      if (left != right) differences += 1;
    }
    const distance = validSites > 0
      ? i32((i64(differences) * 1000000) / i64(validSites))
      : 1000000;
    let insertion = retained;
    for (let index: i32 = 0; index < retained; index += 1) {
      const existingDistance = load<i32>(usize(outDistancesPtr + index * 4));
      const existingIndex = load<i32>(usize(outIndexesPtr + index * 4));
      if (distance < existingDistance || (distance == existingDistance && candidate < existingIndex)) {
        insertion = index;
        break;
      }
    }
    if (insertion >= limit) continue;
    const nextRetained = retained < limit ? retained + 1 : retained;
    for (let move = nextRetained - 1; move > insertion; move -= 1) {
      store<i32>(usize(outIndexesPtr + move * 4), load<i32>(usize(outIndexesPtr + (move - 1) * 4)));
      store<i32>(usize(outDistancesPtr + move * 4), load<i32>(usize(outDistancesPtr + (move - 1) * 4)));
    }
    store<i32>(usize(outIndexesPtr + insertion * 4), candidate);
    store<i32>(usize(outDistancesPtr + insertion * 4), distance);
    retained = nextRetained;
  }
  return retained;
}

@inline
function chiSquare(a: i32, b: i32, c: i32, d: i32): f64 {
  const n: f64 = f64(a + b + c + d);
  if (n <= 0.0) return 0.0;
  const numerator: f64 = n * f64(a * d - b * c) * f64(a * d - b * c);
  const denominator: f64 = f64((a + b) * (c + d) * (a + c) * (b + d));
  return denominator > 0.0 ? numerator / denominator : 0.0;
}

@inline
function clamp(value: i32, lower: i32, upper: i32): i32 {
  return value < lower ? lower : value > upper ? upper : value;
}

@inline
function prefixValue(ptr: i32, position: i32): i32 {
  return load<i32>(usize(ptr + position * 4));
}

@inline
function breakpointChi(
  prefixAPtr: i32,
  prefixBPtr: i32,
  nSites: i32,
  position: i32,
  halfWindow: i32,
): f64 {
  const left = clamp(position - halfWindow, 0, nSites);
  const center = clamp(position, 0, nSites);
  const right = clamp(position + halfWindow, 0, nSites);
  const leftA = prefixValue(prefixAPtr, center) - prefixValue(prefixAPtr, left);
  const leftB = prefixValue(prefixBPtr, center) - prefixValue(prefixBPtr, left);
  const rightA = prefixValue(prefixAPtr, right) - prefixValue(prefixAPtr, center);
  const rightB = prefixValue(prefixBPtr, right) - prefixValue(prefixBPtr, center);
  if (leftA + leftB < 4 || rightA + rightB < 4) return 0.0;
  return chiSquare(leftA, leftB, rightA, rightB);
}

@inline
function informativePrefix(prefixAPtr: i32, prefixBPtr: i32, position: i32): i32 {
  return prefixValue(prefixAPtr, position) + prefixValue(prefixBPtr, position);
}

// RDP5's automated MAXCHI and CHIMAERA scans move in compressed
// informative-site coordinates, not raw alignment coordinates. This lower
// bound maps the boundary after `rank` informative sites back to the original
// alignment without materialising a second coordinate array.
function positionForInformativeRank(
  prefixAPtr: i32,
  prefixBPtr: i32,
  nSites: i32,
  rank: i32,
): i32 {
  if (rank <= 0) return 0;
  const total = informativePrefix(prefixAPtr, prefixBPtr, nSites);
  if (rank >= total) return nSites;
  let lower: i32 = 0;
  let upper: i32 = nSites;
  while (lower < upper) {
    const middle = lower + (upper - lower) / 2;
    if (informativePrefix(prefixAPtr, prefixBPtr, middle) < rank) lower = middle + 1;
    else upper = middle;
  }
  return lower;
}

function variableBreakpointChi(
  prefixAPtr: i32,
  prefixBPtr: i32,
  nSites: i32,
  rank: i32,
  halfWindow: i32,
): f64 {
  const total = informativePrefix(prefixAPtr, prefixBPtr, nSites);
  if (rank < halfWindow || rank + halfWindow > total) return 0.0;
  const left = positionForInformativeRank(prefixAPtr, prefixBPtr, nSites, rank - halfWindow);
  const center = positionForInformativeRank(prefixAPtr, prefixBPtr, nSites, rank);
  const right = positionForInformativeRank(prefixAPtr, prefixBPtr, nSites, rank + halfWindow);
  const leftA = prefixValue(prefixAPtr, center) - prefixValue(prefixAPtr, left);
  const leftB = prefixValue(prefixBPtr, center) - prefixValue(prefixBPtr, left);
  const rightA = prefixValue(prefixAPtr, right) - prefixValue(prefixAPtr, center);
  const rightB = prefixValue(prefixBPtr, right) - prefixValue(prefixBPtr, center);
  return chiSquare(leftA, leftB, rightA, rightB);
}

function polishBreakpoint(
  prefixAPtr: i32,
  prefixBPtr: i32,
  nSites: i32,
  proposed: i32,
  halfWindow: i32,
): i32 {
  const radius = halfWindow > 12 ? halfWindow : 12;
  const first = clamp(proposed - radius, 1, nSites - 1);
  const last = clamp(proposed + radius, 1, nSites - 1);
  let bestPosition = clamp(proposed, 1, nSites - 1);
  let best = breakpointChi(prefixAPtr, prefixBPtr, nSites, bestPosition, halfWindow);
  for (let position = first; position <= last; position += 1) {
    const statistic = breakpointChi(prefixAPtr, prefixBPtr, nSites, position, halfWindow);
    if (statistic > best) {
      best = statistic;
      bestPosition = position;
    }
  }
  return bestPosition;
}

// Evaluates a user-specified recombinant/parent assignment and tract without
// running candidate discovery. This keeps interactive reassignment and manual
// breakpoint edits scientifically honest: the browser can recompute evidence
// for the exact hypothesis instead of displaying statistics from its prior
// assignment.
// Output: informative, inside-minor, inside-major, outside-major,
// outside-minor, chi-square*1000.
export function triplet_counts(
  seqPtr: i32,
  nSites: i32,
  recombinant: i32,
  majorParent: i32,
  minorParent: i32,
  start: i32,
  end: i32,
  outPtr: i32,
): void {
  let informative: i32 = 0;
  let insideMinor: i32 = 0;
  let insideMajor: i32 = 0;
  let outsideMajor: i32 = 0;
  let outsideMinor: i32 = 0;
  for (let site: i32 = 0; site < nSites; site += 1) {
    const r = seqBase(seqPtr, nSites, recombinant, site);
    const major = seqBase(seqPtr, nSites, majorParent, site);
    const minor = seqBase(seqPtr, nSites, minorParent, site);
    if (!valid(r) || !valid(major) || !valid(minor) || major == minor) continue;
    let matchesMajor = r == major;
    let matchesMinor = r == minor;
    if (!matchesMajor && !matchesMinor) continue;
    informative += 1;
    if (site >= start && site < end) {
      if (matchesMinor) insideMinor += 1;
      if (matchesMajor) insideMajor += 1;
    } else {
      if (matchesMajor) outsideMajor += 1;
      if (matchesMinor) outsideMinor += 1;
    }
  }
  const statistic = chiSquare(insideMinor, insideMajor, outsideMinor, outsideMajor);
  store<i32>(usize(outPtr), informative);
  store<i32>(usize(outPtr + 4), insideMinor);
  store<i32>(usize(outPtr + 8), insideMajor);
  store<i32>(usize(outPtr + 12), outsideMajor);
  store<i32>(usize(outPtr + 16), outsideMinor);
  store<i32>(usize(outPtr + 20), i32(statistic * 1000.0));
}

// Calculates method-specific fast statistics for the seven exploratory method
// families after a candidate tract has been localized.  These kernels are
// progressively being replaced by source-compatible method ports; the ABI is
// retained so saved projects and interactive recalculation remain stable.
//
// Output layout (i32 values):
//  0..4  GENECONV G=0 run, eligible/matching sites, run boundaries
//  5..6  BOOTSCAN topology-consistent / decisive windows
//  7..8  MAXCHI and CHIMAERA breakpoint statistics, scaled by 1000
//  9..10 SISCAN-oriented category score and informative-site count
// 11..12 3SEQ maximum HGRW descent and informative-site count
// 13..16 individual MAXCHI/CHIMAERA boundary statistics, scaled by 1000
// 17..18 polished left/right breakpoint positions
// 19..20 3SEQ major-parent and minor-parent informative-site counts
// 21..22 seeded p-distance bootstrap-consistent / decisive replicates
// 23..24 3SEQ maximum-descent start/end genomic positions
// 25..26 independent MAXCHI peak-pair positions
// 27..28 independent CHIMAERA peak-pair positions
// 29..30 independent BOOTSCAN minor-topology run positions
// 31..32 independent SISCAN oriented-category run positions
// 33..34 independent BOOTSCAN/SISCAN run lengths in sampled windows
export function method_stats(
  seqPtr: i32,
  nSites: i32,
  recombinant: i32,
  majorParent: i32,
  minorParent: i32,
  start: i32,
  end: i32,
  window: i32,
  step: i32,
  bootstrapReplicates: i32,
  randomSeed: i32,
  methodMask: i32,
  prefixAPtr: i32,
  prefixBPtr: i32,
  outPtr: i32,
  geneconvGScale: f64,
): void {
  const halfWindow = window / 2 > 8 ? window / 2 : 8;
  const wantGeneconv = (methodMask & 1) != 0;
  const wantBootscan = (methodMask & 2) != 0;
  const wantMaxChi = (methodMask & 4) != 0;
  const wantChimaera = (methodMask & 8) != 0;
  const wantSiScan = (methodMask & 16) != 0;
  const wantThreeSeq = (methodMask & 32) != 0;
  const wantPolish = (methodMask & 64) != 0;

  // RDP5 GENECONV fragment scoring in triplet-polymorphic space. G=0 is the
  // infinite-penalty exact-run special case. Finite G uses the source integer
  // mismatch penalty floor(L*G/mismatches)+1 and maximum local fragment score.
  let eligible: i32 = 0;
  let concordant: i32 = 0;
  let run: i32 = 0;
  let bestRun: i32 = 0;
  let runStart: i32 = start;
  let bestRunStart: i32 = start;
  let bestRunEnd: i32 = start;
  if (wantGeneconv) {
    for (let site: i32 = 0; site < nSites; site += 1) {
      const r = seqBase(seqPtr, nSites, recombinant, site);
      const a = seqBase(seqPtr, nSites, majorParent, site);
      const b = seqBase(seqPtr, nSites, minorParent, site);
      if (!valid(r) || !valid(a) || !valid(b) || (r == a && r == b)) continue;
      eligible += 1;
      if (r == b) concordant += 1;
      if (geneconvGScale <= 0.0) {
        if (r == b) {
          if (run == 0) runStart = site;
          run += 1;
          if (run > bestRun) {
            bestRun = run;
            bestRunStart = runStart;
            bestRunEnd = site + 1;
          }
        } else {
          run = 0;
        }
      }
    }
    if (geneconvGScale > 0.0) {
      const mismatches = eligible - concordant;
      if (mismatches > 0 && concordant > 0) {
        const mismatchPenalty = i32(Math.floor(f64(eligible) * geneconvGScale / f64(mismatches))) + 1;
        let fragmentScore: i32 = 0;
        for (let site: i32 = 0; site < nSites; site += 1) {
          const r = seqBase(seqPtr, nSites, recombinant, site);
          const a = seqBase(seqPtr, nSites, majorParent, site);
          const b = seqBase(seqPtr, nSites, minorParent, site);
          if (!valid(r) || !valid(a) || !valid(b) || (r == a && r == b)) continue;
          const siteScore = r == b ? 1 : -mismatchPenalty;
          if (fragmentScore <= 0) {
            if (siteScore > 0) {
              fragmentScore = siteScore;
              runStart = site;
            } else {
              fragmentScore = 0;
            }
          } else {
            fragmentScore += siteScore;
            if (fragmentScore < 0) fragmentScore = 0;
          }
          if (fragmentScore > bestRun) {
            bestRun = fragmentScore;
            bestRunStart = runStart;
            bestRunEnd = site + 1;
          }
        }
      }
    }
  }

  // Prefix pairwise mismatches for a deterministic, no-resampling RECSCAN
  // topology-switch statistic. Only sites valid in the full triplet are used
  // so both distances have identical denominators.
  let majorDifferences: i32 = 0;
  let minorDifferences: i32 = 0;
  if (wantBootscan) {
    store<i32>(usize(prefixAPtr), 0);
    store<i32>(usize(prefixBPtr), 0);
    for (let site: i32 = 0; site < nSites; site += 1) {
      const r = seqBase(seqPtr, nSites, recombinant, site);
      const a = seqBase(seqPtr, nSites, majorParent, site);
      const b = seqBase(seqPtr, nSites, minorParent, site);
      if (valid(r) && valid(a) && valid(b)) {
        if (r != a) majorDifferences += 1;
        if (r != b) minorDifferences += 1;
      }
      store<i32>(usize(prefixAPtr + (site + 1) * 4), majorDifferences);
      store<i32>(usize(prefixBPtr + (site + 1) * 4), minorDifferences);
    }
  }
  let topologyConsistent: i32 = 0;
  let decisiveWindows: i32 = 0;
  let bootscanRunStart: i32 = start;
  let bootscanRunEnd: i32 = end;
  let currentBootscanStart: i32 = -1;
  let currentBootscanWindows: i32 = 0;
  let bestBootscanWindows: i32 = 0;
  const stride = step > 0 ? step : 1;
  if (wantBootscan) for (let center = halfWindow; center <= nSites - halfWindow; center += stride) {
    const left = center - halfWindow;
    const right = center + halfWindow;
    const majorDistance = prefixValue(prefixAPtr, right) - prefixValue(prefixAPtr, left);
    const minorDistance = prefixValue(prefixBPtr, right) - prefixValue(prefixBPtr, left);
    if (majorDistance == minorDistance) {
      currentBootscanStart = -1;
      currentBootscanWindows = 0;
      continue;
    }
    decisiveWindows += 1;
    const inside = center >= start && center < end;
    if ((inside && minorDistance < majorDistance) || (!inside && majorDistance < minorDistance)) {
      topologyConsistent += 1;
    }
    if (minorDistance < majorDistance) {
      if (currentBootscanStart < 0) currentBootscanStart = center;
      currentBootscanWindows += 1;
      if (currentBootscanWindows > bestBootscanWindows) {
        bestBootscanWindows = currentBootscanWindows;
        bootscanRunStart = clamp(currentBootscanStart - halfWindow, 0, nSites);
        bootscanRunEnd = clamp(center + halfWindow, 0, nSites);
      }
    } else {
      currentBootscanStart = -1;
      currentBootscanWindows = 0;
    }
  }

  // A reproducible triplet bootstrap complements the fast all-window sign
  // scan. Columns are resampled with replacement from one tract window and
  // the available left/right flanks. For a three-sequence candidate, the
  // p-distance ordering is the tree-topology decision relevant to the local
  // parent switch. This avoids JavaScript RNG overhead and keeps the cost
  // bounded at O(3 * replicates * window).
  let bootstrapConsistent: i32 = 0;
  let bootstrapDecisive: i32 = 0;
  const requestedReplicates = clamp(bootstrapReplicates, 0, 1000);
  const targetWindow = clamp(window, 12, nSites);
  let state: u32 = u32(randomSeed) ^ 0x9e3779b9;
  if (state == 0) state = 0x6d2b79f5;
  if (wantBootscan) for (let region: i32 = 0; region < 3; region += 1) {
    let regionStart: i32 = 0;
    let regionLength: i32 = 0;
    let expectMinor = false;
    if (region == 0) {
      regionLength = start < targetWindow ? start : targetWindow;
      regionStart = start - regionLength;
    } else if (region == 1) {
      const tractLength = end - start;
      regionLength = tractLength < targetWindow ? tractLength : targetWindow;
      regionStart = start + (tractLength - regionLength) / 2;
      expectMinor = true;
    } else {
      const available = nSites - end;
      regionLength = available < targetWindow ? available : targetWindow;
      regionStart = end;
    }
    if (regionLength < 8) continue;
    for (let replicate: i32 = 0; replicate < requestedReplicates; replicate += 1) {
      let majorDistance: i32 = 0;
      let minorDistance: i32 = 0;
      let validSites: i32 = 0;
      for (let draw: i32 = 0; draw < regionLength; draw += 1) {
        state ^= state << 13;
        state ^= state >> 17;
        state ^= state << 5;
        const site = regionStart + i32(state % u32(regionLength));
        const r = seqBase(seqPtr, nSites, recombinant, site);
        const a = seqBase(seqPtr, nSites, majorParent, site);
        const b = seqBase(seqPtr, nSites, minorParent, site);
        if (!valid(r) || !valid(a) || !valid(b)) continue;
        validSites += 1;
        if (r != a) majorDistance += 1;
        if (r != b) minorDistance += 1;
      }
      if (validSites < 4 || majorDistance == minorDistance) continue;
      bootstrapDecisive += 1;
      if ((expectMinor && minorDistance < majorDistance)
        || (!expectMinor && majorDistance < minorDistance)) bootstrapConsistent += 1;
    }
  }

  // MAXCHI pairwise variable/non-variable states in triplet-polymorphic space.
  let pairVariable: i32 = 0;
  let pairNonVariable: i32 = 0;
  if (wantMaxChi) {
    store<i32>(usize(prefixAPtr), 0);
    store<i32>(usize(prefixBPtr), 0);
    for (let site: i32 = 0; site < nSites; site += 1) {
      const r = seqBase(seqPtr, nSites, recombinant, site);
      const a = seqBase(seqPtr, nSites, majorParent, site);
      const b = seqBase(seqPtr, nSites, minorParent, site);
      if (valid(r) && valid(a) && valid(b) && (r != a || r != b || a != b)) {
        if (r == b) pairNonVariable += 1;
        else pairVariable += 1;
      }
      store<i32>(usize(prefixAPtr + (site + 1) * 4), pairVariable);
      store<i32>(usize(prefixBPtr + (site + 1) * 4), pairNonVariable);
    }
  }
  let maxChiStart: f64 = 0.0;
  let maxChiEnd: f64 = 0.0;
  let maxChiStatistic: f64 = 0.0;
  let maxChiStartPosition: i32 = start;
  let maxChiEndPosition: i32 = end;
  if (wantMaxChi) {
    const informative = informativePrefix(prefixAPtr, prefixBPtr, nSites);
    let firstPeak: f64 = -1.0;
    let firstPeakRank: i32 = halfWindow;
    for (let rank: i32 = halfWindow; rank + halfWindow <= informative; rank += 1) {
      const value = variableBreakpointChi(prefixAPtr, prefixBPtr, nSites, rank, halfWindow);
      if (value > firstPeak) {
        firstPeak = value;
        firstPeakRank = rank;
        maxChiStartPosition = positionForInformativeRank(prefixAPtr, prefixBPtr, nSites, rank);
      }
    }
    let secondPeak: f64 = -1.0;
    const exclusion = halfWindow > 8 ? halfWindow : 8;
    for (let rank: i32 = halfWindow; rank + halfWindow <= informative; rank += 1) {
      if (Math.abs(rank - firstPeakRank) < exclusion) continue;
      const value = variableBreakpointChi(prefixAPtr, prefixBPtr, nSites, rank, halfWindow);
      if (value > secondPeak) {
        secondPeak = value;
        maxChiEndPosition = positionForInformativeRank(prefixAPtr, prefixBPtr, nSites, rank);
      }
    }
    if (maxChiEndPosition < maxChiStartPosition) {
      const swap = maxChiStartPosition;
      maxChiStartPosition = maxChiEndPosition;
      maxChiEndPosition = swap;
    }
    maxChiStart = firstPeak > 0.0 ? firstPeak : 0.0;
    maxChiEnd = secondPeak > 0.0 ? secondPeak : 0.0;
    maxChiStatistic = maxChiStart < maxChiEnd ? maxChiStart : maxChiEnd;
  }

  // CHIMAERA's compressed binary string and SISCAN-style oriented category
  // score. The same prefixes support an O(window) breakpoint-polishing pass.
  let matchesMajor: i32 = 0;
  let matchesMinor: i32 = 0;
  let sisterScore: i32 = 0;
  let sisterSites: i32 = 0;
  let siScanRunStart: i32 = start;
  let siScanRunEnd: i32 = end;
  let bestSiScanWindows: i32 = 0;
  if (wantChimaera || wantSiScan || wantThreeSeq || wantPolish) {
    store<i32>(usize(prefixAPtr), 0);
    store<i32>(usize(prefixBPtr), 0);
    for (let site: i32 = 0; site < nSites; site += 1) {
      const r = seqBase(seqPtr, nSites, recombinant, site);
      const a = seqBase(seqPtr, nSites, majorParent, site);
      const b = seqBase(seqPtr, nSites, minorParent, site);
      if (valid(r) && valid(a) && valid(b) && a != b) {
        let category: i32 = 0;
        if (r == a) {
          matchesMajor += 1;
          category = -1;
        } else if (r == b) {
          matchesMinor += 1;
          category = 1;
        }
        if (category != 0 && wantSiScan) {
          sisterSites += 1;
          sisterScore += (site >= start && site < end) ? category : -category;
        }
      }
      store<i32>(usize(prefixAPtr + (site + 1) * 4), matchesMajor);
      store<i32>(usize(prefixBPtr + (site + 1) * 4), matchesMinor);
    }
  }
  if (wantSiScan) {
    let currentStart: i32 = -1;
    let currentWindows: i32 = 0;
    for (let center = halfWindow; center <= nSites - halfWindow; center += stride) {
      const left = center - halfWindow;
      const right = center + halfWindow;
      const majorCount = prefixValue(prefixAPtr, right) - prefixValue(prefixAPtr, left);
      const minorCount = prefixValue(prefixBPtr, right) - prefixValue(prefixBPtr, left);
      if (minorCount > majorCount) {
        if (currentStart < 0) currentStart = center;
        currentWindows += 1;
        if (currentWindows > bestSiScanWindows) {
          bestSiScanWindows = currentWindows;
          siScanRunStart = clamp(currentStart - halfWindow, 0, nSites);
          siScanRunEnd = clamp(center + halfWindow, 0, nSites);
        }
      } else {
        currentStart = -1;
        currentWindows = 0;
      }
    }
  }
  let chimaeraStart: f64 = 0.0;
  let chimaeraEnd: f64 = 0.0;
  let chimaeraStatistic: f64 = 0.0;
  let chimaeraStartPosition: i32 = start;
  let chimaeraEndPosition: i32 = end;
  if (wantChimaera) {
    const informative = informativePrefix(prefixAPtr, prefixBPtr, nSites);
    let firstPeak: f64 = -1.0;
    let firstPeakRank: i32 = halfWindow;
    for (let rank: i32 = halfWindow; rank + halfWindow <= informative; rank += 1) {
      const value = variableBreakpointChi(prefixAPtr, prefixBPtr, nSites, rank, halfWindow);
      if (value > firstPeak) {
        firstPeak = value;
        firstPeakRank = rank;
        chimaeraStartPosition = positionForInformativeRank(prefixAPtr, prefixBPtr, nSites, rank);
      }
    }
    let secondPeak: f64 = -1.0;
    const exclusion = halfWindow > 8 ? halfWindow : 8;
    for (let rank: i32 = halfWindow; rank + halfWindow <= informative; rank += 1) {
      if (Math.abs(rank - firstPeakRank) < exclusion) continue;
      const value = variableBreakpointChi(prefixAPtr, prefixBPtr, nSites, rank, halfWindow);
      if (value > secondPeak) {
        secondPeak = value;
        chimaeraEndPosition = positionForInformativeRank(prefixAPtr, prefixBPtr, nSites, rank);
      }
    }
    if (chimaeraEndPosition < chimaeraStartPosition) {
      const swap = chimaeraStartPosition;
      chimaeraStartPosition = chimaeraEndPosition;
      chimaeraEndPosition = swap;
    }
    chimaeraStart = firstPeak > 0.0 ? firstPeak : 0.0;
    chimaeraEnd = secondPeak > 0.0 ? secondPeak : 0.0;
    chimaeraStatistic = chimaeraStart < chimaeraEnd ? chimaeraStart : chimaeraEnd;
  }
  // 3SEQ maps sites matching the major parent to an up-step and sites
  // matching the minor parent to a down-step. Its two-breakpoint statistic is
  // the maximum descent of that hypergeometric random walk from any previous
  // maximum, corresponding to the strongest contiguous minor-parent tract.
  const threeSeqSites = matchesMajor + matchesMinor;
  let threeSeqWalk: i32 = 0;
  let threeSeqMaximum: i32 = 0;
  let threeSeqDescent: i32 = 0;
  let threeSeqMaximumPosition: i32 = 0;
  let threeSeqStart: i32 = start;
  let threeSeqEnd: i32 = end;
  if (wantThreeSeq) for (let position: i32 = 1; position <= nSites; position += 1) {
    const majorStep = prefixValue(prefixAPtr, position) - prefixValue(prefixAPtr, position - 1);
    const minorStep = prefixValue(prefixBPtr, position) - prefixValue(prefixBPtr, position - 1);
    threeSeqWalk += majorStep - minorStep;
    if (threeSeqWalk > threeSeqMaximum) {
      threeSeqMaximum = threeSeqWalk;
      threeSeqMaximumPosition = position;
    }
    const descent = threeSeqMaximum - threeSeqWalk;
    if (descent > threeSeqDescent) {
      threeSeqDescent = descent;
      threeSeqStart = threeSeqMaximumPosition;
      threeSeqEnd = position;
    }
  }
  const polishedStart = wantPolish ? polishBreakpoint(prefixAPtr, prefixBPtr, nSites, start, halfWindow) : start;
  const polishedEnd = wantPolish ? polishBreakpoint(prefixAPtr, prefixBPtr, nSites, end, halfWindow) : end;

  store<i32>(usize(outPtr), bestRun);
  store<i32>(usize(outPtr + 4), eligible);
  store<i32>(usize(outPtr + 8), concordant);
  store<i32>(usize(outPtr + 12), bestRunStart);
  store<i32>(usize(outPtr + 16), bestRunEnd);
  store<i32>(usize(outPtr + 20), topologyConsistent);
  store<i32>(usize(outPtr + 24), decisiveWindows);
  store<i32>(usize(outPtr + 28), i32(maxChiStatistic * 1000.0));
  store<i32>(usize(outPtr + 32), i32(chimaeraStatistic * 1000.0));
  store<i32>(usize(outPtr + 36), sisterScore);
  store<i32>(usize(outPtr + 40), sisterSites);
  store<i32>(usize(outPtr + 44), threeSeqDescent);
  store<i32>(usize(outPtr + 48), threeSeqSites);
  store<i32>(usize(outPtr + 52), i32(maxChiStart * 1000.0));
  store<i32>(usize(outPtr + 56), i32(maxChiEnd * 1000.0));
  store<i32>(usize(outPtr + 60), i32(chimaeraStart * 1000.0));
  store<i32>(usize(outPtr + 64), i32(chimaeraEnd * 1000.0));
  store<i32>(usize(outPtr + 68), polishedStart);
  store<i32>(usize(outPtr + 72), polishedEnd);
  store<i32>(usize(outPtr + 76), matchesMajor);
  store<i32>(usize(outPtr + 80), matchesMinor);
  store<i32>(usize(outPtr + 84), bootstrapConsistent);
  store<i32>(usize(outPtr + 88), bootstrapDecisive);
  store<i32>(usize(outPtr + 92), threeSeqStart);
  store<i32>(usize(outPtr + 96), threeSeqEnd);
  store<i32>(usize(outPtr + 100), maxChiStartPosition);
  store<i32>(usize(outPtr + 104), maxChiEndPosition);
  store<i32>(usize(outPtr + 108), chimaeraStartPosition);
  store<i32>(usize(outPtr + 112), chimaeraEndPosition);
  store<i32>(usize(outPtr + 116), bootscanRunStart);
  store<i32>(usize(outPtr + 120), bootscanRunEnd);
  store<i32>(usize(outPtr + 124), siScanRunStart);
  store<i32>(usize(outPtr + 128), siScanRunEnd);
  store<i32>(usize(outPtr + 132), bestBootscanWindows);
  store<i32>(usize(outPtr + 136), bestSiScanWindows);
}

@inline
function hmmProbability(value: f64): f64 {
  return value < 0.5001 ? 0.5001 : value > 0.9999 ? 0.9999 : value;
}

@inline
function markedPosition(value: i32): i32 {
  return value < 0 ? -value - 1 : value;
}

// Legacy two-state breakpoint-refinement ABI retained for older tests/projects.
// New analyses use the source-compatible BURT engine in public/rdp-burt.js;
// this O(L) kernel is no longer called by the worker.
//
// Output layout (i32 values):
// 0..1 refined start/end; 2 retained informative sites; 3 state switches;
// 4..5 major/minor correct-emission probabilities scaled by 1000;
// 6 informative sites in the selected minor-state run;
// 7..10 start-low/start-high/end-low/end-high informative-site bounds.
export function hmm_polish(
  seqPtr: i32,
  nSites: i32,
  recombinant: i32,
  majorParent: i32,
  minorParent: i32,
  start: i32,
  end: i32,
  positionPtr: i32,
  predecessorPtr: i32,
  outPtr: i32,
): i32 {
  if (nSites < 4 || start < 0 || end > nSites || end <= start) return 0;

  let informative: i32 = 0;
  let insideMajor: i32 = 0;
  let insideMinor: i32 = 0;
  let outsideMajor: i32 = 0;
  let outsideMinor: i32 = 0;
  for (let site: i32 = 0; site < nSites; site += 1) {
    const r = seqBase(seqPtr, nSites, recombinant, site);
    const a = seqBase(seqPtr, nSites, majorParent, site);
    const b = seqBase(seqPtr, nSites, minorParent, site);
    if (!valid(r) || !valid(a) || !valid(b) || a == b || (r != a && r != b)) continue;
    store<i32>(usize(positionPtr + informative * 4), site);
    const inside = site >= start && site < end;
    if (r == a) {
      if (inside) insideMajor += 1;
      else outsideMajor += 1;
    } else {
      if (inside) insideMinor += 1;
      else outsideMinor += 1;
    }
    informative += 1;
  }
  if (informative < 8 || insideMajor + insideMinor < 3 || outsideMajor + outsideMinor < 3) return 0;

  const majorFit = hmmProbability(
    f64(outsideMajor + 1) / f64(outsideMajor + outsideMinor + 2),
  );
  const minorFit = hmmProbability(
    f64(insideMinor + 1) / f64(insideMajor + insideMinor + 2),
  );
  const transition = Math.max(0.000001, Math.min(0.08, 2.0 / f64(informative)));
  const logStay = Math.log(1.0 - transition);
  const logSwitch = Math.log(transition);

  let position = load<i32>(usize(positionPtr));
  let r = seqBase(seqPtr, nSites, recombinant, position);
  let b = seqBase(seqPtr, nSites, minorParent, position);
  let observationMinor = r == b;
  let scoreMajor = Math.log(0.5) + Math.log(observationMinor ? 1.0 - majorFit : majorFit);
  let scoreMinor = Math.log(0.5) + Math.log(observationMinor ? minorFit : 1.0 - minorFit);
  store<i32>(usize(predecessorPtr), 0);

  for (let index: i32 = 1; index < informative; index += 1) {
    position = load<i32>(usize(positionPtr + index * 4));
    r = seqBase(seqPtr, nSites, recombinant, position);
    b = seqBase(seqPtr, nSites, minorParent, position);
    observationMinor = r == b;

    const majorFromMajor = scoreMajor + logStay;
    const majorFromMinor = scoreMinor + logSwitch;
    const minorFromMajor = scoreMajor + logSwitch;
    const minorFromMinor = scoreMinor + logStay;
    const predecessorMajor: i32 = majorFromMinor > majorFromMajor ? 1 : 0;
    const predecessorMinor: i32 = minorFromMinor >= minorFromMajor ? 1 : 0;
    const nextMajor = (predecessorMajor == 1 ? majorFromMinor : majorFromMajor)
      + Math.log(observationMinor ? 1.0 - majorFit : majorFit);
    const nextMinor = (predecessorMinor == 1 ? minorFromMinor : minorFromMajor)
      + Math.log(observationMinor ? minorFit : 1.0 - minorFit);
    store<i32>(usize(predecessorPtr + index * 4), predecessorMajor | (predecessorMinor << 1));
    scoreMajor = nextMajor;
    scoreMinor = nextMinor;
  }

  let state: i32 = scoreMinor > scoreMajor ? 1 : 0;
  let switches: i32 = 0;
  for (let index = informative - 1; index >= 0; index -= 1) {
    const rawPosition = load<i32>(usize(positionPtr + index * 4));
    store<i32>(usize(positionPtr + index * 4), state == 1 ? -rawPosition - 1 : rawPosition);
    if (index == 0) continue;
    const predecessors = load<i32>(usize(predecessorPtr + index * 4));
    const previous = state == 0 ? predecessors & 1 : (predecessors >> 1) & 1;
    if (previous != state) switches += 1;
    state = previous;
  }

  let bestStart: i32 = start;
  let bestEnd: i32 = end;
  let bestOverlap: i32 = -1;
  let bestRunSites: i32 = 0;
  let bestStartLow: i32 = start;
  let bestStartHigh: i32 = start;
  let bestEndLow: i32 = end;
  let bestEndHigh: i32 = end;
  let runStartIndex: i32 = -1;
  for (let index: i32 = 0; index <= informative; index += 1) {
    const isMinor = index < informative && load<i32>(usize(positionPtr + index * 4)) < 0;
    if (isMinor && runStartIndex < 0) runStartIndex = index;
    if (isMinor || runStartIndex < 0) continue;

    const runEndIndex = index - 1;
    const firstPosition = markedPosition(load<i32>(usize(positionPtr + runStartIndex * 4)));
    const lastPosition = markedPosition(load<i32>(usize(positionPtr + runEndIndex * 4)));
    const previousPosition = runStartIndex > 0
      ? markedPosition(load<i32>(usize(positionPtr + (runStartIndex - 1) * 4)))
      : -1;
    const nextPosition = index < informative
      ? markedPosition(load<i32>(usize(positionPtr + index * 4)))
      : nSites;
    const refinedStart = previousPosition >= 0
      ? (previousPosition + firstPosition + 1) / 2
      : 0;
    const refinedEnd = nextPosition < nSites
      ? (lastPosition + nextPosition + 1) / 2
      : nSites;
    const overlapStart = refinedStart > start ? refinedStart : start;
    const overlapEnd = refinedEnd < end ? refinedEnd : end;
    const candidateOverlap = overlapEnd > overlapStart ? overlapEnd - overlapStart : 0;
    const runSites = runEndIndex - runStartIndex + 1;
    if (candidateOverlap > bestOverlap || (candidateOverlap == bestOverlap && runSites > bestRunSites)) {
      bestOverlap = candidateOverlap;
      bestRunSites = runSites;
      bestStart = refinedStart;
      bestEnd = refinedEnd;
      bestStartLow = previousPosition >= 0 ? previousPosition + 1 : 0;
      bestStartHigh = firstPosition;
      bestEndLow = lastPosition + 1;
      bestEndHigh = nextPosition < nSites ? nextPosition : nSites;
    }
    runStartIndex = -1;
  }

  if (bestOverlap <= 0 || bestRunSites < 3 || bestEnd - bestStart < 4) return 0;
  store<i32>(usize(outPtr), bestStart);
  store<i32>(usize(outPtr + 4), bestEnd);
  store<i32>(usize(outPtr + 8), informative);
  store<i32>(usize(outPtr + 12), switches);
  store<i32>(usize(outPtr + 16), i32(majorFit * 1000.0));
  store<i32>(usize(outPtr + 20), i32(minorFit * 1000.0));
  store<i32>(usize(outPtr + 24), bestRunSites);
  store<i32>(usize(outPtr + 28), bestStartLow);
  store<i32>(usize(outPtr + 32), bestStartHigh);
  store<i32>(usize(outPtr + 36), bestEndLow);
  store<i32>(usize(outPtr + 40), bestEndHigh);
  return 1;
}

// Finds the strongest internal CUSUM excursion for a recombinant candidate and
// two candidate parents, then evaluates the resulting inside/outside 2×2 table.
// Prefix sums make each parent-pair scan O(alignment length).
export function scan_pair(
  seqPtr: i32,
  nSites: i32,
  recombinant: i32,
  parentA: i32,
  parentB: i32,
  minimumSegment: i32,
  prefixAPtr: i32,
  prefixBPtr: i32,
  outPtr: i32,
): i32 {
  store<i32>(usize(prefixAPtr), 0);
  store<i32>(usize(prefixBPtr), 0);
  let countA: i32 = 0;
  let countB: i32 = 0;

  for (let site: i32 = 0; site < nSites; site += 1) {
    const r = seqBase(seqPtr, nSites, recombinant, site);
    const a = seqBase(seqPtr, nSites, parentA, site);
    const b = seqBase(seqPtr, nSites, parentB, site);
    if (valid(r) && valid(a) && valid(b) && a != b) {
      if (r == a) countA += 1;
      else if (r == b) countB += 1;
    }
    store<i32>(usize(prefixAPtr + (site + 1) * 4), countA);
    store<i32>(usize(prefixBPtr + (site + 1) * 4), countB);
  }

  const informative: i32 = countA + countB;
  if (informative < 12 || nSites < minimumSegment * 2) return 0;

  const totalDifference: f64 = f64(countB - countA);
  let minimumValue: f64 = 0.0;
  let minimumPosition: i32 = 0;
  let maximumValue: f64 = 0.0;
  let maximumPosition: i32 = 0;
  let bestPositive: f64 = 0.0;
  let positiveStart: i32 = 0;
  let positiveEnd: i32 = 0;
  let bestNegative: f64 = 0.0;
  let negativeStart: i32 = 0;
  let negativeEnd: i32 = 0;

  for (let position: i32 = 1; position < nSites; position += 1) {
    const prefixA = load<i32>(usize(prefixAPtr + position * 4));
    const prefixB = load<i32>(usize(prefixBPtr + position * 4));
    const detrended = f64(prefixB - prefixA) - totalDifference * f64(position) / f64(nSites);

    if (position - minimumPosition >= minimumSegment) {
      const gain = detrended - minimumValue;
      if (gain > bestPositive) {
        bestPositive = gain;
        positiveStart = minimumPosition;
        positiveEnd = position;
      }
    }
    if (position - maximumPosition >= minimumSegment) {
      const loss = maximumValue - detrended;
      if (loss > bestNegative) {
        bestNegative = loss;
        negativeStart = maximumPosition;
        negativeEnd = position;
      }
    }
    if (detrended < minimumValue) {
      minimumValue = detrended;
      minimumPosition = position;
    }
    if (detrended > maximumValue) {
      maximumValue = detrended;
      maximumPosition = position;
    }
  }

  let start: i32 = bestPositive >= bestNegative ? positiveStart : negativeStart;
  let end: i32 = bestPositive >= bestNegative ? positiveEnd : negativeEnd;
  if (end - start < minimumSegment || end - start > nSites - minimumSegment) return 0;

  const insideA = load<i32>(usize(prefixAPtr + end * 4)) - load<i32>(usize(prefixAPtr + start * 4));
  const insideB = load<i32>(usize(prefixBPtr + end * 4)) - load<i32>(usize(prefixBPtr + start * 4));
  const outsideA = countA - insideA;
  const outsideB = countB - insideB;
  if ((insideA - insideB) * (outsideA - outsideB) >= 0) return 0;

  let majorParent = parentA;
  let minorParent = parentB;
  let insideMinor = insideB;
  let insideMajor = insideA;
  let outsideMajor = outsideA;
  let outsideMinor = outsideB;
  if (insideA > insideB) {
    majorParent = parentB;
    minorParent = parentA;
    insideMinor = insideA;
    insideMajor = insideB;
    outsideMajor = outsideB;
    outsideMinor = outsideA;
  }

  if (insideMinor - insideMajor < 3 || outsideMajor - outsideMinor < 3) return 0;
  const statistic = chiSquare(insideMinor, insideMajor, outsideMinor, outsideMajor);
  if (statistic < 6.0) return 0;

  const insideTotal: f64 = f64(insideMinor + insideMajor);
  const outsideTotal: f64 = f64(outsideMinor + outsideMajor);
  const effect = f64(insideMinor) / insideTotal - f64(outsideMinor) / outsideTotal;

  store<i32>(usize(outPtr), start);
  store<i32>(usize(outPtr + 4), end);
  store<i32>(usize(outPtr + 8), majorParent);
  store<i32>(usize(outPtr + 12), minorParent);
  store<i32>(usize(outPtr + 16), i32(statistic * 1000.0));
  store<i32>(usize(outPtr + 20), informative);
  store<i32>(usize(outPtr + 24), insideMinor);
  store<i32>(usize(outPtr + 28), insideMajor);
  store<i32>(usize(outPtr + 32), outsideMajor);
  store<i32>(usize(outPtr + 36), outsideMinor);
  store<i32>(usize(outPtr + 40), i32(effect * 1000000.0));
  store<i32>(usize(outPtr + 44), bestPositive >= bestNegative ? 1 : -1);
  return 1;
}

@inline
function packedTripletCategory(scratchPtr: i32, index: i32): i32 {
  return load<i32>(usize(scratchPtr + index * 4)) & 3;
}

@inline
function packedTripletPosition(scratchPtr: i32, index: i32): i32 {
  return load<i32>(usize(scratchPtr + index * 4)) >> 2;
}

@inline
function pairContains(pair: i32, sequence: i32): bool {
  if (pair == 0) return sequence == 0 || sequence == 1;
  if (pair == 1) return sequence == 0 || sequence == 2;
  return sequence == 1 || sequence == 2;
}

@inline
function sharedPairMember(left: i32, right: i32): i32 {
  for (let sequence: i32 = 0; sequence < 3; sequence += 1) {
    if (pairContains(left, sequence) && pairContains(right, sequence)) return sequence;
  }
  return -1;
}

@inline
function otherPairMember(pair: i32, member: i32): i32 {
  if (pair == 0) return member == 0 ? 1 : 0;
  if (pair == 1) return member == 0 ? 2 : 0;
  return member == 1 ? 2 : 1;
}

@inline
function actualTripletSequence(local: i32, seq1: i32, seq2: i32, seq3: i32): i32 {
  return local == 0 ? seq1 : local == 1 ? seq2 : seq3;
}

@inline
function circularInformativeIndex(index: i32, length: i32): i32 {
  let result = index % length;
  if (result < 0) result += length;
  return result;
}

// Source-compatible RDP triplet detector.  This ports the active
// FindSubSeqPB3 -> XOHomologyP2 -> FindNextP -> DefineEventP2 path in the
// supplied RDP5 DNA5 source.  Only variable positions at which exactly two
// members of the triplet agree are retained.  The three pair-match streams
// are ranked globally, scanned in a 2h+1 VNP circular window, and a local
// medium-pair dominance excursion is delineated with the RDP5 common/different
// rule. The compatibility wrapper returns the best event; the multi-signal
// export retains every distinct excursion up to a caller-selected capacity.
//
// scratchPackedPtr: nSites i32 values (position << 2 | pair category)
// scratchDominancePtr: nSites i32 values
// outPtr layout (i32):
//  0..1 genomic start/end; 2..4 recombinant/major/minor sequence indexes;
//  5 chi-square*1000; 6 informative VNPs; 7..10 inside-minor, inside-major,
//  outside-major, outside-minor; 11 effect*1e6; 12 source common sites;
//  13 source tract VNPs; 14 source medium-pair sites; 15 orientation cycle;
//  16 source full-window size; 17 compatibility marker.
function writeRdpSignal(
  outPtr: i32,
  start: i32,
  end: i32,
  daughter: i32,
  major: i32,
  minor: i32,
  statistic: f64,
  informative: i32,
  insideMinor: i32,
  insideMajor: i32,
  outsideMajor: i32,
  outsideMinor: i32,
  common: i32,
  tractLength: i32,
  mediumCount: i32,
  cycle: i32,
  sourceWindow: i32,
): void {
  const insideTotal = insideMinor + insideMajor;
  const outsideTotal = outsideMinor + outsideMajor;
  const effect = insideTotal > 0 && outsideTotal > 0
    ? f64(insideMinor) / f64(insideTotal) - f64(outsideMinor) / f64(outsideTotal)
    : 0.0;
  store<i32>(usize(outPtr), start);
  store<i32>(usize(outPtr + 4), end);
  store<i32>(usize(outPtr + 8), daughter);
  store<i32>(usize(outPtr + 12), major);
  store<i32>(usize(outPtr + 16), minor);
  store<i32>(usize(outPtr + 20), i32(statistic * 1000.0));
  store<i32>(usize(outPtr + 24), informative);
  store<i32>(usize(outPtr + 28), insideMinor);
  store<i32>(usize(outPtr + 32), insideMajor);
  store<i32>(usize(outPtr + 36), outsideMajor);
  store<i32>(usize(outPtr + 40), outsideMinor);
  store<i32>(usize(outPtr + 44), i32(effect * 1000000.0));
  store<i32>(usize(outPtr + 48), common);
  store<i32>(usize(outPtr + 52), tractLength);
  store<i32>(usize(outPtr + 56), mediumCount);
  store<i32>(usize(outPtr + 60), cycle);
  store<i32>(usize(outPtr + 64), sourceWindow);
  store<i32>(usize(outPtr + 68), 1);
}

function scanRdp5TripletCore(
  seqPtr: i32,
  nSites: i32,
  seq1: i32,
  seq2: i32,
  seq3: i32,
  fullWindow: i32,
  scratchPackedPtr: i32,
  scratchDominancePtr: i32,
  bestOutPtr: i32,
  allOutPtr: i32,
  maximumEvents: i32,
): i32 {
  if (nSites < 4) return 0;
  const sourceWindow = fullWindow < 5 ? 5 : fullWindow;
  const halfWindow = sourceWindow / 2;
  const windowWidth = halfWindow * 2 + 1;
  let informative: i32 = 0;
  let count0: i32 = 0;
  let count1: i32 = 0;
  let count2: i32 = 0;

  for (let site: i32 = 0; site < nSites; site += 1) {
    const a = seqBase(seqPtr, nSites, seq1, site);
    const b = seqBase(seqPtr, nSites, seq2, site);
    const c = seqBase(seqPtr, nSites, seq3, site);
    if (!valid(a) || !valid(b) || !valid(c)) continue;
    let category: i32 = -1;
    if (a == b && a != c) category = 0;
    else if (a == c && a != b) category = 1;
    else if (b == c && b != a) category = 2;
    if (category < 0) continue;
    store<i32>(usize(scratchPackedPtr + informative * 4), (site << 2) | category);
    informative += 1;
    if (category == 0) count0 += 1;
    else if (category == 1) count1 += 1;
    else count2 += 1;
  }

  // FastRecCheckP uses XoverWindow = XOverWindowX / 2 and requires at least
  // twice that many retained sites plus one third-window observation in every
  // pair category.
  if (informative < halfWindow * 2 || informative < windowWidth) return 0;
  const minimumCategory = halfWindow / 3;
  if (count0 < minimumCategory || count1 < minimumCategory || count2 < minimumCategory) return 0;

  let high: i32 = 0;
  if (count1 > count0) high = 1;
  if ((high == 0 ? count0 : count1) < count2) high = 2;
  let medium: i32 = high == 0 ? 1 : 0;
  let low: i32 = high == 2 ? 1 : 2;
  const mediumCountInitial = medium == 0 ? count0 : medium == 1 ? count1 : count2;
  const lowCountInitial = low == 0 ? count0 : low == 1 ? count1 : count2;
  if (lowCountInitial > mediumCountInitial) {
    const swap = medium;
    medium = low;
    low = swap;
  }
  const originalHigh = high;
  const originalMedium = medium;
  const originalLow = low;
  const originalHighCount = high == 0 ? count0 : high == 1 ? count1 : count2;

  let bestStatistic: f64 = 0.0;
  let bestStart: i32 = 0;
  let bestEnd: i32 = 0;
  let bestDaughter: i32 = -1;
  let bestMajor: i32 = -1;
  let bestMinor: i32 = -1;
  let bestInsideMinor: i32 = 0;
  let bestInsideMajor: i32 = 0;
  let bestOutsideMajor: i32 = 0;
  let bestOutsideMinor: i32 = 0;
  let bestCommon: i32 = 0;
  let bestTractLength: i32 = 0;
  let bestMediumCount: i32 = 0;
  let bestCycle: i32 = 0;
  let retainedEvents: i32 = 0;
  let totalEvents: i32 = 0;

  for (let cycle: i32 = 0; cycle < 3; cycle += 1) {
    if (cycle == 0) {
      high = originalHigh;
      medium = originalMedium;
      low = originalLow;
    } else if (cycle == 1) {
      high = originalHigh;
      medium = originalLow;
      low = originalMedium;
    } else {
      // This is the third orientation used by FastRecCheckP after the first
      // medium/low swap.  RDP5 omits it for very high global homology.
      if (f64(originalHighCount) / f64(informative) >= 0.7) break;
      high = originalMedium;
      medium = originalHigh;
      low = originalLow;
    }

    const daughterLocal = sharedPairMember(high, medium);
    if (daughterLocal < 0) continue;
    const majorLocal = otherPairMember(high, daughterLocal);
    const minorLocal = otherPairMember(medium, daughterLocal);

    let local0: i32 = 0;
    let local1: i32 = 0;
    let local2: i32 = 0;
    for (let offset: i32 = -halfWindow; offset <= halfWindow; offset += 1) {
      const category = packedTripletCategory(
        scratchPackedPtr,
        circularInformativeIndex(offset, informative),
      );
      if (category == 0) local0 += 1;
      else if (category == 1) local1 += 1;
      else local2 += 1;
    }
    for (let center: i32 = 0; center < informative; center += 1) {
      if (center > 0) {
        const leaving = packedTripletCategory(
          scratchPackedPtr,
          circularInformativeIndex(center - halfWindow - 1, informative),
        );
        const entering = packedTripletCategory(
          scratchPackedPtr,
          circularInformativeIndex(center + halfWindow, informative),
        );
        if (leaving == 0) local0 -= 1;
        else if (leaving == 1) local1 -= 1;
        else local2 -= 1;
        if (entering == 0) local0 += 1;
        else if (entering == 1) local1 += 1;
        else local2 += 1;
      }
      const highLocal = high == 0 ? local0 : high == 1 ? local1 : local2;
      const mediumLocal = medium == 0 ? local0 : medium == 1 ? local1 : local2;
      const lowLocal = low == 0 ? local0 : low == 1 ? local1 : local2;
      const strictDominance = mediumLocal > highLocal && mediumLocal > lowLocal;
      const continuationDominance = mediumLocal >= highLocal && mediumLocal >= lowLocal;
      store<i32>(usize(scratchDominancePtr + center * 4), strictDominance ? 2 : continuationDominance ? 1 : 0);
    }

    let seed: i32 = 0;
    while (seed < informative) {
      while (seed < informative && load<i32>(usize(scratchDominancePtr + seed * 4)) != 2) seed += 1;
      if (seed >= informative) break;
      let dominanceEnd = seed;
      while (dominanceEnd + 1 < informative
        && load<i32>(usize(scratchDominancePtr + (dominanceEnd + 1) * 4)) == 2) dominanceEnd += 1;

      let eventStart = seed;
      while (eventStart > 0 && packedTripletCategory(scratchPackedPtr, eventStart - 1) == medium) eventStart -= 1;
      while (eventStart < informative && packedTripletCategory(scratchPackedPtr, eventStart) != medium) eventStart += 1;
      if (eventStart >= informative) break;

      let common: i32 = 0;
      let tractSites: i32 = 0;
      let lastCommon: i32 = eventStart;
      let cursor: i32 = eventStart;
      while (cursor < informative) {
        const category = packedTripletCategory(scratchPackedPtr, cursor);
        tractSites += 1;
        if (category == medium) {
          common += 1;
          lastCommon = cursor;
        }
        const next = cursor + 1;
        if (next >= informative) break;
        // DefineEventP2 continues through tied window counts and only stops
        // when medium is strictly below either competitor at a non-common VNP.
        const nextDominant = load<i32>(usize(scratchDominancePtr + next * 4)) != 0;
        const nextCommon = packedTripletCategory(scratchPackedPtr, next) == medium;
        if (!nextDominant && !nextCommon) break;
        cursor = next;
      }
      tractSites -= cursor - lastCommon;
      const different = tractSites - common;
      if (tractSites > 2 && common > f64(different) * 0.8) {
        let insideMajor: i32 = 0;
        let insideMinor: i32 = 0;
        for (let index = eventStart; index <= lastCommon; index += 1) {
          const category = packedTripletCategory(scratchPackedPtr, index);
          if (category == high) insideMajor += 1;
          else if (category == medium) insideMinor += 1;
        }
        const highCount = high == 0 ? count0 : high == 1 ? count1 : count2;
        const mediumCount = medium == 0 ? count0 : medium == 1 ? count1 : count2;
        const outsideMajor = highCount - insideMajor;
        const outsideMinor = mediumCount - insideMinor;
        const statistic = chiSquare(insideMinor, insideMajor, outsideMinor, outsideMajor);
        const signalStart = packedTripletPosition(scratchPackedPtr, eventStart);
        const signalEnd = packedTripletPosition(scratchPackedPtr, lastCommon) + 1;
        const signalDaughter = actualTripletSequence(daughterLocal, seq1, seq2, seq3);
        const signalMajor = actualTripletSequence(majorLocal, seq1, seq2, seq3);
        const signalMinor = actualTripletSequence(minorLocal, seq1, seq2, seq3);
        if (maximumEvents > 0 && statistic > 0.0 && signalEnd > signalStart) {
          let signalSlot: i32 = -1;
          let duplicate: i32 = -1;
          for (let retained: i32 = 0; retained < retainedEvents; retained += 1) {
            const retainedPtr = allOutPtr + retained * 72;
            if (load<i32>(usize(retainedPtr)) == signalStart
              && load<i32>(usize(retainedPtr + 4)) == signalEnd
              && load<i32>(usize(retainedPtr + 8)) == signalDaughter
              && load<i32>(usize(retainedPtr + 12)) == signalMajor
              && load<i32>(usize(retainedPtr + 16)) == signalMinor) {
              duplicate = retained;
              break;
            }
          }
          if (duplicate >= 0) {
            if (i32(statistic * 1000.0) > load<i32>(usize(allOutPtr + duplicate * 72 + 20))) signalSlot = duplicate;
          } else {
            totalEvents += 1;
            if (retainedEvents < maximumEvents) {
              signalSlot = retainedEvents;
              retainedEvents += 1;
            } else {
              let weakestSlot: i32 = 0;
              let weakestStatistic = load<i32>(usize(allOutPtr + 20));
              for (let retained: i32 = 1; retained < retainedEvents; retained += 1) {
                const retainedStatistic = load<i32>(usize(allOutPtr + retained * 72 + 20));
                if (retainedStatistic < weakestStatistic) {
                  weakestStatistic = retainedStatistic;
                  weakestSlot = retained;
                }
              }
              if (i32(statistic * 1000.0) > weakestStatistic) signalSlot = weakestSlot;
            }
          }
          if (signalSlot >= 0) writeRdpSignal(
            allOutPtr + signalSlot * 72,
            signalStart,
            signalEnd,
            signalDaughter,
            signalMajor,
            signalMinor,
            statistic,
            informative,
            insideMinor,
            insideMajor,
            outsideMajor,
            outsideMinor,
            common,
            tractSites,
            mediumCount,
            cycle,
            sourceWindow,
          );
        }
        if (statistic > bestStatistic) {
          bestStatistic = statistic;
          bestStart = signalStart;
          bestEnd = signalEnd;
          bestDaughter = signalDaughter;
          bestMajor = signalMajor;
          bestMinor = signalMinor;
          bestInsideMinor = insideMinor;
          bestInsideMajor = insideMajor;
          bestOutsideMajor = outsideMajor;
          bestOutsideMinor = outsideMinor;
          bestCommon = common;
          bestTractLength = tractSites;
          bestMediumCount = mediumCount;
          bestCycle = cycle;
        }
      }
      seed = dominanceEnd + 1;
    }
  }

  if (bestDaughter < 0 || bestEnd <= bestStart || bestStatistic <= 0.0) return 0;
  writeRdpSignal(
    bestOutPtr,
    bestStart,
    bestEnd,
    bestDaughter,
    bestMajor,
    bestMinor,
    bestStatistic,
    informative,
    bestInsideMinor,
    bestInsideMajor,
    bestOutsideMajor,
    bestOutsideMinor,
    bestCommon,
    bestTractLength,
    bestMediumCount,
    bestCycle,
    sourceWindow,
  );
  return maximumEvents > 0 ? totalEvents : 1;
}

export function scan_rdp5_triplet(
  seqPtr: i32,
  nSites: i32,
  seq1: i32,
  seq2: i32,
  seq3: i32,
  fullWindow: i32,
  scratchPackedPtr: i32,
  scratchDominancePtr: i32,
  outPtr: i32,
): i32 {
  return scanRdp5TripletCore(seqPtr, nSites, seq1, seq2, seq3, fullWindow, scratchPackedPtr, scratchDominancePtr, outPtr, 0, 0);
}

export function scan_rdp5_triplet_all(
  seqPtr: i32,
  nSites: i32,
  seq1: i32,
  seq2: i32,
  seq3: i32,
  fullWindow: i32,
  scratchPackedPtr: i32,
  scratchDominancePtr: i32,
  outPtr: i32,
  maximumEvents: i32,
  bestOutPtr: i32,
): i32 {
  return scanRdp5TripletCore(seqPtr, nSites, seq1, seq2, seq3, fullWindow, scratchPackedPtr, scratchDominancePtr, bestOutPtr, outPtr, maximumEvents);
}
