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

// RDP5 FastRecCheckMC/FastRecCheckChim share a considerably richer peak
// locator than the single strongest-pair shortcut used by early RDP Web
// builds.  The source path below keeps the historical compressed-coordinate
// scan, strict first-maximum ordering, 11-position smoothing basin,
// GrowMChiWin expansion and three triplet tracks in one allocation-free WASM
// kernel.  Scratch buffers are owned by the worker and reused across every
// concrete triplet.

const SOURCE_CHI_ROW_INTS: i32 = 16;
const SOURCE_CHI_PEAK_INTS: i32 = 6;

@inline
function sourceChiScore(scoresPtr: i32, rank: i32, scoreBit: i32): i32 {
  return (i32(load<u8>(usize(scoresPtr + rank))) >> scoreBit) & 1;
}

@inline
function sourceChiScaled(statistic: f64): i32 {
  if (!(statistic > 0.0)) return 0;
  const bounded = statistic > 2_000_000.0 ? 2_000_000.0 : statistic;
  return i32(bounded * 1000.0);
}

@inline
function sourceChiModulo(value: i32, length: i32): i32 {
  let result = value % length;
  if (result < 0) result += length;
  return result;
}

function sourceChiHalfWindow(informative: i32, fullWindow: i32): i32 {
  if (informative < 7) return 0;
  let halfWindow = (fullWindow + 1) / 2;
  const criticalDifference: i32 = 2;
  if (halfWindow * 2 > informative) {
    halfWindow = i32((f64(informative) * 0.75 / 2.0) + 0.51) - 1;
  }
  if (halfWindow <= criticalDifference) {
    halfWindow = i32((f64(informative) / 2.0) + 0.51) - 1;
  }
  return halfWindow >= 6 ? halfWindow : 0;
}

// Builds the shared MAXCHI track and all three CHIMAERA target tracks in one
// nucleotide pass. The desktop routines reuse triplet encodings between
// methods; materializing four compact lanes here gives the browser the same
// triplet-level economy while keeping each downstream source scan unchanged.
function sourceChiBuildAllTracks(
  seqPtr: i32,
  nSites: i32,
  sequence0: i32,
  sequence1: i32,
  sequence2: i32,
  positionsPtr: i32,
  scoresPtr: i32,
  missingPrefixPtr: i32,
  countsPtr: i32,
  methodMask: i32,
): void {
  const positionStride = (nSites + 1) * 4;
  const scoreStride = nSites + 1;
  let maxChiCount: i32 = 0;
  let chimaera0Count: i32 = 0;
  let chimaera1Count: i32 = 0;
  let chimaera2Count: i32 = 0;
  let missing: i32 = 0;
  store<i32>(usize(missingPrefixPtr), 0);
  for (let site: i32 = 0; site < nSites; site += 1) {
    const a = seqBase(seqPtr, nSites, sequence0, site);
    const b = seqBase(seqPtr, nSites, sequence1, site);
    const c = seqBase(seqPtr, nSites, sequence2, site);
    const complete = valid(a) && valid(b) && valid(c);
    if (!complete) {
      missing += 1;
    } else {
      if ((methodMask & 4) != 0 && !(a == b && a == c)) {
        store<i32>(usize(positionsPtr + maxChiCount * 4), site);
        const scoreBits = i32(a == b) | (i32(a == c) << 1) | (i32(b == c) << 2);
        store<u8>(usize(scoresPtr + maxChiCount), u8(scoreBits));
        maxChiCount += 1;
      }
      if ((methodMask & 8) != 0 && b != c && (a == b || a == c)) {
        store<i32>(usize(positionsPtr + positionStride + chimaera0Count * 4), site);
        store<u8>(usize(scoresPtr + scoreStride + chimaera0Count), u8(a == b));
        chimaera0Count += 1;
      }
      if ((methodMask & 8) != 0 && c != a && (b == c || b == a)) {
        store<i32>(usize(positionsPtr + positionStride * 2 + chimaera1Count * 4), site);
        store<u8>(usize(scoresPtr + scoreStride * 2 + chimaera1Count), u8(b == c));
        chimaera1Count += 1;
      }
      if ((methodMask & 8) != 0 && a != b && (c == a || c == b)) {
        store<i32>(usize(positionsPtr + positionStride * 3 + chimaera2Count * 4), site);
        store<u8>(usize(scoresPtr + scoreStride * 3 + chimaera2Count), u8(c == a));
        chimaera2Count += 1;
      }
    }
    store<i32>(usize(missingPrefixPtr + (site + 1) * 4), missing);
  }
  store<i32>(usize(countsPtr), maxChiCount);
  store<i32>(usize(countsPtr + 4), chimaera0Count);
  store<i32>(usize(countsPtr + 8), chimaera1Count);
  store<i32>(usize(countsPtr + 12), chimaera2Count);
  // A negative first prefix entry is the no-missing-data fast-path marker.
  // It lets every downstream window test avoid reading the O(L) prefix lane.
  if (missing == 0) store<i32>(usize(missingPrefixPtr), -1);
}

// The RDP5 desktop first compresses the alignment and then decodes triplet
// states through lookup tables instead of rereading three character strings
// for every triplet. The browser already owns an exact two-bit alignment plus
// one validity bit per nucleotide lane for distance calculations. Reusing it
// here advances sixteen raw sites per word and only materializes informative
// compressed coordinates. This is algebraically identical to
// sourceChiBuildAllTracks; the byte path remains as a validation oracle.
function sourceChiBuildAllTracksPacked(
  packedPtr: i32,
  validityPtr: i32,
  wordsPerSequence: i32,
  nSites: i32,
  sequence0: i32,
  sequence1: i32,
  sequence2: i32,
  positionsPtr: i32,
  scoresPtr: i32,
  missingPrefixPtr: i32,
  countsPtr: i32,
  methodMask: i32,
): void {
  const positionStride = (nSites + 1) * 4;
  const scoreStride = nSites + 1;
  const offset0 = sequence0 * wordsPerSequence;
  const offset1 = sequence1 * wordsPerSequence;
  const offset2 = sequence2 * wordsPerSequence;
  let maxChiCount: i32 = 0;
  let chimaera0Count: i32 = 0;
  let chimaera1Count: i32 = 0;
  let chimaera2Count: i32 = 0;
  let missing: i32 = 0;

  for (let word: i32 = 0; word < wordsPerSequence; word += 1) {
    const a = load<u32>(usize(packedPtr + (offset0 + word) * 4));
    const b = load<u32>(usize(packedPtr + (offset1 + word) * 4));
    const c = load<u32>(usize(packedPtr + (offset2 + word) * 4));
    let rawLaneMask: u32 = 0x55555555;
    const remaining = nSites - word * 16;
    if (remaining < 16) {
      rawLaneMask = 0;
      for (let lane: i32 = 0; lane < remaining; lane += 1) rawLaneMask |= u32(1) << u32(lane * 2);
    }
    const allValid = load<u32>(usize(validityPtr + (offset0 + word) * 4))
      & load<u32>(usize(validityPtr + (offset1 + word) * 4))
      & load<u32>(usize(validityPtr + (offset2 + word) * 4))
      & rawLaneMask;
    missing += popCount32(rawLaneMask & ~allValid);
    const ab = packedEquality(a, b, allValid);
    const ac = packedEquality(a, c, allValid);
    const bc = packedEquality(b, c, allValid);

    if ((methodMask & 4) != 0) {
      let informativeMask = allValid & ~(ab & ac);
      while (informativeMask != 0) {
        const bit = i32(ctz<u32>(informativeMask));
        const lane = u32(1) << u32(bit);
        const site = word * 16 + bit / 2;
        store<i32>(usize(positionsPtr + maxChiCount * 4), site);
        const scoreBits = i32((ab & lane) != 0)
          | (i32((ac & lane) != 0) << 1)
          | (i32((bc & lane) != 0) << 2);
        store<u8>(usize(scoresPtr + maxChiCount), u8(scoreBits));
        maxChiCount += 1;
        informativeMask &= informativeMask - 1;
      }
    }

    if ((methodMask & 8) != 0) {
      let chimaera0 = allValid & ~bc & (ab | ac);
      let chimaera1 = allValid & ~ac & (bc | ab);
      let chimaera2 = allValid & ~ab & (ac | bc);
      while (chimaera0 != 0) {
        const bit = i32(ctz<u32>(chimaera0));
        const lane = u32(1) << u32(bit);
        const site = word * 16 + bit / 2;
        store<i32>(usize(positionsPtr + positionStride + chimaera0Count * 4), site);
        store<u8>(usize(scoresPtr + scoreStride + chimaera0Count), u8((ab & lane) != 0));
        chimaera0Count += 1;
        chimaera0 &= chimaera0 - 1;
      }
      while (chimaera1 != 0) {
        const bit = i32(ctz<u32>(chimaera1));
        const lane = u32(1) << u32(bit);
        const site = word * 16 + bit / 2;
        store<i32>(usize(positionsPtr + positionStride * 2 + chimaera1Count * 4), site);
        store<u8>(usize(scoresPtr + scoreStride * 2 + chimaera1Count), u8((bc & lane) != 0));
        chimaera1Count += 1;
        chimaera1 &= chimaera1 - 1;
      }
      while (chimaera2 != 0) {
        const bit = i32(ctz<u32>(chimaera2));
        const lane = u32(1) << u32(bit);
        const site = word * 16 + bit / 2;
        store<i32>(usize(positionsPtr + positionStride * 3 + chimaera2Count * 4), site);
        store<u8>(usize(scoresPtr + scoreStride * 3 + chimaera2Count), u8((ac & lane) != 0));
        chimaera2Count += 1;
        chimaera2 &= chimaera2 - 1;
      }
    }
  }

  if (missing == 0) {
    store<i32>(usize(missingPrefixPtr), -1);
  } else {
    let missingPrefix: i32 = 0;
    store<i32>(usize(missingPrefixPtr), 0);
    for (let site: i32 = 0; site < nSites; site += 1) {
      const word = site >> 4;
      const lane = u32(1) << u32((site & 15) * 2);
      const complete = (load<u32>(usize(validityPtr + (offset0 + word) * 4))
        & load<u32>(usize(validityPtr + (offset1 + word) * 4))
        & load<u32>(usize(validityPtr + (offset2 + word) * 4))
        & lane) != 0;
      if (!complete) missingPrefix += 1;
      store<i32>(usize(missingPrefixPtr + (site + 1) * 4), missingPrefix);
    }
  }
  store<i32>(usize(countsPtr), maxChiCount);
  store<i32>(usize(countsPtr + 4), chimaera0Count);
  store<i32>(usize(countsPtr + 8), chimaera1Count);
  store<i32>(usize(countsPtr + 12), chimaera2Count);
}

@inline
function sourceChiRawSpanHasMissing(
  positionsPtr: i32,
  missingPrefixPtr: i32,
  firstRank: i32,
  lastRank: i32,
): bool {
  if (prefixValue(missingPrefixPtr, 0) < 0) return false;
  const first = load<i32>(usize(positionsPtr + firstRank * 4));
  const last = load<i32>(usize(positionsPtr + lastRank * 4));
  return prefixValue(missingPrefixPtr, last + 1) - prefixValue(missingPrefixPtr, first) > 0;
}

function sourceChiFillProfile(
  informative: i32,
  halfWindow: i32,
  circular: bool,
  scoreBit: i32,
  positionsPtr: i32,
  scoresPtr: i32,
  missingPrefixPtr: i32,
  chiPtr: i32,
  smoothPtr: i32,
): void {
  for (let rank: i32 = 0; rank < informative; rank += 1) {
    store<f64>(usize(chiPtr + rank * 8), 0.0);
    store<f64>(usize(smoothPtr + rank * 8), 0.0);
  }
  if (halfWindow <= 0 || informative < halfWindow * 2) return;

  let leftOnes: i32 = 0;
  let rightOnes: i32 = 0;
  if (circular) {
    for (let offset: i32 = -halfWindow; offset < 0; offset += 1) {
      leftOnes += sourceChiScore(scoresPtr, sourceChiModulo(offset, informative), scoreBit);
    }
    for (let offset: i32 = 0; offset < halfWindow; offset += 1) {
      rightOnes += sourceChiScore(scoresPtr, offset, scoreBit);
    }
    for (let rank: i32 = 0; rank < informative; rank += 1) {
      const statistic = chiSquare(
        leftOnes,
        halfWindow - leftOnes,
        rightOnes,
        halfWindow - rightOnes,
      );
      const difference = leftOnes - rightOnes;
      store<f64>(usize(chiPtr + rank * 8), Math.abs(difference) > 2 ? statistic : 0.0);
      const leavingLeft = sourceChiModulo(rank - halfWindow, informative);
      const enteringLeft = rank;
      const leavingRight = rank;
      const enteringRight = sourceChiModulo(rank + halfWindow, informative);
      leftOnes += sourceChiScore(scoresPtr, enteringLeft, scoreBit)
        - sourceChiScore(scoresPtr, leavingLeft, scoreBit);
      rightOnes += sourceChiScore(scoresPtr, enteringRight, scoreBit)
        - sourceChiScore(scoresPtr, leavingRight, scoreBit);
    }
  } else {
    for (let rank: i32 = 0; rank < halfWindow; rank += 1) {
      leftOnes += sourceChiScore(scoresPtr, rank, scoreBit);
      rightOnes += sourceChiScore(scoresPtr, rank + halfWindow, scoreBit);
    }
    for (let boundary: i32 = halfWindow; boundary + halfWindow <= informative; boundary += 1) {
      const missing = sourceChiRawSpanHasMissing(
        positionsPtr,
        missingPrefixPtr,
        boundary - halfWindow,
        boundary + halfWindow - 1,
      );
      const difference = leftOnes - rightOnes;
      if (!missing && Math.abs(difference) > 2) {
        store<f64>(usize(chiPtr + boundary * 8), chiSquare(
          leftOnes,
          halfWindow - leftOnes,
          rightOnes,
          halfWindow - rightOnes,
        ));
      }
      if (boundary + halfWindow < informative) {
        leftOnes += sourceChiScore(scoresPtr, boundary, scoreBit)
          - sourceChiScore(scoresPtr, boundary - halfWindow, scoreBit);
        rightOnes += sourceChiScore(scoresPtr, boundary + halfWindow, scoreBit)
          - sourceChiScore(scoresPtr, boundary, scoreBit);
      }
    }
  }

  // SmoothChiValsP/SmoothChiVals3P use a circular 11-position mean even for
  // linear scans; invalid end windows are already zero above.
  let smoothSum: f64 = 0.0;
  for (let offset: i32 = -5; offset <= 5; offset += 1) {
    smoothSum += load<f64>(usize(chiPtr + sourceChiModulo(offset, informative) * 8));
  }
  for (let rank: i32 = 0; rank < informative; rank += 1) {
    store<f64>(usize(smoothPtr + rank * 8), smoothSum / 11.0);
    const leaving = sourceChiModulo(rank - 5, informative);
    const entering = sourceChiModulo(rank + 6, informative);
    smoothSum += load<f64>(usize(chiPtr + entering * 8))
      - load<f64>(usize(chiPtr + leaving * 8));
  }
}

function sourceChiCount(
  scoresPtr: i32,
  informative: i32,
  start: i32,
  count: i32,
  scoreBit: i32,
  circular: bool,
): i32 {
  let total: i32 = 0;
  for (let offset: i32 = 0; offset < count; offset += 1) {
    let rank = start + offset;
    if (circular) rank = sourceChiModulo(rank, informative);
    else if (rank < 0 || rank >= informative) continue;
    total += sourceChiScore(scoresPtr, rank, scoreBit);
  }
  return total;
}

function sourceChiPeakSign(
  scoresPtr: i32,
  informative: i32,
  peak: i32,
  halfWindow: i32,
  scoreBit: i32,
  circular: bool,
): i32 {
  const left = sourceChiCount(scoresPtr, informative, peak - halfWindow, halfWindow, scoreBit, circular);
  const right = sourceChiCount(scoresPtr, informative, peak, halfWindow, scoreBit, circular);
  return left >= right ? 1 : -1;
}

function sourceChiGrowPeak(
  informative: i32,
  peak: i32,
  halfWindow: i32,
  scoreBit: i32,
  circular: bool,
  positionsPtr: i32,
  scoresPtr: i32,
  missingPrefixPtr: i32,
  initialStatistic: f64,
  resultPtr: i32,
): void {
  let tWindow = i32((f64(halfWindow) / 4.0) + 0.51);
  if (tWindow < 6) tWindow = 6;
  if (tWindow > halfWindow) tWindow = halfWindow;
  if (tWindow > informative / 2) tWindow = informative / 2;
  let leftOnes = sourceChiCount(scoresPtr, informative, peak - tWindow, tWindow, scoreBit, circular);
  let rightOnes = sourceChiCount(scoresPtr, informative, peak, tWindow, scoreBit, circular);
  let bestStatistic = initialStatistic;
  let bestWidth = halfWindow;
  let failCount: i32 = 0;
  let maximumFails = halfWindow * 2;
  const availableFails = (informative - tWindow * 2) / 2;
  if (maximumFails > availableFails) maximumFails = availableFails;
  if (maximumFails <= 0) maximumFails = 1;

  tWindow += 1;
  while (failCount <= maximumFails && tWindow <= informative / 2) {
    const leftRank = peak - tWindow;
    const rightRank = peak + tWindow - 1;
    if (!circular) {
      if (leftRank < 0 || rightRank >= informative) break;
      if (sourceChiRawSpanHasMissing(positionsPtr, missingPrefixPtr, leftRank, rightRank)) break;
    }
    leftOnes += sourceChiScore(scoresPtr, circular ? sourceChiModulo(leftRank, informative) : leftRank, scoreBit);
    rightOnes += sourceChiScore(scoresPtr, circular ? sourceChiModulo(rightRank, informative) : rightRank, scoreBit);
    const statistic = chiSquare(
      leftOnes,
      tWindow - leftOnes,
      rightOnes,
      tWindow - rightOnes,
    );
    if (statistic >= bestStatistic) {
      bestStatistic = statistic;
      bestWidth = tWindow;
      failCount = 0;
    } else {
      failCount += 1;
    }
    tWindow += 1;
  }
  store<i32>(usize(resultPtr), sourceChiScaled(bestStatistic));
  store<i32>(usize(resultPtr + 4), bestWidth);
}

function sourceChiSuppressPeak(
  informative: i32,
  peak: i32,
  chiPtr: i32,
  smoothPtr: i32,
): void {
  let right = peak;
  let rightSteps: i32 = 0;
  while (rightSteps < informative) {
    const current = load<f64>(usize(smoothPtr + right * 8));
    const next1 = sourceChiModulo(right + 1, informative);
    const next2 = sourceChiModulo(right + 2, informative);
    if (!(current > 0.0 && (current >= load<f64>(usize(smoothPtr + next1 * 8))
      || current >= load<f64>(usize(smoothPtr + next2 * 8))))) break;
    right = next1;
    rightSteps += 1;
  }
  let left = peak;
  let leftSteps: i32 = 0;
  while (leftSteps < informative) {
    const current = load<f64>(usize(smoothPtr + left * 8));
    const previous1 = sourceChiModulo(left - 1, informative);
    const previous2 = sourceChiModulo(left - 2, informative);
    if (!(current > 0.0 && (current >= load<f64>(usize(smoothPtr + previous1 * 8))
      || current >= load<f64>(usize(smoothPtr + previous2 * 8))))) break;
    left = previous1;
    leftSteps += 1;
  }
  if (leftSteps + rightSteps >= informative - 1) {
    for (let rank: i32 = 0; rank < informative; rank += 1) {
      store<f64>(usize(chiPtr + rank * 8), 0.0);
    }
    return;
  }
  let rank = left;
  for (let step: i32 = 0; step < informative; step += 1) {
    store<f64>(usize(chiPtr + rank * 8), 0.0);
    if (rank == right) break;
    rank = sourceChiModulo(rank + 1, informative);
  }
}

@inline
function sourceChiPeakField(peakPtr: i32, peak: i32, field: i32): i32 {
  return load<i32>(usize(peakPtr + (peak * SOURCE_CHI_PEAK_INTS + field) * 4));
}

@inline
function sourceChiStorePeakField(peakPtr: i32, peak: i32, field: i32, value: i32): void {
  store<i32>(usize(peakPtr + (peak * SOURCE_CHI_PEAK_INTS + field) * 4), value);
}

function sourceChiBoundaryPosition(positionsPtr: i32, informative: i32, rank: i32, nSites: i32): i32 {
  if (rank <= 0) return 0;
  if (rank >= informative) return nSites;
  return load<i32>(usize(positionsPtr + rank * 4));
}

function sourceChiCollectPeaks(
  nSites: i32,
  informative: i32,
  halfWindow: i32,
  circular: bool,
  scoreBit: i32,
  positionsPtr: i32,
  scoresPtr: i32,
  missingPrefixPtr: i32,
  chiPtr: i32,
  smoothPtr: i32,
  peakPtr: i32,
  peakCapacity: i32,
): i32 {
  let peakCount: i32 = 0;
  while (peakCount < peakCapacity) {
    let bestRank: i32 = -1;
    let bestStatistic: f64 = 0.0;
    for (let rank: i32 = 0; rank < informative; rank += 1) {
      const statistic = load<f64>(usize(chiPtr + rank * 8));
      if (statistic > bestStatistic) {
        bestStatistic = statistic;
        bestRank = rank;
      }
    }
    if (bestRank < 0 || !(bestStatistic > 0.0)) break;
    const recordPtr = peakPtr + peakCount * SOURCE_CHI_PEAK_INTS * 4;
    sourceChiGrowPeak(
      informative,
      bestRank,
      halfWindow,
      scoreBit,
      circular,
      positionsPtr,
      scoresPtr,
      missingPrefixPtr,
      bestStatistic,
      recordPtr + 12,
    );
    sourceChiStorePeakField(peakPtr, peakCount, 0, bestRank);
    sourceChiStorePeakField(
      peakPtr,
      peakCount,
      1,
      sourceChiBoundaryPosition(positionsPtr, informative, bestRank, nSites),
    );
    sourceChiStorePeakField(
      peakPtr,
      peakCount,
      2,
      sourceChiPeakSign(scoresPtr, informative, bestRank, halfWindow, scoreBit, circular),
    );
    sourceChiStorePeakField(peakPtr, peakCount, 5, sourceChiScaled(bestStatistic));
    peakCount += 1;
    sourceChiSuppressPeak(informative, bestRank, chiPtr, smoothPtr);
  }

  // FindMChiP keeps the first strict maximum. Restore genomic order after the
  // descending source search so paired state changes form non-overlapping
  // segments.
  for (let index: i32 = 1; index < peakCount; index += 1) {
    const rank = sourceChiPeakField(peakPtr, index, 0);
    const position = sourceChiPeakField(peakPtr, index, 1);
    const sign = sourceChiPeakField(peakPtr, index, 2);
    const statistic = sourceChiPeakField(peakPtr, index, 3);
    const width = sourceChiPeakField(peakPtr, index, 4);
    const rawStatistic = sourceChiPeakField(peakPtr, index, 5);
    let insertion = index;
    while (insertion > 0 && sourceChiPeakField(peakPtr, insertion - 1, 0) > rank) {
      for (let field: i32 = 0; field < SOURCE_CHI_PEAK_INTS; field += 1) {
        sourceChiStorePeakField(
          peakPtr,
          insertion,
          field,
          sourceChiPeakField(peakPtr, insertion - 1, field),
        );
      }
      insertion -= 1;
    }
    sourceChiStorePeakField(peakPtr, insertion, 0, rank);
    sourceChiStorePeakField(peakPtr, insertion, 1, position);
    sourceChiStorePeakField(peakPtr, insertion, 2, sign);
    sourceChiStorePeakField(peakPtr, insertion, 3, statistic);
    sourceChiStorePeakField(peakPtr, insertion, 4, width);
    sourceChiStorePeakField(peakPtr, insertion, 5, rawStatistic);
  }
  return peakCount;
}

function sourceChiWriteRow(
  outPtr: i32,
  outCapacity: i32,
  outputIndex: i32,
  method: i32,
  targetSlot: i32,
  track: i32,
  start: i32,
  end: i32,
  wraps: i32,
  statistic: i32,
  leftStatistic: i32,
  rightStatistic: i32,
  informative: i32,
  halfWindow: i32,
  leftRank: i32,
  rightRank: i32,
  leftWidth: i32,
  rightWidth: i32,
  direction: i32,
): void {
  if (outputIndex >= outCapacity) return;
  const row = outPtr + outputIndex * SOURCE_CHI_ROW_INTS * 4;
  store<i32>(usize(row), method);
  store<i32>(usize(row + 4), targetSlot);
  store<i32>(usize(row + 8), track);
  store<i32>(usize(row + 12), start);
  store<i32>(usize(row + 16), end);
  store<i32>(usize(row + 20), wraps);
  store<i32>(usize(row + 24), statistic);
  store<i32>(usize(row + 28), leftStatistic);
  store<i32>(usize(row + 32), rightStatistic);
  store<i32>(usize(row + 36), informative);
  store<i32>(usize(row + 40), halfWindow);
  store<i32>(usize(row + 44), leftRank);
  store<i32>(usize(row + 48), rightRank);
  store<i32>(usize(row + 52), leftWidth);
  store<i32>(usize(row + 56), rightWidth);
  store<i32>(usize(row + 60), direction);
}

function sourceChiPairPeaks(
  method: i32,
  targetSlot: i32,
  track: i32,
  nSites: i32,
  informative: i32,
  halfWindow: i32,
  circular: bool,
  peakPtr: i32,
  peakCount: i32,
  outPtr: i32,
  outCapacity: i32,
  outputCount: i32,
): i32 {
  let left: i32 = 0;
  while (left + 1 < peakCount) {
    const leftSign = sourceChiPeakField(peakPtr, left, 2);
    let right = left + 1;
    while (right < peakCount && sourceChiPeakField(peakPtr, right, 2) == leftSign) right += 1;
    if (right >= peakCount) break;
    let bestLeft = left;
    for (let candidate = left + 1; candidate < right; candidate += 1) {
      if (sourceChiPeakField(peakPtr, candidate, 3) > sourceChiPeakField(peakPtr, bestLeft, 3)) {
        bestLeft = candidate;
      }
    }
    const start = sourceChiPeakField(peakPtr, bestLeft, 1);
    const end = sourceChiPeakField(peakPtr, right, 1);
    if (end - start >= 4) {
      const leftStatistic = sourceChiPeakField(peakPtr, bestLeft, 3);
      const rightStatistic = sourceChiPeakField(peakPtr, right, 3);
      sourceChiWriteRow(
        outPtr,
        outCapacity,
        outputCount,
        method,
        targetSlot,
        track,
        start,
        end,
        0,
        leftStatistic < rightStatistic ? leftStatistic : rightStatistic,
        leftStatistic,
        rightStatistic,
        informative,
        halfWindow,
        sourceChiPeakField(peakPtr, bestLeft, 0),
        sourceChiPeakField(peakPtr, right, 0),
        sourceChiPeakField(peakPtr, bestLeft, 4),
        sourceChiPeakField(peakPtr, right, 4),
        leftSign,
      );
      outputCount += 1;
    }
    left = right + 1;
  }
  if (circular && peakCount > 1) {
    const last = peakCount - 1;
    const lastSign = sourceChiPeakField(peakPtr, last, 2);
    const firstSign = sourceChiPeakField(peakPtr, 0, 2);
    if (lastSign != firstSign) {
      const start = sourceChiPeakField(peakPtr, last, 1);
      const end = sourceChiPeakField(peakPtr, 0, 1);
      const circularLength = nSites - start + end;
      if (circularLength >= 4 && circularLength < nSites - 4) {
        const leftStatistic = sourceChiPeakField(peakPtr, last, 3);
        const rightStatistic = sourceChiPeakField(peakPtr, 0, 3);
        sourceChiWriteRow(
          outPtr,
          outCapacity,
          outputCount,
          method,
          targetSlot,
          track,
          start,
          end,
          1,
          leftStatistic < rightStatistic ? leftStatistic : rightStatistic,
          leftStatistic,
          rightStatistic,
          informative,
          halfWindow,
          sourceChiPeakField(peakPtr, last, 0),
          sourceChiPeakField(peakPtr, 0, 0),
          sourceChiPeakField(peakPtr, last, 4),
          sourceChiPeakField(peakPtr, 0, 4),
          lastSign,
        );
        outputCount += 1;
      }
    }
  }
  return outputCount;
}

function sourceChiScanTrack(
  method: i32,
  targetSlot: i32,
  track: i32,
  nSites: i32,
  informative: i32,
  fullWindow: i32,
  circular: bool,
  scoreBit: i32,
  positionsPtr: i32,
  scoresPtr: i32,
  missingPrefixPtr: i32,
  chiPtr: i32,
  smoothPtr: i32,
  peakPtr: i32,
  peakCapacity: i32,
  outPtr: i32,
  outCapacity: i32,
  outputCount: i32,
): i32 {
  const halfWindow = sourceChiHalfWindow(informative, fullWindow);
  if (halfWindow <= 0) return outputCount;
  sourceChiFillProfile(
    informative,
    halfWindow,
    circular,
    scoreBit,
    positionsPtr,
    scoresPtr,
    missingPrefixPtr,
    chiPtr,
    smoothPtr,
  );
  const peakCount = sourceChiCollectPeaks(
    nSites,
    informative,
    halfWindow,
    circular,
    scoreBit,
    positionsPtr,
    scoresPtr,
    missingPrefixPtr,
    chiPtr,
    smoothPtr,
    peakPtr,
    peakCapacity,
  );
  return sourceChiPairPeaks(
    method,
    targetSlot,
    track,
    nSites,
    informative,
    halfWindow,
    circular,
    peakPtr,
    peakCount,
    outPtr,
    outCapacity,
    outputCount,
  );
}

function sourceChiScanBuiltTracks(
  nSites: i32,
  fullWindow: i32,
  circularFlag: i32,
  methodMask: i32,
  positionsPtr: i32,
  scoresPtr: i32,
  missingPrefixPtr: i32,
  chiPtr: i32,
  smoothPtr: i32,
  peakPtr: i32,
  peakCapacity: i32,
  outPtr: i32,
  outCapacity: i32,
): i32 {
  let outputCount: i32 = 0;
  const circular = circularFlag != 0;
  const maxChiInformative = load<i32>(usize(peakPtr));
  const chimaera0Informative = load<i32>(usize(peakPtr + 4));
  const chimaera1Informative = load<i32>(usize(peakPtr + 8));
  const chimaera2Informative = load<i32>(usize(peakPtr + 12));
  if ((methodMask & 4) != 0) {
    for (let track: i32 = 0; track < 3; track += 1) {
      outputCount = sourceChiScanTrack(
        3,
        -1,
        track,
        nSites,
        maxChiInformative,
        fullWindow,
        circular,
        track,
        positionsPtr,
        scoresPtr,
        missingPrefixPtr,
        chiPtr,
        smoothPtr,
        peakPtr,
        peakCapacity,
        outPtr,
        outCapacity,
        outputCount,
      );
    }
  }
  if ((methodMask & 8) != 0) {
    const positionStride = (nSites + 1) * 4;
    const scoreStride = nSites + 1;
    for (let targetSlot: i32 = 0; targetSlot < 3; targetSlot += 1) {
      const informative = targetSlot == 0
        ? chimaera0Informative
        : targetSlot == 1 ? chimaera1Informative : chimaera2Informative;
      const trackPositionsPtr = positionsPtr + positionStride * (targetSlot + 1);
      const trackScoresPtr = scoresPtr + scoreStride * (targetSlot + 1);
      outputCount = sourceChiScanTrack(
        4,
        targetSlot,
        targetSlot,
        nSites,
        informative,
        fullWindow,
        circular,
        0,
        trackPositionsPtr,
        trackScoresPtr,
        missingPrefixPtr,
        chiPtr,
        smoothPtr,
        peakPtr,
        peakCapacity,
        outPtr,
        outCapacity,
        outputCount,
      );
    }
  }
  return outputCount;
}

// Output rows (16 i32 each): method (3 MAXCHI, 4 CHIMAERA), target slot
// (-1 for MAXCHI), source track, start/end/wrap, minimum and individual grown
// boundary chi-square values ×1000, compressed-site count, half-window,
// boundary ranks, GrowMChiWin widths, and left-boundary direction. The return
// value is the total row count before output-capacity truncation. positionsPtr
// and scoresPtr each contain four (nSites+1)-element scratch lanes: shared
// MAXCHI followed by CHIMAERA target slots 0, 1 and 2.
export function scan_source_chi_all(
  seqPtr: i32,
  nSites: i32,
  sequence0: i32,
  sequence1: i32,
  sequence2: i32,
  fullWindow: i32,
  circularFlag: i32,
  methodMask: i32,
  positionsPtr: i32,
  scoresPtr: i32,
  missingPrefixPtr: i32,
  chiPtr: i32,
  smoothPtr: i32,
  peakPtr: i32,
  peakCapacity: i32,
  outPtr: i32,
  outCapacity: i32,
): i32 {
  sourceChiBuildAllTracks(
    seqPtr,
    nSites,
    sequence0,
    sequence1,
    sequence2,
    positionsPtr,
    scoresPtr,
    missingPrefixPtr,
    peakPtr,
    methodMask,
  );
  return sourceChiScanBuiltTracks(
    nSites,
    fullWindow,
    circularFlag,
    methodMask,
    positionsPtr,
    scoresPtr,
    missingPrefixPtr,
    chiPtr,
    smoothPtr,
    peakPtr,
    peakCapacity,
    outPtr,
    outCapacity,
  );
}

// Packed production entry point. It is deliberately separate from the byte
// oracle above so tests can enforce bit-for-bit equivalence on every emitted
// source peak while production scans avoid three raw-sequence loads per site.
export function scan_source_chi_all_packed(
  packedPtr: i32,
  validityPtr: i32,
  wordsPerSequence: i32,
  nSites: i32,
  sequence0: i32,
  sequence1: i32,
  sequence2: i32,
  fullWindow: i32,
  circularFlag: i32,
  methodMask: i32,
  positionsPtr: i32,
  scoresPtr: i32,
  missingPrefixPtr: i32,
  chiPtr: i32,
  smoothPtr: i32,
  peakPtr: i32,
  peakCapacity: i32,
  outPtr: i32,
  outCapacity: i32,
): i32 {
  sourceChiBuildAllTracksPacked(
    packedPtr,
    validityPtr,
    wordsPerSequence,
    nSites,
    sequence0,
    sequence1,
    sequence2,
    positionsPtr,
    scoresPtr,
    missingPrefixPtr,
    peakPtr,
    methodMask,
  );
  return sourceChiScanBuiltTracks(
    nSites,
    fullWindow,
    circularFlag,
    methodMask,
    positionsPtr,
    scoresPtr,
    missingPrefixPtr,
    chiPtr,
    smoothPtr,
    peakPtr,
    peakCapacity,
    outPtr,
    outCapacity,
  );
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
  packedPtr: i32,
  validityPtr: i32,
  wordsPerSequence: i32,
  usePacked: bool,
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

  if (usePacked) {
    const offset0 = seq1 * wordsPerSequence;
    const offset1 = seq2 * wordsPerSequence;
    const offset2 = seq3 * wordsPerSequence;
    for (let word: i32 = 0; word < wordsPerSequence; word += 1) {
      const a = load<u32>(usize(packedPtr + (offset0 + word) * 4));
      const b = load<u32>(usize(packedPtr + (offset1 + word) * 4));
      const c = load<u32>(usize(packedPtr + (offset2 + word) * 4));
      const allValid = load<u32>(usize(validityPtr + (offset0 + word) * 4))
        & load<u32>(usize(validityPtr + (offset1 + word) * 4))
        & load<u32>(usize(validityPtr + (offset2 + word) * 4));
      const ab = packedEquality(a, b, allValid);
      const ac = packedEquality(a, c, allValid);
      const bc = packedEquality(b, c, allValid);
      const category0 = ab & ~ac;
      const category1 = ac & ~ab;
      const category2 = bc & ~ab;
      count0 += popCount32(category0);
      count1 += popCount32(category1);
      count2 += popCount32(category2);
      let variable = category0 | category1 | category2;
      while (variable != 0) {
        const bit = i32(ctz<u32>(variable));
        const lane = u32(1) << u32(bit);
        const site = word * 16 + bit / 2;
        if (site >= nSites) break;
        const category = (category0 & lane) != 0 ? 0 : (category1 & lane) != 0 ? 1 : 2;
        store<i32>(usize(scratchPackedPtr + informative * 4), (site << 2) | category);
        informative += 1;
        variable &= variable - 1;
      }
    }
  } else {
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
  return scanRdp5TripletCore(seqPtr, 0, 0, 0, false, nSites, seq1, seq2, seq3, fullWindow, scratchPackedPtr, scratchDominancePtr, outPtr, 0, 0);
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
  return scanRdp5TripletCore(seqPtr, 0, 0, 0, false, nSites, seq1, seq2, seq3, fullWindow, scratchPackedPtr, scratchDominancePtr, bestOutPtr, outPtr, maximumEvents);
}

export function scan_rdp5_triplet_all_packed(
  packedPtr: i32,
  validityPtr: i32,
  wordsPerSequence: i32,
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
  return scanRdp5TripletCore(0, packedPtr, validityPtr, wordsPerSequence, true, nSites, seq1, seq2, seq3, fullWindow, scratchPackedPtr, scratchDominancePtr, bestOutPtr, outPtr, maximumEvents);
}
