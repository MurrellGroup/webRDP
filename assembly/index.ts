// Clean-room WebAssembly kernels for RDP Web.
// The implementation is derived from published triplet-scanning principles,
// not from RDP/OpenRDP source code.

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

// Calculates independent fast statistics for the seven exploratory method
// families after a candidate tract has been localized. These are clean-room
// implementations of the procedures described in the RDP5 manual. Exact
// parity calibrations that are not redistributable (notably GENECONV and
// 3SEQ tables/executables) are deliberately not embedded.
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
  prefixAPtr: i32,
  prefixBPtr: i32,
  outPtr: i32,
): void {
  const halfWindow = window / 2 > 8 ? window / 2 : 8;

  // GENECONV G-scale 0: longest uninterrupted concordant fragment in
  // triplet-polymorphic space. Monomorphic sites are ignored, as documented.
  let eligible: i32 = 0;
  let concordant: i32 = 0;
  let run: i32 = 0;
  let bestRun: i32 = 0;
  let runStart: i32 = start;
  let bestRunStart: i32 = start;
  let bestRunEnd: i32 = start;
  for (let site: i32 = 0; site < nSites; site += 1) {
    const r = seqBase(seqPtr, nSites, recombinant, site);
    const a = seqBase(seqPtr, nSites, majorParent, site);
    const b = seqBase(seqPtr, nSites, minorParent, site);
    if (!valid(r) || !valid(a) || !valid(b) || (r == a && r == b)) continue;
    eligible += 1;
    if (r == b) concordant += 1;
    if (site < start || site >= end) continue;
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

  // Prefix pairwise mismatches for a deterministic, no-resampling RECSCAN
  // topology-switch statistic. Only sites valid in the full triplet are used
  // so both distances have identical denominators.
  store<i32>(usize(prefixAPtr), 0);
  store<i32>(usize(prefixBPtr), 0);
  let majorDifferences: i32 = 0;
  let minorDifferences: i32 = 0;
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
  let topologyConsistent: i32 = 0;
  let decisiveWindows: i32 = 0;
  const stride = step > 0 ? step : 1;
  for (let center = halfWindow; center <= nSites - halfWindow; center += stride) {
    const left = center - halfWindow;
    const right = center + halfWindow;
    const majorDistance = prefixValue(prefixAPtr, right) - prefixValue(prefixAPtr, left);
    const minorDistance = prefixValue(prefixBPtr, right) - prefixValue(prefixBPtr, left);
    if (majorDistance == minorDistance) continue;
    decisiveWindows += 1;
    const inside = center >= start && center < end;
    if ((inside && minorDistance < majorDistance) || (!inside && majorDistance < minorDistance)) {
      topologyConsistent += 1;
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
  for (let region: i32 = 0; region < 3; region += 1) {
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
  store<i32>(usize(prefixAPtr), 0);
  store<i32>(usize(prefixBPtr), 0);
  let pairVariable: i32 = 0;
  let pairNonVariable: i32 = 0;
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
  const maxChiStart = breakpointChi(prefixAPtr, prefixBPtr, nSites, start, halfWindow);
  const maxChiEnd = breakpointChi(prefixAPtr, prefixBPtr, nSites, end, halfWindow);
  const maxChiStatistic = maxChiStart < maxChiEnd ? maxChiStart : maxChiEnd;

  // CHIMAERA's compressed binary string and SISCAN-style oriented category
  // score. The same prefixes support an O(window) breakpoint-polishing pass.
  store<i32>(usize(prefixAPtr), 0);
  store<i32>(usize(prefixBPtr), 0);
  let matchesMajor: i32 = 0;
  let matchesMinor: i32 = 0;
  let sisterScore: i32 = 0;
  let sisterSites: i32 = 0;
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
      if (category != 0) {
        sisterSites += 1;
        sisterScore += (site >= start && site < end) ? category : -category;
      }
    }
    store<i32>(usize(prefixAPtr + (site + 1) * 4), matchesMajor);
    store<i32>(usize(prefixBPtr + (site + 1) * 4), matchesMinor);
  }
  const chimaeraStart = breakpointChi(prefixAPtr, prefixBPtr, nSites, start, halfWindow);
  const chimaeraEnd = breakpointChi(prefixAPtr, prefixBPtr, nSites, end, halfWindow);
  const chimaeraStatistic = chimaeraStart < chimaeraEnd ? chimaeraStart : chimaeraEnd;
  // 3SEQ maps sites matching the major parent to an up-step and sites
  // matching the minor parent to a down-step. Its two-breakpoint statistic is
  // the maximum descent of that hypergeometric random walk from any previous
  // maximum, corresponding to the strongest contiguous minor-parent tract.
  const threeSeqSites = matchesMajor + matchesMinor;
  let threeSeqWalk: i32 = 0;
  let threeSeqMaximum: i32 = 0;
  let threeSeqDescent: i32 = 0;
  for (let position: i32 = 1; position <= nSites; position += 1) {
    const majorStep = prefixValue(prefixAPtr, position) - prefixValue(prefixAPtr, position - 1);
    const minorStep = prefixValue(prefixBPtr, position) - prefixValue(prefixBPtr, position - 1);
    threeSeqWalk += majorStep - minorStep;
    if (threeSeqWalk > threeSeqMaximum) threeSeqMaximum = threeSeqWalk;
    const descent = threeSeqMaximum - threeSeqWalk;
    if (descent > threeSeqDescent) threeSeqDescent = descent;
  }
  const polishedStart = polishBreakpoint(prefixAPtr, prefixBPtr, nSites, start, halfWindow);
  const polishedEnd = polishBreakpoint(prefixAPtr, prefixBPtr, nSites, end, halfWindow);

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
}

@inline
function hmmProbability(value: f64): f64 {
  return value < 0.5001 ? 0.5001 : value > 0.9999 ? 0.9999 : value;
}

@inline
function markedPosition(value: i32): i32 {
  return value < 0 ? -value - 1 : value;
}

// Candidate-seeded, windowless two-state HMM breakpoint refinement. Following
// the public BURT description, only sites where the proposed parents differ
// and the recombinant matches one parent are retained. A Viterbi pass labels
// those sites as major- or minor-parent ancestry and the minor-state run with
// the greatest overlap with the candidate tract supplies the refined bounds.
//
// This is a clean-room two-state precursor, not a claim of RDP5 BURT parity:
// BURT's documented 2..20-state step-up fitting and repeated random-start
// Viterbi training remain validation work. Runtime and scratch space are O(L).
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
