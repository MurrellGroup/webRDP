// RDP5 automated BOOTSCAN/RECSCAN distance-mode batch.
//
// This follows the author-supplied BSXoverR2 -> SEQBOOT2 -> FastBootDistIP ->
// GetPltVal2 path. The unit of work is the active sequence set, not an
// oriented query/reference triplet: one bootstrap table is reused at every
// window, each required pair is evaluated once, then every concrete triplet's
// three pair relationships are interpreted together.

@inline
function valid(base: u8): bool {
  return base < 4;
}

@inline
function byteBase(seqPtr: i32, nSites: i32, sequence: i32, site: i32): u8 {
  return load<u8>(usize(seqPtr + sequence * nSites + site));
}

@inline
function bounded(value: i32, lower: i32, upper: i32): i32 {
  return value < lower ? lower : value > upper ? upper : value;
}

@inline
function sourceBootBase(
  packedMode: i32,
  dataPtr: i32,
  validityPtr: i32,
  wordsPerSequence: i32,
  nSites: i32,
  sequence: i32,
  site: i32,
): u8 {
  if (packedMode == 0) return byteBase(dataPtr, nSites, sequence, site);
  const word = sequence * wordsPerSequence + (site >> 4);
  const shift = (site & 15) << 1;
  const validWord = load<u32>(usize(validityPtr + word * 4));
  if ((validWord & (u32(1) << shift)) == 0) return 4;
  return u8((load<u32>(usize(dataPtr + word * 4)) >> shift) & 3);
}

@inline
function pairIndex(left: i32, right: i32, nSeq: i32): i32 {
  let a = left;
  let b = right;
  if (a > b) {
    const swap = a;
    a = b;
    b = swap;
  }
  return (a * (2 * nSeq - a - 1)) / 2 + b - a - 1;
}

@inline
function quantizedDistance(validSites: i32, differences: i32): u16 {
  if (validSites <= 0) return 32000;
  const identity = f32(validSites - differences) / f32(validSites);
  let distance: f32 = 10.0;
  if (identity > 0.25) {
    const transformed = f32((4.0 * identity - 1.0) / 3.0);
    distance = f32(-0.75 * Math.log(f64(transformed)));
  }
  let scaled = i32(distance * 3200.0);
  if (scaled < 0) scaled = 0;
  if (scaled > 32000) scaled = 32000;
  return u16(scaled);
}

@inline
function globalIdentityGreater(globalPairPtr: i32, leftPair: i32, rightPair: i32): bool {
  const leftValid = load<i32>(usize(globalPairPtr + leftPair * 8));
  const leftDifferences = load<i32>(usize(globalPairPtr + leftPair * 8 + 4));
  const rightValid = load<i32>(usize(globalPairPtr + rightPair * 8));
  const rightDifferences = load<i32>(usize(globalPairPtr + rightPair * 8 + 4));
  const leftMatches = leftValid > 0 ? leftValid - leftDifferences : 0;
  const rightMatches = rightValid > 0 ? rightValid - rightDifferences : 0;
  if (leftValid <= 0) return false;
  if (rightValid <= 0) return leftMatches > 0;
  return i64(leftMatches) * i64(rightValid) > i64(rightMatches) * i64(leftValid);
}

@inline
function globalClosestPair(globalPairPtr: i32, pair0: i32, pair1: i32, pair2: i32): i32 {
  // ScanBSPlots orders relationships with the whole-alignment Distance
  // matrix, whose entries are pairwise identities rather than the window's
  // JC distance. Compare the exact fractions and preserve the source's tie
  // order; JC saturation must not scramble a divergent triplet's baseline.
  let closest: i32 = 0;
  if (globalIdentityGreater(globalPairPtr, pair1, pair0)) {
    closest = 1;
    if (globalIdentityGreater(globalPairPtr, pair2, pair1)) closest = 2;
  } else if (globalIdentityGreater(globalPairPtr, pair2, pair0)) {
    closest = 2;
  }
  return closest;
}

function buildWeights(window: i32, replicates: i32, seed: u32, weightPtr: i32): void {
  for (let position: i32 = 0; position < window; position += 1) {
    const row = weightPtr + position * replicates * 2;
    store<u16>(usize(row), 1);
    for (let replicate: i32 = 1; replicate < replicates; replicate += 1) {
      store<u16>(usize(row + replicate * 2), 0);
    }
  }

  // MSVCRT rand(), used by the supplied Windows build. SEQBOOT2 discards two
  // draws and consumes an additional unused replicate draw per sampled site.
  let state = seed;
  state = state * 214013 + 2531011;
  state = state * 214013 + 2531011;
  const span = window - 1;
  for (let draw: i32 = 0; draw < window; draw += 1) {
    for (let replicate: i32 = 1; replicate <= replicates; replicate += 1) {
      state = state * 214013 + 2531011;
      const randomValue = (state >> 16) & 0x7fff;
      let sampled = i32((f64(randomValue) / 32767.0) * f64(span));
      if (sampled < 0) sampled = 0;
      if (sampled >= window) sampled = window - 1;
      if (replicate < replicates) {
        const address = weightPtr + (sampled * replicates + replicate) * 2;
        store<u16>(usize(address), load<u16>(usize(address)) + 1);
      }
    }
  }
}

function buildDistanceLookup(window: i32, lookupPtr: i32): void {
  for (let validSites: i32 = 0; validSites <= window; validSites += 1) {
    const row = (validSites * (validSites + 1)) / 2;
    for (let differences: i32 = 0; differences <= validSites; differences += 1) {
      store<u16>(
        usize(lookupPtr + (row + differences) * 2),
        quantizedDistance(validSites, differences),
      );
    }
  }
}

@inline
function lookupDistance(lookupPtr: i32, validSites: i32, differences: i32): u16 {
  return load<u16>(usize(
    lookupPtr + (((validSites * (validSites + 1)) / 2 + differences) * 2),
  ));
}

@inline
function pairMatches(topology: i32, first: u8, second: u8, third: u8): bool {
  if (topology == 0) return first == second;
  if (topology == 1) return first == third;
  return second == third;
}

function emitSignal(
  packedMode: i32,
  dataPtr: i32,
  validityPtr: i32,
  wordsPerSequence: i32,
  nSites: i32,
  first: i32,
  second: i32,
  third: i32,
  topology: i32,
  globalClosest: i32,
  startWindow: i32,
  lastHighWindow: i32,
  maximumSupport: i32,
  window: i32,
  step: i32,
  replicates: i32,
  outPtr: i32,
  outCapacity: i32,
  signalCount: i32,
): i32 {
  if (startWindow < 0 || lastHighWindow < startWindow) return signalCount;
  let start = startWindow * step - window / 2;
  let end = lastHighWindow * step + window / 2;
  if (start < 0) start = 0;
  if (end > nSites) end = nSites;
  if (end - start < 4 || (start == 0 && end == nSites)) return signalCount;

  let tractMatches: i32 = 0;
  let backgroundMatches: i32 = 0;
  let tractInformative: i32 = 0;
  let totalInformative: i32 = 0;
  for (let site: i32 = 0; site < nSites; site += 1) {
    const a = sourceBootBase(packedMode, dataPtr, validityPtr, wordsPerSequence, nSites, first, site);
    const b = sourceBootBase(packedMode, dataPtr, validityPtr, wordsPerSequence, nSites, second, site);
    const c = sourceBootBase(packedMode, dataPtr, validityPtr, wordsPerSequence, nSites, third, site);
    if (!valid(a) || !valid(b) || !valid(c)) continue;
    if (a == b && a == c) continue;
    totalInformative += 1;
    const inside = site >= start && site < end;
    const matching = pairMatches(topology, a, b, c);
    if (inside) {
      tractInformative += 1;
      if (matching) tractMatches += 1;
    } else if (matching) {
      backgroundMatches += 1;
    }
  }
  if (tractInformative <= 2 || totalInformative <= tractInformative || tractMatches <= 0) return signalCount;
  // ProbCalc's upper tail cannot be significant unless the selected pair is
  // enriched inside the putative tract. This exact guard prevents null runs
  // from exhausting the finite result buffer.
  if (tractMatches * totalInformative <= (tractMatches + backgroundMatches) * tractInformative) {
    return signalCount;
  }

  if (signalCount < outCapacity) {
    const row = outPtr + signalCount * 16 * 4;
    store<i32>(usize(row), first);
    store<i32>(usize(row + 4), second);
    store<i32>(usize(row + 8), third);
    store<i32>(usize(row + 12), topology);
    store<i32>(usize(row + 16), start);
    store<i32>(usize(row + 20), end);
    store<i32>(usize(row + 24), maximumSupport);
    store<i32>(usize(row + 28), replicates);
    store<i32>(usize(row + 32), lastHighWindow - startWindow + 1);
    store<i32>(usize(row + 36), tractMatches);
    store<i32>(usize(row + 40), backgroundMatches);
    store<i32>(usize(row + 44), tractInformative);
    store<i32>(usize(row + 48), totalInformative);
    store<i32>(usize(row + 52), globalClosest);
    store<i32>(usize(row + 56), window);
    store<i32>(usize(row + 60), step);
  }
  return signalCount + 1;
}

function scanBatchCore(
  packedMode: i32,
  dataPtr: i32,
  validityPtr: i32,
  wordsPerSequence: i32,
  nSeq: i32,
  nSites: i32,
  tripletPtr: i32,
  tripletCount: i32,
  windowInput: i32,
  stepInput: i32,
  replicatesInput: i32,
  cutoffPermilleInput: i32,
  seed: u32,
  pairMapPtr: i32,
  pairListPtr: i32,
  weightPtr: i32,
  pairDistancePtr: i32,
  globalPairPtr: i32,
  statePtr: i32,
  differencePtr: i32,
  validPtr: i32,
  lookupPtr: i32,
  outPtr: i32,
  outCapacity: i32,
): i32 {
  if (nSeq < 3 || nSites < 5 || tripletCount <= 0 || outCapacity <= 0) return 0;
  const maximumWindow = nSites / 2 > 5 ? nSites / 2 : 5;
  let window = bounded(windowInput, 5, maximumWindow);
  if (window > 32767) window = 32767;
  const maximumStep = nSites / 4 > 1 ? nSites / 4 : 1;
  const step = bounded(stepInput, 1, maximumStep);
  const replicates = bounded(replicatesInput, 2, 1000);
  const cutoffPermille = bounded(cutoffPermilleInput, 500, 999);
  const pairCount = (nSeq * (nSeq - 1)) / 2;
  memory.fill(pairMapPtr, 0xff, pairCount * 4);
  for (let triplet: i32 = 0; triplet < tripletCount; triplet += 1) {
    const row = tripletPtr + triplet * 12;
    const first = load<i32>(usize(row));
    const second = load<i32>(usize(row + 4));
    const third = load<i32>(usize(row + 8));
    if (first < 0 || second <= first || third <= second || third >= nSeq) continue;
    store<i32>(usize(pairMapPtr + pairIndex(first, second, nSeq) * 4), 0);
    store<i32>(usize(pairMapPtr + pairIndex(first, third, nSeq) * 4), 0);
    store<i32>(usize(pairMapPtr + pairIndex(second, third, nSeq) * 4), 0);
  }

  // Compact the sparse pair set once. Exact all-triplet scans naturally use
  // every pair; approximate/query scans no longer walk an O(N^2) matrix at
  // every window or reserve bootstrap rows for pairs that were not requested.
  let usedPairCount: i32 = 0;
  let densePair: i32 = 0;
  for (let first: i32 = 0; first < nSeq - 1; first += 1) {
    for (let second: i32 = first + 1; second < nSeq; second += 1) {
      const mapAddress = pairMapPtr + densePair * 4;
      if (load<i32>(usize(mapAddress)) >= 0) {
        store<i32>(usize(mapAddress), usedPairCount);
        store<i32>(usize(pairListPtr + usedPairCount * 8), first);
        store<i32>(usize(pairListPtr + usedPairCount * 8 + 4), second);
        usedPairCount += 1;
      }
      densePair += 1;
    }
  }
  if (usedPairCount <= 0) return 0;

  buildWeights(window, replicates, seed, weightPtr);
  buildDistanceLookup(window, lookupPtr);
  for (let usedPair: i32 = 0; usedPair < usedPairCount; usedPair += 1) {
    const first = load<i32>(usize(pairListPtr + usedPair * 8));
    const second = load<i32>(usize(pairListPtr + usedPair * 8 + 4));
    let validSites: i32 = 0;
    let differences: i32 = 0;
    for (let site: i32 = 0; site < nSites; site += 1) {
      const a = sourceBootBase(packedMode, dataPtr, validityPtr, wordsPerSequence, nSites, first, site);
      const b = sourceBootBase(packedMode, dataPtr, validityPtr, wordsPerSequence, nSites, second, site);
      if (!valid(a) || !valid(b)) continue;
      validSites += 1;
      if (a != b) differences += 1;
    }
    store<i32>(usize(globalPairPtr + usedPair * 8), validSites);
    store<i32>(usize(globalPairPtr + usedPair * 8 + 4), differences);
  }

  // Two alternate tracks per triplet. Each stores start, last full-support
  // window and maximum bootstrap support. The whole-alignment closest pair is
  // the baseline topology and is not emitted as a recombinant tract.
  for (let triplet: i32 = 0; triplet < tripletCount; triplet += 1) {
    const stateRow = statePtr + triplet * 24;
    store<i32>(usize(stateRow), -1);
    store<i32>(usize(stateRow + 4), -1);
    store<i32>(usize(stateRow + 8), 0);
    store<i32>(usize(stateRow + 12), -1);
    store<i32>(usize(stateRow + 16), -1);
    store<i32>(usize(stateRow + 20), 0);
  }

  const overlapWindows = window / step > 1 ? window / step : 1;
  // BSXoverR2 defines NumWins as floor(L / step) + 2 and retains that many
  // distance windows after its on-disk transpose. The final wrapped window is
  // boundary context for the source begin/end scan.
  const windowCount = nSites / step + 2;
  let signalCount: i32 = 0;
  for (let windowIndex: i32 = 0; windowIndex < windowCount; windowIndex += 1) {
    const rawWindowStart = windowIndex * step - window / 2;
    for (let usedPair: i32 = 0; usedPair < usedPairCount; usedPair += 1) {
      const first = load<i32>(usize(pairListPtr + usedPair * 8));
      const second = load<i32>(usize(pairListPtr + usedPair * 8 + 4));
      for (let replicate: i32 = 0; replicate < replicates; replicate += 1) {
        store<i32>(usize(differencePtr + replicate * 4), 0);
        store<i32>(usize(validPtr + replicate * 4), window);
      }
      let completePair = true;
      // Canonical high-identity alignments spend almost all their time on
      // invariant pair positions. Their weights only affect the fixed
      // denominator, so visit the replicate table only at mismatches. If a
      // missing base is encountered, restart this pair on the exact fallback.
      for (let relative: i32 = 0; relative < window; relative += 1) {
        let site = rawWindowStart + relative;
        while (site < 0) site += nSites;
        while (site >= nSites) site -= nSites;
        const a = sourceBootBase(packedMode, dataPtr, validityPtr, wordsPerSequence, nSites, first, site);
        const b = sourceBootBase(packedMode, dataPtr, validityPtr, wordsPerSequence, nSites, second, site);
        if (!valid(a) || !valid(b)) {
          completePair = false;
          break;
        }
        if (a == b) continue;
        const weights = weightPtr + relative * replicates * 2;
        for (let replicate: i32 = 0; replicate < replicates; replicate += 1) {
          const weight = i32(load<u16>(usize(weights + replicate * 2)));
          if (weight == 0) continue;
          const differenceAddress = differencePtr + replicate * 4;
          store<i32>(usize(differenceAddress), load<i32>(usize(differenceAddress)) + weight);
        }
      }
      if (!completePair) {
        for (let replicate: i32 = 0; replicate < replicates; replicate += 1) {
          store<i32>(usize(differencePtr + replicate * 4), 0);
          store<i32>(usize(validPtr + replicate * 4), 0);
        }
        for (let relative: i32 = 0; relative < window; relative += 1) {
          let site = rawWindowStart + relative;
          while (site < 0) site += nSites;
          while (site >= nSites) site -= nSites;
          const a = sourceBootBase(packedMode, dataPtr, validityPtr, wordsPerSequence, nSites, first, site);
          const b = sourceBootBase(packedMode, dataPtr, validityPtr, wordsPerSequence, nSites, second, site);
          if (!valid(a) || !valid(b)) continue;
          const different = a != b;
          const weights = weightPtr + relative * replicates * 2;
          for (let replicate: i32 = 0; replicate < replicates; replicate += 1) {
            const weight = i32(load<u16>(usize(weights + replicate * 2)));
            if (weight == 0) continue;
            const validAddress = validPtr + replicate * 4;
            store<i32>(usize(validAddress), load<i32>(usize(validAddress)) + weight);
            if (different) {
              const differenceAddress = differencePtr + replicate * 4;
              store<i32>(usize(differenceAddress), load<i32>(usize(differenceAddress)) + weight);
            }
          }
        }
      }
      const pairRow = pairDistancePtr + usedPair * replicates * 2;
      for (let replicate: i32 = 0; replicate < replicates; replicate += 1) {
        const validSites = load<i32>(usize(validPtr + replicate * 4));
        const differences = load<i32>(usize(differencePtr + replicate * 4));
        store<u16>(usize(pairRow + replicate * 2), lookupDistance(lookupPtr, validSites, differences));
      }
    }

    for (let triplet: i32 = 0; triplet < tripletCount; triplet += 1) {
      const tripletRow = tripletPtr + triplet * 12;
      const first = load<i32>(usize(tripletRow));
      const second = load<i32>(usize(tripletRow + 4));
      const third = load<i32>(usize(tripletRow + 8));
      if (first < 0 || second <= first || third <= second || third >= nSeq) continue;
      const pair0 = load<i32>(usize(pairMapPtr + pairIndex(first, second, nSeq) * 4));
      const pair1 = load<i32>(usize(pairMapPtr + pairIndex(first, third, nSeq) * 4));
      const pair2 = load<i32>(usize(pairMapPtr + pairIndex(second, third, nSeq) * 4));
      const globalClosest = globalClosestPair(globalPairPtr, pair0, pair1, pair2);

      let support0: i32 = 0;
      let support1: i32 = 0;
      let support2: i32 = 0;
      const distance0 = pairDistancePtr + pair0 * replicates * 2;
      const distance1 = pairDistancePtr + pair1 * replicates * 2;
      const distance2 = pairDistancePtr + pair2 * replicates * 2;
      for (let replicate: i32 = 0; replicate < replicates; replicate += 1) {
        const value0 = load<u16>(usize(distance0 + replicate * 2));
        const value1 = load<u16>(usize(distance1 + replicate * 2));
        const value2 = load<u16>(usize(distance2 + replicate * 2));
        // GetPltVal2's strict comparison/tie order.
        if (value0 < value1) {
          if (value0 < value2) support0 += 1;
          else if (value2 < value0) support2 += 1;
        } else if (value1 < value0) {
          if (value1 < value2) support1 += 1;
          else if (value2 < value1) support2 += 1;
        } else if (value2 < value0) {
          support2 += 1;
        }
      }

      for (let alternate: i32 = 0; alternate < 2; alternate += 1) {
        const topology = alternate == 0 ? (globalClosest + 1) % 3 : (globalClosest + 2) % 3;
        const support = topology == 0 ? support0 : topology == 1 ? support1 : support2;
        const otherA = topology == 0 ? support1 : support0;
        const otherB = topology == 2 ? support1 : support2;
        const stateRow = statePtr + (triplet * 6 + alternate * 3) * 4;
        let startWindow = load<i32>(usize(stateRow));
        let lastHighWindow = load<i32>(usize(stateRow + 4));
        let maximumSupport = load<i32>(usize(stateRow + 8));
        const high = support * 1000 >= cutoffPermille * replicates;
        const dominant = support > otherA && support > otherB && support * 10 > 4 * replicates;
        if (startWindow < 0) {
          if (high) {
            startWindow = windowIndex;
            lastHighWindow = windowIndex;
            maximumSupport = support;
          }
        } else if (high) {
          lastHighWindow = windowIndex;
          if (support > maximumSupport) maximumSupport = support;
        } else if (dominant && windowIndex - lastHighWindow <= overlapWindows) {
          if (support > maximumSupport) maximumSupport = support;
        } else {
          signalCount = emitSignal(
            packedMode, dataPtr, validityPtr, wordsPerSequence, nSites,
            first, second, third, topology, globalClosest, startWindow,
            lastHighWindow, maximumSupport, window, step, replicates,
            outPtr, outCapacity, signalCount,
          );
          startWindow = -1;
          lastHighWindow = -1;
          maximumSupport = 0;
          if (high) {
            startWindow = windowIndex;
            lastHighWindow = windowIndex;
            maximumSupport = support;
          }
        }
        store<i32>(usize(stateRow), startWindow);
        store<i32>(usize(stateRow + 4), lastHighWindow);
        store<i32>(usize(stateRow + 8), maximumSupport);
      }
    }
  }

  for (let triplet: i32 = 0; triplet < tripletCount; triplet += 1) {
    const tripletRow = tripletPtr + triplet * 12;
    const first = load<i32>(usize(tripletRow));
    const second = load<i32>(usize(tripletRow + 4));
    const third = load<i32>(usize(tripletRow + 8));
    if (first < 0 || second <= first || third <= second || third >= nSeq) continue;
    const pair0 = load<i32>(usize(pairMapPtr + pairIndex(first, second, nSeq) * 4));
    const pair1 = load<i32>(usize(pairMapPtr + pairIndex(first, third, nSeq) * 4));
    const pair2 = load<i32>(usize(pairMapPtr + pairIndex(second, third, nSeq) * 4));
    const globalClosest = globalClosestPair(globalPairPtr, pair0, pair1, pair2);
    for (let alternate: i32 = 0; alternate < 2; alternate += 1) {
      const topology = alternate == 0 ? (globalClosest + 1) % 3 : (globalClosest + 2) % 3;
      const stateRow = statePtr + (triplet * 6 + alternate * 3) * 4;
      signalCount = emitSignal(
        packedMode, dataPtr, validityPtr, wordsPerSequence, nSites,
        first, second, third, topology, globalClosest,
        load<i32>(usize(stateRow)), load<i32>(usize(stateRow + 4)),
        load<i32>(usize(stateRow + 8)), window, step, replicates,
        outPtr, outCapacity, signalCount,
      );
    }
  }
  return signalCount;
}

export function scan_source_bootscan_batch(
  seqPtr: i32,
  nSeq: i32,
  nSites: i32,
  tripletPtr: i32,
  tripletCount: i32,
  window: i32,
  step: i32,
  replicates: i32,
  cutoffPermille: i32,
  seed: u32,
  pairMapPtr: i32,
  pairListPtr: i32,
  weightPtr: i32,
  pairDistancePtr: i32,
  globalPairPtr: i32,
  statePtr: i32,
  differencePtr: i32,
  validPtr: i32,
  lookupPtr: i32,
  outPtr: i32,
  outCapacity: i32,
): i32 {
  return scanBatchCore(
    0, seqPtr, 0, 0, nSeq, nSites, tripletPtr, tripletCount,
    window, step, replicates, cutoffPermille, seed, pairMapPtr, pairListPtr, weightPtr,
    pairDistancePtr, globalPairPtr, statePtr, differencePtr, validPtr,
    lookupPtr, outPtr, outCapacity,
  );
}

export function scan_source_bootscan_batch_packed(
  packedPtr: i32,
  validityPtr: i32,
  wordsPerSequence: i32,
  nSeq: i32,
  nSites: i32,
  tripletPtr: i32,
  tripletCount: i32,
  window: i32,
  step: i32,
  replicates: i32,
  cutoffPermille: i32,
  seed: u32,
  pairMapPtr: i32,
  pairListPtr: i32,
  weightPtr: i32,
  pairDistancePtr: i32,
  globalPairPtr: i32,
  statePtr: i32,
  differencePtr: i32,
  validPtr: i32,
  lookupPtr: i32,
  outPtr: i32,
  outCapacity: i32,
): i32 {
  return scanBatchCore(
    1, packedPtr, validityPtr, wordsPerSequence, nSeq, nSites,
    tripletPtr, tripletCount, window, step, replicates, cutoffPermille, seed,
    pairMapPtr, pairListPtr, weightPtr, pairDistancePtr, globalPairPtr, statePtr,
    differencePtr, validPtr, lookupPtr, outPtr, outCapacity,
  );
}
