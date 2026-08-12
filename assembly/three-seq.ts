// RDP5 3SEQ triplet kernel.
//
// This is a direct, fused port of the author-supplied FindSubSeqTS walk used
// by TSXOver.  One concrete unordered triplet is decoded once.  During that
// pass all three possible recombinant assignments are compressed to the
// source +1/-1 walks:
//
//   +1  selected recombinant matches the first parent only
//   -1  selected recombinant matches the second parent only
//
// Invariant sites, incomplete sites, triallelic sites, and sites at which the
// selected recombinant matches neither parent never enter that assignment's
// walk.  The output contains the maximum descent and maximum ascent for each
// recombinant assignment; the worker applies the source probability
// recurrence and retains the lower-p orientation, exactly as TSXOver does.

const THREE_SEQ_ROW_INTS: i32 = 16;
const THREE_SEQ_WORKSPACE_WORDS_PER_SITE: i32 = 6;

// The source CheckwrapC pass needs the compressed coordinates and cumulative
// heights for each of the three target-role walks. Six n-site i32 lanes hold
// those arrays and five trailing words return the chosen interval and
// provenance flags.
export function source_three_seq_workspace_bytes(nSites: i32): i32 {
  return (nSites > 0 ? nSites : 0) * THREE_SEQ_WORKSPACE_WORDS_PER_SITE * 4 + 20;
}

@inline
function threeSeqBase(
  packedMode: bool,
  dataPtr: i32,
  validityPtr: i32,
  wordsPerSequence: i32,
  nSites: i32,
  sequence: i32,
  site: i32,
): u8 {
  if (!packedMode) return load<u8>(usize(dataPtr + sequence * nSites + site));
  const wordIndex = sequence * wordsPerSequence + (site >> 4);
  const shift = (site & 15) << 1;
  const validity = load<u32>(usize(validityPtr + wordIndex * 4));
  if (((validity >> shift) & 1) == 0) return 4;
  return u8((load<u32>(usize(dataPtr + wordIndex * 4)) >> shift) & 3);
}

function writeThreeSeqRow(
  outPtr: i32,
  row: i32,
  target: i32,
  majorParent: i32,
  minorParent: i32,
  direction: i32,
  start: i32,
  end: i32,
  upSteps: i32,
  downSteps: i32,
  descent: i32,
  informative: i32,
  cycle: i32,
  sourceWrap: i32,
  linearComplement: i32,
): void {
  const base = outPtr + row * THREE_SEQ_ROW_INTS * 4;
  store<i32>(usize(base), target);
  store<i32>(usize(base + 4), majorParent);
  store<i32>(usize(base + 8), minorParent);
  store<i32>(usize(base + 12), direction);
  store<i32>(usize(base + 16), start);
  store<i32>(usize(base + 20), end);
  store<i32>(usize(base + 24), upSteps);
  store<i32>(usize(base + 28), downSteps);
  store<i32>(usize(base + 32), descent);
  store<i32>(usize(base + 36), informative);
  store<i32>(usize(base + 40), cycle);
  store<i32>(usize(base + 44), sourceWrap);
  store<i32>(usize(base + 48), linearComplement);
  store<i32>(usize(base + 52), 0); // reserved: split-piece recalibration
  store<i32>(usize(base + 56), 0);
  store<i32>(usize(base + 60), 0);
}

@inline
function threeSeqPrefix(
  heightsPtr: i32,
  count: i32,
  boundary: i32,
  negMod: i32,
): i32 {
  if (boundary <= 0 || count <= 0) return 0;
  const total = load<i32>(usize(heightsPtr + (count - 1) * 4));
  const cycles = boundary / count;
  const rank = boundary % count;
  let value = cycles * total;
  if (rank > 0) value += load<i32>(usize(heightsPtr + (rank - 1) * 4));
  return value * negMod;
}

// Direct translation of CheckwrapC in compressed-rank coordinates. The source
// does not search every rotation: it starts from FindSubSeqTS's final global
// maximum (or minimum for the reverse walk), appends only the prefix through
// the already-known peak/trough boundary, and updates on strict improvements.
// For a linear alignment it reports the complementary non-wrapping interval
// when that extension crosses the origin, deliberately retaining nK.
function sourceThreeSeqExcursion(
  positionsPtr: i32,
  heightsPtr: i32,
  resultPtr: i32,
  count: i32,
  nSites: i32,
  negMod: i32,
  circular: bool,
  peakBoundaryInput: i32,
  endRankInput: i32,
  initialDescent: i32,
): void {
  let peakBoundary = peakBoundaryInput;
  let endRank = endRankInput;
  let best = initialDescent;
  let maximum = threeSeqPrefix(heightsPtr, count, peakBoundary, negMod);
  const total = threeSeqPrefix(heightsPtr, count, count, negMod);
  const peakRank = peakBoundary > 0 ? peakBoundary - 1 : 0;
  const limit = peakRank < endRank ? peakRank : endRank;

  // A source BE of zero is normalized to the first informative coordinate.
  // In boundary form it still represents the pre-walk height, so there is no
  // prefix before it to append. This avoids an artificial one-step increase.
  if (peakBoundary > 0) {
    for (let rank: i32 = 0; rank <= limit; rank += 1) {
      const height = total + load<i32>(usize(heightsPtr + rank * 4)) * negMod;
      if (height > maximum) {
        maximum = height;
        peakBoundary = rank + 1;
      }
      if (maximum - height > best) {
        best = maximum - height;
        endRank = rank;
      }
    }
  }

  let rawStart = load<i32>(usize(positionsPtr + (peakBoundary % count) * 4));
  let rawEnd = load<i32>(usize(positionsPtr + endRank * 4)) + 1;
  if (rawEnd > nSites) rawEnd = nSites;
  const crossedOrigin = rawStart > rawEnd;
  let linearComplement = 0;
  if (crossedOrigin && !circular) {
    const complementStart = (endRank + 1) % count;
    const complementEnd = (peakBoundary + count - 1) % count;
    rawStart = load<i32>(usize(positionsPtr + complementStart * 4));
    rawEnd = load<i32>(usize(positionsPtr + complementEnd * 4)) + 1;
    if (rawEnd > nSites) rawEnd = nSites;
    linearComplement = 1;
  }
  store<i32>(usize(resultPtr), rawStart);
  store<i32>(usize(resultPtr + 4), rawEnd);
  store<i32>(usize(resultPtr + 8), best);
  store<i32>(usize(resultPtr + 12), crossedOrigin ? 1 : 0);
  store<i32>(usize(resultPtr + 16), linearComplement);
}

function scanSourceThreeSeqTripletCore(
  dataPtr: i32,
  validityPtr: i32,
  wordsPerSequence: i32,
  packedMode: bool,
  nSites: i32,
  first: i32,
  second: i32,
  third: i32,
  workspacePtr: i32,
  sourceWrapMode: bool,
  circular: bool,
  outPtr: i32,
): i32 {
  const position0Ptr = workspacePtr;
  const position1Ptr = workspacePtr + nSites * 4;
  const position2Ptr = workspacePtr + nSites * 8;
  const height0Ptr = workspacePtr + nSites * 12;
  const height1Ptr = workspacePtr + nSites * 16;
  const height2Ptr = workspacePtr + nSites * 20;
  const excursionPtr = workspacePtr + nSites * 24;
  // Cycle 0 mirrors TSXOver(first, second, third): third is queried.
  let h0: i32 = 0, maximum0: i32 = 0, minimum0: i32 = 0;
  let descent0: i32 = 0, ascent0: i32 = 0;
  let maxStart0: i32 = 0, minStart0: i32 = 0;
  let descentStart0: i32 = 0, descentEnd0: i32 = 0;
  let ascentStart0: i32 = 0, ascentEnd0: i32 = 0;
  let maxBoundary0: i32 = 0, minBoundary0: i32 = 0;
  let descentEndRank0: i32 = 0, ascentEndRank0: i32 = 0;
  let maxPending0: bool = true, minPending0: bool = true;
  let up0: i32 = 0, down0: i32 = 0;

  // Cycle 1 mirrors TSXOver(second, third, first): first is queried.
  let h1: i32 = 0, maximum1: i32 = 0, minimum1: i32 = 0;
  let descent1: i32 = 0, ascent1: i32 = 0;
  let maxStart1: i32 = 0, minStart1: i32 = 0;
  let descentStart1: i32 = 0, descentEnd1: i32 = 0;
  let ascentStart1: i32 = 0, ascentEnd1: i32 = 0;
  let maxBoundary1: i32 = 0, minBoundary1: i32 = 0;
  let descentEndRank1: i32 = 0, ascentEndRank1: i32 = 0;
  let maxPending1: bool = true, minPending1: bool = true;
  let up1: i32 = 0, down1: i32 = 0;

  // Cycle 2 mirrors TSXOver(third, first, second): second is queried.
  let h2: i32 = 0, maximum2: i32 = 0, minimum2: i32 = 0;
  let descent2: i32 = 0, ascent2: i32 = 0;
  let maxStart2: i32 = 0, minStart2: i32 = 0;
  let descentStart2: i32 = 0, descentEnd2: i32 = 0;
  let ascentStart2: i32 = 0, ascentEnd2: i32 = 0;
  let maxBoundary2: i32 = 0, minBoundary2: i32 = 0;
  let descentEndRank2: i32 = 0, ascentEndRank2: i32 = 0;
  let maxPending2: bool = true, minPending2: bool = true;
  let up2: i32 = 0, down2: i32 = 0;

  for (let site: i32 = 0; site < nSites; site += 1) {
    const a = threeSeqBase(packedMode, dataPtr, validityPtr, wordsPerSequence, nSites, first, site);
    const b = threeSeqBase(packedMode, dataPtr, validityPtr, wordsPerSequence, nSites, second, site);
    const c = threeSeqBase(packedMode, dataPtr, validityPtr, wordsPerSequence, nSites, third, site);

    // A single biallelic classification yields all three target-role walks.
    // This is equivalent to the three FindSubSeqTS2 calls in TSXOverC, but it
    // avoids re-reading and re-comparing the same triplet column three times.
    let code0: i32 = 0;
    let code1: i32 = 0;
    let code2: i32 = 0;
    if (a < 4 && b < 4 && c < 4) {
      if (a == b) {
        if (a != c) { code1 = 1; code2 = -1; }
      } else if (a == c) {
        code0 = 1; code1 = -1;
      } else if (b == c) {
        code0 = -1; code2 = 1;
      }
    }

    if (code0 != 0) {
      const rank = up0 + down0;
      if (maxPending0) { maxStart0 = site; maxPending0 = false; }
      if (minPending0) { minStart0 = site; minPending0 = false; }
      if (code0 > 0) {
        h0 += 1; up0 += 1;
        if (h0 > maximum0) { maximum0 = h0; maxBoundary0 = rank + 1; maxPending0 = true; }
        if (h0 - minimum0 > ascent0) {
          ascent0 = h0 - minimum0; ascentStart0 = minStart0; ascentEnd0 = site + 1; ascentEndRank0 = rank;
        }
      } else {
        h0 -= 1; down0 += 1;
        if (maximum0 - h0 > descent0) {
          descent0 = maximum0 - h0; descentStart0 = maxStart0; descentEnd0 = site + 1; descentEndRank0 = rank;
        }
        if (h0 < minimum0) { minimum0 = h0; minBoundary0 = rank + 1; minPending0 = true; }
      }
      if (sourceWrapMode) {
        store<i32>(usize(position0Ptr + rank * 4), site);
        store<i32>(usize(height0Ptr + rank * 4), h0);
      }
    }

    if (code1 != 0) {
      const rank = up1 + down1;
      if (maxPending1) { maxStart1 = site; maxPending1 = false; }
      if (minPending1) { minStart1 = site; minPending1 = false; }
      if (code1 > 0) {
        h1 += 1; up1 += 1;
        if (h1 > maximum1) { maximum1 = h1; maxBoundary1 = rank + 1; maxPending1 = true; }
        if (h1 - minimum1 > ascent1) {
          ascent1 = h1 - minimum1; ascentStart1 = minStart1; ascentEnd1 = site + 1; ascentEndRank1 = rank;
        }
      } else {
        h1 -= 1; down1 += 1;
        if (maximum1 - h1 > descent1) {
          descent1 = maximum1 - h1; descentStart1 = maxStart1; descentEnd1 = site + 1; descentEndRank1 = rank;
        }
        if (h1 < minimum1) { minimum1 = h1; minBoundary1 = rank + 1; minPending1 = true; }
      }
      if (sourceWrapMode) {
        store<i32>(usize(position1Ptr + rank * 4), site);
        store<i32>(usize(height1Ptr + rank * 4), h1);
      }
    }

    if (code2 != 0) {
      const rank = up2 + down2;
      if (maxPending2) { maxStart2 = site; maxPending2 = false; }
      if (minPending2) { minStart2 = site; minPending2 = false; }
      if (code2 > 0) {
        h2 += 1; up2 += 1;
        if (h2 > maximum2) { maximum2 = h2; maxBoundary2 = rank + 1; maxPending2 = true; }
        if (h2 - minimum2 > ascent2) {
          ascent2 = h2 - minimum2; ascentStart2 = minStart2; ascentEnd2 = site + 1; ascentEndRank2 = rank;
        }
      } else {
        h2 -= 1; down2 += 1;
        if (maximum2 - h2 > descent2) {
          descent2 = maximum2 - h2; descentStart2 = maxStart2; descentEnd2 = site + 1; descentEndRank2 = rank;
        }
        if (h2 < minimum2) { minimum2 = h2; minBoundary2 = rank + 1; minPending2 = true; }
      }
      if (sourceWrapMode) {
        store<i32>(usize(position2Ptr + rank * 4), site);
        store<i32>(usize(height2Ptr + rank * 4), h2);
      }
    }
  }

  let rows: i32 = 0;
  if (sourceWrapMode) {
    const informative0 = up0 + down0;
    if (informative0 >= 4) {
      sourceThreeSeqExcursion(position0Ptr, height0Ptr, excursionPtr, informative0, nSites, 1, circular,
        maxBoundary0, descentEndRank0, descent0);
      const wrappedDescent = load<i32>(usize(excursionPtr + 8));
      if (wrappedDescent > 0) {
        writeThreeSeqRow(outPtr, rows++, third, first, second, 1,
          load<i32>(usize(excursionPtr)), load<i32>(usize(excursionPtr + 4)), up0, down0,
          wrappedDescent, informative0, 0, load<i32>(usize(excursionPtr + 12)), load<i32>(usize(excursionPtr + 16)));
      }
      sourceThreeSeqExcursion(position0Ptr, height0Ptr, excursionPtr, informative0, nSites, -1, circular,
        minBoundary0, ascentEndRank0, ascent0);
      const wrappedAscent = load<i32>(usize(excursionPtr + 8));
      if (wrappedAscent > 0) {
        writeThreeSeqRow(outPtr, rows++, third, second, first, -1,
          load<i32>(usize(excursionPtr)), load<i32>(usize(excursionPtr + 4)), down0, up0,
          wrappedAscent, informative0, 0, load<i32>(usize(excursionPtr + 12)), load<i32>(usize(excursionPtr + 16)));
      }
    }
    const informative1 = up1 + down1;
    if (informative1 >= 4) {
      sourceThreeSeqExcursion(position1Ptr, height1Ptr, excursionPtr, informative1, nSites, 1, circular,
        maxBoundary1, descentEndRank1, descent1);
      const wrappedDescent = load<i32>(usize(excursionPtr + 8));
      if (wrappedDescent > 0) {
        writeThreeSeqRow(outPtr, rows++, first, second, third, 1,
          load<i32>(usize(excursionPtr)), load<i32>(usize(excursionPtr + 4)), up1, down1,
          wrappedDescent, informative1, 1, load<i32>(usize(excursionPtr + 12)), load<i32>(usize(excursionPtr + 16)));
      }
      sourceThreeSeqExcursion(position1Ptr, height1Ptr, excursionPtr, informative1, nSites, -1, circular,
        minBoundary1, ascentEndRank1, ascent1);
      const wrappedAscent = load<i32>(usize(excursionPtr + 8));
      if (wrappedAscent > 0) {
        writeThreeSeqRow(outPtr, rows++, first, third, second, -1,
          load<i32>(usize(excursionPtr)), load<i32>(usize(excursionPtr + 4)), down1, up1,
          wrappedAscent, informative1, 1, load<i32>(usize(excursionPtr + 12)), load<i32>(usize(excursionPtr + 16)));
      }
    }
    const informative2 = up2 + down2;
    if (informative2 >= 4) {
      sourceThreeSeqExcursion(position2Ptr, height2Ptr, excursionPtr, informative2, nSites, 1, circular,
        maxBoundary2, descentEndRank2, descent2);
      const wrappedDescent = load<i32>(usize(excursionPtr + 8));
      if (wrappedDescent > 0) {
        writeThreeSeqRow(outPtr, rows++, second, third, first, 1,
          load<i32>(usize(excursionPtr)), load<i32>(usize(excursionPtr + 4)), up2, down2,
          wrappedDescent, informative2, 2, load<i32>(usize(excursionPtr + 12)), load<i32>(usize(excursionPtr + 16)));
      }
      sourceThreeSeqExcursion(position2Ptr, height2Ptr, excursionPtr, informative2, nSites, -1, circular,
        minBoundary2, ascentEndRank2, ascent2);
      const wrappedAscent = load<i32>(usize(excursionPtr + 8));
      if (wrappedAscent > 0) {
        writeThreeSeqRow(outPtr, rows++, second, first, third, -1,
          load<i32>(usize(excursionPtr)), load<i32>(usize(excursionPtr + 4)), down2, up2,
          wrappedAscent, informative2, 2, load<i32>(usize(excursionPtr + 12)), load<i32>(usize(excursionPtr + 16)));
      }
    }
    return rows;
  }
  if (up0 + down0 >= 4 && descent0 > 0) {
    writeThreeSeqRow(outPtr, rows++, third, first, second, 1, descentStart0, descentEnd0, up0, down0, descent0, up0 + down0, 0, 0, 0);
  }
  if (up0 + down0 >= 4 && ascent0 > 0) {
    writeThreeSeqRow(outPtr, rows++, third, second, first, -1, ascentStart0, ascentEnd0, down0, up0, ascent0, up0 + down0, 0, 0, 0);
  }
  if (up1 + down1 >= 4 && descent1 > 0) {
    writeThreeSeqRow(outPtr, rows++, first, second, third, 1, descentStart1, descentEnd1, up1, down1, descent1, up1 + down1, 1, 0, 0);
  }
  if (up1 + down1 >= 4 && ascent1 > 0) {
    writeThreeSeqRow(outPtr, rows++, first, third, second, -1, ascentStart1, ascentEnd1, down1, up1, ascent1, up1 + down1, 1, 0, 0);
  }
  if (up2 + down2 >= 4 && descent2 > 0) {
    writeThreeSeqRow(outPtr, rows++, second, third, first, 1, descentStart2, descentEnd2, up2, down2, descent2, up2 + down2, 2, 0, 0);
  }
  if (up2 + down2 >= 4 && ascent2 > 0) {
    writeThreeSeqRow(outPtr, rows++, second, first, third, -1, ascentStart2, ascentEnd2, down2, up2, ascent2, up2 + down2, 2, 0, 0);
  }
  return rows;
}

export function scan_source_three_seq_triplet(
  seqPtr: i32,
  nSites: i32,
  first: i32,
  second: i32,
  third: i32,
  outPtr: i32,
): i32 {
  return scanSourceThreeSeqTripletCore(seqPtr, 0, 0, false, nSites, first, second, third, 0, false, false, outPtr);
}

export function scan_source_three_seq_triplet_packed(
  packedPtr: i32,
  validityPtr: i32,
  wordsPerSequence: i32,
  nSites: i32,
  first: i32,
  second: i32,
  third: i32,
  outPtr: i32,
): i32 {
  return scanSourceThreeSeqTripletCore(packedPtr, validityPtr, wordsPerSequence, true, nSites, first, second, third, 0, false, false, outPtr);
}

// Full source path including CheckwrapC. circularFlag=0 preserves the desktop
// linear-complement branch; circularFlag!=0 retains origin-spanning intervals.
export function scan_source_three_seq_triplet_mode(
  seqPtr: i32,
  nSites: i32,
  first: i32,
  second: i32,
  third: i32,
  circularFlag: i32,
  workspacePtr: i32,
  outPtr: i32,
): i32 {
  return scanSourceThreeSeqTripletCore(
    seqPtr, 0, 0, false, nSites, first, second, third,
    workspacePtr, true, circularFlag != 0, outPtr,
  );
}

export function scan_source_three_seq_triplet_packed_mode(
  packedPtr: i32,
  validityPtr: i32,
  wordsPerSequence: i32,
  nSites: i32,
  first: i32,
  second: i32,
  third: i32,
  circularFlag: i32,
  workspacePtr: i32,
  outPtr: i32,
): i32 {
  return scanSourceThreeSeqTripletCore(
    packedPtr, validityPtr, wordsPerSequence, true, nSites, first, second, third,
    workspacePtr, true, circularFlag != 0, outPtr,
  );
}
