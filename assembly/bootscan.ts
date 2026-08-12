// RDP5 automated BOOTSCAN/RECSCAN batch.
//
// This follows the author-supplied BSXoverR2 -> SEQBOOT2 -> FastBootDistIP ->
// GetPltVal2 path. The unit of work is the active sequence set, not an
// oriented query/reference triplet: one bootstrap table is reused at every
// window, each required pair is evaluated once, then every concrete triplet's
// three pair relationships are interpreted together.  Relationship mode 0
// uses those distances directly.  Modes 1 and 2 construct one full-cohort
// UPGMA or neighbour-joining tree for each replicate and replace distances
// with the stored relative tree positions before any triplet is interpreted.

const BOOTSCAN_DISTANCE: i32 = 0;
const BOOTSCAN_UPGMA: i32 = 1;
const BOOTSCAN_NJ: i32 = 2;

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

@inline
function globalTreeCloser(globalPairPtr: i32, leftPair: i32, rightPair: i32): bool {
  return load<i32>(usize(globalPairPtr + leftPair * 8))
    < load<i32>(usize(globalPairPtr + rightPair * 8));
}

@inline
function globalClosestRelationship(
  relationshipMode: i32,
  globalPairPtr: i32,
  pair0: i32,
  pair1: i32,
  pair2: i32,
): i32 {
  if (relationshipMode == BOOTSCAN_DISTANCE) {
    return globalClosestPair(globalPairPtr, pair0, pair1, pair2);
  }
  let closest: i32 = 0;
  if (globalTreeCloser(globalPairPtr, pair1, pair0)) {
    closest = 1;
    if (globalTreeCloser(globalPairPtr, pair2, pair1)) closest = 2;
  } else if (globalTreeCloser(globalPairPtr, pair2, pair0)) {
    closest = 2;
  }
  return closest;
}

@inline
function alignWorkspace(value: i32, alignment: i32): i32 {
  return (value + alignment - 1) & ~(alignment - 1);
}

// One reusable tree is sufficient: RDP's bootstrap matrices are transformed
// serially, then their pair relationships are retained in pairDistancePtr.
export function source_bootscan_tree_workspace_bytes(nSeq: i32): i32 {
  if (nSeq < 3 || nSeq > 4096) return 0;
  const capacity = 2 * nSeq;
  let bytes = capacity * capacity * 8; // mutable f64 distance matrix
  bytes = alignWorkspace(bytes, 8) + capacity * 8; // NJ row sums
  bytes = alignWorkspace(bytes, 4) + capacity * 4; // active nodes
  bytes += capacity * 4; // UPGMA cluster sizes
  bytes += capacity * 3 * 4; // binary-tree adjacency
  bytes += capacity * 4; // node degrees
  bytes += capacity * 4; // BFS queue
  bytes += capacity * 4; // BFS distances
  return alignWorkspace(bytes, 16);
}

@inline
function treeMatrixPtr(treeWorkPtr: i32): i32 {
  return treeWorkPtr;
}

@inline
function treeRowSumPtr(treeWorkPtr: i32, capacity: i32): i32 {
  return alignWorkspace(treeWorkPtr + capacity * capacity * 8, 8);
}

@inline
function treeActivePtr(treeWorkPtr: i32, capacity: i32): i32 {
  return alignWorkspace(treeRowSumPtr(treeWorkPtr, capacity) + capacity * 8, 4);
}

@inline
function treeSizePtr(treeWorkPtr: i32, capacity: i32): i32 {
  return treeActivePtr(treeWorkPtr, capacity) + capacity * 4;
}

@inline
function treeNeighborPtr(treeWorkPtr: i32, capacity: i32): i32 {
  return treeSizePtr(treeWorkPtr, capacity) + capacity * 4;
}

@inline
function treeDegreePtr(treeWorkPtr: i32, capacity: i32): i32 {
  return treeNeighborPtr(treeWorkPtr, capacity) + capacity * 3 * 4;
}

@inline
function treeQueuePtr(treeWorkPtr: i32, capacity: i32): i32 {
  return treeDegreePtr(treeWorkPtr, capacity) + capacity * 4;
}

@inline
function treeGraphDistancePtr(treeWorkPtr: i32, capacity: i32): i32 {
  return treeQueuePtr(treeWorkPtr, capacity) + capacity * 4;
}

@inline
function treeMatrixAddress(matrixPtr: i32, capacity: i32, left: i32, right: i32): i32 {
  return matrixPtr + (left * capacity + right) * 8;
}

@inline
function treeDistance(matrixPtr: i32, capacity: i32, left: i32, right: i32): f64 {
  return load<f64>(usize(treeMatrixAddress(matrixPtr, capacity, left, right)));
}

@inline
function setTreeDistance(
  matrixPtr: i32,
  capacity: i32,
  left: i32,
  right: i32,
  value: f64,
): void {
  store<f64>(usize(treeMatrixAddress(matrixPtr, capacity, left, right)), value);
  store<f64>(usize(treeMatrixAddress(matrixPtr, capacity, right, left)), value);
}

@inline
function addTreeEdge(neighborPtr: i32, degreePtr: i32, left: i32, right: i32): void {
  const leftDegreeAddress = degreePtr + left * 4;
  const rightDegreeAddress = degreePtr + right * 4;
  const leftDegree = load<i32>(usize(leftDegreeAddress));
  const rightDegree = load<i32>(usize(rightDegreeAddress));
  if (leftDegree < 3) {
    store<i32>(usize(neighborPtr + (left * 3 + leftDegree) * 4), right);
    store<i32>(usize(leftDegreeAddress), leftDegree + 1);
  }
  if (rightDegree < 3) {
    store<i32>(usize(neighborPtr + (right * 3 + rightDegree) * 4), left);
    store<i32>(usize(rightDegreeAddress), rightDegree + 1);
  }
}

function replaceActivePair(
  activePtr: i32,
  activeCount: i32,
  firstPosition: i32,
  secondPosition: i32,
  mergedNode: i32,
): i32 {
  let write: i32 = 0;
  for (let position: i32 = 0; position < activeCount; position += 1) {
    if (position == firstPosition || position == secondPosition) continue;
    store<i32>(usize(activePtr + write * 4), load<i32>(usize(activePtr + position * 4)));
    write += 1;
  }
  store<i32>(usize(activePtr + write * 4), mergedNode);
  return write + 1;
}

function buildUpgmaTopology(
  nSeq: i32,
  matrixPtr: i32,
  activePtr: i32,
  sizePtr: i32,
  neighborPtr: i32,
  degreePtr: i32,
): i32 {
  let activeCount = nSeq;
  let nextNode = nSeq;
  while (activeCount > 1) {
    let bestFirst: i32 = 0;
    let bestSecond: i32 = 1;
    let bestDistance = treeDistance(
      matrixPtr, 2 * nSeq,
      load<i32>(usize(activePtr)),
      load<i32>(usize(activePtr + 4)),
    );
    for (let firstPosition: i32 = 0; firstPosition < activeCount - 1; firstPosition += 1) {
      const first = load<i32>(usize(activePtr + firstPosition * 4));
      for (let secondPosition = firstPosition + 1; secondPosition < activeCount; secondPosition += 1) {
        const second = load<i32>(usize(activePtr + secondPosition * 4));
        const distance = treeDistance(matrixPtr, 2 * nSeq, first, second);
        if (distance < bestDistance) {
          bestDistance = distance;
          bestFirst = firstPosition;
          bestSecond = secondPosition;
        }
      }
    }
    const left = load<i32>(usize(activePtr + bestFirst * 4));
    const right = load<i32>(usize(activePtr + bestSecond * 4));
    const leftSize = load<i32>(usize(sizePtr + left * 4));
    const rightSize = load<i32>(usize(sizePtr + right * 4));
    const combinedSize = leftSize + rightSize;
    for (let position: i32 = 0; position < activeCount; position += 1) {
      if (position == bestFirst || position == bestSecond) continue;
      const other = load<i32>(usize(activePtr + position * 4));
      const mergedDistance = (
        treeDistance(matrixPtr, 2 * nSeq, left, other) * f64(leftSize)
        + treeDistance(matrixPtr, 2 * nSeq, right, other) * f64(rightSize)
      ) / f64(combinedSize);
      setTreeDistance(matrixPtr, 2 * nSeq, nextNode, other, mergedDistance);
    }
    store<i32>(usize(sizePtr + nextNode * 4), combinedSize);
    addTreeEdge(neighborPtr, degreePtr, nextNode, left);
    addTreeEdge(neighborPtr, degreePtr, nextNode, right);
    activeCount = replaceActivePair(activePtr, activeCount, bestFirst, bestSecond, nextNode);
    nextNode += 1;
  }
  return nextNode;
}

function buildNeighborJoiningTopology(
  nSeq: i32,
  matrixPtr: i32,
  rowSumPtr: i32,
  activePtr: i32,
  neighborPtr: i32,
  degreePtr: i32,
): i32 {
  let activeCount = nSeq;
  let nextNode = nSeq;
  while (activeCount > 2) {
    for (let position: i32 = 0; position < activeCount; position += 1) {
      const node = load<i32>(usize(activePtr + position * 4));
      let sum: f64 = 0.0;
      for (let otherPosition: i32 = 0; otherPosition < activeCount; otherPosition += 1) {
        if (otherPosition == position) continue;
        const other = load<i32>(usize(activePtr + otherPosition * 4));
        sum += treeDistance(matrixPtr, 2 * nSeq, node, other);
      }
      store<f64>(usize(rowSumPtr + position * 8), sum);
    }
    let bestFirst: i32 = 0;
    let bestSecond: i32 = 1;
    const initialLeft = load<i32>(usize(activePtr));
    const initialRight = load<i32>(usize(activePtr + 4));
    let bestQ = f64(activeCount - 2) * treeDistance(matrixPtr, 2 * nSeq, initialLeft, initialRight)
      - load<f64>(usize(rowSumPtr)) - load<f64>(usize(rowSumPtr + 8));
    for (let firstPosition: i32 = 0; firstPosition < activeCount - 1; firstPosition += 1) {
      const first = load<i32>(usize(activePtr + firstPosition * 4));
      const firstSum = load<f64>(usize(rowSumPtr + firstPosition * 8));
      for (let secondPosition = firstPosition + 1; secondPosition < activeCount; secondPosition += 1) {
        const second = load<i32>(usize(activePtr + secondPosition * 4));
        const q = f64(activeCount - 2) * treeDistance(matrixPtr, 2 * nSeq, first, second)
          - firstSum - load<f64>(usize(rowSumPtr + secondPosition * 8));
        if (q < bestQ) {
          bestQ = q;
          bestFirst = firstPosition;
          bestSecond = secondPosition;
        }
      }
    }
    const left = load<i32>(usize(activePtr + bestFirst * 4));
    const right = load<i32>(usize(activePtr + bestSecond * 4));
    const joinedDistance = treeDistance(matrixPtr, 2 * nSeq, left, right);
    for (let position: i32 = 0; position < activeCount; position += 1) {
      if (position == bestFirst || position == bestSecond) continue;
      const other = load<i32>(usize(activePtr + position * 4));
      let mergedDistance = 0.5 * (
        treeDistance(matrixPtr, 2 * nSeq, left, other)
        + treeDistance(matrixPtr, 2 * nSeq, right, other)
        - joinedDistance
      );
      if (mergedDistance < 0.0) mergedDistance = 0.0;
      setTreeDistance(matrixPtr, 2 * nSeq, nextNode, other, mergedDistance);
    }
    addTreeEdge(neighborPtr, degreePtr, nextNode, left);
    addTreeEdge(neighborPtr, degreePtr, nextNode, right);
    activeCount = replaceActivePair(activePtr, activeCount, bestFirst, bestSecond, nextNode);
    nextNode += 1;
  }
  addTreeEdge(
    neighborPtr,
    degreePtr,
    load<i32>(usize(activePtr)),
    load<i32>(usize(activePtr + 4)),
  );
  return nextNode;
}

function storeTreePathRelationships(
  nSeq: i32,
  nodeCount: i32,
  replicates: i32,
  replicate: i32,
  pairMapPtr: i32,
  pairDistancePtr: i32,
  neighborPtr: i32,
  degreePtr: i32,
  queuePtr: i32,
  graphDistancePtr: i32,
): void {
  for (let source: i32 = 0; source < nSeq - 1; source += 1) {
    memory.fill(graphDistancePtr, 0xff, nodeCount * 4);
    let head: i32 = 0;
    let tail: i32 = 1;
    store<i32>(usize(queuePtr), source);
    store<i32>(usize(graphDistancePtr + source * 4), 0);
    while (head < tail) {
      const node = load<i32>(usize(queuePtr + head * 4));
      head += 1;
      const nextDistance = load<i32>(usize(graphDistancePtr + node * 4)) + 1;
      const degree = load<i32>(usize(degreePtr + node * 4));
      for (let edge: i32 = 0; edge < degree; edge += 1) {
        const neighbor = load<i32>(usize(neighborPtr + (node * 3 + edge) * 4));
        if (load<i32>(usize(graphDistancePtr + neighbor * 4)) >= 0) continue;
        store<i32>(usize(graphDistancePtr + neighbor * 4), nextDistance);
        store<i32>(usize(queuePtr + tail * 4), neighbor);
        tail += 1;
      }
    }
    for (let target = source + 1; target < nSeq; target += 1) {
      const usedPair = load<i32>(usize(pairMapPtr + pairIndex(source, target, nSeq) * 4));
      if (usedPair < 0) continue;
      let path = load<i32>(usize(graphDistancePtr + target * 4));
      if (path < 0) path = 65535;
      if (path > 65535) path = 65535;
      store<u16>(usize(pairDistancePtr + (usedPair * replicates + replicate) * 2), u16(path));
    }
  }
}

function transformTreeRelationships(
  relationshipMode: i32,
  nSeq: i32,
  replicates: i32,
  replicate: i32,
  pairMapPtr: i32,
  pairDistancePtr: i32,
  treeWorkPtr: i32,
): void {
  const capacity = 2 * nSeq;
  const matrixPtr = treeMatrixPtr(treeWorkPtr);
  const rowSumPtr = treeRowSumPtr(treeWorkPtr, capacity);
  const activePtr = treeActivePtr(treeWorkPtr, capacity);
  const sizePtr = treeSizePtr(treeWorkPtr, capacity);
  const neighborPtr = treeNeighborPtr(treeWorkPtr, capacity);
  const degreePtr = treeDegreePtr(treeWorkPtr, capacity);
  const queuePtr = treeQueuePtr(treeWorkPtr, capacity);
  const graphDistancePtr = treeGraphDistancePtr(treeWorkPtr, capacity);

  memory.fill(degreePtr, 0, capacity * 4);
  for (let sequence: i32 = 0; sequence < nSeq; sequence += 1) {
    store<i32>(usize(activePtr + sequence * 4), sequence);
    store<i32>(usize(sizePtr + sequence * 4), 1);
    store<f64>(usize(treeMatrixAddress(matrixPtr, capacity, sequence, sequence)), 0.0);
  }
  for (let first: i32 = 0; first < nSeq - 1; first += 1) {
    for (let second = first + 1; second < nSeq; second += 1) {
      const usedPair = load<i32>(usize(pairMapPtr + pairIndex(first, second, nSeq) * 4));
      const distance = usedPair >= 0
        ? f64(load<u16>(usize(pairDistancePtr + (usedPair * replicates + replicate) * 2)))
        : 32000.0;
      setTreeDistance(matrixPtr, capacity, first, second, distance);
    }
  }
  const nodeCount = relationshipMode == BOOTSCAN_UPGMA
    ? buildUpgmaTopology(nSeq, matrixPtr, activePtr, sizePtr, neighborPtr, degreePtr)
    : buildNeighborJoiningTopology(nSeq, matrixPtr, rowSumPtr, activePtr, neighborPtr, degreePtr);
  storeTreePathRelationships(
    nSeq, nodeCount, replicates, replicate, pairMapPtr, pairDistancePtr,
    neighborPtr, degreePtr, queuePtr, graphDistancePtr,
  );
}

// Differential-test hook for the source tree-position transform. pairMapPtr
// must contain the dense triangular pair mapping and pairDistancePtr the
// quantized distance rows used by the production batch.
export function source_bootscan_transform_tree_relationships(
  relationshipMode: i32,
  nSeq: i32,
  replicates: i32,
  replicate: i32,
  pairMapPtr: i32,
  pairDistancePtr: i32,
  treeWorkPtr: i32,
): void {
  if (nSeq < 3 || replicates < 1 || replicate < 0 || replicate >= replicates) return;
  if (relationshipMode != BOOTSCAN_UPGMA && relationshipMode != BOOTSCAN_NJ) return;
  if (treeWorkPtr <= 0 || source_bootscan_tree_workspace_bytes(nSeq) <= 0) return;
  transformTreeRelationships(
    relationshipMode, nSeq, replicates, replicate,
    pairMapPtr, pairDistancePtr, treeWorkPtr,
  );
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

function accumulatePackedCompletePair(
  packedPtr: i32,
  validityPtr: i32,
  wordsPerSequence: i32,
  nSites: i32,
  first: i32,
  second: i32,
  rawWindowStart: i32,
  window: i32,
  replicates: i32,
  weightPtr: i32,
  differencePtr: i32,
): bool {
  let relative: i32 = 0;
  while (relative < window) {
    let site = rawWindowStart + relative;
    while (site < 0) site += nSites;
    while (site >= nSites) site -= nSites;
    const lane = site & 15;
    let chunk = 16 - lane;
    const sitesToOrigin = nSites - site;
    if (chunk > sitesToOrigin) chunk = sitesToOrigin;
    if (chunk > window - relative) chunk = window - relative;

    const firstWord = first * wordsPerSequence + (site >> 4);
    const secondWord = second * wordsPerSequence + (site >> 4);
    const laneShift = lane << 1;
    let evenMask: u32;
    if (chunk == 16) {
      evenMask = 0x55555555;
    } else {
      evenMask = ((u32(1) << (chunk << 1)) - 1) & 0x55555555;
      evenMask <<= laneShift;
    }
    const firstValid = load<u32>(usize(validityPtr + firstWord * 4));
    const secondValid = load<u32>(usize(validityPtr + secondWord * 4));
    if ((firstValid & secondValid & evenMask) != evenMask) return false;

    const xor = load<u32>(usize(packedPtr + firstWord * 4))
      ^ load<u32>(usize(packedPtr + secondWord * 4));
    const mismatchLanes = (xor | (xor >> 1)) & evenMask;
    if (mismatchLanes != 0) {
      for (let offset: i32 = 0; offset < chunk; offset += 1) {
        const laneBit = u32(1) << (laneShift + (offset << 1));
        if ((mismatchLanes & laneBit) == 0) continue;
        const weights = weightPtr + (relative + offset) * replicates * 2;
        for (let replicate: i32 = 0; replicate < replicates; replicate += 1) {
          const weight = i32(load<u16>(usize(weights + replicate * 2)));
          if (weight == 0) continue;
          const differenceAddress = differencePtr + replicate * 4;
          store<i32>(
            usize(differenceAddress),
            load<i32>(usize(differenceAddress)) + weight,
          );
        }
      }
    }
    relative += chunk;
  }
  return true;
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
  relationshipModeInput: i32,
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
  treeWorkPtr: i32,
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
  const relationshipMode = bounded(relationshipModeInput, BOOTSCAN_DISTANCE, BOOTSCAN_NJ);
  if (relationshipMode != BOOTSCAN_DISTANCE
    && (treeWorkPtr <= 0 || source_bootscan_tree_workspace_bytes(nSeq) <= 0)) return 0;
  const pairCount = (nSeq * (nSeq - 1)) / 2;
  memory.fill(pairMapPtr, 0xff, pairCount * 4);
  if (relationshipMode == BOOTSCAN_DISTANCE) {
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
  } else {
    // Tree relationships are cohort properties.  Even a shortlisted triplet
    // therefore needs every pair in the active alignment when its UPGMA/NJ
    // tree is constructed.
    memory.fill(pairMapPtr, 0, pairCount * 4);
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
  if (relationshipMode != BOOTSCAN_DISTANCE) {
    for (let usedPair: i32 = 0; usedPair < usedPairCount; usedPair += 1) {
      const validSites = load<i32>(usize(globalPairPtr + usedPair * 8));
      const differences = load<i32>(usize(globalPairPtr + usedPair * 8 + 4));
      store<u16>(
        usize(pairDistancePtr + usedPair * replicates * 2),
        quantizedDistance(validSites, differences),
      );
    }
    transformTreeRelationships(
      relationshipMode, nSeq, replicates, 0, pairMapPtr, pairDistancePtr, treeWorkPtr,
    );
    for (let usedPair: i32 = 0; usedPair < usedPairCount; usedPair += 1) {
      store<i32>(
        usize(globalPairPtr + usedPair * 8),
        i32(load<u16>(usize(pairDistancePtr + usedPair * replicates * 2))),
      );
      store<i32>(usize(globalPairPtr + usedPair * 8 + 4), 0);
    }
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
      // denominator, so visit the replicate table only at mismatches. Packed
      // scans reject/compare 16 sites per word; the byte path remains the
      // scalar oracle. If a missing base is encountered, restart this pair on
      // the exact fallback.
      if (packedMode != 0) {
        completePair = accumulatePackedCompletePair(
          dataPtr, validityPtr, wordsPerSequence, nSites, first, second,
          rawWindowStart, window, replicates, weightPtr, differencePtr,
        );
      } else {
        for (let relative: i32 = 0; relative < window; relative += 1) {
          let site = rawWindowStart + relative;
          while (site < 0) site += nSites;
          while (site >= nSites) site -= nSites;
          const a = byteBase(dataPtr, nSites, first, site);
          const b = byteBase(dataPtr, nSites, second, site);
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

    if (relationshipMode != BOOTSCAN_DISTANCE) {
      for (let replicate: i32 = 0; replicate < replicates; replicate += 1) {
        transformTreeRelationships(
          relationshipMode, nSeq, replicates, replicate,
          pairMapPtr, pairDistancePtr, treeWorkPtr,
        );
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
      const globalClosest = globalClosestRelationship(
        relationshipMode, globalPairPtr, pair0, pair1, pair2,
      );

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
    const globalClosest = globalClosestRelationship(
      relationshipMode, globalPairPtr, pair0, pair1, pair2,
    );
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
    window, step, replicates, cutoffPermille, BOOTSCAN_DISTANCE, seed,
    pairMapPtr, pairListPtr, weightPtr,
    pairDistancePtr, globalPairPtr, statePtr, differencePtr, validPtr,
    lookupPtr, 0, outPtr, outCapacity,
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
    tripletPtr, tripletCount, window, step, replicates, cutoffPermille,
    BOOTSCAN_DISTANCE, seed,
    pairMapPtr, pairListPtr, weightPtr, pairDistancePtr, globalPairPtr, statePtr,
    differencePtr, validPtr, lookupPtr, 0, outPtr, outCapacity,
  );
}

export function scan_source_bootscan_batch_mode(
  seqPtr: i32,
  nSeq: i32,
  nSites: i32,
  tripletPtr: i32,
  tripletCount: i32,
  window: i32,
  step: i32,
  replicates: i32,
  cutoffPermille: i32,
  relationshipMode: i32,
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
  treeWorkPtr: i32,
  outPtr: i32,
  outCapacity: i32,
): i32 {
  return scanBatchCore(
    0, seqPtr, 0, 0, nSeq, nSites, tripletPtr, tripletCount,
    window, step, replicates, cutoffPermille, relationshipMode, seed,
    pairMapPtr, pairListPtr, weightPtr, pairDistancePtr, globalPairPtr, statePtr,
    differencePtr, validPtr, lookupPtr, treeWorkPtr, outPtr, outCapacity,
  );
}

export function scan_source_bootscan_batch_mode_packed(
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
  relationshipMode: i32,
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
  treeWorkPtr: i32,
  outPtr: i32,
  outCapacity: i32,
): i32 {
  return scanBatchCore(
    1, packedPtr, validityPtr, wordsPerSequence, nSeq, nSites,
    tripletPtr, tripletCount, window, step, replicates, cutoffPermille,
    relationshipMode, seed, pairMapPtr, pairListPtr, weightPtr, pairDistancePtr,
    globalPairPtr, statePtr, differencePtr, validPtr, lookupPtr, treeWorkPtr,
    outPtr, outCapacity,
  );
}
