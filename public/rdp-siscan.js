// Source-guided Sister-Scanning implementation.
//
// The author-supplied RDP5 sources split this workflow across SetUpSiScan,
// Get3Score/GetPScores2, DoPerms3P, MakeZValue2, DoSums, FindMaxZ,
// ShrinkRegionC and GetSSOL.  Keeping those stages visible here makes the
// historical behavior testable without hiding it behind a generic Z-score.

const PATTERN_TABLE = [
  [0, 0, 0, 0, 0],
  [0, 4, 6, 7, 1],
  [0, 12, 2, 8, 2],
  [0, 13, 9, 3, 3],
  [0, 10, 14, 5, 5],
  [0, 15, 11, 11, 11],
];

const PATTERN_TOPOLOGY_GROUPS = [
  [2, 3, 5],
  [8, 9, 10],
];

const SUM_TOPOLOGY_GROUPS = [
  [1, 2, 3],
  [4, 5, 7],
];

function canonical(base, gapsAsState) {
  return base < 4 || (gapsAsState && base === 5);
}

function tripletState(first, second, third, gapsAsState) {
  if (!canonical(first, gapsAsState) || !canonical(second, gapsAsState) || !canonical(third, gapsAsState)) return 0;
  if (first === second) return first === third ? 5 : 2;
  if (first === third) return 3;
  if (second === third) return 4;
  return 1;
}

export function sourceSiScanPattern(first, second, third, outgroup, options = {}) {
  const gapsAsState = options.gapsAsState === true;
  const state = tripletState(first, second, third, gapsAsState);
  if (!state || !canonical(outgroup, gapsAsState)) return 0;
  let relation = 4;
  if (first === outgroup) relation = 1;
  else if (state === 1) {
    if (second === outgroup) relation = 2;
    else if (third === outgroup) relation = 3;
  } else if (state === 2) {
    if (third === outgroup) relation = 3;
  } else if (state === 3 || state === 4) {
    if (second === outgroup) relation = 2;
  }
  return PATTERN_TABLE[state][relation];
}

function includePattern(category, positionMode) {
  if (!category) return false;
  if (positionMode === "triplet-variable") return category !== 11 && category !== 15;
  if (positionMode === "quartet-variable") return category !== 15;
  return true;
}

// MSVC's 15-bit rand() is the generator used by the supplied DNA.dll source.
// MakeVRand consumes one value immediately after srand(), then fills rows
// 0..N and positions 0..L inclusive.  The otherwise-unused row and position
// zero matter because DoPerms3P addresses the template as a flat byte array.
export function buildSourceSiScanRandomization(length, permutations, seed = 0x5a17c0de) {
  const sites = Math.max(1, Math.trunc(length));
  const replicates = Math.max(2, Math.min(10_000, Math.trunc(permutations)));
  const stride = sites + 1;
  const totalValues = (replicates + 1) * stride;
  // The desktop DLL materialises this entire table.  At bacterial-genome
  // scale that can require several gigabytes even though DoPerms3P walks only
  // a prefix.  Cache modest tables and regenerate the identical MSVC rand()
  // stream in-place for large alignments.  This is bit-for-bit equivalent for
  // every template value actually consumed and bounds peak memory.
  const values = totalValues <= 8_000_000 ? new Uint8Array(totalValues) : null;
  let state = seed >>> 0;
  const next = () => {
    state = (Math.imul(state, 214013) + 2531011) >>> 0;
    return (state >>> 16) & 0x7fff;
  };
  next();
  if (values) {
    for (let replicate = 0; replicate <= replicates; replicate += 1) {
      const offset = replicate * stride;
      for (let site = 0; site <= sites; site += 1) {
        values[offset + site] = Math.floor((next() / 32767) * 11.999) + 1;
      }
    }
  }
  return { values, stride, totalValues, permutations: replicates, seed: seed >>> 0 };
}

function verticalCategory(category, randomClass) {
  if (category >= 2 && category <= 7) return 2 + ((randomClass - 1) % 6);
  if (category >= 8 && category <= 10) return 8 + ((randomClass - 1) % 3);
  if (category >= 11 && category <= 14) return 11 + ((randomClass - 1) % 4);
  return category;
}

function permutationCounts(actual, permutations, randomization) {
  const replicates = Math.max(2, Math.min(permutations, randomization.permutations));
  const counts = new Int32Array((replicates + 1) * 16);
  counts.set(actual.subarray(0, 16));
  let cursor = 0;
  let streamState = randomization.seed >>> 0;
  const nextStreamClass = () => {
    streamState = (Math.imul(streamState, 214013) + 2531011) >>> 0;
    return Math.floor((((streamState >>> 16) & 0x7fff) / 32767) * 11.999) + 1;
  };
  nextStreamClass(); // MakeVRand's post-srand discard.
  // DoPerms3P deliberately permutes only categories 2..14.  Its cursor is
  // advanced by the replicate count for every observed column, producing the
  // unusual flat-template addressing retained below.
  for (let category = 2; category <= 14; category += 1) {
    for (let occurrence = 0; occurrence < actual[category]; occurrence += 1) {
      const base = cursor - 1;
      for (let replicate = 1; replicate <= replicates; replicate += 1) {
        const templateIndex = base + replicate;
        const randomClass = randomization.values?.[templateIndex] ?? nextStreamClass();
        const mapped = verticalCategory(category, randomClass);
        counts[replicate * 16 + mapped] += 1;
      }
      cursor += replicates;
    }
  }
  return { counts, replicates };
}

function zScores(actual, permutations, randomization) {
  const { counts, replicates } = permutationCounts(actual, permutations, randomization);
  const patterns = new Float64Array(16);
  for (let category = 1; category <= 15; category += 1) {
    const observed = counts[category];
    if (observed <= 0) continue;
    let total = 0;
    let totalSquares = 0;
    for (let replicate = 1; replicate <= replicates; replicate += 1) {
      const value = counts[replicate * 16 + category];
      total += value;
      totalSquares += value * value;
    }
    const mean = total / replicates;
    // Automated RDP5 uses MakeZValue2's population variance here.
    const variance = Math.max(0, totalSquares / replicates - mean * mean);
    if (variance > 0) patterns[category] = (observed - mean) / Math.sqrt(variance);
  }

  const sums = new Float64Array(13);
  const sumMembers = [
    [],
    [2, 7, 8],
    [3, 6, 9],
    [4, 5, 10],
    [2, 8, 11, 12],
    [3, 9, 11, 13],
    [],
    [5, 10, 11, 14],
  ];
  for (const group of [1, 2, 3, 4, 5, 7]) {
    const members = sumMembers[group];
    const observed = members.reduce((total, category) => total + counts[category], 0);
    if (observed <= 0) continue;
    let total = 0;
    let totalSquares = 0;
    for (let replicate = 1; replicate <= replicates; replicate += 1) {
      const value = members.reduce((sum, category) => sum + counts[replicate * 16 + category], 0);
      total += value;
      totalSquares += value * value;
    }
    const mean = total / replicates;
    const variance = Math.max(0, totalSquares / replicates - mean * mean);
    if (variance > 0) sums[group] = (observed - mean) / Math.sqrt(variance);
  }
  return { patterns, sums, replicates };
}

// Window scanning calls the same vertical randomization thousands of times.
// DoPerms3P's flat template has a useful invariant: for replicate r, the
// random classes used by occurrence range [a,b) are template cells
// a*P+r..b*P+r.  Prefixing those class counts turns each window from O(P*W)
// random assignments into O(P*61) range lookups without changing one draw.
function buildPermutationPrefix(maxOccurrences, permutations, randomization) {
  const occurrences = Math.max(1, Math.trunc(maxOccurrences));
  const replicates = Math.max(2, Math.min(permutations, randomization.permutations));
  const classes = new Uint8Array(occurrences * replicates);
  if (randomization.values) {
    classes.set(randomization.values.subarray(0, classes.length));
  } else {
    let state = randomization.seed >>> 0;
    const next = () => {
      state = (Math.imul(state, 214013) + 2531011) >>> 0;
      return Math.floor(((((state >>> 16) & 0x7fff) / 32767) * 11.999)) + 1;
    };
    next();
    for (let index = 0; index < classes.length; index += 1) classes[index] = next();
  }
  const Counter = occurrences <= 65_535 ? Uint16Array : Uint32Array;
  const stride = occurrences + 1;
  const prefix6 = new Counter(replicates * stride * 6);
  const prefix3 = new Counter(replicates * stride * 3);
  const prefix4 = new Counter(replicates * stride * 4);
  const fill = (prefix, modulus, replicate, occurrence, randomClass) => {
    const previous = (replicate * stride + occurrence) * modulus;
    const current = previous + modulus;
    for (let category = 0; category < modulus; category += 1) prefix[current + category] = prefix[previous + category];
    prefix[current + ((randomClass - 1) % modulus)] += 1;
  };
  for (let replicate = 0; replicate < replicates; replicate += 1) {
    for (let occurrence = 0; occurrence < occurrences; occurrence += 1) {
      const randomClass = classes[occurrence * replicates + replicate] || 1;
      fill(prefix6, 6, replicate, occurrence, randomClass);
      fill(prefix3, 3, replicate, occurrence, randomClass);
      fill(prefix4, 4, replicate, occurrence, randomClass);
    }
  }
  return { replicates, occurrences, stride, prefix6, prefix3, prefix4 };
}

function reusableRegionPermutationPrefix(actual, permutations, randomization) {
  let occurrences = 0;
  for (let category = 2; category <= 14; category += 1) occurrences += actual[category];
  const replicates = Math.max(2, Math.min(permutations, randomization.permutations));
  const cache = randomization.preparedRegionCache
    ?? (randomization.preparedRegionCache = new Map());
  const existing = cache.get(replicates);
  if (existing?.occurrences >= occurrences) return existing;

  // Region lengths vary from one inferred run to the next. Grow in 256-site
  // quanta (and by at least 1.5× after the first allocation) so later, shorter
  // regions reuse the same exact template prefix without power-of-two waste.
  // The prefix is an accelerator only: very large tracts retain the bounded-
  // memory direct DoPerms3P stream instead of allocating hundreds of MiB.
  const quantum = 256;
  let capacity = Math.max(quantum, Math.ceil(Math.max(1, occurrences) / quantum) * quantum);
  if (existing) {
    capacity = Math.max(capacity, Math.ceil((existing.occurrences * 1.5) / quantum) * quantum);
  }
  const counterBytes = capacity <= 65_535 ? 2 : 4;
  const prefixBytes = replicates * (capacity + 1) * 13 * counterBytes;
  const classBytes = replicates * capacity;
  if (prefixBytes + classBytes > 72 * 1024 * 1024) return null;
  const prepared = buildPermutationPrefix(capacity, replicates, randomization);
  prepared.momentCache = existing?.momentCache ?? new Map();
  cache.set(replicates, prepared);
  return prepared;
}

function prefixedZScores(actual, prepared) {
  // DoPerms3P discards the identity of the source category inside each of its
  // three vertical-randomization bands: 2..7 all map to 2..7, 8..10 all map
  // to 8..10, and 11..14 all map to 11..14.  The old exact prefix path added
  // six adjacent ranges for the first band, three for the second and four for
  // the third.  Adjacent prefix differences telescope, so each destination
  // count is exactly one range query over the band's total occurrence count.
  // This preserves every source template draw and accumulation order while
  // cutting the hot loop from 61 to 13 prefix queries per replicate.
  const end6 = actual[2] + actual[3] + actual[4] + actual[5] + actual[6] + actual[7];
  const end3 = end6 + actual[8] + actual[9] + actual[10];
  const end4 = end3 + actual[11] + actual[12] + actual[13] + actual[14];
  // For a fixed three-band total, the exact randomized mean/variance is also
  // fixed; only the observed category vector changes.  Cache those moments by
  // the three cumulative endpoints so neighbouring windows and other concrete
  // triplets do not replay the same source permutation table.
  const momentKey = `${end6}:${end3}:${end4}`;
  const momentCache = prepared.momentCache ?? (prepared.momentCache = new Map());
  let moments = momentCache.get(momentKey);
  if (!moments) {
    const totals = new Float64Array(16);
    const squares = new Float64Array(16);
    const sumTotals = new Float64Array(8);
    const sumSquares = new Float64Array(8);
    const counts = new Int32Array(16);
    for (let replicate = 0; replicate < prepared.replicates; replicate += 1) {
      const base6 = replicate * prepared.stride * 6;
      const base3 = replicate * prepared.stride * 3;
      const base4 = replicate * prepared.stride * 4;
      const at6 = base6 + end6 * 6;
      const from3 = base3 + end6 * 3;
      const at3 = base3 + end3 * 3;
      const from4 = base4 + end3 * 4;
      const at4 = base4 + end4 * 4;
      for (let category = 0; category < 6; category += 1) {
        counts[2 + category] = prepared.prefix6[at6 + category];
      }
      for (let category = 0; category < 3; category += 1) {
        counts[8 + category] = prepared.prefix3[at3 + category] - prepared.prefix3[from3 + category];
      }
      for (let category = 0; category < 4; category += 1) {
        counts[11 + category] = prepared.prefix4[at4 + category] - prepared.prefix4[from4 + category];
      }
      for (let category = 2; category <= 14; category += 1) {
        totals[category] += counts[category];
        squares[category] += counts[category] * counts[category];
      }
      const sum1 = counts[2] + counts[7] + counts[8];
      const sum2 = counts[3] + counts[6] + counts[9];
      const sum3 = counts[4] + counts[5] + counts[10];
      const sum4 = counts[2] + counts[8] + counts[11] + counts[12];
      const sum5 = counts[3] + counts[9] + counts[11] + counts[13];
      const sum7 = counts[5] + counts[10] + counts[11] + counts[14];
      sumTotals[1] += sum1; sumSquares[1] += sum1 * sum1;
      sumTotals[2] += sum2; sumSquares[2] += sum2 * sum2;
      sumTotals[3] += sum3; sumSquares[3] += sum3 * sum3;
      sumTotals[4] += sum4; sumSquares[4] += sum4 * sum4;
      sumTotals[5] += sum5; sumSquares[5] += sum5 * sum5;
      sumTotals[7] += sum7; sumSquares[7] += sum7 * sum7;
    }
    const patternMeans = new Float64Array(16);
    const patternVariances = new Float64Array(16);
    const sumMeans = new Float64Array(8);
    const sumVariances = new Float64Array(8);
    for (let category = 2; category <= 14; category += 1) {
      const mean = totals[category] / prepared.replicates;
      patternMeans[category] = mean;
      patternVariances[category] = Math.max(0, squares[category] / prepared.replicates - mean * mean);
    }
    for (const group of [1, 2, 3, 4, 5, 7]) {
      const mean = sumTotals[group] / prepared.replicates;
      sumMeans[group] = mean;
      sumVariances[group] = Math.max(0, sumSquares[group] / prepared.replicates - mean * mean);
    }
    moments = { patternMeans, patternVariances, sumMeans, sumVariances };
    momentCache.set(momentKey, moments);
    if (momentCache.size > 8_192) momentCache.delete(momentCache.keys().next().value);
  }
  const patterns = new Float64Array(16);
  const sums = new Float64Array(13);
  for (let category = 1; category <= 15; category += 1) {
    if (actual[category] <= 0) continue;
    const variance = moments.patternVariances[category];
    if (variance > 0) patterns[category] = (actual[category] - moments.patternMeans[category]) / Math.sqrt(variance);
  }
  const observedSums = [
    0,
    actual[2] + actual[7] + actual[8],
    actual[3] + actual[6] + actual[9],
    actual[4] + actual[5] + actual[10],
    actual[2] + actual[8] + actual[11] + actual[12],
    actual[3] + actual[9] + actual[11] + actual[13],
    0,
    actual[5] + actual[10] + actual[11] + actual[14],
  ];
  for (const group of [1, 2, 3, 4, 5, 7]) {
    if (observedSums[group] <= 0) continue;
    const variance = moments.sumVariances[group];
    if (variance > 0) sums[group] = (observedSums[group] - moments.sumMeans[group]) / Math.sqrt(variance);
  }
  return { patterns, sums, replicates: prepared.replicates };
}

function categoryVectorKey(actual) {
  let key = "";
  for (let category = 2; category <= 14; category += 1) {
    if (category > 2) key += ",";
    key += actual[category];
  }
  return key;
}

function retainBounded(cache, key, value, maximum = 16_384) {
  cache.set(key, value);
  if (cache.size > maximum) cache.delete(cache.keys().next().value);
  return value;
}

function pDistance(encoded, length, first, second) {
  let valid = 0;
  let differences = 0;
  const firstOffset = first * length;
  const secondOffset = second * length;
  for (let site = 0; site < length; site += 1) {
    const left = encoded[firstOffset + site];
    const right = encoded[secondOffset + site];
    if (left >= 4 || right >= 4) continue;
    valid += 1;
    if (left !== right) differences += 1;
  }
  return valid ? differences / valid : 1;
}

function matrixDistance(matrix, count, first, second) {
  if (!matrix || matrix.length < count * count) return null;
  const value = matrix[first * count + second];
  return Number.isFinite(value) ? value : null;
}

function directDistance(encoded, length, count, matrix, first, second) {
  return matrixDistance(matrix, count, first, second) ?? pDistance(encoded, length, first, second);
}

export function selectSourceSiScanOutgroup(encoded, length, sequenceCount, triplet, options = {}) {
  if (Number.isInteger(options.outgroupIndex)
      && options.outgroupIndex >= 0
      && options.outgroupIndex < sequenceCount
      && !triplet.includes(options.outgroupIndex)) {
    return { index: options.outgroupIndex, mode: "manual", sampled: false, sourcePath: "analyst-selected fourth sequence" };
  }
  const mode = options.outgroupMode ?? "nearest";
  if (mode === "randomized" || sequenceCount <= 3) {
    return { index: null, mode: "randomized", sampled: false, sourcePath: "horizontal randomization" };
  }
  const allowed = [...new Set((Array.isArray(options.candidatePool)
    ? options.candidatePool
    : Array.from({ length: sequenceCount }, (_, index) => index))
    .filter((index) => index >= 0 && index < sequenceCount))];
  const outside = allowed.filter((index) => !triplet.includes(index));
  if (!outside.length) return { index: null, mode: "randomized", sampled: true, sourcePath: "no fourth sequence available" };
  const distanceMatrix = options.distanceMatrix;
  if (mode === "most-divergent") {
    let best = outside[0];
    let bestDistance = -Infinity;
    for (const candidate of outside) {
      let total = 0;
      let tested = 0;
      const comparisonPool = distanceMatrix?.length >= sequenceCount * sequenceCount
        ? Array.from({ length: sequenceCount }, (_, index) => index)
        : triplet;
      for (const other of comparisonPool) {
        if (other === candidate) continue;
        total += directDistance(encoded, length, sequenceCount, distanceMatrix, candidate, other);
        tested += 1;
      }
      const mean = total / Math.max(1, tested);
      if (mean > bestDistance + 1e-12 || (Math.abs(mean - bestDistance) <= 1e-12 && candidate < best)) {
        best = candidate;
        bestDistance = mean;
      }
    }
    return {
      index: best,
      mode,
      sampled: allowed.length < sequenceCount,
      sourcePath: distanceMatrix?.length >= sequenceCount * sequenceCount
        ? "RDP5 most-divergent whole-matrix rule"
        : "most-divergent bounded candidate rule",
    };
  }

  const [first, second, third] = triplet;
  const tree = options.treeDistanceMatrix;
  if (tree?.length >= sequenceCount * sequenceCount) {
    let inside;
    let outsideLineage;
    const firstSecond = tree[first * sequenceCount + second];
    const firstThird = tree[first * sequenceCount + third];
    if (firstSecond > firstThird) {
      inside = first;
      outsideLineage = third;
    } else if (firstThird > firstSecond) {
      inside = first;
      outsideLineage = second;
    } else {
      inside = second;
      outsideLineage = first;
    }
    const threshold = tree[outsideLineage * sequenceCount + inside];
    let chosen = -1;
    let chosenTreeDistance = 0;
    let chosenDirectDistance = -Infinity;
    for (const candidate of outside) {
      const candidateTreeDistance = tree[outsideLineage * sequenceCount + candidate];
      if (candidate === outsideLineage || !(candidateTreeDistance < threshold) || candidateTreeDistance < chosenTreeDistance) continue;
      const candidateDirectDistance = directDistance(encoded, length, sequenceCount, distanceMatrix, outsideLineage, candidate);
      if (candidateTreeDistance > chosenTreeDistance
          || (candidateTreeDistance === chosenTreeDistance && candidateDirectDistance > chosenDirectDistance)) {
        chosen = candidate;
        chosenTreeDistance = candidateTreeDistance;
        chosenDirectDistance = candidateDirectDistance;
      }
    }
    if (chosen >= 0) {
      return { index: chosen, mode, sampled: allowed.length < sequenceCount, sourcePath: "GetSSOL tree-position rule" };
    }
  }

  // GetSSOL's source fallback is the non-triplet sequence with the smallest
  // mean direct distance to the triplet.  This is also the manual's stated
  // nearest-outlier behavior when a usable tree-position candidate is absent.
  let best = outside[0];
  let bestDistance = Infinity;
  for (const candidate of outside) {
    const mean = triplet.reduce((total, member) => (
      total + directDistance(encoded, length, sequenceCount, distanceMatrix, candidate, member)
    ), 0) / 3;
    if (mean < bestDistance - 1e-12 || (Math.abs(mean - bestDistance) <= 1e-12 && candidate < best)) {
      best = candidate;
      bestDistance = mean;
    }
  }
  return {
    index: best,
    mode,
    sampled: allowed.length < sequenceCount,
    sourcePath: tree ? "GetSSOL direct-distance fallback" : "nearest direct-distance fallback (tree cohort unavailable)",
  };
}

function vb6Random(seed) {
  const buffer = new ArrayBuffer(4);
  const view = new DataView(buffer);
  view.setFloat32(0, -Math.abs(Number(seed) || 1), true);
  const bits = view.getInt32(0, true);
  let state = ((bits & 0xffff) ^ (bits >> 16)) & 0xffffff;
  const next = () => {
    state = (Math.imul(state, 0x43fd43fd) + 0xc39ec3) & 0xffffff;
    return state / 0x1000000;
  };
  next(); // Rnd(-seed) both seeds and returns the first value.
  return next;
}

function horizontalPermutation(size, seed) {
  const random = vb6Random(seed);
  const permutation = new Int32Array(size);
  const taken = new Uint8Array(size);
  for (let index = 0; index < size; index += 1) {
    const initial = Math.min(size - 1, Math.floor(random() * size));
    let selected = initial;
    while (taken[selected]) selected = (selected + 1) % size;
    taken[selected] = 1;
    permutation[index] = selected;
  }
  return permutation;
}

function categoryCounts(encoded, length, triplet, outgroup, start, end, options = {}) {
  const counts = new Int32Array(16);
  const [first, second, third] = triplet;
  const positionMode = options.positionMode ?? "triplet-variable";
  const gapsAsState = options.gapsAsState === true;
  const firstOffset = first * length;
  const secondOffset = second * length;
  const thirdOffset = third * length;
  const outgroupOffset = outgroup === null ? -1 : outgroup * length;
  const regionLength = Math.max(1, end - start);
  const horizontal = outgroup === null
    ? horizontalPermutation(regionLength, (options.seed ?? 0x5a17c0de) ^ Math.imul(start + 1, 0x9e3779b1))
    : null;
  for (let site = start; site < end; site += 1) {
    const relative = site - start;
    const outgroupBase = outgroupOffset >= 0
      ? encoded[outgroupOffset + site]
      : encoded[firstOffset + start + horizontal[relative]];
    const category = sourceSiScanPattern(
      encoded[firstOffset + site],
      encoded[secondOffset + site],
      encoded[thirdOffset + site],
      outgroupBase,
      { gapsAsState },
    );
    if (includePattern(category, positionMode)) counts[category] += 1;
  }
  return counts;
}

function tripletDistances(encoded, length, count, triplet, distanceMatrix) {
  return [
    directDistance(encoded, length, count, distanceMatrix, triplet[0], triplet[1]),
    directDistance(encoded, length, count, distanceMatrix, triplet[0], triplet[2]),
    directDistance(encoded, length, count, distanceMatrix, triplet[1], triplet[2]),
  ];
}

function baselineTopology(distances) {
  if (distances[0] < distances[1] && distances[0] < distances[2]) return 0;
  if (distances[1] < distances[0] && distances[1] < distances[2]) return 1;
  return 2;
}

// Sister-Scanning describes a recombinant as the sequence shared by the
// whole-alignment sister pair and the alternative sister pair inside a run.
// The remaining member of the baseline pair is the major parent and the
// remaining member of the run pair is the minor parent.  SetUpSiScan screens
// one unordered triplet; this conversion is what lets the worker retain that
// single source screen instead of rerunning it under three presumed targets.
export function sourceSiScanRoles(triplet, baseline, inferred) {
  if (!Array.isArray(triplet) || triplet.length !== 3) return null;
  if (![baseline, inferred].every((topology) => Number.isInteger(topology) && topology >= 0 && topology <= 2)) return null;
  if (baseline === inferred) return null;
  const pairSlots = [[0, 1], [0, 2], [1, 2]];
  const baselinePair = pairSlots[baseline];
  const inferredPair = pairSlots[inferred];
  const recombinantSlot = baselinePair.find((slot) => inferredPair.includes(slot));
  if (recombinantSlot === undefined) return null;
  const majorSlot = baselinePair.find((slot) => slot !== recombinantSlot);
  const minorSlot = inferredPair.find((slot) => slot !== recombinantSlot);
  if (majorSlot === undefined || minorSlot === undefined) return null;
  return {
    recombinant: triplet[recombinantSlot],
    majorParent: triplet[majorSlot],
    minorParent: triplet[minorSlot],
  };
}

function winningWindowClass(scores, baseline, extendedGroups) {
  let bestZ = 0;
  let bestIndex = 0;
  let bestFamily = "sum";
  for (let category = 0; category <= 15; category += 1) {
    if (scores.patterns[category] > bestZ) {
      bestZ = scores.patterns[category];
      bestIndex = category;
      bestFamily = "pattern";
    }
  }
  for (let group = 0; group <= 12; group += 1) {
    if (scores.sums[group] > bestZ) {
      bestZ = scores.sums[group];
      bestIndex = group;
      bestFamily = "sum";
    }
  }
  const rows = extendedGroups ? 2 : 1;
  const groups = bestFamily === "pattern" ? PATTERN_TOPOLOGY_GROUPS : SUM_TOPOLOGY_GROUPS;
  for (let row = 0; row < rows; row += 1) {
    for (let topology = 0; topology < 3; topology += 1) {
      if (groups[row][topology] === bestIndex) return { topology, z: bestZ, index: bestIndex, family: bestFamily };
    }
  }
  return { topology: baseline, z: bestZ, index: bestIndex, family: bestFamily };
}

function regionAlternativeScore(scores, baseline, extendedGroups) {
  const alternatives = [0, 1, 2].filter((topology) => topology !== baseline);
  const rows = extendedGroups ? 2 : 1;
  let best = { z: 0, index: 0, family: "pattern", topology: baseline };
  const evaluate = (values, groups, family) => {
    for (let row = 0; row < rows; row += 1) {
      const baselineZ = Math.abs(values[groups[row][baseline]]);
      if (!alternatives.some((topology) => baselineZ < Math.abs(values[groups[row][topology]]))) continue;
      for (const topology of alternatives) {
        const index = groups[row][topology];
        const z = Math.abs(values[index]);
        if (z > best.z) best = { z, index, family, topology };
      }
    }
  };
  evaluate(scores.patterns, PATTERN_TOPOLOGY_GROUPS, "pattern");
  evaluate(scores.sums, SUM_TOPOLOGY_GROUPS, "sum");
  return best;
}

function topologySupport(encoded, length, triplet, topology, site, gapsAsState) {
  const bases = triplet.map((sequence) => encoded[sequence * length + site]);
  if (!bases.every((base) => canonical(base, gapsAsState))) return false;
  if (topology === 0) return bases[0] === bases[1] && bases[0] !== bases[2];
  if (topology === 1) return bases[0] === bases[2] && bases[0] !== bases[1];
  return bases[1] === bases[2] && bases[1] !== bases[0];
}

export function sourceNormalZ(z) {
  const absolute = Math.abs(z);
  if (absolute >= 5.9999999) {
    const exponent = Math.min(170, (absolute - 5.999999) * 10);
    return 1e-9 / Math.pow(1.6, exponent);
  }
  if (absolute === 0) return 1;
  let y = 0.5 * absolute;
  let x;
  if (y >= 3) x = 1;
  else if (y < 1) {
    const w = y * y;
    x = ((((((((0.000124818987 * w - 0.001075204047) * w + 0.005198775019) * w - 0.019198292004) * w + 0.059054035642) * w - 0.151968751364) * w + 0.319152932694) * w - 0.5319230073) * w + 0.797884560593) * y * 2;
  } else {
    y -= 2;
    x = (((((((((((((-0.000045255659 * y + 0.00015252929) * y - 0.000019538132) * y - 0.000676904986) * y + 0.001390604284) * y - 0.00079462082) * y - 0.002034254874) * y + 0.006549791214) * y - 0.010557625006) * y + 0.011630447319) * y - 0.009279453341) * y + 0.005353579108) * y - 0.002141268741) * y + 0.000535310849) * y + 0.999936657524;
  }
  return Math.min(x + 1, 1 - x);
}

export function runSourceSiScan(encoded, length, sequenceCount, triplet, options = {}) {
  const window = Math.max(12, Math.min(length, Math.trunc(options.window ?? 200)));
  const step = Math.max(1, Math.min(window, Math.trunc(options.step ?? 20)));
  const scanPermutations = Math.max(2, Math.trunc(options.scanPermutations ?? 100));
  const pValuePermutations = Math.max(scanPermutations, Math.trunc(options.pValuePermutations ?? 1000));
  const randomization = options.randomization
    ?? buildSourceSiScanRandomization(length, pValuePermutations, options.seed);
  const positionMode = options.positionMode ?? "triplet-variable";
  const gapsAsState = options.gapsAsState === true;
  const extendedGroups = positionMode !== "triplet-variable";
  const outgroup = selectSourceSiScanOutgroup(encoded, length, sequenceCount, triplet, options);
  const distances = tripletDistances(encoded, length, sequenceCount, triplet, options.distanceMatrix);
  const baseline = baselineTopology(distances);
  const windows = [];
  let preparedScanPermutations = null;
  let preparedKey = null;
  if (options.referencePermutationPath !== true) {
    const preparedCache = randomization.preparedScanCache
      ?? (randomization.preparedScanCache = new Map());
    preparedKey = `${window}:${Math.min(scanPermutations, randomization.permutations)}`;
    preparedScanPermutations = preparedCache.get(preparedKey);
    if (!preparedScanPermutations) {
      preparedScanPermutations = buildPermutationPrefix(window, scanPermutations, randomization);
      preparedCache.set(preparedKey, preparedScanPermutations);
    }
  }
  let fixedCategories = null;
  let rollingCounts = null;
  let previousWindowStart = 0;
  if (outgroup.index !== null) {
    fixedCategories = new Uint8Array(length);
    const [first, second, third] = triplet;
    for (let site = 0; site < length; site += 1) {
      const category = sourceSiScanPattern(
        encoded[first * length + site],
        encoded[second * length + site],
        encoded[third * length + site],
        encoded[outgroup.index * length + site],
        { gapsAsState },
      );
      fixedCategories[site] = includePattern(category, positionMode) ? category : 0;
    }
    rollingCounts = new Int32Array(16);
    for (let site = 0; site < window; site += 1) rollingCounts[fixedCategories[site]] += Number(fixedCategories[site] > 0);
  }
  for (let start = 0; start + window <= length; start += step) {
    let counts;
    if (rollingCounts) {
      if (start > 0) {
        const previousEnd = previousWindowStart + window;
        for (let site = previousWindowStart; site < start; site += 1) {
          const category = fixedCategories[site];
          if (category) rollingCounts[category] -= 1;
        }
        for (let site = previousEnd; site < start + window; site += 1) {
          const category = fixedCategories[site];
          if (category) rollingCounts[category] += 1;
        }
        previousWindowStart = start;
      }
      counts = rollingCounts;
    } else {
      counts = categoryCounts(encoded, length, triplet, outgroup.index, start, start + window, {
        positionMode,
        gapsAsState,
        seed: options.seed,
      });
    }
    let scores;
    if (preparedScanPermutations) {
      const scoreCache = randomization.preparedScoreCache
        ?? (randomization.preparedScoreCache = new Map());
      const scoreKey = `${preparedKey}|${categoryVectorKey(counts)}`;
      scores = scoreCache.get(scoreKey)
        ?? retainBounded(scoreCache, scoreKey, prefixedZScores(counts, preparedScanPermutations));
    } else {
      scores = zScores(counts, scanPermutations, randomization);
    }
    const winner = winningWindowClass(scores, baseline, extendedGroups);
    windows.push({ start, end: start + window, center: start + Math.floor(window / 2), ...winner });
  }
  if (!windows.length) return null;

  const runs = [];
  for (let index = 0; index < windows.length;) {
    if (windows[index].topology === baseline) {
      index += 1;
      continue;
    }
    const topology = windows[index].topology;
    const first = index;
    let last = index;
    while (last + 1 < windows.length && windows[last + 1].topology === topology) last += 1;
    let start = windows[first].start;
    let end = windows[last].end;
    while (start < end && !topologySupport(encoded, length, triplet, topology, start, gapsAsState)) start += 1;
    while (end > start && !topologySupport(encoded, length, triplet, topology, end - 1, gapsAsState)) end -= 1;
    if (end > start) runs.push({ first, last, topology, start, end, peakWindowZ: Math.max(...windows.slice(first, last + 1).map((entry) => entry.z)) });
    index = last + 1;
  }

  let best = null;
  const regions = [];
  for (const run of runs) {
    const counts = categoryCounts(encoded, length, triplet, outgroup.index, run.start, run.end, {
      positionMode,
      gapsAsState,
      seed: options.seed,
    });
    const regionScoreCache = randomization.regionScoreCache
      ?? (randomization.regionScoreCache = new Map());
    const regionScoreKey = `${Math.min(pValuePermutations, randomization.permutations)}|${categoryVectorKey(counts)}`;
    let scores = options.referencePermutationPath === true
      ? null
      : regionScoreCache.get(regionScoreKey);
    if (!scores) {
      const preparedRegionPermutations = options.referencePermutationPath === true
        ? null
        : reusableRegionPermutationPrefix(counts, pValuePermutations, randomization);
      scores = preparedRegionPermutations
        ? prefixedZScores(counts, preparedRegionPermutations)
        : zScores(counts, pValuePermutations, randomization);
      if (options.referencePermutationPath !== true) {
        retainBounded(regionScoreCache, regionScoreKey, scores, 4_096);
      }
    }
    const selected = regionAlternativeScore(scores, baseline, extendedGroups);
    if (!(selected.z > 0)) continue;
    const regionLength = run.end - run.start;
    const rawP = Math.min(1, sourceNormalZ(selected.z) * (length / Math.max(1, regionLength)));
    const result = {
      ...run,
      z: selected.z,
      rawP,
      pattern: selected.index,
      scoreFamily: selected.family,
      inferredTopology: selected.topology,
      scanPermutations: Math.min(scanPermutations, randomization.permutations),
      pValuePermutations: Math.min(pValuePermutations, randomization.permutations),
    };
    regions.push(result);
    if (!best || result.rawP < best.rawP - 1e-15 || (Math.abs(result.rawP - best.rawP) <= 1e-15 && result.z > best.z)) best = result;
  }
  if (!best) return null;
  return {
    ...best,
    baselineTopology: baseline,
    distances,
    outgroupIndex: outgroup.index,
    outgroupMode: outgroup.mode,
    outgroupSampled: outgroup.sampled,
    outgroupSourcePath: outgroup.sourcePath,
    positionMode,
    gapsAsState,
    window,
    step,
    windows,
    // A triplet can contain several separate topology runs.  RDP5 queues
    // those as separate event hypotheses; retain every calibrated run so the
    // worker can associate the run co-located with each independently found
    // event rather than incorrectly reusing only the global optimum.
    regions: regions.sort((left, right) => left.rawP - right.rawP || right.z - left.z),
    sourceRoutine: "GetSSOL + Get3Score/GetPScores2 + DoPerms3P + MakeZValue2 + DoSums + FindMaxZ + ShrinkRegionC",
  };
}
