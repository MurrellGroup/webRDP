// RDP5 PHI-test compatibility path.
//
// This is a direct browser translation of the author-supplied RDP5
// `PHITest2`, `PHI`, `pair_score`, `GetFandG`, and
// `AnalyticMeanVariance` routines, themselves described in the desktop source
// as a VB translation of PHIPACK.  The test uses parsimony-informative DNA
// sites, the minimum-reticulation score of each two-site bipartite state graph,
// and the analytic normal approximation from Bruen et al. (2006).
//
// RDP5 randomly thins alignments above 6,000 informative sites.  Browser runs
// use a deterministic, position-balanced subset when a lower work ceiling is
// requested.  The returned provenance makes that distinction impossible to
// confuse with an all-site result.

import { erfc } from "./rdp-statistics.js";

function clamp(value, lower, upper) {
  return Math.max(lower, Math.min(upper, value));
}

// VB6 CLng uses banker's rounding rather than JavaScript's half-up Math.round.
function roundHalfEven(value) {
  if (!Number.isFinite(value)) return 0;
  const floor = Math.floor(value);
  const fraction = value - floor;
  if (fraction < 0.5) return floor;
  if (fraction > 0.5) return floor + 1;
  return floor % 2 === 0 ? floor : floor + 1;
}

function normalCdf(value) {
  if (value === Infinity) return 1;
  if (value === -Infinity) return 0;
  return 0.5 * erfc(-value / Math.SQRT2);
}

function informativeSites(encoded, sequenceCount, length) {
  const sites = [];
  const stateCounts = [];
  for (let site = 0; site < length; site += 1) {
    const counts = new Uint32Array(4);
    for (let sequence = 0; sequence < sequenceCount; sequence += 1) {
      const base = encoded[sequence * length + site];
      if (base < 4) counts[base] += 1;
    }
    let states = 0;
    let repeatedStates = 0;
    for (let base = 0; base < 4; base += 1) {
      if (counts[base] > 0) states += 1;
      if (counts[base] >= 2) repeatedStates += 1;
    }
    if (states >= 2 && repeatedStates >= 2) {
      sites.push(site);
      stateCounts.push(states);
    }
  }
  return { sites, stateCounts };
}

function balancedSubset(values, parallelValues, limit) {
  if (values.length <= limit) return { values, parallelValues, sampled: false };
  const selected = new Array(limit);
  const selectedParallel = new Array(limit);
  for (let index = 0; index < limit; index += 1) {
    const source = Math.floor(index * (values.length - 1) / Math.max(1, limit - 1));
    selected[index] = values[source];
    selectedParallel[index] = parallelValues[source];
  }
  return { values: selected, parallelValues: selectedParallel, sampled: true };
}

// `pair_score` constructs a bipartite graph whose vertices are nucleotide
// states at the two sites and returns edges - vertices + components.  The
// compact implementation below is algebraically identical but avoids a heap
// allocation and linked-list DFS for every site pair.
export function sourcePairIncompatibility(encoded, sequenceCount, length, leftSite, rightSite, leftStates = null, rightStates = null) {
  let edgeMask = 0;
  let leftMask = 0;
  let rightMask = 0;
  for (let sequence = 0; sequence < sequenceCount; sequence += 1) {
    const left = encoded[sequence * length + leftSite];
    const right = encoded[sequence * length + rightSite];
    if (left < 4) leftMask |= 1 << left;
    if (right < 4) rightMask |= 1 << right;
    if (left < 4 && right < 4) edgeMask |= 1 << (left * 4 + right);
  }
  const leftVertices = [];
  const rightVertices = [];
  for (let state = 0; state < 4; state += 1) {
    if (leftMask & (1 << state)) leftVertices.push(state);
    if (rightMask & (1 << state)) rightVertices.push(state);
  }
  const declaredLeftStates = leftStates ?? leftVertices.length;
  const declaredRightStates = rightStates ?? rightVertices.length;
  const vertexCount = declaredLeftStates + declaredRightStates;
  if (vertexCount < 2) return 0;

  const adjacency = new Uint16Array(8);
  let edges = 0;
  for (const left of leftVertices) {
    for (const right of rightVertices) {
      if (!(edgeMask & (1 << (left * 4 + right)))) continue;
      const rightVertex = 4 + right;
      adjacency[left] |= 1 << rightVertex;
      adjacency[rightVertex] |= 1 << left;
      edges += 1;
    }
  }
  let components = 0;
  let remaining = leftMask | (rightMask << 4);
  while (remaining) {
    components += 1;
    const first = remaining & -remaining;
    let frontier = first;
    remaining &= ~first;
    while (frontier) {
      const vertexBit = frontier & -frontier;
      frontier &= ~vertexBit;
      const vertex = 31 - Math.clz32(vertexBit);
      const neighbours = adjacency[vertex] & remaining;
      frontier |= neighbours;
      remaining &= ~neighbours;
    }
  }
  return Math.max(0, edges - vertexCount + components);
}

function incompatibilityMatrix(encoded, sequenceCount, length, sites, stateCounts) {
  const count = sites.length;
  const matrix = new Uint8Array(count * count);
  for (let left = 0; left < count - 1; left += 1) {
    for (let right = left + 1; right < count; right += 1) {
      const score = sourcePairIncompatibility(
        encoded,
        sequenceCount,
        length,
        sites[left],
        sites[right],
        stateCounts[left],
        stateCounts[right],
      );
      matrix[left * count + right] = score;
      matrix[right * count + left] = score;
    }
  }
  return matrix;
}

export function sourcePhiStatistic(matrix, count, k) {
  let score = 0;
  for (let left = 0; left < count - 1; left += 1) {
    for (let distance = 1; distance <= k && left + distance < count; distance += 1) {
      score += matrix[left * count + left + distance];
    }
  }
  const terms = k * (2 * count - k - 1) / 2;
  return terms > 0 ? score / terms : 0;
}

export function sourcePhiAnalyticMeanVariance(matrix, count, k) {
  const fValues = new Float64Array(count);
  const gValues = new Float64Array(count);
  let sumF = 0;
  for (let row = 0; row < count; row += 1) {
    let f = 0;
    let g = 0;
    const offset = row * count;
    for (let column = 0; column < count; column += 1) {
      const value = matrix[offset + column];
      f += value;
      g += value * value;
    }
    fValues[row] = f;
    gValues[row] = g;
    sumF += f;
  }
  const mean = sumF / Math.max(1, count * (count - 1));
  const n = count;
  const top1 = 27 * k * n - 18 * k ** 2 + 28 * k ** 2 * n - 21 * k * n ** 2 - 9 * k + 5 * n
    - 9 * k ** 3 - 11 * n ** 2 + 6 * n ** 3 + 6 * k ** 3 * n - 4 * k ** 2 * n ** 2;
  const bottom1 = k * (k + 1 - 2 * n) ** 2 * (n - 1) ** 2 * (n - 2) * (n - 3) * n ** 2;
  const top2 = 8 * k ** 2 * n - 14 * k ** 2 + 39 * k * n + 19 * n - 21 * k + 3 * k ** 3
    - 15 * k * n ** 2 + 6 * n ** 3 - 21 * n ** 2 - 4;
  const bottom2 = k * (2 * n - k - 1) ** 2 * n * (n - 1) * (n - 2) * (n - 3);
  const top3 = -18 * k * n - 2 * k ** 2 * n + 16 * k ** 2 + 6 * n ** 2 - 10 * n + 2 + 15 * k + 3 * k ** 3;
  const bottom3 = k * (k + 1 - 2 * n) ** 2 * n * (n - 1) * (n - 2) * (n - 3);
  if (!(bottom1 > 0 && bottom2 > 0 && bottom3 > 0)) return { mean, variance: 0 };
  const coefficient1 = (2 / 3) * top1 / bottom1;
  const coefficient2 = (2 / 3) * top2 / bottom2;
  const coefficient3 = (-4 / 3) * top3 / bottom3;
  let sumG = 0;
  let sumFSquared = 0;
  for (let index = 0; index < count; index += 1) {
    sumG += gValues[index];
    sumFSquared += fValues[index] ** 2;
  }
  const variance = coefficient1 * sumF ** 2 + coefficient2 * sumG + coefficient3 * sumFSquared;
  return { mean, variance: Number.isFinite(variance) ? Math.max(0, variance) : 0 };
}

export function sourcePhiTest(encoded, sequenceCount, length, options = {}) {
  if (!(encoded instanceof Uint8Array) || sequenceCount < 4 || length < 1) {
    return {
      pValue: 1,
      statistic: 0,
      mean: 0,
      variance: 0,
      informativeSites: 0,
      totalInformativeSites: 0,
      k: 0,
      window: Math.max(1, Math.trunc(options.window ?? 100)),
      subsampled: false,
      validNormalApproximation: false,
      compatibility: "RDP5 PHITest2/PHI analytic normal approximation",
    };
  }
  const maximumSites = clamp(Math.trunc(options.maxInformativeSites ?? 384), 16, 6000);
  const requestedWindow = clamp(Math.trunc(options.window ?? 100), 1, Math.max(1, length));
  const collected = informativeSites(encoded, sequenceCount, length);
  const subset = balancedSubset(collected.sites, collected.stateCounts, maximumSites);
  const sites = subset.values;
  const stateCounts = subset.parallelValues;
  const count = sites.length;
  if (count < 4) {
    return {
      pValue: 1,
      statistic: 0,
      mean: 0,
      variance: 0,
      informativeSites: count,
      totalInformativeSites: collected.sites.length,
      k: 0,
      window: requestedWindow,
      subsampled: subset.sampled,
      validNormalApproximation: false,
      compatibility: "RDP5 PHITest2/PHI analytic normal approximation",
    };
  }
  let k = roundHalfEven(count / length * requestedWindow);
  if (k < 10) k = 10;
  if (k > count) k = Math.max(1, roundHalfEven(count / 2) - 1);
  k = clamp(k, 1, Math.max(1, count - 1));
  const matrix = incompatibilityMatrix(encoded, sequenceCount, length, sites, stateCounts);
  const statistic = sourcePhiStatistic(matrix, count, k);
  const { mean, variance } = sourcePhiAnalyticMeanVariance(matrix, count, k);
  const validNormalApproximation = count > 2 * k && variance > 0;
  const z = validNormalApproximation ? (statistic - mean) / Math.sqrt(variance) : 0;
  return {
    pValue: validNormalApproximation ? clamp(normalCdf(z), 0, 1) : 1,
    statistic,
    mean,
    variance,
    z,
    informativeSites: count,
    totalInformativeSites: collected.sites.length,
    k,
    window: requestedWindow,
    subsampled: subset.sampled,
    validNormalApproximation,
    compatibility: "RDP5 PHITest2/PHI/pair_score/GetFandG/AnalyticMeanVariance",
  };
}

