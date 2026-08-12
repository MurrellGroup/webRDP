// Seeded bootstrapped neighbor-joining evidence used by the RDP5
// co-recombinant screen. The desktop source constructs six JC/NJ trees per
// event, collapses branches below 50% bootstrap support, and compares tree
// positions across three region pairs. Long regions are represented by
// balanced, interleaved site blocks so every site contributes while keeping
// browser work and memory bounded; setting bootstrapBlocks above the number
// of sites gives ordinary site resampling exactly.

function xorshift32(seed) {
  let state = (seed >>> 0) || 0x9e3779b9;
  return () => {
    state ^= state << 13;
    state ^= state >>> 17;
    state ^= state << 5;
    return state >>> 0;
  };
}

function mixSeed(seed, value) {
  let mixed = (seed ^ Math.imul(value + 1, 0x9e3779b1)) >>> 0;
  mixed ^= mixed >>> 16;
  mixed = Math.imul(mixed, 0x85ebca6b) >>> 0;
  mixed ^= mixed >>> 13;
  return mixed >>> 0;
}

function regionSites(region, length) {
  if (region.kind === "sites") {
    return [...new Set([...region.sites].filter((site) => site >= 0 && site < length))];
  }
  const sites = [];
  for (const [rawStart, rawEnd] of region.segments ?? []) {
    const start = Math.max(0, Math.min(length, Math.trunc(rawStart)));
    const end = Math.max(start, Math.min(length, Math.trunc(rawEnd)));
    for (let site = start; site < end; site += 1) sites.push(site);
  }
  return sites;
}

function canonicalSplit(leaves, taxa) {
  if (leaves.size < 2 || leaves.size > taxa.length - 2) return null;
  const complement = taxa.filter((taxon) => !leaves.has(taxon));
  const selected = [...leaves].sort((left, right) => left - right);
  complement.sort((left, right) => left - right);
  const left = selected.join(",");
  const right = complement.join(",");
  if (selected.length < complement.length) return { signature: left, members: new Set(selected) };
  if (complement.length < selected.length) return { signature: right, members: new Set(complement) };
  return left <= right
    ? { signature: left, members: new Set(selected) }
    : { signature: right, members: new Set(complement) };
}

function neighborJoiningSplits(distances, taxa) {
  const count = taxa.length;
  if (count < 4) return [];
  const capacity = count * 2;
  const matrix = new Float64Array(capacity * capacity);
  for (let left = 0; left < count; left += 1) {
    for (let right = left + 1; right < count; right += 1) {
      const value = Math.max(0, Number.isFinite(distances[left * count + right]) ? distances[left * count + right] : 3);
      matrix[left * capacity + right] = value;
      matrix[right * capacity + left] = value;
    }
  }
  const clusters = new Map(taxa.map((taxon, index) => [index, new Set([taxon])]));
  const splits = new Map();
  let active = Array.from({ length: count }, (_, index) => index);
  let nextId = count;
  while (active.length > 2) {
    const sums = new Float64Array(capacity);
    for (let leftIndex = 0; leftIndex < active.length; leftIndex += 1) {
      const left = active[leftIndex];
      let total = 0;
      for (let rightIndex = 0; rightIndex < active.length; rightIndex += 1) {
        if (rightIndex !== leftIndex) total += matrix[left * capacity + active[rightIndex]];
      }
      sums[left] = total;
    }
    let bestLeft = active[0];
    let bestRight = active[1];
    let bestQ = Infinity;
    for (let leftIndex = 0; leftIndex < active.length - 1; leftIndex += 1) {
      const left = active[leftIndex];
      for (let rightIndex = leftIndex + 1; rightIndex < active.length; rightIndex += 1) {
        const right = active[rightIndex];
        const q = (active.length - 2) * matrix[left * capacity + right] - sums[left] - sums[right];
        if (q < bestQ - 1e-12 || (Math.abs(q - bestQ) <= 1e-12 && (left < bestLeft || (left === bestLeft && right < bestRight)))) {
          bestQ = q;
          bestLeft = left;
          bestRight = right;
        }
      }
    }
    const pairDistance = matrix[bestLeft * capacity + bestRight];
    const merged = nextId++;
    const leaves = new Set([...(clusters.get(bestLeft) ?? []), ...(clusters.get(bestRight) ?? [])]);
    clusters.set(merged, leaves);
    const split = canonicalSplit(leaves, taxa);
    if (split) splits.set(split.signature, split);
    for (const other of active) {
      if (other === bestLeft || other === bestRight) continue;
      const distance = Math.max(0, 0.5 * (
        matrix[bestLeft * capacity + other]
        + matrix[bestRight * capacity + other]
        - pairDistance
      ));
      matrix[merged * capacity + other] = distance;
      matrix[other * capacity + merged] = distance;
    }
    active = active.filter((id) => id !== bestLeft && id !== bestRight);
    active.push(merged);
  }
  return [...splits.values()];
}

function jcDistance(differences, valid) {
  if (valid < 1) return 3;
  const p = differences / valid;
  return p >= 0.749999 ? 3 : Math.max(0, -0.75 * Math.log(1 - 4 * p / 3));
}

function buildBlockCounts(encoded, length, taxa, sites, blockCount) {
  const count = taxa.length;
  const cells = blockCount * count * count;
  const valid = new Uint32Array(cells);
  const differences = new Uint32Array(cells);
  for (let siteIndex = 0; siteIndex < sites.length; siteIndex += 1) {
    const site = sites[siteIndex];
    const blockOffset = (siteIndex % blockCount) * count * count;
    for (let left = 0; left < count - 1; left += 1) {
      const leftBase = encoded[taxa[left] * length + site];
      if (leftBase >= 4) continue;
      for (let right = left + 1; right < count; right += 1) {
        const rightBase = encoded[taxa[right] * length + site];
        if (rightBase >= 4) continue;
        const cell = blockOffset + left * count + right;
        valid[cell] += 1;
        if (leftBase !== rightBase) differences[cell] += 1;
      }
    }
  }
  return { valid, differences };
}

function distanceMatrixFromBlocks(blocks, count, blockCount, weights = null) {
  const matrix = new Float64Array(count * count);
  const blockStride = count * count;
  for (let left = 0; left < count - 1; left += 1) {
    for (let right = left + 1; right < count; right += 1) {
      let valid = 0;
      let differences = 0;
      const pair = left * count + right;
      for (let block = 0; block < blockCount; block += 1) {
        const weight = weights ? weights[block] : 1;
        if (!weight) continue;
        const cell = block * blockStride + pair;
        valid += blocks.valid[cell] * weight;
        differences += blocks.differences[cell] * weight;
      }
      const distance = jcDistance(differences, valid);
      matrix[left * count + right] = distance;
      matrix[right * count + left] = distance;
    }
  }
  return matrix;
}

function pathMatrix(taxa, splits, support, cutoff) {
  const count = taxa.length;
  const matrix = new Float64Array(count * count);
  for (let left = 0; left < count - 1; left += 1) {
    for (let right = left + 1; right < count; right += 1) {
      let distance = 2;
      for (const split of splits) {
        if ((support.get(split.signature) ?? 0) < cutoff) continue;
        if (split.members.has(taxa[left]) !== split.members.has(taxa[right])) distance += 1;
      }
      matrix[left * count + right] = distance;
      matrix[right * count + left] = distance;
    }
  }
  return matrix;
}

// RDP5's nearest-SiScan-outlier rule consumes a whole-alignment tree-distance
// matrix.  Reuse the same deterministic NJ topology representation as the
// event bootstrap code when the worker has an exact all-taxon distance matrix.
export function buildNeighborJoiningPathMatrix(distances, count) {
  const size = Math.max(0, Math.trunc(count));
  const taxa = Array.from({ length: size }, (_, index) => index);
  if (size < 2 || distances.length < size * size) return new Float64Array(size * size);
  const splits = neighborJoiningSplits(distances, taxa);
  return pathMatrix(taxa, splits, new Map(splits.map((split) => [split.signature, 1])), 0);
}

export function buildBootstrapRegionTree(encoded, length, taxa, region, options = {}, seed = 0x5a17c0de) {
  const sites = regionSites(region, length);
  const replicates = Math.max(0, Math.min(1000, Math.trunc(options.replicates ?? 100)));
  const cutoff = Math.max(0, Math.min(1, Number(options.cutoff ?? 0.5)));
  const maximumBlocks = Math.max(8, Math.min(4096, Math.trunc(options.bootstrapBlocks ?? 128)));
  if (taxa.length < 4 || sites.length < 2) {
    return { taxa, index: new Map(taxa.map((taxon, index) => [taxon, index])), sites: sites.length, blocks: 0, replicates, cutoff, exactSiteBootstrap: true, baseDistances: new Float64Array(taxa.length ** 2), collapsed: new Float64Array(taxa.length ** 2), uncollapsed: new Float64Array(taxa.length ** 2), splits: [] };
  }
  const blockCount = Math.min(sites.length, maximumBlocks);
  const blocks = buildBlockCounts(encoded, length, taxa, sites, blockCount);
  const baseDistances = distanceMatrixFromBlocks(blocks, taxa.length, blockCount);
  const baseSplits = neighborJoiningSplits(baseDistances, taxa);
  const counts = new Map(baseSplits.map((split) => [split.signature, 0]));
  const random = xorshift32(seed);
  const weights = new Uint16Array(blockCount);
  for (let replicate = 0; replicate < replicates; replicate += 1) {
    weights.fill(0);
    for (let draw = 0; draw < blockCount; draw += 1) weights[random() % blockCount] += 1;
    const sampled = distanceMatrixFromBlocks(blocks, taxa.length, blockCount, weights);
    const replicateSplits = new Set(neighborJoiningSplits(sampled, taxa).map((split) => split.signature));
    for (const signature of counts.keys()) if (replicateSplits.has(signature)) counts.set(signature, counts.get(signature) + 1);
  }
  const support = new Map([...counts].map(([signature, count]) => [signature, replicates > 0 ? count / replicates : 1]));
  const collapsed = pathMatrix(taxa, baseSplits, support, cutoff);
  const uncollapsed = pathMatrix(taxa, baseSplits, new Map(baseSplits.map((split) => [split.signature, 1])), 0);
  return {
    taxa,
    index: new Map(taxa.map((taxon, index) => [taxon, index])),
    sites: sites.length,
    blocks: blockCount,
    replicates,
    cutoff,
    exactSiteBootstrap: blockCount === sites.length,
    baseDistances,
    collapsed,
    uncollapsed,
    splits: baseSplits.map((split) => ({ signature: split.signature, members: [...split.members], support: support.get(split.signature) ?? 0 })),
  };
}

function selectTreeTaxa(sequenceCount, required, limit, seed) {
  if (sequenceCount <= limit) return Array.from({ length: sequenceCount }, (_, index) => index);
  const selected = new Set(required);
  const candidates = Array.from({ length: sequenceCount }, (_, index) => index).filter((index) => !selected.has(index));
  const random = xorshift32(seed);
  for (let index = candidates.length - 1; index > 0; index -= 1) {
    const swap = random() % (index + 1);
    [candidates[index], candidates[swap]] = [candidates[swap], candidates[index]];
  }
  for (const candidate of candidates) {
    if (selected.size >= limit) break;
    selected.add(candidate);
  }
  return [...selected].sort((left, right) => left - right);
}

export function buildEventBootstrapTrees(encoded, length, sequenceCount, triplet, regionPairs, options = {}, eventIndex = 0) {
  const limit = Math.max(4, Math.min(300, Math.trunc(options.clusterTreeTaxaLimit ?? 32)));
  const seed = mixSeed(options.randomSeed ?? 0x5a17c0de, eventIndex);
  const taxa = selectTreeTaxa(sequenceCount, triplet, Math.max(limit, triplet.length), seed);
  const uniqueRegions = new Map();
  for (const pair of regionPairs) for (const region of pair) uniqueRegions.set(region.id, region);
  const trees = new Map();
  let regionIndex = 0;
  for (const [id, region] of uniqueRegions) {
    trees.set(id, buildBootstrapRegionTree(encoded, length, taxa, region, {
      replicates: options.clusterBootstrapReplicates ?? 100,
      cutoff: options.clusterBootstrapCutoff ?? 0.5,
      bootstrapBlocks: options.clusterBootstrapBlocks ?? 128,
    }, mixSeed(seed, regionIndex++)));
  }
  return {
    taxa,
    sourceSequenceCount: sequenceCount,
    sampled: taxa.length < sequenceCount,
    pairs: regionPairs.map((pair) => pair.map((region) => trees.get(region.id))),
  };
}

// A fixed cohort cap must never turn into a fixed *candidate* cap.  RDP's
// co-recombinant screen considers every sequence in the active alignment.  For
// browser-sized tree work, partition the non-triplet taxa into deterministic
// cohorts while retaining the detecting triplet in every tree.  Each candidate
// is therefore present in one complete six-tree bundle, and all cohorts reuse
// the same region-specific bootstrap stream.  This keeps work approximately
// linear in the number of cohorts instead of silently dropping taxa outside a
// single sampled tree.
export function buildEventBootstrapTreeCohorts(
  encoded,
  length,
  sequenceCount,
  triplet,
  regionPairs,
  options = {},
  eventIndex = 0,
) {
  const limit = Math.max(4, Math.min(300, Math.trunc(options.clusterTreeTaxaLimit ?? 32)));
  const required = [...new Set(triplet)]
    .filter((taxon) => Number.isInteger(taxon) && taxon >= 0 && taxon < sequenceCount)
    .sort((left, right) => left - right);
  if (sequenceCount <= limit || required.length >= sequenceCount) {
    const primary = buildEventBootstrapTrees(
      encoded,
      length,
      sequenceCount,
      required,
      regionPairs,
      options,
      eventIndex,
    );
    return {
      primary,
      cohorts: [primary],
      byTaxon: new Map(primary.taxa.map((taxon) => [taxon, primary])),
      candidateComplete: primary.taxa.length === sequenceCount,
    };
  }

  const requiredSet = new Set(required);
  const candidates = Array.from({ length: sequenceCount }, (_, index) => index)
    .filter((index) => !requiredSet.has(index));
  const random = xorshift32(mixSeed(options.randomSeed ?? 0x5a17c0de, eventIndex));
  for (let index = candidates.length - 1; index > 0; index -= 1) {
    const swap = random() % (index + 1);
    [candidates[index], candidates[swap]] = [candidates[swap], candidates[index]];
  }
  const cohortCapacity = Math.max(1, limit - required.length);
  const cohorts = [];
  const byTaxon = new Map();
  for (let offset = 0; offset < candidates.length; offset += cohortCapacity) {
    const taxa = [...required, ...candidates.slice(offset, offset + cohortCapacity)]
      .sort((left, right) => left - right);
    const bundle = buildEventBootstrapTrees(
      encoded,
      length,
      sequenceCount,
      taxa,
      regionPairs,
      { ...options, clusterTreeTaxaLimit: taxa.length },
      eventIndex,
    );
    cohorts.push(bundle);
    for (const taxon of taxa) if (!requiredSet.has(taxon)) byTaxon.set(taxon, bundle);
  }
  const primary = cohorts[0];
  for (const taxon of required) byTaxon.set(taxon, primary);
  return { primary, cohorts, byTaxon, candidateComplete: byTaxon.size === sequenceCount };
}

export function bootstrapTreeDistance(tree, first, second, collapsed = true) {
  const left = tree?.index.get(first);
  const right = tree?.index.get(second);
  if (left === undefined || right === undefined) return Infinity;
  return (collapsed ? tree.collapsed : tree.uncollapsed)[left * tree.taxa.length + right];
}
