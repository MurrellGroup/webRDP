// RDP5 ancestral-event clustering for RDP Web.
//
// Implemented from the author-supplied RDP5 source and manual with permission.
// It keeps the source workflow's two distinct operations explicit: GetSupers
// collapses similar detections with a 0.1-threshold weighted-average merge;
// the §4.1.4 three-set procedure then identifies co-recombinant descendants.

import { studentTTwoSided } from "./rdp-statistics.js";
import { bootstrapTreeDistance, buildEventBootstrapTreeCohorts } from "./rdp-bootstrap-tree.js";
import { identifyRecombinantRoles } from "./rdp-recombinant-identification.js";

function clamp(value, lower, upper) {
  return Math.max(lower, Math.min(upper, value));
}

function tractSegments(event, length) {
  if (event.wraps && event.start > event.end) {
    return [[event.start, length], ...(event.end > 0 ? [[0, event.end]] : [])];
  }
  return event.end > event.start ? [[event.start, event.end]] : [];
}

function complementSegments(event, length) {
  if (event.wraps && event.start > event.end) return event.end < event.start ? [[event.end, event.start]] : [];
  return [[0, event.start], ...(event.end < length ? [[event.end, length]] : [])].filter(([start, end]) => end > start);
}

function segmentLength(segments) {
  return segments.reduce((total, [start, end]) => total + Math.max(0, end - start), 0);
}

function intersectionBases(left, right, length) {
  let shared = 0;
  for (const [leftStart, leftEnd] of tractSegments(left, length)) {
    for (const [rightStart, rightEnd] of tractSegments(right, length)) {
      shared += Math.max(0, Math.min(leftEnd, rightEnd) - Math.max(leftStart, rightStart));
    }
  }
  return shared;
}

export function reciprocalTractOverlap(left, right, length) {
  const shared = intersectionBases(left, right, length);
  const leftLength = segmentLength(tractSegments(left, length));
  const rightLength = segmentLength(tractSegments(right, length));
  return {
    bases: shared,
    left: shared / Math.max(1, leftLength),
    right: shared / Math.max(1, rightLength),
    minimum: shared / Math.max(1, Math.min(leftLength, rightLength)),
    jaccard: shared / Math.max(1, leftLength + rightLength - shared),
  };
}

function triad(event) {
  return [event.recombinant, event.majorParent, event.minorParent];
}

function sharedMembers(left, right) {
  const rightSet = new Set(right);
  return [...new Set(left)].filter((value) => rightSet.has(value)).length;
}

function tripletVariable(encoded, length, triplet, site) {
  const first = encoded[triplet[0] * length + site];
  const second = encoded[triplet[1] * length + site];
  const third = encoded[triplet[2] * length + site];
  return first < 4 && second < 4 && third < 4 && (first !== second || first !== third);
}

function flankSites(encoded, length, triplet, breakpoint, direction, wanted, circular) {
  const sites = [];
  const maximumSteps = circular ? length : direction < 0 ? breakpoint : length - breakpoint;
  for (let step = 0; step < maximumSteps && sites.length < wanted; step += 1) {
    let site = direction < 0 ? breakpoint - step - 1 : breakpoint + step;
    if (circular) site = (site % length + length) % length;
    if (site < 0 || site >= length) break;
    if (tripletVariable(encoded, length, triplet, site)) sites.push(site);
  }
  sites.sort((left, right) => left - right);
  return { kind: "sites", sites: Int32Array.from(sites) };
}

function makeEventRegions(event, encoded, length, options) {
  const eventTriad = triad(event);
  const flankVnps = Math.max(4, Math.min(200, Math.trunc(options.clusterFlankVnps ?? 60)));
  const circular = options.circular === true || event.wraps === true;
  const startLeft = flankSites(encoded, length, eventTriad, event.start, -1, flankVnps, circular);
  const startRight = flankSites(encoded, length, eventTriad, event.start, 1, flankVnps, circular);
  const endLeft = flankSites(encoded, length, eventTriad, event.end % length, -1, flankVnps, circular);
  const endRight = flankSites(encoded, length, eventTriad, event.end % length, 1, flankVnps, circular);
  const tract = { kind: "segments", segments: tractSegments(event, length) };
  const background = { kind: "segments", segments: complementSegments(event, length) };
  return {
    // RDP5 MakeBPosLR/MakeSDMP2 region numbering:
    // 0/1 flank the beginning, 3/2 flank the ending, and region 4 is the
    // breakpoint-bounded tract. FillRmat's third comparison averages regions
    // 0 and 3 before comparing that vector with region 4.
    sourcePairs: [
      { first: [startLeft], second: [startRight], sdmRegion: startRight },
      { first: [endRight], second: [endLeft], sdmRegion: endLeft },
      { first: [startLeft, endRight], second: [tract], sdmRegion: tract },
    ],
    // The tree set uses the two breakpoint pairs and the full tract versus its
    // complement, as described in manual section 4.1.4.
    phylogeneticPairs: [
      [startLeft, startRight],
      [endRight, endLeft],
      [background, tract],
    ],
  };
}

function visitRegion(region, callback) {
  if (region.kind === "sites") {
    for (const site of region.sites) callback(site);
    return;
  }
  for (const [start, end] of region.segments) for (let site = start; site < end; site += 1) callback(site);
}

function jukesCantorDistance(encoded, length, firstSequence, secondSequence, region, cache) {
  const low = Math.min(firstSequence, secondSequence);
  const high = Math.max(firstSequence, secondSequence);
  const key = `${region.id}:${low}:${high}`;
  const cached = cache.get(key);
  if (cached !== undefined) return cached;
  let valid = 0;
  let differences = 0;
  const firstOffset = firstSequence * length;
  const secondOffset = secondSequence * length;
  visitRegion(region, (site) => {
    const first = encoded[firstOffset + site];
    const second = encoded[secondOffset + site];
    if (first >= 4 || second >= 4) return;
    valid += 1;
    if (first !== second) differences += 1;
  });
  if (valid === 0) {
    cache.set(key, 3);
    return 3;
  }
  const p = differences / valid;
  const distance = p >= 0.749999 ? 3 : Math.max(0, -0.75 * Math.log(1 - 4 * p / 3));
  cache.set(key, distance);
  return distance;
}

function pearson(left, right) {
  const count = Math.min(left.length, right.length);
  let leftMean = 0;
  let rightMean = 0;
  for (let index = 0; index < count; index += 1) {
    leftMean += left[index];
    rightMean += right[index];
  }
  leftMean /= count;
  rightMean /= count;
  let covariance = 0;
  let leftVariance = 0;
  let rightVariance = 0;
  for (let index = 0; index < count; index += 1) {
    const leftDelta = left[index] - leftMean;
    const rightDelta = right[index] - rightMean;
    covariance += leftDelta * rightDelta;
    leftVariance += leftDelta * leftDelta;
    rightVariance += rightDelta * rightDelta;
  }
  if (!(leftVariance > 0) || !(rightVariance > 0)) return 0;
  return clamp(covariance / Math.sqrt(leftVariance * rightVariance), -0.999999999, 0.999999999);
}

function sourceRegionCategoryVector(encoded, length, presumed, parent1, parent2, candidate, region) {
  const counts = [0, 0, 0];
  let total = 0;
  visitRegion(region, (site) => {
    const presumedBase = encoded[presumed * length + site];
    const first = encoded[parent1 * length + site];
    const second = encoded[parent2 * length + site];
    const value = encoded[candidate * length + site];
    if (presumedBase >= 4 || first >= 4 || second >= 4 || value >= 4) return;
    if (value === first && value !== second) {
      counts[0] += 1;
      total += 1;
    } else if (value === second && value !== first) {
      counts[1] += 1;
      total += 1;
    } else if (first === second && value !== first) {
      counts[2] += 1;
      total += 1;
    }
  });
  return total > 0 ? counts.map((count) => count / total) : [10, 10, 10];
}

function sourceAveragedRegionVector(encoded, length, presumed, parent1, parent2, candidate, regions) {
  const vectors = regions.map((region) => sourceRegionCategoryVector(
    encoded,
    length,
    presumed,
    parent1,
    parent2,
    candidate,
    region,
  ));
  return [0, 1, 2].map((category) => vectors.reduce((total, vector) => total + vector[category], 0) / vectors.length);
}

function sourceCategoryVector(encoded, length, presumed, parent1, parent2, candidate, pair) {
  return [
    ...sourceAveragedRegionVector(encoded, length, presumed, parent1, parent2, candidate, pair.first),
    ...sourceAveragedRegionVector(encoded, length, presumed, parent1, parent2, candidate, pair.second),
  ];
}

function sourcePearsonSix(left, right) {
  if (left.some((value) => value > 4) || right.some((value) => value > 4)) return 0;
  const correlation = pearson(left, right);
  const leftConstant = left.every((value) => value === left[0]);
  const rightConstant = right.every((value) => value === right[0]);
  return leftConstant || rightConstant ? 1 : clamp(correlation + 1e-14, -0.999999999, 0.999999999);
}

// CalCR tests the direct six-cell relationship and five category relabelings.
// An inverse code is retained only when direct r is below 0.83 and a relabeled
// vector correlates better. Co-recombinant membership requires inverse code 0.
export function sourceCalCr(presumedVector, candidateVector, threshold = 0.83) {
  const permutations = [
    [0, 1, 2, 3, 4, 5],
    [1, 0, 2, 4, 3, 5],
    [2, 1, 0, 5, 4, 3],
    [0, 2, 1, 3, 5, 4],
    [1, 2, 0, 4, 5, 3],
    [2, 0, 1, 5, 3, 4],
  ];
  const correlations = permutations.map((permutation) => sourcePearsonSix(
    presumedVector,
    permutation.map((index) => candidateVector[index]),
  ));
  let r = correlations[0];
  let inversion = 0;
  if (r < threshold) {
    for (let index = 1; index < correlations.length; index += 1) {
      if (correlations[index] > r) {
        r = correlations[index];
        inversion = index === 5 ? 4 : index;
      }
    }
  }
  const direct = correlations[0];
  const degreesOfFreedom = 4;
  const t = Math.abs(direct) * Math.sqrt(degreesOfFreedom / Math.max(1e-12, 1 - direct * direct));
  return { r: direct, selectedR: r, inversion, pValue: studentTTwoSided(t, degreesOfFreedom), values: 6, permutations: correlations };
}

function distanceCorrelationForPair(encoded, length, presumed, parent1, parent2, candidate, pair, threshold) {
  const presumedVector = sourceCategoryVector(encoded, length, presumed, parent1, parent2, presumed, pair);
  const candidateVector = sourceCategoryVector(encoded, length, presumed, parent1, parent2, candidate, pair);
  return sourceCalCr(presumedVector, candidateVector, threshold);
}

function sourceSdmDistance(encoded, length, triplet, anchor, candidate, region) {
  let valid = 0;
  let differences = 0;
  visitRegion(region, (site) => {
    const first = encoded[triplet[0] * length + site];
    const second = encoded[triplet[1] * length + site];
    const third = encoded[triplet[2] * length + site];
    const anchorBase = encoded[anchor * length + site];
    const candidateBase = encoded[candidate * length + site];
    if (first >= 4 || second >= 4 || third >= 4 || anchorBase >= 4 || candidateBase >= 4) return;
    valid += 1;
    if (anchorBase !== candidateBase) differences += 1;
  });
  return valid > 0 ? differences / valid : 10;
}

// Direct translation of MakeProperRCorr's active SDM post-filter. When only
// one breakpoint pair correlates (or neither does and the combined pair does),
// a moderately high direct correlation is discarded if the candidate is
// farther from the presumed recombinant than from either presumed parent in
// the source comparison region.
function applySourceSdmFilters(encoded, length, focal, candidate, regionBundle, correlations) {
  const sourceThreshold = 0.83;
  const directSupport = correlations.slice(0, 2).filter((entry) => entry.r > sourceThreshold && entry.inversion === 0).length;
  const check = directSupport === 1
    || (directSupport === 0 && correlations[2].r > sourceThreshold && correlations[2].inversion === 0);
  if (!check) return correlations;
  const eventTriplet = triad(focal);
  return correlations.map((entry, index) => {
    if (!(entry.r > sourceThreshold && entry.r < 0.99 && entry.inversion === 0)) return entry;
    const region = regionBundle.sourcePairs[index].sdmRegion;
    const presumedDistance = sourceSdmDistance(encoded, length, eventTriplet, focal.recombinant, candidate, region);
    const firstParentDistance = sourceSdmDistance(encoded, length, eventTriplet, focal.majorParent, candidate, region);
    const secondParentDistance = sourceSdmDistance(encoded, length, eventTriplet, focal.minorParent, candidate, region);
    if (presumedDistance <= firstParentDistance && presumedDistance <= secondParentDistance) return entry;
    return {
      ...entry,
      unfilteredR: entry.r,
      r: 0,
      selectedR: 0,
      pValue: 1,
      sourceSdmFiltered: true,
      sourceSdmDistances: [presumedDistance, firstParentDistance, secondParentDistance],
    };
  });
}

function quartetSplit(encoded, length, candidate, presumed, parent1, parent2, region, cache) {
  const candidatePresumed = jukesCantorDistance(encoded, length, candidate, presumed, region, cache);
  const parentPair = jukesCantorDistance(encoded, length, parent1, parent2, region, cache);
  const candidateParent1 = jukesCantorDistance(encoded, length, candidate, parent1, region, cache);
  const candidateParent2 = jukesCantorDistance(encoded, length, candidate, parent2, region, cache);
  const presumedParent1 = jukesCantorDistance(encoded, length, presumed, parent1, region, cache);
  const presumedParent2 = jukesCantorDistance(encoded, length, presumed, parent2, region, cache);
  const sums = [
    candidatePresumed + parentPair,
    candidateParent1 + presumedParent2,
    candidateParent2 + presumedParent1,
  ];
  const order = [0, 1, 2].sort((left, right) => sums[left] - sums[right]);
  const scale = Math.max(1e-9, sums[order[1]] + sums[order[0]]);
  return { split: order[0], margin: (sums[order[1]] - sums[order[0]]) / scale };
}

function nearestParentMovement(encoded, length, sequence, parent1, parent2, region, cache) {
  const first = jukesCantorDistance(encoded, length, sequence, parent1, region, cache);
  const second = jukesCantorDistance(encoded, length, sequence, parent2, region, cache);
  return { parent: first <= second ? 0 : 1, margin: Math.abs(first - second) / Math.max(1e-9, first + second) };
}

function phylogeneticMovementForPair(encoded, length, presumed, parent1, parent2, candidate, regions, cache, minimumMargin) {
  const presumedLeft = nearestParentMovement(encoded, length, presumed, parent1, parent2, regions[0], cache);
  const presumedRight = nearestParentMovement(encoded, length, presumed, parent1, parent2, regions[1], cache);
  const candidateLeft = nearestParentMovement(encoded, length, candidate, parent1, parent2, regions[0], cache);
  const candidateRight = nearestParentMovement(encoded, length, candidate, parent1, parent2, regions[1], cache);
  const presumedMoves = presumedLeft.parent !== presumedRight.parent;
  const movesTogether = presumedMoves
    && candidateLeft.parent === presumedLeft.parent
    && candidateRight.parent === presumedRight.parent;
  const leftQuartet = quartetSplit(encoded, length, candidate, presumed, parent1, parent2, regions[0], cache);
  const rightQuartet = quartetSplit(encoded, length, candidate, presumed, parent1, parent2, regions[1], cache);
  const sisterTogether = leftQuartet.split === 0 && rightQuartet.split === 0;
  const margin = Math.min(
    presumedLeft.margin,
    presumedRight.margin,
    candidateLeft.margin,
    candidateRight.margin,
    Math.max(leftQuartet.margin, rightQuartet.margin),
  );
  return {
    supported: (movesTogether || sisterTogether) && margin >= minimumMargin,
    movesTogether,
    sisterTogether,
    margin,
  };
}

function sourceTreeComparisonScore(tree, presumed, parent1, parent2, candidate, collapsed) {
  const focalDistance = bootstrapTreeDistance(tree, presumed, candidate, collapsed);
  const presumedParents = [
    bootstrapTreeDistance(tree, presumed, parent1, collapsed),
    bootstrapTreeDistance(tree, presumed, parent2, collapsed),
  ];
  const candidateParents = [
    bootstrapTreeDistance(tree, candidate, parent1, collapsed),
    bootstrapTreeDistance(tree, candidate, parent2, collapsed),
  ];
  if (![focalDistance, ...presumedParents, ...candidateParents].every(Number.isFinite)) return null;
  const scoreComparison = (comparisons) => {
    if (comparisons.every((distance) => focalDistance < distance)) return 4;
    if (comparisons.every((distance) => focalDistance <= distance) && comparisons.some((distance) => focalDistance < distance)) return 2;
    if (comparisons.every((distance) => focalDistance > distance)) return -10;
    if (comparisons.some((distance) => focalDistance > distance)) return -2;
    return 0;
  };
  const comparisons = [...presumedParents, ...candidateParents];
  const scale = Math.max(1, focalDistance, ...comparisons);
  return {
    score: scoreComparison(presumedParents) + scoreComparison(candidateParents),
    close: comparisons.every((distance) => focalDistance <= distance)
      && comparisons.some((distance) => focalDistance < distance),
    margin: (Math.min(...comparisons) - focalDistance) / scale,
  };
}

function sisterSplitSupport(tree, presumed, parent1, parent2, candidate) {
  let maximum = 0;
  for (const split of tree?.splits ?? []) {
    const members = new Set(split.members);
    const focalSide = members.has(presumed) === members.has(candidate);
    const parentsOpposite = members.has(parent1) !== members.has(presumed)
      && members.has(parent2) !== members.has(presumed);
    if (focalSide && parentsOpposite) maximum = Math.max(maximum, split.support);
  }
  return maximum;
}

function treeParentChoice(tree, sequence, parent1, parent2, collapsed) {
  const first = bootstrapTreeDistance(tree, sequence, parent1, collapsed);
  const second = bootstrapTreeDistance(tree, sequence, parent2, collapsed);
  if (!Number.isFinite(first) || !Number.isFinite(second)) return { parent: -1, margin: 0 };
  if (first === second) return { parent: -1, margin: 0 };
  return {
    parent: first < second ? 0 : 1,
    margin: Math.abs(first - second) / Math.max(1, first + second),
  };
}

function treeMovesWithPresumed(treePair, presumed, parent1, parent2, candidate, collapsed) {
  const presumedPositions = treePair.map((tree) => treeParentChoice(tree, presumed, parent1, parent2, collapsed));
  const candidatePositions = treePair.map((tree) => treeParentChoice(tree, candidate, parent1, parent2, collapsed));
  const decisive = [...presumedPositions, ...candidatePositions].every((entry) => entry.parent >= 0);
  const supported = decisive
    && presumedPositions[0].parent !== presumedPositions[1].parent
    && candidatePositions[0].parent === presumedPositions[0].parent
    && candidatePositions[1].parent === presumedPositions[1].parent;
  return {
    supported,
    margin: supported ? Math.min(...presumedPositions.map((entry) => entry.margin), ...candidatePositions.map((entry) => entry.margin)) : 0,
  };
}

function bootstrapPhylogeneticMovementForPair(treePair, presumed, parent1, parent2, candidate, minimumMargin) {
  if (!treePair?.length || treePair.some((tree) => !tree?.index.has(candidate))) {
    return {
      supported: false,
      movesTogether: false,
      sisterTogether: false,
      margin: 0,
      sourceScore: 0,
      bootstrapSupport: 0,
      treeExcluded: true,
      bootstrapReplicates: treePair?.[0]?.replicates ?? 0,
      bootstrapCutoff: treePair?.[0]?.cutoff ?? 0.5,
    };
  }
  const collapsed = treePair.map((tree) => sourceTreeComparisonScore(tree, presumed, parent1, parent2, candidate, true));
  const uncollapsed = treePair.map((tree) => sourceTreeComparisonScore(tree, presumed, parent1, parent2, candidate, false));
  if ([...collapsed, ...uncollapsed].some((entry) => !entry)) {
    return { supported: false, movesTogether: false, sisterTogether: false, margin: 0, sourceScore: 0, bootstrapSupport: 0, treeExcluded: true };
  }
  const collapsedScore = collapsed.reduce((total, entry) => total + entry.score, 0);
  const uncollapsedScore = uncollapsed.reduce((total, entry) => total + entry.score, 0);
  const sourceScore = (collapsedScore + uncollapsedScore) / 2;
  const sisterTogether = collapsed.every((entry) => entry.close);
  const collapsedMovement = treeMovesWithPresumed(treePair, presumed, parent1, parent2, candidate, true);
  const uncollapsedMovement = treeMovesWithPresumed(treePair, presumed, parent1, parent2, candidate, false);
  const movesTogether = collapsedMovement.supported
    || uncollapsedMovement.supported
    || (collapsedScore > 0 && uncollapsedScore > 0);
  const margin = Math.max(
    0,
    collapsedMovement.margin,
    uncollapsedMovement.margin,
    ...collapsed.map((entry) => entry.margin),
    ...uncollapsed.map((entry) => entry.margin),
  );
  const bootstrapSupport = Math.min(...treePair.map((tree) => sisterSplitSupport(tree, presumed, parent1, parent2, candidate)));
  return {
    supported: (movesTogether || sisterTogether) && margin >= minimumMargin,
    movesTogether,
    sisterTogether,
    margin,
    sourceScore,
    collapsedScore,
    uncollapsedScore,
    collapsedMovement: collapsedMovement.supported,
    uncollapsedMovement: uncollapsedMovement.supported,
    bootstrapSupport,
    treeExcluded: false,
    bootstrapReplicates: treePair[0].replicates,
    bootstrapCutoff: treePair[0].cutoff,
    exactSiteBootstrap: treePair.every((tree) => tree.exactSiteBootstrap),
  };
}

function signalMatch(events, focal, candidateSequence, length, minimumOverlap) {
  const focalTriad = triad(focal);
  let best = null;
  for (const signal of events) {
    if (!triad(signal).includes(candidateSequence)) continue;
    const shared = sharedMembers(focalTriad, triad(signal));
    if (shared < 2) continue;
    const overlap = reciprocalTractOverlap(focal, signal, length);
    if (overlap.minimum <= minimumOverlap) continue;
    if (!best || overlap.minimum > best.overlap) best = { eventId: signal.id, sharedTripletMembers: shared, overlap: overlap.minimum };
  }
  return best;
}

function evaluateCandidate(events, focal, candidateSequence, regionBundle, treeBundle, encoded, length, options, cache) {
  const minimumMargin = Math.max(0, Math.min(1, options.clusterTopologyMargin ?? 0.005));
  const correlationAlpha = Math.max(1e-6, Math.min(0.5, options.clusterCorrelationAlpha ?? 0.05));
  const correlationR = Math.max(0, Math.min(0.999999, options.clusterCorrelationR ?? 0.83));
  const signalOverlap = Math.max(0.05, Math.min(1, options.clusterSignalOverlap ?? 0.3));
  const phylogeneticPairs = [];
  let correlationPairs = [];
  for (let pairIndex = 0; pairIndex < regionBundle.phylogeneticPairs.length; pairIndex += 1) {
    const regions = regionBundle.phylogeneticPairs[pairIndex];
    phylogeneticPairs.push(treeBundle
      ? bootstrapPhylogeneticMovementForPair(
          treeBundle.pairs[pairIndex],
          focal.recombinant,
          focal.majorParent,
          focal.minorParent,
          candidateSequence,
          minimumMargin,
        )
      : phylogeneticMovementForPair(
          encoded,
          length,
          focal.recombinant,
          focal.majorParent,
          focal.minorParent,
          candidateSequence,
          regions,
          cache,
          minimumMargin,
        ));
    correlationPairs.push(distanceCorrelationForPair(
      encoded,
      length,
      focal.recombinant,
      focal.majorParent,
      focal.minorParent,
      candidateSequence,
      regionBundle.sourcePairs[pairIndex],
      correlationR,
    ));
  }
  correlationPairs = applySourceSdmFilters(encoded, length, focal, candidateSequence, regionBundle, correlationPairs);
  const phylogenetic = phylogeneticPairs.some((pair) => pair.supported);
  const bestCorrelation = [...correlationPairs].sort((left, right) => left.pValue - right.pValue)[0];
  const distance = Boolean(bestCorrelation && bestCorrelation.inversion === 0 && bestCorrelation.pValue < correlationAlpha && bestCorrelation.r > correlationR);
  const detectableSignal = signalMatch(events, focal, candidateSequence, length, signalOverlap);
  const sets = Number(phylogenetic) + Number(distance) + Number(Boolean(detectableSignal));
  return {
    sequence: candidateSequence,
    supported: sets >= Math.max(1, Math.min(3, Math.trunc(options.clusterMinimumSets ?? 2))),
    sets,
    phylogenetic,
    distance,
    detectableSignal: Boolean(detectableSignal),
    bestCorrelation: bestCorrelation ? { r: bestCorrelation.r, pValue: bestCorrelation.pValue, inversion: bestCorrelation.inversion } : { r: 0, pValue: 1, inversion: 0 },
    topologyMargin: Math.max(0, ...phylogeneticPairs.map((pair) => pair.margin)),
    treeBootstrap: treeBundle ? {
      replicates: Math.max(0, ...phylogeneticPairs.map((pair) => pair.bootstrapReplicates ?? 0)),
      cutoff: phylogeneticPairs[0]?.bootstrapCutoff ?? 0.5,
      cohortTaxa: treeBundle.taxa.length,
      sourceSequenceCount: treeBundle.sourceSequenceCount,
      included: phylogeneticPairs.some((pair) => pair.treeExcluded !== true),
      exactSiteBootstrap: phylogeneticPairs.every((pair) => pair.exactSiteBootstrap !== false),
      sourceScore: Math.max(...phylogeneticPairs.map((pair) => pair.sourceScore ?? 0)),
      cohortCount: treeBundle.cohortCount ?? 1,
      candidateComplete: treeBundle.candidateComplete !== false,
    } : undefined,
    signal: detectableSignal,
    regionEvidence: phylogeneticPairs.map((pair, index) => ({
      pair: index === 0 ? "5-prime breakpoint" : index === 1 ? "3-prime breakpoint" : "tract/background",
      phylogenetic: pair.supported,
      movesTogether: pair.movesTogether,
      sisterTogether: pair.sisterTogether,
      topologyMargin: pair.margin,
      treeSourceScore: pair.sourceScore,
      bootstrapSupport: pair.bootstrapSupport,
      bootstrapReplicates: pair.bootstrapReplicates,
      bootstrapCutoff: pair.bootstrapCutoff,
      treeExcluded: pair.treeExcluded,
      correlationR: correlationPairs[index].r,
      correlationP: correlationPairs[index].pValue,
      correlationInversion: correlationPairs[index].inversion,
      correlationPermutations: correlationPairs[index].permutations,
      correlationSdmFiltered: correlationPairs[index].sourceSdmFiltered === true,
    })),
  };
}

function rawPDistance(encoded, length, firstSequence, secondSequence, regions) {
  let valid = 0;
  let differences = 0;
  const firstOffset = firstSequence * length;
  const secondOffset = secondSequence * length;
  for (const region of regions) visitRegion(region, (site) => {
    const first = encoded[firstOffset + site];
    const second = encoded[secondOffset + site];
    if (first >= 4 || second >= 4) return;
    valid += 1;
    if (first !== second) differences += 1;
  });
  return valid ? differences / valid : 1;
}

function overlapSegments(left, right, length) {
  const segments = [];
  for (const [leftStart, leftEnd] of tractSegments(left, length)) {
    for (const [rightStart, rightEnd] of tractSegments(right, length)) {
      const start = Math.max(leftStart, rightStart);
      const end = Math.min(leftEnd, rightEnd);
      if (end > start) segments.push([start, end]);
    }
  }
  return segments;
}

function circularWindow(center, radius, length) {
  if (radius * 2 >= length) return [{ kind: "segments", segments: [[0, length]] }];
  const start = center - radius;
  const end = center + radius;
  if (start < 0) return [{ kind: "segments", segments: [[0, end], [length + start, length]] }];
  if (end > length) return [{ kind: "segments", segments: [[start, length], [0, end - length]] }];
  return [{ kind: "segments", segments: [[start, end]] }];
}

// Direct port of the event-similarity score used by RDP5 GetSupers: average
// normalised daughter divergence around the overlap boundaries and the lack of
// tract overlap; values below 0.1 are merged agglomeratively in the source.
function sourceSuperEventDistance(left, right, encoded, length, globalMaximumDistance) {
  const sharedSegments = overlapSegments(left, right, length);
  if (!sharedSegments.length || left.recombinant === right.recombinant) return Infinity;
  const sharedBases = segmentLength(sharedSegments);
  const leftLength = segmentLength(tractSegments(left, length));
  const rightLength = segmentLength(tractSegments(right, length));
  const overlapScore = 1 - sharedBases / Math.max(1, (leftLength + rightLength) / 2);
  const first = sharedSegments[0][0];
  const last = sharedSegments.at(-1)[1];
  const boundaryRegions = [
    ...circularWindow(first, 100, length),
    ...circularWindow(last, 100, length),
  ];
  const daughterDistance = rawPDistance(encoded, length, left.recombinant, right.recombinant, boundaryRegions);
  const normalisedDaughterDistance = daughterDistance / Math.max(1e-9, globalMaximumDistance);
  return (normalisedDaughterDistance + overlapScore) / 2;
}

function pairConfidence(left, right, overlap, sourceDistance) {
  const setStrength = (left.sets + right.sets) / 6;
  const correlation = Math.max(0, Math.min(1, (Math.max(left.bestCorrelation.r, right.bestCorrelation.r) + 1) / 2));
  const topology = Math.min(1, Math.max(left.topologyMargin, right.topologyMargin) * 10);
  const sourceSimilarity = Number.isFinite(sourceDistance) ? Math.max(0, 1 - sourceDistance / 0.2) : 0;
  return Math.max(0, Math.min(1, 0.4 * setStrength + 0.2 * overlap.minimum + 0.15 * correlation + 0.1 * topology + 0.15 * sourceSimilarity));
}

function unionFind(size) {
  const parent = Int32Array.from({ length: size }, (_, index) => index);
  const find = (value) => {
    let root = value;
    while (parent[root] !== root) root = parent[root];
    while (parent[value] !== value) {
      const next = parent[value];
      parent[value] = root;
      value = next;
    }
    return root;
  };
  const join = (left, right) => {
    const leftRoot = find(left);
    const rightRoot = find(right);
    if (leftRoot !== rightRoot) parent[rightRoot] = leftRoot;
  };
  return { find, join };
}

// Source-parity translation of Module2.bas GetSupers. The desktop routine is
// WPGMA-like rather than ordinary connected-components clustering: after the
// closest pair is merged, its distance to every remaining cluster is the
// unweighted mean of the two previous cluster distances. That distinction is
// important because it prevents a long chain of weak pairwise matches from
// automatically becoming one event.
export function sourceWeightedEventClusters(distanceMatrix, threshold = 0.1) {
  const size = distanceMatrix.length;
  if (size < 2) return Array.from({ length: size }, (_, index) => [index]);
  const distances = distanceMatrix.map((row, left) => Float64Array.from(
    { length: size },
    (_, right) => left === right ? Infinity : Number.isFinite(row[right]) ? row[right] : 100,
  ));
  const active = new Uint8Array(size);
  active.fill(1);
  const members = Array.from({ length: size }, (_, index) => [index]);
  while (true) {
    let best = threshold;
    let bestLeft = -1;
    let bestRight = -1;
    for (let left = 0; left < size - 1; left += 1) {
      if (!active[left]) continue;
      for (let right = left + 1; right < size; right += 1) {
        if (!active[right]) continue;
        const distance = distances[left][right];
        if (distance < best) {
          best = distance;
          bestLeft = left;
          bestRight = right;
        }
      }
    }
    if (bestLeft < 0) break;
    members[bestLeft].push(...members[bestRight]);
    members[bestRight] = [];
    for (let other = 0; other < size; other += 1) {
      if (!active[other] || other === bestLeft || other === bestRight) continue;
      const averaged = (distances[bestLeft][other] + distances[bestRight][other]) / 2;
      distances[bestLeft][other] = averaged;
      distances[other][bestLeft] = averaged;
      distances[bestRight][other] = 100;
      distances[other][bestRight] = 100;
    }
    active[bestRight] = 0;
  }
  return members.filter((group) => group.length > 0);
}

function stableEventOrder(left, right) {
  const leftP = Math.min(...(left.evidence ?? []).map((item) => item.correctedP), 1);
  const rightP = Math.min(...(right.evidence ?? []).map((item) => item.correctedP), 1);
  return leftP - rightP || right.chiSquare - left.chiSquare || left.start - right.start || left.id.localeCompare(right.id);
}

export function inferAncestralEventClusters(events, encoded, sequenceCount, length, options = {}, sourceMaximumDistance = null) {
  const output = events.map((event) => ({
    ...event,
    groupId: null,
    ancestralCluster: undefined,
    coRecombinantSets: undefined,
    recombinantIdentification: undefined,
  }));
  if (output.length === 0) return { events: output, clusters: [], pairwise: [] };
  if (options.ancestralClustering === false) {
    for (const event of output) {
      event.recombinantIdentification = identifyRecombinantRoles(
        event,
        encoded,
        sequenceCount,
        length,
        null,
        options,
      );
    }
    return { events: output, clusters: [], pairwise: [] };
  }
  const regionCache = new Map();
  const distanceCache = new Map();
  const regionBundles = output.map((event, eventIndex) => {
    const bundle = makeEventRegions(event, encoded, length, options);
    let regionIndex = 0;
    const regions = new Set(bundle.phylogeneticPairs.flat());
    for (const region of regions) region.id = `${eventIndex}:${regionIndex++}`;
    regionCache.set(event.id, bundle);
    return bundle;
  });
  const treeCohorts = sequenceCount > 3
    ? output.map((event, eventIndex) => buildEventBootstrapTreeCohorts(
        encoded,
        length,
        sequenceCount,
        triad(event),
        regionBundles[eventIndex].phylogeneticPairs,
        options,
        eventIndex,
      ))
    : output.map(() => null);
  const treeBundles = treeCohorts.map((cohorts) => cohorts?.primary ?? null);
  const treeBundleForCandidate = (eventIndex, candidate) => {
    const cohorts = treeCohorts[eventIndex];
    if (!cohorts) return null;
    const bundle = cohorts.byTaxon.get(candidate) ?? cohorts.primary;
    return {
      ...bundle,
      cohortCount: cohorts.cohorts.length,
      candidateComplete: cohorts.candidateComplete,
    };
  };

  // RDP5 section 4.1.4 deliberately repeats the co-recombinant search with
  // each member of the detecting triplet treated as the presumed recombinant.
  // Retaining all three sets is important: later role reassignment can select
  // the matching set without rerunning the whole analysis, and descendants
  // without their own retained signal remain visible to the analyst.
  const orientationEvidence = output.map(() => new Map());
  output.forEach((event, eventIndex) => {
    const members = triad(event);
    const sets = members.map((presumedRecombinant, orientation) => {
      const parents = members.filter((_, index) => index !== orientation);
      const focal = {
        ...event,
        recombinant: presumedRecombinant,
        majorParent: parents[0],
        minorParent: parents[1],
      };
      const evidence = [];
      const evaluated = new Map();
      const memberSet = new Set([presumedRecombinant]);
      for (let candidate = 0; candidate < sequenceCount; candidate += 1) {
        if (members.includes(candidate)) continue;
        const result = evaluateCandidate(
          output,
          focal,
          candidate,
          regionBundles[eventIndex],
          treeBundleForCandidate(eventIndex, candidate),
          encoded,
          length,
          options,
          distanceCache,
        );
        evaluated.set(candidate, result);
        if (!result.supported) continue;
        memberSet.add(candidate);
        evidence.push({
          sequence: result.sequence,
          sets: result.sets,
          phylogenetic: result.phylogenetic,
          distance: result.distance,
          detectableSignal: result.detectableSignal,
          bestCorrelation: result.bestCorrelation,
          topologyMargin: result.topologyMargin,
          treeBootstrap: result.treeBootstrap,
          regionEvidence: result.regionEvidence,
        });
      }
      const set = {
        presumedRecombinant,
        parents,
        sequenceMembers: [...memberSet].sort((left, right) => left - right),
        testedSequences: Math.max(0, sequenceCount - members.length),
        requiredEvidenceSets: Math.max(1, Math.min(3, Math.trunc(options.clusterMinimumSets ?? 2))),
        evidence,
      };
      orientationEvidence[eventIndex].set(presumedRecombinant, evaluated);
      return set;
    });
    event.coRecombinantSets = sets;
  });
  output.forEach((event, eventIndex) => {
    const identification = identifyRecombinantRoles(
      event,
      encoded,
      sequenceCount,
      length,
      treeBundles[eventIndex],
      { ...options, roleEventCorpus: output },
    );
    event.recombinantIdentification = identification;
    if (!identification) return;
    if (identification.ambiguous) {
      const warning = `RDP5 source profile consensus is role-ambiguous (${Math.round(identification.confidence * 100)}% relative support); inspect all three recombinant polarities.`;
      if (!event.warnings.includes(warning)) event.warnings = [...event.warnings, warning];
    } else if (identification.recommended !== event.recombinant) {
      const warning = "RDP5 source profile consensus challenges the current recombinant assignment; inspect and, if appropriate, apply the highest-scoring polarity.";
      if (!event.warnings.includes(warning)) event.warnings = [...event.warnings, warning];
    }
  });
  const pairwise = [];
  const fullRegion = { kind: "segments", segments: [[0, length]] };
  let globalMaximumDistance = Number.isFinite(sourceMaximumDistance) ? Number(sourceMaximumDistance) : 0;
  if (!(globalMaximumDistance > 0)) {
    for (let left = 0; left < sequenceCount; left += 1) {
      for (let right = left + 1; right < sequenceCount; right += 1) {
        globalMaximumDistance = Math.max(globalMaximumDistance, rawPDistance(
          encoded,
          length,
          left,
          right,
          [fullRegion],
        ));
      }
    }
  }
  if (!(globalMaximumDistance > 0)) globalMaximumDistance = 1;
  const minimumSignalOverlap = Math.max(0.05, Math.min(1, options.clusterSignalOverlap ?? 0.3));
  const sourceThreshold = Math.max(0, Math.min(1, options.clusterSourceSimilarity ?? 0.1));
  const reciprocal = options.clusterReciprocal === true;
  const sourceDistances = Array.from({ length: output.length }, () => new Float64Array(output.length).fill(100));
  for (let index = 0; index < output.length; index += 1) sourceDistances[index][index] = 0;
  for (let leftIndex = 0; leftIndex < output.length; leftIndex += 1) {
    const leftEvent = output[leftIndex];
    for (let rightIndex = leftIndex + 1; rightIndex < output.length; rightIndex += 1) {
      const rightEvent = output[rightIndex];
      if (leftEvent.recombinant === rightEvent.recombinant) continue;
      const overlap = reciprocalTractOverlap(leftEvent, rightEvent, length);
      // A partial overprint may erase one old breakpoint, so the manual's 30%
      // signal-overlap rule is applied to the shorter tract, not both tracts.
      if (overlap.minimum <= minimumSignalOverlap) continue;
      const leftToRight = orientationEvidence[leftIndex].get(leftEvent.recombinant)?.get(rightEvent.recombinant)
        ?? evaluateCandidate(output, leftEvent, rightEvent.recombinant, regionBundles[leftIndex], treeBundleForCandidate(leftIndex, rightEvent.recombinant), encoded, length, options, distanceCache);
      const rightToLeft = orientationEvidence[rightIndex].get(rightEvent.recombinant)?.get(leftEvent.recombinant)
        ?? evaluateCandidate(output, rightEvent, leftEvent.recombinant, regionBundles[rightIndex], treeBundleForCandidate(rightIndex, leftEvent.recombinant), encoded, length, options, distanceCache);
      const threeSetSupported = reciprocal
        ? leftToRight.supported && rightToLeft.supported
        : leftToRight.supported || rightToLeft.supported;
      const sourceDistance = sourceSuperEventDistance(leftEvent, rightEvent, encoded, length, globalMaximumDistance);
      sourceDistances[leftIndex][rightIndex] = sourceDistance;
      sourceDistances[rightIndex][leftIndex] = sourceDistance;
      const sourceSupported = sourceDistance < sourceThreshold;
      const supported = threeSetSupported || sourceSupported;
      const confidence = pairConfidence(leftToRight, rightToLeft, overlap, sourceDistance);
      pairwise.push({
        leftIndex,
        rightIndex,
        supported,
        confidence,
        overlap,
        partialOverprint: overlap.left < 0.8 || overlap.right < 0.8,
        threeSetSupported,
        sourceSupported,
        sourceDistance,
        leftToRight,
        rightToLeft,
      });
    }
  }

  const sets = unionFind(output.length);
  const minimumConfidence = Math.max(0, Math.min(1, options.clusterMinimumConfidence ?? 0.55));
  const sourceComponents = sourceWeightedEventClusters(sourceDistances, sourceThreshold);
  for (const component of sourceComponents) {
    for (let index = 1; index < component.length; index += 1) sets.join(component[0], component[index]);
  }
  const links = pairwise
    .filter((pair) => pair.threeSetSupported && pair.confidence >= minimumConfidence)
    .sort((left, right) => right.confidence - left.confidence);
  for (const link of links) {
    const leftRoot = sets.find(link.leftIndex);
    const rightRoot = sets.find(link.rightIndex);
    if (leftRoot === rightRoot) continue;
    const leftMembers = output.map((_, index) => index).filter((index) => sets.find(index) === leftRoot);
    const rightMembers = output.map((_, index) => index).filter((index) => sets.find(index) === rightRoot);
    // A co-recombinant assignment is allowed to expand a source-merged event,
    // but every existing member on the opposite side must have either the
    // three-set support or a source-distance relationship. This preserves the
    // source's event merge while avoiding an unrelated one-link bridge.
    const crossLinks = leftMembers.flatMap((left) => rightMembers.map((right) => pairwise.find((pair) =>
      (pair.leftIndex === left && pair.rightIndex === right) || (pair.leftIndex === right && pair.rightIndex === left),
    ))).filter(Boolean);
    const completeLink = crossLinks.every((pair) => pair.threeSetSupported || pair.sourceSupported);
    if (completeLink || leftMembers.length === 1 || rightMembers.length === 1) sets.join(leftRoot, rightRoot);
  }

  const components = new Map();
  output.forEach((_, index) => {
    const root = sets.find(index);
    components.set(root, [...(components.get(root) ?? []), index]);
  });
  const clusters = [...components.values()]
    .filter((members) => members.length > 1)
    .sort((left, right) => stableEventOrder(output[left[0]], output[right[0]]));
  clusters.forEach((members, clusterIndex) => {
    members.sort((left, right) => stableEventOrder(output[left], output[right]));
    const groupId = `ancestry-${String(clusterIndex + 1).padStart(3, "0")}`;
    const memberSet = new Set(members);
    const clusterLinks = pairwise.filter((pair) => memberSet.has(pair.leftIndex) && memberSet.has(pair.rightIndex) && pair.supported);
    const confidence = clusterLinks.length
      ? clusterLinks.reduce((total, pair) => total + pair.confidence, 0) / clusterLinks.length
      : 0;
    const representative = members[0];
    const sequenceMembers = [...new Set(members.flatMap((index) => {
      const event = output[index];
      const sourceSet = event.coRecombinantSets?.find((set) => set.presumedRecombinant === event.recombinant);
      return sourceSet?.sequenceMembers ?? [event.recombinant];
    }))].sort((left, right) => left - right);
    const evidenceCounts = clusterLinks.reduce((counts, pair) => {
      for (const direction of [pair.leftToRight, pair.rightToLeft]) {
        if (direction.phylogenetic) counts.phylogenetic += 1;
        if (direction.distance) counts.distance += 1;
        if (direction.detectableSignal) counts.detectableSignal += 1;
      }
      if (pair.sourceSupported) counts.sourceSimilarity += 1;
      return counts;
    }, { phylogenetic: 0, distance: 0, detectableSignal: 0, sourceSimilarity: 0 });
    for (const member of members) {
      const event = output[member];
      event.groupId = groupId;
      event.ancestralCluster = {
        inference: "rdp5-three-set",
        representativeId: output[representative].id,
        memberEventIds: members.map((index) => output[index].id),
        sequenceMembers,
        confidence,
        evidenceCounts,
        sourceMerge: {
          threshold: sourceThreshold,
          pairDistances: clusterLinks.filter((pair) => Number.isFinite(pair.sourceDistance)).map((pair) => ({
            eventIds: [output[pair.leftIndex].id, output[pair.rightIndex].id],
            distance: pair.sourceDistance,
            belowThreshold: pair.sourceSupported,
          })),
        },
        partialOverprint: clusterLinks.some((pair) => pair.partialOverprint),
        pairwise: clusterLinks.map((pair) => ({
          eventIds: [output[pair.leftIndex].id, output[pair.rightIndex].id],
          confidence: pair.confidence,
          overlap: pair.overlap,
          leftToRight: pair.leftToRight,
          rightToLeft: pair.rightToLeft,
        })),
      };
      const warning = event.ancestralCluster.partialOverprint
        ? `Inferred ${groupId} includes partially overprinted evidence; verify the shared ancestry in tract/background and breakpoint-flank trees.`
        : `Inferred ${groupId} by the RDP5 three-set co-recombinant rule; verify the common ancestral event in local trees.`;
      if (!event.warnings.includes(warning)) event.warnings = [...event.warnings, warning];
    }
  });
  return {
    events: output,
    clusters: clusters.map((members, index) => ({
      id: `ancestry-${String(index + 1).padStart(3, "0")}`,
      eventIndexes: members,
      eventIds: members.map((member) => output[member].id),
      sequenceMembers: [...new Set(members.flatMap((member) => {
        const event = output[member];
        const sourceSet = event.coRecombinantSets?.find((set) => set.presumedRecombinant === event.recombinant);
        return sourceSet?.sequenceMembers ?? [event.recombinant];
      }))].sort((left, right) => left - right),
    })),
    pairwise,
  };
}
