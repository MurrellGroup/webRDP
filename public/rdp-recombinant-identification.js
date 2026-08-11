// Source-guided recombinant identification for an RDP5 event triplet.
//
// The profile calculations below are direct translations of the active
// MakePhPrScore and MakeTrpGroups/MakeTrpScore routines in the author-supplied
// RDP5 sources. The desktop-default MakeConsensusC path is represented as 18
// standalone source statistics, its FinalTrim parsimony penalty, and six
// joint rules. The only unported selector family is named in the returned
// provenance instead of being silently replaced by a generic score.

const IMPLEMENTED_COMPONENTS = [
  "PhPr",
  "TreePhPr",
  "SubPhPr",
  "TreeSubPhPr",
  "SubDist",
  "TreeSubDist",
  "TrpScore",
  "OuCheck",
  "O:E",
  "O:EDist",
  "dMax (VisRD)",
  "ParsimonyO",
  "ParsimonyI",
  "Conflict",
  "SSDist",
  "OUIndex",
  "SetDistT",
  "SetDistP",
  "FinalTrim parsimony penalty",
];

const PENDING_COMPONENTS = [
  "optional logistic / neural-network selector",
];

// VB6 CLng uses round-to-nearest-even, not JavaScript's half-toward-+Infinity
// Math.round rule. MakeConsensusC quantizes several inputs through CLng before
// applying exact-equality/sentinel tests, so this is observable behavior rather
// than display formatting.
export function sourceVbClng(value) {
  if (!Number.isFinite(value)) return value;
  const lower = Math.floor(value);
  const fraction = value - lower;
  if (fraction < 0.5) return lower;
  if (fraction > 0.5) return lower + 1;
  return Math.abs(lower % 2) === 0 ? lower : lower + 1;
}

function sourceQuantize(value, scale) {
  return sourceVbClng(value * scale) / scale;
}

function sourceProfileConsensusInputs(scores, mode) {
  if (mode === "raw") {
    return {
      phPr: scores.phPr.map((value) => sourceQuantize(value, 100_000)),
      // Desktop SubPhPrScore (manual SubDist).
      subDist: scores.subDist.map((value) => sourceQuantize(value, 10_000)),
      // Desktop SubScore (manual SubPhPr); large sentinels are left intact.
      subPhPr: scores.subPhPr.map((value) => Math.abs(value) < 100 ? sourceQuantize(value, 1_000_000) : value),
    };
  }
  return {
    phPr: scores.phPr.map((value) => sourceQuantize(value, 100_000)),
    // Desktop SubPhPrScore2 and SubScore2.
    subDist: scores.subDist.map((value) => sourceQuantize(value, 100_000)),
    subPhPr: scores.subPhPr.map((value) => sourceQuantize(Math.min(value, 100_000), 100_000)),
  };
}

function sourceProfileGate(phPr) {
  // PS1/PS2/PS3 are initialized to three in the source; their only effective
  // disabling branch is a post-quantization absolute sentinel of exactly one.
  return phPr.every((value) => Math.abs(value) !== 1);
}

function tractSegments(event, length) {
  if (event.wraps && event.start > event.end) {
    return [[event.start, length], ...(event.end > 0 ? [[0, event.end]] : [])];
  }
  return event.end > event.start ? [[event.start, event.end]] : [];
}

function backgroundSegments(event, length) {
  if (event.wraps && event.start > event.end) {
    return event.end < event.start ? [[event.end, event.start]] : [];
  }
  return [[0, event.start], ...(event.end < length ? [[event.end, length]] : [])]
    .filter(([start, end]) => end > start);
}

function quartetPatternScores(first, second, third, fourth) {
  if (first >= 4 || second >= 4 || third >= 4 || fourth >= 4) return [0, 0, 0];
  const split = (firstA, firstB, secondA, secondB) => {
    if (firstA === firstB && secondA === secondB && firstA !== secondA) return 2;
    if (firstA === firstB && firstA !== secondA && firstA !== secondB && secondA !== secondB) return 1;
    if (secondA === secondB && secondA !== firstA && secondA !== firstB && firstA !== firstB) return 1;
    return 0;
  };
  return [
    split(first, second, third, fourth),
    split(first, third, second, fourth),
    split(first, fourth, second, third),
  ];
}

// Scalar reference/fallback for CalcMaxD + CMaxD2P3. The analysis worker uses
// the bit-packed WASM implementation of exactly this loop; retaining the
// readable implementation here makes the source behavior testable and keeps
// role identification functional if a browser cannot instantiate WASM.
export function sourceDmaxScores(encoded, length, candidates, cohort, event) {
  const insideSites = new Uint8Array(length);
  for (const [start, end] of tractSegments(event, length)) insideSites.fill(1, start, end);
  const sums = [0, 0, 0];
  const counts = [0, 0, 0];
  for (let ai = 0; ai < cohort.length - 3; ai += 1) {
    for (let bi = ai + 1; bi < cohort.length - 2; bi += 1) {
      for (let ci = bi + 1; ci < cohort.length - 1; ci += 1) {
        for (let di = ci + 1; di < cohort.length; di += 1) {
          const quartet = [cohort[ai], cohort[bi], cohort[ci], cohort[di]];
          const includes = candidates.map((candidate) => quartet.includes(candidate));
          if (!includes.some(Boolean)) continue;
          const inside = [0, 0, 0];
          const outside = [0, 0, 0];
          for (let site = 0; site < length; site += 1) {
            const scores = quartetPatternScores(
              encoded[quartet[0] * length + site],
              encoded[quartet[1] * length + site],
              encoded[quartet[2] * length + site],
              encoded[quartet[3] * length + site],
            );
            const target = insideSites[site] ? inside : outside;
            target[0] += scores[0];
            target[1] += scores[1];
            target[2] += scores[2];
          }
          const insideTotal = inside[0] + inside[1] + inside[2];
          const outsideTotal = outside[0] + outside[1] + outside[2];
          if (!(insideTotal > 0) || !(outsideTotal > 0)) continue;
          const distance = Math.abs(inside[0] / insideTotal - outside[0] / outsideTotal)
            + Math.abs(inside[1] / insideTotal - outside[1] / outsideTotal)
            + Math.abs(inside[2] / insideTotal - outside[2] / outsideTotal);
          includes.forEach((included, index) => {
            if (!included) return;
            sums[index] += distance;
            counts[index] += 1;
          });
        }
      }
    }
  }
  return {
    values: sums.map((sum, index) => counts[index] > 0 ? sum / counts[index] : 0),
    quartetCounts: counts,
    cohortSize: cohort.length,
    sourceRoutine: "CalcMaxD + CMaxD2P3",
    wasmAccelerated: false,
  };
}

function jcDistance(encoded, length, first, second, segments) {
  const firstOffset = first * length;
  const secondOffset = second * length;
  let valid = 0;
  let differences = 0;
  for (const [start, end] of segments) {
    for (let site = start; site < end; site += 1) {
      const left = encoded[firstOffset + site];
      const right = encoded[secondOffset + site];
      if (left >= 4 || right >= 4) continue;
      valid += 1;
      if (left !== right) differences += 1;
    }
  }
  if (valid === 0) return 3;
  const p = differences / valid;
  return p >= 0.749999 ? 3 : Math.max(0, -0.75 * Math.log(1 - 4 * p / 3));
}

function deterministicCohort(sequenceCount, required, maximum) {
  if (sequenceCount <= maximum) return Array.from({ length: sequenceCount }, (_, index) => index);
  const selected = new Set(required);
  const available = Math.max(0, maximum - selected.size);
  for (let draw = 0; draw < available; draw += 1) {
    const index = available <= 1
      ? Math.floor((sequenceCount - 1) / 2)
      : Math.floor(draw * (sequenceCount - 1) / (available - 1));
    selected.add(index);
  }
  if (selected.size < maximum) {
    for (let index = 0; index < sequenceCount && selected.size < maximum; index += 1) selected.add(index);
  }
  return [...selected].sort((left, right) => left - right);
}

function rawProfiles(encoded, length, candidates, cohort, event) {
  const background = backgroundSegments(event, length);
  const tract = tractSegments(event, length);
  const first = new Map();
  const second = new Map();
  for (const candidate of candidates) {
    first.set(candidate, Float64Array.from(cohort, (taxon) => (
      taxon === candidate ? 0 : jcDistance(encoded, length, candidate, taxon, background)
    )));
    second.set(candidate, Float64Array.from(cohort, (taxon) => (
      taxon === candidate ? 0 : jcDistance(encoded, length, candidate, taxon, tract)
    )));
  }
  return { first, second };
}

function treeProfiles(treePair, candidates, cohort, collapsed = true) {
  if (!treePair?.[0] || !treePair?.[1]) return null;
  const first = new Map();
  const second = new Map();
  for (const candidate of candidates) {
    const leftIndex = treePair[0].index.get(candidate);
    const rightIndex = treePair[1].index.get(candidate);
    if (leftIndex === undefined || rightIndex === undefined) return null;
    const leftMatrix = collapsed ? treePair[0].collapsed : treePair[0].uncollapsed;
    const rightMatrix = collapsed ? treePair[1].collapsed : treePair[1].uncollapsed;
    first.set(candidate, Float64Array.from(cohort, (taxon) => {
      const index = treePair[0].index.get(taxon);
      return index === undefined ? Number.NaN : leftMatrix[leftIndex * treePair[0].taxa.length + index];
    }));
    second.set(candidate, Float64Array.from(cohort, (taxon) => {
      const index = treePair[1].index.get(taxon);
      return index === undefined ? Number.NaN : rightMatrix[rightIndex * treePair[1].taxa.length + index];
    }));
  }
  return { first, second };
}

// RDP5 MakePhPrScore uses this algebra rather than a library correlation.
// Constant profiles deliberately return 1, matching the source's sentinel.
function sourcePearson(left, right) {
  const count = Math.min(left.length, right.length);
  if (count <= 1) return 1;
  let sumX = 0;
  let sumY = 0;
  let sumXY = 0;
  let sumX2 = 0;
  let sumY2 = 0;
  for (let index = 0; index < count; index += 1) {
    const x = left[index];
    const y = right[index];
    sumX += x;
    sumY += y;
    sumXY += x * y;
    sumX2 += x * x;
    sumY2 += y * y;
  }
  if (!(sumX2 > 0) || !(sumY2 > 0)) return 1;
  const xVariance = count * sumX2 - sumX * sumX;
  const yVariance = count * sumY2 - sumY * sumY;
  if (!(xVariance > 1e-12) || !(yVariance > 1e-12)) return 1;
  return Math.max(-1, Math.min(1, (count * sumXY - sumX * sumY) / Math.sqrt(xVariance * yVariance)));
}

// Direct MakePhPrScore data flow.  In the source SubPhPrScore is the average
// correlation after removing a candidate (manual: SubDist), while SubScore is
// the distance-change sum (manual: SubPhPr).
export function sourcePhPrScores(candidates, cohort, profiles) {
  const phPr = [1, 1, 1];
  const subDist = [1, 1, 1];
  const subPhPr = [0, 0, 0];
  const cohortPositions = new Map(cohort.map((taxon, index) => [taxon, index]));
  for (let candidateIndex = 0; candidateIndex < candidates.length; candidateIndex += 1) {
    const candidate = candidates[candidateIndex];
    const left = profiles.first.get(candidate);
    const right = profiles.second.get(candidate);
    const valuesLeft = [];
    const valuesRight = [];
    for (let position = 0; position < cohort.length; position += 1) {
      if (cohort[position] === candidate) continue;
      if (!Number.isFinite(left[position]) || !Number.isFinite(right[position])) continue;
      valuesLeft.push(left[position]);
      valuesRight.push(right[position]);
      subPhPr[candidateIndex] += Math.abs(left[position] - right[position]);
    }
    phPr[candidateIndex] = sourcePearson(valuesLeft, valuesRight);
  }
  for (let removedIndex = 0; removedIndex < candidates.length; removedIndex += 1) {
    const correlations = [];
    for (let focalIndex = 0; focalIndex < candidates.length; focalIndex += 1) {
      if (focalIndex === removedIndex) continue;
      const focal = candidates[focalIndex];
      const left = profiles.first.get(focal);
      const right = profiles.second.get(focal);
      const valuesLeft = [];
      const valuesRight = [];
      for (const taxon of cohort) {
        if (taxon === focal || taxon === candidates[removedIndex]) continue;
        const position = cohortPositions.get(taxon);
        if (position === undefined || !Number.isFinite(left[position]) || !Number.isFinite(right[position])) continue;
        valuesLeft.push(left[position]);
        valuesRight.push(right[position]);
      }
      correlations.push(sourcePearson(valuesLeft, valuesRight));
    }
    subDist[removedIndex] = correlations.length
      ? correlations.reduce((total, value) => total + value, 0) / correlations.length
      : 1;
  }
  return { phPr, subDist, subPhPr };
}

function comparisonChanged(firstLeft, firstRight, secondLeft, secondRight) {
  return (firstLeft > firstRight && secondLeft < secondRight)
    || (firstLeft < firstRight && secondLeft > secondRight)
    || (firstLeft === firstRight && secondLeft !== secondRight)
    || (firstLeft !== firstRight && secondLeft === secondRight);
}

// Direct MakeTrpGroups + MakeTrpScore behavior. Exact equal-distance groups
// are important after low-bootstrap branches have been collapsed.
export function sourceTrpScores(candidates, cohort, profiles) {
  const scores = [0, 0, 0];
  for (let candidateIndex = 0; candidateIndex < candidates.length; candidateIndex += 1) {
    const focal = candidates[candidateIndex];
    const otherCandidates = candidates.filter((candidate) => candidate !== focal);
    const first = profiles.first.get(focal);
    const second = profiles.second.get(focal);
    const position = new Map(cohort.map((taxon, index) => [taxon, index]));
    const otherPositions = otherCandidates.map((candidate) => position.get(candidate));
    if (otherPositions.some((index) => index === undefined)) continue;
    const nearThreshold = Math.min(first[otherPositions[0]], first[otherPositions[1]]);
    let nearMaximum = 0;
    const groups = new Int32Array(cohort.length).fill(-1);
    for (let index = 0; index < cohort.length; index += 1) {
      if (first[index] < nearThreshold) {
        groups[index] = 0;
        nearMaximum = Math.max(nearMaximum, first[index]);
      }
    }
    let group = 1;
    const distinct = [...new Set(Array.from(first).filter((value, index) => groups[index] < 0 && value > nearMaximum && Number.isFinite(value)))]
      .sort((left, right) => left - right);
    for (const distance of distinct) {
      for (let index = 0; index < cohort.length; index += 1) {
        if (groups[index] < 0 && first[index] === distance) groups[index] = group;
      }
      group += 1;
    }
    for (let index = 0; index < groups.length; index += 1) if (groups[index] < 0) groups[index] = group++;
    const groupSizes = new Int32Array(group + 1);
    for (const value of groups) groupSizes[value] += 1;
    for (let left = 0; left < cohort.length; left += 1) {
      for (let right = left + 1; right < cohort.length; right += 1) {
        if (!comparisonChanged(first[left], first[right], second[left], second[right])) continue;
        scores[candidateIndex] += 1 / Math.max(1, groupSizes[groups[left]] * groupSizes[groups[right]]);
      }
    }
  }
  return scores;
}

function distanceCategory(value) {
  return Math.round(value * 1_000_000);
}

// Direct MakeLDist + MakeRCompat control flow for one co-recombinant set and
// one tree-distance matrix. Tree-path distances are category-equivalent to
// the desktop's /1000 representation, so category identity rather than its
// arbitrary numeric scale is retained.
function sourceRCompat(matrix, size, recombinantMembers, parentPositions) {
  const members = [...new Set(recombinantMembers)];
  if (members.length === 0) return 0;
  const memberSet = new Set(members);
  let limitingDistance = 0;
  for (let left = 0; left < members.length - 1; left += 1) {
    for (let right = left + 1; right < members.length; right += 1) {
      limitingDistance = Math.max(limitingDistance, matrixDistance(matrix, size, members[left], members[right]));
    }
  }
  const done = new Set(parentPositions);
  const nonRecombinants = [];
  for (const member of members) {
    for (let sequence = 0; sequence < size; sequence += 1) {
      if (done.has(sequence)) continue;
      if (matrixDistance(matrix, size, member, sequence) >= limitingDistance) continue;
      if (memberSet.has(sequence)) continue;
      done.add(sequence);
      nonRecombinants.push(sequence);
    }
  }
  let recombinantCost = 0;
  for (const member of members) {
    const categories = new Set(nonRecombinants.map((sequence) => distanceCategory(
      matrixDistance(matrix, size, member, sequence),
    )));
    for (const parent of parentPositions) {
      if (matrixDistance(matrix, size, member, parent) >= limitingDistance) continue;
      for (const otherMember of members) {
        categories.add(distanceCategory(matrixDistance(matrix, size, otherMember, parent)));
      }
    }
    recombinantCost = Math.max(recombinantCost, categories.size);
  }
  const expandedNonRecombinants = [...nonRecombinants];
  for (const parent of parentPositions) {
    if (members.some((member) => matrixDistance(matrix, size, member, parent) < limitingDistance)) {
      expandedNonRecombinants.push(parent);
    }
  }
  if (expandedNonRecombinants.length > 0) {
    let outsideCost = 0;
    for (const sequence of expandedNonRecombinants) {
      const categories = new Set(members.map((member) => distanceCategory(
        matrixDistance(matrix, size, sequence, member),
      )));
      outsideCost = Math.max(outsideCost, categories.size - 1);
    }
    recombinantCost = Math.min(recombinantCost, outsideCost);
  }
  return Math.min(recombinantCost, Math.max(0, members.length - 1));
}

function segmentsOverlap(left, right, length) {
  for (const [leftStart, leftEnd] of tractSegments(left, length)) {
    for (const [rightStart, rightEnd] of tractSegments(right, length)) {
      if (Math.min(leftEnd, rightEnd) > Math.max(leftStart, rightStart)) return true;
    }
  }
  return false;
}

// Direct FindSets/DoSetsAP closure. Each already reconstructed event acts as a
// three-member hyperedge. A sequence present in exactly two candidate sets is
// added to the third, and the event graph is traversed again until stable.
export function sourceHistoricalSetMembers(candidates, cohort, eventCorpus = [], focalEvent = null, length = 0) {
  const cohortSet = new Set(cohort);
  const focalTriad = focalEvent
    ? [focalEvent.recombinant, focalEvent.majorParent, focalEvent.minorParent]
    : [];
  const edges = eventCorpus.flatMap((event) => {
    if (!event || event.decision === "rejected") return [];
    const edge = [event.recombinant, event.majorParent, event.minorParent];
    if (new Set(edge).size !== 3 || edge.some((taxon) => !cohortSet.has(taxon))) return [];
    const isFocal = event === focalEvent
      || (focalEvent?.id && event.id === focalEvent.id)
      || (focalEvent && edge.every((taxon) => focalTriad.includes(taxon))
        && event.start === focalEvent.start && event.end === focalEvent.end);
    if (isFocal || (focalEvent && length > 0 && !segmentsOverlap(event, focalEvent, length))) return [];
    return [edge];
  });
  const groups = candidates.map((candidate) => new Set([candidate]));
  let changed = true;
  let rounds = 0;
  while (changed && rounds++ <= cohort.length + edges.length + 3) {
    changed = false;
    for (const group of groups) {
      for (const edge of edges) {
        if (!edge.some((taxon) => group.has(taxon))) continue;
        for (const taxon of edge) {
          if (!group.has(taxon)) {
            group.add(taxon);
            changed = true;
          }
        }
      }
    }
    for (const taxon of cohort) {
      const membership = groups.map((group) => group.has(taxon));
      if (membership.filter(Boolean).length !== 2) continue;
      const missing = membership.indexOf(false);
      groups[missing].add(taxon);
      changed = true;
    }
  }
  return groups.map((group) => [...group].sort((left, right) => left - right));
}

export function sourceParsimonyScores(candidates, cohort, treePair, coRecombinantSets = [], eventCorpus = [], focalEvent = null, length = 0) {
  const positions = new Map(cohort.map((taxon, index) => [taxon, index]));
  const ordinaryMembers = candidates.map((candidate) => {
    const set = coRecombinantSets.find((entry) => entry.presumedRecombinant === candidate);
    return [...new Set([candidate, ...(set?.sequenceMembers ?? [])])];
  });
  const historicalMembers = sourceHistoricalSetMembers(candidates, cohort, eventCorpus, focalEvent, length);
  const scoreMatrix = (matrix, memberSets) => candidates.map((candidate, candidateIndex) => {
    const members = memberSets[candidateIndex]
      .map((taxon) => positions.get(taxon))
      .filter((position) => position !== undefined);
    const parents = candidates
      .filter((_, index) => index !== candidateIndex)
      .map((taxon) => positions.get(taxon))
      .filter((position) => position !== undefined);
    return sourceRCompat(matrix, cohort.length, members, parents);
  });
  const informative = (values) => values.some((value) => value !== values[0]);
  const cascade = (tree) => [
    { stage: "ordinary uncollapsed", values: scoreMatrix(tree.uncollapsed, ordinaryMembers) },
    { stage: "ordinary bootstrap-collapsed", values: scoreMatrix(tree.collapsed, ordinaryMembers) },
    { stage: "historical-set uncollapsed", values: scoreMatrix(tree.uncollapsed, historicalMembers) },
    { stage: "historical-set bootstrap-collapsed", values: scoreMatrix(tree.collapsed, historicalMembers) },
  ];
  const outerCascade = cascade(treePair[0]);
  const innerCascade = cascade(treePair[1]);
  const choose = (stages) => stages.find((stage) => informative(stage.values)) ?? stages.at(-1);
  const outerSelected = choose(outerCascade);
  const innerSelected = choose(innerCascade);
  return {
    outer: outerSelected.values,
    inner: innerSelected.values,
    outerStage: outerSelected.stage,
    innerStage: innerSelected.stage,
    outerCascade,
    innerCascade,
    historicalMembers,
  };
}

function closestTripletPair(candidates, cohort, profiles) {
  const positions = new Map(cohort.map((taxon, index) => [taxon, index]));
  const pairs = [[0, 1], [0, 2], [1, 2]];
  let bestPair = 0;
  let bestDistance = Infinity;
  for (let pairIndex = 0; pairIndex < pairs.length; pairIndex += 1) {
    const [left, right] = pairs[pairIndex];
    const rightPosition = positions.get(candidates[right]);
    const distance = rightPosition === undefined ? Infinity : profiles.get(candidates[left])?.[rightPosition];
    if (distance < bestDistance) {
      bestDistance = distance;
      bestPair = pairIndex;
    }
  }
  return bestPair;
}

function sourceInList(firstPair, secondPair) {
  if (firstPair === secondPair) return null;
  return {
    "0:1": [1, 0, 2],
    "0:2": [0, 1, 2],
    "1:0": [2, 0, 1],
    "1:2": [0, 2, 1],
    "2:0": [2, 1, 0],
    "2:1": [1, 2, 0],
  }[`${firstPair}:${secondPair}`] ?? null;
}

function matrixDistance(matrix, size, left, right) {
  return matrix[left * size + right];
}

// Direct translation of Module2.bas SimpleDist. The source names the binary
// result SimScore and the continuous result SimScoreB; the manual calls them
// O:E and O:EDist. The scale adjustment and all three NO/PI/NI formulas are
// intentionally retained, including the source's strict comparisons.
export function sourceSimpleDist(candidates, cohort, firstMatrix, secondMatrix, excludedTaxa = new Set()) {
  const size = cohort.length;
  const positions = new Map(cohort.map((taxon, index) => [taxon, index]));
  const candidatePositions = candidates.map((candidate) => positions.get(candidate));
  if (candidatePositions.some((position) => position === undefined)) {
    return { simScore: [0, 0, 0], simScoreB: [0, 0, 0], rankFirst: [0, 0, 0], rankSecond: [0, 0, 0], adjustment: 1, inList: null };
  }
  let firstTotal = 0;
  let secondTotal = 0;
  for (let left = 0; left < size - 1; left += 1) {
    if (excludedTaxa.has(cohort[left])) continue;
    for (let right = left + 1; right < size; right += 1) {
      if (excludedTaxa.has(cohort[right])) continue;
      const first = matrixDistance(firstMatrix, size, left, right);
      if (!(first < 3)) continue;
      firstTotal += first;
      secondTotal += matrixDistance(secondMatrix, size, left, right);
    }
  }
  const adjustment = firstTotal > 0 && secondTotal > 0 ? firstTotal / secondTotal : 1;
  const moveFirst = new Float64Array(size);
  const moveSecond = new Float64Array(size);
  for (let left = 0; left < size; left += 1) {
    for (let right = 0; right < size; right += 1) {
      moveFirst[left] += matrixDistance(firstMatrix, size, left, right);
      moveSecond[left] += matrixDistance(secondMatrix, size, left, right);
    }
  }
  const rankFirst = candidatePositions.map((position) => {
    let rank = 0;
    for (const value of moveFirst) if (moveFirst[position] > value) rank += 1;
    return rank;
  });
  const rankSecond = candidatePositions.map((position) => {
    let rank = 0;
    for (const value of moveSecond) if (moveSecond[position] > value) rank += 1;
    return rank;
  });
  const pairDistance = (matrix, leftCandidate, rightCandidate) => matrixDistance(
    matrix,
    size,
    candidatePositions[leftCandidate],
    candidatePositions[rightCandidate],
  );
  const firstProfiles = new Map(candidates.map((candidate, index) => [candidate, Float64Array.from(
    { length: 3 },
    (_, right) => pairDistance(firstMatrix, index, right),
  )]));
  const secondProfiles = new Map(candidates.map((candidate, index) => [candidate, Float64Array.from(
    { length: 3 },
    (_, right) => pairDistance(secondMatrix, index, right),
  )]));
  const tripletCohort = [...candidates];
  const firstPair = closestTripletPair(candidates, tripletCohort, firstProfiles);
  const secondPair = closestTripletPair(candidates, tripletCohort, secondProfiles);
  const inList = sourceInList(firstPair, secondPair);
  const simScore = [0, 0, 0];
  const simScoreB = [0, 0, 0];
  if (!inList) return { simScore, simScoreB, rankFirst, rankSecond, adjustment, inList };
  const [no, pi, ni] = inList;
  const f = (left, right) => pairDistance(firstMatrix, left, right);
  const s = (left, right) => pairDistance(secondMatrix, left, right) * adjustment;
  const score = (target, diff0, diff1, diff2) => {
    if (diff0 < diff1 && diff0 < diff2) simScore[target] += 1;
    simScoreB[target] = diff1 + diff2 - diff0;
  };
  score(pi, Math.abs(f(no, ni) - s(no, ni)), s(no, pi) - f(no, pi), f(pi, ni) - s(pi, ni));
  score(no, Math.abs(f(pi, ni) - s(pi, ni)), s(no, pi) - f(no, pi), s(no, ni) - f(no, ni));
  score(ni, Math.abs(f(no, pi) - s(no, pi)), f(pi, ni) - s(pi, ni), f(no, ni) - s(no, ni));
  return { simScore, simScoreB, rankFirst, rankSecond, adjustment, inList };
}

// Direct MakeSSDistB data reduction. RDP5 first rescales the two raw distance
// matrices, then averages squared movement within each background-tree path
// category so large clades do not dominate merely by containing more taxa.
export function sourceSsDistScores(candidates, cohort, firstMatrix, secondMatrix, backgroundTreeMatrix, excludedTaxa = new Set()) {
  const size = cohort.length;
  const positions = new Map(cohort.map((taxon, index) => [taxon, index]));
  let firstTotal = 0;
  let secondTotal = 0;
  for (let left = 0; left < size - 1; left += 1) {
    if (excludedTaxa.has(cohort[left])) continue;
    for (let right = left + 1; right < size; right += 1) {
      if (excludedTaxa.has(cohort[right])) continue;
      firstTotal += matrixDistance(firstMatrix, size, left, right);
      secondTotal += matrixDistance(secondMatrix, size, left, right);
    }
  }
  if (!(secondTotal > 0)) return { values: [0, 0, 0], adjustment: 0 };
  const adjustment = firstTotal / secondTotal;
  const values = candidates.map((candidate) => {
    const candidatePosition = positions.get(candidate);
    if (candidatePosition === undefined) return 0;
    const categorySums = new Map();
    const categoryCounts = new Map();
    for (let taxonPosition = 0; taxonPosition < size; taxonPosition += 1) {
      if (excludedTaxa.has(cohort[taxonPosition])) continue;
      const category = distanceCategory(matrixDistance(backgroundTreeMatrix, size, candidatePosition, taxonPosition));
      const movement = Math.abs(
        matrixDistance(firstMatrix, size, candidatePosition, taxonPosition)
        - matrixDistance(secondMatrix, size, candidatePosition, taxonPosition),
      ) * adjustment;
      categorySums.set(category, (categorySums.get(category) ?? 0) + movement * movement);
      categoryCounts.set(category, (categoryCounts.get(category) ?? 0) + 1);
    }
    let score = 0;
    for (const [category, sum] of categorySums) score += sum / categoryCounts.get(category);
    return score;
  });
  return { values, adjustment };
}

export function sourceOuIndexScores(ssDist, inList) {
  const scores = [0, 0, 0];
  if (!inList) return scores;
  const [no, pi, ni] = inList;
  if (ssDist[no] > ssDist[pi] && ssDist[no] > ssDist[ni]) scores[no] = 1;
  else if (ssDist[no] < ssDist[pi] && ssDist[no] < ssDist[ni]) {
    scores[pi] = 1;
    scores[ni] = 1;
  }
  return scores;
}

function tripletVariable(encoded, length, candidates, site) {
  const values = candidates.map((candidate) => encoded[candidate * length + site]);
  return values.every((value) => value < 4) && (values[0] !== values[1] || values[0] !== values[2]);
}

function sourceFlankSites(encoded, length, candidates, breakpoint, direction, wanted, circular) {
  const sites = [];
  const maximumSteps = circular ? length : direction < 0 ? breakpoint : length - breakpoint;
  for (let step = 0; step < maximumSteps && sites.length < wanted; step += 1) {
    let site = direction < 0 ? breakpoint - step - 1 : breakpoint + step;
    if (circular) site = (site % length + length) % length;
    if (site < 0 || site >= length) break;
    if (tripletVariable(encoded, length, candidates, site)) sites.push(site);
  }
  return sites;
}

function rawSiteDistance(encoded, length, first, second, sites) {
  let valid = 0;
  let differences = 0;
  for (const site of sites) {
    const left = encoded[first * length + site];
    const right = encoded[second * length + site];
    if (left >= 4 || right >= 4) continue;
    valid += 1;
    if (left !== right) differences += 1;
  }
  return valid > 0 ? differences / valid : 1;
}

// Direct GetBadDists decision logic. The upstream source-style co-recombinant
// screen supplies RCorr/RInv-equivalent per-breakpoint evidence; this routine
// applies the original 0.83 gate, asymmetric comparisons, and counts distinct
// background-tree distance categories rather than raw offending sequences.
export function sourceConflictScores(event, encoded, length, candidates, cohort, backgroundTreeMatrix, options = {}) {
  const positions = new Map(cohort.map((taxon, index) => [taxon, index]));
  const wanted = Math.max(4, Math.min(200, Math.trunc(options.clusterFlankVnps ?? 60)));
  const circular = options.circular === true || event.wraps === true;
  const regions = [
    sourceFlankSites(encoded, length, candidates, event.start, -1, wanted, circular),
    sourceFlankSites(encoded, length, candidates, event.start, 1, wanted, circular),
    sourceFlankSites(encoded, length, candidates, event.end % length, -1, wanted, circular),
    sourceFlankSites(encoded, length, candidates, event.end % length, 1, wanted, circular),
  ];
  return candidates.map((candidate, candidateIndex) => {
    const candidatePosition = positions.get(candidate);
    const set = event.coRecombinantSets?.find((entry) => entry.presumedRecombinant === candidate);
    if (candidatePosition === undefined || !set) return 0;
    const parents = candidates.filter((_, index) => index !== candidateIndex);
    const badCategories = new Set();
    for (const evidence of set.evidence ?? []) {
      const sequencePosition = positions.get(evidence.sequence);
      if (sequencePosition === undefined || evidence.sequence === candidate) continue;
      const regional = evidence.regionEvidence ?? [];
      const startR = regional[0]?.correlationSdmFiltered ? 0 : regional[0]?.correlationR ?? 0;
      const endR = regional[1]?.correlationSdmFiltered ? 0 : regional[1]?.correlationR ?? 0;
      const distance = (regionIndex, sequence) => rawSiteDistance(encoded, length, candidate, sequence, regions[regionIndex]);
      let bad = false;
      if (startR > 0.83) {
        const candidateLeft = distance(0, evidence.sequence);
        const candidateRight = distance(1, evidence.sequence);
        bad = parents.some((parent) => candidateLeft > distance(0, parent) || candidateRight >= distance(1, parent));
      }
      if (!bad && endR > 0.83) {
        const candidateLeft = distance(2, evidence.sequence);
        const candidateRight = distance(3, evidence.sequence);
        bad = parents.some((parent) => candidateLeft >= distance(2, parent) || candidateRight > distance(3, parent));
      }
      if (bad) badCategories.add(distanceCategory(
        matrixDistance(backgroundTreeMatrix, cohort.length, candidatePosition, sequencePosition),
      ));
    }
    return badCategories.size;
  });
}

// Direct MakeEList + MakeListCorr reduction for the two set-distance tests.
// SetDistT counts disagreement between observed RCorr inversion classes and
// the classes expected under each NO/PI/NI polarity. SetDistP averages the
// corresponding CalCR permutation correlations with the desktop's exact
// denominator and all-or-zero coverage rule.
export function sourceSetDistanceScores(event, candidates, cohort, firstMatrix, secondMatrix, inList) {
  const size = cohort.length;
  if (!inList) return { setDistT: [0, 0, 0], setDistP: [0, 0, 0], expectedCoverage: [0, 0, 0] };
  const positions = new Map(cohort.map((taxon, index) => [taxon, index]));
  const candidatePositions = candidates.map((candidate) => positions.get(candidate));
  if (candidatePositions.some((position) => position === undefined)) {
    return { setDistT: [0, 0, 0], setDistP: [0, 0, 0], expectedCoverage: [0, 0, 0] };
  }
  const expected = Array.from({ length: 3 }, () => Array.from(
    { length: 3 },
    () => new Int8Array(size).fill(-1),
  ));
  const [no, pi, ni] = inList;
  const value = (matrix, candidateIndex, taxonPosition) => matrixDistance(
    matrix,
    size,
    candidatePositions[candidateIndex],
    taxonPosition,
  );
  const pair = (matrix, leftCandidate, rightCandidate) => value(matrix, leftCandidate, candidatePositions[rightCandidate]);
  for (let taxon = 0; taxon < size; taxon += 1) {
    // Hypothesis 0: NO is recombinant.
    if (value(firstMatrix, no, taxon) < pair(firstMatrix, no, pi)
      && value(secondMatrix, no, taxon) < pair(secondMatrix, no, pi)) expected[0][no][taxon] = 0;
    if (value(firstMatrix, no, taxon) < pair(firstMatrix, no, ni)
      && value(firstMatrix, no, taxon) > 0
      && value(secondMatrix, pi, taxon) < pair(secondMatrix, no, pi)) expected[0][pi][taxon] = 0;
    if (value(firstMatrix, no, taxon) > pair(firstMatrix, no, ni)
      && value(secondMatrix, no, taxon) > value(secondMatrix, ni, taxon)) expected[0][pi][taxon] = 2;
    if (value(firstMatrix, no, taxon) > pair(firstMatrix, no, pi)
      && value(secondMatrix, no, taxon) > value(secondMatrix, pi, taxon)) expected[0][ni][taxon] = 0;
    if (value(firstMatrix, no, taxon) < pair(firstMatrix, no, pi)
      && value(firstMatrix, no, taxon) > 0
      && value(secondMatrix, no, taxon) > value(secondMatrix, pi, taxon)) expected[0][ni][taxon] = 2;

    // Hypothesis 1: PI is recombinant.
    if (value(firstMatrix, pi, taxon) > 0
      && value(firstMatrix, pi, taxon) < pair(firstMatrix, pi, ni)
      && value(secondMatrix, pi, taxon) > value(secondMatrix, no, taxon)) expected[1][no][taxon] = 0;
    if (value(firstMatrix, ni, taxon) > 0
      && value(firstMatrix, ni, taxon) < pair(firstMatrix, pi, ni)
      && value(secondMatrix, ni, taxon) > pair(secondMatrix, pi, ni)) expected[1][no][taxon] = 1;
    if (value(firstMatrix, pi, taxon) < pair(firstMatrix, pi, no)
      && value(secondMatrix, pi, taxon) < pair(secondMatrix, pi, no)) expected[1][pi][taxon] = 0;
    if (value(firstMatrix, ni, taxon) < pair(firstMatrix, no, ni)
      && value(secondMatrix, ni, taxon) < pair(firstMatrix, no, ni)) expected[1][ni][taxon] = 0;
    if (value(firstMatrix, pi, taxon) > pair(firstMatrix, pi, no)
      && value(firstMatrix, pi, taxon) < pair(firstMatrix, pi, ni)) expected[1][ni][taxon] = 1;
    if (value(firstMatrix, pi, taxon) > 0
      && value(firstMatrix, pi, taxon) < pair(firstMatrix, pi, no)
      && value(secondMatrix, pi, taxon) < pair(secondMatrix, pi, no)) expected[1][ni][taxon] = 4;

    // Hypothesis 2: NI is recombinant.
    if (value(firstMatrix, pi, taxon) > 0
      && value(firstMatrix, pi, taxon) < pair(firstMatrix, pi, ni)) expected[2][no][taxon] = 0;
    if (value(firstMatrix, ni, taxon) > 0
      && value(firstMatrix, ni, taxon) < pair(firstMatrix, no, ni)
      && value(secondMatrix, ni, taxon) < value(secondMatrix, no, taxon)) expected[2][no][taxon] = 1;
    if (value(firstMatrix, pi, taxon) < pair(firstMatrix, pi, no)
      && value(secondMatrix, pi, taxon) < pair(secondMatrix, pi, no)) expected[2][pi][taxon] = 0;
    if (value(firstMatrix, ni, taxon) > 0
      && value(firstMatrix, ni, taxon) < pair(firstMatrix, no, ni)
      && value(secondMatrix, ni, taxon) < pair(firstMatrix, no, ni)) expected[2][pi][taxon] = 4;
    if (value(firstMatrix, pi, taxon) > pair(firstMatrix, pi, no)
      && value(firstMatrix, pi, taxon) < pair(firstMatrix, pi, ni)) expected[2][pi][taxon] = 2;
    if (value(firstMatrix, ni, taxon) < pair(firstMatrix, no, ni)
      && value(secondMatrix, ni, taxon) < pair(secondMatrix, no, ni)) expected[2][ni][taxon] = 0;
  }

  const actual = Array.from({ length: 3 }, () => new Int8Array(size).fill(-1));
  const correlations = Array.from({ length: 3 }, () => new Map());
  for (let candidateIndex = 0; candidateIndex < 3; candidateIndex += 1) {
    const candidate = candidates[candidateIndex];
    const set = event.coRecombinantSets?.find((entry) => entry.presumedRecombinant === candidate);
    const candidatePosition = candidatePositions[candidateIndex];
    actual[candidateIndex][candidatePosition] = 0;
    for (const evidence of set?.evidence ?? []) {
      const taxonPosition = positions.get(evidence.sequence);
      if (taxonPosition === undefined) continue;
      const regional = evidence.regionEvidence ?? [];
      const inversions = regional.map((entry) => entry.correlationInversion);
      actual[candidateIndex][taxonPosition] = inversions.includes(1) ? 1
        : inversions.includes(2) ? 2
          : inversions.includes(4) ? 4 : 0;
      correlations[candidateIndex].set(taxonPosition, regional.slice(0, 3).map((entry) => {
        const permutation = entry.correlationPermutations ?? [];
        const categories = [
          permutation[0] ?? 0,
          permutation[1] ?? 0,
          permutation[2] ?? 0,
          permutation[3] ?? 0,
          Math.max(permutation[4] ?? 0, permutation[5] ?? 0),
        ].map((score) => score > 0.5 && score < 1 ? score : 0);
        const shared = Math.max(categories[2], categories[3]);
        categories[2] = shared;
        categories[3] = shared;
        return categories;
      }));
    }
  }

  const setDistT = [0, 0, 0];
  const setDistP = [0, 0, 0];
  const missingScore = [0, 0, 0];
  const expectedCoverage = [0, 0, 0];
  const missingCoverage = [0, 0, 0];
  for (let orientation = 0; orientation < 3; orientation += 1) {
    for (let taxon = 0; taxon < size; taxon += 1) {
      if (actual[orientation][taxon] >= 0) {
        for (let hypothesis = 0; hypothesis < 3; hypothesis += 1) {
          if (actual[orientation][taxon] !== expected[hypothesis][orientation][taxon]) setDistT[inList[hypothesis]] += 1;
        }
      }
      for (let hypothesis = 0; hypothesis < 3; hypothesis += 1) {
        const candidate = inList[hypothesis];
        const expectedCategory = expected[hypothesis][orientation][taxon];
        const regional = correlations[orientation].get(taxon) ?? [];
        if (expectedCategory >= 0) {
          for (let region = 0; region < 3; region += 1) {
            const score = regional[region]?.[expectedCategory] ?? 0;
            if (score > 0) setDistP[candidate] += score;
          }
          expectedCoverage[candidate] += 3;
        } else {
          for (let region = 0; region < 3; region += 1) {
            const score = regional[region]?.[0] ?? 0;
            if (score > 0) missingScore[candidate] += score;
          }
          missingCoverage[candidate] += 3;
        }
      }
    }
  }
  if (expectedCoverage.some((count) => count === 0) || missingCoverage.some((count) => count === 0)) {
    setDistP.fill(0);
  } else {
    for (let index = 0; index < 3; index += 1) setDistP[index] /= expectedCoverage[index];
  }
  return { setDistT, setDistP, expectedCoverage, missingScore, missingCoverage };
}

// Direct MakeINList + MakeOUCheck interpretation. INList describes the
// topology change as NO (outlier in the recombinant region), PI (inlier in
// both regions), and NI (outlier in the background). MakeOUCheck then counts
// which triplet member most often disturbs its relationship with every other
// taxon between the two collapsed trees.
export function sourceOuCheckScores(candidates, cohort, profiles) {
  const firstPair = closestTripletPair(candidates, cohort, profiles.first);
  const secondPair = closestTripletPair(candidates, cohort, profiles.second);
  const inList = sourceInList(firstPair, secondPair);
  if (!inList) return [0, 0, 0];
  const no = inList[0];
  const pi = inList[1];
  const ni = inList[2];
  const positions = new Map(cohort.map((taxon, index) => [taxon, index]));
  const noFirst = profiles.first.get(candidates[no]);
  const piSecond = profiles.second.get(candidates[pi]);
  const noSecond = profiles.second.get(candidates[no]);
  const piPosition = positions.get(candidates[pi]);
  const niPosition = positions.get(candidates[ni]);
  if (piPosition === undefined || niPosition === undefined) return [0, 0, 0];
  const scores = [0, 0, 0];
  for (let index = 0; index < cohort.length; index += 1) {
    if (!(noFirst[index] > noFirst[piPosition] && noFirst[index] < noFirst[niPosition])) continue;
    if (piSecond[index] < noSecond[index]) {
      scores[no] += 1;
      scores[pi] -= 1;
      scores[ni] -= 1;
    } else if (piSecond[index] > noSecond[index]) {
      scores[no] -= 1;
      scores[pi] += 1;
      scores[ni] -= 1;
    } else if (piSecond[index] > noSecond[piPosition]) {
      scores[no] -= 1;
      scores[pi] -= 1;
      scores[ni] += 1;
    }
  }
  return scores;
}

function nearlyEqual(left, right) {
  return Math.abs(left - right) <= 1e-10 * Math.max(1, Math.abs(left), Math.abs(right));
}

function weightedTest(id, label, routine, values, direction, fullWeight, partialWeight) {
  const finite = values.filter(Number.isFinite);
  const best = direction === "lower" ? Math.min(...finite) : Math.max(...finite);
  const worst = direction === "lower" ? Math.max(...finite) : Math.min(...finite);
  const points = values.map((value) => {
    if (!Number.isFinite(value)) return 0;
    if (nearlyEqual(value, best)) return fullWeight;
    const beatsAnother = values.some((other) => Number.isFinite(other) && (
      direction === "lower" ? value < other : value > other
    ));
    return beatsAnother ? partialWeight : 0;
  });
  const winnerIndexes = values.flatMap((value, index) => nearlyEqual(value, best) ? [index] : []);
  return {
    id,
    label,
    sourceRoutine: routine,
    direction,
    values: values.map((value) => Number.isFinite(value) ? value : null),
    points,
    fullWeight,
    partialWeight,
    winnerIndexes,
    decisive: !nearlyEqual(best, worst),
  };
}

function signedOuCheckTest(values) {
  const best = Math.max(...values);
  const worst = Math.min(...values);
  return {
    id: "oucheck",
    label: "OuCheck",
    sourceRoutine: "MakeINList + MakeOUCheck",
    direction: "higher",
    values,
    points: values.map((value) => value > 0 ? 5 : value < 0 ? -5 : 0),
    fullWeight: 5,
    partialWeight: 0,
    winnerIndexes: values.flatMap((value, index) => nearlyEqual(value, best) ? [index] : []),
    decisive: !nearlyEqual(best, worst),
  };
}

function explicitPointTest(id, label, routine, values, points, direction = "higher") {
  const finite = values.filter(Number.isFinite);
  const best = direction === "lower" ? Math.min(...finite) : Math.max(...finite);
  const worst = direction === "lower" ? Math.max(...finite) : Math.min(...finite);
  return {
    id,
    label,
    sourceRoutine: routine,
    direction,
    values: values.map((value) => Number.isFinite(value) ? value : null),
    points,
    fullWeight: Math.max(0, ...points),
    partialWeight: 0,
    winnerIndexes: values.flatMap((value, index) => nearlyEqual(value, best) ? [index] : []),
    decisive: !nearlyEqual(best, worst) || points.some((point) => !nearlyEqual(point, points[0])),
  };
}

// Direct MakeConsensusC dMax block. Besides the normalized 20-point share,
// RDP5 adds tiered bonuses for the largest VisRD displacement. Equalities are
// intentionally inclusive because that is how the source resolves ties.
function sourceDmaxTest(values) {
  const total = values.reduce((sum, value) => sum + Math.max(0, value), 0);
  const points = values.map((value, index) => {
    const others = values.filter((_, otherIndex) => otherIndex !== index);
    let score = total > 0 ? value / total * 20 : 0;
    if (value >= others[0] * 1.1 && value >= others[1] * 1.1) score += 30;
    else if (value >= others[0] && value >= others[1]) score += 20;
    else if (value >= others[0] * 1.1 || value >= others[1] * 1.1) score += 10;
    else if (value >= others[0] || value >= others[1]) score += 5;
    return score;
  });
  return explicitPointTest(
    "dmax",
    "dMax (VisRD)",
    "CalcMaxD + CMaxD2P3 + MakeConsensusC",
    values,
    points,
  );
}

function sourceParsimonyTest(id, label, routine, values) {
  if (values.every((value) => nearlyEqual(value, values[0]))) {
    return explicitPointTest(id, label, routine, values, [0, 0, 0], "lower");
  }
  return weightedTest(id, label, routine, values, "lower", 20, 10);
}

function sourceConflictTest(values) {
  if (values.every((value) => nearlyEqual(value, values[0]))) {
    return explicitPointTest("conflict", "Conflict", "GetBadDists", values, [1, 1, 1], "lower");
  }
  return weightedTest("conflict", "Conflict", "GetBadDists", values, "lower", 10, 5);
}

function sourceSsDistTest(values) {
  const points = values.map((value, index) => {
    const others = values.filter((_, otherIndex) => otherIndex !== index);
    if (nearlyEqual(value, others[0]) || nearlyEqual(value, others[1])) return 0;
    if (value >= others[0] && value >= others[1]) return 5;
    if (value > others[0] || value > others[1]) return 2.5;
    return 0;
  });
  return explicitPointTest("ssdist", "SSDist", "MakeSSDistB", values, points);
}

function sourceOuIndexTest(values) {
  return explicitPointTest("ouindex", "OUIndex", "MakeSSDistB OUIndexA", values, values.map((value) => value === 1 ? 5 : 0));
}

function sourceSetDistPTest(values, sequenceCount) {
  const greatest = Math.max(...values);
  const points = values.map((value, index) => {
    if (!(sequenceCount > 11) || !(value > 0)) return 0;
    if (nearlyEqual(value, greatest)) return 20;
    return values.some((other, otherIndex) => otherIndex !== index && value > other) ? 15 : 0;
  });
  return explicitPointTest("setdist-p", "SetDistP", "MakeEList + MakeListCorr ListCorr2", values, points);
}

export function sourceFinalTrimParsimonyTest(outer, inner) {
  const values = outer.map((value, index) => Math.min(value, inner[index]));
  const points = values.map((value) => value > 1 ? -(value - 1) * 10 : 0);
  return explicitPointTest(
    "final-trim-parsimony",
    "FinalTrim penalty",
    "FinalTrim + RCompatC/RCompatD + MakeConsensusC",
    values,
    points,
    "lower",
  );
}

function gateProfileTest(test, enabled) {
  if (enabled) return test;
  return { ...test, points: [0, 0, 0], decisive: false, sourceRoutine: `${test.sourceRoutine} (PS gate inactive)` };
}

function sourceSynergyTest(id, label, routine, values, qualifies, weight) {
  const points = values.map((_, index) => qualifies(index) ? weight : 0);
  return explicitPointTest(id, label, routine, values, points);
}

function preferredParents(candidates, raw) {
  return candidates.map((recombinant, recombinantIndex) => {
    const parents = candidates.filter((_, index) => index !== recombinantIndex);
    const recombinantPosition = new Map(candidates.map((candidate, index) => [candidate, index]));
    const first = raw.first.get(recombinant);
    const second = raw.second.get(recombinant);
    const cohort = raw.cohort;
    const positions = new Map(cohort.map((taxon, index) => [taxon, index]));
    const scoreAssignment = (majorParent, minorParent) => {
      const majorPosition = positions.get(majorParent);
      const minorPosition = positions.get(minorParent);
      if (majorPosition === undefined || minorPosition === undefined) return -Infinity;
      return (second[majorPosition] - second[minorPosition])
        + (first[minorPosition] - first[majorPosition]);
    };
    const direct = scoreAssignment(parents[0], parents[1]);
    const reverse = scoreAssignment(parents[1], parents[0]);
    const ordered = reverse > direct ? [parents[1], parents[0]] : parents;
    return {
      recombinant,
      majorParent: ordered[0],
      minorParent: ordered[1],
      affinitySwitch: Math.max(direct, reverse),
      candidateIndex: recombinantPosition.get(recombinant),
    };
  });
}

function sourceDistanceTieBreak(candidates, cohort, raw) {
  const positions = new Map(cohort.map((taxon, index) => [taxon, index]));
  const distance = (profiles, left, right) => {
    const position = positions.get(candidates[right]);
    return position === undefined ? Infinity : profiles.get(candidates[left])?.[position] ?? Infinity;
  };
  const f01 = distance(raw.first, 0, 1);
  const f02 = distance(raw.first, 0, 2);
  const f12 = distance(raw.first, 1, 2);
  const s01 = distance(raw.second, 0, 1);
  const s02 = distance(raw.second, 0, 2);
  const s12 = distance(raw.second, 1, 2);
  // Direct GetWinPPfromDists ordering, including its inclusive ties.
  if (f01 <= f02 && f01 <= f12) return s02 <= s12 ? 0 : 1;
  if (f02 <= f01 && f02 <= f12) return s01 <= s12 ? 0 : 2;
  if (f12 <= f01 && f12 <= f02) return s01 <= s02 ? 1 : 2;
  return 2;
}

export function identifyRecombinantRoles(event, encoded, sequenceCount, length, treeBundle = null, options = {}) {
  const candidates = [event.recombinant, event.majorParent, event.minorParent];
  if (new Set(candidates).size !== 3) return null;
  const requestedMaximum = Math.max(3, Math.min(300, Math.trunc(options.roleTaxaLimit ?? options.clusterTreeTaxaLimit ?? 64)));
  const cohort = treeBundle?.taxa?.length
    ? [...treeBundle.taxa]
    : deterministicCohort(sequenceCount, candidates, requestedMaximum);
  const raw = rawProfiles(encoded, length, candidates, cohort, event);
  raw.cohort = cohort;
  const rawScores = sourceProfileConsensusInputs(sourcePhPrScores(candidates, cohort, raw), "raw");
  const rawProfileEnabled = sourceProfileGate(rawScores.phPr);
  const rawSubDistEnabled = rawProfileEnabled
    && rawScores.subDist.some((value) => value !== -1 && value !== 1);
  const tests = [
    gateProfileTest(weightedTest("phpr", "PhPr", "MakePhPrScore(FMat, SMat)", rawScores.phPr, "lower", 8, 4), rawProfileEnabled),
    gateProfileTest(weightedTest("subphpr", "SubPhPr", "MakePhPrScore SubScore", rawScores.subPhPr, "higher", 2, 1), rawProfileEnabled),
    gateProfileTest(weightedTest("subdist", "SubDist", "MakePhPrScore SubPhPrScore", rawScores.subDist, "higher", 10, 5), rawSubDistEnabled),
  ];
  let sourceTieBreakValues = null;
  const quartetLimit = Math.max(4, Math.min(34, Math.trunc(options.roleQuartetTaxaLimit ?? 30)));
  const quartetCohort = deterministicCohort(sequenceCount, candidates, quartetLimit);
  let dmaxEvidence;
  try {
    dmaxEvidence = typeof options.roleDmaxEvaluator === "function"
      ? options.roleDmaxEvaluator({ event, candidates, cohort: quartetCohort })
      : sourceDmaxScores(encoded, length, candidates, quartetCohort, event);
  } catch {
    dmaxEvidence = sourceDmaxScores(encoded, length, candidates, quartetCohort, event);
  }
  if (dmaxEvidence?.values?.length === 3) tests.push(sourceDmaxTest(dmaxEvidence.values));
  const treePair = treeBundle?.pairs?.[2];
  const tree = treeProfiles(treePair, candidates, cohort, true);
  if (tree) {
    const treeScores = sourceProfileConsensusInputs(sourcePhPrScores(candidates, cohort, tree), "tree");
    const trp = sourceTrpScores(candidates, cohort, tree).map((value) => sourceQuantize(value, 100_000));
    const ouCheck = sourceOuCheckScores(candidates, cohort, tree);
    const parsimony = sourceParsimonyScores(
      candidates,
      cohort,
      treePair,
      event.coRecombinantSets ?? [],
      options.roleEventCorpus ?? [],
      event,
      length,
    );
    sourceTieBreakValues = candidates.map((_, index) => rawScores.phPr[index]
      - trp[index]
      + parsimony.outer[index]
      + parsimony.inner[index]);
    const treeProfileEnabled = sourceProfileGate(treeScores.phPr);
    const otherIndexes = (index) => [0, 1, 2].filter((other) => other !== index);
    const greatest = (values, index) => otherIndexes(index).every((other) => values[index] >= values[other]);
    const lowest = (values, index) => otherIndexes(index).every((other) => values[index] <= values[other]);
    const differs = (values) => values.some((value) => !nearlyEqual(value, values[0]));
    tests.push(
      gateProfileTest(weightedTest("tree-phpr", "TreePhPr", "MakePhPrScore(FAMat, SAMat)", treeScores.phPr, "lower", 18, 14), treeProfileEnabled),
      gateProfileTest(weightedTest("tree-subphpr", "TreeSubPhPr", "MakePhPrScore tree SubScore", treeScores.subPhPr, "higher", 10, 5), treeProfileEnabled),
      gateProfileTest(weightedTest("tree-subdist", "TreeSubDist", "MakePhPrScore tree SubPhPrScore", treeScores.subDist, "higher", 8, 4), treeProfileEnabled),
      sequenceCount > 11
        ? weightedTest("trpscore", "TrpScore", "MakeTrpGroups + MakeTrpScore", trp, "higher", 8, 4)
        : { ...weightedTest("trpscore", "TrpScore", "MakeTrpGroups + MakeTrpScore", trp, "higher", 8, 4), points: [0, 0, 0], decisive: false, sourceActive: false },
      signedOuCheckTest(ouCheck),
      sourceParsimonyTest("parsimony-o", "ParsimonyO", `MakeLDist + MakeRCompat (${parsimony.outerStage})`, parsimony.outer),
      sourceParsimonyTest("parsimony-i", "ParsimonyI", `MakeLDist + MakeRCompat (${parsimony.innerStage})`, parsimony.inner),
      sourceFinalTrimParsimonyTest(parsimony.outer, parsimony.inner),
      sourceSynergyTest(
        "phpr-parsimony-o",
        "PhPr × ParsimonyO",
        "MakeConsensusC MaxS joint rule",
        candidates.map((_, index) => lowest(rawScores.phPr, index) && lowest(parsimony.outer, index) ? 1 : 0),
        (index) => rawProfileEnabled
          && lowest(rawScores.phPr, index)
          && differs(parsimony.outer)
          && lowest(parsimony.outer, index),
        10,
      ),
      sourceSynergyTest(
        "parsimony-i-trpscore",
        "ParsimonyI × TrpScore",
        "MakeConsensusC MaxS joint rule",
        candidates.map((_, index) => lowest(parsimony.inner, index) && greatest(trp, index) ? 1 : 0),
        (index) => differs(parsimony.inner) && lowest(parsimony.inner, index) && greatest(trp, index),
        10,
      ),
    );
    const firstDistances = treePair?.[0]?.baseDistances;
    const secondDistances = treePair?.[1]?.baseDistances;
    if (firstDistances?.length === cohort.length ** 2 && secondDistances?.length === cohort.length ** 2) {
      const excludedTaxa = new Set((event.coRecombinantSets ?? []).flatMap((set) => set.sequenceMembers ?? []));
      const simple = sourceSimpleDist(candidates, cohort, firstDistances, secondDistances, excludedTaxa);
      const ssDistRaw = sourceSsDistScores(candidates, cohort, firstDistances, secondDistances, treePair[0].uncollapsed);
      const ssDist = {
        ...ssDistRaw,
        values: ssDistRaw.values.map((value) => value < 10_000 ? sourceQuantize(value, 100_000) : value),
      };
      const ouIndex = sourceOuIndexScores(ssDist.values, simple.inList);
      const conflicts = sourceConflictScores(event, encoded, length, candidates, cohort, treePair[0].uncollapsed, options);
      const setDistancesRaw = sourceSetDistanceScores(event, candidates, cohort, firstDistances, secondDistances, simple.inList);
      const setDistances = {
        ...setDistancesRaw,
        setDistP: setDistancesRaw.setDistP.map((value) => value > 0 ? sourceQuantize(value, 100_000) : value),
      };
      const simpleScoreB = simple.simScoreB.map((value) => sourceQuantize(value, 100_000));
      tests.push(
        explicitPointTest("oe", "O:E", "SimpleDist SimScore", simple.simScore, [0, 0, 0]),
        explicitPointTest("oedist", "O:EDist", "SimpleDist SimScoreB", simpleScoreB, [0, 0, 0]),
        sourceSsDistTest(ssDist.values),
        sourceOuIndexTest(ouIndex),
        sourceConflictTest(conflicts),
        weightedTest("setdist-t", "SetDistT", "MakeEList + MakeListCorr ListCorr", setDistances.setDistT, "lower", 4, 2),
        sourceSetDistPTest(setDistances.setDistP, sequenceCount),
      );
      const uncollapsed = treeProfiles(treePair, candidates, cohort, false);
      const uncollapsedScores = uncollapsed
        ? sourceProfileConsensusInputs(sourcePhPrScores(candidates, cohort, uncollapsed), "tree")
        : null;
      if (uncollapsedScores) {
        const uncollapsedProfileEnabled = sourceProfileGate(uncollapsedScores.phPr);
        tests.push(sourceSynergyTest(
          "subdist-treephpr",
          "SubDist × uncollapsed TreePhPr",
          "MakeConsensusC SMaxS joint rule",
          candidates.map((_, index) => greatest(rawScores.subDist, index) && lowest(uncollapsedScores.phPr, index) ? 1 : 0),
          (index) => rawProfileEnabled && uncollapsedProfileEnabled
            && greatest(rawScores.subDist, index) && lowest(uncollapsedScores.phPr, index),
          20,
        ));
      }
      tests.push(
        sourceSynergyTest(
          "oe-oucheck",
          "O:E × OuCheck",
          "MakeConsensusC MaxS joint rule",
          simple.simScore,
          (index) => {
            const others = otherIndexes(index);
            const topologyDiffers = ouCheck[index] !== ouCheck[others[0]] || ouCheck[index] !== ouCheck[others[1]];
            return topologyDiffers && greatest(ouCheck, index) && simple.simScore[index] === 1;
          },
          10,
        ),
        sourceSynergyTest(
          "oedist-tree-subdist",
          "O:EDist × TreeSubDist",
          "MakeConsensusC MidS joint rule",
          simpleScoreB,
          (index) => treeProfileEnabled && greatest(treeScores.subDist, index) && greatest(simpleScoreB, index),
          5,
        ),
        sourceSynergyTest(
          "tree-subphpr-conflict",
          "TreeSubPhPr × Conflict",
          "MakeConsensusC MidS joint rule",
          treeScores.subPhPr,
          (index) => treeProfileEnabled
            && greatest(treeScores.subPhPr, index)
            && differs(conflicts)
            && lowest(conflicts, index),
          5,
        ),
      );
    }
  }
  const unshiftedScores = candidates.map((_, index) => tests.reduce((total, test) => total + test.points[index], 0));
  const minimumScore = Math.min(0, ...unshiftedScores);
  const scores = unshiftedScores.map((score) => score - minimumScore);
  const maximum = Math.max(...scores);
  const total = scores.reduce((sum, value) => sum + value, 0);
  const maximumIndexes = scores.flatMap((score, index) => nearlyEqual(score, maximum) ? [index] : []);
  // MakeConsensusC's caller walks 0..2 with inclusive >= tests, so the last
  // maximum wins. If all three are tied, RDP5 applies its TBreak statistic or
  // the exact GetWinPPfromDists fallback rather than favoring the input role.
  let recommendedIndex = maximumIndexes.at(-1) ?? 2;
  let sourceTieBreak = "last-inclusive-maximum";
  if (maximumIndexes.length === 3) {
    const profileInformative = rawScores.phPr.some((value) => value !== -1 && value !== 1);
    if (profileInformative && sourceTieBreakValues) {
      if (sourceTieBreakValues[0] < sourceTieBreakValues[1] && sourceTieBreakValues[0] < sourceTieBreakValues[2]) recommendedIndex = 0;
      else if (sourceTieBreakValues[1] < sourceTieBreakValues[0] && sourceTieBreakValues[1] < sourceTieBreakValues[2]) recommendedIndex = 1;
      else recommendedIndex = 2;
      sourceTieBreak = "PhPr-TrpScore+ParsimonyO+ParsimonyI";
    } else {
      recommendedIndex = sourceDistanceTieBreak(candidates, cohort, raw);
      sourceTieBreak = "GetWinPPfromDists";
    }
  }
  const confidence = total > 0 ? maximum / total : 1 / 3;
  const parents = preferredParents(candidates, raw);
  const orientations = candidates.map((candidate, index) => ({
    ...parents[index],
    recombinant: candidate,
    sourcePoints: scores[index],
    sourceScore: maximum > 0 ? Math.round(100 * scores[index] / maximum) : 0,
    sourceShare: total > 0 ? scores[index] / total : 1 / 3,
  }));
  return {
    inference: tree ? "rdp5-source-profile-consensus" : "rdp5-source-distance-consensus",
    candidates,
    recommended: candidates[recommendedIndex],
    recommendedMajorParent: orientations[recommendedIndex].majorParent,
    recommendedMinorParent: orientations[recommendedIndex].minorParent,
    confidence,
    ambiguous: maximumIndexes.length !== 1 || confidence < 0.6,
    sourceThreshold: 0.6,
    orientations,
    tests,
    cohortSize: cohort.length,
    sourceSequenceCount: sequenceCount,
    sampled: cohort.length < sequenceCount,
    treeEvidence: Boolean(tree),
    bootstrapReplicates: treePair?.[0]?.replicates ?? 0,
    bootstrapCutoff: treePair?.[0]?.cutoff ?? 0.5,
    quartetCohortSize: dmaxEvidence?.cohortSize ?? quartetCohort.length,
    quartetCounts: dmaxEvidence?.quartetCounts ?? [0, 0, 0],
    dmaxWasmAccelerated: dmaxEvidence?.wasmAccelerated === true,
    sourceTieBreak,
    sourceTieBreakValues,
    implementedComponents: IMPLEMENTED_COMPONENTS,
    pendingComponents: PENDING_COMPONENTS,
  };
}
