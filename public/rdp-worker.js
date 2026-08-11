import { methodEvidence } from "./rdp-statistics.js";
import { fitBurtTriplet } from "./rdp-burt.js";
import { inferAncestralEventClusters } from "./rdp-clustering.js";
import { buildDisassembledAlignment, candidateComponentProvenance, findComponentIndex, splitCandidateAtStructuralGaps } from "./rdp-disassembly.js";
import { identifyRecombinantRoles } from "./rdp-recombinant-identification.js";
import { buildNeighborJoiningPathMatrix } from "./rdp-bootstrap-tree.js";
import { buildSourceSiScanRandomization, runSourceSiScan } from "./rdp-siscan.js";
import { sourcePhiTest } from "./rdp-phi.js";

let wasmPromise;

const BASES = { A: 0, C: 1, G: 2, T: 3, U: 3, "-": 5 };

function enabledMethodMask(options) {
  const bits = { GENECONV: 1, BootScan: 2, MaxChi: 4, Chimaera: 8, SiScan: 16, "3Seq": 32 };
  let mask = options.polishBreakpoints ? 64 : 0;
  for (const method of options.methods) mask |= bits[method] ?? 0;
  return mask;
}

function align(value, multiple = 8) {
  return Math.ceil(value / multiple) * multiple;
}

async function loadWasm() {
  if (!wasmPromise) {
    wasmPromise = (async () => {
      const imports = {
        env: {
          abort(message, file, line, column) {
            throw new Error(`WebAssembly aborted at ${line}:${column} (${message}/${file})`);
          },
        },
      };
      const workerLocation = typeof self.location?.href === "string"
        ? self.location.href
        : "http://localhost/rdp-worker.js";
      const response = await fetch(new URL("wasm/rdp.wasm", workerLocation));
      if (!response.ok) throw new Error(`Could not load the analysis engine (${response.status}).`);
      try {
        return (await WebAssembly.instantiateStreaming(response.clone(), imports)).instance;
      } catch {
        return (await WebAssembly.instantiate(await response.arrayBuffer(), imports)).instance;
      }
    })();
  }
  return wasmPromise;
}

function encodeSequences(sequences, length) {
  const encoded = new Uint8Array(sequences.length * length);
  encoded.fill(4);
  sequences.forEach((record, sequenceIndex) => {
    for (let site = 0; site < length; site += 1) {
      encoded[sequenceIndex * length + site] = BASES[record.sequence[site]] ?? 4;
    }
  });
  return encoded;
}

function rotateSequences(encoded, sequenceCount, length, rotation) {
  const rotated = new Uint8Array(encoded.length);
  for (let sequence = 0; sequence < sequenceCount; sequence += 1) {
    const source = sequence * length;
    const target = sequence * length;
    rotated.set(encoded.subarray(source + rotation, source + length), target);
    rotated.set(encoded.subarray(source, source + rotation), target + length - rotation);
  }
  return rotated;
}

function packSequences(encoded, sequenceCount, length) {
  const wordsPerSequence = Math.ceil(length / 16);
  const packed = new Uint32Array(sequenceCount * wordsPerSequence);
  const validity = new Uint32Array(sequenceCount * wordsPerSequence);
  for (let sequence = 0; sequence < sequenceCount; sequence += 1) {
    const sequenceOffset = sequence * length;
    const wordOffset = sequence * wordsPerSequence;
    for (let site = 0; site < length; site += 1) {
      const base = encoded[sequenceOffset + site];
      if (base >= 4) continue;
      const word = wordOffset + (site >>> 4);
      const shift = (site & 15) << 1;
      packed[word] |= base << shift;
      validity[word] |= 1 << shift;
    }
  }
  return { packed, validity, wordsPerSequence };
}

function candidateParents(distance, nSeq, target, pool, limit) {
  return pool
    .filter((index) => index !== target)
    .map((index) => ({ index, distance: distance[target * nSeq + index] }))
    .sort((left, right) => left.distance - right.distance)
    .slice(0, Math.max(2, limit))
    .map((entry) => entry.index);
}

function mapInterval(start, end, rotation, length) {
  if (rotation === 0) return { start, end, wraps: start > end };
  const mappedStart = (start + rotation) % length;
  const mappedEnd = (end + rotation) % length;
  if (mappedEnd === 0) return { start: mappedStart, end: length, wraps: false };
  return { start: mappedStart, end: mappedEnd, wraps: mappedStart > mappedEnd };
}

function mapConfidenceInterval(confidence, rotation, length) {
  if (rotation === 0) return confidence;
  const left = (confidence[0] + rotation) % length;
  const right = (confidence[1] + rotation) % length;
  return [left, right];
}

function mapBreakpointModel(model, rotation, length) {
  if (!model || rotation === 0) return model;
  const mapPoint = (position) => (position + rotation) % length;
  return {
    ...model,
    candidateBreakpoints: model.candidateBreakpoints?.map(mapPoint),
    polishedBreakpoints: model.polishedBreakpoints?.map(mapPoint),
    confidence99Start: model.confidence99Start ? mapConfidenceInterval(model.confidence99Start, rotation, length) : undefined,
    confidence99End: model.confidence99End ? mapConfidenceInterval(model.confidence99End, rotation, length) : undefined,
    switches: model.switches?.map((entry) => ({
      ...entry,
      position: mapPoint(entry.position),
      confidence95: mapConfidenceInterval(entry.confidence95, rotation, length),
      confidence99: entry.confidence99 ? mapConfidenceInterval(entry.confidence99, rotation, length) : undefined,
    })),
    posteriorTrace: model.posteriorTrace?.map((entry) => ({ ...entry, position: mapPoint(entry.position) })),
  };
}

function candidateSegments(candidate, length) {
  if (candidate.wraps && candidate.start > candidate.end) {
    return [[candidate.start, length], ...(candidate.end > 0 ? [[0, candidate.end]] : [])];
  }
  return candidate.end > candidate.start ? [[candidate.start, candidate.end]] : [];
}

function tractLength(candidate, length) {
  return candidateSegments(candidate, length)
    .reduce((total, segment) => total + segment[1] - segment[0], 0);
}

function overlap(left, right, length) {
  const leftSegments = candidateSegments(left, length);
  const rightSegments = candidateSegments(right, length);
  let intersection = 0;
  for (const [leftStart, leftEnd] of leftSegments) {
    for (const [rightStart, rightEnd] of rightSegments) {
      intersection += Math.max(0, Math.min(leftEnd, rightEnd) - Math.max(leftStart, rightStart));
    }
  }
  return intersection / Math.max(1, Math.min(tractLength(left, length), tractLength(right, length)));
}

function isCoLocatedSignal(candidate, signal, length, minimumOverlap = 0.25) {
  if (!signal || !(signal.end > signal.start) && !signal.wraps) return false;
  return overlap(candidate, signal, length) >= minimumOverlap;
}

function circularDistance(left, right, length) {
  const distance = Math.abs(left - right);
  return Math.min(distance, length - distance);
}

function sameCircularBreakpoints(left, right, length) {
  const leftParents = [left.majorParent, left.minorParent].sort((a, b) => a - b).join(":");
  const rightParents = [right.majorParent, right.minorParent].sort((a, b) => a - b).join(":");
  if (leftParents !== rightParents) return false;
  const tolerance = Math.max(12, Math.floor(length * 0.02));
  const direct = Math.max(
    circularDistance(left.start, right.start, length),
    circularDistance(left.end % length, right.end % length, length),
  );
  const reversed = Math.max(
    circularDistance(left.start, right.end % length, length),
    circularDistance(left.end % length, right.start, length),
  );
  return Math.min(direct, reversed) <= tolerance;
}

function primaryMethodSignals(signals = []) {
  const best = new Map();
  for (const signal of signals) {
    const previous = best.get(signal.method);
    if (!previous || Math.abs(signal.statistic) > Math.abs(previous.statistic)) best.set(signal.method, signal);
  }
  return [...best.values()];
}

function selectCoLocatedSiScanRegion(result, candidate, rotation, length) {
  if (!result) return null;
  return (result.regions?.length ? result.regions : [result])
    .map((region) => {
      const mapped = mapInterval(region.start, region.end, rotation, length);
      const shared = overlap(candidate, mapped, length);
      return { region, mapped, shared };
    })
    .filter((entry) => entry.shared >= 0.25)
    .sort((left, right) => left.region.rawP - right.region.rawP
      || right.shared - left.shared
      || right.region.z - left.region.z)[0] ?? null;
}

function sourceSiScanProfile(result, rotation, length, maximumPoints = 192) {
  if (!result?.windows?.length) return undefined;
  const stride = Math.max(1, Math.ceil(result.windows.length / maximumPoints));
  const retained = result.windows.filter((_, index) => index % stride === 0 || index === result.windows.length - 1);
  return retained.map((entry) => ({
    position: (entry.center + rotation) % length,
    z: entry.z,
    topology: entry.topology,
    baselineTopology: result.baselineTopology,
    pattern: entry.index,
    scoreFamily: entry.family,
  })).sort((left, right) => left.position - right.position);
}

function addWarnings(candidate, sequences, length, options) {
  const warnings = [];
  if (candidate.structuralUncertainty) {
    warnings.push("RDP5 erased-signal rule split this detection at a prior deleted tract; inspect every gap-adjacent uncertain breakpoint.");
  }
  if (candidate.informative < 25) warnings.push("Low informative-site count");
  const recordIndexes = [candidate.recombinant, candidate.majorParent, candidate.minorParent];
  let ambiguous = 0;
  let examined = 0;
  for (const index of recordIndexes) {
    const sequence = sequences[index].sequence;
    for (const [start, end] of candidateSegments(candidate, length)) {
      for (let site = start; site < end; site += 1) {
        examined += 1;
        if (!(sequence[site] in BASES) || sequence[site] === "-") ambiguous += 1;
      }
    }
  }
  if (options.checkMisalignment && ambiguous / Math.max(1, examined) > 0.035) {
    warnings.push("Gap/ambiguity-rich tract; inspect alignment");
  }
  if (!options.circular && (candidate.start < options.window / 2 || length - candidate.end < options.window / 2)) {
    warnings.push("Breakpoint near an alignment boundary");
  }
  if (candidate.stats.bootscanWindows < 8) warnings.push("Too few decisive windows for stable topology support");
  if (candidate.diagnostics.rateRatio >= 4 || candidate.diagnostics.rateRatio <= 0.25) {
    warnings.push("Strong local variable-site density shift; inspect rate variation or alignment quality");
  }
  if (candidate.diagnostics.parentConflictRate > 0.15) {
    warnings.push("Many tract sites match neither proposed parent; homoplasy or an unsampled parent is plausible");
  }
  if (candidate.diagnostics.diffuseIncompatibility) {
    warnings.push("Widespread four-gamete incompatibility lacks a strong distance trend; homoplasy may contribute");
  }
  return warnings;
}

function bitCount4(value) {
  let count = 0;
  for (let bit = 0; bit < 4; bit += 1) count += (value >>> bit) & 1;
  return count;
}

function firstTwoAlleles(mask) {
  const alleles = [];
  for (let base = 0; base < 4; base += 1) {
    if (mask & (1 << base)) alleles.push(base);
  }
  return alleles;
}

// Bounded, deterministic challenge diagnostics. This is intentionally exposed
// as review evidence rather than used as a hidden filter: four-gamete
// incompatibility can reflect recombination, recurrent mutation, or error.
function alignmentDiagnostics(encoded, sequenceCount, length, window, seed = 0x5a17c0de) {
  const sampleCount = Math.min(256, sequenceCount);
  const sampleIndexes = Array.from({ length: sampleCount }, (_, index) => (
    sampleCount === sequenceCount
      ? index
      : Math.floor(index * (sequenceCount - 1) / Math.max(1, sampleCount - 1))
  ));
  const variablePrefix = new Int32Array(length + 1);
  const biallelic = [];
  let nonCanonical = 0;
  for (let site = 0; site < length; site += 1) {
    let mask = 0;
    for (const sequence of sampleIndexes) {
      const base = encoded[sequence * length + site];
      if (base < 4) mask |= 1 << base;
      else nonCanonical += 1;
    }
    const alleleCount = bitCount4(mask);
    variablePrefix[site + 1] = variablePrefix[site] + (alleleCount > 1 ? 1 : 0);
    if (alleleCount === 2) {
      const [first, second] = firstTwoAlleles(mask);
      biallelic.push({ site, first, second });
    }
  }
  const siteLimit = 96;
  const sampledSites = biallelic.length <= siteLimit
    ? biallelic
    : Array.from({ length: siteLimit }, (_, index) => biallelic[
        Math.floor(index * (biallelic.length - 1) / Math.max(1, siteLimit - 1))
      ]);
  const nearDistance = Math.max(window * 2, Math.floor(length / 50));
  let incompatible = 0;
  let testedPairs = 0;
  let nearIncompatible = 0;
  let nearPairs = 0;
  let farIncompatible = 0;
  let farPairs = 0;
  const pairRows = [];
  for (let left = 0; left < sampledSites.length; left += 1) {
    for (let right = left + 1; right < sampledSites.length; right += 1) {
      const a = sampledSites[left];
      const b = sampledSites[right];
      let gametes = 0;
      for (const sequence of sampleIndexes) {
        const leftBase = encoded[sequence * length + a.site];
        const rightBase = encoded[sequence * length + b.site];
        const leftAllele = leftBase === a.first ? 0 : leftBase === a.second ? 1 : -1;
        const rightAllele = rightBase === b.first ? 0 : rightBase === b.second ? 1 : -1;
        if (leftAllele < 0 || rightAllele < 0) continue;
        gametes |= 1 << (leftAllele * 2 + rightAllele);
        if (gametes === 15) break;
      }
      const isIncompatible = gametes === 15;
      pairRows.push({ left, right, incompatible: isIncompatible });
      testedPairs += 1;
      if (isIncompatible) incompatible += 1;
      if (b.site - a.site <= nearDistance) {
        nearPairs += 1;
        if (isIncompatible) nearIncompatible += 1;
      } else {
        farPairs += 1;
        if (isIncompatible) farIncompatible += 1;
      }
    }
  }
  const nearFraction = nearIncompatible / Math.max(1, nearPairs);
  const farFraction = farIncompatible / Math.max(1, farPairs);
  const proximityRatio = (nearFraction + 1 / Math.max(2, nearPairs))
    / (farFraction + 1 / Math.max(2, farPairs));
  const proximityStatistic = farFraction - nearFraction;
  const permutationReplicates = sampledSites.length >= 4 ? 199 : 0;
  let permutationExtreme = 0;
  let randomState = seed >>> 0;
  const random = () => {
    randomState ^= randomState << 13;
    randomState ^= randomState >>> 17;
    randomState ^= randomState << 5;
    return (randomState >>> 0) / 4294967296;
  };
  const originalPositions = sampledSites.map((site) => site.site);
  for (let replicate = 0; replicate < permutationReplicates; replicate += 1) {
    const permuted = [...originalPositions];
    for (let index = permuted.length - 1; index > 0; index -= 1) {
      const swap = Math.floor(random() * (index + 1));
      [permuted[index], permuted[swap]] = [permuted[swap], permuted[index]];
    }
    let permutedNearPairs = 0;
    let permutedNearIncompatible = 0;
    let permutedFarPairs = 0;
    let permutedFarIncompatible = 0;
    for (const row of pairRows) {
      if (Math.abs(permuted[row.right] - permuted[row.left]) <= nearDistance) {
        permutedNearPairs += 1;
        if (row.incompatible) permutedNearIncompatible += 1;
      } else {
        permutedFarPairs += 1;
        if (row.incompatible) permutedFarIncompatible += 1;
      }
    }
    const statistic = permutedFarIncompatible / Math.max(1, permutedFarPairs)
      - permutedNearIncompatible / Math.max(1, permutedNearPairs);
    if (statistic >= proximityStatistic - 1e-12) permutationExtreme += 1;
  }
  // PHI is calculated separately from the lightweight four-gamete summary:
  // it uses the RDP5/PHIPACK multistate reticulation score and analytic
  // moments. Keep its O(S²N) work bounded for large browser datasets while
  // reporting both retained and total informative-site counts.
  const phiMaximumSites = sequenceCount > 256 ? 160 : sequenceCount > 96 ? 256 : 384;
  const phi = sourcePhiTest(encoded, sequenceCount, length, {
    window: Math.min(length, 100),
    maxInformativeSites: phiMaximumSites,
  });
  const summary = {
    sampledSequences: sampleCount,
    sampledBiallelicSites: sampledSites.length,
    testedSitePairs: testedPairs,
    incompatibleSitePairs: incompatible,
    fourGameteFraction: incompatible / Math.max(1, testedPairs),
    nearIncompatibility: nearFraction,
    farIncompatibility: farFraction,
    proximityRatio,
    proximityStatistic,
    proximityPermutationP: (permutationExtreme + 1) / (permutationReplicates + 1),
    proximityPermutationReplicates: permutationReplicates,
    ambiguityFraction: nonCanonical / Math.max(1, sampleCount * length),
    phiPValue: phi.pValue,
    phiStatistic: phi.statistic,
    phiMean: phi.mean,
    phiVariance: phi.variance,
    phiZ: phi.z ?? 0,
    phiInformativeSites: phi.informativeSites,
    phiTotalInformativeSites: phi.totalInformativeSites,
    phiK: phi.k,
    phiWindow: phi.window,
    phiSubsampled: phi.subsampled,
    phiValidNormalApproximation: phi.validNormalApproximation,
    phiCompatibility: phi.compatibility,
  };
  return { summary, variablePrefix };
}

function candidateDiagnostics(candidate, encoded, length, profile) {
  const segments = candidateSegments(candidate, length);
  const tractSites = tractLength(candidate, length);
  let tractVariable = 0;
  let parentDiscriminating = 0;
  let parentConflicts = 0;
  for (const [start, end] of segments) {
    tractVariable += profile.variablePrefix[end] - profile.variablePrefix[start];
    for (let site = start; site < end; site += 1) {
      const r = encoded[candidate.recombinant * length + site];
      const a = encoded[candidate.majorParent * length + site];
      const b = encoded[candidate.minorParent * length + site];
      if (r >= 4 || a >= 4 || b >= 4 || a === b) continue;
      parentDiscriminating += 1;
      if (r !== a && r !== b) parentConflicts += 1;
    }
  }
  const totalVariable = profile.variablePrefix[length];
  const backgroundSites = Math.max(1, length - tractSites);
  const tractVariableDensity = tractVariable / Math.max(1, tractSites);
  const backgroundVariableDensity = (totalVariable - tractVariable) / backgroundSites;
  const rateRatio = tractVariableDensity / Math.max(1 / Math.max(1, length), backgroundVariableDensity);
  const global = profile.summary;
  const diffuseIncompatibility = global.fourGameteFraction > 0.2
    && global.proximityRatio > 0.75
    && global.proximityRatio < 1.35;
  return {
    tractVariableDensity,
    backgroundVariableDensity,
    rateRatio,
    parentConflictRate: parentConflicts / Math.max(1, parentDiscriminating),
    parentDiscriminatingSites: parentDiscriminating,
    diffuseIncompatibility,
  };
}

function deduplicate(candidates, length, targetCount) {
  const ordered = [...candidates].sort((left, right) =>
    Number(Boolean(right.sourceRdp)) - Number(Boolean(left.sourceRdp))
      || right.chiSquare - left.chiSquare
      || tractLength(left, length) - tractLength(right, length),
  );
  const kept = [];
  for (const candidate of ordered) {
    const duplicate = kept.find((existing) => {
      if (existing.recombinant !== candidate.recombinant) return false;
      const existingLineage = existing.componentProvenance?.recombinant?.lineage?.join("/") ?? "";
      const candidateLineage = candidate.componentProvenance?.recombinant?.lineage?.join("/") ?? "";
      if (existingLineage !== candidateLineage) return false;
      const existingLength = tractLength(existing, length);
      const candidateLength = tractLength(candidate, length);
      const lengthRatio = Math.min(existingLength, candidateLength) / Math.max(1, Math.max(existingLength, candidateLength));
      return sameCircularBreakpoints(existing, candidate, length)
        || (overlap(existing, candidate, length) > 0.75 && lengthRatio > 0.55);
    });
    if (!duplicate) kept.push(candidate);
    else {
      if (!duplicate.alternatives.includes(candidate.minorParent)) duplicate.alternatives.push(candidate.minorParent);
      for (const signal of candidate.methodSignals ?? []) {
        if (!(duplicate.methodSignals ?? []).some((entry) => entry.method === signal.method && entry.start === signal.start && entry.end === signal.end)) {
          duplicate.methodSignals = [...(duplicate.methodSignals ?? []), signal];
        }
      }
    }
  }
  // Keep an adaptive, bounded candidate set so large alignments do not let the
  // first few high-scoring recombinants consume a fixed global quota.
  const maximum = Math.max(500, Math.min(5_000, targetCount * 8));
  const perRecombinant = new Map();
  return kept.filter((candidate) => {
    const count = perRecombinant.get(candidate.recombinant) ?? 0;
    if (count >= 12) return false;
    perRecombinant.set(candidate.recombinant, count + 1);
    return true;
  }).slice(0, maximum);
}

async function analyze(message) {
  const started = performance.now();
  const instance = await loadWasm();
  const { alignment, options, jobId } = message;
  const sequences = alignment.sequences;
  const nSeq = sequences.length;
  const nSites = alignment.length;
  const encoded = encodeSequences(sequences, nSites);
  const disassembly = buildDisassembledAlignment(
    encoded,
    sequences,
    nSites,
    Array.isArray(message.disassemblyEvents) ? message.disassemblyEvents : [],
  );
  const scanEncoded = disassembly.encoded;
  const scanSequences = disassembly.analysisSequences;
  const scanMappings = disassembly.mappings;
  const scanSequenceCount = scanSequences.length;
  const diagnosticsStarted = performance.now();
  const diagnosticProfile = alignmentDiagnostics(encoded, nSeq, nSites, options.window, options.randomSeed);
  const diagnosticsMs = performance.now() - diagnosticsStarted;
  const rotation = options.circular ? Math.floor(nSites / 2) : 0;
  const rotated = rotation > 0 ? rotateSequences(scanEncoded, scanSequenceCount, nSites, rotation) : null;
  const scanPacked = packSequences(scanEncoded, scanSequenceCount, nSites);
  const originalPacked = disassembly.componentCount > 0
    ? packSequences(encoded, nSeq, nSites)
    : scanPacked;
  const wordsPerSequence = scanPacked.wordsPerSequence;
  const exactDistanceMatrix = scanSequenceCount <= 512
    && scanSequenceCount * scanSequenceCount * wordsPerSequence <= 50_000_000;
  const exactDisplayMatrix = nSeq <= 512
    && nSeq * nSeq * originalPacked.wordsPerSequence <= 50_000_000;
  const parentSamples = scanSequenceCount > 2_000 ? 64 : scanSequenceCount > 1_000 ? 128 : 256;
  const scanMatrixCount = exactDistanceMatrix ? scanSequenceCount : Math.min(24, scanSequenceCount);
  const matrixCount = exactDisplayMatrix ? nSeq : Math.min(24, nSeq);
  const seqPtr = 65536;
  const rotatedSeqPtr = rotated ? align(seqPtr + scanEncoded.byteLength, 16) : seqPtr;
  const sequenceEndPtr = rotated ? rotatedSeqPtr + rotated.byteLength : seqPtr + scanEncoded.byteLength;
  const packedPtr = align(sequenceEndPtr, 16);
  const validityPtr = align(packedPtr + scanPacked.packed.byteLength, 16);
  const originalPackedPtr = align(validityPtr + scanPacked.validity.byteLength, 16);
  const originalValidityPtr = align(originalPackedPtr + originalPacked.packed.byteLength, 16);
  const distancePtr = align(originalValidityPtr + originalPacked.validity.byteLength, 16);
  const distanceBytes = Math.max(scanMatrixCount * scanMatrixCount, matrixCount * matrixCount) * 4;
  const prefixAPtr = align(distancePtr + distanceBytes, 16);
  const prefixBPtr = prefixAPtr + (nSites + 1) * 4;
  const outPtr = align(prefixBPtr + (nSites + 1) * 4, 16);
  const rdpSignalCapacity = Math.max(1, Math.min(256, Math.trunc(options.rdpSignalsPerTriplet ?? 128)));
  const rdpBestPtr = outPtr + rdpSignalCapacity * 72;
  const statsPtr = align(rdpBestPtr + 72, 16);
  const poolPtr = align(statsPtr + 160, 16);
  const nearestIndexesPtr = poolPtr + scanSequenceCount * 4;
  const nearestDistancesPtr = nearestIndexesPtr + scanSequenceCount * 4;
  const roleCohortCapacity = Math.max(4, Math.min(34, Math.trunc(options.roleQuartetTaxaLimit ?? 30)));
  const roleCohortPtr = align(nearestDistancesPtr + scanSequenceCount * 4, 16);
  const roleTractMaskPtr = roleCohortPtr + roleCohortCapacity * 4;
  const roleBackgroundMaskPtr = roleTractMaskPtr + wordsPerSequence * 4;
  const roleDmaxOutPtr = align(roleBackgroundMaskPtr + wordsPerSequence * 4, 8);
  const requiredBytes = roleDmaxOutPtr + 40;
  const memory = instance.exports.memory;
  const requiredPages = Math.ceil(requiredBytes / 65536);
  const currentPages = memory.buffer.byteLength / 65536;
  if (requiredPages > currentPages) memory.grow(requiredPages - currentPages);
  new Uint8Array(memory.buffer, seqPtr, scanEncoded.byteLength).set(scanEncoded);
  if (rotated) new Uint8Array(memory.buffer, rotatedSeqPtr, rotated.byteLength).set(rotated);
  new Uint32Array(memory.buffer, packedPtr, scanPacked.packed.length).set(scanPacked.packed);
  new Uint32Array(memory.buffer, validityPtr, scanPacked.validity.length).set(scanPacked.validity);
  new Uint32Array(memory.buffer, originalPackedPtr, originalPacked.packed.length).set(originalPacked.packed);
  new Uint32Array(memory.buffer, originalValidityPtr, originalPacked.validity.length).set(originalPacked.validity);

  const roleDmaxEvaluator = ({ event, candidates: roleCandidates, cohort }) => {
    if (typeof instance.exports.dmax_visrd_packed !== "function" || cohort.length > roleCohortCapacity) return null;
    const tractMask = new Uint32Array(wordsPerSequence);
    const backgroundMask = new Uint32Array(wordsPerSequence);
    for (let site = 0; site < nSites; site += 1) {
      const inside = event.wraps && event.start > event.end
        ? site >= event.start || site < event.end
        : site >= event.start && site < event.end;
      const word = site >>> 4;
      const lane = 1 << ((site & 15) << 1);
      if (inside) tractMask[word] |= lane;
      else backgroundMask[word] |= lane;
    }
    new Int32Array(memory.buffer, roleCohortPtr, cohort.length).set(cohort);
    new Uint32Array(memory.buffer, roleTractMaskPtr, wordsPerSequence).set(tractMask);
    new Uint32Array(memory.buffer, roleBackgroundMaskPtr, wordsPerSequence).set(backgroundMask);
    instance.exports.dmax_visrd_packed(
      packedPtr,
      validityPtr,
      wordsPerSequence,
      roleCohortPtr,
      cohort.length,
      roleCandidates[0],
      roleCandidates[1],
      roleCandidates[2],
      roleTractMaskPtr,
      roleBackgroundMaskPtr,
      roleDmaxOutPtr,
    );
    return {
      values: [...new Float64Array(memory.buffer, roleDmaxOutPtr, 3)],
      quartetCounts: [...new Int32Array(memory.buffer, roleDmaxOutPtr + 24, 3)],
      cohortSize: cohort.length,
      sourceRoutine: "CalcMaxD + CMaxD2P3",
      wasmAccelerated: true,
    };
  };

  const distanceStarted = performance.now();
  instance.exports.distance_matrix_packed(
    packedPtr,
    validityPtr,
    scanMatrixCount,
    wordsPerSequence,
    distancePtr,
  );
  const parentDistance = new Float32Array(memory.buffer, distancePtr, scanMatrixCount * scanMatrixCount).slice();
  instance.exports.distance_matrix_packed(
    originalPackedPtr,
    originalValidityPtr,
    matrixCount,
    originalPacked.wordsPerSequence,
    distancePtr,
  );
  const distance = new Float32Array(memory.buffer, distancePtr, matrixCount * matrixCount).slice();
  let sourceMaximumDistance = exactDisplayMatrix
    ? distance.reduce((maximum, value) => Math.max(maximum, value), 0)
    : instance.exports.maximum_packed_distance(
        originalPackedPtr,
        originalValidityPtr,
        nSeq,
        originalPacked.wordsPerSequence,
      );
  if (!(sourceMaximumDistance > 0)) sourceMaximumDistance = 1;
  // GetSSOL first attempts the desktop program's whole-alignment tree-path
  // outlier rule.  Deterministic NJ is cubic, so retain that exact path for
  // cohorts where it is cheaper than the subsequent triplet scan; larger
  // cohorts use GetSSOL's own direct-distance fallback.
  const sourceSiScanTreeDistance = options.methods.includes("SiScan")
    && exactDistanceMatrix
    && scanSequenceCount <= 128
    ? buildNeighborJoiningPathMatrix(parentDistance, scanSequenceCount)
    : null;
  const distanceMs = performance.now() - distanceStarted;
  const allIndexes = Array.from({ length: scanSequenceCount }, (_, index) => index);
  const excludedTargets = new Set(Array.isArray(message.excludedTargets) ? message.excludedTargets : []);
  const excludedParents = new Set(Array.isArray(message.excludedParents) ? message.excludedParents : []);
  const targetAllowed = (index) => scanMappings[index].kind === "extracted-tract"
    || !excludedTargets.has(scanMappings[index].originIndex);
  const parentAllowed = (index) => scanMappings[index].kind === "extracted-tract"
    || !excludedParents.has(scanMappings[index].originIndex);
  const targets = options.mode === "query-reference"
    ? allIndexes.filter((index) => targetAllowed(index) && (options.testReferences || scanSequences[index].role === "query" || scanSequences[index].role === "both"))
    : allIndexes.filter(targetAllowed);
  const referencePool = options.mode === "query-reference"
    ? allIndexes.filter((index) => parentAllowed(index) && (scanSequences[index].role === "reference" || scanSequences[index].role === "both"))
    : allIndexes.filter(parentAllowed);
  new Int32Array(memory.buffer, poolPtr, referencePool.length).set(referencePool);
  const candidates = [];
  const independentMethods = options.methods.filter((method) => method !== "RDP");
  const independentMethodMask = enabledMethodMask({
    ...options,
    methods: independentMethods,
    polishBreakpoints: false,
  });
  // MAXCHI and CHIMAERA are invariant to swapping the two proposed parents.
  // The remaining locators are directional, so only those bits need the
  // reverse-parent pass.
  const reverseIndependentMethodMask = independentMethodMask & (1 | 2 | 16 | 32);
  const partialBest = new Map();
  const retainCandidate = (candidate) => {
    candidate.structuralUncertaintyVnps = Math.max(1, Math.trunc(options.rdpWindow ?? 30));
    candidate.circular = options.circular === true;
    for (const piece of splitCandidateAtStructuralGaps(candidate, disassembly, nSites)) {
      candidates.push(piece);
      const previousPartial = partialBest.get(piece.analysisRecombinant);
      if (!previousPartial || piece.chiSquare > previousPartial.chiSquare) {
        partialBest.set(piece.analysisRecombinant, piece);
      }
    }
  };
  let comparisons = 0;
  const comparisonKeys = new Set();
  const rdpTripletCache = new Map();
  const rdpTripletCacheLimit = 4_096;
  const countedTruncations = new Set();
  let truncatedRdpSignals = 0;
  const sourceSiScanCache = new Map();
  const sourceSiScanCacheLimit = 512;
  let sourceSiScanRandomization = null;
  const sourceSiScanPermutations = Math.max(
    2,
    Math.trunc(options.siskanPValuePermutations ?? 1000),
    Math.trunc(options.siskanScanPermutations ?? 100),
  );
  const getSourceSiScanRandomization = () => {
    if (!sourceSiScanRandomization) {
      sourceSiScanRandomization = buildSourceSiScanRandomization(
        nSites,
        sourceSiScanPermutations,
        options.randomSeed ?? 0x5a17c0de,
      );
    }
    return sourceSiScanRandomization;
  };
  const scanStarted = performance.now();
  const scanViews = rotated
    ? [{ sequencePtr: seqPtr, rotation: 0 }, { sequencePtr: rotatedSeqPtr, rotation }]
    : [{ sequencePtr: seqPtr, rotation: 0 }];

  for (let targetPosition = 0; targetPosition < targets.length; targetPosition += 1) {
    const analysisRecombinant = targets[targetPosition];
    const recombinant = scanMappings[analysisRecombinant].originIndex;
    const parentLimit = options.exhaustive ? referencePool.length : Math.min(options.candidateParents, referencePool.length);
    let parents;
    if (options.exhaustive) {
      parents = referencePool.filter((index) => index !== analysisRecombinant);
    } else if (exactDistanceMatrix) {
      parents = candidateParents(parentDistance, scanSequenceCount, analysisRecombinant, referencePool, parentLimit);
    } else {
      const nearestLimit = Math.max(2, Math.ceil(parentLimit * 0.625));
      const foundParents = instance.exports.nearest_candidates_sampled(
        seqPtr,
        nSites,
        analysisRecombinant,
        poolPtr,
        referencePool.length,
        Math.min(parentSamples, nSites),
        nearestLimit,
        nearestIndexesPtr,
        nearestDistancesPtr,
      );
      parents = Array.from(new Int32Array(memory.buffer, nearestIndexesPtr, foundParents));
      const stratifiedSlots = Math.max(0, parentLimit - parents.length);
      for (let slot = 0; slot < stratifiedSlots; slot += 1) {
        const position = Math.min(
          referencePool.length - 1,
          Math.floor(((slot + 0.5) * referencePool.length) / Math.max(1, stratifiedSlots)),
        );
        const candidate = referencePool[position];
        if (candidate !== analysisRecombinant && !parents.includes(candidate)) parents.push(candidate);
      }
      for (const candidate of referencePool) {
        if (parents.length >= parentLimit) break;
        if (candidate !== analysisRecombinant && !parents.includes(candidate)) parents.push(candidate);
      }
    }
    // Explicit reference groups reserve a bounded share of the parent shortlist
    // so a dense group cannot crowd every other named lineage out of the scan.
    const groupRepresentatives = new Map();
    for (const candidate of referencePool) {
      const group = scanSequences[candidate].referenceGroup;
      if (!group || candidate === analysisRecombinant || groupRepresentatives.has(group)) continue;
      groupRepresentatives.set(group, candidate);
    }
    const reserve = Math.min(groupRepresentatives.size, Math.max(0, Math.floor(parentLimit / 3)));
    let groupSlot = 0;
    for (const representative of groupRepresentatives.values()) {
      if (groupSlot >= reserve) break;
      if (!parents.includes(representative)) parents[Math.max(0, parents.length - 1 - groupSlot)] = representative;
      groupSlot += 1;
    }
    parents = [...new Set(parents)].filter((index) => index !== analysisRecombinant);
    for (let left = 0; left < parents.length; left += 1) {
      for (let right = left + 1; right < parents.length; right += 1) {
        const inputMajorParent = parents[left];
        const inputMinorParent = parents[right];
        if (options.mode === "query-reference") {
          const firstGroup = scanSequences[inputMajorParent].referenceGroup?.trim();
          const secondGroup = scanSequences[inputMinorParent].referenceGroup?.trim();
          if (firstGroup && secondGroup && firstGroup === secondGroup) continue;
        }
        const inputOrigins = [
          recombinant,
          scanMappings[inputMajorParent].originIndex,
          scanMappings[inputMinorParent].originIndex,
        ];
        if (new Set(inputOrigins).size < 3) continue;
        const comparisonKey = options.mode !== "query-reference"
          ? [analysisRecombinant, inputMajorParent, inputMinorParent].sort((a, b) => a - b).join(":")
          : `${analysisRecombinant}:${[inputMajorParent, inputMinorParent].sort((a, b) => a - b).join(":")}`;
        if (!comparisonKeys.has(comparisonKey)) {
          comparisonKeys.add(comparisonKey);
          comparisons += 1;
        }
        for (const view of scanViews) {
          if (options.methods.includes("RDP")) {
            const rdpCacheKey = `${[analysisRecombinant, inputMajorParent, inputMinorParent].sort((a, b) => a - b).join(":")}@${view.rotation}`;
            let cachedRdp = rdpTripletCache.get(rdpCacheKey);
            if (!cachedRdp) {
              const totalRdpSignals = instance.exports.scan_rdp5_triplet_all(
                view.sequencePtr,
                nSites,
                analysisRecombinant,
                inputMajorParent,
                inputMinorParent,
                Math.max(5, options.rdpWindow ?? 30),
                prefixBPtr,
                prefixAPtr,
                outPtr,
                rdpSignalCapacity,
                rdpBestPtr,
              );
              const retainedRdpSignals = Math.min(rdpSignalCapacity, Math.max(0, totalRdpSignals));
              cachedRdp = {
                total: totalRdpSignals,
                outputs: Array.from({ length: retainedRdpSignals }, (_, signalIndex) => (
                  Array.from(new Int32Array(memory.buffer, outPtr + signalIndex * 72, 18))
                )),
              };
              rdpTripletCache.set(rdpCacheKey, cachedRdp);
              if (rdpTripletCache.size > rdpTripletCacheLimit) {
                rdpTripletCache.delete(rdpTripletCache.keys().next().value);
              }
            }
            if (cachedRdp.total > rdpSignalCapacity && !countedTruncations.has(rdpCacheKey)) {
              truncatedRdpSignals += cachedRdp.total - rdpSignalCapacity;
              countedTruncations.add(rdpCacheKey);
            }
            for (const output of cachedRdp.outputs) {
              // RDP can assign any triplet member as the daughter. This pass is
              // scoped to one requested target, so defer the other polarities
              // to their own target passes.
              if (output[2] === analysisRecombinant) {
                const rawStart = output[0];
                const rawEnd = output[1];
                const mapped = mapInterval(rawStart, rawEnd, view.rotation, nSites);
                const majorParent = scanMappings[output[3]].originIndex;
                const minorParent = scanMappings[output[4]].originIndex;
                const candidate = {
                  recombinant,
                  ...mapped,
                  rawStart,
                  rawEnd,
                  sequencePtr: view.sequencePtr,
                  rotation: view.rotation,
                  majorParent,
                  minorParent,
                  analysisRecombinant: output[2],
                  analysisMajorParent: output[3],
                  analysisMinorParent: output[4],
                  siskanCandidatePool: exactDistanceMatrix ? undefined : [...parents],
                  componentProvenance: candidateComponentProvenance(disassembly, output[2], output[3], output[4]),
                  chiSquare: output[5] / 1000,
                  informative: output[6],
                  insideMinor: output[7],
                  insideMajor: output[8],
                  outsideMajor: output[9],
                  outsideMinor: output[10],
                  effect: output[11] / 1e6,
                  alternatives: [],
                  sourceRdp: {
                    common: output[12],
                    tractSites: output[13],
                    mediumSites: output[14],
                    probabilitySites: Math.max(1, output[6] - 1),
                    orientationCycle: output[15],
                    window: output[16],
                    compatibility: "RDP5 FindSubSeqPB3/XOHomologyP2/FindNextP/DefineEventP2",
                  },
                  methodSignals: [{ method: "RDP", start: mapped.start, end: mapped.end, wraps: mapped.wraps, statistic: output[5] / 1000, locator: "RDP5 source" }],
                };
                retainCandidate(candidate);
              }
            }
          }

          // Each enabled non-RDP family now contributes its own discovery
          // interval. This replaces the former shared CUSUM seed, which could
          // not discover a signal seen only by one method and could attach a
          // method's strongest peak to an unrelated RDP event.
          if (independentMethodMask !== 0) {
            const parentOrientations = [
              [inputMajorParent, inputMinorParent, independentMethodMask],
              [inputMinorParent, inputMajorParent, reverseIndependentMethodMask],
            ];
            for (let orientation = 0; orientation < parentOrientations.length; orientation += 1) {
              const [analysisMajorParent, analysisMinorParent, methodMask] = parentOrientations[orientation];
              if (methodMask === 0) continue;
              instance.exports.method_stats(
                view.sequencePtr,
                nSites,
                analysisRecombinant,
                analysisMajorParent,
                analysisMinorParent,
                0,
                nSites,
                Math.max(20, options.window),
                Math.max(1, options.step),
                0,
                (
                  (options.randomSeed ?? 0x5a17c0de)
                  ^ Math.imul(analysisRecombinant + 1, 0x9e3779b1)
                  ^ Math.imul(analysisMajorParent + 1, 0x85ebca6b)
                  ^ Math.imul(analysisMinorParent + 1, 0xc2b2ae35)
                ) | 0,
                methodMask,
                prefixAPtr,
                prefixBPtr,
                statsPtr,
                options.geneconvGScale ?? 1,
              );
              const output = new Int32Array(memory.buffer, statsPtr, 35).slice();
              const signals = [
                { method: "GENECONV", start: output[3], end: output[4], statistic: output[0], present: output[0] > 0, locator: `independent RDP5 G=${options.geneconvGScale ?? 1} fragment` },
                { method: "BootScan", start: output[29], end: output[30], statistic: output[33], present: output[33] > 0, locator: "independent topology-window run" },
                { method: "MaxChi", start: output[25], end: output[26], statistic: output[7] / 1000, present: output[7] > 0, locator: "independent paired χ² peaks" },
                { method: "Chimaera", start: output[27], end: output[28], statistic: output[8] / 1000, present: output[8] > 0, locator: "independent binary χ² peaks" },
                { method: "SiScan", start: output[31], end: output[32], statistic: output[34], present: output[34] > 0, locator: "independent oriented-category run" },
                { method: "3Seq", start: output[23], end: output[24], statistic: output[11], present: output[11] > 0, locator: "independent maximum HGRW descent" },
              ];
              for (const signal of signals) {
                if (!signal.present || !independentMethods.includes(signal.method)) continue;
                const rawStart = signal.start;
                const rawEnd = signal.end;
                if (!(rawEnd - rawStart >= 4 && rawStart >= 0 && rawEnd <= nSites)) continue;
                instance.exports.triplet_counts(
                  view.sequencePtr,
                  nSites,
                  analysisRecombinant,
                  analysisMajorParent,
                  analysisMinorParent,
                  rawStart,
                  rawEnd,
                  rdpBestPtr,
                );
                const counts = new Int32Array(memory.buffer, rdpBestPtr, 6).slice();
                const mapped = mapInterval(rawStart, rawEnd, view.rotation, nSites);
                const insideTotal = counts[1] + counts[2];
                const outsideTotal = counts[3] + counts[4];
                const effect = counts[1] / Math.max(1, insideTotal) - counts[4] / Math.max(1, outsideTotal);
                retainCandidate({
                  recombinant,
                  ...mapped,
                  rawStart,
                  rawEnd,
                  sequencePtr: view.sequencePtr,
                  rotation: view.rotation,
                  majorParent: scanMappings[analysisMajorParent].originIndex,
                  minorParent: scanMappings[analysisMinorParent].originIndex,
                  analysisRecombinant,
                  analysisMajorParent,
                  analysisMinorParent,
                  siskanCandidatePool: exactDistanceMatrix ? undefined : [...parents],
                  componentProvenance: candidateComponentProvenance(disassembly, analysisRecombinant, analysisMajorParent, analysisMinorParent),
                  chiSquare: counts[5] / 1000,
                  informative: counts[0],
                  insideMinor: counts[1],
                  insideMajor: counts[2],
                  outsideMajor: counts[3],
                  outsideMinor: counts[4],
                  effect,
                  alternatives: [],
                  methodSignals: [{
                    method: signal.method,
                    start: mapped.start,
                    end: mapped.end,
                    wraps: mapped.wraps,
                    statistic: signal.statistic,
                    locator: signal.locator,
                  }],
                });
              }
            }
          }
        }
      }
    }
    postMessage({
      type: "progress",
      jobId,
      progress: 0.85 * (targetPosition + 1) / targets.length,
      phase: `Scanning ${scanSequences[analysisRecombinant].name}`,
    });
    const partialInterval = Math.max(1, Math.floor(targets.length / 25));
    if ((targetPosition + 1) % partialInterval === 0 || targetPosition === targets.length - 1) {
      const confidence = Math.max(4, Math.floor(options.window / 8));
      const partialEvents = [...partialBest.values()]
        .sort((left, right) => right.chiSquare - left.chiSquare)
        .slice(0, 100)
        .map((candidate, index) => ({
          id: `partial-${jobId}-${candidate.recombinant}-${index + 1}`,
          recombinant: candidate.recombinant,
          majorParent: candidate.majorParent,
          minorParent: candidate.minorParent,
          start: candidate.start,
          end: candidate.end,
          wraps: candidate.wraps,
          confidenceStart: [Math.max(0, candidate.start - confidence), Math.min(nSites, candidate.start + confidence)],
          confidenceEnd: [Math.max(0, candidate.end - confidence), Math.min(nSites, candidate.end + confidence)],
          breakpointModel: { method: "local-chi-square", informativeSites: candidate.informative },
          evidence: [],
          chiSquare: candidate.chiSquare,
          informativeSites: candidate.informative,
          decision: "unreviewed",
          warnings: ["Recovered partial candidate: method evidence and experiment-wide correction are incomplete. Rerun or recalculate before review."],
          note: "Checkpointed while the scan was still running.",
          source: "wasm",
          supportedCount: 0,
          diagnostics: { tractVariableDensity: 0, backgroundVariableDensity: 0, rateRatio: 1, parentConflictRate: 0, parentDiscriminatingSites: 0, diffuseIncompatibility: false },
          groupId: null,
          alternativeParents: [],
          hypothesisTests: Math.max(1, comparisons),
          history: [{ id: `partial-history-${jobId}-${index + 1}`, timestamp: new Date().toISOString(), action: "Checkpointed partial candidate", summary: "Candidate discovery completed for a subset of target sequences; evidence calibration did not complete." }],
          evidenceStale: true,
          componentProvenance: candidate.componentProvenance,
          structuralUncertainty: candidate.structuralUncertainty,
        }));
      postMessage({ type: "partial", jobId, events: partialEvents, comparisons });
    }
  }
  const scanMs = performance.now() - scanStarted;

  const unique = deduplicate(candidates, nSites, targets.length);
  const calibratedCandidateLimit = Math.max(500, Math.min(5_000, targets.length * 8));
  const statisticsStarted = performance.now();
  for (let candidateIndex = 0; candidateIndex < unique.length; candidateIndex += 1) {
    const candidate = unique[candidateIndex];
    if (candidate.structuralUncertainty || candidate.needsTripletRecount) {
      instance.exports.triplet_counts(
        candidate.sequencePtr,
        nSites,
        candidate.analysisRecombinant,
        candidate.analysisMajorParent,
        candidate.analysisMinorParent,
        candidate.rawStart,
        candidate.rawEnd,
        outPtr,
      );
      const counts = new Int32Array(memory.buffer, outPtr, 6);
      candidate.informative = counts[0];
      candidate.insideMinor = counts[1];
      candidate.insideMajor = counts[2];
      candidate.outsideMajor = counts[3];
      candidate.outsideMinor = counts[4];
      candidate.chiSquare = counts[5] / 1000;
      candidate.effect = candidate.insideMinor / Math.max(1, candidate.insideMinor + candidate.insideMajor)
        - candidate.outsideMinor / Math.max(1, candidate.outsideMinor + candidate.outsideMajor);
    }
    instance.exports.method_stats(
      candidate.sequencePtr,
      nSites,
      candidate.analysisRecombinant,
      candidate.analysisMajorParent,
      candidate.analysisMinorParent,
      candidate.rawStart,
      candidate.rawEnd,
      Math.max(20, options.window),
      Math.max(1, options.step),
      options.methods.includes("BootScan") ? Math.max(0, options.bootstrapReplicates ?? 100) : 0,
      (
        (options.randomSeed ?? 0x5a17c0de)
        ^ Math.imul(candidate.analysisRecombinant + 1, 0x9e3779b1)
        ^ Math.imul(candidate.analysisMajorParent + 1, 0x85ebca6b)
        ^ Math.imul(candidate.analysisMinorParent + 1, 0xc2b2ae35)
        ^ candidate.rawStart
        ^ Math.imul(candidate.rawEnd, 31)
      ) | 0,
      enabledMethodMask(options),
      prefixAPtr,
      prefixBPtr,
      statsPtr,
      options.geneconvGScale ?? 1,
    );
    const methodOutput = new Int32Array(memory.buffer, statsPtr, 35);
    candidate.stats = {
      genconvRun: methodOutput[0],
      genconvEligible: methodOutput[1],
      genconvMatches: methodOutput[2],
      genconvStart: methodOutput[3],
      genconvEnd: methodOutput[4],
      bootscanConsistent: methodOutput[5],
      bootscanWindows: methodOutput[6],
      maxChi: methodOutput[7] / 1000,
      chimaera: methodOutput[8] / 1000,
      siskanScore: methodOutput[9],
      siskanSites: methodOutput[10],
      threeSeqDescent: methodOutput[11],
      threeSeqSites: methodOutput[12],
      maxChiBoundaries: [methodOutput[13] / 1000, methodOutput[14] / 1000],
      chimaeraBoundaries: [methodOutput[15] / 1000, methodOutput[16] / 1000],
      threeSeqMajorSites: methodOutput[19],
      threeSeqMinorSites: methodOutput[20],
      bootscanBootstrapConsistent: methodOutput[21],
      bootscanBootstrapReplicates: methodOutput[22],
      threeSeqStart: methodOutput[23],
      threeSeqEnd: methodOutput[24],
      maxChiStart: methodOutput[25],
      maxChiEnd: methodOutput[26],
      chimaeraStart: methodOutput[27],
      chimaeraEnd: methodOutput[28],
      bootscanStart: methodOutput[29],
      bootscanEnd: methodOutput[30],
      siskanStart: methodOutput[31],
      siskanEnd: methodOutput[32],
      bootscanRunWindows: methodOutput[33],
      siskanRunWindows: methodOutput[34],
      rdpSource: candidate.sourceRdp,
    };
    const mapSignal = (method, start, end, statistic, locator) => {
      if (!(end > start)) return;
      const mapped = mapInterval(start, end, candidate.rotation, nSites);
      const signal = { method, ...mapped, statistic, locator };
      if (!isCoLocatedSignal(candidate, signal, nSites)) return;
      if (!(candidate.methodSignals ?? []).some((entry) => entry.method === method && entry.start === signal.start && entry.end === signal.end)) {
        candidate.methodSignals = [...(candidate.methodSignals ?? []), signal];
      }
    };
    if (options.methods.includes("GENECONV")) mapSignal("GENECONV", methodOutput[3], methodOutput[4], methodOutput[0], "maximum concordant fragment");
    if (options.methods.includes("BootScan")) mapSignal("BootScan", methodOutput[29], methodOutput[30], methodOutput[21] / Math.max(1, methodOutput[22]), "minor-topology window run");
    if (options.methods.includes("MaxChi")) mapSignal("MaxChi", methodOutput[25], methodOutput[26], methodOutput[7] / 1000, "independent χ² peak pair");
    if (options.methods.includes("Chimaera")) mapSignal("Chimaera", methodOutput[27], methodOutput[28], methodOutput[8] / 1000, "independent binary χ² peak pair");
    if (options.methods.includes("SiScan")) {
      const preliminary = methodOutput[32] > methodOutput[31]
        ? { method: "SiScan", ...mapInterval(methodOutput[31], methodOutput[32], candidate.rotation, nSites) }
        : null;
      const sourceRequested = (candidate.methodSignals ?? []).some((signal) => signal.method === "SiScan")
        || isCoLocatedSignal(candidate, preliminary, nSites);
      // The WASM category run is a deliberately cheap locator.  It may seed a
      // candidate, but it never supplies final SiScan evidence: the source
      // 15-category/outgroup/permutation workflow must confirm a co-located
      // topology run first.
      candidate.methodSignals = (candidate.methodSignals ?? []).filter((signal) => signal.method !== "SiScan");
      candidate.stats.siskanScore = 0;
      candidate.stats.siskanSites = 0;
      candidate.stats.siskanSourceP = 1;
      candidate.stats.siskanSourceZ = 0;
      if (sourceRequested) {
        const triplet = [candidate.analysisRecombinant, candidate.analysisMajorParent, candidate.analysisMinorParent];
        const tripletOrigins = new Set(triplet.map((index) => scanMappings[index].originIndex));
        const candidatePool = (candidate.siskanCandidatePool ?? allIndexes)
          .filter((index) => !tripletOrigins.has(scanMappings[index].originIndex));
        const manualOutgroup = options.siskanOutgroupMode === "manual"
          ? allIndexes.find((index) => scanMappings[index].originIndex === options.siskanOutgroupSequence
            && !triplet.includes(index)
            && !tripletOrigins.has(scanMappings[index].originIndex))
          : undefined;
        const cacheKey = `${triplet.join(":")}@${candidate.rotation}:${candidatePool.join(",")}`;
        let sourceResult = sourceSiScanCache.get(cacheKey);
        if (sourceResult === undefined) {
          const seed = (
            (options.randomSeed ?? 0x5a17c0de)
            ^ Math.imul(candidate.analysisRecombinant + 1, 0x9e3779b1)
            ^ Math.imul(candidate.analysisMajorParent + 1, 0x85ebca6b)
            ^ Math.imul(candidate.analysisMinorParent + 1, 0xc2b2ae35)
          ) >>> 0;
          sourceResult = runSourceSiScan(
            candidate.rotation === 0 ? scanEncoded : rotated,
            nSites,
            scanSequenceCount,
            triplet,
            {
              window: Math.max(20, options.window),
              step: Math.max(1, options.step),
              scanPermutations: options.siskanScanPermutations ?? 100,
              pValuePermutations: options.siskanPValuePermutations ?? 1000,
              seed,
              outgroupMode: options.siskanOutgroupMode === "manual" && manualOutgroup === undefined
                ? "nearest"
                : options.siskanOutgroupMode ?? "nearest",
              outgroupIndex: manualOutgroup,
              positionMode: options.siskanPositionMode ?? "triplet-variable",
              gapsAsState: options.siskanGapMode === "fifth-state",
              candidatePool,
              distanceMatrix: exactDistanceMatrix ? parentDistance : null,
              treeDistanceMatrix: sourceSiScanTreeDistance,
              randomization: getSourceSiScanRandomization(),
            },
          );
          sourceSiScanCache.set(cacheKey, sourceResult);
          if (sourceSiScanCache.size > sourceSiScanCacheLimit) {
            sourceSiScanCache.delete(sourceSiScanCache.keys().next().value);
          }
        }
        // The source routine can return several disjoint topology runs for one
        // ordered triplet.  Queue every locally significant run once, bounded
        // by the same per-recombinant/global retention limits as discovery.
        // This preserves events that a single "best run" locator would hide.
        if (sourceResult?.regions?.length && unique.length < calibratedCandidateLimit) {
          let recombinantCount = unique.filter((entry) => entry.recombinant === candidate.recombinant).length;
          for (const region of sourceResult.regions) {
            if (region.rawP > Math.max(1e-300, options.alpha ?? 0.05) || recombinantCount >= 12 || unique.length >= calibratedCandidateLimit) continue;
            const mapped = mapInterval(region.start, region.end, candidate.rotation, nSites);
            const duplicateRegion = unique.some((entry) => (
              entry.rotation === candidate.rotation
              && entry.analysisRecombinant === candidate.analysisRecombinant
              && entry.analysisMajorParent === candidate.analysisMajorParent
              && entry.analysisMinorParent === candidate.analysisMinorParent
              && overlap(entry, mapped, nSites) > 0.75
            ));
            if (duplicateRegion) continue;
            const sourceSignal = {
              method: "SiScan",
              ...mapped,
              statistic: region.z,
              locator: `RDP5 Sister-Scanning ${region.scoreFamily} ${region.pattern} topology run`,
            };
            const expanded = splitCandidateAtStructuralGaps({
              ...candidate,
              ...mapped,
              rawStart: region.start,
              rawEnd: region.end,
              sourceRdp: undefined,
              breakpointModel: undefined,
              confidenceStart: undefined,
              confidenceEnd: undefined,
              structuralUncertainty: undefined,
              needsTripletRecount: true,
              methodSignals: [sourceSignal],
              alternatives: [...(candidate.alternatives ?? [])],
            }, disassembly, nSites);
            for (const piece of expanded) {
              if (recombinantCount >= 12 || unique.length >= calibratedCandidateLimit) break;
              unique.push(piece);
              recombinantCount += 1;
            }
          }
        }
        const selected = selectCoLocatedSiScanRegion(sourceResult, candidate, candidate.rotation, nSites);
        if (selected) {
          const { region, mapped } = selected;
          candidate.stats.siskanScore = region.z;
          candidate.stats.siskanSites = region.end - region.start;
          candidate.stats.siskanStart = region.start;
          candidate.stats.siskanEnd = region.end;
          candidate.stats.siskanSourceP = region.rawP;
          candidate.stats.siskanSourceZ = region.z;
          candidate.stats.siskanOutgroupIndex = sourceResult.outgroupIndex;
          candidate.stats.siskanOutgroupMode = sourceResult.outgroupMode;
          candidate.stats.siskanOutgroupSampled = sourceResult.outgroupSampled;
          candidate.stats.siskanOutgroupSourcePath = sourceResult.outgroupSourcePath;
          candidate.stats.siskanPositionMode = sourceResult.positionMode;
          candidate.stats.siskanGapMode = sourceResult.gapsAsState ? "fifth-state" : "strip";
          candidate.stats.siskanScanPermutations = region.scanPermutations;
          candidate.stats.siskanPValuePermutations = region.pValuePermutations;
          candidate.stats.siskanPattern = region.pattern;
          candidate.stats.siskanScoreFamily = region.scoreFamily;
          candidate.stats.siskanBaselineTopology = sourceResult.baselineTopology;
          candidate.stats.siskanInferredTopology = region.inferredTopology;
          candidate.stats.siskanSourceRoutine = sourceResult.sourceRoutine;
          candidate.methodSignals.push({
            method: "SiScan",
            ...mapped,
            statistic: region.z,
            locator: `RDP5 Sister-Scanning ${region.scoreFamily} ${region.pattern} topology run`,
            sourceRoutine: sourceResult.sourceRoutine,
            outgroup: sourceResult.outgroupIndex === null ? null : scanMappings[sourceResult.outgroupIndex].originIndex,
            outgroupMode: sourceResult.outgroupMode,
            outgroupSampled: sourceResult.outgroupSampled,
            permutations: region.pValuePermutations,
            scanPermutations: region.scanPermutations,
            pattern: region.pattern,
            scoreFamily: region.scoreFamily,
            baselineTopology: sourceResult.baselineTopology,
            inferredTopology: region.inferredTopology,
            profile: sourceSiScanProfile(
              sourceResult,
              candidate.rotation,
              nSites,
              unique.length > 1000 ? 48 : unique.length > 250 ? 96 : 192,
            ),
          });
        }
      }
    }
    if (options.methods.includes("3Seq")) mapSignal("3Seq", methodOutput[23], methodOutput[24], methodOutput[11], "maximum HGRW descent");
    if (options.polishBreakpoints) {
      let polishedStart = methodOutput[17];
      let polishedEnd = methodOutput[18];
      candidate.breakpointModel = {
        method: "local-chi-square",
        informativeSites: candidate.informative,
      };
      const burt = fitBurtTriplet(
        candidate.rotation === 0 ? scanEncoded : rotated,
        nSites,
        candidate.analysisRecombinant,
        candidate.analysisMajorParent,
        candidate.analysisMinorParent,
        candidate.rawStart,
        candidate.rawEnd,
        {
          sourceParity: options.burtMode !== "manual-step-up",
          randomStarts: options.burtRandomStarts,
          maxIterations: options.burtMaxIterations,
          maxStates: options.burtMaxStates,
          exhaustiveModels: options.burtExhaustiveModels,
          posteriorThreshold: options.burtPosteriorThreshold,
          circular: options.circular,
          seed: options.burtMode !== "manual-step-up"
            ? (options.randomSeed ?? 0x5a17c0de) >>> 0
            : (
                (options.randomSeed ?? 0x5a17c0de)
                ^ Math.imul(candidate.analysisRecombinant + 1, 0x9e3779b1)
                ^ Math.imul(candidate.analysisMajorParent + 1, 0x85ebca6b)
                ^ Math.imul(candidate.analysisMinorParent + 1, 0xc2b2ae35)
              ) >>> 0,
        },
      );
      if (burt) {
        polishedStart = burt.start;
        polishedEnd = burt.end;
        candidate.breakpointModel = mapBreakpointModel(burt.model, candidate.rotation, nSites);
        candidate.confidenceStart = mapConfidenceInterval(burt.confidenceStart, candidate.rotation, nSites);
        candidate.confidenceEnd = mapConfidenceInterval(burt.confidenceEnd, candidate.rotation, nSites);
      }
      const polishedLength = polishedStart <= polishedEnd ? polishedEnd - polishedStart : nSites + polishedEnd - polishedStart;
      if (polishedStart >= 0 && polishedEnd <= nSites && polishedLength >= 12) {
        Object.assign(candidate, mapInterval(polishedStart, polishedEnd, candidate.rotation, nSites));
      }
    }
    candidate.methodSignals = (candidate.methodSignals ?? []).filter((signal) => isCoLocatedSignal(candidate, signal, nSites));
    candidate.diagnostics = candidateDiagnostics(candidate, encoded, nSites, diagnosticProfile);
    if (candidateIndex % 16 === 0 || candidateIndex === unique.length - 1) {
      postMessage({
        type: "progress",
        jobId,
        progress: 0.85 + 0.15 * (candidateIndex + 1) / Math.max(1, unique.length),
        phase: "Calibrating method evidence",
      });
    }
  }
  let statisticsMs = performance.now() - statisticsStarted;
  let exactThreeSeqBudget = options.methods.includes("3Seq") ? 20_000_000 : 0;
  let events = unique.map((candidate, index) => {
    const threeSeqOperations = (candidate.stats.threeSeqMajorSites + 1)
      * (candidate.stats.threeSeqMinorSites + 1)
      * Math.max(1, candidate.stats.threeSeqDescent);
    const exactOperations = Math.min(4_000_000, exactThreeSeqBudget);
    const evidence = methodEvidence(candidate, candidate.stats, {
      ...options,
      threeSeqMaxOperations: exactOperations,
    }, Math.max(1, comparisons), nSites);
    if (options.methods.includes("3Seq") && threeSeqOperations <= exactOperations) exactThreeSeqBudget -= threeSeqOperations;
    const supportedCount = evidence.filter((item) => item.supported).length;
    const confidence = Math.max(4, Math.floor(options.window / Math.max(6, 2 + Math.sqrt(candidate.informative))));
    return {
      id: `wasm-${jobId}-${index + 1}`,
      recombinant: candidate.recombinant,
      majorParent: candidate.majorParent,
      minorParent: candidate.minorParent,
      start: candidate.start,
      end: candidate.end,
      wraps: candidate.wraps,
      confidenceStart: candidate.confidenceStart ?? [Math.max(0, candidate.start - confidence), Math.min(nSites, candidate.start + confidence)],
      confidenceEnd: candidate.confidenceEnd ?? [Math.max(0, candidate.end - confidence), Math.min(nSites, candidate.end + confidence)],
      breakpointModel: candidate.breakpointModel,
      methodSignals: primaryMethodSignals(candidate.methodSignals),
      evidence,
      chiSquare: candidate.chiSquare,
      informativeSites: candidate.informative,
      decision: "unreviewed",
      warnings: addWarnings(candidate, sequences, nSites, options),
      note: [
        candidate.wraps ? "Origin-spanning circular tract." : "",
        candidate.componentProvenance
          ? `Detected after RDP5 signal disassembly of ${candidate.componentProvenance.appliedEventIds.length} accepted event${candidate.componentProvenance.appliedEventIds.length === 1 ? "" : "s"}.`
          : "",
        candidate.structuralUncertainty
          ? `Continuous piece ${candidate.structuralUncertainty.piece}/${candidate.structuralUncertainty.pieces}; gap-adjacent breakpoint${candidate.structuralUncertainty.uncertainStart && candidate.structuralUncertainty.uncertainEnd ? "s are" : " is"} uncertain.`
          : "",
        candidate.alternatives.length
          ? `${candidate.alternatives.length} alternative minor-parent candidate${candidate.alternatives.length === 1 ? "" : "s"} grouped with this signal.`
          : "",
      ].filter(Boolean).join(" "),
      source: "wasm",
      supportedCount,
      diagnostics: candidate.diagnostics,
      groupId: null,
      alternativeParents: candidate.alternatives.filter((parent) => parent !== candidate.recombinant && parent !== candidate.majorParent && parent !== candidate.minorParent),
      hypothesisTests: Math.max(1, comparisons),
      history: [{
        id: `history-${jobId}-${index + 1}-1`,
        timestamp: new Date().toISOString(),
        action: "Detected by scan",
        summary: `${supportedCount} method families supported the initial hypothesis${candidate.componentProvenance ? " on the disassembled component alignment" : ""}.`,
      }],
      evidenceStale: false,
      componentProvenance: candidate.componentProvenance,
      structuralUncertainty: candidate.structuralUncertainty,
      analysisRecombinant: candidate.analysisRecombinant,
      analysisMajorParent: candidate.analysisMajorParent,
      analysisMinorParent: candidate.analysisMinorParent,
    };
  });

  if (options.correction === "holm") {
    for (const method of options.methods) {
      const methodRows = events
        .map((event) => ({ event, evidence: event.evidence.find((item) => item.method === method) }))
        .filter((row) => row.evidence)
        .sort((left, right) => left.evidence.pValue - right.evidence.pValue);
      let previous = 0;
      methodRows.forEach((row, rank) => {
        const corrected = Math.max(previous, Math.min(1, row.evidence.pValue * (methodRows.length - rank)));
        row.evidence.correctedP = corrected;
        row.evidence.supported = corrected <= options.alpha;
        row.evidence.correctionScope = `Holm family across ${methodRows.length.toLocaleString()} retained candidate hypotheses`;
        previous = corrected;
      });
    }
    events.forEach((event) => { event.supportedCount = event.evidence.filter((item) => item.supported).length; });
  }

  events = events
    .filter((event) => event.supportedCount >= options.minMethods)
    .sort((left, right) => {
      const leftP = Math.min(...left.evidence.map((item) => item.correctedP));
      const rightP = Math.min(...right.evidence.map((item) => item.correctedP));
      return leftP - rightP;
    });

  const clusteringStarted = performance.now();
  const analysisEvents = events.map((event) => ({
    ...event,
    displayRoles: {
      recombinant: event.recombinant,
      majorParent: event.majorParent,
      minorParent: event.minorParent,
    },
    recombinant: event.analysisRecombinant,
    majorParent: event.analysisMajorParent,
    minorParent: event.analysisMinorParent,
  }));
  const rawClustered = inferAncestralEventClusters(
    analysisEvents,
    scanEncoded,
    scanSequenceCount,
    nSites,
    { ...options, roleDmaxEvaluator },
    sourceMaximumDistance,
  );
  const originIndex = (index) => scanMappings[index]?.originIndex ?? index;
  const mapEvidenceRows = (rows = []) => {
    const byOrigin = new Map();
    for (const row of rows) {
      const mapped = { ...row, sequence: originIndex(row.sequence) };
      const previous = byOrigin.get(mapped.sequence);
      if (!previous || mapped.sets > previous.sets || (mapped.sets === previous.sets && mapped.topologyMargin > previous.topologyMargin)) {
        byOrigin.set(mapped.sequence, mapped);
      }
    }
    return [...byOrigin.values()].sort((left, right) => left.sequence - right.sequence);
  };
  const mapDirection = (direction) => direction && typeof direction === "object"
    ? { ...direction, sequence: originIndex(direction.sequence) }
    : direction;
  const mapClusterEvent = (event) => {
    const displayRoles = event.displayRoles ?? {
      recombinant: originIndex(event.recombinant),
      majorParent: originIndex(event.majorParent),
      minorParent: originIndex(event.minorParent),
    };
    const mappedSets = event.coRecombinantSets?.map((set) => ({
      ...set,
      presumedRecombinant: originIndex(set.presumedRecombinant),
      parents: [...new Set(set.parents.map(originIndex))],
      sequenceMembers: [...new Set(set.sequenceMembers.map(originIndex))].sort((left, right) => left - right),
      evidence: mapEvidenceRows(set.evidence),
    }));
    const mappedCluster = event.ancestralCluster
      ? {
          ...event.ancestralCluster,
          sequenceMembers: [...new Set(event.ancestralCluster.sequenceMembers.map(originIndex))].sort((left, right) => left - right),
          pairwise: event.ancestralCluster.pairwise.map((pair) => ({
            ...pair,
            leftToRight: mapDirection(pair.leftToRight),
            rightToLeft: mapDirection(pair.rightToLeft),
          })),
        }
      : undefined;
    const mappedIdentification = event.recombinantIdentification
      ? {
          ...event.recombinantIdentification,
          candidates: event.recombinantIdentification.candidates.map(originIndex),
          recommended: originIndex(event.recombinantIdentification.recommended),
          recommendedMajorParent: originIndex(event.recombinantIdentification.recommendedMajorParent),
          recommendedMinorParent: originIndex(event.recombinantIdentification.recommendedMinorParent),
          orientations: event.recombinantIdentification.orientations.map((orientation) => ({
            ...orientation,
            recombinant: originIndex(orientation.recombinant),
            majorParent: originIndex(orientation.majorParent),
            minorParent: originIndex(orientation.minorParent),
          })),
        }
      : undefined;
    const rest = { ...event };
    delete rest.displayRoles;
    delete rest.analysisRecombinant;
    delete rest.analysisMajorParent;
    delete rest.analysisMinorParent;
    return {
      ...rest,
      ...displayRoles,
      coRecombinantSets: mappedSets,
      ancestralCluster: mappedCluster,
      recombinantIdentification: mappedIdentification,
    };
  };
  const clustered = {
    ...rawClustered,
    events: rawClustered.events.map(mapClusterEvent),
    clusters: rawClustered.clusters.map((cluster) => ({
      ...cluster,
      sequenceMembers: [...new Set(cluster.sequenceMembers.map(originIndex))].sort((left, right) => left - right),
    })),
  };
  events = clustered.events;
  const clusteringMs = performance.now() - clusteringStarted;
  statisticsMs += clusteringMs;

  const eventsByRecombinant = new Map();
  events.forEach((event) => eventsByRecombinant.set(event.recombinant, [...(eventsByRecombinant.get(event.recombinant) ?? []), event]));
  for (const group of eventsByRecombinant.values()) {
    if (group.length < 2) continue;
    for (let left = 0; left < group.length; left += 1) {
      for (let right = left + 1; right < group.length; right += 1) {
        const shared = overlap(group[left], group[right], nSites);
        if (shared <= 0.05) continue;
        const warning = "Multiple retained signals overlap in this recombinant; treat them as nested/overprinted hypotheses, not as co-recombinant descendants.";
        if (!group[left].warnings.includes(warning)) group[left].warnings.push(warning);
        if (!group[right].warnings.includes(warning)) group[right].warnings.push(warning);
      }
    }
  }

  postMessage({
    type: "result",
    jobId,
    events,
    distance: Array.from(distance),
    comparisons,
    elapsedMs: performance.now() - started,
    timing: { distanceMs, scanMs, statisticsMs, diagnosticsMs, clusteringMs },
    ancestralClusters: clustered.clusters,
    diagnostics: diagnosticProfile.summary,
    disassembly: {
      appliedEvents: disassembly.appliedEventIds.length,
      components: disassembly.componentCount,
      erasedCanonicalBases: disassembly.erasedCanonicalBases,
    },
    rdpSignalTruncations: truncatedRdpSignals,
    matrixCount,
    parentSamples: exactDistanceMatrix ? nSites : Math.min(parentSamples, nSites),
    matrixMode: exactDistanceMatrix && exactDisplayMatrix
      ? "exact"
      : `${exactDisplayMatrix ? "exact display" : "24-sequence display"} + ${exactDistanceMatrix ? "exact component parent search" : "sampled/stratified component parent search"}`,
    tripletMode: options.exhaustive ? "all-concrete-triplets" : "approximate-parent-shortlist",
    // Every method invocation above receives three explicit indexes into the
    // working alignment. No alignment consensus or rest-of-alignment proxy is
    // ever substituted for a triplet member. SiScan may separately select a
    // real fourth outgroup, as required by that method.
    concreteTripletInputs: true,
    engine: [
      "WebAssembly",
      options.exhaustive ? "all concrete sequence triplets" : "approximate parent-shortlist triplets",
      exactDistanceMatrix ? "packed distance" : "sampled parent search",
      options.circular ? "dual-origin circular scan" : "linear scan",
      options.polishBreakpoints ? (options.burtMode === "manual-step-up" ? "BURT 2–20-state step-up" : "RDP5-source BURT") : "raw breakpoints",
      options.ancestralClustering === false ? "event clustering off" : `${clustered.clusters.length} ancestral clusters`,
      disassembly.componentCount > 0 ? `${disassembly.componentCount} extracted analysis components` : "intact alignment",
    ].join(" · "),
  });
}

async function recalculate(message) {
  const started = performance.now();
  const instance = await loadWasm();
  const { alignment, options, event, jobId } = message;
  const sequences = alignment.sequences;
  const nSeq = sequences.length;
  const nSites = alignment.length;
  const encoded = encodeSequences(sequences, nSites);
  const profile = alignmentDiagnostics(encoded, nSeq, nSites, options.window, options.randomSeed);
  const disassembly = buildDisassembledAlignment(
    encoded,
    sequences,
    nSites,
    Array.isArray(message.disassemblyEvents) ? message.disassemblyEvents : [],
  );
  let analysisRecombinant = event.recombinant;
  let analysisMajorParent = event.majorParent;
  let analysisMinorParent = event.minorParent;
  if (event.componentProvenance) {
    analysisRecombinant = findComponentIndex(disassembly, event.componentProvenance.recombinant);
    analysisMajorParent = findComponentIndex(disassembly, event.componentProvenance.majorParent);
    analysisMinorParent = findComponentIndex(disassembly, event.componentProvenance.minorParent);
    if ([analysisRecombinant, analysisMajorParent, analysisMinorParent].some((index) => index < 0)) {
      throw new Error("Could not rebuild this event's signal-disassembly lineage. Restore its accepted predecessor events, then recalculate again.");
    }
  }
  let working = disassembly.encoded;
  const workingSequenceCount = disassembly.analysisSequences.length;
  let rawStart = event.start;
  let rawEnd = event.end;
  let rotation = 0;
  if (event.wraps && event.start > event.end) {
    const backgroundLength = event.start - event.end;
    rotation = (event.end + Math.floor(backgroundLength / 2)) % nSites;
    working = rotateSequences(disassembly.encoded, workingSequenceCount, nSites, rotation);
    rawStart = (event.start - rotation + nSites) % nSites;
    rawEnd = (event.end - rotation + nSites) % nSites;
  }
  if (rawEnd <= rawStart) throw new Error("The edited event does not define a valid tract.");

  let mappedInterval = { start: event.start, end: event.end, wraps: event.wraps };
  let confidenceStart = [event.start, event.start];
  let confidenceEnd = [event.end, event.end];
  let breakpointModel = event.breakpointModel?.method === "manual"
    ? event.breakpointModel
    : { method: "local-chi-square", informativeSites: event.informativeSites };
  if (options.polishBreakpoints && event.breakpointModel?.method !== "manual") {
    const burt = fitBurtTriplet(
      working,
      nSites,
      analysisRecombinant,
      analysisMajorParent,
      analysisMinorParent,
      rawStart,
      rawEnd,
      {
        sourceParity: options.burtMode !== "manual-step-up",
        randomStarts: options.burtRandomStarts,
        maxIterations: options.burtMaxIterations,
        maxStates: options.burtMaxStates,
        exhaustiveModels: options.burtExhaustiveModels,
        posteriorThreshold: options.burtPosteriorThreshold,
        circular: options.circular,
        seed: options.burtMode !== "manual-step-up"
          ? (options.randomSeed ?? 0x5a17c0de) >>> 0
          : (
              (options.randomSeed ?? 0x5a17c0de)
              ^ Math.imul(event.recombinant + 1, 0x9e3779b1)
              ^ Math.imul(event.majorParent + 1, 0x85ebca6b)
              ^ Math.imul(event.minorParent + 1, 0xc2b2ae35)
            ) >>> 0,
      },
    );
    const burtLength = burt ? (burt.start <= burt.end ? burt.end - burt.start : nSites + burt.end - burt.start) : 0;
    if (burt && burtLength >= 4) {
      rawStart = burt.start;
      rawEnd = burt.end;
      mappedInterval = mapInterval(rawStart, rawEnd, rotation, nSites);
      confidenceStart = mapConfidenceInterval(burt.confidenceStart, rotation, nSites);
      confidenceEnd = mapConfidenceInterval(burt.confidenceEnd, rotation, nSites);
      breakpointModel = mapBreakpointModel(burt.model, rotation, nSites);
    }
  }
  const seqPtr = 65536;
  const prefixAPtr = align(seqPtr + working.byteLength, 16);
  const prefixBPtr = prefixAPtr + (nSites + 1) * 4;
  const outPtr = align(prefixBPtr + (nSites + 1) * 4, 16);
  const statsPtr = align(outPtr + 96, 16);
  const workingPacked = packSequences(working, workingSequenceCount, nSites);
  const packedPtr = align(statsPtr + 160, 16);
  const validityPtr = packedPtr + workingPacked.packed.byteLength;
  const sourceSiScanExactDistance = options.methods.includes("SiScan")
    && workingSequenceCount <= 512
    && workingSequenceCount * workingSequenceCount * workingPacked.wordsPerSequence <= 50_000_000;
  const sourceSiScanDistancePtr = align(validityPtr + workingPacked.validity.byteLength, 16);
  const sourceSiScanDistanceBytes = sourceSiScanExactDistance ? workingSequenceCount * workingSequenceCount * 4 : 0;
  const roleCohortCapacity = Math.max(4, Math.min(34, Math.trunc(options.roleQuartetTaxaLimit ?? 30)));
  const roleCohortPtr = align(sourceSiScanDistancePtr + sourceSiScanDistanceBytes, 16);
  const roleTractMaskPtr = roleCohortPtr + roleCohortCapacity * 4;
  const roleBackgroundMaskPtr = roleTractMaskPtr + workingPacked.wordsPerSequence * 4;
  const roleDmaxOutPtr = align(roleBackgroundMaskPtr + workingPacked.wordsPerSequence * 4, 8);
  const requiredPages = Math.ceil((roleDmaxOutPtr + 40) / 65536);
  const currentPages = instance.exports.memory.buffer.byteLength / 65536;
  if (requiredPages > currentPages) instance.exports.memory.grow(requiredPages - currentPages);
  const memory = instance.exports.memory;
  new Uint8Array(memory.buffer, seqPtr, working.byteLength).set(working);
  new Uint32Array(memory.buffer, packedPtr, workingPacked.packed.length).set(workingPacked.packed);
  new Uint32Array(memory.buffer, validityPtr, workingPacked.validity.length).set(workingPacked.validity);
  let sourceSiScanDistance = null;
  let sourceSiScanTreeDistance = null;
  if (sourceSiScanExactDistance) {
    instance.exports.distance_matrix_packed(
      packedPtr,
      validityPtr,
      workingSequenceCount,
      workingPacked.wordsPerSequence,
      sourceSiScanDistancePtr,
    );
    sourceSiScanDistance = new Float32Array(
      memory.buffer,
      sourceSiScanDistancePtr,
      workingSequenceCount * workingSequenceCount,
    ).slice();
    if (workingSequenceCount <= 128) {
      sourceSiScanTreeDistance = buildNeighborJoiningPathMatrix(sourceSiScanDistance, workingSequenceCount);
    }
  }
  const roleDmaxEvaluator = ({ event: roleEvent, candidates: roleCandidates, cohort }) => {
    if (typeof instance.exports.dmax_visrd_packed !== "function" || cohort.length > roleCohortCapacity) return null;
    const tractMask = new Uint32Array(workingPacked.wordsPerSequence);
    const backgroundMask = new Uint32Array(workingPacked.wordsPerSequence);
    for (let site = 0; site < nSites; site += 1) {
      const inside = roleEvent.wraps && roleEvent.start > roleEvent.end
        ? site >= roleEvent.start || site < roleEvent.end
        : site >= roleEvent.start && site < roleEvent.end;
      const word = site >>> 4;
      const lane = 1 << ((site & 15) << 1);
      if (inside) tractMask[word] |= lane;
      else backgroundMask[word] |= lane;
    }
    new Int32Array(memory.buffer, roleCohortPtr, cohort.length).set(cohort);
    new Uint32Array(memory.buffer, roleTractMaskPtr, workingPacked.wordsPerSequence).set(tractMask);
    new Uint32Array(memory.buffer, roleBackgroundMaskPtr, workingPacked.wordsPerSequence).set(backgroundMask);
    instance.exports.dmax_visrd_packed(
      packedPtr,
      validityPtr,
      workingPacked.wordsPerSequence,
      roleCohortPtr,
      cohort.length,
      roleCandidates[0],
      roleCandidates[1],
      roleCandidates[2],
      roleTractMaskPtr,
      roleBackgroundMaskPtr,
      roleDmaxOutPtr,
    );
    return {
      values: [...new Float64Array(memory.buffer, roleDmaxOutPtr, 3)],
      quartetCounts: [...new Int32Array(memory.buffer, roleDmaxOutPtr + 24, 3)],
      cohortSize: cohort.length,
      sourceRoutine: "CalcMaxD + CMaxD2P3",
      wasmAccelerated: true,
    };
  };

  instance.exports.triplet_counts(
    seqPtr,
    nSites,
    analysisRecombinant,
    analysisMajorParent,
    analysisMinorParent,
    rawStart,
    rawEnd,
    outPtr,
  );
  const counts = new Int32Array(instance.exports.memory.buffer, outPtr, 6).slice();
  instance.exports.method_stats(
    seqPtr,
    nSites,
    analysisRecombinant,
    analysisMajorParent,
    analysisMinorParent,
    rawStart,
    rawEnd,
    Math.max(20, options.window),
    Math.max(1, options.step),
    options.methods.includes("BootScan") ? Math.max(0, options.bootstrapReplicates ?? 100) : 0,
    (options.randomSeed ?? 0x5a17c0de) | 0,
    enabledMethodMask(options),
    prefixAPtr,
    prefixBPtr,
    statsPtr,
    options.geneconvGScale ?? 1,
  );
  const output = new Int32Array(instance.exports.memory.buffer, statsPtr, 35);
  const stats = {
    genconvRun: output[0],
    genconvEligible: output[1],
    genconvMatches: output[2],
    genconvStart: output[3],
    genconvEnd: output[4],
    bootscanConsistent: output[5],
    bootscanWindows: output[6],
    maxChi: output[7] / 1000,
    chimaera: output[8] / 1000,
    siskanScore: output[9],
    siskanSites: output[10],
    threeSeqDescent: output[11],
    threeSeqSites: output[12],
    maxChiBoundaries: [output[13] / 1000, output[14] / 1000],
    chimaeraBoundaries: [output[15] / 1000, output[16] / 1000],
    threeSeqMajorSites: output[19],
    threeSeqMinorSites: output[20],
    bootscanBootstrapConsistent: output[21],
    bootscanBootstrapReplicates: output[22],
    threeSeqStart: output[23],
    threeSeqEnd: output[24],
    maxChiStart: output[25],
    maxChiEnd: output[26],
    chimaeraStart: output[27],
    chimaeraEnd: output[28],
    bootscanStart: output[29],
    bootscanEnd: output[30],
    siskanStart: output[31],
    siskanEnd: output[32],
    bootscanRunWindows: output[33],
    siskanRunWindows: output[34],
  };
  const candidate = {
    recombinant: event.recombinant,
    majorParent: event.majorParent,
    minorParent: event.minorParent,
    start: mappedInterval.start,
    end: mappedInterval.end,
    wraps: mappedInterval.wraps,
    rawStart,
    rawEnd,
    rotation,
    sequencePtr: seqPtr,
    analysisRecombinant,
    analysisMajorParent,
    analysisMinorParent,
    componentProvenance: event.componentProvenance
      ? candidateComponentProvenance(disassembly, analysisRecombinant, analysisMajorParent, analysisMinorParent)
      : undefined,
    insideMinor: counts[1],
    insideMajor: counts[2],
    outsideMajor: counts[3],
    outsideMinor: counts[4],
    informative: counts[0],
    chiSquare: counts[5] / 1000,
    stats,
    diagnostics: null,
  };
  let recalculatedSiScanSignal = null;
  if (options.methods.includes("SiScan")) {
    const preliminary = output[32] > output[31]
      ? { method: "SiScan", ...mapInterval(output[31], output[32], rotation, nSites) }
      : null;
    const sourceRequested = (event.methodSignals ?? []).some((signal) => signal.method === "SiScan")
      || isCoLocatedSignal(candidate, preliminary, nSites);
    stats.siskanScore = 0;
    stats.siskanSites = 0;
    stats.siskanSourceP = 1;
    stats.siskanSourceZ = 0;
    if (sourceRequested) {
      const triplet = [analysisRecombinant, analysisMajorParent, analysisMinorParent];
      const tripletOrigins = new Set(triplet.map((index) => disassembly.mappings[index].originIndex));
      const candidatePool = Array.from({ length: workingSequenceCount }, (_, index) => index)
        .filter((index) => !tripletOrigins.has(disassembly.mappings[index].originIndex));
      const manualOutgroup = options.siskanOutgroupMode === "manual"
        ? Array.from({ length: workingSequenceCount }, (_, index) => index).find((index) => (
            disassembly.mappings[index].originIndex === options.siskanOutgroupSequence
            && !tripletOrigins.has(disassembly.mappings[index].originIndex)
          ))
        : undefined;
      const sourceResult = runSourceSiScan(
        working,
        nSites,
        workingSequenceCount,
        triplet,
        {
          window: Math.max(20, options.window),
          step: Math.max(1, options.step),
          scanPermutations: options.siskanScanPermutations ?? 100,
          pValuePermutations: options.siskanPValuePermutations ?? 1000,
          seed: (options.randomSeed ?? 0x5a17c0de) >>> 0,
          outgroupMode: options.siskanOutgroupMode === "manual" && manualOutgroup === undefined
            ? "nearest"
            : options.siskanOutgroupMode ?? "nearest",
          outgroupIndex: manualOutgroup,
          positionMode: options.siskanPositionMode ?? "triplet-variable",
          gapsAsState: options.siskanGapMode === "fifth-state",
          candidatePool,
          distanceMatrix: sourceSiScanDistance,
          treeDistanceMatrix: sourceSiScanTreeDistance,
          randomization: buildSourceSiScanRandomization(
            nSites,
            Math.max(options.siskanScanPermutations ?? 100, options.siskanPValuePermutations ?? 1000),
            options.randomSeed ?? 0x5a17c0de,
          ),
        },
      );
      const selected = selectCoLocatedSiScanRegion(sourceResult, candidate, rotation, nSites);
      if (selected) {
        const { region, mapped } = selected;
        Object.assign(stats, {
          siskanScore: region.z,
          siskanSites: region.end - region.start,
          siskanStart: region.start,
          siskanEnd: region.end,
          siskanSourceP: region.rawP,
          siskanSourceZ: region.z,
          siskanOutgroupIndex: sourceResult.outgroupIndex,
          siskanOutgroupMode: sourceResult.outgroupMode,
          siskanOutgroupSampled: sourceResult.outgroupSampled,
          siskanOutgroupSourcePath: sourceResult.outgroupSourcePath,
          siskanPositionMode: sourceResult.positionMode,
          siskanGapMode: sourceResult.gapsAsState ? "fifth-state" : "strip",
          siskanScanPermutations: region.scanPermutations,
          siskanPValuePermutations: region.pValuePermutations,
          siskanPattern: region.pattern,
          siskanScoreFamily: region.scoreFamily,
          siskanBaselineTopology: sourceResult.baselineTopology,
          siskanInferredTopology: region.inferredTopology,
          siskanSourceRoutine: sourceResult.sourceRoutine,
        });
        recalculatedSiScanSignal = {
          method: "SiScan",
          ...mapped,
          statistic: region.z,
          locator: `RDP5 Sister-Scanning ${region.scoreFamily} ${region.pattern} topology run`,
          sourceRoutine: sourceResult.sourceRoutine,
          outgroup: sourceResult.outgroupIndex === null ? null : disassembly.mappings[sourceResult.outgroupIndex].originIndex,
          outgroupMode: sourceResult.outgroupMode,
          outgroupSampled: sourceResult.outgroupSampled,
          permutations: region.pValuePermutations,
          scanPermutations: region.scanPermutations,
          pattern: region.pattern,
          scoreFamily: region.scoreFamily,
          baselineTopology: sourceResult.baselineTopology,
          inferredTopology: region.inferredTopology,
          profile: sourceSiScanProfile(sourceResult, rotation, nSites, 192),
        };
      }
    }
  }
  candidate.diagnostics = candidateDiagnostics(candidate, encoded, nSites, profile);
  const recalculatedMethodSignals = [
    ...(options.methods.includes("RDP") ? [{ method: "RDP", ...mappedInterval, statistic: candidate.chiSquare, locator: "edited hypothesis recalculation" }] : []),
    ...(options.methods.includes("GENECONV") && output[4] > output[3] ? [{ method: "GENECONV", ...mapInterval(output[3], output[4], rotation, nSites), statistic: output[0], locator: "maximum concordant fragment" }] : []),
    ...(options.methods.includes("BootScan") && output[30] > output[29] ? [{ method: "BootScan", ...mapInterval(output[29], output[30], rotation, nSites), statistic: output[21] / Math.max(1, output[22]), locator: "minor-topology window run" }] : []),
    ...(options.methods.includes("MaxChi") && output[26] > output[25] ? [{ method: "MaxChi", ...mapInterval(output[25], output[26], rotation, nSites), statistic: output[7] / 1000, locator: "independent χ² peak pair" }] : []),
    ...(options.methods.includes("Chimaera") && output[28] > output[27] ? [{ method: "Chimaera", ...mapInterval(output[27], output[28], rotation, nSites), statistic: output[8] / 1000, locator: "independent binary χ² peak pair" }] : []),
    ...(recalculatedSiScanSignal ? [recalculatedSiScanSignal] : []),
    ...(options.methods.includes("3Seq") && output[24] > output[23] ? [{ method: "3Seq", ...mapInterval(output[23], output[24], rotation, nSites), statistic: output[11], locator: "maximum HGRW descent" }] : []),
  ].filter((signal) => isCoLocatedSignal(candidate, signal, nSites));
  candidate.methodSignals = primaryMethodSignals(recalculatedMethodSignals);
  const familyComparisons = Math.max(1, Math.trunc(message.comparisons ?? event.hypothesisTests ?? 1));
  // A single edited hypothesis cannot be inserted into the original Holm rank
  // ordering without rescanning the whole family. Use conservative Bonferroni
  // for this local recalculation instead of silently treating it as one test.
  const recalculationOptions = options.correction === "holm"
    ? { ...options, correction: "bonferroni" }
    : options;
  const evidence = methodEvidence(candidate, stats, recalculationOptions, familyComparisons, nSites);
  const recalculationScope = candidate.componentProvenance
    ? ` The exact signal-disassembly lineage was rebuilt from ${candidate.componentProvenance.appliedEventIds.length} accepted predecessor event${candidate.componentProvenance.appliedEventIds.length === 1 ? "" : "s"}.`
    : "";
  const recalculationNote = options.correction === "holm"
    ? `Edited-event recalculation uses conservative Bonferroni across ${familyComparisons.toLocaleString()} original scan triplets; rerun the full scan for exact Holm ranks and ancestral-event clustering.`
    : `Multiplicity correction retained the original scan scope of ${familyComparisons.toLocaleString()} triplets; rerun the full scan to refresh ancestral-event clustering.`;
  const rawIdentification = identifyRecombinantRoles({
    recombinant: analysisRecombinant,
    majorParent: analysisMajorParent,
    minorParent: analysisMinorParent,
    start: rawStart,
    end: rawEnd,
    wraps: false,
  }, working, workingSequenceCount, nSites, null, { ...options, roleDmaxEvaluator });
  const mapComponentOrigin = (index) => disassembly.mappings[index]?.originIndex ?? index;
  const recombinantIdentification = rawIdentification
    ? {
        ...rawIdentification,
        candidates: rawIdentification.candidates.map(mapComponentOrigin),
        recommended: mapComponentOrigin(rawIdentification.recommended),
        recommendedMajorParent: mapComponentOrigin(rawIdentification.recommendedMajorParent),
        recommendedMinorParent: mapComponentOrigin(rawIdentification.recommendedMinorParent),
        orientations: rawIdentification.orientations.map((orientation) => ({
          ...orientation,
          recombinant: mapComponentOrigin(orientation.recombinant),
          majorParent: mapComponentOrigin(orientation.majorParent),
          minorParent: mapComponentOrigin(orientation.minorParent),
        })),
      }
    : undefined;
  const roleWarnings = recombinantIdentification?.ambiguous
    ? [`RDP5 source profile consensus is role-ambiguous (${Math.round(recombinantIdentification.confidence * 100)}% relative support); inspect all three recombinant polarities.`]
    : recombinantIdentification && recombinantIdentification.recommended !== event.recombinant
      ? ["RDP5 source profile consensus challenges the current recombinant assignment; inspect and, if appropriate, apply the highest-scoring polarity."]
      : [];
  postMessage({
    type: "recalculated",
    jobId,
    patch: {
      evidence,
      chiSquare: candidate.chiSquare,
      informativeSites: candidate.informative,
      warnings: [...addWarnings(candidate, sequences, nSites, options), ...roleWarnings, ...(options.correction === "holm" ? [recalculationNote] : [])],
      diagnostics: candidate.diagnostics,
      evidenceStale: false,
      start: mappedInterval.start,
      end: mappedInterval.end,
      wraps: mappedInterval.wraps,
      confidenceStart,
      confidenceEnd,
      breakpointModel,
      methodSignals: candidate.methodSignals,
      ancestralCluster: undefined,
      componentProvenance: candidate.componentProvenance,
      structuralUncertainty: event.structuralUncertainty,
      recombinantIdentification,
      hypothesisTests: familyComparisons,
      recalculationNote: `${recalculationNote}${recalculationScope}`,
    },
    diagnostics: profile.summary,
    elapsedMs: performance.now() - started,
  });
}

self.onmessage = (event) => {
  if (event.data?.type !== "analyze" && event.data?.type !== "recalculate") return;
  const operation = event.data.type === "recalculate" ? recalculate : analyze;
  operation(event.data).catch((error) => {
    postMessage({
      type: "error",
      jobId: event.data.jobId,
      message: error instanceof Error ? error.message : String(error),
    });
  });
};
