import { methodEvidence, rdp5SourceProbability, threeSeqSiegmundP, threeSeqSourceP } from "./rdp-statistics.js";
import { fitBurtTriplet } from "./rdp-burt.js";
import { inferAncestralEventClusters } from "./rdp-clustering.js";
import { buildDisassembledAlignment, candidateComponentProvenance, findComponentIndex, sourceThreeSeqSubPValExcursion, splitCandidateAtStructuralGaps } from "./rdp-disassembly.js";
import { identifyRecombinantRoles } from "./rdp-recombinant-identification.js";
import { buildNeighborJoiningPathMatrix } from "./rdp-bootstrap-tree.js";
import { buildSourceSiScanRandomization, runSourceSiScan, sourceSiScanRoles } from "./rdp-siscan.js";
import { sourcePhiTest } from "./rdp-phi.js";

let wasmPromise;

const BASES = { A: 0, C: 1, G: 2, T: 3, U: 3, "-": 5 };

function enabledMethodMask(options) {
  // Production discovery is source-only. BootScan has its own whole-cohort
  // batch below; 3Seq has its own fused concrete-triplet kernel.
  const bits = { GENECONV: 2, MaxChi: 4, Chimaera: 8, SiScan: 16 };
  let mask = 0;
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

function rawIntervalLength(start, end, length) {
  return end >= start ? end - start : length - start + end;
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

const SOURCE_CHI_ROW_INTS = 16;
const SOURCE_CHI_PEAK_INTS = 6;
const SOURCE_GENECONV_ROW_INTS = 16;
const SOURCE_BOOTSCAN_ROW_INTS = 16;
const SOURCE_THREE_SEQ_ROW_INTS = 16;

function sourceChiRoutine(method) {
  return method === "MaxChi"
    ? "FindSubSeqMCPB2 → WinScoreCalcP → CalcChiVals4P3 → SmoothChiValsP → FindMChiP → GrowMChiWinP2"
    : "FindSubSeqDP3 → WinScoreCalc4P2 → CalcChiVals5P → SmoothChiVals3P → FindMChi3P → GrowMChiWinP";
}

function sourceChiSignal(row, rotation, length) {
  const method = row[0] === 3 ? "MaxChi" : "Chimaera";
  const mapped = mapInterval(row[3], row[4], rotation, length);
  return {
    method,
    ...mapped,
    statistic: row[6] / 1000,
    locator: method === "MaxChi"
      ? `RDP5 pair-equality track ${row[2] + 1} · paired source peak basin`
      : `RDP5 recombinant track ${row[2] + 1} · paired source peak basin`,
    sourceRoutine: sourceChiRoutine(method),
    sourceChi: {
      track: row[2],
      targetSlot: row[1] < 0 ? null : row[1],
      informativeSites: row[9],
      halfWindow: row[10],
      boundaryStatistics: [row[7] / 1000, row[8] / 1000],
      boundaryRanks: [row[11], row[12]],
      growthWidths: [row[13], row[14]],
      direction: row[15] < 0 ? -1 : 1,
    },
  };
}

function sourceGeneconvProbability(row) {
  const negativeLog = Math.max(0, row[12]) / 1_000_000;
  const probability = Math.exp(-negativeLog);
  return probability > 0 ? probability : Number.MIN_VALUE;
}

function sourceGeneconvSignal(row, rotation, length) {
  const mapped = mapInterval(row[4], row[5], rotation, length);
  return {
    method: "GENECONV",
    ...mapped,
    statistic: row[7],
    locator: `RDP5 six-track fragment queue · track ${row[0] + 1}`,
    sourceRoutine: "FindSubSeqGCAP6/7 → GetFragsP → GetMaxFragScoreP → CalcKMaxP → GCCalcPValP → GCGetHiPValP/DelPValsP",
    sourceGeneconv: {
      track: row[0],
      targetSlot: row[1],
      minorSlot: row[2],
      majorSlot: row[3],
      fragmentScore: row[7],
      informativeSites: row[8],
      matchingSites: row[9],
      mismatchSites: row[10],
      mismatchPenalty: row[11],
      rawP: sourceGeneconvProbability(row),
      startRank: row[13],
      endRank: row[14],
    },
  };
}

function sourceBootscanSignal(row, rotation, length) {
  const mapped = mapInterval(row[4], row[5], rotation, length);
  const rawP = rdp5SourceProbability({
    common: row[9],
    tractSites: row[11],
    mediumSites: row[9] + row[10],
    informativeSites: row[12],
    probabilitySites: row[12],
  }) ?? 1;
  return {
    method: "BootScan",
    ...mapped,
    statistic: row[6] / Math.max(1, row[7]),
    locator: `RDP5 RecScan distance topology ${row[3] + 1} · ${row[8]} supported windows`,
    sourceRoutine: "BSXoverR2 → SEQBOOT2 → FastBootDistIP → GetPltVal2 → FindBeginBS/FindEndBS → MakeScoresBS/ProbCalc",
    sourceBootscan: {
      topology: row[3],
      baselineTopology: row[13],
      bootstrapSupport: row[6] / Math.max(1, row[7]),
      bootstrapReplicates: row[7],
      runWindows: row[8],
      tractPairMatches: row[9],
      backgroundPairMatches: row[10],
      tractInformativeSites: row[11],
      informativeSites: row[12],
      rawP,
      window: row[14],
      step: row[15],
      relationshipMode: "distance",
    },
  };
}

function sourceBootscanRoles(row) {
  const sequences = [row[0], row[1], row[2]];
  const pairSlots = [[0, 1], [0, 2], [1, 2]];
  const tractPair = pairSlots[row[3]];
  const baselinePair = pairSlots[row[13]];
  const recombinantSlot = tractPair.find((slot) => baselinePair.includes(slot));
  if (recombinantSlot === undefined) return null;
  const majorSlot = baselinePair.find((slot) => slot !== recombinantSlot);
  const minorSlot = tractPair.find((slot) => slot !== recombinantSlot);
  if (majorSlot === undefined || minorSlot === undefined) return null;
  return {
    recombinant: sequences[recombinantSlot],
    majorParent: sequences[majorSlot],
    minorParent: sequences[minorSlot],
  };
}

function sourceThreeSeqSignal(row, rotation, length, exactOperations = 1_000_000) {
  const probability = threeSeqSourceP(row[6], row[7], row[8], exactOperations);
  const mapped = mapInterval(row[4], row[5], rotation, length);
  return {
    method: "3Seq",
    ...mapped,
    statistic: row[8],
    locator: `RDP5 3Seq target walk · cycle ${row[10] + 1} · ${row[3] > 0 ? "descent" : "ascent"}`,
    sourceRoutine: "FindSubSeqTS/FindSubSeqTS2 → CheckwrapC → GetTSPVal → Seq3PVals/Get3SeqPvalC",
    sourceThreeSeq: {
      target: row[0],
      majorParent: row[1],
      minorParent: row[2],
      direction: row[3] > 0 ? 1 : -1,
      upSteps: row[6],
      downSteps: row[7],
      descent: row[8],
      informativeSites: row[9],
      cycle: row[10],
      rawStart: row[4],
      rawEnd: row[5],
      rawP: probability.p,
      probabilityMode: probability.mode,
      sourceWrap: row[11] === 1,
      linearComplement: row[12] === 1,
    },
  };
}

function refineSourceThreeSeqSignalForPiece(
  signal,
  candidate,
  encoded,
  length,
  exactOperations = 1_000_000,
) {
  if (!signal?.sourceThreeSeq || !candidate?.structuralUncertainty) return signal;
  const source = signal.sourceThreeSeq;
  if (candidate.rawStart === source.rawStart && candidate.rawEnd === source.rawEnd) return signal;
  const piece = sourceThreeSeqSubPValExcursion(
    encoded,
    length,
    source.target,
    source.majorParent,
    source.minorParent,
    candidate.rawStart,
    candidate.rawEnd,
  );
  if (piece.informativeSites < 1 || piece.excursion < 1) return null;
  // SubPVal deliberately uses the original full-walk nM/nN counts and only
  // substitutes the continuously observed piece's height range for nK.
  const probability = threeSeqSourceP(
    source.upSteps,
    source.downSteps,
    piece.excursion,
    exactOperations,
  );
  return {
    ...signal,
    start: candidate.start,
    end: candidate.end,
    wraps: candidate.wraps === true,
    statistic: piece.excursion,
    locator: `${signal.locator} · CheckSplit3Seq piece ${candidate.structuralUncertainty.piece}/${candidate.structuralUncertainty.pieces}`,
    sourceRoutine: `${signal.sourceRoutine} → CheckSplit3Seq/SubPVal`,
    sourceThreeSeq: {
      ...source,
      fullDescent: source.fullDescent ?? source.descent,
      descent: piece.excursion,
      splitInformativeSites: piece.informativeSites,
      splitRefined: true,
      rawStart: candidate.rawStart,
      rawEnd: candidate.rawEnd,
      rawP: probability.p,
      probabilityMode: probability.mode,
    },
  };
}

// TSXOver evaluates both directions, then keeps the smaller raw probability
// (descent wins exact ties). A cheap source Siegmund screen prevents exact
// table-equivalent DPs from dominating large all-triplet analyses; plausible
// signals are then recalculated through the exact/source dispatcher.
function selectSourceThreeSeqSignals(rows, rotation, length, alpha, exactOperations = 1_000_000) {
  const grouped = new Map();
  for (const row of rows) grouped.set(row[0], [...(grouped.get(row[0]) ?? []), row]);
  const selected = [];
  const prefilter = Math.min(0.5, Math.max(0.1, alpha * 20));
  for (const targetRows of grouped.values()) {
    let approximateMinimum = 1;
    for (const row of targetRows) {
      const approximate = threeSeqSiegmundP(row[6], row[7], row[8]);
      if (approximate !== null) approximateMinimum = Math.min(approximateMinimum, approximate);
      else if ((row[6] + 1) * (row[7] + 1) * Math.max(1, row[8]) <= exactOperations) approximateMinimum = 0;
    }
    if (approximateMinimum > prefilter) continue;
    let best = null;
    for (const row of targetRows) {
      // TSXOver rejects two degenerate walk configurations before recording.
      if ((row[7] > 0 && row[8] === 1) || row[7] - row[6] === row[8]) continue;
      const signal = sourceThreeSeqSignal(row, rotation, length, exactOperations);
      if (!best || signal.sourceThreeSeq.rawP < best.sourceThreeSeq.rawP) best = signal;
    }
    if (best) selected.push(best);
  }
  return selected;
}

function selectCoLocatedSiScanRegion(result, candidate, rotation, length, analysisTriplet = null) {
  if (!result) return null;
  return (result.regions?.length ? result.regions : [result])
    .filter((region) => {
      if (!analysisTriplet) return true;
      const roles = sourceSiScanRoles(analysisTriplet, result.baselineTopology, region.inferredTopology);
      return roles?.recombinant === candidate.analysisRecombinant
        && roles?.majorParent === candidate.analysisMajorParent
        && roles?.minorParent === candidate.analysisMinorParent;
    })
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

function sourceSiScanMethodSignal(result, region, analysisTriplet, rotation, length, mappings, maximumPoints = 192) {
  const mapped = mapInterval(region.start, region.end, rotation, length);
  const topologyTriplet = analysisTriplet.map((index) => mappings[index].originIndex);
  const analysisRoles = sourceSiScanRoles(analysisTriplet, result.baselineTopology, region.inferredTopology);
  const resolvedRoles = analysisRoles ? {
    recombinant: mappings[analysisRoles.recombinant].originIndex,
    majorParent: mappings[analysisRoles.majorParent].originIndex,
    minorParent: mappings[analysisRoles.minorParent].originIndex,
  } : null;
  const signal = {
    method: "SiScan",
    ...mapped,
    statistic: region.z,
    locator: `RDP5 Sister-Scanning ${region.scoreFamily} ${region.pattern} topology run`,
    sourceRoutine: result.sourceRoutine,
    outgroup: result.outgroupIndex === null ? null : mappings[result.outgroupIndex].originIndex,
    outgroupMode: result.outgroupMode,
    outgroupSampled: result.outgroupSampled,
    permutations: region.pValuePermutations,
    scanPermutations: region.scanPermutations,
    pattern: region.pattern,
    scoreFamily: region.scoreFamily,
    baselineTopology: result.baselineTopology,
    inferredTopology: region.inferredTopology,
    profile: sourceSiScanProfile(result, rotation, length, maximumPoints),
    sourceSiScan: {
      rawP: region.rawP,
      rawStart: region.start,
      rawEnd: region.end,
      runWindows: Math.max(1, region.last - region.first + 1),
      outgroupSourcePath: result.outgroupSourcePath,
      positionMode: result.positionMode,
      gapMode: result.gapsAsState ? "fifth-state" : "strip",
      window: result.window,
      step: result.step,
      topologyTriplet,
      recombinant: resolvedRoles?.recombinant,
      majorParent: resolvedRoles?.majorParent,
      minorParent: resolvedRoles?.minorParent,
    },
  };
  // Analysis component indexes are intentionally transient. Imported
  // projects retain topologyTriplet above, while a live scan uses this
  // non-enumerable key to reuse the exact once-per-triplet result during
  // characterization without leaking internal component indexes to output.
  Object.defineProperty(signal, "sourceAnalysisTripletKey", {
    value: [...analysisTriplet].sort((left, right) => left - right).join(":"),
    enumerable: false,
  });
  return signal;
}

function applySourceSiScanStats(stats, signal) {
  const source = signal?.sourceSiScan;
  if (!source) return;
  stats.siskanScore = signal.statistic;
  stats.siskanSites = Math.max(0, source.rawEnd - source.rawStart);
  stats.siskanStart = source.rawStart;
  stats.siskanEnd = source.rawEnd;
  stats.siskanRunWindows = source.runWindows;
  stats.siskanSourceP = source.rawP;
  stats.siskanSourceZ = signal.statistic;
  stats.siskanOutgroupIndex = signal.outgroup;
  stats.siskanOutgroupMode = signal.outgroupMode;
  stats.siskanOutgroupSampled = signal.outgroupSampled;
  stats.siskanOutgroupSourcePath = source.outgroupSourcePath;
  stats.siskanPositionMode = source.positionMode;
  stats.siskanGapMode = source.gapMode;
  stats.siskanScanPermutations = signal.scanPermutations;
  stats.siskanPValuePermutations = signal.permutations;
  stats.siskanPattern = signal.pattern;
  stats.siskanScoreFamily = signal.scoreFamily;
  stats.siskanBaselineTopology = signal.baselineTopology;
  stats.siskanInferredTopology = signal.inferredTopology;
  stats.siskanSourceRoutine = signal.sourceRoutine;
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
  if (options.methods.includes("BootScan") && candidate.stats.bootscanWindows > 0 && candidate.stats.bootscanWindows < 8) {
    warnings.push("BootScan topology run contains fewer than eight windows; inspect support stability");
  }
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

async function analyze(message, emit = postMessage) {
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
  const affectedOrigins = new Set(
    Array.isArray(message.affectedOrigins)
      ? message.affectedOrigins.map((value) => Math.trunc(value)).filter((value) => value >= 0 && value < nSeq)
      : [],
  );
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
  const chiSignalCapacity = Math.max(1, Math.min(256, Math.trunc(options.chiSignalsPerTriplet ?? 24)));
  const geneconvSignalCapacity = Math.max(1, Math.min(256, Math.trunc(options.geneconvSignalsPerTriplet ?? 64)));
  const chiTrackCount = 3 * (
    Number(options.methods.includes("MaxChi")) + Number(options.methods.includes("Chimaera"))
  );
  const chiPeakCapacity = Math.max(8, Math.min(256, Math.ceil(
    (chiSignalCapacity * 2) / Math.max(3, chiTrackCount),
  )));
  const rdpBestPtr = outPtr + rdpSignalCapacity * 72;
  const chiPositionsPtr = align(rdpBestPtr + 72, 16);
  const chiScoresPtr = chiPositionsPtr + 4 * (nSites + 1) * 4;
  const chiMissingPrefixPtr = align(chiScoresPtr + 4 * (nSites + 1), 4);
  const chiProfilePtr = align(chiMissingPrefixPtr + (nSites + 1) * 4, 8);
  const chiSmoothPtr = chiProfilePtr + (nSites + 1) * 8;
  const chiPeakPtr = align(chiSmoothPtr + (nSites + 1) * 8, 4);
  const chiOutPtr = chiPeakPtr + chiPeakCapacity * SOURCE_CHI_PEAK_INTS * 4;
  const geneconvPositionsPtr = align(chiOutPtr + chiSignalCapacity * SOURCE_CHI_ROW_INTS * 4, 16);
  const geneconvCategoriesPtr = geneconvPositionsPtr + nSites * 4;
  const geneconvRunStartPtr = align(geneconvCategoriesPtr + nSites, 4);
  const geneconvRunEndPtr = geneconvRunStartPtr + nSites * 4;
  const geneconvRunScorePtr = geneconvRunEndPtr + nSites * 4;
  const geneconvPrefixPtr = align(geneconvRunScorePtr + nSites * 4, 8);
  const geneconvTreePtr = geneconvPrefixPtr + (nSites + 1) * 8;
  // Monotone excursion workspace: i32 stack + f64 maximum + i32 endpoint,
  // with up to seven alignment bytes before the f64 block.
  const geneconvWorkspaceBytes = (nSites + 1) * 16 + 8;
  const geneconvCalibrationPtr = align(geneconvTreePtr + geneconvWorkspaceBytes, 8);
  const geneconvCandidatePtr = align(geneconvCalibrationPtr + 6 * 40, 8);
  const geneconvCandidateCapacity = 3 * (nSites + 1);
  const geneconvDeletePtr = align(geneconvCandidatePtr + geneconvCandidateCapacity * 24, 4);
  const geneconvOutPtr = align(geneconvDeletePtr + nSites * 4, 8);
  const poolPtr = align(geneconvOutPtr + geneconvSignalCapacity * SOURCE_GENECONV_ROW_INTS * 4, 16);
  const nearestIndexesPtr = poolPtr + scanSequenceCount * 4;
  const nearestDistancesPtr = nearestIndexesPtr + scanSequenceCount * 4;
  const roleCohortCapacity = Math.max(4, Math.min(34, Math.trunc(options.roleQuartetTaxaLimit ?? 30)));
  const roleCohortPtr = align(nearestDistancesPtr + scanSequenceCount * 4, 16);
  const roleTractMaskPtr = roleCohortPtr + roleCohortCapacity * 4;
  const roleBackgroundMaskPtr = roleTractMaskPtr + wordsPerSequence * 4;
  const roleDmaxOutPtr = align(roleBackgroundMaskPtr + wordsPerSequence * 4, 8);
  const sourceThreeSeqOutPtr = align(roleDmaxOutPtr + 40, 16);
  const sourceThreeSeqWorkspacePtr = align(sourceThreeSeqOutPtr + 6 * SOURCE_THREE_SEQ_ROW_INTS * 4, 16);
  const sourceThreeSeqWorkspaceBytes = typeof instance.exports.source_three_seq_workspace_bytes === "function"
    ? instance.exports.source_three_seq_workspace_bytes(nSites)
    : 0;
  const requiredBytes = sourceThreeSeqWorkspacePtr + sourceThreeSeqWorkspaceBytes;
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
  const targetSet = new Set(targets);
  const referenceSet = new Set(referencePool);
  new Int32Array(memory.buffer, poolPtr, referencePool.length).set(referencePool);
  const candidates = [];
  const independentMethods = options.methods.filter((method) => method !== "RDP");
  const independentMethodMask = enabledMethodMask({
    ...options,
    methods: independentMethods,
    polishBreakpoints: false,
  });
  const sourceGeneconvEnabled = (independentMethodMask & 2) !== 0;
  const sourceChiMethodMask = independentMethodMask & (4 | 8);
  const sourceBootscanEnabled = options.methods.includes("BootScan");
  const sourceThreeSeqEnabled = options.methods.includes("3Seq");
  const sourceSiScanEnabled = options.methods.includes("SiScan");
  // Every source family below consumes one explicit unordered triplet and
  // resolves its own topology/roles. No outer presumed-target orientation is
  // allowed to turn one biological triplet into three detector calls.
  const sourceOnlyUnorderedPass = options.mode === "exploratory"
    && options.exhaustive;
  const partialBest = new Map();
  const retainCandidate = (candidate) => {
    candidate.structuralUncertaintyVnps = Math.max(1, Math.trunc(options.rdpWindow ?? 30));
    candidate.circular = options.circular === true;
    let retainedPieces = splitCandidateAtStructuralGaps(candidate, disassembly, nSites);
    const fullThreeSeqSignal = (candidate.methodSignals ?? []).find((signal) => signal.sourceThreeSeq);
    if (fullThreeSeqSignal && retainedPieces.length > 1) {
      // CheckSplit3Seq does not emit every side of an interrupted 3Seq walk as
      // a new event. It calls SubPVal for the continuously observed pieces and
      // keeps the lower-p piece. Preserve that source distinction here; other
      // detector families may still retain both structural pieces.
      const sourceSplitPieceCandidates = retainedPieces
        .map((piece) => {
          const signal = refineSourceThreeSeqSignalForPiece(
            fullThreeSeqSignal,
            piece,
            candidate.rotation === 0 ? scanEncoded : rotated,
            nSites,
            Math.max(10_000, Math.trunc(options.threeSeqExactOperations ?? 1_000_000)),
          );
          return signal ? { piece, signal } : null;
        })
        .filter(Boolean);
      let winner = null;
      for (const entry of sourceSplitPieceCandidates) {
        if (!winner || entry.signal.sourceThreeSeq.rawP < winner.signal.sourceThreeSeq.rawP) winner = entry;
      }
      retainedPieces = winner ? [{
        ...winner.piece,
        methodSignals: [
          ...(winner.piece.methodSignals ?? []).filter((signal) => !signal.sourceThreeSeq),
          winner.signal,
        ],
      }] : [];
    }
    for (const piece of retainedPieces) {
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
  const processedRdpTriplets = new Set();
  const countedTruncations = new Set();
  let truncatedRdpSignals = 0;
  const sourceGeneconvTripletCache = new Map();
  const sourceGeneconvTripletCacheLimit = 4_096;
  const processedSourceGeneconvTriplets = new Set();
  const countedGeneconvTruncations = new Set();
  let truncatedGeneconvSignals = 0;
  const sourceChiTripletCache = new Map();
  const sourceChiTripletCacheLimit = 4_096;
  const processedSourceChiTriplets = new Set();
  const processedSourceBootscanTriplets = new Set();
  const processedSourceThreeSeqTriplets = new Set();
  const processedSourceSiScanTriplets = new Set();
  const countedChiTruncations = new Set();
  let truncatedChiSignals = 0;
  let rdpTripletKernelCalls = 0;
  let sourceGeneconvTripletKernelCalls = 0;
  let sourceChiTripletKernelCalls = 0;
  let sourceBootscanBatchCalls = 0;
  let sourceThreeSeqTripletKernelCalls = 0;
  let sourceSiScanTripletCalls = 0;
  let truncatedBootscanSignals = 0;
  let sourceBootscanWorkspaceBytes = 0;
  let sourceBootscanTripletCount = 0;
  let sourceBootscanUsedPairCount = 0;
  const sourceSiScanCache = new Map();
  const sourceThreeSeqTripletCache = new Map();
  const sourceThreeSeqTripletCacheLimit = 4_096;
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
  const getSourceSiScanResult = (triplet, view, candidatePoolInput = null) => {
    if (!sourceSiScanEnabled) return { result: null, triplet: [], cacheKey: "" };
    const ordered = [...triplet].sort((left, right) => left - right);
    const tripletOrigins = new Set(ordered.map((index) => scanMappings[index].originIndex));
    if (tripletOrigins.size !== 3) return { result: null, triplet: ordered, cacheKey: "" };
    const candidatePool = (candidatePoolInput ?? (exactDistanceMatrix ? allIndexes : referencePool))
      .filter((index) => !tripletOrigins.has(scanMappings[index].originIndex));
    const manualOutgroup = options.siskanOutgroupMode === "manual"
      ? allIndexes.find((index) => scanMappings[index].originIndex === options.siskanOutgroupSequence
        && !ordered.includes(index)
        && !tripletOrigins.has(scanMappings[index].originIndex))
      : undefined;
    const cacheKey = `${ordered.join(":")}@${view.rotation}`;
    if (sourceSiScanCache.has(cacheKey)) {
      return { result: sourceSiScanCache.get(cacheKey), triplet: ordered, cacheKey };
    }
    sourceSiScanTripletCalls += 1;
    const seed = (
      (options.randomSeed ?? 0x5a17c0de)
      ^ Math.imul(ordered[0] + 1, 0x9e3779b1)
      ^ Math.imul(ordered[1] + 1, 0x85ebca6b)
      ^ Math.imul(ordered[2] + 1, 0xc2b2ae35)
      ^ Math.imul(view.rotation + 1, 0x27d4eb2f)
    ) >>> 0;
    const result = runSourceSiScan(
      view.rotation === 0 ? scanEncoded : rotated,
      nSites,
      scanSequenceCount,
      ordered,
      {
        window: Math.max(12, options.siskanWindow ?? options.window),
        step: Math.max(1, options.siskanStep ?? options.step),
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
    sourceSiScanCache.set(cacheKey, result);
    if (sourceSiScanCache.size > sourceSiScanCacheLimit) {
      sourceSiScanCache.delete(sourceSiScanCache.keys().next().value);
    }
    return { result, triplet: ordered, cacheKey };
  };
  const sourceBootscanRowsByRotation = new Map();
  let sourceBootscanWindow = Math.max(5, Math.min(32_767, Math.min(Math.floor(nSites / 2), Math.trunc(options.bootscanWindow ?? 200))));
  if (sourceBootscanWindow > nSites / 2) sourceBootscanWindow = Math.max(5, Math.floor(nSites / 2));
  const sourceBootscanStep = Math.max(1, Math.min(Math.max(1, Math.floor(nSites / 4)), Math.trunc(options.bootscanStep ?? 20)));
  const sourceBootscanReplicates = Math.max(0, Math.min(1000, Math.trunc(options.bootstrapReplicates ?? 100)));
  const sourceBootscanTriplets = [];
  const sourceBootscanTripletKeys = new Set();
  const sourceBootscanPairKeys = new Set();
  const sourceBootscanPairIndex = (left, right) => {
    const first = Math.min(left, right);
    const second = Math.max(left, right);
    return first * (2 * scanSequenceCount - first - 1) / 2 + second - first - 1;
  };
  const addSourceBootscanTriplet = (target, firstParent, secondParent) => {
    const ordered = [target, firstParent, secondParent].sort((left, right) => left - right);
    if (ordered[0] === ordered[1] || ordered[1] === ordered[2]) return;
    const origins = ordered.map((index) => scanMappings[index].originIndex);
    if (origins[0] === origins[1] || origins[0] === origins[2] || origins[1] === origins[2]) return;
    if (affectedOrigins.size > 0 && !origins.some((origin) => affectedOrigins.has(origin))) return;
    const key = (ordered[0] * scanSequenceCount + ordered[1]) * scanSequenceCount + ordered[2];
    if (sourceBootscanTripletKeys.has(key)) return;
    sourceBootscanTripletKeys.add(key);
    sourceBootscanTriplets.push(ordered);
    sourceBootscanPairKeys.add(sourceBootscanPairIndex(ordered[0], ordered[1]));
    sourceBootscanPairKeys.add(sourceBootscanPairIndex(ordered[0], ordered[2]));
    sourceBootscanPairKeys.add(sourceBootscanPairIndex(ordered[1], ordered[2]));
  };
  if (sourceBootscanEnabled && sourceBootscanReplicates >= 2) {
    if (sourceOnlyUnorderedPass) {
      // The source primary scan consumes each unordered concrete triplet once.
      // Enumerating it directly avoids generating the same triple under three
      // temporary target orientations only to deduplicate it afterwards.
      for (let first = 0; first < scanSequenceCount - 2; first += 1) {
        for (let second = first + 1; second < scanSequenceCount - 1; second += 1) {
          for (let third = second + 1; third < scanSequenceCount; third += 1) {
            const triplet = [first, second, third];
            const hasEligibleOrientation = triplet.some((candidateTarget) => (
              targetSet.has(candidateTarget)
              && triplet.every((member) => member === candidateTarget || referenceSet.has(member))
            ));
            if (hasEligibleOrientation) addSourceBootscanTriplet(first, second, third);
          }
        }
      }
    } else {
      for (const target of targets) {
        const parentLimit = options.exhaustive
          ? referencePool.length
          : Math.min(Math.max(2, options.candidateParents), referencePool.length);
        let bootParents;
        if (options.exhaustive) {
          bootParents = referencePool.filter((index) => index !== target);
        } else if (exactDistanceMatrix) {
          bootParents = candidateParents(parentDistance, scanSequenceCount, target, referencePool, parentLimit);
        } else {
          bootParents = [];
          for (let slot = 0; slot < parentLimit; slot += 1) {
            const position = Math.min(referencePool.length - 1, Math.floor(((slot + 0.5) * referencePool.length) / Math.max(1, parentLimit)));
            const candidate = referencePool[position];
            if (candidate !== target && !bootParents.includes(candidate)) bootParents.push(candidate);
          }
        }
        for (let left = 0; left < bootParents.length; left += 1) {
          for (let right = left + 1; right < bootParents.length; right += 1) {
            const firstParent = bootParents[left];
            const secondParent = bootParents[right];
            if (options.mode === "query-reference") {
              const firstGroup = scanSequences[firstParent].referenceGroup?.trim();
              const secondGroup = scanSequences[secondParent].referenceGroup?.trim();
              if (firstGroup && secondGroup && firstGroup === secondGroup) continue;
            }
            addSourceBootscanTriplet(target, firstParent, secondParent);
          }
        }
      }
    }
  }
  sourceBootscanTripletCount = sourceBootscanTriplets.length;
  sourceBootscanUsedPairCount = sourceBootscanPairKeys.size;
  if (sourceBootscanEnabled && sourceBootscanReplicates < 2) {
    throw new Error("BootScan requires at least two replicates. Increase Bootstrap replicates or disable BootScan.");
  }
  if (sourceBootscanTripletCount > 0) {
    const pairCount = scanSequenceCount * (scanSequenceCount - 1) / 2;
    const outputCapacity = Math.max(128, Math.min(50_000, Math.trunc(options.bootscanSignals ?? 20_000)));
    const bootTripletPtr = align(roleDmaxOutPtr + 40, 16);
    const bootPairMapPtr = bootTripletPtr + sourceBootscanTripletCount * 12;
    const bootPairListPtr = bootPairMapPtr + pairCount * 4;
    const bootWeightPtr = align(bootPairListPtr + sourceBootscanUsedPairCount * 8, 2);
    const bootPairDistancePtr = bootWeightPtr + sourceBootscanWindow * sourceBootscanReplicates * 2;
    const bootGlobalPairPtr = bootPairDistancePtr + sourceBootscanUsedPairCount * sourceBootscanReplicates * 2;
    const bootStatePtr = align(bootGlobalPairPtr + sourceBootscanUsedPairCount * 8, 4);
    const bootDifferencePtr = bootStatePtr + sourceBootscanTripletCount * 24;
    const bootValidPtr = bootDifferencePtr + sourceBootscanReplicates * 4;
    const bootLookupPtr = align(bootValidPtr + sourceBootscanReplicates * 4, 2);
    const bootLookupEntries = (sourceBootscanWindow + 1) * (sourceBootscanWindow + 2) / 2;
    const bootOutPtr = align(bootLookupPtr + bootLookupEntries * 2, 4);
    const bootRequiredBytes = bootOutPtr + outputCapacity * SOURCE_BOOTSCAN_ROW_INTS * 4;
    sourceBootscanWorkspaceBytes = bootRequiredBytes - bootTripletPtr;
    if (sourceBootscanWorkspaceBytes > 512 * 1024 * 1024) {
      throw new Error(`BootScan batch needs ${(sourceBootscanWorkspaceBytes / (1024 * 1024)).toFixed(0)} MiB for ${sourceBootscanTripletCount.toLocaleString()} concrete triplets. Use the approximate parent shortlist or reduce the active sequence set.`);
    }
    const bootRequiredPages = Math.ceil(bootRequiredBytes / 65536);
    const bootCurrentPages = memory.buffer.byteLength / 65536;
    if (bootRequiredPages > bootCurrentPages) memory.grow(bootRequiredPages - bootCurrentPages);
    new Int32Array(memory.buffer, bootTripletPtr, sourceBootscanTripletCount * 3).set(sourceBootscanTriplets.flat());
    for (const view of scanViews) {
      sourceBootscanBatchCalls += 1;
      const sharedArguments = [
        scanSequenceCount,
        nSites,
        bootTripletPtr,
        sourceBootscanTripletCount,
        sourceBootscanWindow,
        sourceBootscanStep,
        sourceBootscanReplicates,
        Math.round(Math.max(0.5, Math.min(0.999, options.bootscanCutoff ?? 0.7)) * 1000),
        (options.randomSeed ?? 0x5a17c0de) >>> 0,
        bootPairMapPtr,
        bootPairListPtr,
        bootWeightPtr,
        bootPairDistancePtr,
        bootGlobalPairPtr,
        bootStatePtr,
        bootDifferencePtr,
        bootValidPtr,
        bootLookupPtr,
        bootOutPtr,
        outputCapacity,
      ];
      const total = view.rotation === 0
        ? instance.exports.scan_source_bootscan_batch_packed(
            packedPtr,
            validityPtr,
            wordsPerSequence,
            ...sharedArguments,
          )
        : instance.exports.scan_source_bootscan_batch(view.sequencePtr, ...sharedArguments);
      const retained = Math.min(outputCapacity, Math.max(0, total));
      if (total > outputCapacity) truncatedBootscanSignals += total - outputCapacity;
      const rowsByTriplet = new Map();
      for (let rowIndex = 0; rowIndex < retained; rowIndex += 1) {
        const row = Array.from(new Int32Array(
          memory.buffer,
          bootOutPtr + rowIndex * SOURCE_BOOTSCAN_ROW_INTS * 4,
          SOURCE_BOOTSCAN_ROW_INTS,
        ));
        const key = `${row[0]}:${row[1]}:${row[2]}`;
        rowsByTriplet.set(key, [...(rowsByTriplet.get(key) ?? []), row]);
      }
      sourceBootscanRowsByRotation.set(view.rotation, rowsByTriplet);
    }
  }
  const getSourceBootscanRows = (triplet, view) => {
    if (!sourceBootscanEnabled) return [];
    const key = [...triplet].sort((left, right) => left - right).join(":");
    return sourceBootscanRowsByRotation.get(view.rotation)?.get(key) ?? [];
  };
  const getSourceThreeSeqRows = (triplet, view) => {
    if (!sourceThreeSeqEnabled) return [];
    const ordered = [...triplet].sort((left, right) => left - right);
    const cacheKey = `${ordered.join(":")}@${view.rotation}`;
    let cached = sourceThreeSeqTripletCache.get(cacheKey);
    if (!cached) {
      sourceThreeSeqTripletKernelCalls += 1;
      // FindSubSeqTS visits essentially every column, so the contiguous byte
      // layout is measurably faster than repeated packed lane extraction on
      // current browser/Node WASM engines. Rotated and ordinary views both
      // already live in this layout; retain the packed export as an exact
      // oracle and future SIMD target.
      const total = typeof instance.exports.scan_source_three_seq_triplet_mode === "function"
        ? instance.exports.scan_source_three_seq_triplet_mode(
            view.sequencePtr,
            nSites,
            ordered[0],
            ordered[1],
            ordered[2],
            options.circular ? 1 : 0,
            sourceThreeSeqWorkspacePtr,
            sourceThreeSeqOutPtr,
          )
        : instance.exports.scan_source_three_seq_triplet(
            view.sequencePtr,
            nSites,
            ordered[0],
            ordered[1],
            ordered[2],
            sourceThreeSeqOutPtr,
          );
      cached = Array.from({ length: Math.min(6, Math.max(0, total)) }, (_, index) => (
        Array.from(new Int32Array(
          memory.buffer,
          sourceThreeSeqOutPtr + index * SOURCE_THREE_SEQ_ROW_INTS * 4,
          SOURCE_THREE_SEQ_ROW_INTS,
        ))
      ));
      sourceThreeSeqTripletCache.set(cacheKey, cached);
      if (sourceThreeSeqTripletCache.size > sourceThreeSeqTripletCacheLimit) {
        sourceThreeSeqTripletCache.delete(sourceThreeSeqTripletCache.keys().next().value);
      }
    }
    return cached;
  };
  const getSourceGeneconvRows = (triplet, view) => {
    if (!sourceGeneconvEnabled) return [];
    const ordered = [...triplet].sort((left, right) => left - right);
    const cacheKey = `${ordered.join(":")}@${view.rotation}`;
    let cached = sourceGeneconvTripletCache.get(cacheKey);
    if (!cached) {
      sourceGeneconvTripletKernelCalls += 1;
      const sharedArguments = [
        nSites,
        ordered[0],
        ordered[1],
        ordered[2],
        Math.max(0, Math.round(options.geneconvGScale ?? 1)),
        Math.max(Number.MIN_VALUE, Math.min(1, options.alpha ?? 0.05)),
        geneconvPositionsPtr,
        geneconvCategoriesPtr,
        geneconvRunStartPtr,
        geneconvRunEndPtr,
        geneconvRunScorePtr,
        geneconvPrefixPtr,
        geneconvTreePtr,
        geneconvCalibrationPtr,
        geneconvCandidatePtr,
        geneconvCandidateCapacity,
        geneconvDeletePtr,
        geneconvOutPtr,
        geneconvSignalCapacity,
      ];
      const total = view.rotation === 0 && typeof instance.exports.scan_source_geneconv_all_packed === "function"
        ? instance.exports.scan_source_geneconv_all_packed(
            packedPtr,
            validityPtr,
            wordsPerSequence,
            ...sharedArguments,
          )
        : instance.exports.scan_source_geneconv_all(
            view.sequencePtr,
            ...sharedArguments,
          );
      const retained = Math.min(geneconvSignalCapacity, Math.max(0, total));
      cached = {
        ordered,
        total,
        rows: Array.from({ length: retained }, (_, signalIndex) => (
          Array.from(new Int32Array(
            memory.buffer,
            geneconvOutPtr + signalIndex * SOURCE_GENECONV_ROW_INTS * 4,
            SOURCE_GENECONV_ROW_INTS,
          ))
        )),
      };
      sourceGeneconvTripletCache.set(cacheKey, cached);
      if (sourceGeneconvTripletCache.size > sourceGeneconvTripletCacheLimit) {
        sourceGeneconvTripletCache.delete(sourceGeneconvTripletCache.keys().next().value);
      }
    }
    if (cached.total > geneconvSignalCapacity && !countedGeneconvTruncations.has(cacheKey)) {
      truncatedGeneconvSignals += cached.total - geneconvSignalCapacity;
      countedGeneconvTruncations.add(cacheKey);
    }
    return cached.rows;
  };
  const getSourceChiRows = (triplet, view) => {
    if (sourceChiMethodMask === 0) return [];
    const ordered = [...triplet].sort((left, right) => left - right);
    const cacheKey = `${ordered.join(":")}@${view.rotation}`;
    let cached = sourceChiTripletCache.get(cacheKey);
    if (!cached) {
      sourceChiTripletKernelCalls += 1;
      const total = view.rotation === 0 && typeof instance.exports.scan_source_chi_all_packed === "function"
        ? instance.exports.scan_source_chi_all_packed(
            packedPtr,
            validityPtr,
            wordsPerSequence,
            nSites,
            ordered[0],
            ordered[1],
            ordered[2],
            Math.max(20, options.window),
            0,
            sourceChiMethodMask,
            chiPositionsPtr,
            chiScoresPtr,
            chiMissingPrefixPtr,
            chiProfilePtr,
            chiSmoothPtr,
            chiPeakPtr,
            chiPeakCapacity,
            chiOutPtr,
            chiSignalCapacity,
          )
        : instance.exports.scan_source_chi_all(
            view.sequencePtr,
            nSites,
            ordered[0],
            ordered[1],
            ordered[2],
            Math.max(20, options.window),
            0,
            sourceChiMethodMask,
            chiPositionsPtr,
            chiScoresPtr,
            chiMissingPrefixPtr,
            chiProfilePtr,
            chiSmoothPtr,
            chiPeakPtr,
            chiPeakCapacity,
            chiOutPtr,
            chiSignalCapacity,
          );
      const retained = Math.min(chiSignalCapacity, Math.max(0, total));
      cached = {
        ordered,
        total,
        rows: Array.from({ length: retained }, (_, signalIndex) => (
          Array.from(new Int32Array(
            memory.buffer,
            chiOutPtr + signalIndex * SOURCE_CHI_ROW_INTS * 4,
            SOURCE_CHI_ROW_INTS,
          ))
        )),
      };
      sourceChiTripletCache.set(cacheKey, cached);
      if (sourceChiTripletCache.size > sourceChiTripletCacheLimit) {
        sourceChiTripletCache.delete(sourceChiTripletCache.keys().next().value);
      }
    }
    if (cached.total > chiSignalCapacity && !countedChiTruncations.has(cacheKey)) {
      truncatedChiSignals += cached.total - chiSignalCapacity;
      countedChiTruncations.add(cacheKey);
    }
    return cached.rows;
  };

  for (let targetPosition = 0; targetPosition < targets.length; targetPosition += 1) {
    const analysisRecombinant = targets[targetPosition];
    const recombinant = scanMappings[analysisRecombinant].originIndex;
    const parentLimit = options.exhaustive ? referencePool.length : Math.min(options.candidateParents, referencePool.length);
    let parents;
    if (options.exhaustive) {
      // AlistRDP/AlistMC/AlistChi consume one explicit unordered triplet and
      // resolve their internal pair tracks/orientations themselves. In the
      // role-agnostic source-only path, a<b<c therefore enumerates the full
      // search exactly once instead of revisiting it under three outer roles.
      parents = referencePool.filter((index) => sourceOnlyUnorderedPass
        ? index > analysisRecombinant
        : index !== analysisRecombinant);
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
    const reserve = sourceOnlyUnorderedPass
      ? 0
      : Math.min(groupRepresentatives.size, Math.max(0, Math.floor(parentLimit / 3)));
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
        // RDP5's redo list re-examines only triplets containing a sequence
        // changed by the latest erase/extract operation. An empty filter is
        // the initial complete a<b<c screen.
        if (affectedOrigins.size > 0 && !inputOrigins.some((origin) => affectedOrigins.has(origin))) continue;
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
            const globalConcretePass = options.mode === "exploratory" && options.exhaustive;
            if (!globalConcretePass || !processedRdpTriplets.has(rdpCacheKey)) {
              if (globalConcretePass) processedRdpTriplets.add(rdpCacheKey);
              let cachedRdp = rdpTripletCache.get(rdpCacheKey);
              if (!cachedRdp) {
                rdpTripletKernelCalls += 1;
                const totalRdpSignals = view.rotation === 0 && typeof instance.exports.scan_rdp5_triplet_all_packed === "function"
                  ? instance.exports.scan_rdp5_triplet_all_packed(
                      packedPtr,
                      validityPtr,
                      wordsPerSequence,
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
                    )
                  : instance.exports.scan_rdp5_triplet_all(
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
                const targetEligible = globalConcretePass
                  ? targetSet.has(output[2]) && referenceSet.has(output[3]) && referenceSet.has(output[4])
                  : output[2] === analysisRecombinant;
                if (targetEligible) {
                  const rawStart = output[0];
                  const rawEnd = output[1];
                  const mapped = mapInterval(rawStart, rawEnd, view.rotation, nSites);
                  const majorParent = scanMappings[output[3]].originIndex;
                  const minorParent = scanMappings[output[4]].originIndex;
                  const candidate = {
                    recombinant: scanMappings[output[2]].originIndex,
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
                    siskanCandidatePool: exactDistanceMatrix ? undefined : [...referencePool],
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
          }

          if (sourceGeneconvEnabled) {
            const triplet = [analysisRecombinant, inputMajorParent, inputMinorParent];
            const orderedTriplet = [...triplet].sort((left, right) => left - right);
            const globalConcretePass = options.mode === "exploratory" && options.exhaustive;
            const processedKey = `${orderedTriplet.join(":")}@${view.rotation}`;
            if (!globalConcretePass || !processedSourceGeneconvTriplets.has(processedKey)) {
              if (globalConcretePass) processedSourceGeneconvTriplets.add(processedKey);
              for (const row of getSourceGeneconvRows(triplet, view)) {
                const sourceTarget = orderedTriplet[row[1]];
                const sourceMinor = orderedTriplet[row[2]];
                const sourceMajor = orderedTriplet[row[3]];
                const targetEligible = globalConcretePass
                  ? targetSet.has(sourceTarget) && referenceSet.has(sourceMinor) && referenceSet.has(sourceMajor)
                  : sourceTarget === analysisRecombinant;
                if (!targetEligible) continue;
                const rawStart = row[4];
                const rawEnd = row[5];
                if (!(rawEnd - rawStart >= 4 && rawStart >= 0 && rawEnd <= nSites)) continue;
                instance.exports.triplet_counts(
                  view.sequencePtr,
                  nSites,
                  sourceTarget,
                  sourceMajor,
                  sourceMinor,
                  rawStart,
                  rawEnd,
                  rdpBestPtr,
                );
                const counts = Array.from(new Int32Array(memory.buffer, rdpBestPtr, 6));
                if (counts[0] < 4) continue;
                const mapped = mapInterval(rawStart, rawEnd, view.rotation, nSites);
                const insideTotal = counts[1] + counts[2];
                const outsideTotal = counts[3] + counts[4];
                const effect = counts[1] / Math.max(1, insideTotal)
                  - counts[4] / Math.max(1, outsideTotal);
                retainCandidate({
                  recombinant: scanMappings[sourceTarget].originIndex,
                  ...mapped,
                  rawStart,
                  rawEnd,
                  sequencePtr: view.sequencePtr,
                  rotation: view.rotation,
                  majorParent: scanMappings[sourceMajor].originIndex,
                  minorParent: scanMappings[sourceMinor].originIndex,
                  analysisRecombinant: sourceTarget,
                  analysisMajorParent: sourceMajor,
                  analysisMinorParent: sourceMinor,
                  siskanCandidatePool: exactDistanceMatrix ? undefined : [...referencePool],
                  componentProvenance: candidateComponentProvenance(
                    disassembly,
                    sourceTarget,
                    sourceMajor,
                    sourceMinor,
                  ),
                  chiSquare: counts[5] / 1000,
                  informative: counts[0],
                  insideMinor: counts[1],
                  insideMajor: counts[2],
                  outsideMajor: counts[3],
                  outsideMinor: counts[4],
                  effect,
                  alternatives: [],
                  methodSignals: [sourceGeneconvSignal(row, view.rotation, nSites)],
                });
              }
            }
          }

          if (sourceBootscanEnabled) {
            const triplet = [analysisRecombinant, inputMajorParent, inputMinorParent];
            const orderedTriplet = [...triplet].sort((left, right) => left - right);
            const globalConcretePass = options.mode === "exploratory" && options.exhaustive;
            const processedKey = `${orderedTriplet.join(":")}@${view.rotation}`;
            if (!globalConcretePass || !processedSourceBootscanTriplets.has(processedKey)) {
              if (globalConcretePass) processedSourceBootscanTriplets.add(processedKey);
              const sourceTargets = globalConcretePass
                ? orderedTriplet.filter((target) => (
                    targetSet.has(target)
                    && orderedTriplet.every((member) => member === target || referenceSet.has(member))
                  ))
                : [analysisRecombinant];
              for (const row of getSourceBootscanRows(triplet, view)) {
                const rawStart = row[4];
                const rawEnd = row[5];
                if (!(rawEnd - rawStart >= 4 && rawStart >= 0 && rawEnd <= nSites)) continue;
                const sourceSignal = sourceBootscanSignal(row, view.rotation, nSites);
                if (sourceSignal.sourceBootscan.rawP > Math.max(Number.MIN_VALUE, options.alpha ?? 0.05)) continue;
                const sourceRoles = sourceBootscanRoles(row);
                if (!sourceRoles
                  || !sourceTargets.includes(sourceRoles.recombinant)
                  || !referenceSet.has(sourceRoles.majorParent)
                  || !referenceSet.has(sourceRoles.minorParent)) continue;
                instance.exports.triplet_counts(
                  view.sequencePtr,
                  nSites,
                  sourceRoles.recombinant,
                  sourceRoles.majorParent,
                  sourceRoles.minorParent,
                  rawStart,
                  rawEnd,
                  rdpBestPtr,
                );
                const sourceCounts = Array.from(new Int32Array(memory.buffer, rdpBestPtr, 6));
                const insideTotal = sourceCounts[1] + sourceCounts[2];
                const outsideTotal = sourceCounts[3] + sourceCounts[4];
                const effect = sourceCounts[1] / Math.max(1, insideTotal)
                  - sourceCounts[4] / Math.max(1, outsideTotal);
                if (sourceCounts[0] < 4 || !(effect > 0)) continue;
                const mapped = mapInterval(rawStart, rawEnd, view.rotation, nSites);
                retainCandidate({
                  recombinant: scanMappings[sourceRoles.recombinant].originIndex,
                  ...mapped,
                  rawStart,
                  rawEnd,
                  sequencePtr: view.sequencePtr,
                  rotation: view.rotation,
                  majorParent: scanMappings[sourceRoles.majorParent].originIndex,
                  minorParent: scanMappings[sourceRoles.minorParent].originIndex,
                  analysisRecombinant: sourceRoles.recombinant,
                  analysisMajorParent: sourceRoles.majorParent,
                  analysisMinorParent: sourceRoles.minorParent,
                  siskanCandidatePool: exactDistanceMatrix ? undefined : [...referencePool],
                  componentProvenance: candidateComponentProvenance(
                    disassembly,
                    sourceRoles.recombinant,
                    sourceRoles.majorParent,
                    sourceRoles.minorParent,
                  ),
                  chiSquare: sourceCounts[5] / 1000,
                  informative: sourceCounts[0],
                  insideMinor: sourceCounts[1],
                  insideMajor: sourceCounts[2],
                  outsideMajor: sourceCounts[3],
                  outsideMinor: sourceCounts[4],
                  effect,
                  alternatives: [],
                  methodSignals: [sourceSignal],
                });
              }
            }
          }

          // CheckwrapC already performs the source's bounded origin-crossing
          // extension. Do not feed 3Seq the synthetic half-genome origin as a
          // second biological test: it can change which complementary walk
          // wins and duplicates work. Keep the legacy rotated fallback only
          // for older WASM builds that do not expose the CheckwrapC mode.
          if (sourceThreeSeqEnabled && (
            !options.circular
            || view.rotation === 0
            || typeof instance.exports.scan_source_three_seq_triplet_mode !== "function"
          )) {
            const triplet = [analysisRecombinant, inputMajorParent, inputMinorParent];
            const orderedTriplet = [...triplet].sort((left, right) => left - right);
            const globalConcretePass = options.mode === "exploratory" && options.exhaustive;
            const processedKey = `${orderedTriplet.join(":")}@${view.rotation}`;
            if (!globalConcretePass || !processedSourceThreeSeqTriplets.has(processedKey)) {
              if (globalConcretePass) processedSourceThreeSeqTriplets.add(processedKey);
              const sourceTargets = globalConcretePass
                ? orderedTriplet.filter((target) => (
                    targetSet.has(target)
                    && orderedTriplet.every((member) => member === target || referenceSet.has(member))
                  ))
                : [analysisRecombinant];
              const selectedSignals = selectSourceThreeSeqSignals(
                getSourceThreeSeqRows(triplet, view),
                view.rotation,
                nSites,
                Math.max(Number.MIN_VALUE, options.alpha ?? 0.05),
                Math.max(10_000, Math.trunc(options.threeSeqExactOperations ?? 1_000_000)),
              );
              for (const sourceSignal of selectedSignals) {
                const sourceRoles = sourceSignal.sourceThreeSeq;
                if (sourceSignal.sourceThreeSeq.rawP > Math.max(Number.MIN_VALUE, options.alpha ?? 0.05)) continue;
                if (!sourceTargets.includes(sourceRoles.target)
                  || !referenceSet.has(sourceRoles.majorParent)
                  || !referenceSet.has(sourceRoles.minorParent)) continue;
                const kernelStart = sourceRoles.rawStart;
                const kernelEnd = sourceRoles.rawEnd;
                const kernelLength = rawIntervalLength(kernelStart, kernelEnd, nSites);
                if (!(kernelLength >= 4 && kernelStart >= 0 && kernelStart < nSites && kernelEnd >= 0 && kernelEnd <= nSites)) continue;
                let bestOrientation = null;
                for (const [majorParent, minorParent] of [
                  [sourceRoles.majorParent, sourceRoles.minorParent],
                  [sourceRoles.minorParent, sourceRoles.majorParent],
                ]) {
                  instance.exports.triplet_counts(
                    view.sequencePtr,
                    nSites,
                    sourceRoles.target,
                    majorParent,
                    minorParent,
                    kernelStart,
                    kernelEnd,
                    rdpBestPtr,
                  );
                  const counts = Array.from(new Int32Array(memory.buffer, rdpBestPtr, 6));
                  const insideTotal = counts[1] + counts[2];
                  const outsideTotal = counts[3] + counts[4];
                  const effect = counts[1] / Math.max(1, insideTotal)
                    - counts[4] / Math.max(1, outsideTotal);
                  if (!bestOrientation || effect > bestOrientation.effect) {
                    bestOrientation = { majorParent, minorParent, counts, effect };
                  }
                }
                if (!bestOrientation || bestOrientation.counts[0] < 4 || !(bestOrientation.effect > 0)) continue;
                const mapped = mapInterval(kernelStart, kernelEnd, view.rotation, nSites);
                const mappedSignal = {
                  ...sourceSignal,
                  ...mapped,
                  sourceThreeSeq: {
                    ...sourceSignal.sourceThreeSeq,
                    target: scanMappings[sourceRoles.target].originIndex,
                    // Preserve the source walk's +1/-1 parent orientation in
                    // its ledger. Event major/minor roles are inferred below
                    // from the tract contrast and may legitimately be swapped.
                    majorParent: scanMappings[sourceRoles.majorParent].originIndex,
                    minorParent: scanMappings[sourceRoles.minorParent].originIndex,
                  },
                };
                retainCandidate({
                  recombinant: scanMappings[sourceRoles.target].originIndex,
                  ...mapped,
                  rawStart: kernelStart,
                  rawEnd: kernelEnd,
                  sequencePtr: view.sequencePtr,
                  rotation: view.rotation,
                  majorParent: scanMappings[bestOrientation.majorParent].originIndex,
                  minorParent: scanMappings[bestOrientation.minorParent].originIndex,
                  analysisRecombinant: sourceRoles.target,
                  analysisMajorParent: bestOrientation.majorParent,
                  analysisMinorParent: bestOrientation.minorParent,
                  siskanCandidatePool: exactDistanceMatrix ? undefined : [...referencePool],
                  componentProvenance: candidateComponentProvenance(
                    disassembly,
                    sourceRoles.target,
                    bestOrientation.majorParent,
                    bestOrientation.minorParent,
                  ),
                  chiSquare: bestOrientation.counts[5] / 1000,
                  informative: bestOrientation.counts[0],
                  insideMinor: bestOrientation.counts[1],
                  insideMajor: bestOrientation.counts[2],
                  outsideMajor: bestOrientation.counts[3],
                  outsideMinor: bestOrientation.counts[4],
                  effect: bestOrientation.effect,
                  alternatives: [],
                  methodSignals: [mappedSignal],
                });
              }
            }
          }

          if (sourceChiMethodMask !== 0) {
            const triplet = [analysisRecombinant, inputMajorParent, inputMinorParent];
            const orderedTriplet = [...triplet].sort((left, right) => left - right);
            const globalConcretePass = options.mode === "exploratory" && options.exhaustive;
            const processedKey = `${orderedTriplet.join(":")}@${view.rotation}`;
            if (!globalConcretePass || !processedSourceChiTriplets.has(processedKey)) {
              if (globalConcretePass) processedSourceChiTriplets.add(processedKey);
              const sourceRows = getSourceChiRows(triplet, view);
              const sourceTargets = globalConcretePass
                ? orderedTriplet.filter((target) => (
                    targetSet.has(target)
                    && orderedTriplet.every((member) => member === target || referenceSet.has(member))
                  ))
                : [analysisRecombinant];
              for (const sourceTarget of sourceTargets) {
                const targetSlot = orderedTriplet.indexOf(sourceTarget);
                const applicableRows = new Map();
                for (const row of sourceRows) {
                  if (row[0] === 4 && row[1] !== targetSlot) continue;
                  const key = `${row[0]}:${row[3]}:${row[4]}:${row[5]}`;
                  const previous = applicableRows.get(key);
                  if (!previous || row[6] > previous[6]) applicableRows.set(key, row);
                }
                const sourceParents = orderedTriplet.filter((index) => index !== sourceTarget);
                for (const row of applicableRows.values()) {
                  const rawStart = row[3];
                  const rawEnd = row[4];
                  if (!(rawEnd - rawStart >= 4 && rawStart >= 0 && rawEnd <= nSites)) continue;
                  let bestOrientation = null;
                  for (const [analysisMajorParent, analysisMinorParent] of [
                    [sourceParents[0], sourceParents[1]],
                    [sourceParents[1], sourceParents[0]],
                  ]) {
                    instance.exports.triplet_counts(
                      view.sequencePtr,
                      nSites,
                      sourceTarget,
                      analysisMajorParent,
                      analysisMinorParent,
                      rawStart,
                      rawEnd,
                      rdpBestPtr,
                    );
                    const counts = Array.from(new Int32Array(memory.buffer, rdpBestPtr, 6));
                    const insideTotal = counts[1] + counts[2];
                    const outsideTotal = counts[3] + counts[4];
                    const effect = counts[1] / Math.max(1, insideTotal) - counts[4] / Math.max(1, outsideTotal);
                    if (!bestOrientation || effect > bestOrientation.effect) {
                      bestOrientation = { analysisMajorParent, analysisMinorParent, counts, effect };
                    }
                  }
                  if (!bestOrientation || bestOrientation.counts[0] < 4) continue;
                  const mapped = mapInterval(rawStart, rawEnd, view.rotation, nSites);
                  const signal = sourceChiSignal(row, view.rotation, nSites);
                  retainCandidate({
                    recombinant: scanMappings[sourceTarget].originIndex,
                    ...mapped,
                    rawStart,
                    rawEnd,
                    sequencePtr: view.sequencePtr,
                    rotation: view.rotation,
                    majorParent: scanMappings[bestOrientation.analysisMajorParent].originIndex,
                    minorParent: scanMappings[bestOrientation.analysisMinorParent].originIndex,
                    analysisRecombinant: sourceTarget,
                    analysisMajorParent: bestOrientation.analysisMajorParent,
                    analysisMinorParent: bestOrientation.analysisMinorParent,
                    siskanCandidatePool: exactDistanceMatrix ? undefined : [...referencePool],
                    componentProvenance: candidateComponentProvenance(
                      disassembly,
                      sourceTarget,
                      bestOrientation.analysisMajorParent,
                      bestOrientation.analysisMinorParent,
                    ),
                    chiSquare: bestOrientation.counts[5] / 1000,
                    informative: bestOrientation.counts[0],
                    insideMinor: bestOrientation.counts[1],
                    insideMajor: bestOrientation.counts[2],
                    outsideMajor: bestOrientation.counts[3],
                    outsideMinor: bestOrientation.counts[4],
                    effect: bestOrientation.effect,
                    alternatives: [],
                    methodSignals: [signal],
                  });
                }
              }
            }
          }

          if (sourceSiScanEnabled) {
            const triplet = [analysisRecombinant, inputMajorParent, inputMinorParent];
            const orderedTriplet = [...triplet].sort((left, right) => left - right);
            const processedKey = `${orderedTriplet.join(":")}@${view.rotation}`;
            if (!processedSourceSiScanTriplets.has(processedKey)) {
              processedSourceSiScanTriplets.add(processedKey);
              const source = getSourceSiScanResult(
                orderedTriplet,
                view,
                exactDistanceMatrix ? allIndexes : referencePool,
              );
              for (const region of source.result?.regions ?? []) {
                if (region.rawP > Math.max(Number.MIN_VALUE, options.alpha ?? 0.05)) continue;
                const sourceRoles = sourceSiScanRoles(
                  source.triplet,
                  source.result.baselineTopology,
                  region.inferredTopology,
                );
                if (!sourceRoles
                  || !targetSet.has(sourceRoles.recombinant)
                  || !referenceSet.has(sourceRoles.majorParent)
                  || !referenceSet.has(sourceRoles.minorParent)) continue;
                const rawStart = region.start;
                const rawEnd = region.end;
                if (!(rawEnd - rawStart >= 4 && rawStart >= 0 && rawEnd <= nSites)) continue;
                instance.exports.triplet_counts(
                  view.sequencePtr,
                  nSites,
                  sourceRoles.recombinant,
                  sourceRoles.majorParent,
                  sourceRoles.minorParent,
                  rawStart,
                  rawEnd,
                  rdpBestPtr,
                );
                const counts = Array.from(new Int32Array(memory.buffer, rdpBestPtr, 6));
                const insideTotal = counts[1] + counts[2];
                const outsideTotal = counts[3] + counts[4];
                const effect = counts[1] / Math.max(1, insideTotal)
                  - counts[4] / Math.max(1, outsideTotal);
                if (counts[0] < 4 || !(effect > 0)) continue;
                const mapped = mapInterval(rawStart, rawEnd, view.rotation, nSites);
                const sourceSignal = sourceSiScanMethodSignal(
                  source.result,
                  region,
                  source.triplet,
                  view.rotation,
                  nSites,
                  scanMappings,
                  candidates.length > 1_000 ? 48 : candidates.length > 250 ? 96 : 192,
                );
                retainCandidate({
                  recombinant: scanMappings[sourceRoles.recombinant].originIndex,
                  ...mapped,
                  rawStart,
                  rawEnd,
                  sequencePtr: view.sequencePtr,
                  rotation: view.rotation,
                  majorParent: scanMappings[sourceRoles.majorParent].originIndex,
                  minorParent: scanMappings[sourceRoles.minorParent].originIndex,
                  analysisRecombinant: sourceRoles.recombinant,
                  analysisMajorParent: sourceRoles.majorParent,
                  analysisMinorParent: sourceRoles.minorParent,
                  siskanCandidatePool: exactDistanceMatrix ? undefined : [...referencePool],
                  componentProvenance: candidateComponentProvenance(
                    disassembly,
                    sourceRoles.recombinant,
                    sourceRoles.majorParent,
                    sourceRoles.minorParent,
                  ),
                  chiSquare: counts[5] / 1000,
                  informative: counts[0],
                  insideMinor: counts[1],
                  insideMajor: counts[2],
                  outsideMajor: counts[3],
                  outsideMinor: counts[4],
                  effect,
                  alternatives: [],
                  methodSignals: [sourceSignal],
                });
              }
            }
          }

        }
      }
    }
    emit({
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
      emit({ type: "partial", jobId, events: partialEvents, comparisons });
    }
  }
  const scanMs = performance.now() - scanStarted;

  const unique = deduplicate(candidates, nSites, targets.length);
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
    const candidateTriplet = [candidate.analysisRecombinant, candidate.analysisMajorParent, candidate.analysisMinorParent];
    const orderedCandidateTriplet = [...candidateTriplet].sort((left, right) => left - right);
    const candidateTargetSlot = orderedCandidateTriplet.indexOf(candidate.analysisRecombinant);
    let sourceGeneconv = null;
    for (const row of getSourceGeneconvRows(candidateTriplet, {
      sequencePtr: candidate.sequencePtr,
      rotation: candidate.rotation,
    })) {
      if (row[1] !== candidateTargetSlot) continue;
      const signal = sourceGeneconvSignal(row, candidate.rotation, nSites);
      if (!isCoLocatedSignal(candidate, signal, nSites)) continue;
      if (!sourceGeneconv || signal.sourceGeneconv.rawP < sourceGeneconv.sourceGeneconv.rawP) {
        sourceGeneconv = signal;
      }
    }
    if (sourceGeneconv) {
      candidate.methodSignals = [
        ...(candidate.methodSignals ?? []).filter((entry) => entry.method !== "GENECONV"),
        sourceGeneconv,
      ];
    }
    let sourceBootscan = null;
    for (const row of getSourceBootscanRows(candidateTriplet, {
      sequencePtr: candidate.sequencePtr,
      rotation: candidate.rotation,
    })) {
      const signal = sourceBootscanSignal(row, candidate.rotation, nSites);
      if (!isCoLocatedSignal(candidate, signal, nSites)) continue;
      if (!sourceBootscan || signal.sourceBootscan.rawP < sourceBootscan.sourceBootscan.rawP) {
        sourceBootscan = signal;
      }
    }
    if (sourceBootscan) {
      candidate.methodSignals = [
        ...(candidate.methodSignals ?? []).filter((entry) => entry.method !== "BootScan"),
        sourceBootscan,
      ];
    }
    const sourceChiByMethod = new Map();
    for (const row of getSourceChiRows(candidateTriplet, {
      sequencePtr: candidate.sequencePtr,
      rotation: candidate.rotation,
    })) {
      if (row[0] === 4 && row[1] !== candidateTargetSlot) continue;
      const signal = sourceChiSignal(row, candidate.rotation, nSites);
      if (!isCoLocatedSignal(candidate, signal, nSites)) continue;
      const previous = sourceChiByMethod.get(signal.method);
      if (!previous || signal.statistic > previous.statistic) sourceChiByMethod.set(signal.method, signal);
    }
    for (const signal of sourceChiByMethod.values()) {
      if (!(candidate.methodSignals ?? []).some((entry) => (
        entry.method === signal.method && entry.start === signal.start && entry.end === signal.end
      ))) {
        candidate.methodSignals = [...(candidate.methodSignals ?? []), signal];
      }
    }
    const sourceMaxChi = sourceChiByMethod.get("MaxChi")
      ?? (candidate.methodSignals ?? []).find((signal) => signal.method === "MaxChi" && signal.sourceChi);
    const sourceChimaera = sourceChiByMethod.get("Chimaera")
      ?? (candidate.methodSignals ?? []).find((signal) => signal.method === "Chimaera" && signal.sourceChi);
    let sourceThreeSeq = null;
    for (const fullSignal of selectSourceThreeSeqSignals(
      getSourceThreeSeqRows(candidateTriplet, {
        sequencePtr: candidate.sequencePtr,
        rotation: candidate.rotation,
      }),
      candidate.rotation,
      nSites,
      Math.max(Number.MIN_VALUE, options.alpha ?? 0.05),
      Math.max(10_000, Math.trunc(options.threeSeqExactOperations ?? 1_000_000)),
    )) {
      const signal = refineSourceThreeSeqSignalForPiece(
        fullSignal,
        candidate,
        candidate.rotation === 0 ? scanEncoded : rotated,
        nSites,
        Math.max(10_000, Math.trunc(options.threeSeqExactOperations ?? 1_000_000)),
      );
      if (!signal) continue;
      if (signal.sourceThreeSeq.target !== candidate.analysisRecombinant) continue;
      if (!isCoLocatedSignal(candidate, signal, nSites)) continue;
      if (!sourceThreeSeq || signal.sourceThreeSeq.rawP < sourceThreeSeq.sourceThreeSeq.rawP) {
        sourceThreeSeq = {
          ...signal,
          sourceThreeSeq: {
            ...signal.sourceThreeSeq,
            target: scanMappings[signal.sourceThreeSeq.target].originIndex,
            majorParent: scanMappings[signal.sourceThreeSeq.majorParent].originIndex,
            minorParent: scanMappings[signal.sourceThreeSeq.minorParent].originIndex,
          },
        };
      }
    }
    if (sourceThreeSeq) {
      candidate.methodSignals = [
        ...(candidate.methodSignals ?? []).filter((entry) => entry.method !== "3Seq"),
        sourceThreeSeq,
      ];
    }
    candidate.stats = {
      genconvRun: sourceGeneconv?.sourceGeneconv.fragmentScore ?? sourceGeneconv?.statistic ?? 0,
      genconvEligible: sourceGeneconv?.sourceGeneconv.informativeSites ?? 0,
      genconvMatches: sourceGeneconv?.sourceGeneconv.matchingSites ?? 0,
      genconvStart: sourceGeneconv ? sourceGeneconv.start : 0,
      genconvEnd: sourceGeneconv ? sourceGeneconv.end : 0,
      bootscanConsistent: sourceBootscan?.sourceBootscan.runWindows ?? 0,
      bootscanWindows: sourceBootscan?.sourceBootscan.runWindows ?? 0,
      maxChi: sourceMaxChi?.statistic ?? 0,
      chimaera: sourceChimaera?.statistic ?? 0,
      siskanScore: 0,
      siskanSites: 0,
      threeSeqDescent: sourceThreeSeq?.sourceThreeSeq.descent ?? 0,
      threeSeqSites: sourceThreeSeq?.sourceThreeSeq.informativeSites ?? 0,
      maxChiBoundaries: sourceMaxChi?.sourceChi?.boundaryStatistics ?? [0, 0],
      chimaeraBoundaries: sourceChimaera?.sourceChi?.boundaryStatistics ?? [0, 0],
      maxChiInformative: sourceMaxChi?.sourceChi?.informativeSites,
      maxChiHalfWindow: sourceMaxChi?.sourceChi?.halfWindow,
      chimaeraInformative: sourceChimaera?.sourceChi?.informativeSites,
      chimaeraHalfWindow: sourceChimaera?.sourceChi?.halfWindow,
      threeSeqMajorSites: sourceThreeSeq?.sourceThreeSeq.upSteps ?? 0,
      threeSeqMinorSites: sourceThreeSeq?.sourceThreeSeq.downSteps ?? 0,
      bootscanBootstrapConsistent: sourceBootscan
        ? Math.round(sourceBootscan.sourceBootscan.bootstrapSupport * sourceBootscan.sourceBootscan.bootstrapReplicates)
        : 0,
      bootscanBootstrapReplicates: sourceBootscan?.sourceBootscan.bootstrapReplicates ?? 0,
      threeSeqStart: sourceThreeSeq?.start ?? 0,
      threeSeqEnd: sourceThreeSeq?.end ?? 0,
      maxChiStart: 0,
      maxChiEnd: 0,
      chimaeraStart: 0,
      chimaeraEnd: 0,
      bootscanStart: sourceBootscan?.start ?? 0,
      bootscanEnd: sourceBootscan?.end ?? 0,
      siskanStart: 0,
      siskanEnd: 0,
      bootscanRunWindows: sourceBootscan?.sourceBootscan.runWindows ?? 0,
      siskanRunWindows: 0,
      rdpSource: candidate.sourceRdp,
      bootscanSource: sourceBootscan?.sourceBootscan,
    };
    if (sourceSiScanEnabled) {
      const triplet = [candidate.analysisRecombinant, candidate.analysisMajorParent, candidate.analysisMinorParent];
      const analysisTripletKey = [...triplet].sort((left, right) => left - right).join(":");
      const existingSourceSignal = (candidate.methodSignals ?? [])
        .filter((signal) => signal.method === "SiScan"
          && signal.sourceAnalysisTripletKey === analysisTripletKey
          && signal.sourceSiScan?.recombinant === candidate.recombinant
          && signal.sourceSiScan?.majorParent === candidate.majorParent
          && signal.sourceSiScan?.minorParent === candidate.minorParent
          && isCoLocatedSignal(candidate, signal, nSites))
        .sort((left, right) => (left.sourceSiScan?.rawP ?? 1) - (right.sourceSiScan?.rawP ?? 1)
          || Math.abs(right.statistic) - Math.abs(left.statistic))[0];
      candidate.methodSignals = (candidate.methodSignals ?? []).filter((signal) => signal.method !== "SiScan");
      candidate.stats.siskanScore = 0;
      candidate.stats.siskanSites = 0;
      candidate.stats.siskanSourceP = 1;
      candidate.stats.siskanSourceZ = 0;
      let sourceSignal = existingSourceSignal;
      if (!sourceSignal) {
        const source = getSourceSiScanResult(
          triplet,
          { sequencePtr: candidate.sequencePtr, rotation: candidate.rotation },
          candidate.siskanCandidatePool ?? (exactDistanceMatrix ? allIndexes : referencePool),
        );
        const selected = selectCoLocatedSiScanRegion(
          source.result,
          candidate,
          candidate.rotation,
          nSites,
          source.triplet,
        );
        if (selected) {
          sourceSignal = sourceSiScanMethodSignal(
            source.result,
            selected.region,
            source.triplet,
            candidate.rotation,
            nSites,
            scanMappings,
            unique.length > 1_000 ? 48 : unique.length > 250 ? 96 : 192,
          );
        }
      }
      if (sourceSignal) {
        candidate.methodSignals.push(sourceSignal);
        applySourceSiScanStats(candidate.stats, sourceSignal);
      }
    }
    if (options.polishBreakpoints) {
      let polishedStart = candidate.rawStart;
      let polishedEnd = candidate.rawEnd;
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
      emit({
        type: "progress",
        jobId,
        progress: 0.85 + 0.15 * (candidateIndex + 1) / Math.max(1, unique.length),
        phase: "Calibrating method evidence",
      });
    }
  }
  let statisticsMs = performance.now() - statisticsStarted;
  let events = unique.map((candidate, index) => {
    const evidence = methodEvidence(candidate, candidate.stats, options, Math.max(1, comparisons), nSites);
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

  const result = {
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
    geneconvSignalTruncations: truncatedGeneconvSignals,
    chiSignalTruncations: truncatedChiSignals,
    bootscanSignalTruncations: truncatedBootscanSignals,
    bootscanBatch: sourceBootscanEnabled ? {
      calls: sourceBootscanBatchCalls,
      triplets: sourceBootscanTripletCount,
      usedPairs: sourceBootscanUsedPairCount,
      windows: Math.floor(nSites / Math.max(1, sourceBootscanStep)) + 2,
      replicates: sourceBootscanReplicates,
      workspaceBytes: sourceBootscanWorkspaceBytes,
      relationshipMode: "distance",
    } : undefined,
    tripletKernelCalls: {
      rdp: rdpTripletKernelCalls,
      geneconv: sourceGeneconvTripletKernelCalls,
      sourceChi: sourceChiTripletKernelCalls,
      threeSeq: sourceThreeSeqTripletKernelCalls,
      siscan: sourceSiScanTripletCalls,
    },
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
      `source RDP/GENECONV${sourceBootscanEnabled ? "/BOOTSCAN" : ""}/MAXCHI/CHIMAERA${sourceSiScanEnabled ? "/SISCAN" : ""}${sourceThreeSeqEnabled ? "/3SEQ" : ""}`,
      options.circular ? "dual-origin circular scan" : "linear scan",
      options.polishBreakpoints ? (options.burtMode === "manual-step-up" ? "BURT 2–20-state step-up" : "RDP5-source BURT") : "raw breakpoints",
      options.ancestralClustering === false ? "event clustering off" : `${clustered.clusters.length} ancestral clusters`,
      disassembly.componentCount > 0 ? `${disassembly.componentCount} extracted analysis components` : "intact alignment",
    ].join(" · "),
  };
  emit(result);
  return result;
}

function cycleEventP(event) {
  return Math.min(1, ...(event.evidence ?? []).map((item) => item.correctedP));
}

function cycleEventOrigins(event) {
  return [event.recombinant, event.majorParent, event.minorParent];
}

function cycleRecombinantOrigins(event) {
  const selected = event.coRecombinantSets?.find((set) => set.presumedRecombinant === event.recombinant);
  const members = selected?.sequenceMembers?.length
    ? selected.sequenceMembers
    : event.ancestralCluster?.sequenceMembers?.length
      ? event.ancestralCluster.sequenceMembers
      : [event.recombinant];
  return [...new Set([event.recombinant, ...members])].sort((left, right) => left - right);
}

function sameCycleEvent(left, right, length) {
  if (left.recombinant !== right.recombinant) return false;
  const leftParents = [left.majorParent, left.minorParent].sort((a, b) => a - b).join(":");
  const rightParents = [right.majorParent, right.minorParent].sort((a, b) => a - b).join(":");
  if (leftParents !== rightParents) return false;
  return sameCircularBreakpoints(left, right, length) || overlap(left, right, length) > 0.75;
}

function mergeCyclePool(current, incoming, selected, invalidatedOrigins, length) {
  const retained = current.filter((event) => (
    !cycleEventOrigins(event).some((origin) => invalidatedOrigins.has(origin))
    && !selected.some((accepted) => sameCycleEvent(accepted, event, length))
  ));
  for (const event of incoming) {
    if (selected.some((accepted) => sameCycleEvent(accepted, event, length))) continue;
    const duplicate = retained.find((candidate) => sameCycleEvent(candidate, event, length));
    if (!duplicate) retained.push(event);
    else if (cycleEventP(event) < cycleEventP(duplicate)) retained[retained.indexOf(duplicate)] = event;
  }
  return retained.sort((left, right) => cycleEventP(left) - cycleEventP(right));
}

// RDP5 manual section 4.1.6: screen, take the best signal, identify its
// recombinant/co-recombinant set, erase and extract those tracts, then repeat
// until no signal remains. Unaffected detections stay in the pool; only
// triplets containing an origin modified by the latest event are put on the
// redo list. This is the source workflow's sequential decomposition rather
// than a post-hoc reconstruction of a one-pass result list.
async function analyzeCyclic(message) {
  const started = performance.now();
  const selected = [];
  const applied = Array.isArray(message.disassemblyEvents)
    ? message.disassemblyEvents.filter((event) => event?.decision === "accepted" && event.evidenceStale !== true)
    : [];
  const preexistingAppliedCount = applied.length;
  let pool = [];
  let affectedOrigins = [];
  let executedPasses = 0;
  let lastResult = null;
  let stoppedBecause = "no-detectable-signals";
  const maximumCycles = Math.max(1, Math.min(1000, Math.trunc(message.options.maximumDetectionCycles ?? 250)));
  const aggregate = {
    comparisons: 0,
    initialComparisons: 0,
    distanceMs: 0,
    scanMs: 0,
    statisticsMs: 0,
    diagnosticsMs: 0,
    clusteringMs: 0,
    rdpSignalTruncations: 0,
    geneconvSignalTruncations: 0,
    chiSignalTruncations: 0,
    bootscanSignalTruncations: 0,
    bootscanBatchCalls: 0,
    bootscanTriplets: 0,
    bootscanUsedPairs: 0,
    bootscanWorkspaceBytes: 0,
    rdpCalls: 0,
    geneconvCalls: 0,
    sourceChiCalls: 0,
    threeSeqCalls: 0,
    siscanCalls: 0,
  };

  while (selected.length < maximumCycles && executedPasses < maximumCycles * 4 + 4) {
    let passResult = null;
    const passNumber = executedPasses + 1;
    const passJobId = `${message.jobId}-cycle-${passNumber}`;
    await analyze({
      ...message,
      jobId: passJobId,
      disassemblyEvents: applied,
      affectedOrigins: executedPasses === 0 ? [] : affectedOrigins,
    }, (payload) => {
      if (payload.type === "result") {
        passResult = payload;
        return;
      }
      if (payload.type === "progress") {
        postMessage({
          ...payload,
          jobId: message.jobId,
          progress: Math.min(0.98, 1 - 1 / (passNumber + Math.max(0, Math.min(1, payload.progress)))),
          phase: `Detection pass ${passNumber} · ${payload.phase}`,
        });
        return;
      }
      if (payload.type === "partial") {
        postMessage({ ...payload, jobId: message.jobId, events: [...selected, ...payload.events] });
      }
    });
    if (!passResult) throw new Error(`Detection pass ${passNumber} did not return a result.`);
    lastResult = passResult;
    executedPasses += 1;
    aggregate.comparisons += passResult.comparisons;
    if (executedPasses === 1) aggregate.initialComparisons = passResult.comparisons;
    aggregate.distanceMs += passResult.timing?.distanceMs ?? 0;
    aggregate.scanMs += passResult.timing?.scanMs ?? 0;
    aggregate.statisticsMs += passResult.timing?.statisticsMs ?? 0;
    aggregate.diagnosticsMs += passResult.timing?.diagnosticsMs ?? 0;
    aggregate.clusteringMs += passResult.timing?.clusteringMs ?? 0;
    aggregate.rdpSignalTruncations += passResult.rdpSignalTruncations ?? 0;
    aggregate.geneconvSignalTruncations += passResult.geneconvSignalTruncations ?? 0;
    aggregate.chiSignalTruncations += passResult.chiSignalTruncations ?? 0;
    aggregate.bootscanSignalTruncations += passResult.bootscanSignalTruncations ?? 0;
    aggregate.bootscanBatchCalls += passResult.bootscanBatch?.calls ?? 0;
    aggregate.bootscanTriplets += passResult.bootscanBatch?.triplets ?? 0;
    aggregate.bootscanUsedPairs += passResult.bootscanBatch?.usedPairs ?? 0;
    aggregate.bootscanWorkspaceBytes = Math.max(aggregate.bootscanWorkspaceBytes, passResult.bootscanBatch?.workspaceBytes ?? 0);
    aggregate.rdpCalls += passResult.tripletKernelCalls?.rdp ?? 0;
    aggregate.geneconvCalls += passResult.tripletKernelCalls?.geneconv ?? 0;
    aggregate.sourceChiCalls += passResult.tripletKernelCalls?.sourceChi ?? 0;
    aggregate.threeSeqCalls += passResult.tripletKernelCalls?.threeSeq ?? 0;
    aggregate.siscanCalls += passResult.tripletKernelCalls?.siscan ?? 0;

    pool = mergeCyclePool(
      pool,
      passResult.events,
      selected,
      new Set(executedPasses === 1 ? [] : affectedOrigins),
      message.alignment.length,
    );
    const best = pool[0];
    if (!best) break;
    // Detection signals that survived from an earlier pool are still valid,
    // but RDP5 characterizes the chosen signal against the *current*
    // component alignment. Refresh every concrete triplet involving that
    // candidate before applying it if its role/group ledger predates an erase.
    const appliedIds = applied.map((event) => event.id);
    const characterizedIds = new Set(best.componentProvenance?.appliedEventIds ?? []);
    if (appliedIds.some((id) => !characterizedIds.has(id))) {
      affectedOrigins = cycleEventOrigins(best);
      continue;
    }
    pool.shift();
    const round = selected.length + 1;
    const members = cycleRecombinantOrigins(best);
    const cycleSummary = `Selected as the strongest remaining signal in RDP5 erase/extract cycle ${round}; ${members.length} recombinant lineage${members.length === 1 ? "" : "s"} split before the redo-list scan.`;
    const recorded = {
      ...best,
      decision: "unreviewed",
      note: [best.note, cycleSummary].filter(Boolean).join(" "),
      history: [
        ...(best.history ?? []),
        {
          id: `history-${message.jobId}-cycle-${round}`,
          timestamp: new Date().toISOString(),
          action: `Cyclic detection round ${round}`,
          summary: cycleSummary,
        },
      ],
    };
    selected.push(recorded);
    applied.push({ ...recorded, decision: "accepted", evidenceStale: false });
    affectedOrigins = members;
  }

  if (!lastResult) throw new Error("Cyclic detection did not start.");
  if (selected.length >= maximumCycles || executedPasses >= maximumCycles * 4 + 4) {
    stoppedBecause = "cycle-cap";
    const final = selected[selected.length - 1];
    const warning = `Cyclical detection stopped at the configured ${maximumCycles}-event safety cap while detectable signals remained.`;
    if (!final.warnings.includes(warning)) final.warnings.push(warning);
    final.note = [final.note, warning].filter(Boolean).join(" ");
  }
  const totalMs = performance.now() - started;
  postMessage({
    ...lastResult,
    jobId: message.jobId,
    events: selected,
    comparisons: aggregate.comparisons,
    elapsedMs: totalMs,
    timing: {
      distanceMs: aggregate.distanceMs,
      scanMs: aggregate.scanMs,
      statisticsMs: aggregate.statisticsMs,
      diagnosticsMs: aggregate.diagnosticsMs,
      clusteringMs: aggregate.clusteringMs,
    },
    rdpSignalTruncations: aggregate.rdpSignalTruncations,
    geneconvSignalTruncations: aggregate.geneconvSignalTruncations,
    chiSignalTruncations: aggregate.chiSignalTruncations,
    bootscanSignalTruncations: aggregate.bootscanSignalTruncations,
    bootscanBatch: lastResult.bootscanBatch ? {
      ...lastResult.bootscanBatch,
      calls: aggregate.bootscanBatchCalls,
      triplets: aggregate.bootscanTriplets,
      usedPairs: aggregate.bootscanUsedPairs,
      workspaceBytes: aggregate.bootscanWorkspaceBytes,
    } : undefined,
    tripletKernelCalls: {
      rdp: aggregate.rdpCalls,
      geneconv: aggregate.geneconvCalls,
      sourceChi: aggregate.sourceChiCalls,
      threeSeq: aggregate.threeSeqCalls,
      siscan: aggregate.siscanCalls,
    },
    detectionCycle: {
      enabled: true,
      eventsApplied: selected.length,
      passes: executedPasses,
      initialComparisons: aggregate.initialComparisons,
      redoComparisons: Math.max(0, aggregate.comparisons - aggregate.initialComparisons),
      stoppedBecause,
      maximumCycles,
    },
    disassembly: {
      ...(lastResult.disassembly ?? {}),
      appliedEvents: applied.length,
    },
    engine: `${lastResult.engine} · sequential erase/extract redo queue (${selected.length} new events; ${preexistingAppliedCount} pre-applied)`,
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
  const chiSignalCapacity = Math.max(1, Math.min(256, Math.trunc(options.chiSignalsPerTriplet ?? 24)));
  const geneconvSignalCapacity = Math.max(1, Math.min(256, Math.trunc(options.geneconvSignalsPerTriplet ?? 64)));
  const chiTrackCount = 3 * (
    Number(options.methods.includes("MaxChi")) + Number(options.methods.includes("Chimaera"))
  );
  const chiPeakCapacity = Math.max(8, Math.min(256, Math.ceil(
    (chiSignalCapacity * 2) / Math.max(3, chiTrackCount),
  )));
  const chiPositionsPtr = align(outPtr + 96, 16);
  const chiScoresPtr = chiPositionsPtr + 4 * (nSites + 1) * 4;
  const chiMissingPrefixPtr = align(chiScoresPtr + 4 * (nSites + 1), 4);
  const chiProfilePtr = align(chiMissingPrefixPtr + (nSites + 1) * 4, 8);
  const chiSmoothPtr = chiProfilePtr + (nSites + 1) * 8;
  const chiPeakPtr = align(chiSmoothPtr + (nSites + 1) * 8, 4);
  const chiOutPtr = chiPeakPtr + chiPeakCapacity * SOURCE_CHI_PEAK_INTS * 4;
  const geneconvPositionsPtr = align(chiOutPtr + chiSignalCapacity * SOURCE_CHI_ROW_INTS * 4, 16);
  const geneconvCategoriesPtr = geneconvPositionsPtr + nSites * 4;
  const geneconvRunStartPtr = align(geneconvCategoriesPtr + nSites, 4);
  const geneconvRunEndPtr = geneconvRunStartPtr + nSites * 4;
  const geneconvRunScorePtr = geneconvRunEndPtr + nSites * 4;
  const geneconvPrefixPtr = align(geneconvRunScorePtr + nSites * 4, 8);
  const geneconvTreePtr = geneconvPrefixPtr + (nSites + 1) * 8;
  const geneconvWorkspaceBytes = (nSites + 1) * 16 + 8;
  const geneconvCalibrationPtr = align(geneconvTreePtr + geneconvWorkspaceBytes, 8);
  const geneconvCandidatePtr = align(geneconvCalibrationPtr + 6 * 40, 8);
  const geneconvCandidateCapacity = 3 * (nSites + 1);
  const geneconvDeletePtr = align(geneconvCandidatePtr + geneconvCandidateCapacity * 24, 4);
  const geneconvOutPtr = align(geneconvDeletePtr + nSites * 4, 8);
  const workingPacked = packSequences(working, workingSequenceCount, nSites);
  const packedPtr = align(geneconvOutPtr + geneconvSignalCapacity * SOURCE_GENECONV_ROW_INTS * 4, 16);
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
  const recalcBootscanEnabled = options.methods.includes("BootScan");
  const recalcBootscanWindow = Math.max(5, Math.min(32_767, Math.min(Math.floor(nSites / 2), Math.trunc(options.bootscanWindow ?? 200))));
  const recalcBootscanStep = Math.max(1, Math.min(Math.max(1, Math.floor(nSites / 4)), Math.trunc(options.bootscanStep ?? 20)));
  const recalcBootscanReplicates = Math.max(0, Math.min(1000, Math.trunc(options.bootstrapReplicates ?? 100)));
  if (recalcBootscanEnabled && recalcBootscanReplicates < 2) {
    throw new Error("BootScan requires at least two replicates. Increase Bootstrap replicates or disable BootScan.");
  }
  const recalcBootscanPairCount = workingSequenceCount * (workingSequenceCount - 1) / 2;
  const recalcThreeSeqOutPtr = align(roleDmaxOutPtr + 40, 16);
  const recalcThreeSeqWorkspacePtr = align(recalcThreeSeqOutPtr + 6 * SOURCE_THREE_SEQ_ROW_INTS * 4, 16);
  const recalcThreeSeqWorkspaceBytes = typeof instance.exports.source_three_seq_workspace_bytes === "function"
    ? instance.exports.source_three_seq_workspace_bytes(nSites)
    : 0;
  const recalcBootscanTripletPtr = align(recalcThreeSeqWorkspacePtr + recalcThreeSeqWorkspaceBytes, 16);
  const recalcBootscanPairMapPtr = recalcBootscanTripletPtr + 12;
  const recalcBootscanPairListPtr = recalcBootscanPairMapPtr + recalcBootscanPairCount * 4;
  const recalcBootscanWeightPtr = align(recalcBootscanPairListPtr + 3 * 8, 2);
  const recalcBootscanPairDistancePtr = recalcBootscanWeightPtr + recalcBootscanWindow * Math.max(2, recalcBootscanReplicates) * 2;
  const recalcBootscanGlobalPairPtr = recalcBootscanPairDistancePtr + 3 * Math.max(2, recalcBootscanReplicates) * 2;
  const recalcBootscanStatePtr = align(recalcBootscanGlobalPairPtr + 3 * 8, 4);
  const recalcBootscanDifferencePtr = recalcBootscanStatePtr + 24;
  const recalcBootscanValidPtr = recalcBootscanDifferencePtr + Math.max(2, recalcBootscanReplicates) * 4;
  const recalcBootscanLookupPtr = align(recalcBootscanValidPtr + Math.max(2, recalcBootscanReplicates) * 4, 2);
  // A noisy or long triplet can contain many source runs.  Keep enough rows to
  // find the event being recalculated instead of silently treating a truncated
  // prefix as the complete RDP5 signal ledger.
  const recalcBootscanOutCapacity = 4_096;
  const recalcBootscanLookupEntries = (recalcBootscanWindow + 1) * (recalcBootscanWindow + 2) / 2;
  const recalcBootscanOutPtr = align(recalcBootscanLookupPtr + recalcBootscanLookupEntries * 2, 4);
  const recalcRequiredBytes = recalcBootscanEnabled
    ? recalcBootscanOutPtr + recalcBootscanOutCapacity * SOURCE_BOOTSCAN_ROW_INTS * 4
    : recalcThreeSeqWorkspacePtr + recalcThreeSeqWorkspaceBytes;
  const requiredPages = Math.ceil(recalcRequiredBytes / 65536);
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
  const stats = {
    genconvRun: 0,
    genconvEligible: 0,
    genconvMatches: 0,
    genconvStart: 0,
    genconvEnd: 0,
    bootscanConsistent: 0,
    bootscanWindows: 0,
    maxChi: 0,
    chimaera: 0,
    siskanScore: 0,
    siskanSites: 0,
    threeSeqDescent: 0,
    threeSeqSites: 0,
    maxChiBoundaries: [0, 0],
    chimaeraBoundaries: [0, 0],
    threeSeqMajorSites: 0,
    threeSeqMinorSites: 0,
    bootscanBootstrapConsistent: 0,
    bootscanBootstrapReplicates: 0,
    threeSeqStart: 0,
    threeSeqEnd: 0,
    maxChiStart: 0,
    maxChiEnd: 0,
    chimaeraStart: 0,
    chimaeraEnd: 0,
    bootscanStart: 0,
    bootscanEnd: 0,
    siskanStart: 0,
    siskanEnd: 0,
    bootscanRunWindows: 0,
    siskanRunWindows: 0,
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
    structuralUncertainty: event.structuralUncertainty,
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
  let recalculatedGeneconvSignal = null;
  if (options.methods.includes("GENECONV")) {
    const orderedTriplet = [analysisRecombinant, analysisMajorParent, analysisMinorParent]
      .sort((left, right) => left - right);
    const targetSlot = orderedTriplet.indexOf(analysisRecombinant);
    const total = instance.exports.scan_source_geneconv_all_packed(
      packedPtr,
      validityPtr,
      workingPacked.wordsPerSequence,
      nSites,
      orderedTriplet[0],
      orderedTriplet[1],
      orderedTriplet[2],
      Math.max(0, Math.round(options.geneconvGScale ?? 1)),
      Math.max(Number.MIN_VALUE, Math.min(1, options.alpha ?? 0.05)),
      geneconvPositionsPtr,
      geneconvCategoriesPtr,
      geneconvRunStartPtr,
      geneconvRunEndPtr,
      geneconvRunScorePtr,
      geneconvPrefixPtr,
      geneconvTreePtr,
      geneconvCalibrationPtr,
      geneconvCandidatePtr,
      geneconvCandidateCapacity,
      geneconvDeletePtr,
      geneconvOutPtr,
      geneconvSignalCapacity,
    );
    for (let index = 0; index < Math.min(total, geneconvSignalCapacity); index += 1) {
      const row = Array.from(new Int32Array(
        memory.buffer,
        geneconvOutPtr + index * SOURCE_GENECONV_ROW_INTS * 4,
        SOURCE_GENECONV_ROW_INTS,
      ));
      if (row[1] !== targetSlot) continue;
      const signal = sourceGeneconvSignal(row, rotation, nSites);
      if (!isCoLocatedSignal(candidate, signal, nSites)) continue;
      if (!recalculatedGeneconvSignal
        || signal.sourceGeneconv.rawP < recalculatedGeneconvSignal.sourceGeneconv.rawP) {
        recalculatedGeneconvSignal = signal;
      }
    }
    if (recalculatedGeneconvSignal) {
      stats.genconvRun = recalculatedGeneconvSignal.sourceGeneconv.fragmentScore;
      stats.genconvEligible = recalculatedGeneconvSignal.sourceGeneconv.informativeSites;
      stats.genconvMatches = recalculatedGeneconvSignal.sourceGeneconv.matchingSites;
      stats.genconvStart = recalculatedGeneconvSignal.start;
      stats.genconvEnd = recalculatedGeneconvSignal.end;
    }
  }
  let recalculatedBootscanSignal = null;
  if (recalcBootscanEnabled) {
    const orderedTriplet = [analysisRecombinant, analysisMajorParent, analysisMinorParent]
      .sort((left, right) => left - right);
    new Int32Array(memory.buffer, recalcBootscanTripletPtr, 3).set(orderedTriplet);
    const total = instance.exports.scan_source_bootscan_batch_packed(
      packedPtr,
      validityPtr,
      workingPacked.wordsPerSequence,
      workingSequenceCount,
      nSites,
      recalcBootscanTripletPtr,
      1,
      recalcBootscanWindow,
      recalcBootscanStep,
      recalcBootscanReplicates,
      Math.round(Math.max(0.5, Math.min(0.999, options.bootscanCutoff ?? 0.7)) * 1000),
      (options.randomSeed ?? 0x5a17c0de) >>> 0,
      recalcBootscanPairMapPtr,
      recalcBootscanPairListPtr,
      recalcBootscanWeightPtr,
      recalcBootscanPairDistancePtr,
      recalcBootscanGlobalPairPtr,
      recalcBootscanStatePtr,
      recalcBootscanDifferencePtr,
      recalcBootscanValidPtr,
      recalcBootscanLookupPtr,
      recalcBootscanOutPtr,
      recalcBootscanOutCapacity,
    );
    for (let index = 0; index < Math.min(total, recalcBootscanOutCapacity); index += 1) {
      const row = Array.from(new Int32Array(
        memory.buffer,
        recalcBootscanOutPtr + index * SOURCE_BOOTSCAN_ROW_INTS * 4,
        SOURCE_BOOTSCAN_ROW_INTS,
      ));
      const signal = sourceBootscanSignal(row, rotation, nSites);
      if (!isCoLocatedSignal(candidate, signal, nSites)) continue;
      if (!recalculatedBootscanSignal
        || signal.sourceBootscan.rawP < recalculatedBootscanSignal.sourceBootscan.rawP) {
        recalculatedBootscanSignal = signal;
      }
    }
    if (recalculatedBootscanSignal) {
      stats.bootscanConsistent = recalculatedBootscanSignal.sourceBootscan.runWindows;
      stats.bootscanWindows = recalculatedBootscanSignal.sourceBootscan.runWindows;
      stats.bootscanBootstrapConsistent = Math.round(
        recalculatedBootscanSignal.sourceBootscan.bootstrapSupport
          * recalculatedBootscanSignal.sourceBootscan.bootstrapReplicates,
      );
      stats.bootscanBootstrapReplicates = recalculatedBootscanSignal.sourceBootscan.bootstrapReplicates;
      stats.bootscanStart = recalculatedBootscanSignal.start;
      stats.bootscanEnd = recalculatedBootscanSignal.end;
      stats.bootscanRunWindows = recalculatedBootscanSignal.sourceBootscan.runWindows;
      stats.bootscanSource = recalculatedBootscanSignal.sourceBootscan;
    }
  }
  const recalculatedSourceChiSignals = [];
  const recalculatedChiMask = (options.methods.includes("MaxChi") ? 4 : 0)
    | (options.methods.includes("Chimaera") ? 8 : 0);
  if (recalculatedChiMask !== 0) {
    const orderedTriplet = [analysisRecombinant, analysisMajorParent, analysisMinorParent]
      .sort((left, right) => left - right);
    const targetSlot = orderedTriplet.indexOf(analysisRecombinant);
    const total = instance.exports.scan_source_chi_all(
      seqPtr,
      nSites,
      orderedTriplet[0],
      orderedTriplet[1],
      orderedTriplet[2],
      Math.max(20, options.window),
      0,
      recalculatedChiMask,
      chiPositionsPtr,
      chiScoresPtr,
      chiMissingPrefixPtr,
      chiProfilePtr,
      chiSmoothPtr,
      chiPeakPtr,
      chiPeakCapacity,
      chiOutPtr,
      chiSignalCapacity,
    );
    const bestByMethod = new Map();
    for (let index = 0; index < Math.min(total, chiSignalCapacity); index += 1) {
      const row = Array.from(new Int32Array(
        memory.buffer,
        chiOutPtr + index * SOURCE_CHI_ROW_INTS * 4,
        SOURCE_CHI_ROW_INTS,
      ));
      if (row[0] === 4 && row[1] !== targetSlot) continue;
      const signal = sourceChiSignal(row, rotation, nSites);
      if (!isCoLocatedSignal(candidate, signal, nSites)) continue;
      const previous = bestByMethod.get(signal.method);
      if (!previous || signal.statistic > previous.statistic) bestByMethod.set(signal.method, signal);
    }
    recalculatedSourceChiSignals.push(...bestByMethod.values());
    const maxChi = bestByMethod.get("MaxChi");
    const chimaera = bestByMethod.get("Chimaera");
    if (maxChi) {
      stats.maxChi = maxChi.statistic;
      stats.maxChiBoundaries = maxChi.sourceChi.boundaryStatistics;
      stats.maxChiInformative = maxChi.sourceChi.informativeSites;
      stats.maxChiHalfWindow = maxChi.sourceChi.halfWindow;
    } else {
      stats.maxChi = 0;
      stats.maxChiBoundaries = [0, 0];
    }
    if (chimaera) {
      stats.chimaera = chimaera.statistic;
      stats.chimaeraBoundaries = chimaera.sourceChi.boundaryStatistics;
      stats.chimaeraInformative = chimaera.sourceChi.informativeSites;
      stats.chimaeraHalfWindow = chimaera.sourceChi.halfWindow;
    } else {
      stats.chimaera = 0;
      stats.chimaeraBoundaries = [0, 0];
    }
  }
  let recalculatedThreeSeqSignal = null;
  if (options.methods.includes("3Seq")) {
    const orderedTriplet = [analysisRecombinant, analysisMajorParent, analysisMinorParent]
      .sort((left, right) => left - right);
    const total = typeof instance.exports.scan_source_three_seq_triplet_mode === "function"
      ? instance.exports.scan_source_three_seq_triplet_mode(
          seqPtr,
          nSites,
          orderedTriplet[0],
          orderedTriplet[1],
          orderedTriplet[2],
          options.circular ? 1 : 0,
          recalcThreeSeqWorkspacePtr,
          recalcThreeSeqOutPtr,
        )
      : instance.exports.scan_source_three_seq_triplet_packed(
          packedPtr,
          validityPtr,
          workingPacked.wordsPerSequence,
          nSites,
          orderedTriplet[0],
          orderedTriplet[1],
          orderedTriplet[2],
          recalcThreeSeqOutPtr,
        );
    const rows = Array.from({ length: Math.min(6, Math.max(0, total)) }, (_, index) => (
      Array.from(new Int32Array(
        memory.buffer,
        recalcThreeSeqOutPtr + index * SOURCE_THREE_SEQ_ROW_INTS * 4,
        SOURCE_THREE_SEQ_ROW_INTS,
      ))
    ));
    for (const fullSignal of selectSourceThreeSeqSignals(
      rows,
      rotation,
      nSites,
      Math.max(Number.MIN_VALUE, options.alpha ?? 0.05),
      Math.max(10_000, Math.trunc(options.threeSeqExactOperations ?? 1_000_000)),
    )) {
      const signal = refineSourceThreeSeqSignalForPiece(
        fullSignal,
        candidate,
        working,
        nSites,
        Math.max(10_000, Math.trunc(options.threeSeqExactOperations ?? 1_000_000)),
      );
      if (!signal) continue;
      if (signal.sourceThreeSeq.target !== analysisRecombinant) continue;
      if (!isCoLocatedSignal(candidate, signal, nSites)) continue;
      if (!recalculatedThreeSeqSignal
        || signal.sourceThreeSeq.rawP < recalculatedThreeSeqSignal.sourceThreeSeq.rawP) {
        recalculatedThreeSeqSignal = {
          ...signal,
          sourceThreeSeq: {
            ...signal.sourceThreeSeq,
            target: disassembly.mappings[signal.sourceThreeSeq.target].originIndex,
            majorParent: disassembly.mappings[signal.sourceThreeSeq.majorParent].originIndex,
            minorParent: disassembly.mappings[signal.sourceThreeSeq.minorParent].originIndex,
          },
        };
      }
    }
    if (recalculatedThreeSeqSignal) {
      stats.threeSeqDescent = recalculatedThreeSeqSignal.sourceThreeSeq.descent;
      stats.threeSeqSites = recalculatedThreeSeqSignal.sourceThreeSeq.informativeSites;
      stats.threeSeqMajorSites = recalculatedThreeSeqSignal.sourceThreeSeq.upSteps;
      stats.threeSeqMinorSites = recalculatedThreeSeqSignal.sourceThreeSeq.downSteps;
      stats.threeSeqStart = recalculatedThreeSeqSignal.start;
      stats.threeSeqEnd = recalculatedThreeSeqSignal.end;
    }
  }
  let recalculatedSiScanSignal = null;
  if (options.methods.includes("SiScan")) {
    const sourceRequested = true;
    stats.siskanScore = 0;
    stats.siskanSites = 0;
    stats.siskanSourceP = 1;
    stats.siskanSourceZ = 0;
    if (sourceRequested) {
      const triplet = [analysisRecombinant, analysisMajorParent, analysisMinorParent]
        .sort((left, right) => left - right);
      const tripletOrigins = new Set(triplet.map((index) => disassembly.mappings[index].originIndex));
      const candidatePool = Array.from({ length: workingSequenceCount }, (_, index) => index)
        .filter((index) => !tripletOrigins.has(disassembly.mappings[index].originIndex));
      const manualOutgroup = options.siskanOutgroupMode === "manual"
        ? Array.from({ length: workingSequenceCount }, (_, index) => index).find((index) => (
            disassembly.mappings[index].originIndex === options.siskanOutgroupSequence
            && !tripletOrigins.has(disassembly.mappings[index].originIndex)
          ))
        : undefined;
      const sourceSeed = (
        (options.randomSeed ?? 0x5a17c0de)
        ^ Math.imul(triplet[0] + 1, 0x9e3779b1)
        ^ Math.imul(triplet[1] + 1, 0x85ebca6b)
        ^ Math.imul(triplet[2] + 1, 0xc2b2ae35)
        ^ Math.imul(rotation + 1, 0x27d4eb2f)
      ) >>> 0;
      const sourceResult = runSourceSiScan(
        working,
        nSites,
        workingSequenceCount,
        triplet,
        {
          window: Math.max(12, options.siskanWindow ?? options.window),
          step: Math.max(1, options.siskanStep ?? options.step),
          scanPermutations: options.siskanScanPermutations ?? 100,
          pValuePermutations: options.siskanPValuePermutations ?? 1000,
          seed: sourceSeed,
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
      const selected = selectCoLocatedSiScanRegion(sourceResult, candidate, rotation, nSites, triplet);
      if (selected) {
        recalculatedSiScanSignal = sourceSiScanMethodSignal(
          sourceResult,
          selected.region,
          triplet,
          rotation,
          nSites,
          disassembly.mappings,
          192,
        );
        applySourceSiScanStats(stats, recalculatedSiScanSignal);
      }
    }
  }
  candidate.diagnostics = candidateDiagnostics(candidate, encoded, nSites, profile);
  const recalculatedMethodSignals = [
    ...(options.methods.includes("RDP") ? [{ method: "RDP", ...mappedInterval, statistic: candidate.chiSquare, locator: "edited hypothesis recalculation" }] : []),
    ...(recalculatedGeneconvSignal ? [recalculatedGeneconvSignal] : []),
    ...(recalculatedBootscanSignal ? [recalculatedBootscanSignal] : []),
    ...recalculatedSourceChiSignals,
    ...(recalculatedThreeSeqSignal ? [recalculatedThreeSeqSignal] : []),
    ...(recalculatedSiScanSignal ? [recalculatedSiScanSignal] : []),
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
  const operation = event.data.type === "recalculate"
    ? recalculate
    : event.data.cyclicDetection === true
      ? analyzeCyclic
      : analyze;
  operation(event.data).catch((error) => {
    postMessage({
      type: "error",
      jobId: event.data.jobId,
      message: error instanceof Error ? error.message : String(error),
    });
  });
};
