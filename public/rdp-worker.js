import { methodEvidence } from "./rdp-statistics.js";

let wasmPromise;

const BASES = { A: 0, C: 1, G: 2, T: 3, U: 3, "-": 5 };

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
  if (rotation === 0) return { start, end, wraps: false };
  const mappedStart = (start + rotation) % length;
  const mappedEnd = (end + rotation) % length;
  if (mappedEnd === 0) return { start: mappedStart, end: length, wraps: false };
  return { start: mappedStart, end: mappedEnd, wraps: mappedStart > mappedEnd };
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

function addWarnings(candidate, sequences, length, options) {
  const warnings = [];
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
function alignmentDiagnostics(encoded, sequenceCount, length, window) {
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
  const summary = {
    sampledSequences: sampleCount,
    sampledBiallelicSites: sampledSites.length,
    testedSitePairs: testedPairs,
    incompatibleSitePairs: incompatible,
    fourGameteFraction: incompatible / Math.max(1, testedPairs),
    nearIncompatibility: nearFraction,
    farIncompatibility: farFraction,
    proximityRatio,
    ambiguityFraction: nonCanonical / Math.max(1, sampleCount * length),
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

function deduplicate(candidates, length) {
  const ordered = [...candidates].sort((left, right) =>
    right.chiSquare - left.chiSquare || tractLength(left, length) - tractLength(right, length),
  );
  const kept = [];
  for (const candidate of ordered) {
    const duplicate = kept.find((existing) =>
      existing.recombinant === candidate.recombinant &&
      (overlap(existing, candidate, length) > 0.62 || sameCircularBreakpoints(existing, candidate, length)),
    );
    if (!duplicate) kept.push(candidate);
    else if (!duplicate.alternatives.includes(candidate.minorParent)) duplicate.alternatives.push(candidate.minorParent);
  }
  return kept.slice(0, 500);
}

async function analyze(message) {
  const started = performance.now();
  const instance = await loadWasm();
  const { alignment, options, jobId } = message;
  const sequences = alignment.sequences;
  const nSeq = sequences.length;
  const nSites = alignment.length;
  const encoded = encodeSequences(sequences, nSites);
  const diagnosticsStarted = performance.now();
  const diagnosticProfile = alignmentDiagnostics(encoded, nSeq, nSites, options.window);
  const diagnosticsMs = performance.now() - diagnosticsStarted;
  const rotation = options.circular ? Math.floor(nSites / 2) : 0;
  const rotated = rotation > 0 ? rotateSequences(encoded, nSeq, nSites, rotation) : null;
  const { packed, validity, wordsPerSequence } = packSequences(encoded, nSeq, nSites);
  const exactDistanceMatrix = nSeq <= 512
    && nSeq * nSeq * wordsPerSequence <= 50_000_000;
  const parentSamples = nSeq > 2_000 ? 64 : nSeq > 1_000 ? 128 : 256;
  const matrixCount = exactDistanceMatrix ? nSeq : Math.min(24, nSeq);
  const seqPtr = 65536;
  const rotatedSeqPtr = rotated ? align(seqPtr + encoded.byteLength, 16) : seqPtr;
  const sequenceEndPtr = rotated ? rotatedSeqPtr + rotated.byteLength : seqPtr + encoded.byteLength;
  const packedPtr = align(sequenceEndPtr, 16);
  const validityPtr = align(packedPtr + packed.byteLength, 16);
  const distancePtr = align(validityPtr + validity.byteLength, 16);
  const distanceBytes = matrixCount * matrixCount * 4;
  const prefixAPtr = align(distancePtr + distanceBytes, 16);
  const prefixBPtr = prefixAPtr + (nSites + 1) * 4;
  const outPtr = align(prefixBPtr + (nSites + 1) * 4, 16);
  const statsPtr = align(outPtr + 64, 16);
  const poolPtr = align(statsPtr + 96, 16);
  const nearestIndexesPtr = poolPtr + nSeq * 4;
  const nearestDistancesPtr = nearestIndexesPtr + nSeq * 4;
  const requiredBytes = nearestDistancesPtr + nSeq * 4;
  const memory = instance.exports.memory;
  const requiredPages = Math.ceil(requiredBytes / 65536);
  const currentPages = memory.buffer.byteLength / 65536;
  if (requiredPages > currentPages) memory.grow(requiredPages - currentPages);
  new Uint8Array(memory.buffer, seqPtr, encoded.byteLength).set(encoded);
  if (rotated) new Uint8Array(memory.buffer, rotatedSeqPtr, rotated.byteLength).set(rotated);
  new Uint32Array(memory.buffer, packedPtr, packed.length).set(packed);
  new Uint32Array(memory.buffer, validityPtr, validity.length).set(validity);

  const distanceStarted = performance.now();
  instance.exports.distance_matrix_packed(
    packedPtr,
    validityPtr,
    matrixCount,
    wordsPerSequence,
    distancePtr,
  );
  const distanceMs = performance.now() - distanceStarted;
  const distance = new Float32Array(memory.buffer, distancePtr, matrixCount * matrixCount).slice();
  const allIndexes = Array.from({ length: nSeq }, (_, index) => index);
  const targets = options.mode === "query-reference"
    ? allIndexes.filter((index) => sequences[index].role === "query" || sequences[index].role === "both")
    : allIndexes;
  const referencePool = options.mode === "query-reference"
    ? allIndexes.filter((index) => sequences[index].role === "reference" || sequences[index].role === "both")
    : allIndexes;
  new Int32Array(memory.buffer, poolPtr, referencePool.length).set(referencePool);
  const candidates = [];
  let comparisons = 0;
  const scanStarted = performance.now();
  const scanViews = rotated
    ? [{ sequencePtr: seqPtr, rotation: 0 }, { sequencePtr: rotatedSeqPtr, rotation }]
    : [{ sequencePtr: seqPtr, rotation: 0 }];

  for (let targetPosition = 0; targetPosition < targets.length; targetPosition += 1) {
    const recombinant = targets[targetPosition];
    const parentLimit = options.exhaustive ? referencePool.length : Math.min(options.candidateParents, referencePool.length);
    let parents;
    if (options.exhaustive) {
      parents = referencePool.filter((index) => index !== recombinant);
    } else if (exactDistanceMatrix) {
      parents = candidateParents(distance, nSeq, recombinant, referencePool, parentLimit);
    } else {
      const nearestLimit = Math.max(2, Math.ceil(parentLimit * 0.625));
      const foundParents = instance.exports.nearest_candidates_sampled(
        seqPtr,
        nSites,
        recombinant,
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
        if (candidate !== recombinant && !parents.includes(candidate)) parents.push(candidate);
      }
      for (const candidate of referencePool) {
        if (parents.length >= parentLimit) break;
        if (candidate !== recombinant && !parents.includes(candidate)) parents.push(candidate);
      }
    }
    for (let left = 0; left < parents.length; left += 1) {
      for (let right = left + 1; right < parents.length; right += 1) {
        for (const view of scanViews) {
          comparisons += 1;
          const found = instance.exports.scan_pair(
            view.sequencePtr,
            nSites,
            recombinant,
            parents[left],
            parents[right],
            Math.max(20, Math.floor(options.window / 2)),
            prefixAPtr,
            prefixBPtr,
            outPtr,
          );
          if (!found) continue;
          const output = new Int32Array(memory.buffer, outPtr, 12);
          const rawStart = output[0];
          const rawEnd = output[1];
          const mapped = mapInterval(rawStart, rawEnd, view.rotation, nSites);
          candidates.push({
            recombinant,
            ...mapped,
            rawStart,
            rawEnd,
            sequencePtr: view.sequencePtr,
            rotation: view.rotation,
            majorParent: output[2],
            minorParent: output[3],
            chiSquare: output[4] / 1000,
            informative: output[5],
            insideMinor: output[6],
            insideMajor: output[7],
            outsideMajor: output[8],
            outsideMinor: output[9],
            effect: output[10] / 1e6,
            alternatives: [],
          });
        }
      }
    }
    postMessage({
      type: "progress",
      jobId,
      progress: 0.85 * (targetPosition + 1) / targets.length,
      phase: `Scanning ${sequences[recombinant].name}`,
    });
  }
  const scanMs = performance.now() - scanStarted;

  const unique = deduplicate(candidates, nSites);
  const statisticsStarted = performance.now();
  unique.forEach((candidate, candidateIndex) => {
    instance.exports.method_stats(
      candidate.sequencePtr,
      nSites,
      candidate.recombinant,
      candidate.majorParent,
      candidate.minorParent,
      candidate.rawStart,
      candidate.rawEnd,
      Math.max(20, options.window),
      Math.max(1, options.step),
      Math.max(0, options.bootstrapReplicates ?? 100),
      (
        (options.randomSeed ?? 0x5a17c0de)
        ^ Math.imul(candidate.recombinant + 1, 0x9e3779b1)
        ^ Math.imul(candidate.majorParent + 1, 0x85ebca6b)
        ^ Math.imul(candidate.minorParent + 1, 0xc2b2ae35)
        ^ candidate.rawStart
        ^ Math.imul(candidate.rawEnd, 31)
      ) | 0,
      prefixAPtr,
      prefixBPtr,
      statsPtr,
    );
    const methodOutput = new Int32Array(memory.buffer, statsPtr, 23);
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
    };
    if (options.polishBreakpoints) {
      let polishedStart = methodOutput[17];
      let polishedEnd = methodOutput[18];
      candidate.breakpointModel = {
        method: "local-chi-square",
        informativeSites: candidate.informative,
      };
      const hmmFound = instance.exports.hmm_polish(
        candidate.sequencePtr,
        nSites,
        candidate.recombinant,
        candidate.majorParent,
        candidate.minorParent,
        candidate.rawStart,
        candidate.rawEnd,
        prefixAPtr,
        prefixBPtr,
        outPtr,
      );
      if (hmmFound) {
        const hmmOutput = new Int32Array(memory.buffer, outPtr, 11);
        polishedStart = hmmOutput[0];
        polishedEnd = hmmOutput[1];
        candidate.breakpointModel = {
          method: "two-state-hmm",
          informativeSites: hmmOutput[2],
          stateSwitches: hmmOutput[3],
          majorFit: hmmOutput[4] / 1000,
          minorFit: hmmOutput[5] / 1000,
        };
      }
      if (polishedStart >= 0 && polishedEnd <= nSites && polishedEnd - polishedStart >= 12) {
        Object.assign(candidate, mapInterval(polishedStart, polishedEnd, candidate.rotation, nSites));
      }
    }
    candidate.diagnostics = candidateDiagnostics(candidate, encoded, nSites, diagnosticProfile);
    if (candidateIndex % 16 === 0 || candidateIndex === unique.length - 1) {
      postMessage({
        type: "progress",
        jobId,
        progress: 0.85 + 0.15 * (candidateIndex + 1) / Math.max(1, unique.length),
        phase: "Calibrating independent method evidence",
      });
    }
  });
  const statisticsMs = performance.now() - statisticsStarted;
  let exactThreeSeqBudget = 20_000_000;
  let events = unique.map((candidate, index) => {
    const threeSeqOperations = (candidate.stats.threeSeqMajorSites + 1)
      * (candidate.stats.threeSeqMinorSites + 1)
      * Math.max(1, candidate.stats.threeSeqDescent);
    const exactOperations = Math.min(4_000_000, exactThreeSeqBudget);
    const evidence = methodEvidence(candidate, candidate.stats, {
      ...options,
      threeSeqMaxOperations: exactOperations,
    }, Math.max(1, comparisons), nSites);
    if (threeSeqOperations <= exactOperations) exactThreeSeqBudget -= threeSeqOperations;
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
      confidenceStart: [Math.max(0, candidate.start - confidence), Math.min(nSites, candidate.start + confidence)],
      confidenceEnd: [Math.max(0, candidate.end - confidence), Math.min(nSites, candidate.end + confidence)],
      breakpointModel: candidate.breakpointModel,
      evidence,
      chiSquare: candidate.chiSquare,
      informativeSites: candidate.informative,
      decision: "unreviewed",
      warnings: addWarnings(candidate, sequences, nSites, options),
      note: [
        candidate.wraps ? "Origin-spanning circular tract." : "",
        candidate.alternatives.length
          ? `${candidate.alternatives.length} alternative minor-parent candidate${candidate.alternatives.length === 1 ? "" : "s"} grouped with this signal.`
          : "",
      ].filter(Boolean).join(" "),
      source: "wasm",
      supportedCount,
      diagnostics: candidate.diagnostics,
      groupId: null,
      history: [{
        id: `history-${jobId}-${index + 1}-1`,
        timestamp: new Date().toISOString(),
        action: "Detected by scan",
        summary: `${supportedCount} method families supported the initial hypothesis.`,
      }],
      evidenceStale: false,
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

  postMessage({
    type: "result",
    jobId,
    events,
    distance: Array.from(distance),
    comparisons,
    elapsedMs: performance.now() - started,
    timing: { distanceMs, scanMs, statisticsMs, diagnosticsMs },
    diagnostics: diagnosticProfile.summary,
    matrixCount,
    parentSamples: exactDistanceMatrix ? nSites : Math.min(parentSamples, nSites),
    matrixMode: exactDistanceMatrix ? "exact" : "24-sequence view + sampled/stratified parent search",
    engine: [
      "WebAssembly",
      exactDistanceMatrix ? "packed distance" : "sampled parent search",
      options.circular ? "dual-origin circular scan" : "linear scan",
      options.polishBreakpoints ? "two-state HMM polish" : "raw breakpoints",
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
  const profile = alignmentDiagnostics(encoded, nSeq, nSites, options.window);
  let working = encoded;
  let rawStart = event.start;
  let rawEnd = event.end;
  let rotation = 0;
  if (event.wraps && event.start > event.end) {
    const backgroundLength = event.start - event.end;
    rotation = (event.end + Math.floor(backgroundLength / 2)) % nSites;
    working = rotateSequences(encoded, nSeq, nSites, rotation);
    rawStart = (event.start - rotation + nSites) % nSites;
    rawEnd = (event.end - rotation + nSites) % nSites;
  }
  if (rawEnd <= rawStart) throw new Error("The edited event does not define a valid tract.");
  const seqPtr = 65536;
  const prefixAPtr = align(seqPtr + working.byteLength, 16);
  const prefixBPtr = prefixAPtr + (nSites + 1) * 4;
  const outPtr = align(prefixBPtr + (nSites + 1) * 4, 16);
  const statsPtr = align(outPtr + 64, 16);
  const requiredPages = Math.ceil((statsPtr + 96) / 65536);
  const currentPages = instance.exports.memory.buffer.byteLength / 65536;
  if (requiredPages > currentPages) instance.exports.memory.grow(requiredPages - currentPages);
  new Uint8Array(instance.exports.memory.buffer, seqPtr, working.byteLength).set(working);

  instance.exports.triplet_counts(
    seqPtr,
    nSites,
    event.recombinant,
    event.majorParent,
    event.minorParent,
    rawStart,
    rawEnd,
    outPtr,
  );
  const counts = new Int32Array(instance.exports.memory.buffer, outPtr, 6).slice();
  instance.exports.method_stats(
    seqPtr,
    nSites,
    event.recombinant,
    event.majorParent,
    event.minorParent,
    rawStart,
    rawEnd,
    Math.max(20, options.window),
    Math.max(1, options.step),
    Math.max(0, options.bootstrapReplicates ?? 100),
    (options.randomSeed ?? 0x5a17c0de) | 0,
    prefixAPtr,
    prefixBPtr,
    statsPtr,
  );
  const output = new Int32Array(instance.exports.memory.buffer, statsPtr, 23);
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
  };
  const candidate = {
    recombinant: event.recombinant,
    majorParent: event.majorParent,
    minorParent: event.minorParent,
    start: event.start,
    end: event.end,
    wraps: event.wraps,
    rawStart,
    rawEnd,
    rotation,
    sequencePtr: seqPtr,
    insideMinor: counts[1],
    insideMajor: counts[2],
    outsideMajor: counts[3],
    outsideMinor: counts[4],
    informative: counts[0],
    chiSquare: counts[5] / 1000,
    stats,
    diagnostics: null,
  };
  candidate.diagnostics = candidateDiagnostics(candidate, encoded, nSites, profile);
  const evidence = methodEvidence(candidate, stats, options, 1, nSites);
  postMessage({
    type: "recalculated",
    jobId,
    patch: {
      evidence,
      chiSquare: candidate.chiSquare,
      informativeSites: candidate.informative,
      warnings: addWarnings(candidate, sequences, nSites, options),
      diagnostics: candidate.diagnostics,
      evidenceStale: false,
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
