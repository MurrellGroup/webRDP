import { performance } from "node:perf_hooks";
import fs from "node:fs";

const nSeq = Number(process.argv[2] ?? 100);
const nSites = Number(process.argv[3] ?? 10_000);
const parentCount = Math.min(8, nSeq - 1);
const bytes = fs.readFileSync(new URL("../public/wasm/rdp.wasm", import.meta.url));
const { instance } = await WebAssembly.instantiate(bytes);
const align = (value, multiple = 16) => Math.ceil(value / multiple) * multiple;
const seqPtr = 65536;
const wordsPerSequence = Math.ceil(nSites / 16);
const packedPtr = align(seqPtr + nSeq * nSites);
const validityPtr = packedPtr + nSeq * wordsPerSequence * 4;
const distancePtr = align(validityPtr + nSeq * wordsPerSequence * 4);
const prefixAPtr = align(distancePtr + nSeq * nSeq * 4);
const prefixBPtr = prefixAPtr + (nSites + 1) * 4;
const outPtr = align(prefixBPtr + (nSites + 1) * 4);
const statsPtr = outPtr + 64;
const required = statsPtr + 96;
const missing = Math.ceil((required - instance.exports.memory.buffer.byteLength) / 65536);
if (missing > 0) instance.exports.memory.grow(missing);

const sequences = new Uint8Array(instance.exports.memory.buffer, seqPtr, nSeq * nSites);
let state = 0x9e3779b9;
for (let site = 0; site < nSites; site += 1) {
  state ^= state << 13; state ^= state >>> 17; state ^= state << 5;
  sequences[site] = state & 3;
}
for (let sequence = 1; sequence < nSeq; sequence += 1) {
  const source = sequences.subarray(0, nSites);
  sequences.set(source, sequence * nSites);
  for (let site = 0; site < nSites; site += 1) {
    state ^= state << 13; state ^= state >>> 17; state ^= state << 5;
    if ((state >>> 0) % 100 < 8) sequences[sequence * nSites + site] = (sequences[sequence * nSites + site] + 1 + (state & 1)) & 3;
  }
}

const packed = new Uint32Array(instance.exports.memory.buffer, packedPtr, nSeq * wordsPerSequence);
const validity = new Uint32Array(instance.exports.memory.buffer, validityPtr, nSeq * wordsPerSequence);
for (let sequence = 0; sequence < nSeq; sequence += 1) {
  for (let site = 0; site < nSites; site += 1) {
    const base = sequences[sequence * nSites + site];
    if (base >= 4) continue;
    const word = sequence * wordsPerSequence + (site >>> 4);
    const shift = (site & 15) * 2;
    packed[word] |= base << shift;
    validity[word] |= 1 << shift;
  }
}

const scalarDistanceStart = performance.now();
instance.exports.distance_matrix(seqPtr, nSeq, nSites, distancePtr);
const scalarDistanceMs = performance.now() - scalarDistanceStart;
const packedDistanceStart = performance.now();
instance.exports.distance_matrix_packed(
  packedPtr,
  validityPtr,
  nSeq,
  wordsPerSequence,
  distancePtr,
);
const distanceMs = performance.now() - packedDistanceStart;
let comparisons = 0;
let signals = 0;
const candidates = [];
const scanStart = performance.now();
for (let recombinant = 0; recombinant < nSeq; recombinant += 1) {
  const parents = [];
  for (let candidate = 0; candidate < nSeq && parents.length < parentCount; candidate += 1) {
    if (candidate !== recombinant) parents.push(candidate);
  }
  for (let left = 0; left < parents.length; left += 1) {
    for (let right = left + 1; right < parents.length; right += 1) {
      comparisons += 1;
      const found = instance.exports.scan_pair(
        seqPtr,
        nSites,
        recombinant,
        parents[left],
        parents[right],
        60,
        prefixAPtr,
        prefixBPtr,
        outPtr,
      );
      signals += found;
      if (found) {
        const result = new Int32Array(instance.exports.memory.buffer, outPtr, 12);
        candidates.push({
          recombinant,
          start: result[0],
          end: result[1],
          majorParent: result[2],
          minorParent: result[3],
          chiSquare: result[4],
        });
      }
    }
  }
}
const scanMs = performance.now() - scanStart;
const retained = candidates.sort((left, right) => right.chiSquare - left.chiSquare).slice(0, 500);
const statisticsStart = performance.now();
for (const candidate of retained) {
  instance.exports.method_stats(
    seqPtr,
    nSites,
    candidate.recombinant,
    candidate.majorParent,
    candidate.minorParent,
    candidate.start,
    candidate.end,
    120,
    5,
    100,
    1511506142,
    127,
    prefixAPtr,
    prefixBPtr,
    statsPtr,
  );
}
const statisticsMs = performance.now() - statisticsStart;
let hmmPolished = 0;
const hmmStart = performance.now();
for (const candidate of retained) {
  hmmPolished += instance.exports.hmm_polish(
    seqPtr,
    nSites,
    candidate.recombinant,
    candidate.majorParent,
    candidate.minorParent,
    candidate.start,
    candidate.end,
    prefixAPtr,
    prefixBPtr,
    outPtr,
  );
}
const hmmMs = performance.now() - hmmStart;

const report = {
  dataset: `${nSeq} × ${nSites}`,
  nucleotides: nSeq * nSites,
  candidateParents: parentCount,
  comparisons,
  signals,
  retainedSignals: retained.length,
  hmmPolished,
  scalarDistanceMs: Number(scalarDistanceMs.toFixed(2)),
  distanceMs: Number(distanceMs.toFixed(2)),
  packedDistanceSpeedup: Number((scalarDistanceMs / Math.max(distanceMs, 0.001)).toFixed(2)),
  scanMs: Number(scanMs.toFixed(2)),
  statisticsMs: Number(statisticsMs.toFixed(2)),
  hmmMs: Number(hmmMs.toFixed(2)),
  totalMs: Number((distanceMs + scanMs + statisticsMs + hmmMs).toFixed(2)),
  millionSiteComparisonsPerSecond: Number(((comparisons * nSites) / (scanMs * 1000)).toFixed(1)),
};

console.log(JSON.stringify(report, null, 2));

if (process.env.RDP_PERFORMANCE_GATE === "1") {
  const failures = [];
  if (report.distanceMs > 150) failures.push(`packed distance ${report.distanceMs} ms > 150 ms`);
  if (report.totalMs > 1_500) failures.push(`production total ${report.totalMs} ms > 1500 ms`);
  if (report.millionSiteComparisonsPerSecond < 30) failures.push(`scan throughput ${report.millionSiteComparisonsPerSecond} M/s < 30 M/s`);
  if (failures.length) throw new Error(`Performance regression gate failed: ${failures.join("; ")}`);
}
