import { performance } from "node:perf_hooks";
import fs from "node:fs";
import { buildSourceSiScanRandomization, runSourceSiScan } from "../public/rdp-siscan.js";
import { sourcePhiTest } from "../public/rdp-phi.js";

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
const rdpSignalCapacity = 32;
const rdpBestPtr = outPtr + rdpSignalCapacity * 72;
const chiSignalCapacity = 24;
const chiPeakCapacity = 8;
const chiPositionsPtr = align(rdpBestPtr + 72);
const chiScoresPtr = chiPositionsPtr + 4 * (nSites + 1) * 4;
const chiMissingPtr = align(chiScoresPtr + 4 * (nSites + 1), 4);
const chiProfilePtr = align(chiMissingPtr + (nSites + 1) * 4, 8);
const chiSmoothPtr = chiProfilePtr + (nSites + 1) * 8;
const chiPeakPtr = align(chiSmoothPtr + (nSites + 1) * 8, 4);
const chiOutPtr = chiPeakPtr + chiPeakCapacity * 6 * 4;
const geneconvSignalCapacity = 64;
const geneconvPositionsPtr = align(chiOutPtr + chiSignalCapacity * 16 * 4);
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
const roleCohortCount = Math.min(30, nSeq);
const roleCohortPtr = align(geneconvOutPtr + geneconvSignalCapacity * 16 * 4);
const tractMaskPtr = align(roleCohortPtr + roleCohortCount * 4);
const backgroundMaskPtr = tractMaskPtr + wordsPerSequence * 4;
const dmaxOutPtr = align(backgroundMaskPtr + wordsPerSequence * 4, 8);
const required = dmaxOutPtr + 40;
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
let retainedRdpSignals = 0;
let geneconvSignals = 0;
let geneconvMs = 0;
const concreteTriplets = new Set();
const scanStart = performance.now();
for (let recombinant = 0; recombinant < nSeq; recombinant += 1) {
  const parents = [];
  for (let candidate = 0; candidate < nSeq && parents.length < parentCount; candidate += 1) {
    if (candidate !== recombinant) parents.push(candidate);
  }
  for (let left = 0; left < parents.length; left += 1) {
    for (let right = left + 1; right < parents.length; right += 1) {
      const triplet = [recombinant, parents[left], parents[right]].sort((a, b) => a - b);
      const key = triplet.join(":");
      if (concreteTriplets.has(key)) continue;
      concreteTriplets.add(key);
      comparisons += 1;
      const rdpSignals = instance.exports.scan_rdp5_triplet_all_packed(
        packedPtr,
        validityPtr,
        wordsPerSequence,
        nSites,
        triplet[0],
        triplet[1],
        triplet[2],
        30,
        prefixBPtr,
        prefixAPtr,
        outPtr,
        rdpSignalCapacity,
        rdpBestPtr,
      );
      const retainedRdp = Math.min(rdpSignalCapacity, Math.max(0, rdpSignals));
      signals += rdpSignals;
      retainedRdpSignals += retainedRdp;
      const geneconvStarted = performance.now();
      const geneconv = instance.exports.scan_source_geneconv_all_packed(
        packedPtr,
        validityPtr,
        wordsPerSequence,
        nSites,
        triplet[0],
        triplet[1],
        triplet[2],
        1,
        0.05,
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
      geneconvMs += performance.now() - geneconvStarted;
      geneconvSignals += geneconv;
      signals += geneconv;
      signals += instance.exports.scan_source_chi_all_packed(
        packedPtr,
        validityPtr,
        wordsPerSequence,
        nSites,
        triplet[0],
        triplet[1],
        triplet[2],
        120,
        0,
        12,
        chiPositionsPtr,
        chiScoresPtr,
        chiMissingPtr,
        chiProfilePtr,
        chiSmoothPtr,
        chiPeakPtr,
        chiPeakCapacity,
        chiOutPtr,
        chiSignalCapacity,
      );
    }
  }
}
const scanMs = performance.now() - scanStart;
// VisRD is the most expensive term in the source recombinant-role consensus.
// Keep its production cohort and bit-packed path in the performance ledger so
// a superficially faster detector cannot hide a role-classification regression.
const roleCohort = new Int32Array(instance.exports.memory.buffer, roleCohortPtr, roleCohortCount);
for (let index = 0; index < roleCohortCount; index += 1) roleCohort[index] = index;
const tractMask = new Uint32Array(instance.exports.memory.buffer, tractMaskPtr, wordsPerSequence);
const backgroundMask = new Uint32Array(instance.exports.memory.buffer, backgroundMaskPtr, wordsPerSequence);
for (let site = 0; site < nSites; site += 1) {
  const word = site >>> 4;
  const bit = 1 << ((site & 15) * 2);
  if (site >= Math.floor(nSites / 3) && site < Math.floor(2 * nSites / 3)) tractMask[word] |= bit;
  else backgroundMask[word] |= bit;
}
const dmaxStart = performance.now();
if (roleCohortCount >= 4) instance.exports.dmax_visrd_packed(
  packedPtr,
  validityPtr,
  wordsPerSequence,
  roleCohortPtr,
  roleCohortCount,
  0,
  1,
  2,
  tractMaskPtr,
  backgroundMaskPtr,
  dmaxOutPtr,
);
const dmaxMs = performance.now() - dmaxStart;

// The source PHI graph is quadratic in retained informative sites. Exercise
// the worker's adaptive site ceiling against the full production taxon count
// so the browser safeguard cannot regress into an unnoticed stall.
const phiStart = performance.now();
const phiResult = sourcePhiTest(sequences, nSeq, nSites, {
  window: Math.min(100, nSites),
  maxInformativeSites: nSeq > 256 ? 160 : nSeq > 96 ? 256 : 384,
});
const sourcePhiMs = performance.now() - phiStart;

// Long-genome Sister-Scanning exercises the source-compatible 15-category
// vertical-permutation path and its bounded-memory random stream.  This is a
// separate worker-JavaScript hot path, not hidden inside the native total.
const sisterLength = 80_000;
const sisterEncoded = new Uint8Array(sisterLength * 4);
for (let site = 0; site < sisterLength; site += 1) {
  const inside = site >= 30_000 && site < 45_000;
  sisterEncoded[site] = inside ? 1 : 0;
  sisterEncoded[sisterLength + site] = 0;
  sisterEncoded[sisterLength * 2 + site] = 1;
  sisterEncoded[sisterLength * 3 + site] = 2;
}
const sisterRandomization = buildSourceSiScanRandomization(sisterLength, 1000, 1511506142);
const sisterStart = performance.now();
const sisterResult = runSourceSiScan(sisterEncoded, sisterLength, 4, [0, 1, 2], {
  window: 600,
  step: 100,
  scanPermutations: 100,
  pValuePermutations: 1000,
  candidatePool: [3],
  randomization: sisterRandomization,
  seed: 1511506142,
});
const sourceSiScanMs = performance.now() - sisterStart;
const choose = (count, size) => {
  if (count < size) return 0;
  let value = 1;
  for (let index = 1; index <= size; index += 1) value = value * (count - size + index) / index;
  return value;
};
const dmaxQuartets = choose(roleCohortCount, 4) - choose(Math.max(0, roleCohortCount - 3), 4);

const report = {
  dataset: `${nSeq} × ${nSites}`,
  nucleotides: nSeq * nSites,
  candidateParents: parentCount,
  comparisons,
  signals,
  geneconvSignals,
  geneconvMs: Number(geneconvMs.toFixed(2)),
  geneconvMillionTripletSitesPerSecond: Number(((comparisons * nSites) / (geneconvMs * 1000)).toFixed(1)),
  retainedSignals: retainedRdpSignals,
  scalarDistanceMs: Number(scalarDistanceMs.toFixed(2)),
  distanceMs: Number(distanceMs.toFixed(2)),
  packedDistanceSpeedup: Number((scalarDistanceMs / Math.max(distanceMs, 0.001)).toFixed(2)),
  scanMs: Number(scanMs.toFixed(2)),
  dmaxCohort: roleCohortCount,
  dmaxQuartets,
  dmaxMs: Number(dmaxMs.toFixed(2)),
  dmaxMillionSiteQuartetsPerSecond: Number(((dmaxQuartets * nSites) / (dmaxMs * 1000)).toFixed(1)),
  sourcePhiMs: Number(sourcePhiMs.toFixed(2)),
  sourcePhiInformativeSites: phiResult.informativeSites,
  sourcePhiTotalInformativeSites: phiResult.totalInformativeSites,
  sourceSiScanDataset: `4 × ${sisterLength}`,
  sourceSiScanMs: Number(sourceSiScanMs.toFixed(2)),
  sourceSiScanRecovered: sisterResult ? [sisterResult.start, sisterResult.end] : null,
  sourceSiScanMaterializedTable: Boolean(sisterRandomization.values),
  totalMs: Number((distanceMs + scanMs + dmaxMs + sourcePhiMs).toFixed(2)),
  millionSiteComparisonsPerSecond: Number(((comparisons * nSites) / (scanMs * 1000)).toFixed(1)),
};

console.log(JSON.stringify(report, null, 2));

if (process.env.RDP_PERFORMANCE_GATE === "1") {
  const failures = [];
  if (report.distanceMs > 150) failures.push(`packed distance ${report.distanceMs} ms > 150 ms`);
  if (report.dmaxMs > 500) failures.push(`packed VisRD dMax ${report.dmaxMs} ms > 500 ms`);
  if (report.sourcePhiMs > 500) failures.push(`bounded source PHI ${report.sourcePhiMs} ms > 500 ms`);
  if (report.sourceSiScanMs > 2_000) failures.push(`80 kb source SiScan ${report.sourceSiScanMs} ms > 2000 ms`);
  if (String(report.sourceSiScanRecovered) !== "30000,45000") failures.push(`source SiScan recovered ${report.sourceSiScanRecovered} instead of 30000,45000`);
  if (report.geneconvMs > 500) failures.push(`six-track source GENECONV ${report.geneconvMs} ms > 500 ms`);
  if (report.totalMs > 2_000) failures.push(`production total ${report.totalMs} ms > 2000 ms`);
  // The aggregate now includes three complete source detector families per
  // concrete triplet (RDP, six-track GENECONV, and MAXCHI/CHIMAERA).  Keep a
  // strict total-time gate and independently require the new kernel to exceed
  // 50 million full six-track triplet-sites/s.
  if (report.geneconvMillionTripletSitesPerSecond < 50) failures.push(`six-track source GENECONV throughput ${report.geneconvMillionTripletSitesPerSecond} M/s < 50 M/s`);
  if (failures.length) throw new Error(`Performance regression gate failed: ${failures.join("; ")}`);
}
