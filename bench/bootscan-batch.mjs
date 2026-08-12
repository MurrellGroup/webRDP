import fs from "node:fs";
import { performance } from "node:perf_hooks";

const nSeq = Math.max(3, Number(process.argv[2] ?? 24));
const nSites = Math.max(400, Number(process.argv[3] ?? 2_000));
const window = Math.min(200, Math.floor(nSites / 2));
const step = 20;
const replicates = 100;
const cutoffPermille = 700;
const pairCount = nSeq * (nSeq - 1) / 2;
const tripletCount = nSeq * (nSeq - 1) * (nSeq - 2) / 6;
const windowCount = Math.floor(nSites / step) + 2;
const align = (value, multiple = 16) => Math.ceil(value / multiple) * multiple;

const wasm = fs.readFileSync(new URL("../public/wasm/rdp.wasm", import.meta.url));
const { instance } = await WebAssembly.instantiate(wasm);
const wordsPerSequence = Math.ceil(nSites / 16);
const packedPtr = 65_536;
const validityPtr = packedPtr + nSeq * wordsPerSequence * 4;
const tripletPtr = align(validityPtr + nSeq * wordsPerSequence * 4);
const pairMapPtr = tripletPtr + tripletCount * 12;
const pairListPtr = pairMapPtr + pairCount * 4;
const weightPtr = align(pairListPtr + pairCount * 8, 2);
const pairDistancePtr = weightPtr + window * replicates * 2;
const globalPairPtr = pairDistancePtr + pairCount * replicates * 2;
const statePtr = align(globalPairPtr + pairCount * 8, 4);
const differencePtr = statePtr + tripletCount * 24;
const validPtr = differencePtr + replicates * 4;
const lookupPtr = align(validPtr + replicates * 4, 2);
const lookupEntries = (window + 1) * (window + 2) / 2;
const outputCapacity = 8_192;
const outPtr = align(lookupPtr + lookupEntries * 2, 4);
const requiredBytes = outPtr + outputCapacity * 16 * 4;
const missingPages = Math.ceil((requiredBytes - instance.exports.memory.buffer.byteLength) / 65_536);
if (missingPages > 0) instance.exports.memory.grow(missingPages);

const packed = new Uint32Array(instance.exports.memory.buffer, packedPtr, nSeq * wordsPerSequence);
const validity = new Uint32Array(instance.exports.memory.buffer, validityPtr, nSeq * wordsPerSequence);
let randomState = 0x9e3779b9;
const nextRandom = () => {
  randomState ^= randomState << 13;
  randomState ^= randomState >>> 17;
  randomState ^= randomState << 5;
  return randomState >>> 0;
};
const ancestor = new Uint8Array(nSites);
for (let site = 0; site < nSites; site += 1) ancestor[site] = nextRandom() & 3;
for (let sequence = 0; sequence < nSeq; sequence += 1) {
  for (let site = 0; site < nSites; site += 1) {
    let base = ancestor[site];
    if (sequence > 0 && nextRandom() % 100 < 6) base = (base + 1 + (nextRandom() & 1)) & 3;
    const word = sequence * wordsPerSequence + (site >>> 4);
    const shift = (site & 15) * 2;
    packed[word] |= base << shift;
    validity[word] |= 1 << shift;
  }
}

const triplets = new Int32Array(instance.exports.memory.buffer, tripletPtr, tripletCount * 3);
let triplet = 0;
for (let first = 0; first < nSeq - 2; first += 1) {
  for (let second = first + 1; second < nSeq - 1; second += 1) {
    for (let third = second + 1; third < nSeq; third += 1) {
      triplets[triplet * 3] = first;
      triplets[triplet * 3 + 1] = second;
      triplets[triplet * 3 + 2] = third;
      triplet += 1;
    }
  }
}

const args = [
  packedPtr,
  validityPtr,
  wordsPerSequence,
  nSeq,
  nSites,
  tripletPtr,
  tripletCount,
  window,
  step,
  replicates,
  cutoffPermille,
  0x5a17c0de,
  pairMapPtr,
  pairListPtr,
  weightPtr,
  pairDistancePtr,
  globalPairPtr,
  statePtr,
  differencePtr,
  validPtr,
  lookupPtr,
  outPtr,
  outputCapacity,
];

instance.exports.scan_source_bootscan_batch_packed(...args);
const started = performance.now();
const signals = instance.exports.scan_source_bootscan_batch_packed(...args);
const elapsedMs = performance.now() - started;
const naivePairRowsPerWindow = tripletCount * 3;
const report = {
  dataset: `${nSeq} × ${nSites}`,
  triplets: tripletCount,
  usedPairs: pairCount,
  windows: windowCount,
  replicates,
  elapsedMs: Number(elapsedMs.toFixed(2)),
  signals,
  pairReuseFactor: Number((naivePairRowsPerWindow / pairCount).toFixed(1)),
  tripletReplicateEvaluationsPerSecond: Math.round(tripletCount * windowCount * replicates / (elapsedMs / 1000)),
  workspaceMiB: Number(((requiredBytes - tripletPtr) / (1024 * 1024)).toFixed(2)),
};
console.log(JSON.stringify(report, null, 2));

if (process.env.RDP_BOOTSCAN_GATE === "1" && elapsedMs > 2_000) {
  throw new Error(`BootScan batch exceeded 2 s gate: ${elapsedMs.toFixed(1)} ms`);
}
