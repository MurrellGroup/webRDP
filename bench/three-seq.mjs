import fs from "node:fs";
import { performance } from "node:perf_hooks";

const nSeq = Math.max(3, Number(process.argv[2] ?? 32));
const nSites = Math.max(100, Number(process.argv[3] ?? 5_000));
const bytes = fs.readFileSync(new URL("../public/wasm/rdp.wasm", import.meta.url));
const instance = (await WebAssembly.instantiate(bytes)).instance;
const align = (value, multiple = 16) => Math.ceil(value / multiple) * multiple;
const seqPtr = 65_536;
const encoded = new Uint8Array(nSeq * nSites);
let state = 0x5a17c0de;
for (let index = 0; index < encoded.length; index += 1) {
  state = (Math.imul(state, 1_664_525) + 1_013_904_223) >>> 0;
  encoded[index] = state >>> 30;
}
const wordsPerSequence = Math.ceil(nSites / 16);
const packedPtr = align(seqPtr + encoded.byteLength);
const validityPtr = packedPtr + nSeq * wordsPerSequence * 4;
const workspacePtr = align(validityPtr + nSeq * wordsPerSequence * 4);
const workspaceBytes = instance.exports.source_three_seq_workspace_bytes(nSites);
const outPtr = align(workspacePtr + workspaceBytes);
const required = outPtr + 6 * 16 * 4;
const missingPages = Math.ceil((required - instance.exports.memory.buffer.byteLength) / 65_536);
if (missingPages > 0) instance.exports.memory.grow(missingPages);
new Uint8Array(instance.exports.memory.buffer, seqPtr, encoded.length).set(encoded);
const packed = new Uint32Array(instance.exports.memory.buffer, packedPtr, nSeq * wordsPerSequence);
const validity = new Uint32Array(instance.exports.memory.buffer, validityPtr, nSeq * wordsPerSequence);
for (let sequence = 0; sequence < nSeq; sequence += 1) {
  for (let site = 0; site < nSites; site += 1) {
    const word = sequence * wordsPerSequence + (site >>> 4);
    const shift = (site & 15) * 2;
    packed[word] |= encoded[sequence * nSites + site] << shift;
    validity[word] |= 1 << shift;
  }
}
const triplets = [];
for (let first = 0; first < nSeq - 2; first += 1) {
  for (let second = first + 1; second < nSeq - 1; second += 1) {
    for (let third = second + 1; third < nSeq; third += 1) triplets.push([first, second, third]);
  }
}

function run(packedMode) {
  let rows = 0;
  const started = performance.now();
  for (const [first, second, third] of triplets) {
    rows += packedMode
      ? instance.exports.scan_source_three_seq_triplet_packed_mode(
          packedPtr, validityPtr, wordsPerSequence, nSites, first, second, third, 0, workspacePtr, outPtr,
        )
      : instance.exports.scan_source_three_seq_triplet_mode(
          seqPtr, nSites, first, second, third, 0, workspacePtr, outPtr,
        );
  }
  return { milliseconds: performance.now() - started, rows };
}

run(true);
run(false);
const packedResult = run(true);
const byteResult = run(false);
const fastest = packedResult.milliseconds <= byteResult.milliseconds ? "packed" : "byte";
const bestMs = Math.min(packedResult.milliseconds, byteResult.milliseconds);
console.log(JSON.stringify({
  sequences: nSeq,
  sites: nSites,
  concreteTriplets: triplets.length,
  tripletSites: triplets.length * nSites,
  packedMs: Number(packedResult.milliseconds.toFixed(2)),
  byteMs: Number(byteResult.milliseconds.toFixed(2)),
  fastest,
  millionTripletSitesPerSecond: Number((triplets.length * nSites / bestMs / 1_000).toFixed(1)),
  outputRows: packedResult.rows,
  checkwrapWorkspaceBytes: workspaceBytes,
}, null, 2));
