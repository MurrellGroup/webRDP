import assert from "node:assert/strict";
import fs from "node:fs";
import test from "node:test";

async function engine() {
  const bytes = fs.readFileSync(new URL("../public/wasm/rdp.wasm", import.meta.url));
  return (await WebAssembly.instantiate(bytes)).instance;
}

function reserve(instance, nSeq, nSites) {
  const align = (value, multiple = 16) => Math.ceil(value / multiple) * multiple;
  const seqPtr = 65536;
  const distancePtr = align(seqPtr + nSeq * nSites);
  const prefixAPtr = align(distancePtr + nSeq * nSeq * 4);
  const prefixBPtr = prefixAPtr + (nSites + 1) * 4;
  const outPtr = align(prefixBPtr + (nSites + 1) * 4);
  const statsPtr = outPtr + 64;
  const wordsPerSequence = Math.ceil(nSites / 16);
  const packedPtr = align(statsPtr + 96);
  const validityPtr = packedPtr + nSeq * wordsPerSequence * 4;
  const poolPtr = align(validityPtr + nSeq * wordsPerSequence * 4);
  const nearestIndexesPtr = poolPtr + nSeq * 4;
  const nearestDistancesPtr = nearestIndexesPtr + nSeq * 4;
  const required = nearestDistancesPtr + nSeq * 4;
  const missing = Math.ceil((required - instance.exports.memory.buffer.byteLength) / 65536);
  if (missing > 0) instance.exports.memory.grow(missing);
  return { seqPtr, distancePtr, prefixAPtr, prefixBPtr, outPtr, statsPtr, packedPtr, validityPtr, poolPtr, nearestIndexesPtr, nearestDistancesPtr, wordsPerSequence };
}

function pack(bytes, nSeq, nSites) {
  const wordsPerSequence = Math.ceil(nSites / 16);
  const packed = new Uint32Array(nSeq * wordsPerSequence);
  const validity = new Uint32Array(nSeq * wordsPerSequence);
  for (let sequence = 0; sequence < nSeq; sequence += 1) {
    for (let site = 0; site < nSites; site += 1) {
      const base = bytes[sequence * nSites + site];
      if (base >= 4) continue;
      const word = sequence * wordsPerSequence + (site >>> 4);
      const shift = (site & 15) * 2;
      packed[word] |= base << shift;
      validity[word] |= 1 << shift;
    }
  }
  return { packed, validity };
}

test("distance kernel returns canonical p-distances", async () => {
  const instance = await engine();
  const nSeq = 3;
  const nSites = 100;
  const pointers = reserve(instance, nSeq, nSites);
  const sequences = new Uint8Array(instance.exports.memory.buffer, pointers.seqPtr, nSeq * nSites);
  sequences.fill(0, 0, nSites);
  sequences.fill(1, nSites, nSites * 2);
  sequences.fill(0, nSites * 2, nSites * 3);
  sequences.fill(1, nSites * 2 + 40, nSites * 2 + 60);
  instance.exports.distance_matrix(pointers.seqPtr, nSeq, nSites, pointers.distancePtr);
  const matrix = new Float32Array(instance.exports.memory.buffer, pointers.distancePtr, 9);
  assert.equal(matrix[0], 0);
  assert.equal(matrix[1], 1);
  assert.ok(Math.abs(matrix[2] - 0.2) < 1e-6);
  assert.ok(Math.abs(matrix[5] - 0.8) < 1e-6);
});

test("packed p-distance kernel is bit-for-bit equivalent to the scalar kernel", async () => {
  const instance = await engine();
  const nSeq = 3;
  const nSites = 101;
  const pointers = reserve(instance, nSeq, nSites);
  const sequences = new Uint8Array(instance.exports.memory.buffer, pointers.seqPtr, nSeq * nSites);
  for (let index = 0; index < sequences.length; index += 1) sequences[index] = (index * 17 + (index >>> 2)) % 6;
  instance.exports.distance_matrix(pointers.seqPtr, nSeq, nSites, pointers.distancePtr);
  const scalar = new Float32Array(instance.exports.memory.buffer, pointers.distancePtr, nSeq * nSeq).slice();
  const { packed, validity } = pack(sequences, nSeq, nSites);
  new Uint32Array(instance.exports.memory.buffer, pointers.packedPtr, packed.length).set(packed);
  new Uint32Array(instance.exports.memory.buffer, pointers.validityPtr, validity.length).set(validity);
  instance.exports.distance_matrix_packed(
    pointers.packedPtr,
    pointers.validityPtr,
    nSeq,
    pointers.wordsPerSequence,
    pointers.distancePtr,
  );
  const optimized = new Float32Array(instance.exports.memory.buffer, pointers.distancePtr, nSeq * nSeq);
  assert.deepEqual([...optimized], [...scalar]);
});

test("sampled nearest-candidate kernel retains the closest references", async () => {
  const instance = await engine();
  const nSeq = 5;
  const nSites = 160;
  const pointers = reserve(instance, nSeq, nSites);
  const sequences = new Uint8Array(instance.exports.memory.buffer, pointers.seqPtr, nSeq * nSites);
  sequences.fill(0);
  sequences.fill(1, nSites * 3, nSites * 4);
  for (let site = 0; site < 5; site += 1) sequences[nSites + site] = 1;
  for (let site = 0; site < 10; site += 1) sequences[nSites * 2 + site] = 1;
  for (let site = 0; site < 2; site += 1) sequences[nSites * 4 + site] = 1;
  new Int32Array(instance.exports.memory.buffer, pointers.poolPtr, 4).set([1, 2, 3, 4]);
  const found = instance.exports.nearest_candidates_sampled(
    pointers.seqPtr,
    nSites,
    0,
    pointers.poolPtr,
    4,
    nSites,
    2,
    pointers.nearestIndexesPtr,
    pointers.nearestDistancesPtr,
  );
  assert.equal(found, 2);
  assert.deepEqual(
    [...new Int32Array(instance.exports.memory.buffer, pointers.nearestIndexesPtr, 2)],
    [4, 1],
  );
});

test("triplet kernel localizes a known internal mosaic tract", async () => {
  const instance = await engine();
  const nSeq = 3;
  const nSites = 300;
  const pointers = reserve(instance, nSeq, nSites);
  const sequences = new Uint8Array(instance.exports.memory.buffer, pointers.seqPtr, nSeq * nSites);
  sequences.fill(0, 0, nSites);
  sequences.fill(1, nSites, nSites * 2);
  sequences.fill(0, nSites * 2, nSites * 3);
  sequences.fill(1, nSites * 2 + 90, nSites * 2 + 210);
  const found = instance.exports.scan_pair(
    pointers.seqPtr,
    nSites,
    2,
    0,
    1,
    30,
    pointers.prefixAPtr,
    pointers.prefixBPtr,
    pointers.outPtr,
  );
  assert.equal(found, 1);
  const result = new Int32Array(instance.exports.memory.buffer, pointers.outPtr, 12);
  assert.equal(result[0], 90);
  assert.equal(result[1], 210);
  assert.equal(result[2], 0);
  assert.equal(result[3], 1);
  assert.ok(result[4] / 1000 > 250);

  instance.exports.method_stats(
    pointers.seqPtr,
    nSites,
    2,
    0,
    1,
    result[0],
    result[1],
    60,
    5,
    100,
    12345,
    pointers.prefixAPtr,
    pointers.prefixBPtr,
    pointers.statsPtr,
  );
  const stats = new Int32Array(instance.exports.memory.buffer, pointers.statsPtr, 23);
  assert.equal(stats[0], 120, "GENECONV G=0 run should span the inserted tract");
  assert.ok(stats[5] / stats[6] > 0.9, "window topology should switch with the known tract");
  assert.ok(Math.abs(stats[17] - 90) <= 2, "left breakpoint polishing should remain at the truth");
  assert.ok(Math.abs(stats[18] - 210) <= 2, "right breakpoint polishing should remain at the truth");
  assert.equal(stats[11], 120, "3Seq maximum HGRW descent should span the inserted tract");
  assert.equal(stats[19], 180, "3Seq should count major-parent matches");
  assert.equal(stats[20], 120, "3Seq should count minor-parent matches");
  assert.equal(stats[22], 300, "three representative windows should each contribute 100 decisive bootstraps");
  assert.equal(stats[21], 300, "all bootstrap topologies should match the known mosaic");

  const hmmFound = instance.exports.hmm_polish(
    pointers.seqPtr,
    nSites,
    2,
    0,
    1,
    82,
    218,
    pointers.prefixAPtr,
    pointers.prefixBPtr,
    pointers.outPtr,
  );
  assert.equal(hmmFound, 1);
  const hmm = new Int32Array(instance.exports.memory.buffer, pointers.outPtr, 11);
  assert.ok(Math.abs(hmm[0] - 90) <= 1, `HMM left breakpoint ${hmm[0]} should recover 90`);
  assert.ok(Math.abs(hmm[1] - 210) <= 1, `HMM right breakpoint ${hmm[1]} should recover 210`);
  assert.equal(hmm[2], 300, "all sites are informative in the synthetic triplet");
  assert.equal(hmm[3], 2, "the two-state path should contain the two known ancestry switches");
  assert.ok(hmm[4] >= 950 && hmm[5] >= 850, "candidate-seeded emissions should fit both ancestry states");
});
