import assert from "node:assert/strict";
import fs from "node:fs";
import test from "node:test";

const ROW_INTS = 16;

async function engine() {
  const bytes = fs.readFileSync(new URL("../public/wasm/rdp.wasm", import.meta.url));
  return (await WebAssembly.instantiate(bytes)).instance;
}

function pack(sequences) {
  const nSites = sequences[0].length;
  const wordsPerSequence = Math.ceil(nSites / 16);
  const packed = new Uint32Array(sequences.length * wordsPerSequence);
  const validity = new Uint32Array(sequences.length * wordsPerSequence);
  for (let sequence = 0; sequence < sequences.length; sequence += 1) {
    for (let site = 0; site < nSites; site += 1) {
      const base = sequences[sequence][site];
      if (base >= 4) continue;
      const word = sequence * wordsPerSequence + (site >>> 4);
      const shift = (site & 15) * 2;
      packed[word] |= base << shift;
      validity[word] |= 1 << shift;
    }
  }
  return { packed, validity, wordsPerSequence };
}

function scan(instance, sequences, triplets, { packedScan = true, seed = 1511506142 } = {}) {
  const align = (value, multiple = 16) => Math.ceil(value / multiple) * multiple;
  const nSeq = sequences.length;
  const nSites = sequences[0].length;
  const window = 80;
  const step = 10;
  const replicates = 50;
  const pairCount = nSeq * (nSeq - 1) / 2;
  const packedData = pack(sequences);
  const seqPtr = 65_536;
  const packedPtr = align(seqPtr + nSeq * nSites);
  const validityPtr = packedPtr + packedData.packed.byteLength;
  const tripletPtr = align(validityPtr + packedData.validity.byteLength);
  const pairMapPtr = tripletPtr + triplets.length * 12;
  const pairListPtr = pairMapPtr + pairCount * 4;
  const weightPtr = align(pairListPtr + pairCount * 8, 2);
  const pairDistancePtr = weightPtr + window * replicates * 2;
  const globalPairPtr = pairDistancePtr + pairCount * replicates * 2;
  const statePtr = align(globalPairPtr + pairCount * 8, 4);
  const differencePtr = statePtr + triplets.length * 24;
  const validPtr = differencePtr + replicates * 4;
  const lookupPtr = align(validPtr + replicates * 4, 2);
  const outputCapacity = 128;
  const lookupEntries = (window + 1) * (window + 2) / 2;
  const outPtr = align(lookupPtr + lookupEntries * 2, 4);
  const required = outPtr + outputCapacity * ROW_INTS * 4;
  const missing = Math.ceil((required - instance.exports.memory.buffer.byteLength) / 65_536);
  if (missing > 0) instance.exports.memory.grow(missing);

  const encoded = new Uint8Array(nSeq * nSites);
  sequences.forEach((sequence, index) => encoded.set(sequence, index * nSites));
  new Uint8Array(instance.exports.memory.buffer, seqPtr, encoded.length).set(encoded);
  new Uint32Array(instance.exports.memory.buffer, packedPtr, packedData.packed.length).set(packedData.packed);
  new Uint32Array(instance.exports.memory.buffer, validityPtr, packedData.validity.length).set(packedData.validity);
  new Int32Array(instance.exports.memory.buffer, tripletPtr, triplets.length * 3).set(triplets.flat());
  const shared = [
    nSeq, nSites, tripletPtr, triplets.length, window, step, replicates, 700,
    seed, pairMapPtr, pairListPtr, weightPtr, pairDistancePtr, globalPairPtr, statePtr,
    differencePtr, validPtr, lookupPtr, outPtr, outputCapacity,
  ];
  const total = packedScan
    ? instance.exports.scan_source_bootscan_batch_packed(
        packedPtr, validityPtr, packedData.wordsPerSequence, ...shared,
      )
    : instance.exports.scan_source_bootscan_batch(seqPtr, ...shared);
  const rows = Array.from({ length: Math.min(total, outputCapacity) }, (_, index) => (
    Array.from(new Int32Array(instance.exports.memory.buffer, outPtr + index * ROW_INTS * 4, ROW_INTS))
  ));
  return { rows, total };
}

function mosaic(length = 480, decoyBase = 2) {
  const recombinant = new Uint8Array(length);
  const major = new Uint8Array(length);
  const minor = new Uint8Array(length);
  const decoy = new Uint8Array(length);
  for (let site = 0; site < length; site += 1) {
    major[site] = site & 3;
    minor[site] = (site + 1) & 3;
    recombinant[site] = site >= 150 && site < 300 ? minor[site] : major[site];
    decoy[site] = (site + decoyBase) & 3;
  }
  return [recombinant, major, minor, decoy];
}

function divergentMosaic(length = 480) {
  const first = new Uint8Array(length);
  const second = new Uint8Array(length);
  const third = new Uint8Array(length);
  const decoy = new Uint8Array(length);
  for (let site = 0; site < length; site += 1) {
    const inside = site >= 160 && site < 280;
    if (inside && site % 4 !== 0) {
      // Topology 2 (second/third) dominates the implanted tract.
      first[site] = 0;
      second[site] = 1;
      third[site] = 1;
    } else if (!inside && site % 10 < 3) {
      // Topology 1 (first/third) is the strongest whole-alignment identity.
      first[site] = 0;
      second[site] = 1;
      third[site] = 0;
    } else {
      first[site] = 0;
      second[site] = 1;
      third[site] = 2;
    }
    decoy[site] = 3;
  }
  return [first, second, third, decoy];
}

test("source RecScan batch recovers the alternate topology from one concrete triplet", async () => {
  const instance = await engine();
  const result = scan(instance, mosaic(), [[0, 1, 2]]);
  const signal = result.rows.find((row) => row[3] === 1);
  assert.ok(signal, "recombinant/minor topology must replace the whole-alignment recombinant/major topology");
  assert.equal(signal[13], 0, "the whole-alignment closest pair is the baseline topology");
  assert.ok(signal[4] <= 150 && signal[5] >= 300, "sliding-window bounds must contain the implanted tract");
  assert.ok(signal[6] >= 35, "the alternate topology must exceed the 70% source cutoff");
  assert.equal(signal[7], 50);
  assert.ok(signal[9] > 0 && signal[10] === 0, "BSSubSeq/MakeScoresBS pair matches are triplet-local");
});

test("packed production RecScan equals the byte oracle, is seeded, and ignores unrelated decoys", async () => {
  const instance = await engine();
  const byte = scan(instance, mosaic(480, 2), [[0, 1, 2]], { packedScan: false });
  const packed = scan(instance, mosaic(480, 2), [[0, 1, 2]], { packedScan: true });
  const changedDecoy = scan(instance, mosaic(480, 3), [[0, 1, 2]], { packedScan: true });
  const repeated = scan(instance, mosaic(480, 2), [[0, 1, 2]], { packedScan: true });
  assert.deepEqual(packed, byte);
  assert.deepEqual(changedDecoy, packed, "no rest-of-alignment proxy may enter the concrete triplet");
  assert.deepEqual(repeated, packed, "the RDP seed must reproduce the same bootstrap table");
});

test("the batch accepts several triplets while sharing one pair-distance pass", async () => {
  const instance = await engine();
  const result = scan(instance, mosaic(), [[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]]);
  assert.ok(result.rows.some((row) => row[0] === 0 && row[1] === 1 && row[2] === 2));
  assert.ok(result.rows.every((row) => row[0] < row[1] && row[1] < row[2]));
});

test("whole-alignment identity retains the RDP5 baseline beyond JC saturation", async () => {
  const instance = await engine();
  const result = scan(instance, divergentMosaic(), [[0, 1, 2]]);
  const signal = result.rows.find((row) => row[3] === 2);
  assert.ok(signal, "the second/third tract topology must be recovered in a divergent triplet");
  assert.equal(signal[13], 1, "global first/third identity, not a saturated JC tie, defines the baseline");
});
