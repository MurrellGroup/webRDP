import assert from "node:assert/strict";
import fs from "node:fs";
import test from "node:test";
import { buildNeighborJoiningPathMatrix } from "../public/rdp-bootstrap-tree.js";

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

function scan(instance, sequences, triplets, {
  packedScan = true,
  seed = 1511506142,
  mode = "distance",
  legacy = false,
} = {}) {
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
  const modeCode = mode === "upgma" ? 1 : mode === "neighbor-joining" ? 2 : 0;
  const treeWorkPtr = align(lookupPtr + lookupEntries * 2, 16);
  const treeWorkspaceBytes = modeCode === 0
    ? 0
    : instance.exports.source_bootscan_tree_workspace_bytes(nSeq);
  const outPtr = align(treeWorkPtr + treeWorkspaceBytes, 4);
  const required = outPtr + outputCapacity * ROW_INTS * 4;
  const missing = Math.ceil((required - instance.exports.memory.buffer.byteLength) / 65_536);
  if (missing > 0) instance.exports.memory.grow(missing);

  const encoded = new Uint8Array(nSeq * nSites);
  sequences.forEach((sequence, index) => encoded.set(sequence, index * nSites));
  new Uint8Array(instance.exports.memory.buffer, seqPtr, encoded.length).set(encoded);
  new Uint32Array(instance.exports.memory.buffer, packedPtr, packedData.packed.length).set(packedData.packed);
  new Uint32Array(instance.exports.memory.buffer, validityPtr, packedData.validity.length).set(packedData.validity);
  new Int32Array(instance.exports.memory.buffer, tripletPtr, triplets.length * 3).set(triplets.flat());
  const legacyShared = [
    nSeq, nSites, tripletPtr, triplets.length, window, step, replicates, 700,
    seed, pairMapPtr, pairListPtr, weightPtr, pairDistancePtr, globalPairPtr, statePtr,
    differencePtr, validPtr, lookupPtr, outPtr, outputCapacity,
  ];
  const modeShared = [
    nSeq, nSites, tripletPtr, triplets.length, window, step, replicates, 700,
    modeCode, seed, pairMapPtr, pairListPtr, weightPtr, pairDistancePtr, globalPairPtr,
    statePtr, differencePtr, validPtr, lookupPtr, treeWorkPtr, outPtr, outputCapacity,
  ];
  const total = legacy
    ? packedScan
      ? instance.exports.scan_source_bootscan_batch_packed(
          packedPtr, validityPtr, packedData.wordsPerSequence, ...legacyShared,
        )
      : instance.exports.scan_source_bootscan_batch(seqPtr, ...legacyShared)
    : packedScan
      ? instance.exports.scan_source_bootscan_batch_mode_packed(
          packedPtr, validityPtr, packedData.wordsPerSequence, ...modeShared,
        )
      : instance.exports.scan_source_bootscan_batch_mode(seqPtr, ...modeShared);
  const rows = Array.from({ length: Math.min(total, outputCapacity) }, (_, index) => (
    Array.from(new Int32Array(instance.exports.memory.buffer, outPtr + index * ROW_INTS * 4, ROW_INTS))
  ));
  const pairMap = new Int32Array(instance.exports.memory.buffer, pairMapPtr, pairCount);
  const usedPairs = Math.max(-1, ...pairMap) + 1;
  return { rows, total, usedPairs };
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

test("the mode-aware distance ABI is exactly backward compatible", async () => {
  const instance = await engine();
  const current = scan(instance, mosaic(), [[0, 1, 2]], { mode: "distance" });
  const legacy = scan(instance, mosaic(), [[0, 1, 2]], { legacy: true });
  assert.deepEqual(current, legacy);
});

test("UPGMA and neighbor-joining transform one full-cohort matrix into tree paths", async () => {
  const instance = await engine();
  const align = (value, multiple = 16) => Math.ceil(value / multiple) * multiple;
  const nSeq = 4;
  const replicates = 1;
  const pairCount = 6;
  const pairMapPtr = 65_536;
  const pairDistancePtr = align(pairMapPtr + pairCount * 4, 2);
  const treeWorkPtr = align(pairDistancePtr + pairCount * 2, 16);
  const required = treeWorkPtr + instance.exports.source_bootscan_tree_workspace_bytes(nSeq);
  const missing = Math.ceil((required - instance.exports.memory.buffer.byteLength) / 65_536);
  if (missing > 0) instance.exports.memory.grow(missing);
  new Int32Array(instance.exports.memory.buffer, pairMapPtr, pairCount).set([0, 1, 2, 3, 4, 5]);
  const sourceDistances = [1, 10, 10, 10, 10, 2]; // 01 and 23 are the two cherries
  const transform = (modeCode) => {
    new Uint16Array(instance.exports.memory.buffer, pairDistancePtr, pairCount).set(sourceDistances);
    instance.exports.source_bootscan_transform_tree_relationships(
      modeCode, nSeq, replicates, 0, pairMapPtr, pairDistancePtr, treeWorkPtr,
    );
    return [...new Uint16Array(instance.exports.memory.buffer, pairDistancePtr, pairCount)];
  };
  assert.deepEqual(transform(1), [2, 4, 4, 4, 4, 2], "UPGMA paths must include the midpoint root");
  assert.deepEqual(transform(2), [2, 3, 3, 3, 3, 2], "NJ paths must retain the cohort split");
});

test("the WASM NJ transform matches the independent split-path implementation", async () => {
  const instance = await engine();
  const align = (value, multiple = 16) => Math.ceil(value / multiple) * multiple;
  let randomState = 0x4d595df4;
  const random = () => {
    randomState ^= randomState << 13;
    randomState ^= randomState >>> 17;
    randomState ^= randomState << 5;
    return randomState >>> 0;
  };
  for (let nSeq = 4; nSeq <= 10; nSeq += 1) {
    const pairCount = nSeq * (nSeq - 1) / 2;
    const pairMapPtr = 65_536;
    const pairDistancePtr = align(pairMapPtr + pairCount * 4, 2);
    const treeWorkPtr = align(pairDistancePtr + pairCount * 2, 16);
    const required = treeWorkPtr + instance.exports.source_bootscan_tree_workspace_bytes(nSeq);
    const missing = Math.ceil((required - instance.exports.memory.buffer.byteLength) / 65_536);
    if (missing > 0) instance.exports.memory.grow(missing);
    new Int32Array(instance.exports.memory.buffer, pairMapPtr, pairCount)
      .set(Array.from({ length: pairCount }, (_, index) => index));
    const triangular = new Uint16Array(pairCount);
    const square = new Float64Array(nSeq * nSeq);
    let pair = 0;
    for (let left = 0; left < nSeq - 1; left += 1) {
      for (let right = left + 1; right < nSeq; right += 1) {
        const distance = 1 + (random() % 30_000);
        triangular[pair++] = distance;
        square[left * nSeq + right] = distance;
        square[right * nSeq + left] = distance;
      }
    }
    new Uint16Array(instance.exports.memory.buffer, pairDistancePtr, pairCount).set(triangular);
    instance.exports.source_bootscan_transform_tree_relationships(
      2, nSeq, 1, 0, pairMapPtr, pairDistancePtr, treeWorkPtr,
    );
    const actual = [...new Uint16Array(instance.exports.memory.buffer, pairDistancePtr, pairCount)];
    const expectedSquare = buildNeighborJoiningPathMatrix(square, nSeq);
    const expected = [];
    for (let left = 0; left < nSeq - 1; left += 1) {
      for (let right = left + 1; right < nSeq; right += 1) {
        expected.push(expectedSquare[left * nSeq + right]);
      }
    }
    assert.deepEqual(actual, expected, `NJ split paths must agree for ${nSeq} taxa`);
  }
});

test("NJ relationships are inferred in full-cohort context, not from three isolated pairs", async () => {
  const instance = await engine();
  const align = (value, multiple = 16) => Math.ceil(value / multiple) * multiple;
  const nSeq = 4;
  const pairCount = 6;
  const pairMapPtr = 65_536;
  const pairDistancePtr = align(pairMapPtr + pairCount * 4, 2);
  const treeWorkPtr = align(pairDistancePtr + pairCount * 2, 16);
  const required = treeWorkPtr + instance.exports.source_bootscan_tree_workspace_bytes(nSeq);
  const missing = Math.ceil((required - instance.exports.memory.buffer.byteLength) / 65_536);
  if (missing > 0) instance.exports.memory.grow(missing);
  new Int32Array(instance.exports.memory.buffer, pairMapPtr, pairCount).set([0, 1, 2, 3, 4, 5]);
  // Within triplet 0/1/2 the raw distances prefer 0/1 (2,4,4). Taxon 3's
  // distances alter the four-taxon NJ split, making 0/2 the closest stored
  // tree relationship. A per-triplet tree or three independent pair screens
  // cannot produce this result.
  new Uint16Array(instance.exports.memory.buffer, pairDistancePtr, pairCount)
    .set([2, 4, 1, 4, 1, 4]);
  instance.exports.source_bootscan_transform_tree_relationships(
    2, nSeq, 1, 0, pairMapPtr, pairDistancePtr, treeWorkPtr,
  );
  const paths = [...new Uint16Array(instance.exports.memory.buffer, pairDistancePtr, pairCount)];
  assert.deepEqual(paths, [3, 2, 3, 3, 2, 3]);
  assert.ok(paths[1] < paths[0] && paths[1] < paths[3]);
});

test("tree-mode RecScan calculates every cohort pair even for one shortlisted triplet", async () => {
  const instance = await engine();
  for (const mode of ["upgma", "neighbor-joining"]) {
    const packed = scan(instance, mosaic(), [[0, 1, 2]], { mode });
    const byte = scan(instance, mosaic(), [[0, 1, 2]], { mode, packedScan: false });
    assert.equal(packed.usedPairs, 6, `${mode} must include all pairs of the four-sequence cohort`);
    assert.deepEqual(packed, byte, `${mode} packed and byte paths must agree`);
    assert.ok(packed.rows.some((row) => row[3] === 1), `${mode} must recover the implanted topology`);
  }
});

test("16-site packed blocks retain byte-exact missing and circular-boundary behavior", async () => {
  const instance = await engine();
  const sequences = mosaic(487).map((sequence) => sequence.slice());
  const missingSites = [0, 7, 15, 16, 79, 80, 159, 320, 479, 486];
  for (let sequence = 0; sequence < sequences.length; sequence += 1) {
    for (let index = sequence; index < missingSites.length; index += 3) {
      sequences[sequence][missingSites[index]] = 4;
    }
  }
  const triplets = [[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]];
  for (const mode of ["distance", "upgma", "neighbor-joining"]) {
    assert.deepEqual(
      scan(instance, sequences, triplets, { mode }),
      scan(instance, sequences, triplets, { mode, packedScan: false }),
      `${mode} packed blocks must match scalar extraction across word/origin boundaries`,
    );
  }
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
