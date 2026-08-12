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

function scan(instance, sequences, { gScale = 1, packedScan = false, pCutoff = 1 } = {}) {
  const nSites = sequences[0].length;
  const align = (value, multiple = 16) => Math.ceil(value / multiple) * multiple;
  const seqPtr = 65_536;
  const encodedBytes = sequences.length * nSites;
  const packedData = pack(sequences);
  const packedPtr = align(seqPtr + encodedBytes);
  const validityPtr = packedPtr + packedData.packed.byteLength;
  const positionsPtr = align(validityPtr + packedData.validity.byteLength);
  const categoriesPtr = positionsPtr + nSites * 4;
  const runStartPtr = align(categoriesPtr + nSites, 4);
  const runEndPtr = runStartPtr + nSites * 4;
  const runScorePtr = runEndPtr + nSites * 4;
  const prefixPtr = align(runScorePtr + nSites * 4, 8);
  const treePtr = prefixPtr + (nSites + 1) * 8;
  const workspaceBytes = (nSites + 1) * 16 + 8;
  const calibrationPtr = align(treePtr + workspaceBytes, 8);
  const candidatePtr = align(calibrationPtr + 6 * 40, 8);
  const candidateCapacity = 3 * (nSites + 1);
  const deletePtr = align(candidatePtr + candidateCapacity * 24, 4);
  const outCapacity = 128;
  const outPtr = align(deletePtr + nSites * 4, 8);
  const required = outPtr + outCapacity * ROW_INTS * 4;
  const missingPages = Math.ceil((required - instance.exports.memory.buffer.byteLength) / 65_536);
  if (missingPages > 0) instance.exports.memory.grow(missingPages);

  const encoded = new Uint8Array(sequences.length * nSites);
  sequences.forEach((sequence, index) => encoded.set(sequence, index * nSites));
  new Uint8Array(instance.exports.memory.buffer, seqPtr, encoded.length).set(encoded);
  new Uint32Array(instance.exports.memory.buffer, packedPtr, packedData.packed.length).set(packedData.packed);
  new Uint32Array(instance.exports.memory.buffer, validityPtr, packedData.validity.length).set(packedData.validity);
  const scratch = [
    positionsPtr,
    categoriesPtr,
    runStartPtr,
    runEndPtr,
    runScorePtr,
    prefixPtr,
    treePtr,
    calibrationPtr,
    candidatePtr,
    candidateCapacity,
    deletePtr,
    outPtr,
    outCapacity,
  ];
  const total = packedScan
    ? instance.exports.scan_source_geneconv_all_packed(
        packedPtr,
        validityPtr,
        packedData.wordsPerSequence,
        nSites,
        0,
        1,
        2,
        gScale,
        pCutoff,
        ...scratch,
      )
    : instance.exports.scan_source_geneconv_all(
        seqPtr,
        nSites,
        0,
        1,
        2,
        gScale,
        pCutoff,
        ...scratch,
      );
  const rows = Array.from({ length: Math.min(total, outCapacity) }, (_, index) => (
    Array.from(new Int32Array(instance.exports.memory.buffer, outPtr + index * ROW_INTS * 4, ROW_INTS))
  ));
  return { rows, total };
}

function sourceTriplet(length = 600, holes = []) {
  const sequences = [new Uint8Array(length), new Uint8Array(length), new Uint8Array(length)];
  const holeSet = new Set(holes);
  for (let site = 0; site < length; site += 1) {
    let category;
    if (site >= 100 && site < 180 && !holeSet.has(site)) category = 0;
    else category = site % 2 ? 1 : 2;
    if (category === 0) {
      sequences[0][site] = 0;
      sequences[1][site] = 0;
      sequences[2][site] = 1;
    } else if (category === 1) {
      sequences[0][site] = 0;
      sequences[1][site] = 1;
      sequences[2][site] = 0;
    } else {
      sequences[0][site] = 1;
      sequences[1][site] = 0;
      sequences[2][site] = 0;
    }
  }
  return sequences;
}

function tripletFromCategories(categories) {
  const sequences = [new Uint8Array(categories.length), new Uint8Array(categories.length), new Uint8Array(categories.length)];
  categories.forEach((category, site) => {
    const bases = category === 0 ? [0, 0, 1] : category === 1 ? [0, 1, 0] : [1, 0, 0];
    for (let sequence = 0; sequence < 3; sequence += 1) sequences[sequence][site] = bases[sequence];
  });
  return sequences;
}

function naiveTrackZero(categories, gScale) {
  const mismatchSites = categories.filter((category) => category !== 0).length;
  const mismatchPenalty = gScale > 0 ? Math.floor(categories.length * gScale / mismatchSites) + 1 : 1;
  const runs = [];
  let start = 0;
  let positive = categories[0] === 0;
  for (let rank = 1; rank <= categories.length; rank += 1) {
    const nextPositive = rank < categories.length ? categories[rank] === 0 : !positive;
    if (rank < categories.length && nextPositive === positive) continue;
    runs.push({ start, end: rank - 1, score: positive ? rank - start : -(rank - start) });
    start = rank;
    positive = nextPositive;
  }
  let best = null;
  runs.forEach((run, runIndex) => {
    if (run.score <= 0) return;
    let score = 0;
    let highScore = 0;
    let highEnd = run.end;
    for (let extension = runIndex; extension < runs.length; extension += 1) {
      score += runs[extension].score > 0
        ? runs[extension].score
        : runs[extension].score * mismatchPenalty;
      if (score < 0) break;
      if (score >= highScore) {
        highScore = score;
        highEnd = runs[extension].end;
      }
    }
    const candidate = { start: run.start, end: highEnd + 1, score: highScore, runIndex };
    if (!best || candidate.score > best.score || (candidate.score === best.score && candidate.runIndex < best.runIndex)) {
      best = candidate;
    }
  });
  return { ...best, mismatchPenalty };
}

test("source GENECONV emits the six-track queue's strongest exact fragment", async () => {
  const instance = await engine();
  const result = scan(instance, sourceTriplet(), { gScale: 0 });
  assert.ok(result.total > 0);
  const fragment = result.rows.find((row) => row[0] === 0);
  assert.ok(fragment, "pair 0/1 should win the source tie order");
  assert.deepEqual(fragment.slice(1, 4), [0, 1, 2], "GCXoverD track 0 role mapping must be retained");
  assert.deepEqual(fragment.slice(4, 6), [100, 180]);
  assert.equal(fragment[7], 80);
  assert.equal(fragment[8], 600, "only this concrete triplet's variable complete sites enter the batch");
  assert.equal(fragment[9], 80);
  assert.equal(fragment[10], 520);
  assert.equal(fragment[11], 1);
  assert.ok(fragment[12] > 0, "the WASM row must carry the calibrated raw probability");
});

test("finite-G source scoring bridges short mismatch runs with the native penalty", async () => {
  const instance = await engine();
  const result = scan(instance, sourceTriplet(600, [139, 140]), { gScale: 1 });
  const fragment = result.rows.find((row) => row[0] === 0 && row[4] === 100 && row[5] === 180);
  assert.ok(fragment, "the two positive fragments should be joined across the short mismatch run");
  assert.equal(fragment[11], 2, "Int(L*G/NDiff)+1 must set the RDP5 mismatch penalty");
  assert.equal(fragment[7], 74, "78 matches - 2 mismatches * penalty 2");
});

test("packed GENECONV is exactly equivalent to the byte oracle and ignores decoys", async () => {
  const instance = await engine();
  const triplet = sourceTriplet(617, [139, 140]);
  const decoyA = Uint8Array.from({ length: 617 }, (_, site) => (site * 3 + 1) % 4);
  const decoyB = new Uint8Array(617).fill(3);
  const byte = scan(instance, [...triplet, decoyA], { gScale: 1, packedScan: false });
  const packedA = scan(instance, [...triplet, decoyA], { gScale: 1, packedScan: true });
  const packedB = scan(instance, [...triplet, decoyB], { gScale: 1, packedScan: true });
  assert.deepEqual(packedA, byte);
  assert.deepEqual(packedB, byte, "alignment-global variable columns must not change a concrete triplet scan");
});

test("all-different triplet sites remain in GENECONV compression", async () => {
  const instance = await engine();
  const sequences = sourceTriplet();
  for (let site = 220; site < 270; site += 1) {
    sequences[0][site] = 0;
    sequences[1][site] = 1;
    sequences[2][site] = 2;
  }
  const result = scan(instance, sequences, { gScale: 0, packedScan: true });
  assert.ok(result.rows.length > 0);
  assert.ok(result.rows.every((row) => row[8] === 600));
});

test("linear excursion index is differential-tested against the native quadratic extension rule", async () => {
  const instance = await engine();
  let state = 0x5a17c0de;
  for (let fixture = 0; fixture < 32; fixture += 1) {
    const categories = Array.from({ length: 600 }, (_, site) => {
      state ^= state << 13; state ^= state >>> 17; state ^= state << 5;
      if (site >= 100 && site < 300) return (state >>> 0) % 5 === 0 ? 1 + ((state >>> 8) & 1) : 0;
      return site & 1 ? 1 : 2;
    });
    const expected = naiveTrackZero(categories, 1);
    const result = scan(instance, tripletFromCategories(categories), { gScale: 1, packedScan: true });
    const actual = result.rows.find((row) => row[0] === 0);
    assert.ok(actual, `fixture ${fixture} must retain the dominant track-zero fragment`);
    assert.deepEqual(
      [actual[4], actual[5], actual[7], actual[11]],
      [expected.start, expected.end, expected.score, expected.mismatchPenalty],
      `fixture ${fixture} must preserve GetMaxFragScoreP endpoint and score ties`,
    );
  }
});
