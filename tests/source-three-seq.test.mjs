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

function scan(instance, sequences, packedMode = false, sourceMode = false, circular = false) {
  const nSites = sequences[0].length;
  const seqPtr = 65_536;
  const align = (value, multiple = 16) => Math.ceil(value / multiple) * multiple;
  const packedData = pack(sequences);
  const packedPtr = align(seqPtr + sequences.length * nSites);
  const validityPtr = packedPtr + packedData.packed.byteLength;
  const workspacePtr = align(validityPtr + packedData.validity.byteLength);
  const workspaceBytes = sourceMode ? instance.exports.source_three_seq_workspace_bytes(nSites) : 0;
  const outPtr = align(workspacePtr + workspaceBytes);
  const required = outPtr + 6 * ROW_INTS * 4;
  const missingPages = Math.ceil((required - instance.exports.memory.buffer.byteLength) / 65_536);
  if (missingPages > 0) instance.exports.memory.grow(missingPages);
  const encoded = new Uint8Array(sequences.length * nSites);
  sequences.forEach((sequence, index) => encoded.set(sequence, index * nSites));
  new Uint8Array(instance.exports.memory.buffer, seqPtr, encoded.length).set(encoded);
  new Uint32Array(instance.exports.memory.buffer, packedPtr, packedData.packed.length).set(packedData.packed);
  new Uint32Array(instance.exports.memory.buffer, validityPtr, packedData.validity.length).set(packedData.validity);
  const total = sourceMode
    ? packedMode
      ? instance.exports.scan_source_three_seq_triplet_packed_mode(
          packedPtr,
          validityPtr,
          packedData.wordsPerSequence,
          nSites,
          0,
          1,
          2,
          circular ? 1 : 0,
          workspacePtr,
          outPtr,
        )
      : instance.exports.scan_source_three_seq_triplet_mode(
          seqPtr,
          nSites,
          0,
          1,
          2,
          circular ? 1 : 0,
          workspacePtr,
          outPtr,
        )
    : packedMode
      ? instance.exports.scan_source_three_seq_triplet_packed(
        packedPtr,
        validityPtr,
        packedData.wordsPerSequence,
        nSites,
        0,
        1,
        2,
        outPtr,
        )
      : instance.exports.scan_source_three_seq_triplet(seqPtr, nSites, 0, 1, 2, outPtr);
  return Array.from({ length: total }, (_, index) => (
    Array.from(new Int32Array(instance.exports.memory.buffer, outPtr + index * ROW_INTS * 4, ROW_INTS))
  ));
}

function mosaicTriplet(length = 240) {
  const recombinant = new Uint8Array(length);
  const major = new Uint8Array(length);
  const minor = new Uint8Array(length);
  for (let site = 0; site < length; site += 1) {
    major[site] = site % 4;
    minor[site] = (site + 1) % 4;
    recombinant[site] = site >= 80 && site < 160 ? minor[site] : major[site];
  }
  return [recombinant, major, minor];
}

// Independent, deliberately scalar transcription of the author source's
// FindSubSeqTS + CheckwrapC control flow. CheckwrapC is not an exhaustive
// circular maximum-subarray search: it extends the source-selected excursion
// only through the earlier of its peak and trough ranks, with strict ties.
function sourceCheckwrapOracle(sequences, row, circular) {
  const [target, plusParent, minusParent] = row;
  const positions = [];
  const steps = [];
  for (let site = 0; site < sequences[0].length; site += 1) {
    const targetBase = sequences[target][site];
    const plus = sequences[plusParent][site];
    const minus = sequences[minusParent][site];
    if (targetBase >= 4 || plus >= 4 || minus >= 4 || plus === minus) continue;
    if (targetBase === plus) { positions.push(site); steps.push(1); }
    else if (targetBase === minus) { positions.push(site); steps.push(-1); }
  }
  const heights = [];
  let height = 0;
  let maximum = 0;
  let best = 0;
  let peakBoundary = 0;
  let endRank = 0;
  for (let rank = 0; rank < steps.length; rank += 1) {
    height += steps[rank];
    heights.push(height);
    if (height > maximum) {
      maximum = height;
      peakBoundary = rank + 1;
    }
    if (maximum - height > best) {
      best = maximum - height;
      endRank = rank;
    }
  }

  // CheckwrapC appends this bounded prefix to the end of the same walk.
  maximum = peakBoundary > 0 ? heights[peakBoundary - 1] : 0;
  const total = heights.at(-1);
  const limit = Math.min(peakBoundary > 0 ? peakBoundary - 1 : 0, endRank);
  if (peakBoundary > 0) {
    for (let rank = 0; rank <= limit; rank += 1) {
      const appendedHeight = total + heights[rank];
      if (appendedHeight > maximum) {
        maximum = appendedHeight;
        peakBoundary = rank + 1;
      }
      if (maximum - appendedHeight > best) {
        best = maximum - appendedHeight;
        endRank = rank;
      }
    }
  }

  let start = positions[peakBoundary % steps.length];
  let end = positions[endRank] + 1;
  const sourceWrap = start > end;
  if (sourceWrap && !circular) {
    start = positions[(endRank + 1) % steps.length];
    end = positions[(peakBoundary + steps.length - 1) % steps.length] + 1;
  }
  return { start, end, best, sourceWrap, linearComplement: sourceWrap && !circular };
}

test("source 3Seq screens all three recombinant roles in one concrete-triplet call", async () => {
  const instance = await engine();
  const rows = scan(instance, mosaicTriplet());
  assert.ok(new Set(rows.map((row) => row[0])).size === 3, "all three TSXOver target cycles must be represented");
  const forward = rows.find((row) => row[0] === 0 && row[3] === 1);
  const reverse = rows.find((row) => row[0] === 0 && row[3] === -1);
  assert.ok(forward && reverse);
  assert.deepEqual(forward.slice(1, 10), [1, 2, 1, 80, 160, 160, 80, 80, 240]);
  assert.deepEqual(reverse.slice(1, 10), [2, 1, -1, 0, 80, 80, 160, 80, 240], "strict source tie order retains the first maximum ascent");
});

test("packed source 3Seq is exactly equivalent to the byte oracle and ignores decoys", async () => {
  const instance = await engine();
  const triplet = mosaicTriplet(257);
  const decoyA = Uint8Array.from({ length: 257 }, (_, site) => (site * 3 + 1) % 4);
  const decoyB = new Uint8Array(257).fill(3);
  const byteRows = scan(instance, [...triplet, decoyA], false);
  assert.deepEqual(scan(instance, [...triplet, decoyA], true), byteRows);
  assert.deepEqual(scan(instance, [...triplet, decoyB], true), byteRows, "unrelated alignment rows must not define informative sites");
});

test("source 3Seq strips invariant, incomplete, and neither-parent-match sites locally", async () => {
  const instance = await engine();
  const triplet = mosaicTriplet();
  for (let site = 20; site < 30; site += 1) {
    triplet[0][site] = 2;
    triplet[1][site] = 2;
    triplet[2][site] = 2;
  }
  for (let site = 100; site < 110; site += 1) triplet[0][site] = (triplet[2][site] + 1) % 4;
  for (let site = 180; site < 190; site += 1) triplet[1][site] = 4;
  const rows = scan(instance, triplet, true);
  const forward = rows.find((row) => row[0] === 0 && row[3] === 1);
  assert.ok(forward);
  assert.equal(forward[9], 210, "only sites where the target matches exactly one valid parent enter its walk");
  assert.equal(forward[7], 70, "neither-parent-match sites inside the tract must not count as -1 steps");
});

test("CheckwrapC mode retains a circular tract and reports its linear complement", async () => {
  const instance = await engine();
  const [recombinant, major, minor] = mosaicTriplet();
  for (let site = 0; site < recombinant.length; site += 1) {
    recombinant[site] = site < 40 || site >= 200 ? minor[site] : major[site];
  }
  const circularRows = scan(instance, [recombinant, major, minor], false, true, true);
  const circular = circularRows.find((row) => row[0] === 0 && row[1] === 1 && row[2] === 2);
  assert.ok(circular);
  assert.deepEqual(circular.slice(3, 13), [1, 200, 40, 160, 80, 80, 240, 1, 1, 0]);

  const linearRows = scan(instance, [recombinant, major, minor], false, true, false);
  const linear = linearRows.find((row) => row[0] === 0 && row[1] === 1 && row[2] === 2);
  assert.ok(linear);
  assert.deepEqual(linear.slice(3, 13), [1, 40, 200, 160, 80, 80, 240, 1, 1, 1]);
  assert.deepEqual(
    scan(instance, [recombinant, major, minor], true, true, true),
    circularRows,
    "packed and byte CheckwrapC paths must be exact",
  );
});

test("optimized CheckwrapC agrees with the scalar source-control-flow oracle", async () => {
  const instance = await engine();
  let state = 0x73e2a91d;
  const random = () => {
    state ^= state << 13; state ^= state >>> 17; state ^= state << 5;
    return state >>> 0;
  };
  for (let replicate = 0; replicate < 80; replicate += 1) {
    const sequences = Array.from({ length: 3 }, () => Uint8Array.from(
      { length: 23 },
      () => random() % 19 === 0 ? 4 : random() % 4,
    ));
    for (const circular of [false, true]) {
      const rows = scan(instance, sequences, false, true, circular);
      assert.deepEqual(scan(instance, sequences, true, true, circular), rows);
      for (const row of rows) {
        const expected = sourceCheckwrapOracle(sequences, row, circular);
        assert.deepEqual(
          [row[4], row[5], row[8], row[11] === 1, row[12] === 1],
          [expected.start, expected.end, expected.best, expected.sourceWrap, expected.linearComplement],
        );
      }
    }
  }
});
