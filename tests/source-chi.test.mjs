import assert from "node:assert/strict";
import fs from "node:fs";
import test from "node:test";

const ROW_INTS = 16;

async function engine() {
  const bytes = fs.readFileSync(new URL("../public/wasm/rdp.wasm", import.meta.url));
  return (await WebAssembly.instantiate(bytes)).instance;
}

function scan(instance, sequences, fullWindow = 40, methodMask = 12, packedScan = false) {
  const nSites = sequences[0].length;
  const align = (value, multiple = 16) => Math.ceil(value / multiple) * multiple;
  const seqPtr = 65_536;
  const wordsPerSequence = Math.ceil(nSites / 16);
  const packedPtr = align(seqPtr + sequences.length * nSites);
  const validityPtr = packedPtr + sequences.length * wordsPerSequence * 4;
  const positionsPtr = align(validityPtr + sequences.length * wordsPerSequence * 4);
  const scoresPtr = positionsPtr + 4 * (nSites + 1) * 4;
  const missingPtr = align(scoresPtr + 4 * (nSites + 1), 4);
  const chiPtr = align(missingPtr + (nSites + 1) * 4, 8);
  const smoothPtr = chiPtr + (nSites + 1) * 8;
  const peakCapacity = 64;
  const peakPtr = align(smoothPtr + (nSites + 1) * 8, 4);
  const outputCapacity = 128;
  const outPtr = peakPtr + peakCapacity * 6 * 4;
  const required = outPtr + outputCapacity * ROW_INTS * 4;
  const missingPages = Math.ceil((required - instance.exports.memory.buffer.byteLength) / 65_536);
  if (missingPages > 0) instance.exports.memory.grow(missingPages);
  const encoded = new Uint8Array(sequences.length * nSites);
  sequences.forEach((sequence, index) => encoded.set(sequence, index * nSites));
  new Uint8Array(instance.exports.memory.buffer, seqPtr, encoded.length).set(encoded);
  const packed = new Uint32Array(instance.exports.memory.buffer, packedPtr, sequences.length * wordsPerSequence);
  const validity = new Uint32Array(instance.exports.memory.buffer, validityPtr, sequences.length * wordsPerSequence);
  packed.fill(0);
  validity.fill(0);
  for (let sequence = 0; sequence < sequences.length; sequence += 1) {
    for (let site = 0; site < nSites; site += 1) {
      const base = encoded[sequence * nSites + site];
      if (base >= 4) continue;
      const word = sequence * wordsPerSequence + (site >>> 4);
      const shift = (site & 15) * 2;
      packed[word] |= base << shift;
      validity[word] |= 1 << shift;
    }
  }
  const sharedArguments = [
    nSites, 0, 1, 2, fullWindow, 0, methodMask, positionsPtr, scoresPtr,
    missingPtr, chiPtr, smoothPtr, peakPtr, peakCapacity, outPtr, outputCapacity,
  ];
  const total = packedScan
    ? instance.exports.scan_source_chi_all_packed(packedPtr, validityPtr, wordsPerSequence, ...sharedArguments)
    : instance.exports.scan_source_chi_all(seqPtr, ...sharedArguments);
  const rows = Array.from({ length: Math.min(total, outputCapacity) }, (_, index) => (
    Array.from(new Int32Array(instance.exports.memory.buffer, outPtr + index * ROW_INTS * 4, ROW_INTS))
  ));
  return { rows, total };
}

function twoTractTriplet(length = 600) {
  const recombinant = new Uint8Array(length);
  const major = new Uint8Array(length);
  const minor = new Uint8Array(length);
  for (let site = 0; site < length; site += 1) {
    major[site] = site % 4;
    minor[site] = (site + 1) % 4;
    const inside = (site >= 100 && site < 180) || (site >= 350 && site < 450);
    recombinant[site] = inside ? minor[site] : major[site];
  }
  return [recombinant, major, minor];
}

test("source MAXCHI/CHIMAERA scanner retains two disjoint peak pairs", async () => {
  const instance = await engine();
  const first = scan(instance, twoTractTriplet());
  const second = scan(instance, twoTractTriplet());
  assert.deepEqual(second, first, "source peak discovery must be deterministic");

  const maxChi = first.rows.filter((row) => row[0] === 3);
  const chimaera = first.rows.filter((row) => row[0] === 4 && row[1] === 0);
  for (const rows of [maxChi, chimaera]) {
    assert.ok(rows.some((row) => row[3] === 100 && row[4] === 180), "first mosaic tract must be retained");
    assert.ok(rows.some((row) => row[3] === 350 && row[4] === 450), "second mosaic tract must be retained");
  }
  assert.ok(new Set(maxChi.map((row) => row[2])).size >= 2, "MAXCHI must evaluate separate pair-equality tracks");
  assert.ok(chimaera.every((row) => row[9] === 600 && row[10] === 20), "CHIMAERA must report compressed-site and source half-window metadata");
  assert.ok(first.rows.every((row) => row[6] === Math.min(row[7], row[8])), "a paired call uses its weaker grown boundary statistic");
  assert.ok(first.rows.every((row) => row[13] >= row[10] && row[14] >= row[10]), "GrowMChiWin must never shrink below the scan half-window");
});

test("MAXCHI keeps all-different triplet sites while CHIMAERA strips them", async () => {
  const instance = await engine();
  const sequences = twoTractTriplet();
  for (let site = 220; site < 270; site += 1) sequences[0][site] = (sequences[2][site] + 1) % 4;
  const { rows } = scan(instance, sequences);
  const maxChiInformative = rows.find((row) => row[0] === 3)?.[9];
  const chimaeraInformative = rows.find((row) => row[0] === 4 && row[1] === 0)?.[9];
  assert.equal(maxChiInformative, 600, "FSSMC-compatible compression includes all-different sites");
  assert.equal(chimaeraInformative, 550, "FSSRDP-compatible compression retains only exact binary parent matches");
});

test("packed production scan exactly matches the byte validation oracle", async () => {
  const instance = await engine();
  const sequences = twoTractTriplet(617);
  for (let site = 17; site < sequences[0].length; site += 97) sequences[site % 3][site] = 4;
  const byteResult = scan(instance, sequences, 42, 12, false);
  const packedResult = scan(instance, sequences, 42, 12, true);
  assert.deepEqual(packedResult, byteResult);
});

test("source chi compression is local to the concrete triplet", async () => {
  const instance = await engine();
  const triplet = twoTractTriplet(617);
  const firstDecoy = Uint8Array.from({ length: 617 }, (_, site) => (site * 3 + 1) % 4);
  const secondDecoy = new Uint8Array(617).fill(3);
  const first = scan(instance, [...triplet, firstDecoy], 42, 12, true);
  const second = scan(instance, [...triplet, secondDecoy], 42, 12, true);
  assert.deepEqual(second, first, "alignment-global variable sites from an unrelated sequence must not enter MAXCHI/CHIMAERA tracks");
});
