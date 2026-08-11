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
  const wordsPerSequence = Math.ceil(nSites / 16);
  const packedPtr = align(outPtr + 4096);
  const validityPtr = packedPtr + nSeq * wordsPerSequence * 4;
  const poolPtr = align(validityPtr + nSeq * wordsPerSequence * 4);
  const nearestIndexesPtr = poolPtr + nSeq * 4;
  const nearestDistancesPtr = nearestIndexesPtr + nSeq * 4;
  const cohortPtr = align(nearestDistancesPtr + nSeq * 4);
  const tractMaskPtr = cohortPtr + nSeq * 4;
  const backgroundMaskPtr = tractMaskPtr + wordsPerSequence * 4;
  const dmaxOutPtr = align(backgroundMaskPtr + wordsPerSequence * 4, 8);
  const required = dmaxOutPtr + 40;
  const missing = Math.ceil((required - instance.exports.memory.buffer.byteLength) / 65536);
  if (missing > 0) instance.exports.memory.grow(missing);
  return { seqPtr, distancePtr, prefixAPtr, prefixBPtr, outPtr, packedPtr, validityPtr, poolPtr, nearestIndexesPtr, nearestDistancesPtr, cohortPtr, tractMaskPtr, backgroundMaskPtr, dmaxOutPtr, wordsPerSequence };
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

function quartetPatternScores(a, b, c, d) {
  if ([a, b, c, d].some((base) => base >= 4)) return [0, 0, 0];
  const split = (firstA, firstB, secondA, secondB) => {
    if (firstA === firstB && secondA === secondB && firstA !== secondA) return 2;
    if (firstA === firstB && firstA !== secondA && firstA !== secondB && secondA !== secondB) return 1;
    if (secondA === secondB && secondA !== firstA && secondA !== firstB && firstA !== firstB) return 1;
    return 0;
  };
  return [split(a, b, c, d), split(a, c, b, d), split(a, d, b, c)];
}

function scalarDmax(bytes, nSites, cohort, candidates, start, end) {
  const sums = [0, 0, 0];
  const counts = [0, 0, 0];
  for (let ai = 0; ai < cohort.length - 3; ai += 1) {
    for (let bi = ai + 1; bi < cohort.length - 2; bi += 1) {
      for (let ci = bi + 1; ci < cohort.length - 1; ci += 1) {
        for (let di = ci + 1; di < cohort.length; di += 1) {
          const quartet = [cohort[ai], cohort[bi], cohort[ci], cohort[di]];
          const included = candidates.map((candidate) => quartet.includes(candidate));
          if (!included.some(Boolean)) continue;
          const inside = [0, 0, 0];
          const outside = [0, 0, 0];
          for (let site = 0; site < nSites; site += 1) {
            const scores = quartetPatternScores(...quartet.map((sequence) => bytes[sequence * nSites + site]));
            const target = site >= start && site < end ? inside : outside;
            for (let split = 0; split < 3; split += 1) target[split] += scores[split];
          }
          const insideTotal = inside.reduce((total, value) => total + value, 0);
          const outsideTotal = outside.reduce((total, value) => total + value, 0);
          if (insideTotal === 0 || outsideTotal === 0) continue;
          const distance = inside.reduce((total, value, split) => (
            total + Math.abs(value / insideTotal - outside[split] / outsideTotal)
          ), 0);
          included.forEach((value, index) => {
            if (!value) return;
            sums[index] += distance;
            counts[index] += 1;
          });
        }
      }
    }
  }
  return { values: sums.map((sum, index) => counts[index] ? sum / counts[index] : 0), counts };
}

test("packed VisRD dMax kernel matches the source quartet table", async () => {
  const instance = await engine();
  const nSeq = 7;
  const nSites = 73;
  const pointers = reserve(instance, nSeq, nSites);
  const bytes = new Uint8Array(nSeq * nSites);
  for (let sequence = 0; sequence < nSeq; sequence += 1) {
    for (let site = 0; site < nSites; site += 1) {
      bytes[sequence * nSites + site] = (sequence * 3 + site * 5 + (site >>> 2)) % 5;
    }
  }
  const { packed, validity } = pack(bytes, nSeq, nSites);
  new Uint32Array(instance.exports.memory.buffer, pointers.packedPtr, packed.length).set(packed);
  new Uint32Array(instance.exports.memory.buffer, pointers.validityPtr, validity.length).set(validity);
  const cohort = [0, 1, 2, 3, 4, 5, 6];
  const candidates = [1, 3, 5];
  const start = 19;
  const end = 57;
  const tract = new Uint32Array(pointers.wordsPerSequence);
  const background = new Uint32Array(pointers.wordsPerSequence);
  for (let site = 0; site < nSites; site += 1) {
    const word = site >>> 4;
    const lane = 1 << ((site & 15) * 2);
    if (site >= start && site < end) tract[word] |= lane;
    else background[word] |= lane;
  }
  new Int32Array(instance.exports.memory.buffer, pointers.cohortPtr, cohort.length).set(cohort);
  new Uint32Array(instance.exports.memory.buffer, pointers.tractMaskPtr, tract.length).set(tract);
  new Uint32Array(instance.exports.memory.buffer, pointers.backgroundMaskPtr, background.length).set(background);
  instance.exports.dmax_visrd_packed(
    pointers.packedPtr,
    pointers.validityPtr,
    pointers.wordsPerSequence,
    pointers.cohortPtr,
    cohort.length,
    ...candidates,
    pointers.tractMaskPtr,
    pointers.backgroundMaskPtr,
    pointers.dmaxOutPtr,
  );
  const actual = [...new Float64Array(instance.exports.memory.buffer, pointers.dmaxOutPtr, 3)];
  const actualCounts = [...new Int32Array(instance.exports.memory.buffer, pointers.dmaxOutPtr + 24, 3)];
  const expected = scalarDmax(bytes, nSites, cohort, candidates, start, end);
  actual.forEach((value, index) => assert.ok(Math.abs(value - expected.values[index]) < 1e-12));
  assert.deepEqual(actualCounts, expected.counts);
});

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
  const maximum = instance.exports.maximum_packed_distance(
    pointers.packedPtr,
    pointers.validityPtr,
    nSeq,
    pointers.wordsPerSequence,
  );
  assert.equal(maximum, Math.max(...scalar));
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

test("triplet count kernel evaluates a known internal mosaic tract", async () => {
  const instance = await engine();
  const nSeq = 3;
  const nSites = 300;
  const pointers = reserve(instance, nSeq, nSites);
  const sequences = new Uint8Array(instance.exports.memory.buffer, pointers.seqPtr, nSeq * nSites);
  sequences.fill(0, 0, nSites);
  sequences.fill(1, nSites, nSites * 2);
  sequences.fill(0, nSites * 2, nSites * 3);
  sequences.fill(1, nSites * 2 + 90, nSites * 2 + 210);
  instance.exports.triplet_counts(
    pointers.seqPtr, nSites, 2, 0, 1, 90, 210, pointers.outPtr,
  );
  const counts = new Int32Array(instance.exports.memory.buffer, pointers.outPtr, 6);
  assert.deepEqual([...counts.slice(0, 5)], [300, 120, 0, 180, 0]);
  assert.ok(counts[5] > 250_000, `expected a strong tract/background contrast, got ${counts[5] / 1000}`);
});

test("RDP5 source detector localizes its own VNP-window event and polarity", async () => {
  const instance = await engine();
  const nSeq = 4;
  const nSites = 300;
  const pointers = reserve(instance, nSeq, nSites);
  const sequences = new Uint8Array(instance.exports.memory.buffer, pointers.seqPtr, nSeq * nSites);
  sequences.fill(0, 0, nSites); // major parent
  sequences.fill(1, nSites, nSites * 2); // minor parent
  sequences.fill(0, nSites * 2, nSites * 3); // recombinant background
  sequences.fill(1, nSites * 2 + 90, nSites * 2 + 210);
  // A fourth sequence makes many otherwise invariant alignment columns
  // globally variable. It must have no influence on the concrete 0/1/2
  // triplet's compressed VNP stream.
  for (let site = 0; site < nSites; site += 1) sequences[nSites * 3 + site] = (site * 3 + 1) % 4;
  // Supply the third pair-match category required by FastRecCheckP without
  // changing the true ancestry switch.
  for (let site = 12; site < nSites; site += 17) {
    if (site === 97 || site === 199) continue;
    sequences[site] = 2;
    sequences[nSites + site] = 2;
    sequences[nSites * 2 + site] = 3;
  }
  const found = instance.exports.scan_rdp5_triplet(
    pointers.seqPtr,
    nSites,
    2,
    0,
    1,
    30,
    pointers.prefixBPtr,
    pointers.prefixAPtr,
    pointers.outPtr,
  );
  assert.equal(found, 1);
  const result = new Int32Array(instance.exports.memory.buffer, pointers.outPtr, 18);
  let expectedTripletVnps = 0;
  for (let site = 0; site < nSites; site += 1) {
    const values = [sequences[site], sequences[nSites + site], sequences[nSites * 2 + site]];
    const pairMatches = Number(values[0] === values[1]) + Number(values[0] === values[2]) + Number(values[1] === values[2]);
    if (pairMatches === 1) expectedTripletVnps += 1;
  }
  assert.equal(result[2], 2, "source polarity should identify the requested recombinant");
  assert.equal(result[3], 0);
  assert.equal(result[4], 1);
  assert.ok(Math.abs(result[0] - 90) <= 20, `unexpected source start ${result[0]}`);
  assert.ok(Math.abs(result[1] - 210) <= 20, `unexpected source end ${result[1]}`);
  assert.equal(result[16], 30);
  assert.equal(result[17], 1);
  assert.equal(result[6], expectedTripletVnps, "RDP must compress to this triplet's VNPs and ignore its invariant columns");
  const capacity = 8;
  const bestPtr = pointers.outPtr + capacity * 72;
  const totalSignals = instance.exports.scan_rdp5_triplet_all(
    pointers.seqPtr,
    nSites,
    2,
    0,
    1,
    30,
    pointers.prefixBPtr,
    pointers.prefixAPtr,
    pointers.outPtr,
    capacity,
    bestPtr,
  );
  assert.ok(totalSignals >= 1, "multi-signal RDP export should retain the source excursion");
  const retained = Array.from({ length: Math.min(capacity, totalSignals) }, (_, index) => (
    new Int32Array(instance.exports.memory.buffer, pointers.outPtr + index * 72, 18).slice()
  ));
  assert.ok(retained.some((signal) => signal[2] === 2 && Math.abs(signal[0] - 90) <= 20 && Math.abs(signal[1] - 210) <= 20));
  const best = new Int32Array(instance.exports.memory.buffer, bestPtr, 18);
  assert.equal(best[5], Math.max(...retained.map((signal) => signal[5])));
  const permutedTotal = instance.exports.scan_rdp5_triplet_all(
    pointers.seqPtr,
    nSites,
    0,
    1,
    2,
    30,
    pointers.prefixBPtr,
    pointers.prefixAPtr,
    pointers.outPtr,
    capacity,
    bestPtr,
  );
  const permuted = Array.from({ length: Math.min(capacity, permutedTotal) }, (_, index) => (
    new Int32Array(instance.exports.memory.buffer, pointers.outPtr + index * 72, 18).slice()
  ));
  const normalize = (signals) => signals.map((signal) => [
    signal[2], signal[3], signal[4], signal[0], signal[1], signal[5], signal[12], signal[13], signal[14],
  ].join(":")).sort();
  assert.deepEqual(normalize(permuted), normalize(retained), "RDP all-signal output must be invariant to input triplet order for worker caching");
  const packedTriplet = pack(sequences, nSeq, nSites);
  new Uint32Array(instance.exports.memory.buffer, pointers.packedPtr, packedTriplet.packed.length).set(packedTriplet.packed);
  new Uint32Array(instance.exports.memory.buffer, pointers.validityPtr, packedTriplet.validity.length).set(packedTriplet.validity);
  const packedTotal = instance.exports.scan_rdp5_triplet_all_packed(
    pointers.packedPtr,
    pointers.validityPtr,
    pointers.wordsPerSequence,
    nSites,
    2,
    0,
    1,
    30,
    pointers.prefixBPtr,
    pointers.prefixAPtr,
    pointers.outPtr,
    capacity,
    bestPtr,
  );
  const packedSignals = Array.from({ length: Math.min(capacity, packedTotal) }, (_, index) => (
    new Int32Array(instance.exports.memory.buffer, pointers.outPtr + index * 72, 18).slice()
  ));
  assert.equal(packedTotal, totalSignals);
  assert.deepEqual(normalize(packedSignals), normalize(retained), "packed production RDP scan must exactly match its byte oracle");

  sequences.fill(3, nSites * 3, nSites * 4);
  const repacked = pack(sequences, nSeq, nSites);
  new Uint32Array(instance.exports.memory.buffer, pointers.packedPtr, repacked.packed.length).set(repacked.packed);
  new Uint32Array(instance.exports.memory.buffer, pointers.validityPtr, repacked.validity.length).set(repacked.validity);
  const decoyChangedTotal = instance.exports.scan_rdp5_triplet_all_packed(
    pointers.packedPtr,
    pointers.validityPtr,
    pointers.wordsPerSequence,
    nSites,
    2,
    0,
    1,
    30,
    pointers.prefixBPtr,
    pointers.prefixAPtr,
    pointers.outPtr,
    capacity,
    bestPtr,
  );
  const decoyChanged = Array.from({ length: Math.min(capacity, decoyChangedTotal) }, (_, index) => (
    new Int32Array(instance.exports.memory.buffer, pointers.outPtr + index * 72, 18).slice()
  ));
  assert.equal(decoyChangedTotal, packedTotal);
  assert.deepEqual(normalize(decoyChanged), normalize(packedSignals), "a non-triplet sequence must never leak its variable columns into RDP's VNP stream");
});
