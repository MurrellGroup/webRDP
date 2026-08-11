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
  const statsPtr = outPtr + 96;
  const wordsPerSequence = Math.ceil(nSites / 16);
  const packedPtr = align(statsPtr + 160);
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
  return { seqPtr, distancePtr, prefixAPtr, prefixBPtr, outPtr, statsPtr, packedPtr, validityPtr, poolPtr, nearestIndexesPtr, nearestDistancesPtr, cohortPtr, tractMaskPtr, backgroundMaskPtr, dmaxOutPtr, wordsPerSequence };
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
    127,
    pointers.prefixAPtr,
    pointers.prefixBPtr,
    pointers.statsPtr,
    0,
  );
  const stats = new Int32Array(instance.exports.memory.buffer, pointers.statsPtr, 33);
  assert.equal(stats[0], 120, "GENECONV G=0 run should span the inserted tract");
  assert.ok(stats[5] / stats[6] > 0.9, "window topology should switch with the known tract");
  assert.ok(Math.abs(stats[17] - 90) <= 2, "left breakpoint polishing should remain at the truth");
  assert.ok(Math.abs(stats[18] - 210) <= 2, "right breakpoint polishing should remain at the truth");
  assert.equal(stats[11], 120, "3Seq maximum HGRW descent should span the inserted tract");
  assert.equal(stats[19], 180, "3Seq should count major-parent matches");
  assert.equal(stats[20], 120, "3Seq should count minor-parent matches");
  assert.equal(stats[22], 300, "three representative windows should each contribute 100 decisive bootstraps");
  assert.equal(stats[21], 300, "all bootstrap topologies should match the known mosaic");
  assert.ok(Math.abs(stats[25] - 90) <= 10, `MAXCHI left peak should localize the true boundary, got ${stats[25]}`);
  assert.ok(Math.abs(stats[26] - 210) <= 10, `MAXCHI right peak should localize the true boundary, got ${stats[26]}`);
  assert.ok(Math.abs(stats[27] - 90) <= 10, `CHIMAERA left peak should localize the true boundary, got ${stats[27]}`);
  assert.ok(Math.abs(stats[28] - 210) <= 10, `CHIMAERA right peak should localize the true boundary, got ${stats[28]}`);
  assert.ok(Math.abs(stats[29] - 90) <= 60, `BOOTSCAN run should cover the left boundary, got ${stats[29]}`);
  assert.ok(Math.abs(stats[30] - 210) <= 60, `BOOTSCAN run should cover the right boundary, got ${stats[30]}`);
  assert.ok(Math.abs(stats[31] - 90) <= 60, `SISCAN run should cover the left boundary, got ${stats[31]}`);
  assert.ok(Math.abs(stats[32] - 210) <= 60, `SISCAN run should cover the right boundary, got ${stats[32]}`);

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
    25,
    12345,
    2,
    pointers.prefixAPtr,
    pointers.prefixBPtr,
    pointers.statsPtr,
    0,
  );
  const bootscanOnly = new Int32Array(instance.exports.memory.buffer, pointers.statsPtr, 23);
  assert.equal(bootscanOnly[0], 0, "disabled GENECONV kernel should not run");
  assert.ok(bootscanOnly[6] > 0, "enabled BootScan kernel should run");
  assert.equal(bootscanOnly[7], 0, "disabled MaxChi kernel should not run");
  assert.equal(bootscanOnly[11], 0, "disabled 3SEQ kernel should not run");
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

test("MAXCHI and CHIMAERA windows remain in informative-site coordinates", async () => {
  const instance = await engine();
  const nSeq = 3;
  const nSites = 600;
  const pointers = reserve(instance, nSeq, nSites);
  const sequences = new Uint8Array(instance.exports.memory.buffer, pointers.seqPtr, nSeq * nSites);
  sequences.fill(0);
  // Only every tenth nucleotide is polymorphic. The recombinant changes from
  // the major to the minor state for informative ranks 20..39. Adding the 540
  // monomorphic columns must not change a 20-VNP MAXCHI/CHIMAERA window.
  for (let rank = 0; rank < 60; rank += 1) {
    const site = rank * 10 + 5;
    sequences[nSites + site] = 1;
    sequences[nSites * 2 + site] = rank >= 20 && rank < 40 ? 1 : 0;
  }
  instance.exports.method_stats(
    pointers.seqPtr,
    nSites,
    2,
    0,
    1,
    205,
    405,
    20,
    1,
    0,
    99,
    4 | 8,
    pointers.prefixAPtr,
    pointers.prefixBPtr,
    pointers.statsPtr,
    0,
  );
  const stats = new Int32Array(instance.exports.memory.buffer, pointers.statsPtr, 35);
  assert.ok(Math.abs(stats[25] - 205) <= 12, `MAXCHI left VNP peak should map near 205, got ${stats[25]}`);
  assert.ok(Math.abs(stats[26] - 405) <= 12, `MAXCHI right VNP peak should map near 405, got ${stats[26]}`);
  assert.ok(Math.abs(stats[27] - 205) <= 12, `CHIMAERA left VNP peak should map near 205, got ${stats[27]}`);
  assert.ok(Math.abs(stats[28] - 405) <= 12, `CHIMAERA right VNP peak should map near 405, got ${stats[28]}`);
});

test("GENECONV finite G-scale bridges source-scored internal mismatches", async () => {
  const instance = await engine();
  const nSeq = 3;
  const nSites = 100;
  const pointers = reserve(instance, nSeq, nSites);
  const sequences = new Uint8Array(instance.exports.memory.buffer, pointers.seqPtr, nSeq * nSites);
  sequences.fill(0, 0, nSites);
  sequences.fill(1, nSites, nSites * 2);
  sequences.fill(0, nSites * 2, nSites * 3);
  sequences.fill(1, nSites * 2 + 30, nSites * 2 + 70);
  sequences[nSites * 2 + 50] = 2;
  const run = (gScale) => {
    instance.exports.method_stats(
      pointers.seqPtr,
      nSites,
      2,
      0,
      1,
      30,
      70,
      20,
      1,
      0,
      7,
      1,
      pointers.prefixAPtr,
      pointers.prefixBPtr,
      pointers.statsPtr,
      gScale,
    );
    return new Int32Array(instance.exports.memory.buffer, pointers.statsPtr, 35).slice();
  };
  const exactRuns = run(0);
  const mismatchTolerant = run(1);
  assert.equal(exactRuns[0], 20, "G=0 must stop at the internal discordance");
  assert.equal(mismatchTolerant[0], 37, "G=1 must apply the source integer mismatch penalty");
  assert.equal(mismatchTolerant[3], 30);
  assert.equal(mismatchTolerant[4], 70);
});

test("RDP5 source detector localizes its own VNP-window event and polarity", async () => {
  const instance = await engine();
  const nSeq = 3;
  const nSites = 300;
  const pointers = reserve(instance, nSeq, nSites);
  const sequences = new Uint8Array(instance.exports.memory.buffer, pointers.seqPtr, nSeq * nSites);
  sequences.fill(0, 0, nSites); // major parent
  sequences.fill(1, nSites, nSites * 2); // minor parent
  sequences.fill(0, nSites * 2, nSites * 3); // recombinant background
  sequences.fill(1, nSites * 2 + 90, nSites * 2 + 210);
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
  assert.equal(result[2], 2, "source polarity should identify the requested recombinant");
  assert.equal(result[3], 0);
  assert.equal(result[4], 1);
  assert.ok(Math.abs(result[0] - 90) <= 20, `unexpected source start ${result[0]}`);
  assert.ok(Math.abs(result[1] - 210) <= 20, `unexpected source end ${result[1]}`);
  assert.equal(result[16], 30);
  assert.equal(result[17], 1);
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
});
