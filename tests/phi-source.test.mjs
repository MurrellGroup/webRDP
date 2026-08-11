import assert from "node:assert/strict";
import test from "node:test";

import {
  sourcePairIncompatibility,
  sourcePhiAnalyticMeanVariance,
  sourcePhiStatistic,
  sourcePhiTest,
} from "../public/rdp-phi.js";

function encodeRows(rows) {
  const map = new Map([["A", 0], ["C", 1], ["G", 2], ["T", 3]]);
  return Uint8Array.from(rows.flatMap((row) => [...row].map((base) => map.get(base) ?? 4)));
}

test("RDP5 pair_score graph distinguishes a four-gamete cycle", () => {
  const incompatible = encodeRows(["AA", "AC", "CA", "CC"]);
  assert.equal(sourcePairIncompatibility(incompatible, 4, 2, 0, 1), 1);

  const compatible = encodeRows(["AA", "AA", "CA", "CC"]);
  assert.equal(sourcePairIncompatibility(compatible, 4, 2, 0, 1), 0);
});

test("PHI and analytic moments follow the supplied RDP5 matrix formulas", () => {
  const matrix = Uint8Array.from([
    0, 0, 1, 1, 0, 1,
    0, 0, 0, 1, 1, 0,
    1, 0, 0, 0, 1, 1,
    1, 1, 0, 0, 0, 1,
    0, 1, 1, 0, 0, 0,
    1, 0, 1, 1, 0, 0,
  ]);
  assert.equal(sourcePhiStatistic(matrix, 6, 2), 4 / 9);
  const moments = sourcePhiAnalyticMeanVariance(matrix, 6, 2);
  assert.ok(Math.abs(moments.mean - 8 / 15) < 1e-12);
  assert.ok(Math.abs(moments.variance - 0.012949245541838135) < 1e-15);
  assert.ok(Number.isFinite(moments.variance));
  assert.ok(moments.variance >= 0);
});

test("source PHI detects spatially clustered compatible genealogies reproducibly", () => {
  const partitions = [
    [0, 0, 0, 0, 1, 1, 1, 1],
    [0, 0, 1, 1, 0, 0, 1, 1],
    [0, 1, 0, 1, 0, 1, 0, 1],
  ];
  const sequenceCount = 8;
  const length = 180;
  const encoded = new Uint8Array(sequenceCount * length);
  for (let site = 0; site < length; site += 1) {
    const partition = partitions[Math.floor(site / 60)];
    for (let sequence = 0; sequence < sequenceCount; sequence += 1) {
      encoded[sequence * length + site] = partition[sequence];
    }
  }
  const full = sourcePhiTest(encoded, sequenceCount, length, { window: 20, maxInformativeSites: 384 });
  assert.equal(full.informativeSites, 180);
  assert.equal(full.subsampled, false);
  assert.equal(full.validNormalApproximation, true);
  assert.ok(full.pValue < 0.01, `expected clustered incompatibilities to be significant, got ${full.pValue}`);
  assert.deepEqual(sourcePhiTest(encoded, sequenceCount, length, { window: 20, maxInformativeSites: 384 }), full);

  const bounded = sourcePhiTest(encoded, sequenceCount, length, { window: 20, maxInformativeSites: 64 });
  assert.equal(bounded.informativeSites, 64);
  assert.equal(bounded.totalInformativeSites, 180);
  assert.equal(bounded.subsampled, true);
  assert.match(bounded.compatibility, /PHITest2/);
});
