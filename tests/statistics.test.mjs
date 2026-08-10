import assert from "node:assert/strict";
import test from "node:test";
import { binomialUpper, methodEvidence, threeSeqExactP } from "../public/rdp-statistics.js";

test("exact binomial upper tail matches closed-form small cases", () => {
  assert.ok(Math.abs(binomialUpper(5, 10, 0.5) - 0.623046875) < 1e-12);
  assert.ok(Math.abs(binomialUpper(10, 10, 0.5) - 0.0009765625) < 1e-12);
  assert.equal(binomialUpper(0, 10, 0.5), 1);
  assert.ok(binomialUpper(900, 1_000, 0.5) > 0);
});

test("method families retain independent statistics and calibrations", () => {
  const candidate = {
    insideMinor: 95,
    insideMajor: 5,
    outsideMinor: 8,
    outsideMajor: 92,
  };
  const stats = {
    genconvRun: 28,
    genconvEligible: 200,
    genconvMatches: 100,
    bootscanConsistent: 38,
    bootscanWindows: 40,
    bootscanBootstrapConsistent: 286,
    bootscanBootstrapReplicates: 300,
    maxChi: 31,
    chimaera: 28,
    siskanScore: 72,
    siskanSites: 200,
    threeSeqDescent: 48,
    threeSeqSites: 200,
    threeSeqMajorSites: 100,
    threeSeqMinorSites: 100,
  };
  const evidence = methodEvidence(candidate, stats, {
    methods: ["RDP", "GENECONV", "BootScan", "MaxChi", "Chimaera", "SiScan", "3Seq"],
    window: 120,
    step: 5,
    correction: "none",
    alpha: 0.05,
  }, 1, 2_400);
  assert.equal(evidence.length, 7);
  assert.equal(new Set(evidence.map((item) => item.statisticLabel)).size, 6);
  assert.equal(new Set(evidence.map((item) => item.calibration)).size, 7);
  assert.ok(new Set(evidence.map((item) => item.pValue.toPrecision(8))).size >= 5);
});

function bruteForceThreeSeq(up, down, threshold) {
  let paths = 0;
  let hits = 0;
  function visit(usedUp, usedDown, walk, maximum, descent) {
    if (usedUp === up && usedDown === down) {
      paths += 1;
      if (descent >= threshold) hits += 1;
      return;
    }
    if (usedUp < up) {
      const nextWalk = walk + 1;
      visit(usedUp + 1, usedDown, nextWalk, Math.max(maximum, nextWalk), descent);
    }
    if (usedDown < down) {
      const nextWalk = walk - 1;
      visit(usedUp, usedDown + 1, nextWalk, maximum, Math.max(descent, maximum - nextWalk));
    }
  }
  visit(0, 0, 0, 0, 0);
  return hits / paths;
}

test("3SEQ exact HGRW tail agrees with exhaustive small-path enumeration", () => {
  for (let up = 1; up <= 4; up += 1) {
    for (let down = 1; down <= 4; down += 1) {
      for (let threshold = 1; threshold <= down; threshold += 1) {
        const result = threeSeqExactP(up, down, threshold);
        assert.equal(result.exact, true);
        assert.ok(Math.abs(result.p - bruteForceThreeSeq(up, down, threshold)) < 1e-12);
      }
    }
  }
});

test("3SEQ exact HGRW tail reproduces the published 30/30 example", () => {
  const result = threeSeqExactP(30, 30, 18);
  assert.equal(result.exact, true);
  assert.ok(result.p > 1.6e-4 && result.p < 1.9e-4, `unexpected p=${result.p}`);
});

test("3SEQ calibration falls back safely when the bounded DP would be too large", () => {
  assert.deepEqual(threeSeqExactP(1_000, 1_000, 400, 10), { p: null, exact: false });
});
