import assert from "node:assert/strict";
import test from "node:test";
import { binomialUpper, chiSquareP, geneconvSourceG0Probability, geneconvSourceProbability, methodEvidence, rdp5SourceProbability, sourceChiWindowProbability, threeSeqExactP } from "../public/rdp-statistics.js";

test("exact binomial upper tail matches closed-form small cases", () => {
  assert.ok(Math.abs(binomialUpper(5, 10, 0.5) - 0.623046875) < 1e-12);
  assert.ok(Math.abs(binomialUpper(10, 10, 0.5) - 0.0009765625) < 1e-12);
  assert.equal(binomialUpper(0, 10, 0.5), 1);
  assert.ok(binomialUpper(900, 1_000, 0.5) > 0);
});

test("RDP5 source calibration retains its y-1 scan multiplier", () => {
  const source = {
    tractSites: 10,
    common: 8,
    mediumSites: 20,
    informativeSites: 100,
    probabilitySites: 99,
  };
  const expected = binomialUpper(8, 10, 0.2) * 9.9;
  assert.ok(Math.abs(rdp5SourceProbability(source) - expected) < 1e-14);
  assert.notEqual(rdp5SourceProbability(source), binomialUpper(8, 10, 0.2) * 10);
});

test("GENECONV G=0 calibration is the RDP5 CalcKMax/GCCalc specialization", () => {
  const run = 12;
  const sites = 40;
  const matches = 30;
  const expected = 1 - Math.exp(-10 * Math.pow(0.75, run));
  assert.ok(Math.abs(geneconvSourceG0Probability(run, sites, matches) - expected) < 1e-14);
  assert.equal(geneconvSourceG0Probability(3, sites, matches), 1, "RDP5 requires a score above three");
  assert.equal(geneconvSourceG0Probability(run, sites, sites), 1, "an identical pair is not a valid source fragment family");
});

test("GENECONV finite G-scale follows the RDP5 lambda/K root calculation", () => {
  assert.ok(Math.abs(geneconvSourceProbability(30, 100, 70, 1) - 0.00853636746733932) < 1e-14);
  assert.equal(geneconvSourceProbability(30, 100, 70, 0), geneconvSourceG0Probability(30, 100, 70));
  assert.equal(geneconvSourceProbability(20, 100, 90, 0.1), 1, "a non-negative expected fragment score has no finite KA root");
});

test("MAXCHI/CHIMAERA calibration uses the RDP5 informative half-window multiplier", () => {
  const statistic = 18.5;
  const expected = chiSquareP(statistic) * (240 / 50) * 3;
  assert.ok(Math.abs(sourceChiWindowProbability(statistic, 240, 100) - expected) < 1e-14);
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

test("a method cannot confirm a distant candidate without a co-located signal", () => {
  const evidence = methodEvidence({
    insideMinor: 20,
    insideMajor: 2,
    outsideMinor: 2,
    outsideMajor: 20,
    methodSignals: [{ method: "RDP" }],
  }, {
    genconvRun: 20,
    genconvEligible: 40,
    genconvMatches: 30,
    bootscanConsistent: 0,
    bootscanWindows: 0,
    bootscanBootstrapConsistent: 0,
    bootscanBootstrapReplicates: 0,
    maxChi: 0,
    chimaera: 0,
    siskanScore: 0,
    siskanSites: 0,
    threeSeqDescent: 0,
    threeSeqSites: 0,
    threeSeqMajorSites: 0,
    threeSeqMinorSites: 0,
  }, {
    methods: ["RDP", "GENECONV"],
    window: 20,
    step: 1,
    correction: "none",
    alpha: 0.05,
  }, 1, 100);
  const geneconv = evidence.find((row) => row.method === "GENECONV");
  assert.equal(geneconv?.supported, false);
  assert.equal(geneconv?.pValue, 1);
  assert.equal(geneconv?.calibration, "no co-located discovery signal");
});

test("BootScan never falls back to the retired pairwise/window approximation", () => {
  const evidence = methodEvidence({
    insideMinor: 95,
    insideMajor: 5,
    outsideMinor: 8,
    outsideMajor: 92,
  }, {
    genconvRun: 0,
    genconvEligible: 200,
    genconvMatches: 100,
    bootscanConsistent: 38,
    bootscanWindows: 40,
    bootscanBootstrapConsistent: 286,
    bootscanBootstrapReplicates: 300,
    maxChi: 0,
    chimaera: 0,
    siskanScore: 0,
    siskanSites: 0,
    threeSeqDescent: 0,
    threeSeqSites: 0,
    threeSeqMajorSites: 0,
    threeSeqMinorSites: 0,
  }, {
    methods: ["BootScan"],
    window: 120,
    step: 5,
    correction: "none",
    alpha: 0.05,
  }, 1, 2_400);
  assert.equal(evidence[0]?.pValue, 1);
  assert.equal(evidence[0]?.statistic, 0);
  assert.match(evidence[0]?.calibration ?? "", /source BootScan batch signal unavailable/);
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
