import assert from "node:assert/strict";
import test from "node:test";

import {
  identifyRecombinantRoles,
  sourceVbClng,
  sourceDmaxScores,
  sourceFinalTrimParsimonyTest,
  sourceConflictScores,
  sourceHistoricalSetMembers,
  sourceOuCheckScores,
  sourceOuIndexScores,
  sourceParsimonyScores,
  sourcePhPrScores,
  sourceSimpleDist,
  sourceSetDistanceScores,
  sourceSsDistScores,
  sourceTrpScores,
} from "../public/rdp-recombinant-identification.js";

test("MakeConsensusC quantization uses VB6 banker's rounding", () => {
  assert.equal(sourceVbClng(2.5), 2);
  assert.equal(sourceVbClng(3.5), 4);
  assert.equal(sourceVbClng(-1.5), -2);
  assert.equal(sourceVbClng(-2.5), -2);
});

test("FinalTrim parsimony penalties are applied independently by polarity", () => {
  const result = sourceFinalTrimParsimonyTest([2, 0, 3], [4, 0, 1]);
  assert.deepEqual(result.values, [2, 0, 1]);
  assert.deepEqual(result.points, [-10, 0, 0]);
});

const encode = (sequences) => {
  const bases = { A: 0, C: 1, G: 2, T: 3 };
  const length = sequences[0].length;
  const encoded = new Uint8Array(sequences.length * length);
  sequences.forEach((sequence, index) => {
    encoded.set(Uint8Array.from(sequence, (base) => bases[base] ?? 4), index * length);
  });
  return { encoded, length };
};

function syntheticMosaic() {
  const length = 600;
  const major = Array(length).fill("A");
  const minor = Array(length).fill("A");
  for (let site = 0; site < length; site += 5) major[site] = "C";
  for (let site = 2; site < length; site += 5) minor[site] = "G";
  const recombinant = [...major];
  for (let site = 200; site < 400; site += 1) recombinant[site] = minor[site];
  const majorRelative = major.map((base, site) => site % 17 === 0 ? "T" : base);
  const minorRelative = minor.map((base, site) => site % 19 === 0 ? "C" : base);
  return [recombinant, major, minor, majorRelative, minorRelative].map((sequence) => sequence.join(""));
}

test("MakePhPrScore translation preserves the source directions", () => {
  const candidates = [0, 1, 2];
  const cohort = [0, 1, 2, 3, 4];
  const profiles = {
    first: new Map([
      [0, Float64Array.from([0, 1, 2, 3, 4])],
      [1, Float64Array.from([1, 0, 1, 2, 3])],
      [2, Float64Array.from([2, 1, 0, 1, 2])],
    ]),
    second: new Map([
      [0, Float64Array.from([0, 4, 3, 2, 1])],
      [1, Float64Array.from([1, 0, 1.1, 2.1, 3.1])],
      [2, Float64Array.from([2, 1.1, 0, 1.1, 2.1])],
    ]),
  };
  const scores = sourcePhPrScores(candidates, cohort, profiles);
  assert.equal(scores.phPr.indexOf(Math.min(...scores.phPr)), 0, "the recombinant has the lowest profile correlation");
  assert.equal(scores.subDist.indexOf(Math.max(...scores.subDist)), 0, "removing the recombinant maximizes the remaining average correlation");
  assert.equal(scores.subPhPr.indexOf(Math.max(...scores.subPhPr)), 0, "the recombinant has the largest distance-profile change");
});

test("MakeTrpScore translation counts reordered tree positions", () => {
  const candidates = [0, 1, 2];
  const cohort = [0, 1, 2, 3, 4];
  const profiles = {
    first: new Map([
      [0, Float64Array.from([0, 1, 2, 3, 4])],
      [1, Float64Array.from([1, 0, 1, 2, 3])],
      [2, Float64Array.from([2, 1, 0, 1, 2])],
    ]),
    second: new Map([
      [0, Float64Array.from([0, 4, 3, 2, 1])],
      [1, Float64Array.from([1, 0, 1, 2, 3])],
      [2, Float64Array.from([2, 1, 0, 1, 2])],
    ]),
  };
  const scores = sourceTrpScores(candidates, cohort, profiles);
  assert.ok(scores[0] > scores[1] && scores[0] > scores[2]);
});

test("MakeOUCheck translation follows the source NO/PI/NI topology mapping", () => {
  const candidates = [0, 1, 2];
  const cohort = [0, 1, 2, 3, 4];
  const profiles = {
    first: new Map([
      [0, Float64Array.from([0, 1, 4, 2, 3])],
      [1, Float64Array.from([1, 0, 4, 2, 3])],
      [2, Float64Array.from([4, 4, 0, 3, 2])],
    ]),
    second: new Map([
      [0, Float64Array.from([0, 4, 1, 1, 2])],
      [1, Float64Array.from([4, 0, 4, 3, 3])],
      [2, Float64Array.from([1, 4, 0, 2, 1])],
    ]),
  };
  const scores = sourceOuCheckScores(candidates, cohort, profiles);
  assert.ok(scores[1] > scores[0] && scores[1] > scores[2], "the NO member should gain support when its relationship is the disturbed one");
});

test("SimpleDist exposes the source O:E and O:EDist formulas", () => {
  const candidates = [0, 1, 2];
  const cohort = [0, 1, 2, 3];
  const matrix = (rows) => Float64Array.from(rows.flat());
  const background = matrix([
    [0, 1, 4, 1], [1, 0, 4, 1], [4, 4, 0, 3], [1, 1, 3, 0],
  ]);
  const tract = matrix([
    [0, 4, 1, 2], [4, 0, 4, 1], [1, 4, 0, 1], [2, 1, 1, 0],
  ]);
  const result = sourceSimpleDist(candidates, cohort, background, tract);
  assert.equal(result.inList.length, 3);
  assert.equal(result.simScore.length, 3);
  assert.equal(result.simScoreB.length, 3);
  assert.ok(result.simScore.some((value) => value === 1), "one observed/expected topology should satisfy the strict source rule");
});

test("CalcMaxD scalar reference gives the known mosaic the largest quartet displacement", () => {
  const sequences = syntheticMosaic();
  const { encoded, length } = encode(sequences);
  const result = sourceDmaxScores(encoded, length, [0, 1, 2], [0, 1, 2, 3, 4], { start: 200, end: 400, wraps: false });
  assert.equal(result.values.indexOf(Math.max(...result.values)), 0);
  assert.ok(result.quartetCounts.every((count) => count > 0));
});

test("MakeRCompat parsimony port uses the matching co-recombinant set for every polarity", () => {
  const candidates = [0, 1, 2];
  const cohort = [0, 1, 2, 3, 4];
  const distances = Float64Array.from([
    0, 2, 4, 1, 5,
    2, 0, 4, 3, 4,
    4, 4, 0, 5, 1,
    1, 3, 5, 0, 6,
    5, 4, 1, 6, 0,
  ]);
  const tree = { collapsed: distances, uncollapsed: distances };
  const result = sourceParsimonyScores(candidates, cohort, [tree, tree], [
    { presumedRecombinant: 0, sequenceMembers: [0, 3] },
    { presumedRecombinant: 1, sequenceMembers: [1, 4] },
    { presumedRecombinant: 2, sequenceMembers: [2] },
  ]);
  assert.deepEqual(result.outer, result.inner);
  assert.ok(result.outer.every((value) => Number.isInteger(value) && value >= 0));
  assert.equal(result.outer.length, 3);
});

test("FindSets closure propagates overlapping historical recombinant hyperedges", () => {
  const members = sourceHistoricalSetMembers(
    [0, 1, 2],
    [0, 1, 2, 3, 4, 5],
    [
      { id: "older-a", recombinant: 0, majorParent: 3, minorParent: 4, start: 10, end: 90 },
      { id: "older-b", recombinant: 1, majorParent: 3, minorParent: 5, start: 20, end: 80 },
    ],
    { id: "focal", recombinant: 0, majorParent: 1, minorParent: 2, start: 30, end: 70 },
    100,
  );
  assert.ok(members[0].includes(3));
  assert.ok(members[1].includes(3));
  assert.ok(members[2].includes(3), "membership in exactly two sets must add the taxon to the third");
});

test("MakeSSDistB and OUIndex preserve tree-category weighting", () => {
  const candidates = [0, 1, 2];
  const cohort = [0, 1, 2, 3];
  const first = Float64Array.from([
    0, 1, 3, 1, 1, 0, 2, 1, 3, 2, 0, 2, 1, 1, 2, 0,
  ]);
  const second = Float64Array.from([
    0, 3, 1, 3, 3, 0, 2, 1, 1, 2, 0, 1, 3, 1, 1, 0,
  ]);
  const tree = Float64Array.from([
    0, 2, 4, 2, 2, 0, 4, 2, 4, 4, 0, 3, 2, 2, 3, 0,
  ]);
  const result = sourceSsDistScores(candidates, cohort, first, second, tree);
  assert.equal(result.values.length, 3);
  assert.ok(result.values.every(Number.isFinite));
  assert.equal(sourceOuIndexScores([5, 2, 3], [0, 1, 2])[0], 1);
  assert.deepEqual(sourceOuIndexScores([1, 3, 2], [0, 1, 2]), [0, 1, 1]);
});

test("GetBadDists counts distinct conflicting background-tree categories", () => {
  const { encoded, length } = encode([
    "A".repeat(100),
    "A".repeat(100),
    "T".repeat(100),
    "C".repeat(100),
  ]);
  const tree = Float64Array.from([
    0, 2, 4, 3, 2, 0, 4, 3, 4, 4, 0, 2, 3, 3, 2, 0,
  ]);
  const event = {
    recombinant: 0,
    majorParent: 1,
    minorParent: 2,
    start: 30,
    end: 70,
    coRecombinantSets: [{
      presumedRecombinant: 0,
      evidence: [{
        sequence: 3,
        regionEvidence: [
          { correlationR: 0.9, correlationSdmFiltered: false },
          { correlationR: 0, correlationSdmFiltered: false },
        ],
      }],
    }],
  };
  assert.deepEqual(sourceConflictScores(event, encoded, length, [0, 1, 2], [0, 1, 2, 3], tree), [1, 0, 0]);
});

test("MakeEList and MakeListCorr expose SetDistT and SetDistP", () => {
  const candidates = [0, 1, 2];
  const cohort = [0, 1, 2, 3];
  const background = Float64Array.from([
    0, 1, 4, 1, 1, 0, 4, 1, 4, 4, 0, 3, 1, 1, 3, 0,
  ]);
  const tract = Float64Array.from([
    0, 4, 1, 2, 4, 0, 4, 1, 1, 4, 0, 1, 2, 1, 1, 0,
  ]);
  const inList = sourceSimpleDist(candidates, cohort, background, tract).inList;
  const event = {
    coRecombinantSets: [{
      presumedRecombinant: 0,
      evidence: [{
        sequence: 3,
        regionEvidence: Array.from({ length: 3 }, () => ({
          correlationInversion: 1,
          correlationPermutations: [0.7, 0.9, 0.6, 0.65, 0.55, 0.58],
        })),
      }],
    }],
  };
  const scores = sourceSetDistanceScores(event, candidates, cohort, background, tract, inList);
  assert.equal(scores.setDistT.length, 3);
  assert.equal(scores.setDistP.length, 3);
  assert.ok(scores.setDistT.every((value) => Number.isInteger(value) && value >= 0));
  assert.ok(scores.expectedCoverage.every((value) => value >= 0));
});

test("source profile consensus identifies a known mosaic independent of input polarity", () => {
  const sequences = syntheticMosaic();
  const { encoded, length } = encode(sequences);
  const correct = identifyRecombinantRoles({ recombinant: 0, majorParent: 1, minorParent: 2, start: 200, end: 400, wraps: false }, encoded, sequences.length, length);
  assert.equal(correct.recommended, 0);
  assert.ok(correct.confidence >= 0.6);
  assert.equal(correct.inference, "rdp5-source-distance-consensus");
  assert.deepEqual(correct.implementedComponents.slice(0, 3), ["PhPr", "TreePhPr", "SubPhPr"]);

  const wrongInitialPolarity = identifyRecombinantRoles({ recombinant: 1, majorParent: 0, minorParent: 2, start: 200, end: 400, wraps: false }, encoded, sequences.length, length);
  assert.equal(wrongInitialPolarity.recommended, 0);
  assert.equal(wrongInitialPolarity.orientations.find((orientation) => orientation.recombinant === 0)?.minorParent, 2);
});

test("tree-path evidence activates the expanded source component ledger", () => {
  const sequences = syntheticMosaic();
  const { encoded, length } = encode(sequences);
  const taxa = [0, 1, 2, 3, 4];
  const matrix = (rows) => Float64Array.from(rows.flat());
  const background = matrix([
    [0, 1, 4, 1, 4], [1, 0, 4, 1, 4], [4, 4, 0, 4, 1], [1, 1, 4, 0, 4], [4, 4, 1, 4, 0],
  ]);
  const tract = matrix([
    [0, 4, 1, 4, 1], [4, 0, 4, 1, 4], [1, 4, 0, 4, 1], [4, 1, 4, 0, 4], [1, 4, 1, 4, 0],
  ]);
  const tree = (distances) => ({ taxa, index: new Map(taxa.map((taxon, index) => [taxon, index])), baseDistances: distances, collapsed: distances, uncollapsed: distances, replicates: 100, cutoff: 0.5 });
  const treeBundle = { taxa, pairs: [[null, null], [null, null], [tree(background), tree(tract)]] };
  const result = identifyRecombinantRoles({ recombinant: 0, majorParent: 1, minorParent: 2, start: 200, end: 400, wraps: false }, encoded, sequences.length, length, treeBundle);
  assert.equal(result.treeEvidence, true);
  assert.equal(result.tests.length, 25);
  assert.ok(result.tests.some((entry) => entry.id === "tree-subphpr-conflict"));
  assert.ok(result.tests.some((entry) => entry.sourceRoutine === "MakeTrpGroups + MakeTrpScore"));
  assert.ok(result.tests.some((entry) => entry.sourceRoutine === "SimpleDist SimScore"));
  assert.ok(result.tests.some((entry) => entry.label === "ParsimonyO"));
  assert.ok(result.tests.some((entry) => entry.sourceRoutine === "MakeSSDistB"));
  assert.ok(result.tests.some((entry) => entry.sourceRoutine === "GetBadDists"));
  assert.ok(result.tests.some((entry) => entry.label === "SetDistT"));
  assert.ok(result.tests.some((entry) => entry.label === "SetDistP"));
  assert.ok(result.tests.some((entry) => entry.sourceRoutine === "FinalTrim + RCompatC/RCompatD + MakeConsensusC"));
  assert.ok(result.tests.some((entry) => entry.sourceRoutine === "CalcMaxD + CMaxD2P3 + MakeConsensusC"));
  assert.ok(result.tests.some((entry) => entry.sourceRoutine === "MakeConsensusC MaxS joint rule"));
});

test("an uninformative triplet is flagged role-ambiguous", () => {
  const sequences = ["A".repeat(300), "A".repeat(300), "A".repeat(300)];
  const { encoded, length } = encode(sequences);
  const result = identifyRecombinantRoles({ recombinant: 0, majorParent: 1, minorParent: 2, start: 100, end: 200, wraps: false }, encoded, sequences.length, length);
  assert.equal(result.ambiguous, true);
  assert.equal(result.confidence, 1 / 3);
  assert.equal(result.sourceTieBreak, "GetWinPPfromDists");
  assert.equal(result.recommended, 0);
});
