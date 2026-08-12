import assert from "node:assert/strict";
import test from "node:test";
import { inferAncestralEventClusters, sourceCalCr, sourceWeightedEventClusters } from "../public/rdp-clustering.js";
import {
  bootstrapTreeDistance,
  buildBootstrapRegionTree,
  buildEventBootstrapTreeCohorts,
} from "../public/rdp-bootstrap-tree.js";

test("seeded JC/NJ bootstrap retains a strongly supported sister split", () => {
  const length = 200;
  const encoded = new Uint8Array(4 * length);
  for (let site = 0; site < length; site += 1) {
    encoded[site] = site < 100 ? 0 : 2;
    encoded[length + site] = site < 100 ? 0 : 2;
    encoded[2 * length + site] = site < 100 ? 1 : 3;
    encoded[3 * length + site] = site < 100 ? 1 : 3;
  }
  const region = { kind: "segments", segments: [[0, length]], id: "all" };
  const first = buildBootstrapRegionTree(encoded, length, [0, 1, 2, 3], region, { replicates: 100, cutoff: 0.5, bootstrapBlocks: 256 }, 91);
  const second = buildBootstrapRegionTree(encoded, length, [0, 1, 2, 3], region, { replicates: 100, cutoff: 0.5, bootstrapBlocks: 256 }, 91);
  assert.equal(first.exactSiteBootstrap, true);
  assert.deepEqual(first.splits, second.splits);
  assert.ok(first.splits.some((split) => split.support >= 0.99));
  assert.ok(bootstrapTreeDistance(first, 0, 1) < bootstrapTreeDistance(first, 0, 2));
});

test("CalCR distinguishes direct ancestry pattern from a category inversion", () => {
  const presumed = [0.8, 0.1, 0.1, 0.1, 0.8, 0.1];
  const direct = sourceCalCr(presumed, [...presumed]);
  assert.equal(direct.inversion, 0);
  assert.ok(direct.r > 0.99);
  const inverse = sourceCalCr(presumed, [0.1, 0.8, 0.1, 0.8, 0.1, 0.1]);
  assert.equal(inverse.inversion, 1);
  assert.ok(inverse.selectedR > 0.99);
  assert.ok(inverse.r < 0.83);
});

test("GetSupers-compatible weighted merging does not single-link a weak chain", () => {
  const distances = [
    [0, 0.05, 0.30],
    [0.05, 0, 0.05],
    [0.30, 0.05, 0],
  ];
  assert.deepEqual(sourceWeightedEventClusters(distances, 0.1), [[0, 1], [2]]);
});

test("bounded bootstrap cohorts retain the detecting triplet while covering every candidate", () => {
  const length = 16;
  const sequenceCount = 11;
  const encoded = new Uint8Array(sequenceCount * length);
  const left = { id: "left", kind: "segments", segments: [[0, 8]] };
  const right = { id: "right", kind: "segments", segments: [[8, 16]] };
  const result = buildEventBootstrapTreeCohorts(
    encoded,
    length,
    sequenceCount,
    [0, 1, 2],
    [[left, right]],
    { clusterTreeTaxaLimit: 6, clusterBootstrapReplicates: 0, clusterBootstrapBlocks: 16 },
    7,
  );
  assert.equal(result.candidateComplete, true);
  assert.equal(result.cohorts.length, 3);
  assert.equal(result.byTaxon.size, sequenceCount);
  for (let candidate = 3; candidate < sequenceCount; candidate += 1) {
    const bundle = result.byTaxon.get(candidate);
    assert.ok(bundle.taxa.includes(candidate));
    assert.ok([0, 1, 2].every((member) => bundle.taxa.includes(member)));
    assert.ok(bundle.taxa.length <= 6);
  }
});

function event(id, recombinant, start, end) {
  return {
    id,
    recombinant,
    majorParent: 3,
    minorParent: 4,
    start,
    end,
    wraps: false,
    evidence: [{ correctedP: 1e-8 }],
    chiSquare: 40,
    warnings: [],
  };
}

test("signals in distinct descendant sequences cluster as one ancestral recombinant", () => {
  const length = 1_000;
  const encoded = new Uint8Array(6 * length);
  for (let site = 0; site < length; site += 1) {
    encoded[3 * length + site] = 0;
    encoded[4 * length + site] = 1;
    encoded[0 * length + site] = site >= 300 && site < 600 ? 1 : 0;
    encoded[1 * length + site] = site >= 300 && site < 600 ? 1 : 0;
    encoded[2 * length + site] = site >= 720 && site < 860 ? 1 : 0;
    encoded[5 * length + site] = site >= 300 && site < 600 ? 1 : 0;
    if (site % 137 === 0) encoded[1 * length + site] = 2;
  }
  const result = inferAncestralEventClusters([
    event("e1", 0, 300, 600),
    event("e2", 1, 300, 600),
    event("e3", 2, 720, 860),
  ], encoded, 6, length, {
    ancestralClustering: true,
    clusterFlankVnps: 60,
    clusterMinimumSets: 2,
    clusterCorrelationAlpha: 0.05,
    clusterCorrelationR: 0.83,
    clusterSignalOverlap: 0.3,
    clusterTopologyMargin: 0,
    clusterMinimumConfidence: 0,
    clusterSourceSimilarity: 0.1,
  });
  assert.equal(result.events[0].groupId, "ancestry-001");
  assert.equal(result.events[1].groupId, "ancestry-001");
  assert.equal(result.events[2].groupId, null);
  assert.deepEqual(result.events[0].ancestralCluster.sequenceMembers, [0, 1, 5]);
  assert.equal(result.events[0].coRecombinantSets.length, 3);
  const actualOrientation = result.events[0].coRecombinantSets.find((set) => set.presumedRecombinant === 0);
  assert.deepEqual(actualOrientation.sequenceMembers, [0, 1, 5]);
  const unsignalledDescendant = actualOrientation.evidence.find((entry) => entry.sequence === 5);
  assert.ok(unsignalledDescendant && unsignalledDescendant.sets >= 2);
  assert.equal(unsignalledDescendant.treeBootstrap.included, true);
  assert.equal(unsignalledDescendant.treeBootstrap.candidateComplete, true);
  assert.equal(unsignalledDescendant.treeBootstrap.cohortCount, 1);
  assert.ok(result.events[0].ancestralCluster.sourceMerge.pairDistances.some((pair) => pair.belowThreshold));
});

test("two signals in the same daughter are never called co-recombinant descendants", () => {
  const length = 500;
  const encoded = new Uint8Array(3 * length);
  const result = inferAncestralEventClusters([
    { ...event("e1", 0, 100, 250), majorParent: 1, minorParent: 2 },
    { ...event("e2", 0, 120, 260), majorParent: 1, minorParent: 2 },
  ], encoded, 3, length, { ancestralClustering: true, clusterSourceSimilarity: 0.1 });
  assert.equal(result.clusters.length, 0);
  assert.ok(result.events.every((candidate) => candidate.groupId === null));
});
