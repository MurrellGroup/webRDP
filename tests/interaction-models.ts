import assert from "node:assert/strict";
import { affinityDescription, classifyParentAffinity, parentInformativeSites } from "../app/alignment-highlighter";
import { buildReconstructionModel, eventOverlapBases } from "../app/reconstruction";
import type { RdpEvent } from "../app/rdp-core";
import { layoutNeighborJoiningTree } from "../app/tree-layout";

assert.deepEqual(classifyParentAffinity("A", ["A", "C"]), { kind: "unique", parentSlots: [0] });
assert.deepEqual(classifyParentAffinity("A", ["A", "A", "G"]), { kind: "shared", parentSlots: [0, 1] });
assert.deepEqual(classifyParentAffinity("T", ["A", "C"]), { kind: "novel", parentSlots: [] });
assert.deepEqual(classifyParentAffinity("N", ["A", "C"]), { kind: "missing", parentSlots: [] });
assert.equal(affinityDescription({ kind: "shared", parentSlots: [0, 2] }, ["P1", "P2", "P3"]), "shared by P1, P3");
assert.deepEqual(parentInformativeSites(["ACGT", "AGGT", "ACGA"], [0, 1], 4), [1]);

const treeLayout = layoutNeighborJoiningTree({
  length: 0,
  children: [
    { name: "A", length: 0.1 },
    { length: 0.05, children: [{ name: "B", length: 0.02 }, { name: "C", length: 0 }] },
  ],
});
assert.equal(treeLayout.labels.length, 3);
assert.equal(treeLayout.edges.length, 4);
assert.equal(treeLayout.joints.length, 2);
assert.equal(treeLayout.zeroLengthBranches, 1);
assert.ok(treeLayout.path.includes("V") && treeLayout.path.includes("H"));
assert.doesNotMatch(treeLayout.path, /NaN|Infinity/);
treeLayout.edges.forEach((edge) => {
  assert.ok(Number.isFinite(edge.parentX) && Number.isFinite(edge.childX) && Number.isFinite(edge.childY));
  assert.ok(edge.childX >= edge.parentX, "non-negative tree branches must not draw backwards");
});

function makeEvent(index: number, patch: Partial<RdpEvent> = {}): RdpEvent {
  return {
    id: `event-${index}`,
    recombinant: index,
    majorParent: (index + 1) % 4,
    minorParent: (index + 2) % 4,
    start: 100,
    end: 500,
    wraps: false,
    confidenceStart: [95, 105],
    confidenceEnd: [495, 505],
    evidence: [],
    chiSquare: 0,
    informativeSites: 0,
    decision: "accepted",
    warnings: [],
    note: "",
    source: "manual",
    groupId: null,
    history: [],
    evidenceStale: false,
    diagnostics: { tractVariableDensity: 0, backgroundVariableDensity: 0, rateRatio: 0, parentConflictRate: 0, parentDiscriminatingSites: 0, diffuseIncompatibility: false },
    ...patch,
  };
}

const reconstructionEvents = [
  makeEvent(0, { recombinant: 0, groupId: "ancestry-a", evidenceStale: true }),
  makeEvent(1, { recombinant: 0, start: 300, end: 420, groupId: "ancestry-a", decision: "unreviewed" }),
  makeEvent(2, { recombinant: 2, majorParent: 0, start: 600, end: 800 }),
];
const reconstruction = buildReconstructionModel(reconstructionEvents, 1_000);
assert.equal(reconstruction.staleFromIndex, 0);
assert.deepEqual(reconstruction.downstreamIndexes, [1, 2]);
assert.equal(reconstruction.nextReviewIndex, 0);
assert.ok(reconstruction.relationships.some((relationship) => relationship.kind === "possible-overprint" && relationship.overlapBases === 120));
assert.ok(reconstruction.relationships.some((relationship) => relationship.kind === "event-group"));
assert.ok(reconstruction.relationships.some((relationship) => relationship.kind === "recombinant-parent"));

const circular = makeEvent(3, { start: 900, end: 120, wraps: true });
const origin = makeEvent(4, { start: 20, end: 80 });
assert.equal(eventOverlapBases(circular, origin, 1_000), 60);

console.log("interaction model checks passed");
