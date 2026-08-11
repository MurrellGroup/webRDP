import assert from "node:assert/strict";
import { affinityDescription, classifyParentAffinity, parentInformativeSites } from "../app/alignment-highlighter";
import { buildReconstructionModel, eventOverlapBases } from "../app/reconstruction";
import { computeBreakpointPairDensity, computeLocalDiscordanceMatrices, computeRegionSeparationMatrices, eventContainsPosition } from "../app/pattern-matrices";
import { AUTO_RESOLVE_PRESETS, applyAutoResolutionPlan, filterResolvedEventDuplicates, planAutoResolution, rescanTargetsForBarrier } from "../app/auto-resolve";
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

const evidenceMethods = ["RDP", "GENECONV", "BootScan", "MaxChi", "Chimaera", "SiScan", "3Seq"] as const;
const evidence = (supporting: number, correctedP: number): RdpEvent["evidence"] => evidenceMethods.map((method, index) => ({
  method,
  pValue: correctedP / 2,
  correctedP,
  score: index < supporting ? 0.9 : 0.1,
  supported: index < supporting,
  statistic: index < supporting ? 12 : 1,
  statisticLabel: "test statistic",
  calibration: "deterministic fixture",
}));
const strongEvent = makeEvent(20, {
  recombinant: 0,
  start: 100,
  end: 500,
  confidenceStart: [98, 102],
  confidenceEnd: [498, 502],
  evidence: evidence(5, 1e-8),
  informativeSites: 140,
  decision: "unreviewed",
  diagnostics: { tractVariableDensity: 0.2, backgroundVariableDensity: 0.18, rateRatio: 1.1, parentConflictRate: 0.02, parentDiscriminatingSites: 140, diffuseIncompatibility: false },
});
const dependentEvent = makeEvent(21, {
  recombinant: 3,
  majorParent: 0,
  minorParent: 2,
  start: 600,
  end: 850,
  confidenceStart: [598, 602],
  confidenceEnd: [848, 852],
  evidence: evidence(5, 1e-7),
  informativeSites: 120,
  decision: "unreviewed",
  diagnostics: { tractVariableDensity: 0.2, backgroundVariableDensity: 0.19, rateRatio: 1.05, parentConflictRate: 0.01, parentDiscriminatingSites: 120, diffuseIncompatibility: false },
});
const weakEvent = makeEvent(22, {
  recombinant: 2,
  evidence: evidence(0, 0.8),
  informativeSites: 3,
  decision: "unreviewed",
  warnings: ["Possible misalignment artefact"],
  diagnostics: { tractVariableDensity: 0.9, backgroundVariableDensity: 0.05, rateRatio: 18, parentConflictRate: 0.8, parentDiscriminatingSites: 3, diffuseIncompatibility: true },
});
const conservativePlan = planAutoResolution([strongEvent, dependentEvent, weakEvent], 1_000, AUTO_RESOLVE_PRESETS.conservative);
assert.equal(conservativePlan.entries[0].recommendation, "accept");
assert.equal(conservativePlan.entries[1].recommendation, "accept");
assert.equal(conservativePlan.entries[2].recommendation, "reject");
assert.equal(conservativePlan.barriers[0].afterEventIndex, 0, "recombinant-parent use must stop the ordered queue after its causal event");
assert.deepEqual(conservativePlan.barriers[0].impactedTargetIndexes, [2, 3]);
const firstAutoBatch = applyAutoResolutionPlan([strongEvent, dependentEvent, weakEvent], conservativePlan, conservativePlan.barriers[0].afterEventIndex, "conservative", "2026-08-11T00:00:00.000Z");
assert.equal(firstAutoBatch.events[0].decision, "accepted");
assert.equal(firstAutoBatch.events[1].decision, "unreviewed", "events beyond a rescan barrier must remain untouched");
assert.match(firstAutoBatch.events[0].history.at(-1)?.summary ?? "", /score/);
assert.deepEqual(rescanTargetsForBarrier(conservativePlan.barriers[0], 10, AUTO_RESOLVE_PRESETS.conservative), { targetIndexes: [2, 3], scope: "targeted" });
const filteredRescan = filterResolvedEventDuplicates([
  makeEvent(30, { recombinant: 0, majorParent: 1, minorParent: 2, start: 102, end: 498 }),
  makeEvent(31, { recombinant: 0, majorParent: 1, minorParent: 2, start: 220, end: 320 }),
  makeEvent(32, { recombinant: 0, majorParent: 1, minorParent: 3, start: 102, end: 498 }),
], [strongEvent], 1_000);
assert.deepEqual(filteredRescan.map((event) => event.id), ["event-31", "event-32"], "rescans suppress only close same-parent duplicates, not nested or alternative-parent hypotheses");
const noRescanPlan = planAutoResolution([strongEvent, dependentEvent], 1_000, { ...AUTO_RESOLVE_PRESETS.conservative, rescanStrategy: "off" });
assert.equal(noRescanPlan.barriers.length, 0);

const moderateEvent = makeEvent(23, {
  evidence: evidence(2, 0.05),
  informativeSites: 20,
  decision: "unreviewed",
  diagnostics: { tractVariableDensity: 0.2, backgroundVariableDensity: 0.18, rateRatio: 1.15, parentConflictRate: 0.04, parentDiscriminatingSites: 20, diffuseIncompatibility: false },
});
assert.notEqual(planAutoResolution([moderateEvent], 1_000, AUTO_RESOLVE_PRESETS.conservative).entries[0].recommendation, "accept");
assert.equal(planAutoResolution([moderateEvent], 1_000, AUTO_RESOLVE_PRESETS.aggressive).entries[0].recommendation, "accept");
assert.equal(planAutoResolution([{ ...strongEvent, evidenceStale: true }], 1_000, AUTO_RESOLVE_PRESETS.aggressive).entries[0].recommendation, "review");
const blockedDependentPlan = planAutoResolution([{ ...strongEvent, evidenceStale: true }, dependentEvent], 1_000, AUTO_RESOLVE_PRESETS.conservative);
assert.equal(blockedDependentPlan.entries[1].recommendation, "review");
assert.match(blockedDependentPlan.entries[1].reasons.at(-1) ?? "", /blocked because unresolved E1/);
assert.equal(planAutoResolution([{ ...strongEvent, decision: "accepted" }], 1_000, AUTO_RESOLVE_PRESETS.conservative).entries[0].recommendation, "keep");

const circular = makeEvent(3, { start: 900, end: 120, wraps: true });
const origin = makeEvent(4, { start: 20, end: 80 });
assert.equal(eventOverlapBases(circular, origin, 1_000), 60);
assert.equal(eventContainsPosition(circular, 950), true);
assert.equal(eventContainsPosition(circular, 60), true);
assert.equal(eventContainsPosition(circular, 500), false);

const pairDensity = computeBreakpointPairDensity([
  makeEvent(5, { start: 100, end: 300 }),
  makeEvent(6, { start: 100, end: 300 }),
], 1_000, 10);
assert.equal(pairDensity.values[1 * 10 + 3], 2);
assert.equal(pairDensity.values[3 * 10 + 1], 2);
assert.equal(pairDensity.maximum, 2);

const separation = computeRegionSeparationMatrices([makeEvent(7, { start: 200, end: 600 })], 1_000, 10);
assert.equal(separation.observed[0 * 10 + 2], 1);
assert.equal(separation.observed[2 * 10 + 4], 0);
assert.ok([...separation.standardizedResidual].every(Number.isFinite));

const discordance = computeLocalDiscordanceMatrices([
  "A".repeat(80),
  `${"A".repeat(40)}${"C".repeat(40)}`,
  `${"C".repeat(40)}${"A".repeat(40)}`,
  "C".repeat(80),
], 80, 4, 4, 20);
assert.equal(discordance.pairCount, 6);
assert.ok(discordance.maximumRmsDeviation > 0);
assert.ok(discordance.maximumCorrelationLoss > 0);
assert.deepEqual([...discordance.rmsDeviation], [...discordance.rmsDeviation].map((value, index, values) => values[(index % 4) * 4 + Math.floor(index / 4)]));

console.log("interaction model checks passed");
