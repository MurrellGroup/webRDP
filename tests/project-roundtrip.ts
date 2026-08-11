import assert from "node:assert/strict";
import {
  DEFAULT_OPTIONS,
  demoEvent,
  eventLength,
  exportRecombinationFree,
  makeDemoAlignment,
  parseGenomeAnnotations,
  parseAlignment,
  parseProject,
  serializeProject,
} from "../app/rdp-core";

const event = {
  ...demoEvent(),
  start: 2_200,
  end: 100,
  wraps: true,
  breakpointModel: {
    method: "burt-hmm" as const,
    informativeSites: 284,
    states: 3,
    stateSwitches: 2,
    majorFit: 0.96,
    minorFit: 0.94,
    sourceParity: true,
    sourceCompatibility: "RDP5 BenHMM + DoHMMCyclesSerial",
    posteriorThreshold: 0.995,
    confidence99Start: [2_190, 2_210] as [number, number],
    confidence99End: [90, 110] as [number, number],
    switches: [{ position: 2_200, fromState: 0, toState: 1, confidence95: [2_195, 2_205] as [number, number], confidence99: [2_190, 2_210] as [number, number] }],
  },
  methodSignals: [
    { method: "RDP" as const, start: 2_200, end: 100, wraps: true, statistic: 91.2, locator: "RDP5 source" },
    { method: "SiScan" as const, start: 2_210, end: 90, wraps: true, statistic: 12.5, locator: "RDP5 Sister-Scanning pattern 3 topology run", sourceRoutine: "GetSSOL + DoPerms3P + ShrinkRegionC", outgroup: 10, outgroupMode: "manual" as const, permutations: 1000, scanPermutations: 100, pattern: 3, scoreFamily: "pattern" as const, baselineTopology: 0, inferredTopology: 1, profile: [{ position: 10, z: 3.2, topology: 1, baselineTopology: 0, pattern: 3, scoreFamily: "pattern" as const }] },
  ],
  coRecombinantSets: [{
    presumedRecombinant: 0,
    parents: [1, 2],
    sequenceMembers: [0, 3],
    testedSequences: 5,
    requiredEvidenceSets: 2,
    evidence: [{
      sequence: 3,
      sets: 2,
      phylogenetic: true,
      distance: true,
      detectableSignal: false,
      bestCorrelation: { r: 0.94, pValue: 0.005, inversion: 0 },
      topologyMargin: 0.12,
      regionEvidence: [{ pair: "5-prime breakpoint", phylogenetic: true, movesTogether: true, sisterTogether: false, topologyMargin: 0.12, treeSourceScore: 8, bootstrapSupport: 0.93, bootstrapReplicates: 100, bootstrapCutoff: 0.5, treeExcluded: false, correlationR: 0.94, correlationP: 0.005, correlationInversion: 0, correlationPermutations: [0.94, 0.3, 0.2, 0.1, 0.4, 0.35], correlationSdmFiltered: false }],
    }],
  }],
  componentProvenance: {
    reconstruction: "rdp5-signal-disassembly" as const,
    appliedEventIds: ["accepted-predecessor"],
    recombinant: { originIndex: 0, kind: "extracted-tract" as const, lineage: ["accepted-predecessor"], sourceEventId: "accepted-predecessor", parentLineage: [], start: 1_800, end: 2_300, wraps: false, erasedEventIds: [] },
    majorParent: { originIndex: 1, kind: "remainder" as const, lineage: [], erasedEventIds: [] },
    minorParent: { originIndex: 2, kind: "remainder" as const, lineage: [], erasedEventIds: [] },
  },
  structuralUncertainty: {
    source: "rdp5-erased-signal-boundary" as const,
    originalStart: 2_000,
    originalEnd: 300,
    originalWraps: true,
    piece: 2,
    pieces: 3,
    uncertainStart: true,
    uncertainEnd: true,
    adjacentEventIds: ["accepted-predecessor"],
  },
  recombinantIdentification: {
    inference: "rdp5-source-profile-consensus" as const,
    candidates: [0, 1, 2],
    recommended: 0,
    recommendedMajorParent: 1,
    recommendedMinorParent: 2,
    confidence: 0.72,
    ambiguous: false,
    sourceThreshold: 0.6,
    orientations: [
      { recombinant: 0, majorParent: 1, minorParent: 2, affinitySwitch: 0.2, candidateIndex: 0, sourcePoints: 54, sourceScore: 100, sourceShare: 0.72 },
      { recombinant: 1, majorParent: 0, minorParent: 2, affinitySwitch: 0.01, candidateIndex: 1, sourcePoints: 14, sourceScore: 26, sourceShare: 0.19 },
      { recombinant: 2, majorParent: 1, minorParent: 0, affinitySwitch: -0.02, candidateIndex: 2, sourcePoints: 7, sourceScore: 13, sourceShare: 0.09 },
    ],
    tests: [
      { id: "phpr", label: "PhPr", sourceRoutine: "MakePhPrScore(FMat, SMat)", direction: "lower" as const, values: [-0.8, 0.7, 0.8], points: [8, 4, 0], fullWeight: 8, partialWeight: 4, winnerIndexes: [0], decisive: true },
      { id: "oucheck", label: "OuCheck", sourceRoutine: "MakeINList + MakeOUCheck", direction: "higher" as const, values: [2, -1, -1], points: [5, -5, -5], fullWeight: 5, partialWeight: 0, winnerIndexes: [0], decisive: true },
    ],
    cohortSize: 12,
    sourceSequenceCount: 12,
    sampled: false,
    treeEvidence: true,
    bootstrapReplicates: 100,
    bootstrapCutoff: 0.5,
    quartetCohortSize: 12,
    quartetCounts: [165, 165, 165],
    dmaxWasmAccelerated: true,
    sourceTieBreak: "last-inclusive-maximum",
    sourceTieBreakValues: [1, 2, 3],
    implementedComponents: ["PhPr", "TreePhPr", "dMax (VisRD)", "ParsimonyO", "Conflict"],
    pendingComponents: ["SetDistT / SetDistP"],
  },
  decision: "accepted" as const,
  note: "reviewed circular positive control",
};
const serialized = serializeProject({
  alignment: {
    ...makeDemoAlignment(),
    features: parseGenomeAnnotations("genome\ttest\tCDS\t101\t900\t.\t+\t0\tID=cds1;Name=rep", "fixture.gff3", 2_400),
  },
  options: { ...DEFAULT_OPTIONS, candidateParents: 12, siskanOutgroupMode: "manual", siskanOutgroupSequence: 10, siskanPositionMode: "quartet-variable", siskanGapMode: "fifth-state", siskanScanPermutations: 200, siskanPValuePermutations: 2000 },
  events: [event],
  metrics: { elapsedMs: 12.5, comparisons: 56, engine: "test", rdpSignalTruncations: 3, disassembly: { appliedEvents: 1, components: 2, erasedCanonicalBases: 1_000 } },
  distance: [0, 0.1, 0.1, 0],
  auditLog: [{ id: "audit-1", timestamp: "2026-08-11T00:00:00.000Z", action: "Accepted event", summary: "Reviewed fixture.", eventId: event.id }],
});
const restored = parseProject(serialized);

assert.equal(restored.schema, "rdp-web/0.5");
assert.equal(restored.alignment.sequences[0].sequence, makeDemoAlignment().sequences[0].sequence);
assert.equal(restored.options.candidateParents, 12);
assert.equal(restored.options.siskanOutgroupMode, "manual");
assert.equal(restored.options.siskanOutgroupSequence, 10);
assert.equal(restored.options.siskanPositionMode, "quartet-variable");
assert.equal(restored.options.siskanGapMode, "fifth-state");
assert.equal(restored.options.siskanPValuePermutations, 2000);
assert.equal(restored.events[0].decision, "accepted");
assert.equal(restored.events[0].note, "reviewed circular positive control");
assert.equal(restored.events[0].wraps, true);
assert.equal(restored.events[0].breakpointModel?.method, "burt-hmm");
assert.equal(restored.events[0].breakpointModel?.sourceCompatibility, "RDP5 BenHMM + DoHMMCyclesSerial");
assert.deepEqual(restored.events[0].breakpointModel?.confidence99Start, [2_190, 2_210]);
assert.deepEqual(restored.events[0].breakpointModel?.switches?.[0].confidence99, [2_190, 2_210]);
assert.equal(restored.events[0].methodSignals?.[0].locator, "RDP5 source");
assert.equal(restored.events[0].methodSignals?.[1].outgroup, 10);
assert.equal(restored.events[0].methodSignals?.[1].profile?.[0].z, 3.2);
assert.deepEqual(restored.events[0].coRecombinantSets?.[0].sequenceMembers, [0, 3]);
assert.equal(restored.events[0].coRecombinantSets?.[0].evidence[0].bestCorrelation?.r, 0.94);
assert.equal(restored.events[0].coRecombinantSets?.[0].evidence[0].regionEvidence?.[0].bootstrapSupport, 0.93);
assert.equal(restored.events[0].coRecombinantSets?.[0].evidence[0].regionEvidence?.[0].correlationPermutations?.[0], 0.94);
assert.deepEqual(restored.events[0].componentProvenance?.recombinant.lineage, ["accepted-predecessor"]);
assert.equal(restored.events[0].componentProvenance?.recombinant.kind, "extracted-tract");
assert.equal(restored.events[0].structuralUncertainty?.piece, 2);
assert.equal(restored.events[0].structuralUncertainty?.originalWraps, true);
assert.deepEqual(restored.events[0].structuralUncertainty?.adjacentEventIds, ["accepted-predecessor"]);
assert.equal(restored.events[0].recombinantIdentification?.recommended, 0);
assert.equal(restored.events[0].recombinantIdentification?.tests[0].sourceRoutine, "MakePhPrScore(FMat, SMat)");
assert.equal(restored.events[0].recombinantIdentification?.tests[1].points[1], -5);
assert.equal(restored.events[0].recombinantIdentification?.orientations[0].sourceShare, 0.72);
assert.equal(restored.events[0].recombinantIdentification?.quartetCohortSize, 12);
assert.deepEqual(restored.events[0].recombinantIdentification?.quartetCounts, [165, 165, 165]);
assert.equal(restored.events[0].recombinantIdentification?.dmaxWasmAccelerated, true);
assert.equal(restored.events[0].recombinantIdentification?.sourceTieBreak, "last-inclusive-maximum");
assert.deepEqual(restored.events[0].recombinantIdentification?.sourceTieBreakValues, [1, 2, 3]);
assert.equal(restored.events[0].history[0].action, event.history[0].action);
assert.equal(restored.events[0].evidenceStale, false);
assert.equal(restored.alignment.features?.[0].name, "rep");
assert.equal(eventLength(restored.events[0], restored.alignment.length), 300);
assert.equal(restored.events[0].evidence[0].calibration, event.evidence[0].calibration);
assert.deepEqual(restored.distance, [0, 0.1, 0.1, 0]);
assert.equal(restored.auditLog[0].eventId, event.id);
assert.equal(restored.metrics?.disassembly?.components, 2);
assert.equal(restored.metrics?.rdpSignalTruncations, 3);

const masked = exportRecombinationFree(restored.alignment, restored.events, "mask")[0];
const maskedAlignment = parseAlignment(masked.content, masked.filename);
assert.equal(maskedAlignment.sequences[0].sequence.match(/N/g)?.length, 300);
assert.equal(maskedAlignment.sequences[1].sequence.match(/N/g)?.length ?? 0, 0);
