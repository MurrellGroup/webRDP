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
    sourceCompatibility: "RDP5 BenHMM + DoHMMCyclesSerial + MatchBPtoCI + PolishBP",
    sourceRoutines: ["BenHMM", "DoHMMCyclesSerial", "MatchBPtoCI", "PolishBP"],
    sequenceOrder: [0, 1, 2],
    stateDominantCategories: [0, 2, 1],
    viterbiLogLikelihood: -83.5,
    winningRestart: 3,
    circularPadding: { offset: 142, fittedSites: 569, croppedSites: 285 },
    candidateBreakpoints: [2_180, 120] as [number, number],
    polishedBreakpoints: [2_200, 100] as [number, number],
    polishDecision: { startAdopted: true, endAdopted: true, sameSwitchResolved: false, startMissingBoundary: false, endMissingBoundary: true, revertedForInformation: false, insideVariableSites: 42, outsideVariableSites: 88, startWithin99: true, endWithin99: false, startVariableSiteDistance: 2, endVariableSiteDistance: 4 },
    posteriorThreshold: 0.995,
    confidence99Start: [2_190, 2_210] as [number, number],
    confidence99End: [90, 110] as [number, number],
    switches: [{ position: 2_200, informativeIndex: 141, fromState: 0, toState: 1, confidence95: [2_195, 2_205] as [number, number], confidence99: [2_190, 2_210] as [number, number], sourceCoordinates: [2_190, 2_210, 2_200, 2_195, 2_205], matchedStart: true, matchedEnd: false }],
    posteriorTrace: [{ position: 2_200, informativeIndex: 141, state: 1, probabilities: [0.01, 0.98, 0.01] }],
  },
  methodSignals: [
    { method: "RDP" as const, start: 2_200, end: 100, wraps: true, statistic: 91.2, locator: "RDP5 source" },
    { method: "SiScan" as const, start: 2_210, end: 90, wraps: true, statistic: 12.5, locator: "RDP5 Sister-Scanning pattern 3 topology run", sourceRoutine: "GetSSOL + DoPerms3P + ShrinkRegionC", sourceSiScan: { rawP: 1e-12, rawStart: 2_210, rawEnd: 2_300, runWindows: 7, outgroupSourcePath: "analyst-selected fourth sequence", positionMode: "triplet-variable" as const, gapMode: "strip" as const, window: 200, step: 20, topologyTriplet: [0, 1, 2] as [number, number, number], recombinant: 0, majorParent: 1, minorParent: 2 }, outgroup: 10, outgroupMode: "manual" as const, permutations: 1000, scanPermutations: 100, pattern: 3, scoreFamily: "pattern" as const, baselineTopology: 0, inferredTopology: 1, profile: [{ position: 10, z: 3.2, topology: 1, baselineTopology: 0, pattern: 3, scoreFamily: "pattern" as const }] },
    { method: "MaxChi" as const, start: 2_200, end: 100, wraps: true, statistic: 18.2, locator: "RDP5 pair-equality track 2 · paired source peak basin", sourceRoutine: "FindSubSeqMCPB2 → GrowMChiWinP2", sourceChi: { track: 1, targetSlot: null, informativeSites: 284, halfWindow: 60, boundaryStatistics: [18.2, 22.1] as [number, number], boundaryRanks: [140, 230] as [number, number], growthWidths: [74, 81] as [number, number], direction: -1 as const } },
    { method: "GENECONV" as const, start: 2_205, end: 95, wraps: true, statistic: 86, locator: "RDP5 six-track fragment queue · track 4", sourceRoutine: "FindSubSeqGCAP6/7 → GetFragsP → CalcKMaxP", sourceGeneconv: { track: 3, targetSlot: 0, minorSlot: 1, majorSlot: 2, fragmentScore: 86, informativeSites: 284, matchingSites: 170, mismatchSites: 114, mismatchPenalty: 3, rawP: 1.2e-9, startRank: 41, endRank: 198 } },
    { method: "BootScan" as const, start: 2_190, end: 110, wraps: true, statistic: 0.92, locator: "RDP5 RecScan distance topology 2", sourceRoutine: "BSXoverR2 → SEQBOOT2 → FastBootDistIP → GetPltVal2 → MakeScoresBS/ProbCalc", sourceBootscan: { topology: 1, baselineTopology: 0, bootstrapSupport: 0.92, bootstrapReplicates: 100, runWindows: 17, tractPairMatches: 94, backgroundPairMatches: 12, tractInformativeSites: 108, informativeSites: 284, rawP: 3.4e-12, window: 200, step: 20, relationshipMode: "distance" as const } },
    { method: "3Seq" as const, start: 2_200, end: 100, wraps: true, statistic: 71, locator: "RDP5 3Seq target walk · cycle 2 · descent", sourceRoutine: "FindSubSeqTS → CheckwrapC → GetTSPVal → CheckSplit3Seq/SubPVal", sourceThreeSeq: { target: 0, majorParent: 1, minorParent: 2, direction: 1 as const, upSteps: 180, downSteps: 104, descent: 71, informativeSites: 284, cycle: 1, rawStart: 12, rawEnd: 312, rawP: 8.2e-14, probabilityMode: "exact-table" as const, sourceWrap: true, linearComplement: false, splitRefined: true, fullDescent: 80, splitInformativeSites: 240 } },
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
  options: { ...DEFAULT_OPTIONS, candidateParents: 12, chiSignalsPerTriplet: 96, geneconvSignalsPerTriplet: 112, bootscanWindow: 240, bootscanStep: 12, bootscanCutoff: 0.82, bootscanSignals: 8192, threeSeqExactOperations: 2_000_000, siskanOutgroupMode: "manual", siskanOutgroupSequence: 10, siskanPositionMode: "quartet-variable", siskanGapMode: "fifth-state", siskanScanPermutations: 200, siskanPValuePermutations: 2000 },
  events: [event],
  metrics: { elapsedMs: 12.5, comparisons: 56, engine: "test", tripletMode: "all-concrete-triplets", concreteTripletInputs: true, rdpSignalTruncations: 3, geneconvSignalTruncations: 4, chiSignalTruncations: 2, bootscanSignalTruncations: 5, bootscanBatch: { calls: 2, triplets: 56, usedPairs: 28, windows: 202, replicates: 100, workspaceBytes: 65_536, relationshipMode: "distance" as const }, tripletKernelCalls: { rdp: 56, geneconv: 56, sourceChi: 56, threeSeq: 56, siscan: 56 }, disassembly: { appliedEvents: 1, components: 2, erasedCanonicalBases: 1_000 } },
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
assert.equal(restored.options.chiSignalsPerTriplet, 96);
assert.equal(restored.options.geneconvSignalsPerTriplet, 112);
assert.equal(restored.options.bootscanWindow, 240);
assert.equal(restored.options.bootscanStep, 12);
assert.equal(restored.options.bootscanCutoff, 0.82);
assert.equal(restored.options.bootscanSignals, 8192);
assert.equal(restored.options.threeSeqExactOperations, 2_000_000);
assert.equal(restored.events[0].decision, "accepted");
assert.equal(restored.events[0].note, "reviewed circular positive control");
assert.equal(restored.events[0].wraps, true);
assert.equal(restored.events[0].breakpointModel?.method, "burt-hmm");
assert.equal(restored.events[0].breakpointModel?.sourceCompatibility, "RDP5 BenHMM + DoHMMCyclesSerial + MatchBPtoCI + PolishBP");
assert.equal(restored.events[0].breakpointModel?.viterbiLogLikelihood, -83.5);
assert.equal(restored.events[0].breakpointModel?.winningRestart, 3);
assert.deepEqual(restored.events[0].breakpointModel?.sourceRoutines, ["BenHMM", "DoHMMCyclesSerial", "MatchBPtoCI", "PolishBP"]);
assert.deepEqual(restored.events[0].breakpointModel?.circularPadding, { offset: 142, fittedSites: 569, croppedSites: 285 });
assert.deepEqual(restored.events[0].breakpointModel?.candidateBreakpoints, [2_180, 120]);
assert.equal(restored.events[0].breakpointModel?.polishDecision?.endMissingBoundary, true);
assert.deepEqual(restored.events[0].breakpointModel?.confidence99Start, [2_190, 2_210]);
assert.deepEqual(restored.events[0].breakpointModel?.switches?.[0].confidence99, [2_190, 2_210]);
assert.equal(restored.events[0].breakpointModel?.switches?.[0].informativeIndex, 141);
assert.equal(restored.events[0].breakpointModel?.switches?.[0].matchedStart, true);
assert.equal(restored.events[0].breakpointModel?.posteriorTrace?.[0].informativeIndex, 141);
assert.equal(restored.events[0].methodSignals?.[0].locator, "RDP5 source");
assert.equal(restored.events[0].methodSignals?.[1].outgroup, 10);
assert.equal(restored.events[0].methodSignals?.[1].profile?.[0].z, 3.2);
assert.deepEqual(restored.events[0].methodSignals?.[2].sourceChi?.growthWidths, [74, 81]);
assert.equal(restored.events[0].methodSignals?.[3].sourceGeneconv?.rawP, 1.2e-9);
assert.deepEqual(restored.events[0].methodSignals?.[3].sourceGeneconv && [
  restored.events[0].methodSignals?.[3].sourceGeneconv?.track,
  restored.events[0].methodSignals?.[3].sourceGeneconv?.targetSlot,
  restored.events[0].methodSignals?.[3].sourceGeneconv?.minorSlot,
  restored.events[0].methodSignals?.[3].sourceGeneconv?.majorSlot,
], [3, 0, 1, 2]);
assert.equal(restored.events[0].methodSignals?.[1].sourceSiScan?.rawP, 1e-12);
assert.equal(restored.events[0].methodSignals?.[1].sourceSiScan?.outgroupSourcePath, "analyst-selected fourth sequence");
assert.deepEqual(restored.events[0].methodSignals?.[1].sourceSiScan?.topologyTriplet, [0, 1, 2]);
assert.equal(restored.events[0].methodSignals?.[4].sourceBootscan?.bootstrapSupport, 0.92);
assert.equal(restored.events[0].methodSignals?.[4].sourceBootscan?.baselineTopology, 0);
assert.equal(restored.events[0].methodSignals?.[4].sourceBootscan?.relationshipMode, "distance");
assert.equal(restored.events[0].methodSignals?.[5].sourceThreeSeq?.rawP, 8.2e-14);
assert.equal(restored.events[0].methodSignals?.[5].sourceThreeSeq?.probabilityMode, "exact-table");
assert.equal(restored.events[0].methodSignals?.[5].sourceThreeSeq?.sourceWrap, true);
assert.equal(restored.events[0].methodSignals?.[5].sourceThreeSeq?.linearComplement, false);
assert.equal(restored.events[0].methodSignals?.[5].sourceThreeSeq?.splitRefined, true);
assert.equal(restored.events[0].methodSignals?.[5].sourceThreeSeq?.fullDescent, 80);
assert.equal(restored.events[0].methodSignals?.[5].sourceThreeSeq?.splitInformativeSites, 240);
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
assert.equal(restored.metrics?.geneconvSignalTruncations, 4);
assert.equal(restored.metrics?.chiSignalTruncations, 2);
assert.equal(restored.metrics?.bootscanSignalTruncations, 5);
assert.deepEqual(restored.metrics?.bootscanBatch, { calls: 2, triplets: 56, usedPairs: 28, windows: 202, replicates: 100, workspaceBytes: 65_536, relationshipMode: "distance" });
assert.deepEqual(restored.metrics?.tripletKernelCalls, { rdp: 56, geneconv: 56, sourceChi: 56, threeSeq: 56, siscan: 56 });
assert.equal(restored.metrics?.tripletMode, "all-concrete-triplets");
assert.equal(restored.metrics?.concreteTripletInputs, true);

const masked = exportRecombinationFree(restored.alignment, restored.events, "mask")[0];
const maskedAlignment = parseAlignment(masked.content, masked.filename);
assert.equal(maskedAlignment.sequences[0].sequence.match(/N/g)?.length, 300);
assert.equal(maskedAlignment.sequences[1].sequence.match(/N/g)?.length ?? 0, 0);
