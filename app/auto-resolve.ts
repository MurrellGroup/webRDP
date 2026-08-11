import { eventLength, type EventDecision, type RdpEvent } from "./rdp-core";
import { eventOverlapBases } from "./reconstruction";

export type AutoResolvePresetName = "conservative" | "balanced" | "aggressive";
export type AutoResolveRecommendation = "accept" | "reject" | "review" | "keep";
export type AutoResolveRescanStrategy = "off" | "targeted" | "adaptive" | "full";

export interface AutoResolveSettings {
  acceptScore: number;
  rejectScore: number;
  minimumSupportingMethods: number;
  maximumCorrectedP: number;
  minimumInformativeSites: number;
  maximumBreakpointUncertainty: number;
  maximumParentConflict: number;
  maximumRateFold: number;
  minimumRoleConfidence: number;
  blockSevereWarnings: boolean;
  revisitReviewed: boolean;
  consensusWeight: number;
  significanceWeight: number;
  informationWeight: number;
  breakpointWeight: number;
  diagnosticsWeight: number;
  roleWeight: number;
  rescanStrategy: AutoResolveRescanStrategy;
  rescanRiskThreshold: number;
  overlapTriggerFraction: number;
  minimumEventsBetweenRescans: number;
  maximumRescanRounds: number;
  adaptiveFullTargetFraction: number;
  sameRecombinantRisk: number;
  recombinantParentRisk: number;
  groupedEventRisk: number;
  acceptedWithdrawalRisk: number;
}

export const AUTO_RESOLVE_PRESETS: Record<AutoResolvePresetName, AutoResolveSettings> = {
  conservative: {
    acceptScore: 86,
    rejectScore: 18,
    minimumSupportingMethods: 4,
    maximumCorrectedP: 0.01,
    minimumInformativeSites: 40,
    maximumBreakpointUncertainty: 0.15,
    maximumParentConflict: 0.1,
    maximumRateFold: 2,
    minimumRoleConfidence: 0.7,
    blockSevereWarnings: true,
    revisitReviewed: false,
    consensusWeight: 30,
    significanceWeight: 25,
    informationWeight: 15,
    breakpointWeight: 15,
    diagnosticsWeight: 15,
    roleWeight: 15,
    rescanStrategy: "targeted",
    rescanRiskThreshold: 35,
    overlapTriggerFraction: 0.15,
    minimumEventsBetweenRescans: 1,
    maximumRescanRounds: 3,
    adaptiveFullTargetFraction: 0.25,
    sameRecombinantRisk: 52,
    recombinantParentRisk: 58,
    groupedEventRisk: 28,
    acceptedWithdrawalRisk: 75,
  },
  balanced: {
    acceptScore: 74,
    rejectScore: 30,
    minimumSupportingMethods: 3,
    maximumCorrectedP: 0.05,
    minimumInformativeSites: 24,
    maximumBreakpointUncertainty: 0.3,
    maximumParentConflict: 0.2,
    maximumRateFold: 3,
    minimumRoleConfidence: 0.6,
    blockSevereWarnings: true,
    revisitReviewed: false,
    consensusWeight: 28,
    significanceWeight: 22,
    informationWeight: 16,
    breakpointWeight: 14,
    diagnosticsWeight: 20,
    roleWeight: 12,
    rescanStrategy: "adaptive",
    rescanRiskThreshold: 48,
    overlapTriggerFraction: 0.25,
    minimumEventsBetweenRescans: 3,
    maximumRescanRounds: 2,
    adaptiveFullTargetFraction: 0.3,
    sameRecombinantRisk: 46,
    recombinantParentRisk: 54,
    groupedEventRisk: 24,
    acceptedWithdrawalRisk: 70,
  },
  aggressive: {
    acceptScore: 62,
    rejectScore: 42,
    minimumSupportingMethods: 2,
    maximumCorrectedP: 0.1,
    minimumInformativeSites: 12,
    maximumBreakpointUncertainty: 0.5,
    maximumParentConflict: 0.35,
    maximumRateFold: 5,
    minimumRoleConfidence: 0.5,
    blockSevereWarnings: false,
    revisitReviewed: false,
    consensusWeight: 25,
    significanceWeight: 20,
    informationWeight: 15,
    breakpointWeight: 10,
    diagnosticsWeight: 30,
    roleWeight: 8,
    rescanStrategy: "adaptive",
    rescanRiskThreshold: 65,
    overlapTriggerFraction: 0.4,
    minimumEventsBetweenRescans: 8,
    maximumRescanRounds: 1,
    adaptiveFullTargetFraction: 0.4,
    sameRecombinantRisk: 40,
    recombinantParentRisk: 48,
    groupedEventRisk: 18,
    acceptedWithdrawalRisk: 65,
  },
};

export interface AutoResolveMetrics {
  supportingMethods: number;
  testedMethods: number;
  bestCorrectedP: number;
  informativeSites: number;
  breakpointUncertainty: number;
  parentConflict: number;
  rateFold: number;
  roleConfidence: number;
  roleMatchesCurrent: boolean;
  roleAmbiguous: boolean;
  severeWarnings: number;
}

export interface AutoResolveEntry {
  eventId: string;
  eventIndex: number;
  previousDecision: EventDecision;
  recommendation: AutoResolveRecommendation;
  score: number;
  metrics: AutoResolveMetrics;
  reasons: string[];
  impactRisk: number;
  impactedEventIndexes: number[];
  impactedTargetIndexes: number[];
  impactReasons: string[];
}

export interface AutoResolveRescanBarrier {
  afterEventIndex: number;
  risk: number;
  triggerEventIndexes: number[];
  impactedEventIndexes: number[];
  impactedTargetIndexes: number[];
  reasons: string[];
}

export interface AutoResolvePlan {
  entries: AutoResolveEntry[];
  barriers: AutoResolveRescanBarrier[];
  acceptCount: number;
  rejectCount: number;
  reviewCount: number;
  keepCount: number;
}

export interface AppliedAutoResolvePlan {
  events: RdpEvent[];
  changedIndexes: number[];
  accepted: number;
  rejected: number;
  reviewed: number;
}

const clamp = (value: number, minimum = 0, maximum = 1) => Math.max(minimum, Math.min(maximum, value));

function rateFold(rateRatio: number): number {
  if (!Number.isFinite(rateRatio) || rateRatio <= 0) return 1;
  return Math.max(rateRatio, 1 / rateRatio);
}

function breakpointUncertainty(event: RdpEvent, length: number): number {
  const uncertainty = Math.max(0, event.confidenceStart[1] - event.confidenceStart[0])
    + Math.max(0, event.confidenceEnd[1] - event.confidenceEnd[0]);
  return uncertainty / Math.max(1, 2 * eventLength(event, length));
}

function severeWarningCount(event: RdpEvent): number {
  const pattern = /misalign|artefact|artifact|false[- ]positive|rate[- ]?(?:shift|variation)|diffuse|homoplas|parent conflict|gap block|incomplete|partial candidate/i;
  return event.warnings.filter((warning) => pattern.test(warning)).length;
}

function significanceQuality(bestP: number, maximumP: number): number {
  const boundedP = Math.max(1e-300, Math.min(1, bestP));
  const thresholdStrength = Math.max(1, -Math.log10(Math.max(1e-300, maximumP)));
  return clamp(-Math.log10(boundedP) / (thresholdStrength * 2));
}

function formatP(value: number): string {
  if (!Number.isFinite(value)) return "not available";
  if (value === 0) return "<1×10⁻³⁰⁰";
  return value < 0.001 ? value.toExponential(1) : value.toPrecision(2);
}

function scoreEvent(event: RdpEvent, length: number, settings: AutoResolveSettings): Omit<AutoResolveEntry, "eventIndex" | "impactRisk" | "impactedEventIndexes" | "impactedTargetIndexes" | "impactReasons"> {
  const supportedEvidence = event.evidence.filter((evidence) => evidence.supported && evidence.correctedP <= settings.maximumCorrectedP);
  const supportingMethods = supportedEvidence.length;
  const testedMethods = event.evidence.length;
  const bestCorrectedP = event.evidence.length
    ? Math.min(...event.evidence.map((evidence) => evidence.correctedP))
    : Number.POSITIVE_INFINITY;
  const uncertainty = breakpointUncertainty(event, length);
  const parentConflict = Math.max(0, event.diagnostics.parentConflictRate || 0);
  const fold = rateFold(event.diagnostics.rateRatio);
  const roleIdentification = event.recombinantIdentification;
  const roleConfidence = roleIdentification?.confidence ?? 1;
  const roleMatchesCurrent = !roleIdentification || roleIdentification.recommended === event.recombinant;
  const roleAmbiguous = roleIdentification?.ambiguous ?? false;
  const severeWarnings = severeWarningCount(event);
  const metrics: AutoResolveMetrics = {
    supportingMethods,
    testedMethods,
    bestCorrectedP,
    informativeSites: event.informativeSites,
    breakpointUncertainty: uncertainty,
    parentConflict,
    rateFold: fold,
    roleConfidence,
    roleMatchesCurrent,
    roleAmbiguous,
    severeWarnings,
  };

  if (!settings.revisitReviewed && event.decision !== "unreviewed" && !event.evidenceStale) {
    return {
      eventId: event.id,
      previousDecision: event.decision,
      recommendation: "keep",
      score: event.decision === "accepted" ? 100 : 0,
      metrics,
      reasons: [`Existing ${event.decision} decision is locked by this profile.`],
    };
  }
  if (event.evidenceStale) {
    return {
      eventId: event.id,
      previousDecision: event.decision,
      recommendation: "review",
      score: 0,
      metrics,
      reasons: ["Evidence is stale; automatic acceptance is blocked until recalculation or rescan."],
    };
  }
  if (testedMethods === 0) {
    return {
      eventId: event.id,
      previousDecision: event.decision,
      recommendation: "review",
      score: 0,
      metrics,
      reasons: ["No calibrated method evidence is available for this manual or known-truth hypothesis."],
    };
  }

  const consensusQuality = (clamp(supportingMethods / Math.max(1, settings.minimumSupportingMethods))
    + clamp(supportingMethods / Math.max(1, testedMethods))) / 2;
  const informationQuality = clamp(event.informativeSites / Math.max(1, settings.minimumInformativeSites * 2));
  const precisionQuality = clamp(1 - uncertainty / Math.max(0.0001, settings.maximumBreakpointUncertainty * 2));
  const conflictQuality = clamp(1 - parentConflict / Math.max(0.0001, settings.maximumParentConflict * 2));
  const foldQuality = clamp(1 - (fold - 1) / Math.max(0.0001, (settings.maximumRateFold - 1) * 2));
  const warningQuality = 1 / (1 + severeWarnings);
  const diagnosticQuality = (conflictQuality + foldQuality + (event.diagnostics.diffuseIncompatibility ? 0 : 1) + warningQuality) / 4;
  const roleQuality = roleMatchesCurrent && !roleAmbiguous ? roleConfidence : Math.max(0, 1 - roleConfidence);
  const effectiveRoleWeight = roleIdentification ? settings.roleWeight : 0;
  const weightedTotal = settings.consensusWeight + settings.significanceWeight + settings.informationWeight + settings.breakpointWeight + settings.diagnosticsWeight + effectiveRoleWeight;
  const score = Math.round(100 * (
    consensusQuality * settings.consensusWeight
    + significanceQuality(bestCorrectedP, settings.maximumCorrectedP) * settings.significanceWeight
    + informationQuality * settings.informationWeight
    + precisionQuality * settings.breakpointWeight
    + diagnosticQuality * settings.diagnosticsWeight
    + roleQuality * effectiveRoleWeight
  ) / Math.max(1, weightedTotal));

  const roleGate = !roleIdentification || (
    roleMatchesCurrent
    && !roleAmbiguous
    && roleConfidence >= settings.minimumRoleConfidence
  );

  const acceptanceGates = supportingMethods >= settings.minimumSupportingMethods
    && bestCorrectedP <= settings.maximumCorrectedP
    && event.informativeSites >= settings.minimumInformativeSites
    && uncertainty <= settings.maximumBreakpointUncertainty
    && parentConflict <= settings.maximumParentConflict
    && fold <= settings.maximumRateFold
    && roleGate
    && !event.diagnostics.diffuseIncompatibility
    && (!settings.blockSevereWarnings || severeWarnings === 0);
  const rejectionEvidence = supportingMethods < settings.minimumSupportingMethods
    || bestCorrectedP > settings.maximumCorrectedP
    || parentConflict > settings.maximumParentConflict
    || fold > settings.maximumRateFold
    || event.diagnostics.diffuseIncompatibility
    || (settings.blockSevereWarnings && severeWarnings > 0);
  const recommendation: AutoResolveRecommendation = !roleGate
    ? "review"
    : acceptanceGates && score >= settings.acceptScore
      ? "accept"
      : rejectionEvidence && score <= settings.rejectScore
        ? "reject"
        : "review";
  const reasons = [
    `${supportingMethods}/${testedMethods} methods pass adjusted P ≤ ${settings.maximumCorrectedP}.`,
    `Best adjusted P is ${formatP(bestCorrectedP)}; ${event.informativeSites} informative sites.`,
    `Breakpoint uncertainty is ${(uncertainty * 100).toFixed(1)}% of tract length.`,
    `Parent conflict ${(parentConflict * 100).toFixed(1)}%; rate-density deviation ${fold.toFixed(2)}×.`,
  ];
  if (roleIdentification) reasons.push(`Source role consensus is ${(roleConfidence * 100).toFixed(1)}% for ${roleMatchesCurrent ? "the current recombinant" : "another triplet member"}${roleAmbiguous ? " and is flagged ambiguous" : ""}.`);
  if (event.diagnostics.diffuseIncompatibility) reasons.push("Diffuse incompatibility diagnostic is active.");
  if (severeWarnings) reasons.push(`${severeWarnings} high-risk warning${severeWarnings === 1 ? "" : "s"} matched the configured blocker.`);
  if (recommendation === "review") reasons.push("The score falls inside the analyst-review band or a hard acceptance gate failed.");
  return { eventId: event.id, previousDecision: event.decision, recommendation, score, metrics, reasons };
}

interface DependencyIndex {
  byRecombinant: Map<number, number[]>;
  byParent: Map<number, number[]>;
  byGroup: Map<string, number[]>;
  intervalBuckets: Map<number, Map<number, number[]>>;
  intervalBinCount: number;
  length: number;
}

function appendIndex<K>(map: Map<K, number[]>, key: K, eventIndex: number): void {
  const indexes = map.get(key);
  if (indexes) indexes.push(eventIndex);
  else map.set(key, [eventIndex]);
}

function eventIntervalBins(event: RdpEvent, length: number, binCount: number): number[] {
  if (length <= 0) return [0];
  const intervals: Array<[number, number]> = event.wraps ? [[event.start, length], [0, event.end]] : [[event.start, event.end]];
  const bins = new Set<number>();
  intervals.forEach(([start, end]) => {
    if (end <= start) return;
    const first = Math.max(0, Math.min(binCount - 1, Math.floor(start / length * binCount)));
    const last = Math.max(first, Math.min(binCount - 1, Math.ceil(end / length * binCount) - 1));
    for (let bin = first; bin <= last; bin += 1) bins.add(bin);
  });
  return [...bins];
}

function buildDependencyIndex(events: RdpEvent[], length: number): DependencyIndex {
  const intervalBinCount = Math.max(32, Math.min(256, Math.ceil(Math.sqrt(Math.max(1, events.length)) * 2)));
  const index: DependencyIndex = { byRecombinant: new Map(), byParent: new Map(), byGroup: new Map(), intervalBuckets: new Map(), intervalBinCount, length };
  events.forEach((event, eventIndex) => {
    appendIndex(index.byRecombinant, event.recombinant, eventIndex);
    appendIndex(index.byParent, event.majorParent, eventIndex);
    appendIndex(index.byParent, event.minorParent, eventIndex);
    if (event.groupId) appendIndex(index.byGroup, event.groupId, eventIndex);
    let sequenceBuckets = index.intervalBuckets.get(event.recombinant);
    if (!sequenceBuckets) {
      sequenceBuckets = new Map();
      index.intervalBuckets.set(event.recombinant, sequenceBuckets);
    }
    eventIntervalBins(event, length, intervalBinCount).forEach((bin) => appendIndex(sequenceBuckets!, bin, eventIndex));
  });
  return index;
}

function overlapCandidates(event: RdpEvent, dependencyIndex: DependencyIndex): number[] {
  const buckets = dependencyIndex.intervalBuckets.get(event.recombinant);
  if (!buckets) return [];
  const candidates = new Set<number>();
  eventIntervalBins(event, dependencyIndex.length, dependencyIndex.intervalBinCount).forEach((bin) => {
    (buckets.get(bin) ?? []).forEach((eventIndex) => candidates.add(eventIndex));
  });
  return [...candidates];
}

function downstreamImpact(
  event: RdpEvent,
  eventIndex: number,
  recommendation: AutoResolveRecommendation,
  events: RdpEvent[],
  length: number,
  settings: AutoResolveSettings,
  dependencyIndex: DependencyIndex,
): Pick<AutoResolveEntry, "impactRisk" | "impactedEventIndexes" | "impactedTargetIndexes" | "impactReasons"> {
  const becomesAccepted = recommendation === "accept" && event.decision !== "accepted";
  const withdrawsAcceptance = recommendation === "reject" && event.decision === "accepted";
  if (!becomesAccepted && !withdrawsAcceptance) return { impactRisk: 0, impactedEventIndexes: [], impactedTargetIndexes: [], impactReasons: [] };
  const impactedEvents = new Set<number>();
  const impactedTargets = new Set<number>();
  const reasons = new Set<string>();
  let sameRecombinantFraction = 0;
  let recombinantParent = false;
  let reciprocalNestedParent = false;
  let grouped = false;
  const candidateIndexes = new Set<number>([
    ...overlapCandidates(event, dependencyIndex),
    ...(dependencyIndex.byParent.get(event.recombinant) ?? []),
    ...(dependencyIndex.byRecombinant.get(event.majorParent) ?? []),
    ...(dependencyIndex.byRecombinant.get(event.minorParent) ?? []),
    ...(event.groupId ? dependencyIndex.byGroup.get(event.groupId) ?? [] : []),
  ]);
  for (const laterIndex of candidateIndexes) {
    if (laterIndex <= eventIndex) continue;
    const later = events[laterIndex];
    if (later.decision === "rejected" && !settings.revisitReviewed) continue;
    let affected = false;
    if (later.recombinant === event.recombinant) {
      const overlap = eventOverlapBases(event, later, length);
      const overlapFraction = overlap / Math.max(1, Math.min(eventLength(event, length), eventLength(later, length)));
      if (overlapFraction >= settings.overlapTriggerFraction) {
        sameRecombinantFraction = Math.max(sameRecombinantFraction, overlapFraction);
        reasons.add(`overlapping later tract in the same recombinant (${Math.round(overlapFraction * 100)}%)`);
        affected = true;
      }
    }
    if (later.majorParent === event.recombinant || later.minorParent === event.recombinant) {
      recombinantParent = true;
      reasons.add("the resolved recombinant is used as a parent proxy later");
      affected = true;
    }
    if (event.majorParent === later.recombinant || event.minorParent === later.recombinant) {
      reciprocalNestedParent = true;
      reasons.add("a later recombinant was used as a parent proxy in the earlier event");
      affected = true;
    }
    if (event.groupId && event.groupId === later.groupId) {
      grouped = true;
      reasons.add("a later hypothesis shares the same co-recombinant group");
      affected = true;
    }
    if (affected) {
      impactedEvents.add(laterIndex);
      impactedTargets.add(later.recombinant);
    }
  }
  let risk = withdrawsAcceptance ? settings.acceptedWithdrawalRisk : 0;
  if (sameRecombinantFraction > 0) risk += settings.sameRecombinantRisk * Math.max(0.35, sameRecombinantFraction);
  if (recombinantParent) risk += settings.recombinantParentRisk;
  if (reciprocalNestedParent) risk += settings.recombinantParentRisk * 0.55;
  if (grouped) risk += settings.groupedEventRisk;
  return {
    impactRisk: Math.round(clamp(risk, 0, 100)),
    impactedEventIndexes: [...impactedEvents].sort((left, right) => left - right),
    impactedTargetIndexes: [...impactedTargets].sort((left, right) => left - right),
    impactReasons: [...reasons],
  };
}

function combineRisk(left: number, right: number): number {
  return Math.round(100 * (1 - (1 - clamp(left / 100)) * (1 - clamp(right / 100))));
}

export function planAutoResolution(events: RdpEvent[], length: number, settings: AutoResolveSettings): AutoResolvePlan {
  const dependencyIndex = buildDependencyIndex(events, length);
  let entries: AutoResolveEntry[] = events.map((event, eventIndex) => {
    const scored = scoreEvent(event, length, settings);
    return {
      ...scored,
      eventIndex,
      impactRisk: 0,
      impactedEventIndexes: [],
      impactedTargetIndexes: [],
      impactReasons: [],
    };
  });
  const blockerQueue = entries
    .filter((entry) => entry.recommendation === "review" && events[entry.eventIndex].decision !== "rejected")
    .map((entry) => entry.eventIndex);
  const propagatedBlockers = new Set<number>();
  while (blockerQueue.length) {
    const blockerIndex = blockerQueue.shift()!;
    if (propagatedBlockers.has(blockerIndex)) continue;
    propagatedBlockers.add(blockerIndex);
    const potential = downstreamImpact({ ...events[blockerIndex], decision: "unreviewed" }, blockerIndex, "accept", events, length, settings, dependencyIndex);
    if (potential.impactRisk < settings.rescanRiskThreshold) continue;
    const impacted = new Set(potential.impactedEventIndexes);
    entries = entries.map((candidate, candidateIndex) => {
      if (!impacted.has(candidateIndex) || candidate.recommendation === "keep") return candidate;
      if (candidate.recommendation !== "review") blockerQueue.push(candidateIndex);
      return {
        ...candidate,
        recommendation: "review",
        impactRisk: 0,
        impactedEventIndexes: [],
        impactedTargetIndexes: [],
        impactReasons: [],
        reasons: [...candidate.reasons, `Automatic decision blocked because unresolved E${blockerIndex + 1} has ${potential.impactRisk}/100 downstream dependency risk.`],
      };
    });
  }
  const barriers: AutoResolveRescanBarrier[] = [];
  if (settings.rescanStrategy !== "off") {
    let pendingRisk = 0;
    let pendingFirstIndex = -1;
    let pendingLastIndex = -1;
    let pendingTriggers = new Set<number>();
    let pendingEvents = new Set<number>();
    let pendingTargets = new Set<number>();
    let pendingReasons = new Set<string>();
    const flush = (afterEventIndex: number) => {
      barriers.push({
        afterEventIndex,
        risk: pendingRisk,
        triggerEventIndexes: [...pendingTriggers],
        impactedEventIndexes: [...pendingEvents].sort((left, right) => left - right),
        impactedTargetIndexes: [...pendingTargets].sort((left, right) => left - right),
        reasons: [...pendingReasons],
      });
      pendingRisk = 0;
      pendingFirstIndex = -1;
      pendingLastIndex = -1;
      pendingTriggers = new Set();
      pendingEvents = new Set();
      pendingTargets = new Set();
      pendingReasons = new Set();
    };
    for (let entryIndex = 0; entryIndex < entries.length; entryIndex += 1) {
      const entry = entries[entryIndex];
      if (entry.recommendation === "accept" || entry.recommendation === "reject") {
        const impact = downstreamImpact(events[entryIndex], entryIndex, entry.recommendation, events, length, settings, dependencyIndex);
        entries[entryIndex] = { ...entry, ...impact };
        if (impact.impactRisk > 0 && impact.impactedTargetIndexes.length > 0) {
          if (pendingFirstIndex < 0) pendingFirstIndex = entryIndex;
          pendingRisk = combineRisk(pendingRisk, impact.impactRisk);
          pendingLastIndex = entryIndex;
          pendingTriggers.add(entryIndex);
          impact.impactedEventIndexes.forEach((index) => pendingEvents.add(index));
          impact.impactedTargetIndexes.forEach((index) => pendingTargets.add(index));
          impact.impactReasons.forEach((reason) => pendingReasons.add(reason));
        }
      }
      const enoughBatching = pendingFirstIndex >= 0 && entryIndex - pendingFirstIndex + 1 >= settings.minimumEventsBetweenRescans;
      if (pendingRisk >= settings.rescanRiskThreshold && enoughBatching) {
        flush(entryIndex);
        break;
      }
    }
    if (barriers.length === 0 && pendingRisk >= settings.rescanRiskThreshold && pendingLastIndex >= 0) flush(pendingLastIndex);
  }
  return {
    entries,
    barriers,
    acceptCount: entries.filter((entry) => entry.recommendation === "accept").length,
    rejectCount: entries.filter((entry) => entry.recommendation === "reject").length,
    reviewCount: entries.filter((entry) => entry.recommendation === "review").length,
    keepCount: entries.filter((entry) => entry.recommendation === "keep").length,
  };
}

export function applyAutoResolutionPlan(
  events: RdpEvent[],
  plan: AutoResolvePlan,
  throughEventIndex: number,
  profileLabel: string,
  timestamp: string,
): AppliedAutoResolvePlan {
  const entries = new Map(plan.entries.map((entry) => [entry.eventId, entry]));
  const changedIndexes: number[] = [];
  let accepted = 0;
  let rejected = 0;
  let reviewed = 0;
  const nextEvents = events.map((event, eventIndex) => {
    if (eventIndex > throughEventIndex) return event;
    const entry = entries.get(event.id);
    if (!entry || entry.recommendation === "keep") return event;
    if (entry.recommendation === "review") {
      return event;
    }
    const decision: EventDecision = entry.recommendation === "accept" ? "accepted" : "rejected";
    if (event.decision === decision) return event;
    if (decision === "accepted") accepted += 1;
    else rejected += 1;
    changedIndexes.push(eventIndex);
    return {
      ...event,
      decision,
      history: [...event.history, {
        id: `auto-resolve-${timestamp}-${eventIndex}`,
        timestamp,
        action: `Auto-${entry.recommendation}ed event`,
        summary: `${profileLabel} heuristic score ${entry.score}/100. ${entry.reasons.slice(0, 2).join(" ")}`,
      }],
    };
  });
  reviewed = plan.entries.filter((entry) => entry.eventIndex <= throughEventIndex && entry.recommendation === "review").length;
  return { events: nextEvents, changedIndexes, accepted, rejected, reviewed };
}

export function rescanTargetsForBarrier(
  barrier: AutoResolveRescanBarrier,
  sequenceCount: number,
  settings: AutoResolveSettings,
): { targetIndexes: number[]; scope: "targeted" | "full" } {
  const impacted = [...new Set(barrier.impactedTargetIndexes.filter((index) => index >= 0 && index < sequenceCount))];
  const useFull = settings.rescanStrategy === "full"
    || (settings.rescanStrategy === "adaptive" && impacted.length / Math.max(1, sequenceCount) >= settings.adaptiveFullTargetFraction);
  return useFull
    ? { targetIndexes: Array.from({ length: sequenceCount }, (_, index) => index), scope: "full" }
    : { targetIndexes: impacted, scope: "targeted" };
}

function circularBoundaryDistance(left: number, right: number, length: number): number {
  const direct = Math.abs(left - right);
  return length > 0 ? Math.min(direct, Math.max(0, length - direct)) : direct;
}

function sameParentPair(left: RdpEvent, right: RdpEvent): boolean {
  return (left.majorParent === right.majorParent && left.minorParent === right.minorParent)
    || (left.majorParent === right.minorParent && left.minorParent === right.majorParent);
}

export function filterResolvedEventDuplicates(candidates: RdpEvent[], resolved: RdpEvent[], length: number): RdpEvent[] {
  const byRecombinant = new Map<number, RdpEvent[]>();
  resolved.forEach((event) => appendIndexEvent(byRecombinant, event.recombinant, event));
  const tolerance = Math.max(4, Math.round(length * 0.01));
  return candidates.filter((candidate) => !(byRecombinant.get(candidate.recombinant) ?? []).some((event) => {
    if (!sameParentPair(candidate, event)) return false;
    const closeBoundaries = circularBoundaryDistance(candidate.start, event.start, length) <= tolerance
      && circularBoundaryDistance(candidate.end, event.end, length) <= tolerance;
    const candidateLength = eventLength(candidate, length);
    const resolvedLength = eventLength(event, length);
    const reciprocalOverlap = eventOverlapBases(candidate, event, length) / Math.max(1, Math.max(candidateLength, resolvedLength));
    const comparableLengths = Math.min(candidateLength, resolvedLength) / Math.max(1, Math.max(candidateLength, resolvedLength)) >= 0.7;
    return closeBoundaries || (comparableLengths && reciprocalOverlap >= 0.85);
  }));
}

function appendIndexEvent<K, V>(map: Map<K, V[]>, key: K, value: V): void {
  const values = map.get(key);
  if (values) values.push(value);
  else map.set(key, [value]);
}
