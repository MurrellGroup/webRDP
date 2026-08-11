import {
  AlignmentData,
  EventDecision,
  RdpEvent,
  eventLength,
  eventSegments,
} from "./rdp-core";

export interface ReviewQueueFilters {
  skipAccepted: boolean;
  skipRejected: boolean;
  minimumMethods: number;
  warningsOnly: boolean;
  staleOnly: boolean;
}

export interface ReviewChecklistItem {
  key: "methods" | "breakpoints" | "roles" | "artefacts" | "grouping";
  label: string;
  detail: string;
  state: "pass" | "review" | "fail";
}

export interface RoleAssignmentTrial {
  key: "current" | "swap-major" | "swap-minor";
  recombinant: number;
  majorParent: number;
  minorParent: number;
  tractMinorIdentity: number;
  tractMajorIdentity: number;
  backgroundMajorIdentity: number;
  backgroundMinorIdentity: number;
  switchSignal: number;
  score: number;
  informativeSites: number;
  sourcePoints?: number;
  sourceShare?: number;
  sourceRecommended?: boolean;
  sourceAmbiguous?: boolean;
  sourceTests?: Array<{
    id: string;
    label: string;
    sourceRoutine: string;
    value: number | null;
    points: number;
    winner: boolean;
    decisive: boolean;
  }>;
}

export const DEFAULT_REVIEW_FILTERS: ReviewQueueFilters = {
  skipAccepted: true,
  skipRejected: true,
  minimumMethods: 0,
  warningsOnly: false,
  staleOnly: false,
};

export function supportingMethodCount(event: RdpEvent): number {
  return event.evidence.reduce((count, item) => count + Number(item.supported), 0);
}

export function bestCorrectedP(event: RdpEvent): number {
  return event.evidence.length
    ? Math.min(...event.evidence.map((item) => item.correctedP))
    : Number.POSITIVE_INFINITY;
}

export function eventMatchesReviewFilters(event: RdpEvent, filters: ReviewQueueFilters): boolean {
  if (filters.skipAccepted && event.decision === "accepted") return false;
  if (filters.skipRejected && event.decision === "rejected") return false;
  if (supportingMethodCount(event) < filters.minimumMethods) return false;
  if (filters.warningsOnly && event.warnings.length === 0) return false;
  if (filters.staleOnly && !event.evidenceStale) return false;
  return true;
}

export function filteredReviewIndexes(events: RdpEvent[], filters: ReviewQueueFilters): number[] {
  const indexes: number[] = [];
  events.forEach((event, index) => {
    if (eventMatchesReviewFilters(event, filters)) indexes.push(index);
  });
  return indexes;
}

export function navigateReviewEvent(
  events: RdpEvent[],
  selectedId: string | null,
  direction: -1 | 1,
  filters: ReviewQueueFilters,
): string | null {
  const candidates = filteredReviewIndexes(events, filters);
  if (!candidates.length) return null;
  const selectedIndex = events.findIndex((event) => event.id === selectedId);
  if (selectedIndex < 0) return events[candidates[direction > 0 ? 0 : candidates.length - 1]].id;
  if (direction > 0) {
    const next = candidates.find((index) => index > selectedIndex) ?? candidates[0];
    return events[next].id;
  }
  const previous = [...candidates].reverse().find((index) => index < selectedIndex) ?? candidates.at(-1)!;
  return events[previous].id;
}

export function bestUnresolvedEventId(events: RdpEvent[], filters: ReviewQueueFilters): string | null {
  let best: RdpEvent | null = null;
  let bestP = Number.POSITIVE_INFINITY;
  events.forEach((event) => {
    if (event.decision !== "unreviewed") return;
    if (!eventMatchesReviewFilters(event, filters)) return;
    const correctedP = bestCorrectedP(event);
    if (!best || correctedP < bestP) {
      best = event;
      bestP = correctedP;
    }
  });
  return best ? (best as RdpEvent).id : null;
}

export function eventGroupIndexes(events: RdpEvent[], event: RdpEvent): number[] {
  if (!event.groupId) {
    const index = events.findIndex((candidate) => candidate.id === event.id);
    return index < 0 ? [] : [index];
  }
  return events.flatMap((candidate, index) => candidate.groupId === event.groupId ? [index] : []);
}

function warningMatches(event: RdpEvent, expression: RegExp): boolean {
  return event.warnings.some((warning) => expression.test(warning));
}

export function buildReviewChecklist(event: RdpEvent, alignmentLength: number): ReviewChecklistItem[] {
  const supported = supportingMethodCount(event);
  const tractLength = Math.max(1, eventLength(event, alignmentLength));
  const uncertainty = (
    Math.max(0, event.confidenceStart[1] - event.confidenceStart[0])
    + Math.max(0, event.confidenceEnd[1] - event.confidenceEnd[0])
  ) / tractLength;
  const roleWarning = warningMatches(event, /wrong (?:sequence|recombinant)|recombinant.*uncertain|parent.*(?:uncertain|unsampled)/i);
  const artefactWarning = warningMatches(event, /misalign|artefact|rate shift|homoplas|diffuse incompat|false positive/i);
  return [
    {
      key: "methods",
      label: "Method confirmation",
      detail: supported >= 2 ? `${supported} independent method families support the event.` : `${supported} supporting method${supported === 1 ? "" : "s"}; seek a second line of evidence.`,
      state: supported >= 2 ? "pass" : supported === 1 ? "review" : "fail",
    },
    {
      key: "breakpoints",
      label: "Breakpoint placement",
      detail: event.evidenceStale ? "The event was edited and its saved breakpoint evidence is stale." : uncertainty <= 0.12 ? "Confidence bounds are narrow relative to the proposed tract." : `Combined uncertainty spans ${Math.round(uncertainty * 100)}% of the tract.`,
      state: event.evidenceStale ? "fail" : uncertainty <= 0.12 ? "pass" : "review",
    },
    {
      key: "roles",
      label: "Recombinant polarity",
      detail: roleWarning ? "The source diagnostics warn that the recombinant or a parent proxy may be misassigned." : "No saved diagnostic specifically challenges the current role assignment.",
      state: roleWarning ? "review" : "pass",
    },
    {
      key: "artefacts",
      label: "False-positive challenge",
      detail: artefactWarning ? "At least one alignment, rate-variation, homoplasy, or diffuse-incompatibility warning needs review." : "No high-risk artefact warning is attached to this event.",
      state: artefactWarning ? "review" : "pass",
    },
    {
      key: "grouping",
      label: "Co-recombinant grouping",
      detail: event.groupId ? `Assigned to ${event.groupId}; verify that its descendants move together across local trees.` : "No common-ancestor group is asserted; leave ungrouped unless tree and plot comparisons support one.",
      state: "review",
    },
  ];
}

function complementarySegments(event: RdpEvent, length: number): [number, number][] {
  if (event.wraps) return event.end < event.start ? [[event.end, event.start]] : [];
  const segments: [number, number][] = [];
  if (event.start > 0) segments.push([0, event.start]);
  if (event.end < length) segments.push([event.end, length]);
  return segments;
}

function sampledIdentity(
  left: string,
  right: string,
  segments: [number, number][],
  maximumSites = 4096,
): { identity: number; callable: number } {
  const totalLength = segments.reduce((sum, [start, end]) => sum + Math.max(0, end - start), 0);
  const stride = Math.max(1, Math.ceil(totalLength / maximumSites));
  let callable = 0;
  let matches = 0;
  for (const [start, end] of segments) {
    for (let site = start; site < end; site += stride) {
      const leftBase = left[site];
      const rightBase = right[site];
      if (!"ACGT".includes(leftBase) || !"ACGT".includes(rightBase)) continue;
      callable += 1;
      if (leftBase === rightBase) matches += 1;
    }
  }
  return { identity: callable ? matches / callable : 0, callable };
}

function informativeContrastCount(
  recombinant: string,
  major: string,
  minor: string,
  segments: [number, number][],
  maximumSites = 4096,
): number {
  const totalLength = segments.reduce((sum, [start, end]) => sum + Math.max(0, end - start), 0);
  const stride = Math.max(1, Math.ceil(totalLength / maximumSites));
  let count = 0;
  for (const [start, end] of segments) {
    for (let site = start; site < end; site += stride) {
      const recombinantBase = recombinant[site];
      const majorBase = major[site];
      const minorBase = minor[site];
      if (![recombinantBase, majorBase, minorBase].every((base) => "ACGT".includes(base))) continue;
      if (majorBase !== minorBase && (recombinantBase === majorBase || recombinantBase === minorBase)) count += 1;
    }
  }
  return count;
}

export function roleAssignmentTrials(alignment: AlignmentData, event: RdpEvent): RoleAssignmentTrial[] {
  const source = event.recombinantIdentification;
  const sourceOrientation = (recombinant: number) => source?.orientations.find((orientation) => orientation.recombinant === recombinant);
  const majorOrientation = sourceOrientation(event.majorParent);
  const minorOrientation = sourceOrientation(event.minorParent);
  const assignments: Array<Pick<RoleAssignmentTrial, "key" | "recombinant" | "majorParent" | "minorParent">> = [
    { key: "current", recombinant: event.recombinant, majorParent: event.majorParent, minorParent: event.minorParent },
    { key: "swap-major", recombinant: event.majorParent, majorParent: majorOrientation?.majorParent ?? event.recombinant, minorParent: majorOrientation?.minorParent ?? event.minorParent },
    { key: "swap-minor", recombinant: event.minorParent, majorParent: minorOrientation?.majorParent ?? event.majorParent, minorParent: minorOrientation?.minorParent ?? event.recombinant },
  ];
  const tract = eventSegments(event, alignment.length);
  const background = complementarySegments(event, alignment.length);
  return assignments.map((assignment) => {
    const recombinant = alignment.sequences[assignment.recombinant]?.sequence ?? "";
    const major = alignment.sequences[assignment.majorParent]?.sequence ?? "";
    const minor = alignment.sequences[assignment.minorParent]?.sequence ?? "";
    const tractMinor = sampledIdentity(recombinant, minor, tract);
    const tractMajor = sampledIdentity(recombinant, major, tract);
    const backgroundMajor = sampledIdentity(recombinant, major, background);
    const backgroundMinor = sampledIdentity(recombinant, minor, background);
    const switchSignal = (
      tractMinor.identity - tractMajor.identity
      + backgroundMajor.identity - backgroundMinor.identity
    ) / 2;
    const orientation = sourceOrientation(assignment.recombinant);
    const candidateIndex = source?.candidates.indexOf(assignment.recombinant) ?? -1;
    return {
      ...assignment,
      tractMinorIdentity: tractMinor.identity,
      tractMajorIdentity: tractMajor.identity,
      backgroundMajorIdentity: backgroundMajor.identity,
      backgroundMinorIdentity: backgroundMinor.identity,
      switchSignal,
      score: orientation?.sourceScore ?? Math.round(Math.max(0, Math.min(100, 50 + switchSignal * 250))),
      informativeSites: informativeContrastCount(recombinant, major, minor, [...tract, ...background]),
      sourcePoints: orientation?.sourcePoints,
      sourceShare: orientation?.sourceShare,
      sourceRecommended: source?.recommended === assignment.recombinant,
      sourceAmbiguous: source?.ambiguous,
      sourceTests: candidateIndex >= 0 ? source?.tests.map((test) => ({
        id: test.id,
        label: test.label,
        sourceRoutine: test.sourceRoutine,
        value: test.values[candidateIndex] ?? null,
        points: test.points[candidateIndex] ?? 0,
        winner: test.winnerIndexes.includes(candidateIndex),
        decisive: test.decisive,
      })) : undefined,
    };
  });
}

export function applyDecisionToIndexes(events: RdpEvent[], indexes: number[], decision: EventDecision): RdpEvent[] {
  const targets = new Set(indexes);
  return events.map((event, index) => targets.has(index) ? { ...event, decision } : event);
}
