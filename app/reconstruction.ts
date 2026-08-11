import type { RdpEvent } from "./rdp-core";

export type ReconstructionRelationshipKind = "possible-overprint" | "recombinant-parent" | "event-group";

export interface ReconstructionRelationship {
  kind: ReconstructionRelationshipKind;
  fromIndex: number;
  toIndex: number;
  overlapBases: number;
}

export interface ReconstructionSequenceRow {
  sequenceIndex: number;
  eventIndexes: number[];
}

export interface ReconstructionModel {
  retainedIndexes: number[];
  relationships: ReconstructionRelationship[];
  relationshipsTruncated: boolean;
  sequenceRows: ReconstructionSequenceRow[];
  nextReviewIndex: number | null;
  staleFromIndex: number | null;
  downstreamIndexes: number[];
}

function segments(event: RdpEvent, length: number): Array<[number, number]> {
  if (!event.wraps) return [[event.start, event.end]];
  return [[event.start, length], [0, event.end]];
}

export function eventOverlapBases(left: RdpEvent, right: RdpEvent, length: number): number {
  let overlap = 0;
  for (const [leftStart, leftEnd] of segments(left, length)) {
    for (const [rightStart, rightEnd] of segments(right, length)) {
      overlap += Math.max(0, Math.min(leftEnd, rightEnd) - Math.max(leftStart, rightStart));
    }
  }
  return overlap;
}

export function buildReconstructionModel(
  events: RdpEvent[],
  length: number,
): ReconstructionModel {
  const retainedIndexes = events
    .map((event, index) => ({ event, index }))
    .filter(({ event }) => event.decision !== "rejected")
    .map(({ index }) => index);
  const relationships: ReconstructionRelationship[] = [];
  const maximumRelationships = 20_000;
  let relationshipsTruncated = false;
  const relationshipKeys = new Set<string>();
  const byRecombinant = new Map<number, number[]>();
  const byGroup = new Map<string, number[]>();
  retainedIndexes.forEach((eventIndex) => {
    const event = events[eventIndex];
    byRecombinant.set(event.recombinant, [...(byRecombinant.get(event.recombinant) ?? []), eventIndex]);
    if (event.groupId) byGroup.set(event.groupId, [...(byGroup.get(event.groupId) ?? []), eventIndex]);
  });
  const addRelationship = (kind: ReconstructionRelationshipKind, fromIndex: number, toIndex: number, knownOverlap?: number) => {
    if (fromIndex === toIndex) return;
    if (relationships.length >= maximumRelationships) {
      relationshipsTruncated = true;
      return;
    }
    const key = `${kind}:${fromIndex}:${toIndex}`;
    if (relationshipKeys.has(key)) return;
    relationshipKeys.add(key);
    relationships.push({ kind, fromIndex, toIndex, overlapBases: knownOverlap ?? eventOverlapBases(events[fromIndex], events[toIndex], length) });
  };
  byRecombinant.forEach((eventIndexes) => {
    if (relationships.length >= maximumRelationships) { relationshipsTruncated = true; return; }
    for (let left = 0; left < eventIndexes.length && relationships.length < maximumRelationships; left += 1) {
      for (let right = left + 1; right < eventIndexes.length && relationships.length < maximumRelationships; right += 1) {
        const leftIndex = eventIndexes[left];
        const rightIndex = eventIndexes[right];
        const overlap = eventOverlapBases(events[leftIndex], events[rightIndex], length);
        if (overlap > 0) addRelationship("possible-overprint", leftIndex, rightIndex, overlap);
      }
    }
  });
  byGroup.forEach((eventIndexes) => {
    if (relationships.length >= maximumRelationships) { relationshipsTruncated = true; return; }
    for (let left = 0; left < eventIndexes.length && relationships.length < maximumRelationships; left += 1) {
      for (let right = left + 1; right < eventIndexes.length && relationships.length < maximumRelationships; right += 1) addRelationship("event-group", eventIndexes[left], eventIndexes[right]);
    }
  });
  retainedIndexes.forEach((eventIndex) => {
    if (relationships.length >= maximumRelationships) { relationshipsTruncated = true; return; }
    const event = events[eventIndex];
    const parentSequences = new Set([event.majorParent, event.minorParent]);
    parentSequences.forEach((parentSequence) => {
      (byRecombinant.get(parentSequence) ?? []).forEach((carrierIndex) => addRelationship("recombinant-parent", carrierIndex, eventIndex));
    });
  });
  const sequenceRows = [...byRecombinant]
    .map(([sequenceIndex, eventIndexes]) => ({ sequenceIndex, eventIndexes }))
    .sort((left, right) => left.sequenceIndex - right.sequenceIndex);
  const staleFromIndex = events.findIndex((event) => event.decision !== "rejected" && event.evidenceStale);
  const normalizedStaleFrom = staleFromIndex >= 0 ? staleFromIndex : null;
  const downstreamIndexes = normalizedStaleFrom === null
    ? []
    : retainedIndexes.filter((index) => index > normalizedStaleFrom);
  const nextReviewIndex = events.findIndex((event) => event.decision === "unreviewed" || event.evidenceStale);
  return {
    retainedIndexes,
    relationships,
    relationshipsTruncated,
    sequenceRows,
    nextReviewIndex: nextReviewIndex >= 0 ? nextReviewIndex : null,
    staleFromIndex: normalizedStaleFrom,
    downstreamIndexes,
  };
}
