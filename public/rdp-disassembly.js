function eventSegments(event, length) {
  if (event.wraps === true && event.start > event.end) {
    return [[event.start, length], ...(event.end > 0 ? [[0, event.end]] : [])];
  }
  return event.end > event.start ? [[event.start, event.end]] : [];
}

function sameLineage(left, right) {
  if (!Array.isArray(left) || !Array.isArray(right) || left.length !== right.length) return false;
  return left.every((value, index) => value === right[index]);
}

function selectedCoRecombinants(event, sequenceCount) {
  const selected = event.coRecombinantSets?.find((set) => set.presumedRecombinant === event.recombinant);
  const members = selected?.sequenceMembers?.length
    ? selected.sequenceMembers
    : event.ancestralCluster?.sequenceMembers?.length
      ? event.ancestralCluster.sequenceMembers
      : [event.recombinant];
  return [...new Set(members)]
    .map((value) => Math.trunc(value))
    .filter((value) => value >= 0 && value < sequenceCount);
}

// RDP5 manual section 4.1.6 erases an accepted recombinant tract from the
// working sequence and adds that tract as a separate, gap-padded sequence.
// Repeating the operation against the event's recorded component lineage lets
// later scans expose nested and overprinted signals without presenting private
// component rows as extra user sequences.
export function buildDisassembledAlignment(originalEncoded, sequences, length, acceptedEvents = []) {
  const originalCount = sequences.length;
  const rows = Array.from({ length: originalCount }, (_, originIndex) => ({
    data: originalEncoded.slice(originIndex * length, (originIndex + 1) * length),
    structural: new Uint8Array(length),
    mapping: { originIndex, kind: "remainder", lineage: [], erasedEventIds: [] },
  }));
  const appliedEventIds = [];
  let erasedCanonicalBases = 0;

  for (const event of acceptedEvents) {
    if (!event || event.decision !== "accepted" || event.evidenceStale === true || typeof event.id !== "string") continue;
    const segments = eventSegments(event, length);
    if (!segments.length) continue;
    const requestedLineage = Array.isArray(event.componentProvenance?.recombinant?.lineage)
      ? event.componentProvenance.recombinant.lineage.filter((value) => typeof value === "string")
      : [];
    const members = selectedCoRecombinants(event, originalCount);
    let applied = false;
    for (const originIndex of members) {
      let target = rows.find((row) => row.mapping.originIndex === originIndex
        && sameLineage(row.mapping.lineage, requestedLineage));
      if (!target) {
        target = rows.find((row) => row.mapping.originIndex === originIndex
          && row.mapping.kind === "remainder");
      }
      if (!target) continue;

      const component = new Uint8Array(length);
      component.fill(4);
      const componentStructural = new Uint8Array(length);
      componentStructural.fill(1);
      let copied = 0;
      for (const [start, end] of segments) {
        for (let site = start; site < end; site += 1) {
          const base = target.data[site];
          component[site] = base;
          componentStructural[site] = target.structural[site];
          if (base < 4) {
            copied += 1;
            erasedCanonicalBases += 1;
          }
          target.data[site] = 4;
          target.structural[site] = 1;
        }
      }
      if (!target.mapping.erasedEventIds.includes(event.id)) target.mapping.erasedEventIds.push(event.id);
      if (copied < 1) continue;
      rows.push({
        data: component,
        structural: componentStructural,
        mapping: {
          originIndex,
          kind: "extracted-tract",
          lineage: [...target.mapping.lineage, event.id],
          sourceEventId: event.id,
          parentLineage: [...target.mapping.lineage],
          start: event.start,
          end: event.end,
          wraps: event.wraps === true,
          erasedEventIds: [],
        },
      });
      applied = true;
    }
    if (applied) appliedEventIds.push(event.id);
  }

  const encoded = new Uint8Array(rows.length * length);
  const structuralMasks = new Uint8Array(rows.length * length);
  const analysisSequences = rows.map((row, index) => {
    encoded.set(row.data, index * length);
    structuralMasks.set(row.structural, index * length);
    const original = sequences[row.mapping.originIndex];
    return row.mapping.kind === "remainder"
      ? original
      : { ...original, name: `${original.name} · component ${row.mapping.lineage.length}` };
  });
  return {
    encoded,
    analysisSequences,
    mappings: rows.map((row) => row.mapping),
    structuralMasks,
    appliedEventIds,
    componentCount: rows.length - originalCount,
    erasedCanonicalBases,
  };
}

export function componentReference(mapping) {
  if (!mapping) return undefined;
  return {
    originIndex: mapping.originIndex,
    kind: mapping.kind,
    lineage: [...mapping.lineage],
    sourceEventId: mapping.sourceEventId,
    parentLineage: mapping.parentLineage ? [...mapping.parentLineage] : undefined,
    start: mapping.start,
    end: mapping.end,
    wraps: mapping.wraps,
    erasedEventIds: [...(mapping.erasedEventIds ?? [])],
  };
}

export function findComponentIndex(disassembly, reference) {
  if (!reference) return -1;
  const requestedLineage = Array.isArray(reference.lineage) ? reference.lineage : [];
  return disassembly.mappings.findIndex((mapping) => mapping.originIndex === reference.originIndex
    && mapping.kind === reference.kind
    && sameLineage(mapping.lineage, requestedLineage));
}

export function candidateComponentProvenance(disassembly, recombinant, majorParent, minorParent) {
  if (!disassembly.appliedEventIds.length) return undefined;
  return {
    reconstruction: "rdp5-signal-disassembly",
    appliedEventIds: [...disassembly.appliedEventIds],
    recombinant: componentReference(disassembly.mappings[recombinant]),
    majorParent: componentReference(disassembly.mappings[majorParent]),
    minorParent: componentReference(disassembly.mappings[minorParent]),
  };
}

function mapRawInterval(start, end, rotation, length) {
  if (rotation === 0) return { start, end, wraps: false };
  const mappedStart = (start + rotation) % length;
  const mappedEnd = (end + rotation) % length;
  if (mappedEnd === 0) return { start: mappedStart, end: length, wraps: false };
  return { start: mappedStart, end: mappedEnd, wraps: mappedStart > mappedEnd };
}

// After signal erasure, RDP5 breaks any newly detected signal spanning a
// deleted tract into continuously observed pieces and labels gap-adjacent
// breakpoints uncertain. Work in the scan's (possibly rotated) coordinates so
// origin-spanning candidates are split without inventing a linear cut.
export function splitCandidateAtStructuralGaps(candidate, disassembly, length) {
  if (!disassembly.componentCount || !disassembly.structuralMasks) return [candidate];
  const rotation = Number.isInteger(candidate.rotation) ? candidate.rotation : 0;
  const roles = [candidate.analysisRecombinant, candidate.analysisMajorParent, candidate.analysisMinorParent];
  const structuralAt = (rawSite) => {
    const site = (rawSite + rotation) % length;
    return roles.some((role) => disassembly.structuralMasks[role * length + site] !== 0);
  };
  const available = (rawSite) => !structuralAt(rawSite);
  const runs = [];
  let runStart = -1;
  let interrupted = false;
  for (let site = candidate.rawStart; site <= candidate.rawEnd; site += 1) {
    const observed = site < candidate.rawEnd && available(site);
    if (observed && runStart < 0) runStart = site;
    if (observed) continue;
    if (site < candidate.rawEnd) interrupted = true;
    if (runStart >= 0) {
      runs.push([runStart, site]);
      runStart = -1;
    }
  }
  const sourceEventIds = [...new Set(roles.flatMap((role) => {
    const mapping = disassembly.mappings[role];
    return [...(mapping.erasedEventIds ?? []), ...(mapping.sourceEventId ? [mapping.sourceEventId] : [])];
  }))];
  const vnpMargin = Math.max(1, Math.trunc(candidate.structuralUncertaintyVnps ?? 30));
  const variableAt = (rawSite) => {
    const site = (rawSite + rotation) % length;
    const values = roles.map((role) => disassembly.encoded[role * length + site]);
    return values.every((value) => value < 4) && (values[0] !== values[1] || values[0] !== values[2]);
  };
  const nearDeletedTract = (boundary, circular = false) => [-1, 1].some((direction) => {
    let variableSites = 0;
    for (let step = 0; step < length; step += 1) {
      let rawSite = direction < 0 ? boundary - step - 1 : boundary + step;
      if (!circular && (rawSite < 0 || rawSite >= length)) return false;
      rawSite = (rawSite % length + length) % length;
      if (structuralAt(rawSite)) return true;
      if (variableAt(rawSite)) variableSites += 1;
      if (variableSites > vnpMargin) return false;
    }
    return false;
  });
  const originalStartNear = nearDeletedTract(candidate.rawStart, candidate.circular === true);
  const originalEndNear = nearDeletedTract(candidate.rawEnd, candidate.circular === true);
  if (!interrupted && !originalStartNear && !originalEndNear) return [candidate];
  if (!interrupted) {
    return [{
      ...candidate,
      structuralUncertainty: {
        source: "rdp5-erased-signal-boundary",
        originalStart: candidate.start,
        originalEnd: candidate.end,
        originalWraps: candidate.wraps === true,
        piece: 1,
        pieces: 1,
        uncertainStart: originalStartNear,
        uncertainEnd: originalEndNear,
        adjacentEventIds: sourceEventIds,
      },
    }];
  }
  const retained = runs.filter(([start, end]) => end - start >= 4);
  return retained.map(([rawStart, rawEnd], index) => {
    const mapped = mapRawInterval(rawStart, rawEnd, rotation, length);
    return {
      ...candidate,
      ...mapped,
      rawStart,
      rawEnd,
      sourceRdp: undefined,
      methodSignals: (candidate.methodSignals ?? []).map((signal) => ({ ...signal, ...mapped })),
      structuralUncertainty: {
        source: "rdp5-erased-signal-boundary",
        originalStart: candidate.start,
        originalEnd: candidate.end,
        originalWraps: candidate.wraps === true,
        piece: index + 1,
        pieces: retained.length,
        uncertainStart: rawStart > candidate.rawStart || nearDeletedTract(rawStart, candidate.circular === true),
        uncertainEnd: rawEnd < candidate.rawEnd || nearDeletedTract(rawEnd, candidate.circular === true),
        adjacentEventIds: sourceEventIds,
      },
    };
  });
}
