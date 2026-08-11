import assert from "node:assert/strict";
import test from "node:test";
import {
  buildDisassembledAlignment,
  candidateComponentProvenance,
  findComponentIndex,
  splitCandidateAtStructuralGaps,
} from "../public/rdp-disassembly.js";

function row(disassembly, index, length) {
  return disassembly.encoded.slice(index * length, (index + 1) * length);
}

test("accepted co-recombinant tracts are erased, extracted, and recursively disassembled", () => {
  const length = 20;
  const sequences = Array.from({ length: 4 }, (_, index) => ({ name: `s${index + 1}`, sequence: "A".repeat(length) }));
  const encoded = new Uint8Array(sequences.length * length);
  for (let sequence = 0; sequence < sequences.length; sequence += 1) {
    for (let site = 0; site < length; site += 1) encoded[sequence * length + site] = (sequence + site) % 4;
  }
  const first = {
    id: "event-a",
    decision: "accepted",
    evidenceStale: false,
    recombinant: 0,
    start: 5,
    end: 15,
    wraps: false,
    coRecombinantSets: [{ presumedRecombinant: 0, sequenceMembers: [0, 3] }],
  };
  const nested = {
    id: "event-b",
    decision: "accepted",
    evidenceStale: false,
    recombinant: 0,
    start: 8,
    end: 12,
    wraps: false,
    componentProvenance: { recombinant: { lineage: ["event-a"] } },
    coRecombinantSets: [{ presumedRecombinant: 0, sequenceMembers: [0, 3] }],
  };

  const disassembly = buildDisassembledAlignment(encoded, sequences, length, [first, nested]);
  assert.deepEqual(disassembly.appliedEventIds, ["event-a", "event-b"]);
  assert.equal(disassembly.componentCount, 4);
  assert.equal(disassembly.erasedCanonicalBases, 28);

  const remainder0 = findComponentIndex(disassembly, { originIndex: 0, kind: "remainder", lineage: [] });
  const componentA0 = findComponentIndex(disassembly, { originIndex: 0, kind: "extracted-tract", lineage: ["event-a"] });
  const componentB0 = findComponentIndex(disassembly, { originIndex: 0, kind: "extracted-tract", lineage: ["event-a", "event-b"] });
  const componentB3 = findComponentIndex(disassembly, { originIndex: 3, kind: "extracted-tract", lineage: ["event-a", "event-b"] });
  assert.ok([remainder0, componentA0, componentB0, componentB3].every((index) => index >= 0));
  assert.ok([...row(disassembly, remainder0, length).slice(5, 15)].every((base) => base === 4));
  assert.ok([...row(disassembly, componentA0, length).slice(8, 12)].every((base) => base === 4));
  assert.ok([...row(disassembly, componentA0, length).slice(5, 8)].every((base) => base < 4));
  assert.ok([...row(disassembly, componentB0, length).slice(8, 12)].every((base) => base < 4));
  assert.ok([...row(disassembly, componentB0, length).slice(0, 8)].every((base) => base === 4));
  assert.ok([...row({ encoded: disassembly.structuralMasks }, remainder0, length).slice(5, 15)].every((value) => value === 1));
  assert.ok([...row({ encoded: disassembly.structuralMasks }, componentA0, length).slice(5, 8)].every((value) => value === 0));
  assert.ok([...row({ encoded: disassembly.structuralMasks }, componentA0, length).slice(8, 12)].every((value) => value === 1));
  assert.ok([...row({ encoded: disassembly.structuralMasks }, componentB0, length).slice(8, 12)].every((value) => value === 0));

  const provenance = candidateComponentProvenance(disassembly, componentB0, 1, 2);
  assert.equal(provenance.recombinant.kind, "extracted-tract");
  assert.deepEqual(provenance.recombinant.lineage, ["event-a", "event-b"]);
  assert.equal(provenance.majorParent.originIndex, 1);
  assert.deepEqual(provenance.appliedEventIds, ["event-a", "event-b"]);
});

test("signals crossing erased tracts split into continuous pieces with uncertain adjacent breakpoints", () => {
  const length = 30;
  const sequences = Array.from({ length: 3 }, (_, index) => ({ name: `s${index + 1}`, sequence: "A".repeat(length) }));
  const encoded = new Uint8Array(sequences.length * length);
  const outer = {
    id: "outer",
    decision: "accepted",
    evidenceStale: false,
    recombinant: 0,
    start: 5,
    end: 25,
  };
  const nested = {
    id: "nested",
    decision: "accepted",
    evidenceStale: false,
    recombinant: 0,
    start: 12,
    end: 16,
    componentProvenance: { recombinant: { lineage: ["outer"] } },
  };
  const disassembly = buildDisassembledAlignment(encoded, sequences, length, [outer, nested]);
  const recombinant = findComponentIndex(disassembly, { originIndex: 0, kind: "extracted-tract", lineage: ["outer"] });
  assert.ok(recombinant >= 0);
  const pieces = splitCandidateAtStructuralGaps({
    recombinant: 0,
    majorParent: 1,
    minorParent: 2,
    analysisRecombinant: recombinant,
    analysisMajorParent: 1,
    analysisMinorParent: 2,
    rotation: 0,
    rawStart: 5,
    rawEnd: 25,
    start: 5,
    end: 25,
    wraps: false,
    sourceRdp: { pValue: 1e-12 },
    methodSignals: [{ method: "RDP", start: 5, end: 25, wraps: false }],
  }, disassembly, length);

  assert.deepEqual(pieces.map((piece) => [piece.start, piece.end]), [[5, 12], [16, 25]]);
  assert.equal(pieces[0].sourceRdp, undefined);
  assert.deepEqual(pieces.map((piece) => [piece.structuralUncertainty.uncertainStart, piece.structuralUncertainty.uncertainEnd]), [[true, true], [true, true]]);
  assert.deepEqual(pieces[0].structuralUncertainty.adjacentEventIds.sort(), ["nested", "outer"]);
  assert.deepEqual(pieces[0].methodSignals.map((signal) => [signal.start, signal.end]), [[5, 12]]);
});

test("RDP-window proximity to an erased tract marks an unsplit breakpoint uncertain", () => {
  const length = 40;
  const sequences = Array.from({ length: 3 }, (_, index) => ({ name: `s${index + 1}`, sequence: "A".repeat(length) }));
  const encoded = new Uint8Array(sequences.length * length);
  const disassembly = buildDisassembledAlignment(encoded, sequences, length, [{
    id: "prior",
    decision: "accepted",
    evidenceStale: false,
    recombinant: 0,
    start: 10,
    end: 20,
  }]);
  const candidate = {
    analysisRecombinant: 0,
    analysisMajorParent: 1,
    analysisMinorParent: 2,
    rotation: 0,
    rawStart: 4,
    rawEnd: 9,
    start: 4,
    end: 9,
    wraps: false,
    circular: false,
    structuralUncertaintyVnps: 30,
  };
  const [marked] = splitCandidateAtStructuralGaps(candidate, disassembly, length);
  assert.equal(marked.structuralUncertainty.pieces, 1);
  assert.equal(marked.structuralUncertainty.uncertainEnd, true);
  assert.ok(marked.structuralUncertainty.adjacentEventIds.includes("prior"));
});

test("stale and rejected events never alter the component alignment", () => {
  const length = 12;
  const sequences = [{ name: "one", sequence: "A".repeat(length) }];
  const encoded = new Uint8Array(length);
  const disassembly = buildDisassembledAlignment(encoded, sequences, length, [
    { id: "rejected", decision: "rejected", recombinant: 0, start: 2, end: 8 },
    { id: "stale", decision: "accepted", evidenceStale: true, recombinant: 0, start: 2, end: 8 },
  ]);
  assert.equal(disassembly.componentCount, 0);
  assert.deepEqual(disassembly.appliedEventIds, []);
  assert.deepEqual([...disassembly.encoded], [...encoded]);
});
