import assert from "node:assert/strict";
import fs from "node:fs";
import { DEFAULT_OPTIONS, makeDemoAlignment } from "../app/rdp-core";

const wasm = fs.readFileSync(new URL("../public/wasm/rdp.wasm", import.meta.url));
const runtime = globalThis as typeof globalThis & {
  self: typeof globalThis & { onmessage: (event: { data: unknown }) => void };
  postMessage: (payload: unknown) => void;
};

runtime.self = runtime;
runtime.fetch = async () => new Response(wasm, {
  headers: { "content-type": "application/wasm" },
});
let partialCheckpoints = 0;
const alignment = {
  ...makeDemoAlignment(),
  sequences: makeDemoAlignment().sequences.map((record) => ({
    ...record,
    referenceGroup: record.role === "reference" ? record.name.split("-")[0] : undefined,
  })),
};

const result = new Promise<Record<string, unknown>>((resolve, reject) => {
  runtime.postMessage = (payload: unknown) => {
    const message = payload as Record<string, unknown>;
    if (message.type === "partial") partialCheckpoints += 1;
    if (message.type === "result") resolve(message);
    if (message.type === "error") reject(new Error(String(message.message)));
  };
});

await import("../public/rdp-worker.js");
runtime.self.onmessage({
  data: {
    type: "analyze",
    jobId: 1,
    alignment,
    options: { ...DEFAULT_OPTIONS, mode: "query-reference" },
  },
});

const message = await result;
const events = message.events as Array<{
  id: string;
  recombinant: number;
  start: number;
  end: number;
  majorParent: number;
  minorParent: number;
  coRecombinantSets?: Array<{ presumedRecombinant: number; sequenceMembers: number[] }>;
  recombinantIdentification?: { recommended: number; confidence: number; tests: Array<{ sourceRoutine: string }> };
  methodSignals?: Array<{ method: string; sourceRoutine?: string; outgroup?: number | null; permutations?: number; profile?: Array<{ position: number; z: number; topology: number }> }>;
  evidence?: Array<{ method: string; calibration: string; statisticLabel: string }>;
}>;
assert.ok(events.length > 0, "the positive-control alignment should yield events");
assert.ok(partialCheckpoints > 0, "the worker should checkpoint recoverable partial candidates during discovery");
const mosaic = events.find((event) => event.recombinant === 0);
assert.ok(mosaic, "Mosaic-X should be identified as a recombinant candidate");
assert.ok(Math.abs(mosaic.start - 782) <= 80, `left breakpoint ${mosaic.start} should localize near 782`);
assert.ok(Math.abs(mosaic.end - 1538) <= 80, `right breakpoint ${mosaic.end} should localize near 1538`);
assert.ok([2, 3, 4].includes(mosaic.majorParent), "an Alpha sequence should be the major parent");
assert.ok([5, 6, 7, 11].includes(mosaic.minorParent), "a Beta sequence should be the minor parent");
assert.notEqual(alignment.sequences[mosaic.majorParent].referenceGroup, alignment.sequences[mosaic.minorParent].referenceGroup, "grouped query/reference scans must draw parents from distinct reference groups");
assert.ok(mosaic.recombinantIdentification, "a source recombinant-identification ledger must be persisted for every characterized event");
assert.equal(mosaic.recombinantIdentification?.recommended, 0, "the source profile consensus should recover the known mosaic polarity");
assert.ok(mosaic.recombinantIdentification?.tests.some((test) => /MakePhPrScore/.test(test.sourceRoutine)));
const sisterSignal = mosaic.methodSignals?.find((signal) => signal.method === "SiScan");
assert.ok(sisterSignal, "the known mosaic should retain a source-confirmed SiScan topology run");
assert.match(sisterSignal?.sourceRoutine ?? "", /GetSSOL.*DoPerms3P.*ShrinkRegionC/);
assert.equal(sisterSignal?.permutations, DEFAULT_OPTIONS.siskanPValuePermutations);
assert.ok((sisterSignal?.profile?.length ?? 0) > 10, "the source topology-window trace must remain available to the interactive method plot");
assert.match(mosaic.evidence?.find((row) => row.method === "SiScan")?.calibration ?? "", /vertical permutation Z/);
const parentFilteredResult = new Promise<Record<string, unknown>>((resolve, reject) => {
  runtime.postMessage = (payload: unknown) => {
    const filteredMessage = payload as Record<string, unknown>;
    if (filteredMessage.type === "result") resolve(filteredMessage);
    if (filteredMessage.type === "error") reject(new Error(String(filteredMessage.message)));
  };
});
runtime.self.onmessage({
  data: {
    type: "analyze",
    jobId: 2,
    alignment,
    options: { ...DEFAULT_OPTIONS, mode: "query-reference" },
    excludedParents: [2, 3, 4],
  },
});
const parentFilteredEvents = (await parentFilteredResult).events as Array<{ majorParent: number; minorParent: number }>;
assert.ok(parentFilteredEvents.every((event) => ![2, 3, 4].includes(event.majorParent) && ![2, 3, 4].includes(event.minorParent)), "accepted mosaic parent proxies must be absent from an event-aware rescan parent pool");

const accepted = { ...mosaic, decision: "accepted", evidenceStale: false };
const acceptedMembers = accepted.coRecombinantSets?.find((set) => set.presumedRecombinant === accepted.recombinant)?.sequenceMembers ?? [accepted.recombinant];
const disassembledResult = new Promise<Record<string, unknown>>((resolve, reject) => {
  runtime.postMessage = (payload: unknown) => {
    const disassembledMessage = payload as Record<string, unknown>;
    if (disassembledMessage.type === "result") resolve(disassembledMessage);
    if (disassembledMessage.type === "error") reject(new Error(String(disassembledMessage.message)));
  };
});
runtime.self.onmessage({
  data: {
    type: "analyze",
    jobId: 3,
    alignment,
    options: { ...DEFAULT_OPTIONS, mode: "query-reference" },
    disassemblyEvents: [accepted],
    excludedTargets: acceptedMembers,
    excludedParents: acceptedMembers,
  },
});
const disassembledMessage = await disassembledResult;
const disassembly = disassembledMessage.disassembly as { appliedEvents: number; components: number; erasedCanonicalBases: number };
const disassembledEvents = disassembledMessage.events as Array<{ recombinant: number; majorParent: number; minorParent: number; componentProvenance?: { appliedEventIds: string[] } }>;
assert.equal(disassembly.appliedEvents, 1);
assert.ok(disassembly.components >= 1);
assert.ok(disassembly.erasedCanonicalBases > 0);
assert.ok(disassembledEvents.every((event) => [event.recombinant, event.majorParent, event.minorParent].every((index) => index < alignment.sequences.length)), "private component indexes must never leak into user-facing event roles");
assert.ok(disassembledEvents.every((event) => event.componentProvenance?.appliedEventIds.includes(accepted.id)), "rescanned events must retain their signal-disassembly provenance");
console.log(JSON.stringify({
  events: events.length,
  mosaic: {
    recombinant: mosaic.recombinant,
    majorParent: mosaic.majorParent,
    minorParent: mosaic.minorParent,
    start: mosaic.start,
    end: mosaic.end,
  },
  elapsedMs: message.elapsedMs,
}, null, 2));
