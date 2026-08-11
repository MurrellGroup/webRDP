import assert from "node:assert/strict";
import fs from "node:fs";
import { DEFAULT_OPTIONS, SOURCE_READY_METHODS, demoEvent, makeDemoAlignment } from "../app/rdp-core";

const wasm = fs.readFileSync(new URL("../public/wasm/rdp.wasm", import.meta.url));
const runtime = globalThis as typeof globalThis & {
  self: typeof globalThis & { onmessage: (event: { data: unknown }) => void };
  postMessage: (payload: unknown) => void;
};
runtime.self = runtime;
runtime.fetch = async () => new Response(wasm, { headers: { "content-type": "application/wasm" } });

const result = new Promise<Record<string, unknown>>((resolve, reject) => {
  runtime.postMessage = (payload: unknown) => {
    const message = payload as Record<string, unknown>;
    if (message.type === "recalculated") resolve(message);
    if (message.type === "error") reject(new Error(String(message.message)));
  };
});

await import("../public/rdp-worker.js");
runtime.self.onmessage({
  data: {
    type: "recalculate",
    jobId: 11,
    alignment: makeDemoAlignment(),
    options: { ...DEFAULT_OPTIONS, mode: "query-reference", correction: "holm" },
    event: { ...demoEvent(), start: 790, end: 1_530, evidenceStale: true },
    comparisons: 56,
  },
});

const message = await result;
const patch = message.patch as {
  evidence: Array<{ method: string; calibration: string; correctionScope: string }>;
  informativeSites: number;
  evidenceStale: boolean;
  diagnostics: { rateRatio: number };
  hypothesisTests: number;
  recalculationNote: string;
};
assert.deepEqual(
  patch.evidence.map((item) => item.method),
  SOURCE_READY_METHODS,
  "edited-event recalculation must expose only production source ports",
);
assert.equal(patch.evidenceStale, false);
assert.ok(patch.informativeSites > 250);
assert.ok(!patch.evidence.some((item) => ["GENECONV", "BootScan", "3Seq"].includes(item.method)));
assert.match(patch.evidence[0].correctionScope, /56 scanned triplets/);
assert.equal(patch.hypothesisTests, 56);
assert.match(patch.recalculationNote, /conservative Bonferroni/);
assert.ok(patch.diagnostics.rateRatio > 0);

const predecessor = { ...demoEvent(), id: "accepted-predecessor", decision: "accepted" as const, evidenceStale: false, start: 782, end: 1_538, wraps: false };
const derived = {
  ...demoEvent(),
  id: "derived-component-event",
  start: 900,
  end: 1_250,
  evidenceStale: true,
  componentProvenance: {
    reconstruction: "rdp5-signal-disassembly" as const,
    appliedEventIds: [predecessor.id],
    recombinant: { originIndex: predecessor.recombinant, kind: "extracted-tract" as const, lineage: [predecessor.id], sourceEventId: predecessor.id, parentLineage: [], start: predecessor.start, end: predecessor.end, wraps: false, erasedEventIds: [] },
    majorParent: { originIndex: predecessor.majorParent, kind: "remainder" as const, lineage: [], erasedEventIds: [] },
    minorParent: { originIndex: predecessor.minorParent, kind: "remainder" as const, lineage: [], erasedEventIds: [] },
  },
};
const componentResult = new Promise<Record<string, unknown>>((resolve, reject) => {
  runtime.postMessage = (payload: unknown) => {
    const componentMessage = payload as Record<string, unknown>;
    if (componentMessage.type === "recalculated") resolve(componentMessage);
    if (componentMessage.type === "error") reject(new Error(String(componentMessage.message)));
  };
});
runtime.self.onmessage({
  data: {
    type: "recalculate",
    jobId: 12,
    alignment: makeDemoAlignment(),
    options: { ...DEFAULT_OPTIONS, mode: "query-reference" },
    event: derived,
    disassemblyEvents: [predecessor],
    comparisons: 56,
  },
});
const componentMessage = await componentResult;
const componentPatch = componentMessage.patch as { componentProvenance: { recombinant: { kind: string; lineage: string[] } }; recalculationNote: string };
assert.equal(componentPatch.componentProvenance.recombinant.kind, "extracted-tract");
assert.deepEqual(componentPatch.componentProvenance.recombinant.lineage, [predecessor.id]);
assert.match(componentPatch.recalculationNote, /signal-disassembly lineage was rebuilt/);
