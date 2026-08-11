import assert from "node:assert/strict";
import fs from "node:fs";
import { DEFAULT_OPTIONS, makeDemoAlignment } from "../app/rdp-core";

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
    if (message.type === "result") resolve(message);
    if (message.type === "error") reject(new Error(String(message.message)));
  };
});

await import("../public/rdp-worker.js");
runtime.self.onmessage({
  data: {
    type: "analyze",
    jobId: 881,
    alignment: makeDemoAlignment(),
    options: {
      ...DEFAULT_OPTIONS,
      mode: "query-reference",
      methods: ["RDP", "MaxChi", "Chimaera"],
      minMethods: 1,
      correction: "none",
      polishBreakpoints: false,
      ancestralClustering: false,
      checkMisalignment: false,
      maximumDetectionCycles: 8,
    },
    cyclicDetection: true,
  },
});

const message = await result;
const cycle = message.detectionCycle as {
  eventsApplied: number;
  passes: number;
  initialComparisons: number;
  redoComparisons: number;
  stoppedBecause: string;
};
const events = message.events as Array<{
  decision: string;
  note: string;
  history: Array<{ action: string }>;
  componentProvenance?: { appliedEventIds: string[] };
}>;

assert.ok(events.length >= 1, "the positive control must enter the cyclic event ledger");
assert.equal(cycle.eventsApplied, events.length);
assert.ok(cycle.passes >= 2, "the strongest event must be erased/extracted before at least one redo pass");
assert.ok(cycle.initialComparisons > 0);
assert.ok(cycle.redoComparisons > 0, "changed-origin triplets must be put on the redo list");
assert.ok(events.every((event) => event.decision === "unreviewed"), "automatic signal disassembly must not masquerade as analyst acceptance");
assert.ok(events.every((event) => /erase\/extract cycle/.test(event.note)));
assert.ok(events.every((event) => event.history.some((entry) => /Cyclic detection round/.test(entry.action))));
for (let index = 1; index < events.length; index += 1) {
  assert.ok(
    (events[index].componentProvenance?.appliedEventIds.length ?? 0) >= 1,
    "later detections must retain their component-lineage provenance",
  );
}
