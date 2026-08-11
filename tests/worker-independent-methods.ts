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

const source = makeDemoAlignment();
const alignment = {
  ...source,
  sequences: source.sequences.map((record) => ({
    ...record,
    referenceGroup: record.role === "reference" ? record.name.split("-")[0] : undefined,
  })),
};

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
    jobId: 71,
    alignment,
    options: {
      ...DEFAULT_OPTIONS,
      mode: "query-reference",
      methods: ["GENECONV"],
      correction: "none",
      alpha: 1,
      minMethods: 1,
      polishBreakpoints: false,
      ancestralClustering: false,
    },
  },
});

const message = await result;
const events = message.events as Array<{
  recombinant: number;
  start: number;
  end: number;
  methodSignals: Array<{ method: string; locator: string }>;
  evidence: Array<{ method: string; calibration: string; supported: boolean }>;
}>;
assert.ok(events.length > 0, "GENECONV must discover candidates when RDP is disabled");
assert.ok(events.every((event) => event.methodSignals.every((signal) => signal.method !== "shared-screen")), "the retired shared CUSUM seed must not reappear");
assert.ok(events.every((event) => event.methodSignals.some((signal) => signal.method === "GENECONV")), "retained GENECONV candidates must carry their own locator provenance");
assert.ok(events.every((event) => event.evidence.length === 1 && event.evidence[0].method === "GENECONV" && event.evidence[0].supported), "single-method evidence must be calibrated from the co-located GENECONV signal");
const mosaic = events.find((event) => event.recombinant === 0);
assert.ok(mosaic, "the independent GENECONV scan must recover the Mosaic-X positive control");
assert.ok(mosaic.end > mosaic.start, "the independently discovered fragment must define a tract");

