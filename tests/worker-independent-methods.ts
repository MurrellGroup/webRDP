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
const events = message.events as Array<unknown>;
assert.equal(events.length, 0, "a method awaiting its complete author-source port must not silently run a simplified discovery kernel");
assert.deepEqual(message.tripletKernelCalls, { rdp: 0, sourceChi: 0 });
