import assert from "node:assert/strict";
import fs from "node:fs";
import { DEFAULT_OPTIONS, demoEvent, makeDemoAlignment } from "../app/rdp-core";

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
    options: { ...DEFAULT_OPTIONS, mode: "query-reference", correction: "none" },
    event: { ...demoEvent(), start: 790, end: 1_530, evidenceStale: true },
  },
});

const message = await result;
const patch = message.patch as {
  evidence: Array<{ method: string; calibration: string }>;
  informativeSites: number;
  evidenceStale: boolean;
  diagnostics: { rateRatio: number };
};
assert.equal(patch.evidence.length, 7);
assert.equal(patch.evidenceStale, false);
assert.ok(patch.informativeSites > 250);
assert.match(patch.evidence.find((item) => item.method === "3Seq")?.calibration ?? "", /HGRW/);
assert.ok(patch.diagnostics.rateRatio > 0);
