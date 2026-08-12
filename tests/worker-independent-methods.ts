import assert from "node:assert/strict";
import fs from "node:fs";
import { DEFAULT_OPTIONS, parseAlignment } from "../app/rdp-core";

const wasm = fs.readFileSync(new URL("../public/wasm/rdp.wasm", import.meta.url));
const runtime = globalThis as typeof globalThis & {
  self: typeof globalThis & { onmessage: (event: { data: unknown }) => void };
  postMessage: (payload: unknown) => void;
};

runtime.self = runtime;
runtime.fetch = async () => new Response(wasm, {
  headers: { "content-type": "application/wasm" },
});

const records = ["", "", ""];
for (let site = 0; site < 600; site += 1) {
  const category = site >= 100 && site < 180 ? 0 : site % 2 ? 1 : 2;
  const bases = category === 0 ? ["A", "A", "C"] : category === 1 ? ["A", "C", "A"] : ["C", "A", "A"];
  for (let sequence = 0; sequence < 3; sequence += 1) records[sequence] += bases[sequence];
}
const alignment = parseAlignment(records.map((sequence, index) => `>sequence_${index + 1}\n${sequence}`).join("\n"));

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
      mode: "exploratory",
      methods: ["GENECONV"],
      correction: "none",
      alpha: 0.05,
      minMethods: 1,
      cyclicDetection: false,
      polishBreakpoints: false,
      ancestralClustering: false,
      checkMisalignment: false,
    },
  },
});

const message = await result;
const events = message.events as Array<{
  start: number;
  end: number;
  methodSignals?: Array<{ method: string; sourceGeneconv?: { track: number; informativeSites: number; rawP: number } }>;
}>;
assert.ok(events.some((event) => (
  event.start === 100
  && event.end === 180
  && event.methodSignals?.some((signal) => (
    signal.method === "GENECONV"
    && signal.sourceGeneconv?.track === 0
    && signal.sourceGeneconv.informativeSites === 600
    && signal.sourceGeneconv.rawP < 0.05
  ))
)), "GENECONV must independently discover the source fragment without an RDP seed");
assert.equal(message.comparisons, 1);
assert.deepEqual(message.tripletKernelCalls, { rdp: 0, geneconv: 1, sourceChi: 0, threeSeq: 0, siscan: 0 });
