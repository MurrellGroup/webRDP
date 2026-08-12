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
for (let site = 0; site < 480; site += 1) {
  const major = "ACGT"[site & 3];
  const minor = "ACGT"[(site + 1) & 3];
  records[0] += site >= 150 && site < 300 ? minor : major;
  records[1] += major;
  records[2] += minor;
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
    jobId: 73,
    alignment,
    options: {
      ...DEFAULT_OPTIONS,
      mode: "exploratory",
      methods: ["BootScan"],
      correction: "none",
      alpha: 0.05,
      minMethods: 1,
      bootstrapReplicates: 50,
      bootscanWindow: 80,
      bootscanStep: 10,
      bootscanCutoff: 0.7,
      cyclicDetection: false,
      polishBreakpoints: false,
      ancestralClustering: false,
      checkMisalignment: false,
    },
  },
});

const message = await result;
const events = message.events as Array<{
  recombinant: number;
  start: number;
  end: number;
  evidence: Array<{ method: string; pValue: number; calibration: string }>;
  methodSignals?: Array<{
    method: string;
    sourceBootscan?: {
      topology: number;
      baselineTopology: number;
      bootstrapReplicates: number;
      rawP: number;
    };
  }>;
}>;
const detected = events.find((event) => (
  event.recombinant === 0
  && event.start <= 150
  && event.end >= 300
  && event.methodSignals?.some((signal) => signal.method === "BootScan" && signal.sourceBootscan?.rawP < 0.05)
));
assert.ok(detected, "BootScan must independently seed the mosaic event without an RDP/GENECONV locator");
const evidence = detected.evidence.find((item) => item.method === "BootScan");
assert.ok(evidence && evidence.pValue < 0.05);
assert.match(evidence.calibration, /RDP5 RecScan distance batch/);
assert.equal(message.comparisons, 1);
assert.deepEqual(message.bootscanBatch, {
  calls: 1,
  triplets: 1,
  usedPairs: 3,
  windows: 50,
  replicates: 50,
  workspaceBytes: (message.bootscanBatch as { workspaceBytes: number }).workspaceBytes,
  relationshipMode: "distance",
});
