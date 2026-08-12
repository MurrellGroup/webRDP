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

const records = ["", "", "", ""];
for (let site = 0; site < 480; site += 1) {
  const major = "ACGT"[site & 3];
  const minor = "ACGT"[(site + 1) & 3];
  records[0] += site >= 150 && site < 300 ? minor : major;
  records[1] += major;
  records[2] += minor;
  records[3] += "ACGT"[(site + 2) & 3];
}
const alignment = parseAlignment(records.map((sequence, index) => `>sequence_${index + 1}\n${sequence}`).join("\n"));

await import("../public/rdp-worker.js");
const analyze = (jobId: number, bootscanRelationshipMode: "distance" | "upgma" | "neighbor-joining") => (
  new Promise<Record<string, unknown>>((resolve, reject) => {
    runtime.postMessage = (payload: unknown) => {
      const message = payload as Record<string, unknown>;
      if (message.type === "result") resolve(message);
      if (message.type === "error") reject(new Error(String(message.message)));
    };
    runtime.self.onmessage({
      data: {
        type: "analyze",
        jobId,
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
          bootscanRelationshipMode,
          cyclicDetection: false,
          polishBreakpoints: false,
          ancestralClustering: false,
          checkMisalignment: false,
        },
      },
    });
  })
);

type BootscanEvent = {
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
      relationshipMode: "distance" | "upgma" | "neighbor-joining";
    };
  }>;
};

for (const [index, relationshipMode] of (["distance", "upgma", "neighbor-joining"] as const).entries()) {
  const message = await analyze(73 + index, relationshipMode);
  const events = message.events as BootscanEvent[];
  const detected = events.find((event) => (
    event.recombinant === 0
    && event.start <= 150
    && event.end >= 300
    && event.methodSignals?.some((signal) => (
      signal.method === "BootScan"
      && signal.sourceBootscan?.rawP < 0.05
      && signal.sourceBootscan.relationshipMode === relationshipMode
    ))
  ));
  assert.ok(detected, `${relationshipMode} BootScan must independently seed the mosaic event`);
  const evidence = detected.evidence.find((item) => item.method === "BootScan");
  assert.ok(evidence && evidence.pValue < 0.05);
  assert.match(evidence.calibration, new RegExp(`RDP5 RecScan ${relationshipMode} batch`));
  assert.equal(message.comparisons, 4);
  const batch = message.bootscanBatch as {
    calls: number;
    triplets: number;
    usedPairs: number;
    windows: number;
    replicates: number;
    workspaceBytes: number;
    relationshipMode: string;
  };
  assert.equal(batch.calls, 1);
  assert.equal(batch.triplets, 4);
  assert.equal(batch.usedPairs, 6);
  assert.equal(batch.windows, 50);
  assert.equal(batch.replicates, 50);
  assert.ok(batch.workspaceBytes > 0);
  assert.equal(batch.relationshipMode, relationshipMode);
}
