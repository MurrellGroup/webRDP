import assert from "node:assert/strict";
import fs from "node:fs";
import { DEFAULT_OPTIONS, type AlignmentData } from "../app/rdp-core";

const length = 600;
const bases = ["A", "C", "G", "T"];
const sequences = ["", "", "", ""];
for (let site = 0; site < length; site += 1) {
  const major = bases[site % 4];
  const minor = bases[(site + 1) % 4];
  const outgroup = bases[(site + 2) % 4];
  const inside = (site >= 100 && site < 180) || (site >= 350 && site < 450);
  sequences[0] += inside ? minor : major;
  sequences[1] += major;
  sequences[2] += minor;
  sequences[3] += outgroup;
}
const alignment: AlignmentData = {
  name: "two SiScan topology runs",
  format: "generated",
  length,
  createdAt: 0,
  sequences: [
    { name: "mosaic", sequence: sequences[0], role: "query" },
    { name: "major", sequence: sequences[1], role: "reference", referenceGroup: "major" },
    { name: "minor", sequence: sequences[2], role: "reference", referenceGroup: "minor" },
    { name: "outgroup", sequence: sequences[3], role: "reference", referenceGroup: "outgroup" },
  ],
};

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
    jobId: 71,
    alignment,
    options: {
      ...DEFAULT_OPTIONS,
      mode: "query-reference",
      // Source SiScan currently confirms/expands candidates from a source
      // triplet detector; its independent desktop batch orchestration remains
      // a disclosed parity task.
      methods: ["RDP", "MaxChi", "Chimaera", "SiScan"],
      minMethods: 1,
      candidateParents: 3,
      window: 40,
      step: 10,
      polishBreakpoints: false,
      ancestralClustering: false,
      siskanScanPermutations: 100,
      siskanPValuePermutations: 500,
      randomSeed: 91,
    },
  },
});

const message = await result;
const events = message.events as Array<{
  recombinant: number;
  majorParent: number;
  minorParent: number;
  start: number;
  end: number;
  methodSignals?: Array<{ method: string; sourceRoutine?: string }>;
}>;
const calls = events.filter((event) => (
  event.recombinant === 0
  && new Set([event.majorParent, event.minorParent]).has(1)
  && new Set([event.majorParent, event.minorParent]).has(2)
  && event.methodSignals?.some((signal) => signal.method === "SiScan" && /DoPerms3P/.test(signal.sourceRoutine ?? ""))
));
assert.ok(calls.some((event) => Math.abs(event.start - 100) <= 2 && Math.abs(event.end - 180) <= 2), `first source topology run must become an event: calls=${JSON.stringify(calls)} all=${JSON.stringify(events)}`);
assert.ok(calls.some((event) => Math.abs(event.start - 350) <= 2 && Math.abs(event.end - 450) <= 2), "second source topology run must become a separate event");
