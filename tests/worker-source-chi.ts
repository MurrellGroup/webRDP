import assert from "node:assert/strict";
import fs from "node:fs";
import { DEFAULT_OPTIONS, type AlignmentData } from "../app/rdp-core";

const length = 600;
const bases = ["A", "C", "G", "T"];
let recombinant = "";
let major = "";
let minor = "";
for (let site = 0; site < length; site += 1) {
  const majorBase = bases[site % 4];
  const minorBase = bases[(site + 1) % 4];
  const inside = (site >= 100 && site < 180) || (site >= 350 && site < 450);
  recombinant += inside ? minorBase : majorBase;
  major += majorBase;
  minor += minorBase;
}
const alignment: AlignmentData = {
  name: "two source chi peak pairs",
  format: "generated",
  length,
  createdAt: 0,
  sequences: [
    { name: "mosaic", sequence: recombinant, role: "query" },
    { name: "major", sequence: major, role: "reference", referenceGroup: "major" },
    { name: "minor", sequence: minor, role: "reference", referenceGroup: "minor" },
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
    jobId: 93,
    alignment,
    options: {
      ...DEFAULT_OPTIONS,
      mode: "query-reference",
      methods: ["MaxChi", "Chimaera"],
      minMethods: 1,
      correction: "none",
      window: 40,
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
  methodSignals?: Array<{
    method: string;
    locator: string;
    sourceRoutine?: string;
    sourceChi?: {
      informativeSites: number;
      halfWindow: number;
      boundaryStatistics: [number, number];
      growthWidths: [number, number];
    };
  }>;
}>;
const calls = events.filter((event) => event.recombinant === 0);
assert.ok(calls.some((event) => event.start === 100 && event.end === 180), `first source peak pair must become a candidate: ${JSON.stringify(calls)}`);
assert.ok(calls.some((event) => event.start === 350 && event.end === 450), "second source peak pair must become a separate candidate");
for (const event of calls) {
  const signals = event.methodSignals?.filter((signal) => signal.method === "MaxChi" || signal.method === "Chimaera") ?? [];
  assert.ok(signals.length > 0, "retained chi-square events need method-specific source calls");
  assert.ok(signals.every((signal) => signal.sourceChi?.informativeSites === 600 && signal.sourceChi.halfWindow === 20));
  assert.ok(signals.every((signal) => /GrowMChiWin/.test(signal.sourceRoutine ?? "")), "source routine provenance must reach the event ledger");
  assert.ok(signals.every((signal) => /paired source peak basin/.test(signal.locator)), "the retired generic peak-pair locator must not be used");
}
assert.equal(message.chiSignalTruncations, 0);
