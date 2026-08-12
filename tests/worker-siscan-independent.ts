import assert from "node:assert/strict";
import fs from "node:fs";
import { DEFAULT_OPTIONS, type AlignmentData } from "../app/rdp-core";

const length = 240;
const records = ["", "", "", ""];
for (let site = 0; site < length; site += 1) {
  records[0] += site >= 80 && site < 160 ? "C" : "A";
  records[1] += "A";
  records[2] += "C";
  records[3] += "G";
}

const alignment: AlignmentData = {
  name: "role-agnostic standalone SiScan",
  format: "generated",
  length,
  createdAt: 0,
  // Deliberately designate no references. Exploratory RDP analysis must infer
  // the recombinant and both parents from the alignment itself.
  sequences: records.map((sequence, index) => ({
    name: ["mosaic", "lineage A", "lineage C", "fourth sequence"][index],
    sequence,
    role: "query" as const,
  })),
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
    jobId: 72,
    alignment,
    options: {
      ...DEFAULT_OPTIONS,
      mode: "exploratory",
      exhaustive: true,
      methods: ["SiScan"],
      minMethods: 1,
      correction: "none",
      window: 40,
      step: 10,
      siskanWindow: 40,
      siskanStep: 10,
      polishBreakpoints: false,
      ancestralClustering: false,
      cyclicalDetection: false,
      siskanScanPermutations: 100,
      siskanPValuePermutations: 500,
      randomSeed: 91,
    },
  },
});

const message = await result;
assert.equal(message.comparisons, 4, "four sequences contain exactly four unordered triplets");
assert.deepEqual(message.tripletKernelCalls, {
  rdp: 0,
  geneconv: 0,
  sourceChi: 0,
  threeSeq: 0,
  siscan: 4,
});
const events = message.events as Array<{
  recombinant: number;
  majorParent: number;
  minorParent: number;
  start: number;
  end: number;
  methodSignals?: Array<{
    method: string;
    sourceRoutine?: string;
    sourceSiScan?: { rawP: number; topologyTriplet: [number, number, number] };
  }>;
}>;
const event = events.find((candidate) => (
  candidate.recombinant === 0
  && candidate.majorParent === 1
  && candidate.minorParent === 2
  && Math.abs(candidate.start - 80) <= 2
  && Math.abs(candidate.end - 160) <= 2
));
assert.ok(event, `standalone SiScan must discover the known mosaic without reference designations: ${JSON.stringify(events)}`);
const signal = event.methodSignals?.find((entry) => entry.method === "SiScan");
assert.match(signal?.sourceRoutine ?? "", /GetSSOL.*DoPerms3P.*ShrinkRegionC/);
assert.deepEqual(signal?.sourceSiScan?.topologyTriplet, [0, 1, 2]);
assert.ok((signal?.sourceSiScan?.rawP ?? 1) < 1e-10);
for (const candidate of events) {
  for (const sourceSignal of candidate.methodSignals?.filter((entry) => entry.method === "SiScan") ?? []) {
    const source = sourceSignal.sourceSiScan as typeof sourceSignal.sourceSiScan & {
      recombinant: number;
      majorParent: number;
      minorParent: number;
    };
    assert.deepEqual(
      [source.recombinant, source.majorParent, source.minorParent],
      [candidate.recombinant, candidate.majorParent, candidate.minorParent],
      "SiScan evidence must never be borrowed from another topology orientation of the same triplet",
    );
  }
}
