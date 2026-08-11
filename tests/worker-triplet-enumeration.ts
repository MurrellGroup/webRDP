import assert from "node:assert/strict";
import fs from "node:fs";
import { DEFAULT_OPTIONS, parseAlignment } from "../app/rdp-core";

assert.equal(DEFAULT_OPTIONS.exhaustive, true, "new analyses must default to full RDP5 triplet enumeration");

const wasm = fs.readFileSync(new URL("../public/wasm/rdp.wasm", import.meta.url));
const runtime = globalThis as typeof globalThis & {
  self: typeof globalThis & { onmessage: (event: { data: unknown }) => void };
  postMessage: (payload: unknown) => void;
};
runtime.self = runtime;
runtime.fetch = async () => new Response(wasm, { headers: { "content-type": "application/wasm" } });

const bases = ["A", "C", "G", "T"];
const fasta = Array.from({ length: 5 }, (_, sequence) => {
  const value = Array.from({ length: 180 }, (_, site) => bases[(site * 7 + sequence * 11 + Math.floor(site / 37)) % 4]).join("");
  return `>sequence_${sequence + 1}\n${value}`;
}).join("\n");
const alignment = parseAlignment(fasta);

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
    jobId: 711,
    alignment,
    options: {
      ...DEFAULT_OPTIONS,
      methods: ["RDP", "MaxChi", "Chimaera"],
      correction: "none",
      alpha: 1,
      minMethods: 1,
      polishBreakpoints: false,
      ancestralClustering: false,
      checkMisalignment: false,
    },
  },
});

const message = await result;
assert.equal(message.comparisons, 10, "five sequences must produce C(5,3)=10 unordered concrete triplets");
assert.equal(message.tripletMode, "all-concrete-triplets");
assert.equal(message.concreteTripletInputs, true);
assert.match(String(message.engine), /all concrete sequence triplets/);
assert.deepEqual(message.tripletKernelCalls, { rdp: 10, sourceChi: 10 }, "source kernels must run once—not three times—for each unordered triplet");
