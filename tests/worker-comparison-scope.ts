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

const major = Array.from({ length: 360 }, () => "A");
const minor = Array.from({ length: 360 }, () => "G");
const recombinant = Array.from({ length: 360 }, (_, site) => site >= 120 && site < 240 ? "G" : "A");
// RDP5's triplet screen deliberately requires all three pair-match classes.
// Sparse parent-parent sites keep this fixture inside that real source gate.
for (let site = 7; site < 360; site += 17) {
  major[site] = "C";
  minor[site] = "C";
  recombinant[site] = "T";
}
const alignment = parseAlignment(`>major\n${major.join("")}\n>minor\n${minor.join("")}\n>recombinant\n${recombinant.join("")}\n`);

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
    jobId: 81,
    alignment,
    options: {
      ...DEFAULT_OPTIONS,
      mode: "exploratory",
      exhaustive: true,
      methods: ["RDP"],
      correction: "none",
      alpha: 1,
      minMethods: 1,
      polishBreakpoints: false,
      ancestralClustering: false,
    },
  },
});

const message = await result;
assert.equal(
  message.comparisons,
  1,
  "three exploratory target passes are one unordered triplet, not three multiplicity tests",
);
assert.ok((message.events as unknown[]).length > 0, "the scope fixture must exercise an actual retained event");
