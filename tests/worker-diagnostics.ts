import assert from "node:assert/strict";
import fs from "node:fs";
import { DEFAULT_OPTIONS } from "../app/rdp-core";

const length = 600;
const major = "A".repeat(length);
const minor = Array.from({ length }, (_, site) => (
  (site >= 200 && site < 400) || site % 20 === 0 ? "C" : "A"
)).join("");
const recombinant = major.slice(0, 200) + minor.slice(200, 400) + major.slice(400);
const alignment = {
  name: "Rate-variation challenge",
  format: "generated" as const,
  length,
  createdAt: 1,
  sequences: [
    { name: "Mosaic", sequence: recombinant, role: "query" as const },
    { name: "Major", sequence: major, role: "reference" as const },
    { name: "Minor", sequence: minor, role: "reference" as const },
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
    jobId: 7,
    alignment,
    options: {
      ...DEFAULT_OPTIONS,
      mode: "query-reference",
      correction: "none",
      minMethods: 1,
      candidateParents: 3,
      window: 60,
      bootstrapReplicates: 25,
    },
  },
});

const message = await result;
const events = message.events as Array<{
  recombinant: number;
  diagnostics: { rateRatio: number };
  warnings: string[];
}>;
const event = events.find((candidate) => candidate.recombinant === 0);
assert.ok(event, "the challenge should retain its strong mosaic-like signal");
assert.ok(event.diagnostics.rateRatio > 4, `expected a strong rate shift, got ${event.diagnostics.rateRatio}`);
assert.ok(event.warnings.some((warning) => warning.includes("variable-site density")));

const diagnostics = message.diagnostics as {
  sampledSequences: number;
  testedSitePairs: number;
  fourGameteFraction: number;
  proximityPermutationP: number;
  proximityPermutationReplicates: number;
};
assert.equal(diagnostics.sampledSequences, 3);
assert.ok(diagnostics.testedSitePairs > 0);
assert.ok(diagnostics.fourGameteFraction >= 0 && diagnostics.fourGameteFraction <= 1);
assert.ok(diagnostics.proximityPermutationP > 0 && diagnostics.proximityPermutationP <= 1);
assert.equal(diagnostics.proximityPermutationReplicates, 199);
