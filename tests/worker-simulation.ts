import assert from "node:assert/strict";
import fs from "node:fs";
import { DEFAULT_OPTIONS, type AlignmentData } from "../app/rdp-core";

function randomGenerator(seed: number) {
  let state = seed >>> 0;
  return () => {
    state ^= state << 13;
    state ^= state >>> 17;
    state ^= state << 5;
    return (state >>> 0) / 4294967296;
  };
}

function mutate(sequence: string, rate: number, seed: number) {
  const random = randomGenerator(seed);
  const bases = "ACGT";
  return [...sequence].map((base) => {
    if (random() >= rate) return base;
    const alternatives = bases.replace(base, "");
    return alternatives[Math.floor(random() * alternatives.length)];
  }).join("");
}

const random = randomGenerator(0x1234abcd);
const ancestor = Array.from({ length: 1_200 }, () => "ACGT"[Math.floor(random() * 4)]).join("");
const nullAlignment: AlignmentData = {
  name: "Independent-mutation null",
  format: "generated",
  length: ancestor.length,
  createdAt: 1,
  sequences: Array.from({ length: 12 }, (_, index) => ({
    name: `Null-${index + 1}`,
    sequence: mutate(ancestor, 0.018, 100 + index),
    role: "both" as const,
  })),
};
const gappedAlignment: AlignmentData = {
  ...nullAlignment,
  name: "Gap-block null",
  sequences: nullAlignment.sequences.map((record, index) => ({
    ...record,
    sequence: index === 0 ? `${record.sequence.slice(0, 450)}${"-".repeat(120)}${record.sequence.slice(570)}` : record.sequence,
  })),
};

const wasm = fs.readFileSync(new URL("../public/wasm/rdp.wasm", import.meta.url));
const runtime = globalThis as typeof globalThis & {
  self: typeof globalThis & { onmessage: (event: { data: unknown }) => void };
  postMessage: (payload: unknown) => void;
};
runtime.self = runtime;
runtime.fetch = async () => new Response(wasm, { headers: { "content-type": "application/wasm" } });
await import("../public/rdp-worker.js");

function analyze(alignment: AlignmentData, jobId: number) {
  return new Promise<Record<string, unknown>>((resolve, reject) => {
    runtime.postMessage = (payload: unknown) => {
      const message = payload as Record<string, unknown>;
      if (message.jobId !== jobId) return;
      if (message.type === "result") resolve(message);
      if (message.type === "error") reject(new Error(String(message.message)));
    };
    runtime.self.onmessage({
      data: {
        type: "analyze",
        jobId,
        alignment,
        options: { ...DEFAULT_OPTIONS, candidateParents: 8, bootstrapReplicates: 25 },
      },
    });
  });
}

const nullResult = await analyze(nullAlignment, 21);
assert.equal((nullResult.events as unknown[]).length, 0, "independent mutation should not create a consensus event");
const gapResult = await analyze(gappedAlignment, 22);
assert.equal((gapResult.events as unknown[]).length, 0, "a gap block alone should not create a consensus event");
