import assert from "node:assert/strict";
import fs from "node:fs";
import { DEFAULT_OPTIONS } from "../app/rdp-core";

const length = 600;
const major = "A".repeat(length);
const minor = "C".repeat(length);
const recombinant = minor.slice(0, 100) + major.slice(100, 520) + minor.slice(520);
const alignment = {
  name: "Circular origin-spanning positive control",
  format: "generated" as const,
  length,
  createdAt: 1,
  sequences: [
    { name: "Circular-mosaic", sequence: recombinant, role: "query" as const },
    { name: "Major-parent", sequence: major, role: "reference" as const },
    { name: "Minor-parent", sequence: minor, role: "reference" as const },
  ],
};

const wasm = fs.readFileSync(new URL("../public/wasm/rdp.wasm", import.meta.url));
const runtime = globalThis as typeof globalThis & {
  self: typeof globalThis & { onmessage: (event: { data: unknown }) => void };
  postMessage: (payload: unknown) => void;
};

runtime.self = runtime;
runtime.fetch = async () => new Response(wasm, {
  headers: { "content-type": "application/wasm" },
});

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
    jobId: 3,
    alignment,
    options: {
      ...DEFAULT_OPTIONS,
      mode: "query-reference",
      circular: true,
      correction: "none",
      minMethods: 1,
      window: 60,
      candidateParents: 3,
    },
  },
});

const message = await result;
const events = message.events as Array<{
  recombinant: number;
  start: number;
  end: number;
  wraps: boolean;
  majorParent: number;
  minorParent: number;
  breakpointModel?: { method: string; stateSwitches?: number };
}>;
const event = events.find((candidate) => candidate.recombinant === 0 && candidate.wraps);
assert.ok(event, "the origin-spanning tract should be retained as a wrapping event");
assert.ok(Math.abs(event.start - 520) <= 2, `circular start ${event.start} should recover 520`);
assert.ok(Math.abs(event.end - 100) <= 2, `circular end ${event.end} should recover 100`);
assert.equal(event.majorParent, 1);
assert.equal(event.minorParent, 2);
assert.equal(event.breakpointModel?.method, "two-state-hmm");
assert.equal(event.breakpointModel?.stateSwitches, 2);
assert.equal(message.comparisons, 2, "circular mode should correct for both tested origins");
assert.match(String(message.engine), /dual-origin circular scan/);
