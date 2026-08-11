import assert from "node:assert/strict";
import fs from "node:fs";
import { DEFAULT_OPTIONS, type AlignmentData } from "../app/rdp-core";

const length = 900;
const alpha = Array.from({ length }, (_, index) => "ACGT"[index % 4]).join("");
const beta = [...alpha].map((base, index) => index % 5 === 0 ? ({ A: "C", C: "G", G: "T", T: "A" }[base] ?? base) : base).join("");
const query = alpha.slice(0, 300) + beta.slice(300, 600) + alpha.slice(600);
const variant = (sequence: string, seed: number) => {
  const bases = [...sequence];
  for (let index = seed % 97; index < bases.length; index += 211) {
    bases[index] = ({ A: "C", C: "G", G: "T", T: "A" }[bases[index]] ?? bases[index]);
  }
  return bases.join("");
};
const alignment: AlignmentData = {
  name: "large-path regression",
  format: "generated",
  length,
  createdAt: 1,
  sequences: [
    { name: "Mosaic", sequence: query, role: "query" },
    ...Array.from({ length: 256 }, (_, index) => ({ name: `Alpha-${index}`, sequence: variant(alpha, index), role: "reference" as const })),
    ...Array.from({ length: 256 }, (_, index) => ({ name: `Beta-${index}`, sequence: variant(beta, index + 17), role: "reference" as const })),
  ],
};

const wasm = fs.readFileSync(new URL("../public/wasm/rdp.wasm", import.meta.url));
const runtime = globalThis as typeof globalThis & {
  self: typeof globalThis & { onmessage: (event: { data: unknown }) => void; location: { href: string } };
  postMessage: (payload: unknown) => void;
};
runtime.self = runtime;
runtime.location = { href: "https://example.test/rdp-web/rdp-worker.js" } as Location;
runtime.fetch = async (input) => {
  assert.equal(String(input), "https://example.test/rdp-web/wasm/rdp.wasm");
  return new Response(wasm, { headers: { "content-type": "application/wasm" } });
};

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
    jobId: 2,
    alignment,
    options: {
      ...DEFAULT_OPTIONS,
      mode: "query-reference",
      // This fixture exercises the explicitly approximate large-cohort path;
      // full-triplet enumeration has its own combinatorial regression.
      exhaustive: false,
      correction: "none",
      minMethods: 2,
    },
  },
});

const message = await result;
assert.equal(message.matrixCount, 24);
assert.match(String(message.matrixMode), /sampled\/stratified/);
assert.equal((message.distance as number[]).length, 24 ** 2);
assert.equal(message.comparisons, 28);
const events = message.events as Array<{ recombinant: number; start: number; end: number }>;
assert.ok(events.some((event) => event.recombinant === 0 && event.start < 360 && event.end > 540));
