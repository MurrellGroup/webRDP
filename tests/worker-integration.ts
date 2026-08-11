import assert from "node:assert/strict";
import fs from "node:fs";
import { DEFAULT_OPTIONS, makeDemoAlignment } from "../app/rdp-core";

const wasm = fs.readFileSync(new URL("../public/wasm/rdp.wasm", import.meta.url));
const runtime = globalThis as typeof globalThis & {
  self: typeof globalThis & { onmessage: (event: { data: unknown }) => void };
  postMessage: (payload: unknown) => void;
};

runtime.self = runtime;
runtime.fetch = async () => new Response(wasm, {
  headers: { "content-type": "application/wasm" },
});
let partialCheckpoints = 0;

const result = new Promise<Record<string, unknown>>((resolve, reject) => {
  runtime.postMessage = (payload: unknown) => {
    const message = payload as Record<string, unknown>;
    if (message.type === "partial") partialCheckpoints += 1;
    if (message.type === "result") resolve(message);
    if (message.type === "error") reject(new Error(String(message.message)));
  };
});

await import("../public/rdp-worker.js");
runtime.self.onmessage({
  data: {
    type: "analyze",
    jobId: 1,
    alignment: makeDemoAlignment(),
    options: { ...DEFAULT_OPTIONS, mode: "query-reference" },
  },
});

const message = await result;
const events = message.events as Array<{
  recombinant: number;
  start: number;
  end: number;
  majorParent: number;
  minorParent: number;
}>;
assert.ok(events.length > 0, "the positive-control alignment should yield events");
assert.ok(partialCheckpoints > 0, "the worker should checkpoint recoverable partial candidates during discovery");
const mosaic = events.find((event) => event.recombinant === 0);
assert.ok(mosaic, "Mosaic-X should be identified as a recombinant candidate");
assert.ok(Math.abs(mosaic.start - 782) <= 80, `left breakpoint ${mosaic.start} should localize near 782`);
assert.ok(Math.abs(mosaic.end - 1538) <= 80, `right breakpoint ${mosaic.end} should localize near 1538`);
assert.ok([2, 3, 4].includes(mosaic.majorParent), "an Alpha sequence should be the major parent");
assert.ok([5, 6, 7, 11].includes(mosaic.minorParent), "a Beta sequence should be the minor parent");
const parentFilteredResult = new Promise<Record<string, unknown>>((resolve, reject) => {
  runtime.postMessage = (payload: unknown) => {
    const filteredMessage = payload as Record<string, unknown>;
    if (filteredMessage.type === "result") resolve(filteredMessage);
    if (filteredMessage.type === "error") reject(new Error(String(filteredMessage.message)));
  };
});
runtime.self.onmessage({
  data: {
    type: "analyze",
    jobId: 2,
    alignment: makeDemoAlignment(),
    options: { ...DEFAULT_OPTIONS, mode: "query-reference" },
    excludedParents: [2, 3, 4],
  },
});
const parentFilteredEvents = (await parentFilteredResult).events as Array<{ majorParent: number; minorParent: number }>;
assert.ok(parentFilteredEvents.every((event) => ![2, 3, 4].includes(event.majorParent) && ![2, 3, 4].includes(event.minorParent)), "accepted mosaic parent proxies must be absent from an event-aware rescan parent pool");
console.log(JSON.stringify({ events: events.length, mosaic, elapsedMs: message.elapsedMs }, null, 2));
