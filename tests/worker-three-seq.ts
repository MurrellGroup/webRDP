import assert from "node:assert/strict";
import fs from "node:fs";
import { DEFAULT_OPTIONS, parseAlignment } from "../app/rdp-core";

const wasm = fs.readFileSync(new URL("../public/wasm/rdp.wasm", import.meta.url));
const runtime = globalThis as typeof globalThis & {
  self: typeof globalThis & { onmessage: (event: { data: unknown }) => void };
  postMessage: (payload: unknown) => void;
};
runtime.self = runtime;
runtime.fetch = async () => new Response(wasm, { headers: { "content-type": "application/wasm" } });

const records = ["", "", ""];
const bases = ["A", "C", "G", "T"];
for (let site = 0; site < 240; site += 1) {
  const major = bases[site % 4];
  const minor = bases[(site + 1) % 4];
  records[0] += site >= 80 && site < 160 ? minor : major;
  records[1] += major;
  records[2] += minor;
}
const alignment = parseAlignment(records.map((sequence, index) => `>sequence_${index + 1}\n${sequence}`).join("\n"));

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
    jobId: 712,
    alignment,
    options: {
      ...DEFAULT_OPTIONS,
      methods: ["3Seq"],
      correction: "none",
      alpha: 0.05,
      minMethods: 1,
      cyclicDetection: false,
      polishBreakpoints: false,
      ancestralClustering: false,
      checkMisalignment: false,
    },
  },
});

const message = await result;
const events = message.events as Array<{
  recombinant: number;
  start: number;
  end: number;
  wraps?: boolean;
  evidence: Array<{ method: string; pValue: number; calibration: string }>;
  methodSignals?: Array<{
    method: string;
    sourceThreeSeq?: {
      upSteps: number;
      downSteps: number;
      descent: number;
      informativeSites: number;
      rawP: number;
      sourceWrap?: boolean;
      linearComplement?: boolean;
      splitRefined?: boolean;
      fullDescent?: number;
      splitInformativeSites?: number;
    };
  }>;
  structuralUncertainty?: { piece: number; pieces: number; adjacentEventIds: string[] };
}>;
const event = events.find((candidate) => (
  candidate.recombinant === 0
  && candidate.start === 80
  && candidate.end === 160
  && candidate.methodSignals?.some((signal) => signal.method === "3Seq")
));
assert.ok(event, "3Seq must independently discover the mosaic without an RDP or consensus seed");
const signal = event.methodSignals?.find((entry) => entry.method === "3Seq")?.sourceThreeSeq;
assert.deepEqual(
  signal && [signal.upSteps, signal.downSteps, signal.descent, signal.informativeSites],
  [160, 80, 80, 240],
);
assert.ok((signal?.rawP ?? 1) < 0.05);
assert.match(event.evidence.find((entry) => entry.method === "3Seq")?.calibration ?? "", /FindSubSeqTS\/Seq3PVals/);
assert.equal(message.comparisons, 1);
assert.deepEqual(message.tripletKernelCalls, { rdp: 0, geneconv: 0, sourceChi: 0, threeSeq: 1, siscan: 0 });

// The exploratory path has no designated references. CheckwrapC must retain
// the origin-spanning interval itself on a circular alignment, not silently
// replace it with the linear complement.
const circularRecords = ["", "", ""];
for (let site = 0; site < 240; site += 1) {
  const major = bases[site % 4];
  const minor = bases[(site + 1) % 4];
  circularRecords[0] += site < 70 || site >= 170 ? minor : major;
  circularRecords[1] += major;
  circularRecords[2] += minor;
}
const circularAlignment = parseAlignment(
  circularRecords.map((sequence, index) => `>circular_${index + 1}\n${sequence}`).join("\n"),
);
const circularResult = new Promise<Record<string, unknown>>((resolve, reject) => {
  runtime.postMessage = (payload: unknown) => {
    const next = payload as Record<string, unknown>;
    if (next.type === "result") resolve(next);
    if (next.type === "error") reject(new Error(String(next.message)));
  };
});
runtime.self.onmessage({
  data: {
    type: "analyze",
    jobId: 713,
    alignment: circularAlignment,
    options: {
      ...DEFAULT_OPTIONS,
      methods: ["3Seq"],
      circular: true,
      correction: "none",
      alpha: 0.05,
      minMethods: 1,
      cyclicDetection: false,
      polishBreakpoints: false,
      ancestralClustering: false,
      checkMisalignment: false,
    },
  },
});
const circularMessage = await circularResult;
const circularEvents = circularMessage.events as typeof events;
const circularEvent = circularEvents.find((candidate) => (
  candidate.recombinant === 0
  && candidate.start === 170
  && candidate.end === 70
  && candidate.wraps === true
));
assert.ok(circularEvent, "3Seq CheckwrapC must preserve an origin-spanning event without reference designations");
const circularSignal = circularEvent.methodSignals?.find((entry) => entry.method === "3Seq")?.sourceThreeSeq;
assert.equal(circularSignal?.sourceWrap, true);
assert.equal(circularSignal?.linearComplement, false);
assert.equal(circularMessage.comparisons, 1);
assert.deepEqual(circularMessage.tripletKernelCalls, { rdp: 0, geneconv: 0, sourceChi: 0, threeSeq: 1, siscan: 0 });

// CheckSplit3Seq/SubPVal: erase the middle of one otherwise continuous 3Seq
// signal. The two equal continuous sides must not become two events; the
// deterministic first-on-tie rule retains one piece and calibrates it with the
// original full-walk up/down counts.
const splitRecords = ["", "", ""];
for (let site = 0; site < 200; site += 1) {
  const major = bases[site % 4];
  const minor = bases[(site + 1) % 4];
  splitRecords[0] += site >= 50 && site < 150 ? minor : major;
  splitRecords[1] += major;
  splitRecords[2] += minor;
}
const splitAlignment = parseAlignment(
  splitRecords.map((sequence, index) => `>split_${index + 1}\n${sequence}`).join("\n"),
);
const splitResult = new Promise<Record<string, unknown>>((resolve, reject) => {
  runtime.postMessage = (payload: unknown) => {
    const next = payload as Record<string, unknown>;
    if (next.type === "result") resolve(next);
    if (next.type === "error") reject(new Error(String(next.message)));
  };
});
runtime.self.onmessage({
  data: {
    type: "analyze",
    jobId: 714,
    alignment: splitAlignment,
    disassemblyEvents: [{
      id: "erased-middle",
      decision: "accepted",
      evidenceStale: false,
      recombinant: 0,
      start: 90,
      end: 110,
      wraps: false,
    }],
    options: {
      ...DEFAULT_OPTIONS,
      methods: ["3Seq"],
      correction: "none",
      alpha: 0.05,
      minMethods: 1,
      cyclicDetection: false,
      polishBreakpoints: false,
      ancestralClustering: false,
      checkMisalignment: false,
    },
  },
});
const splitMessage = await splitResult;
const splitEvents = splitMessage.events as typeof events;
const splitEvent = splitEvents.find((candidate) => (
  candidate.recombinant === 0
  && candidate.start === 50
  && candidate.end === 90
  && candidate.structuralUncertainty?.pieces === 2
));
assert.ok(splitEvent, "CheckSplit3Seq must keep one lowest-p continuous side of an interrupted signal");
assert.equal(
  splitEvents.some((candidate) => candidate.recombinant === 0 && candidate.start === 110 && candidate.end === 150),
  false,
  "the other side is a SubPVal alternative, not a second 3Seq event",
);
const splitSignal = splitEvent.methodSignals?.find((entry) => entry.method === "3Seq")?.sourceThreeSeq;
assert.deepEqual(
  splitSignal && [splitSignal.splitRefined, splitSignal.fullDescent, splitSignal.descent, splitSignal.splitInformativeSites],
  [true, 100, 40, 40],
);
assert.ok(splitEvent.structuralUncertainty?.adjacentEventIds.includes("erased-middle"));
