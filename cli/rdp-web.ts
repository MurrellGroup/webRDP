import fs from "node:fs";
import path from "node:path";
import { fileURLToPath } from "node:url";
import {
  DEFAULT_OPTIONS,
  MethodName,
  PRIMARY_METHODS,
  parseAlignment,
  serializeProject,
} from "../app/rdp-core";

const HELP = `RDP Web batch runner

Usage:
  npm run cli -- input.fasta [output.rdpweb] [options]

Options:
  --mode exploratory|query-reference
  --circular
  --methods RDP,GENECONV,BootScan,MaxChi,Chimaera,SiScan,3Seq
  --min-methods N
  --candidate-parents N
  --bootstrap N
  --seed N
  --exhaustive
  --help

The CLI runs the same checked-in WebAssembly and probability code as the
GitHub Pages application. It writes a restorable .rdpweb 0.5 project.
`;

const arguments_ = process.argv.slice(2);
if (!arguments_.length || arguments_.includes("--help") || arguments_.includes("-h")) {
  process.stdout.write(HELP);
  process.exit(0);
}

function option(name: string): string | undefined {
  const index = arguments_.indexOf(name);
  return index >= 0 ? arguments_[index + 1] : undefined;
}

const positional = arguments_.filter((value, index) => !value.startsWith("--") && (index === 0 || !arguments_[index - 1].startsWith("--")));
const inputPath = path.resolve(positional[0]);
const outputPath = path.resolve(positional[1] ?? `${inputPath.replace(/\.[^.]+$/u, "")}.rdpweb`);
const input = fs.readFileSync(inputPath, "utf8");
const alignment = parseAlignment(input, path.basename(inputPath));
const requestedMethods = option("--methods")?.split(",").filter((method): method is MethodName => PRIMARY_METHODS.includes(method as MethodName));
const options = {
  ...DEFAULT_OPTIONS,
  mode: option("--mode") === "query-reference" ? "query-reference" as const : "exploratory" as const,
  circular: arguments_.includes("--circular"),
  exhaustive: arguments_.includes("--exhaustive"),
  methods: requestedMethods?.length ? requestedMethods : [...PRIMARY_METHODS],
  minMethods: Math.max(1, Number(option("--min-methods") ?? DEFAULT_OPTIONS.minMethods)),
  candidateParents: Math.max(3, Number(option("--candidate-parents") ?? DEFAULT_OPTIONS.candidateParents)),
  bootstrapReplicates: Math.max(0, Number(option("--bootstrap") ?? DEFAULT_OPTIONS.bootstrapReplicates)),
  randomSeed: Number(option("--seed") ?? DEFAULT_OPTIONS.randomSeed) >>> 0,
};

const wasmPath = fileURLToPath(new URL("../public/wasm/rdp.wasm", import.meta.url));
const wasm = fs.readFileSync(wasmPath);
const runtime = globalThis as typeof globalThis & {
  self: typeof globalThis & { onmessage: (event: { data: unknown }) => void };
  postMessage: (payload: unknown) => void;
};
runtime.self = runtime;
runtime.fetch = async () => new Response(wasm, { headers: { "content-type": "application/wasm" } });

const result = new Promise<Record<string, unknown>>((resolve, reject) => {
  runtime.postMessage = (payload: unknown) => {
    const message = payload as Record<string, unknown>;
    if (message.type === "progress") process.stderr.write(`\r${String(message.phase).padEnd(58)} ${Math.round(Number(message.progress) * 100)}%`);
    if (message.type === "result") resolve(message);
    if (message.type === "error") reject(new Error(String(message.message)));
  };
});

await import("../public/rdp-worker.js");
runtime.self.onmessage({ data: { type: "analyze", jobId: 1, alignment, options } });
const message = await result;
process.stderr.write("\n");
const events = message.events as ReturnType<typeof JSON.parse>;
const metrics = {
  elapsedMs: Number(message.elapsedMs),
  comparisons: Number(message.comparisons),
  engine: String(message.engine),
  matrixMode: String(message.matrixMode),
  parentSamples: Number(message.parentSamples),
  timing: message.timing as { distanceMs: number; scanMs: number; statisticsMs: number; diagnosticsMs?: number },
  diagnostics: message.diagnostics as NonNullable<Parameters<typeof serializeProject>[0]["metrics"]>["diagnostics"],
};
fs.writeFileSync(outputPath, serializeProject({
  alignment,
  options,
  events,
  metrics,
  distance: message.distance as number[],
  auditLog: [{
    id: `cli-${Date.now()}`,
    timestamp: new Date().toISOString(),
    action: "Completed CLI batch scan",
    summary: `${metrics.comparisons.toLocaleString()} triplets tested; ${events.length} hypotheses retained.`,
  }],
}));
process.stdout.write(`${outputPath}\n${events.length} hypotheses · ${metrics.comparisons.toLocaleString()} triplets · ${metrics.elapsedMs.toFixed(1)} ms\n`);
