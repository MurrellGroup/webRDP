import fs from "node:fs";
import path from "node:path";
import { fileURLToPath } from "node:url";
import {
  DEFAULT_OPTIONS,
  MethodName,
  PRIMARY_METHODS,
  SOURCE_READY_METHODS,
  parseAlignment,
  serializeProject,
} from "../app/rdp-core";

const HELP = `RDP Web batch runner

Usage:
  npm run cli -- input.fasta [output.rdpweb] [options]

Options:
  --mode exploratory|query-reference
  --circular
  --methods RDP,MaxChi,Chimaera,SiScan
  --min-methods N
  --approximate-parent-shortlist
  --candidate-parents N          (only with approximate shortlist)
  --one-pass                     disable RDP5 erase/extract detection cycles
  --max-detection-cycles N       safety cap; default 250
  --chi-signals N                MAXCHI/CHIMAERA peak pairs retained/triplet
  --bootstrap N
  --seed N
  --siscan-outgroup nearest|most-divergent|randomized
  --siscan-outgroup-sequence N   (1-based; selects manual mode)
  --siscan-positions triplet-variable|quartet-variable|all
  --siscan-gaps strip|fifth-state
  --siscan-window N             source default 200 alignment positions
  --siscan-step N               source default 20 alignment positions
  --siscan-scan-permutations N
  --siscan-final-permutations N
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
const requestedNames = option("--methods")?.split(",").filter(Boolean);
const unavailable = requestedNames?.filter((method) => PRIMARY_METHODS.includes(method as MethodName) && !SOURCE_READY_METHODS.includes(method as MethodName)) ?? [];
if (unavailable.length) throw new Error(`${unavailable.join(", ")} ${unavailable.length === 1 ? "is" : "are"} disabled until the complete author-source batch port is available.`);
const requestedMethods = requestedNames?.filter((method): method is MethodName => SOURCE_READY_METHODS.includes(method as MethodName));
const methods = requestedMethods?.length ? requestedMethods : [...SOURCE_READY_METHODS];
const requestedOutgroup = option("--siscan-outgroup");
const manualOutgroup = Number(option("--siscan-outgroup-sequence"));
const requestedPositions = option("--siscan-positions");
const requestedGapMode = option("--siscan-gaps");
const siskanScanPermutations = Math.max(2, Math.min(1000, Number(option("--siscan-scan-permutations") ?? DEFAULT_OPTIONS.siskanScanPermutations)));
const options = {
  ...DEFAULT_OPTIONS,
  mode: option("--mode") === "query-reference" ? "query-reference" as const : "exploratory" as const,
  circular: arguments_.includes("--circular"),
  exhaustive: !arguments_.includes("--approximate-parent-shortlist"),
  cyclicDetection: !arguments_.includes("--one-pass"),
  maximumDetectionCycles: Math.max(1, Math.min(1000, Number(option("--max-detection-cycles") ?? DEFAULT_OPTIONS.maximumDetectionCycles))),
  methods,
  minMethods: Math.min(methods.length, Math.max(1, Number(option("--min-methods") ?? DEFAULT_OPTIONS.minMethods))),
  candidateParents: Math.max(3, Number(option("--candidate-parents") ?? DEFAULT_OPTIONS.candidateParents)),
  chiSignalsPerTriplet: Math.max(1, Math.min(256, Number(option("--chi-signals") ?? DEFAULT_OPTIONS.chiSignalsPerTriplet))),
  bootstrapReplicates: Math.max(0, Number(option("--bootstrap") ?? DEFAULT_OPTIONS.bootstrapReplicates)),
  randomSeed: Number(option("--seed") ?? DEFAULT_OPTIONS.randomSeed) >>> 0,
  siskanOutgroupMode: Number.isFinite(manualOutgroup) && manualOutgroup >= 1
    ? "manual" as const
    : requestedOutgroup === "most-divergent" || requestedOutgroup === "randomized"
      ? requestedOutgroup
      : "nearest" as const,
  siskanOutgroupSequence: Number.isFinite(manualOutgroup) && manualOutgroup >= 1 && manualOutgroup <= alignment.sequences.length
    ? Math.trunc(manualOutgroup - 1)
    : null,
  siskanPositionMode: requestedPositions === "quartet-variable" || requestedPositions === "all"
    ? requestedPositions
    : "triplet-variable" as const,
  siskanGapMode: requestedGapMode === "fifth-state" ? "fifth-state" as const : "strip" as const,
  siskanWindow: Math.max(12, Math.min(alignment.length, Number(option("--siscan-window") ?? DEFAULT_OPTIONS.siskanWindow))),
  siskanStep: Math.max(1, Math.min(alignment.length, Number(option("--siscan-step") ?? DEFAULT_OPTIONS.siskanStep))),
  siskanScanPermutations,
  siskanPValuePermutations: Math.max(
    siskanScanPermutations,
    Math.min(10_000, Number(option("--siscan-final-permutations") ?? DEFAULT_OPTIONS.siskanPValuePermutations)),
  ),
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
runtime.self.onmessage({ data: { type: "analyze", jobId: 1, alignment, options, cyclicDetection: options.cyclicDetection } });
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
  disassembly: message.disassembly as NonNullable<Parameters<typeof serializeProject>[0]["metrics"]>["disassembly"],
  rdpSignalTruncations: Number(message.rdpSignalTruncations) || undefined,
  geneconvSignalTruncations: Number(message.geneconvSignalTruncations) || undefined,
  chiSignalTruncations: Number(message.chiSignalTruncations) || undefined,
  tripletKernelCalls: message.tripletKernelCalls as NonNullable<Parameters<typeof serializeProject>[0]["metrics"]>["tripletKernelCalls"],
  detectionCycle: message.detectionCycle as NonNullable<Parameters<typeof serializeProject>[0]["metrics"]>["detectionCycle"],
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
