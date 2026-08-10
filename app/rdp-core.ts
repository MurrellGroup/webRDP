export const PRIMARY_METHODS = [
  "RDP",
  "GENECONV",
  "BootScan",
  "MaxChi",
  "Chimaera",
  "SiScan",
  "3Seq",
] as const;

export type MethodName = (typeof PRIMARY_METHODS)[number];
export type EventDecision = "unreviewed" | "accepted" | "rejected";

export interface SequenceRecord {
  name: string;
  sequence: string;
  role?: "query" | "reference" | "both";
}

export interface GenomeFeature {
  id: string;
  type: string;
  start: number;
  end: number;
  strand: "+" | "-" | ".";
  phase?: 0 | 1 | 2;
  name: string;
  source: string;
  attributes: Record<string, string>;
}

export interface AlignmentData {
  name: string;
  format: "FASTA" | "CLUSTAL" | "PHYLIP" | "NEXUS" | "generated";
  sequences: SequenceRecord[];
  length: number;
  createdAt: number;
  features?: GenomeFeature[];
}

export interface MethodEvidence {
  method: MethodName;
  pValue: number;
  correctedP: number;
  score: number;
  supported: boolean;
  statistic: number;
  statisticLabel: string;
  calibration: string;
}

export interface BreakpointModel {
  method: "two-state-hmm" | "local-chi-square" | "manual";
  informativeSites: number;
  stateSwitches?: number;
  majorFit?: number;
  minorFit?: number;
}

export interface EventAuditEntry {
  id: string;
  timestamp: string;
  action: string;
  summary: string;
}

export interface EventDiagnostics {
  tractVariableDensity: number;
  backgroundVariableDensity: number;
  rateRatio: number;
  parentConflictRate: number;
  parentDiscriminatingSites: number;
  diffuseIncompatibility: boolean;
}

export interface AlignmentDiagnostics {
  sampledSequences: number;
  sampledBiallelicSites: number;
  testedSitePairs: number;
  incompatibleSitePairs: number;
  fourGameteFraction: number;
  nearIncompatibility: number;
  farIncompatibility: number;
  proximityRatio: number;
  ambiguityFraction: number;
}

export interface RdpEvent {
  id: string;
  recombinant: number;
  majorParent: number;
  minorParent: number;
  start: number;
  end: number;
  wraps: boolean;
  confidenceStart: [number, number];
  confidenceEnd: [number, number];
  breakpointModel?: BreakpointModel;
  evidence: MethodEvidence[];
  chiSquare: number;
  informativeSites: number;
  decision: EventDecision;
  warnings: string[];
  note: string;
  source: "wasm" | "example" | "manual";
  groupId: string | null;
  history: EventAuditEntry[];
  evidenceStale: boolean;
  diagnostics: EventDiagnostics;
}

export interface AnalysisOptions {
  mode: "exploratory" | "query-reference";
  circular: boolean;
  window: number;
  step: number;
  alpha: number;
  correction: "bonferroni" | "holm" | "none";
  minMethods: number;
  candidateParents: number;
  methods: MethodName[];
  exhaustive: boolean;
  polishBreakpoints: boolean;
  checkMisalignment: boolean;
  bootstrapReplicates: number;
  randomSeed: number;
}

export interface EvidencePoint {
  position: number;
  recombinantMajor: number;
  recombinantMinor: number;
  parentParent: number;
}

export interface NeighborJoiningNode {
  name?: string;
  length: number;
  children?: NeighborJoiningNode[];
}

export interface NeighborJoiningTree {
  root: NeighborJoiningNode;
  newick: string;
}

export interface HotspotTest {
  observedMaximum: number;
  expectedPerBin: number;
  empiricalP: number;
  bins: number;
  replicates: number;
}

const VALID = new Set(["A", "C", "G", "T", "U", "R", "Y", "S", "W", "K", "M", "B", "D", "H", "V", "N", "-", "?", ".", "!"]);

function cleanSequence(value: string): string {
  const sequence = value.toUpperCase().replace(/\s+/g, "").replace(/U/g, "T");
  for (const symbol of sequence) {
    if (!VALID.has(symbol)) {
      throw new Error(`Unsupported sequence symbol “${symbol}”.`);
    }
  }
  return sequence;
}

function uniqueNames(records: SequenceRecord[]): SequenceRecord[] {
  const seen = new Map<string, number>();
  return records.map((record, index) => {
    const base = record.name.trim() || `Sequence_${index + 1}`;
    const count = seen.get(base) ?? 0;
    seen.set(base, count + 1);
    return { ...record, name: count === 0 ? base : `${base}_${count + 1}` };
  });
}

function finalizeAlignment(
  name: string,
  format: AlignmentData["format"],
  records: SequenceRecord[],
): AlignmentData {
  const sequences = uniqueNames(
    records
      .map((record) => ({ ...record, sequence: cleanSequence(record.sequence) }))
      .filter((record) => record.sequence.length > 0),
  );
  if (sequences.length < 3) {
    throw new Error("An analysis needs at least three aligned sequences.");
  }
  const length = sequences[0].sequence.length;
  const mismatch = sequences.find((record) => record.sequence.length !== length);
  if (mismatch) {
    throw new Error(
      `Sequences must already be aligned. “${mismatch.name}” has ${mismatch.sequence.length.toLocaleString()} sites; expected ${length.toLocaleString()}.`,
    );
  }
  return { name, format, sequences, length, createdAt: Date.now(), features: [] };
}

function parseAttributes(value: string): Record<string, string> {
  const attributes: Record<string, string> = {};
  for (const field of value.split(";")) {
    const [rawKey, ...rawValue] = field.trim().split("=");
    if (!rawKey) continue;
    attributes[decodeURIComponent(rawKey)] = decodeURIComponent(rawValue.join("=") || "");
  }
  return attributes;
}

export function parseGenomeAnnotations(text: string, filename: string, length: number): GenomeFeature[] {
  const features: GenomeFeature[] = [];
  const lines = text.split(/\r?\n/);
  const isGff = lines.some((line) => line.split("\t").length >= 8);
  if (isGff) {
    for (const line of lines) {
      if (!line || line.startsWith("#")) continue;
      const columns = line.split("\t");
      if (columns.length < 8) continue;
      const start = Math.max(0, Number(columns[3]) - 1);
      const end = Math.min(length, Number(columns[4]));
      if (!Number.isFinite(start) || !Number.isFinite(end) || end <= start) continue;
      const attributes = parseAttributes(columns[8] ?? "");
      const phase = columns[7] === "0" || columns[7] === "1" || columns[7] === "2" ? Number(columns[7]) as 0 | 1 | 2 : undefined;
      features.push({
        id: attributes.ID || `${filename}-${features.length + 1}`,
        type: columns[2] || "feature",
        start,
        end,
        strand: columns[6] === "+" || columns[6] === "-" ? columns[6] : ".",
        phase,
        name: attributes.Name || attributes.gene || attributes.product || attributes.ID || `${columns[2]} ${features.length + 1}`,
        source: columns[1] || filename,
        attributes,
      });
    }
  } else if (/^\s*FEATURES\s+Location\/Qualifiers/m.test(text)) {
    let current: GenomeFeature | null = null;
    for (const line of lines) {
      const feature = line.match(/^\s{5}(\S+)\s+(.+)$/);
      if (feature) {
        const coordinates = [...feature[2].matchAll(/\d+/g)].map((match) => Number(match[0]));
        if (!coordinates.length) continue;
        const start = Math.max(0, Math.min(...coordinates) - 1);
        const end = Math.min(length, Math.max(...coordinates));
        if (end <= start) continue;
        current = {
          id: `${filename}-${features.length + 1}`,
          type: feature[1],
          start,
          end,
          strand: /complement\(/i.test(feature[2]) ? "-" : "+",
          name: `${feature[1]} ${features.length + 1}`,
          source: filename,
          attributes: { location: feature[2].trim() },
        };
        features.push(current);
        continue;
      }
      const qualifier = line.match(/^\s+\/(\w+)=(?:"([^"]*)"|(\S+))/);
      if (qualifier && current) {
        current.attributes[qualifier[1]] = qualifier[2] ?? qualifier[3] ?? "";
        if (["label", "gene", "product", "locus_tag"].includes(qualifier[1])) current.name = qualifier[2] ?? qualifier[3] ?? current.name;
      }
    }
  } else {
    for (const line of lines) {
      if (!line || line.startsWith("#")) continue;
      const columns = line.split(/\t|\s+/);
      if (columns.length < 3) continue;
      const start = Math.max(0, Number(columns[1]));
      const end = Math.min(length, Number(columns[2]));
      if (!Number.isFinite(start) || !Number.isFinite(end) || end <= start) continue;
      features.push({
        id: `${filename}-${features.length + 1}`,
        type: columns[4] || "region",
        start,
        end,
        strand: columns[5] === "+" || columns[5] === "-" ? columns[5] : ".",
        name: columns[3] || `region ${features.length + 1}`,
        source: filename,
        attributes: {},
      });
    }
  }
  if (!features.length) throw new Error("No GFF3, GenBank FEATURES, or BED annotations were recognized.");
  return features.sort((left, right) => left.start - right.start || left.end - right.end);
}

export function toGff3(features: GenomeFeature[], sequenceName = "alignment"): string {
  const lines = ["##gff-version 3"];
  for (const feature of features) {
    const attributes = {
      ...feature.attributes,
      ID: feature.id,
      Name: feature.name,
    };
    const encoded = Object.entries(attributes)
      .map(([key, value]) => `${encodeURIComponent(key)}=${encodeURIComponent(value)}`)
      .join(";");
    lines.push([
      sequenceName,
      feature.source || "RDP-Web",
      feature.type,
      feature.start + 1,
      feature.end,
      ".",
      feature.strand,
      feature.phase ?? ".",
      encoded,
    ].join("\t"));
  }
  return `${lines.join("\n")}\n`;
}

function parseFasta(text: string): SequenceRecord[] {
  const records: SequenceRecord[] = [];
  let current: SequenceRecord | null = null;
  for (const rawLine of text.split(/\r?\n/)) {
    const line = rawLine.trim();
    if (!line || line.startsWith(";")) continue;
    if (line.startsWith(">")) {
      if (current) records.push(current);
      current = { name: line.slice(1).trim(), sequence: "", role: "both" };
    } else {
      if (!current) throw new Error("FASTA sequence data appeared before its >header.");
      current.sequence += line;
    }
  }
  if (current) records.push(current);
  return records;
}

function parseClustal(text: string): SequenceRecord[] {
  const chunks = new Map<string, string>();
  for (const rawLine of text.split(/\r?\n/).slice(1)) {
    if (!rawLine.trim() || /^\s/.test(rawLine)) continue;
    const match = rawLine.trim().match(/^(\S+)\s+([A-Za-z?.!*\-]+)(?:\s+\d+)?$/);
    if (!match) continue;
    chunks.set(match[1], (chunks.get(match[1]) ?? "") + match[2]);
  }
  return [...chunks].map(([name, sequence]) => ({ name, sequence, role: "both" }));
}

function parsePhylip(text: string): SequenceRecord[] {
  const lines = text.split(/\r?\n/).filter((line) => line.trim());
  const header = lines.shift()?.trim().match(/^(\d+)\s+(\d+)/);
  if (!header) throw new Error("The PHYLIP header must contain sequence and site counts.");
  const expected = Number(header[1]);
  const chunks = new Map<string, string>();
  const order: string[] = [];
  let continuationIndex = 0;
  for (const rawLine of lines) {
    const line = rawLine.trim();
    const named = line.match(/^(\S+)\s+([A-Za-z?.!*\-\s]+)$/);
    if (named && (!chunks.has(named[1]) || chunks.size < expected)) {
      const name = named[1];
      if (!chunks.has(name)) order.push(name);
      chunks.set(name, (chunks.get(name) ?? "") + named[2].replace(/\s/g, ""));
      continue;
    }
    if (order.length) {
      const name = order[continuationIndex % order.length];
      chunks.set(name, (chunks.get(name) ?? "") + line.replace(/\s/g, ""));
      continuationIndex += 1;
    }
  }
  if (chunks.size !== expected) {
    throw new Error(`PHYLIP declares ${expected} sequences but ${chunks.size} were parsed.`);
  }
  return order.map((name) => ({ name, sequence: chunks.get(name) ?? "", role: "both" }));
}

function parseNexus(text: string): SequenceRecord[] {
  const matrix = text.match(/\bmatrix\b([\s\S]*?);/i)?.[1];
  if (!matrix) throw new Error("No MATRIX block was found in the NEXUS file.");
  const chunks = new Map<string, string>();
  for (const rawLine of matrix.split(/\r?\n/)) {
    const line = rawLine.replace(/\[[^\]]*\]/g, "").trim();
    if (!line) continue;
    const match = line.match(/^(?:'([^']+)'|"([^"]+)"|(\S+))\s+([A-Za-z?.!*\-\s]+)$/);
    if (!match) continue;
    const name = match[1] ?? match[2] ?? match[3];
    chunks.set(name, (chunks.get(name) ?? "") + match[4].replace(/\s/g, ""));
  }
  return [...chunks].map(([name, sequence]) => ({ name, sequence, role: "both" }));
}

export function parseAlignment(text: string, filename = "Pasted alignment"): AlignmentData {
  const trimmed = text.trim();
  if (!trimmed) throw new Error("The alignment is empty.");
  if (trimmed.startsWith(">")) return finalizeAlignment(filename, "FASTA", parseFasta(trimmed));
  if (/^(CLUSTAL|MUSCLE)/i.test(trimmed)) {
    return finalizeAlignment(filename, "CLUSTAL", parseClustal(trimmed));
  }
  if (/^#NEXUS/i.test(trimmed)) return finalizeAlignment(filename, "NEXUS", parseNexus(trimmed));
  if (/^\s*\d+\s+\d+/.test(trimmed)) {
    return finalizeAlignment(filename, "PHYLIP", parsePhylip(trimmed));
  }
  throw new Error("Format not recognized. Use aligned FASTA, CLUSTAL, PHYLIP, or NEXUS.");
}

function xorshift(seed: number): () => number {
  let value = seed >>> 0;
  return () => {
    value ^= value << 13;
    value ^= value >>> 17;
    value ^= value << 5;
    return (value >>> 0) / 4294967296;
  };
}

function mutate(sequence: string, rate: number, seed: number): string {
  const random = xorshift(seed);
  const bases = "ACGT";
  return [...sequence]
    .map((base) => {
      if (random() >= rate) return base;
      const alternatives = bases.replace(base, "");
      return alternatives[Math.floor(random() * alternatives.length)];
    })
    .join("");
}

export function makeDemoAlignment(): AlignmentData {
  const random = xorshift(0x5a17c0de);
  const ancestor = Array.from({ length: 2400 }, () => "ACGT"[Math.floor(random() * 4)]).join("");
  const alpha = mutate(ancestor, 0.072, 17);
  const beta = mutate(ancestor, 0.075, 29);
  const gamma = mutate(ancestor, 0.11, 31);
  const mosaic = mutate(alpha.slice(0, 782) + beta.slice(782, 1538) + alpha.slice(1538), 0.006, 43);
  const mosaic2 = mutate(alpha.slice(0, 795) + beta.slice(795, 1531) + alpha.slice(1531), 0.011, 47);
  return finalizeAlignment("Mosaic virus tutorial · 12 × 2,400 nt", "generated", [
    { name: "Mosaic-X", sequence: mosaic, role: "query" },
    { name: "Mosaic-Y", sequence: mosaic2, role: "query" },
    { name: "Alpha-01", sequence: alpha, role: "reference" },
    { name: "Alpha-02", sequence: mutate(alpha, 0.009, 53), role: "reference" },
    { name: "Alpha-03", sequence: mutate(alpha, 0.014, 59), role: "reference" },
    { name: "Beta-01", sequence: beta, role: "reference" },
    { name: "Beta-02", sequence: mutate(beta, 0.01, 61), role: "reference" },
    { name: "Beta-03", sequence: mutate(beta, 0.016, 67), role: "reference" },
    { name: "Gamma-01", sequence: gamma, role: "reference" },
    { name: "Gamma-02", sequence: mutate(gamma, 0.015, 71), role: "reference" },
    { name: "Alpha-04", sequence: mutate(alpha, 0.021, 73), role: "reference" },
    { name: "Beta-04", sequence: mutate(beta, 0.022, 79), role: "reference" },
  ]);
}

export function demoEvent(): RdpEvent {
  const pValues = [2.1e-31, 7.2e-24, 4.8e-27, 1.3e-29, 5.9e-26, 8.1e-19, 3.6e-22];
  return {
    id: "example-1",
    recombinant: 0,
    majorParent: 2,
    minorParent: 5,
    start: 782,
    end: 1538,
    wraps: false,
    confidenceStart: [770, 796],
    confidenceEnd: [1525, 1550],
    evidence: PRIMARY_METHODS.map((method, index) => ({
      method,
      pValue: pValues[index],
      correctedP: pValues[index] * 198,
      score: 30 - index * 1.6,
      supported: true,
      statistic: [0.82, 31, 0.97, 42.1, 38.6, 9.8, 13.4][index],
      statisticLabel: ["identity shift", "concordant run", "bootstrap topology support", "boundary χ²", "boundary χ²", "category Z", "maximum HGRW descent"][index],
      calibration: ["binomial", "G-scale 0 run", "seeded p-distance bootstrap + window sign", "χ²", "binary-triplet χ²", "category Z", "exact HGRW first-passage DP"][index],
    })),
    chiSquare: 132.4,
    informativeSites: 284,
    decision: "unreviewed",
    warnings: [],
    note: "Synthetic positive control: the known Beta-derived tract spans sites 783–1,538.",
    source: "example",
    groupId: null,
    history: [{
      id: "example-history-1",
      timestamp: "2026-01-01T00:00:00.000Z",
      action: "Loaded positive control",
      summary: "Known synthetic Alpha/Beta mosaic event.",
    }],
    evidenceStale: false,
    diagnostics: {
      tractVariableDensity: 0.13,
      backgroundVariableDensity: 0.12,
      rateRatio: 1.08,
      parentConflictRate: 0.01,
      parentDiscriminatingSites: 284,
      diffuseIncompatibility: false,
    },
  };
}

function canonicalBase(base: string): boolean {
  return base === "A" || base === "C" || base === "G" || base === "T";
}

export function eventSegments(
  event: Pick<RdpEvent, "start" | "end" | "wraps">,
  length: number,
): [number, number][] {
  const start = Math.max(0, Math.min(length, Math.trunc(event.start)));
  const end = Math.max(0, Math.min(length, Math.trunc(event.end)));
  if (event.wraps && start > end) {
    return [
      ...(start < length ? [[start, length] as [number, number]] : []),
      ...(end > 0 ? [[0, end] as [number, number]] : []),
    ];
  }
  return end > start ? [[start, end]] : [];
}

export function eventLength(
  event: Pick<RdpEvent, "start" | "end" | "wraps">,
  length: number,
): number {
  return eventSegments(event, length).reduce((total, [start, end]) => total + end - start, 0);
}

export function formatEventRegion(
  event: Pick<RdpEvent, "start" | "end" | "wraps">,
  length: number,
): string {
  if (event.wraps && event.start > event.end) {
    const tail = `${(event.start + 1).toLocaleString()}–${length.toLocaleString()}`;
    return event.end > 0
      ? `${tail} ↻ 1–${event.end.toLocaleString()}`
      : `${tail} ↻ origin`;
  }
  return `${(event.start + 1).toLocaleString()}–${event.end.toLocaleString()}`;
}

function identityCounts(
  a: string,
  b: string,
  start: number,
  end: number,
): { matches: number; valid: number } {
  let matches = 0;
  let valid = 0;
  for (let index = start; index < Math.min(end, a.length, b.length); index += 1) {
    if (!canonicalBase(a[index]) || !canonicalBase(b[index])) continue;
    valid += 1;
    if (a[index] === b[index]) matches += 1;
  }
  return { matches, valid };
}

export function pairwiseIdentity(a: string, b: string, start = 0, end = a.length): number {
  const { matches, valid } = identityCounts(a, b, start, end);
  return valid ? matches / valid : 0;
}

export function pairwiseIdentitySampled(a: string, b: string, maxSites = 4096): number {
  const length = Math.min(a.length, b.length);
  const samples = Math.min(length, maxSites);
  let matches = 0;
  let valid = 0;
  for (let sample = 0; sample < samples; sample += 1) {
    const site = samples <= 1 ? 0 : Math.floor(sample * (length - 1) / (samples - 1));
    if (!canonicalBase(a[site]) || !canonicalBase(b[site])) continue;
    valid += 1;
    if (a[site] === b[site]) matches += 1;
  }
  return valid ? matches / valid : 0;
}

function pairwiseIdentitySegmentsSampled(
  a: string,
  b: string,
  segments: [number, number][],
  maxSites = 8192,
): number {
  const totalSites = segments.reduce((total, [start, end]) => total + Math.max(0, end - start), 0);
  const samples = Math.min(totalSites, maxSites);
  let matches = 0;
  let valid = 0;
  for (let sample = 0; sample < samples; sample += 1) {
    let offset = samples <= 1 ? 0 : Math.floor(sample * (totalSites - 1) / (samples - 1));
    let site = 0;
    for (const [start, end] of segments) {
      const segmentLength = Math.max(0, end - start);
      if (offset < segmentLength) {
        site = start + offset;
        break;
      }
      offset -= segmentLength;
    }
    if (!canonicalBase(a[site]) || !canonicalBase(b[site])) continue;
    valid += 1;
    if (a[site] === b[site]) matches += 1;
  }
  return valid ? matches / valid : 0;
}

export function pairwiseIdentitySegments(
  a: string,
  b: string,
  segments: [number, number][],
): number {
  let matches = 0;
  let valid = 0;
  for (const [start, end] of segments) {
    const counts = identityCounts(a, b, start, end);
    matches += counts.matches;
    valid += counts.valid;
  }
  return valid ? matches / valid : 0;
}

function newickLabel(name: string): string {
  return /^[A-Za-z0-9_.-]+$/.test(name) ? name : `'${name.replaceAll("'", "''")}'`;
}

function nodeNewick(node: NeighborJoiningNode, root = false): string {
  const label = node.children?.length
    ? `(${node.children.map((child) => nodeNewick(child)).join(",")})`
    : newickLabel(node.name ?? "unnamed");
  return root ? `${label};` : `${label}:${Math.max(0, node.length).toFixed(8)}`;
}

export function neighborJoining(names: string[], distances: number[]): NeighborJoiningTree {
  if (names.length === 0) throw new Error("A tree needs at least one sequence.");
  if (distances.length !== names.length ** 2) throw new Error("The tree distance matrix has the wrong size.");
  if (names.length === 1) {
    const root = { name: names[0], length: 0 };
    return { root, newick: `${newickLabel(names[0])};` };
  }
  const nodes = new Map<number, NeighborJoiningNode>(names.map((name, index) => [index, { name, length: 0 }]));
  const matrix = new Map<string, number>();
  const key = (left: number, right: number) => left < right ? `${left}:${right}` : `${right}:${left}`;
  for (let left = 0; left < names.length; left += 1) {
    for (let right = left + 1; right < names.length; right += 1) {
      matrix.set(key(left, right), Math.max(0, distances[left * names.length + right]));
    }
  }
  let active = names.map((_, index) => index);
  let nextId = names.length;
  while (active.length > 2) {
    const sums = new Map<number, number>();
    for (const left of active) {
      let sum = 0;
      for (const right of active) if (left !== right) sum += matrix.get(key(left, right)) ?? 0;
      sums.set(left, sum);
    }
    let bestLeft = active[0];
    let bestRight = active[1];
    let bestQ = Number.POSITIVE_INFINITY;
    for (let leftIndex = 0; leftIndex < active.length; leftIndex += 1) {
      for (let rightIndex = leftIndex + 1; rightIndex < active.length; rightIndex += 1) {
        const left = active[leftIndex];
        const right = active[rightIndex];
        const q = (active.length - 2) * (matrix.get(key(left, right)) ?? 0)
          - (sums.get(left) ?? 0) - (sums.get(right) ?? 0);
        if (q < bestQ) {
          bestQ = q;
          bestLeft = left;
          bestRight = right;
        }
      }
    }
    const pairDistance = matrix.get(key(bestLeft, bestRight)) ?? 0;
    const denominator = 2 * (active.length - 2);
    const leftLength = Math.max(0, 0.5 * pairDistance + ((sums.get(bestLeft) ?? 0) - (sums.get(bestRight) ?? 0)) / denominator);
    const rightLength = Math.max(0, pairDistance - leftLength);
    const mergedId = nextId;
    nextId += 1;
    nodes.set(mergedId, {
      length: 0,
      children: [
        { ...(nodes.get(bestLeft) as NeighborJoiningNode), length: leftLength },
        { ...(nodes.get(bestRight) as NeighborJoiningNode), length: rightLength },
      ],
    });
    for (const other of active) {
      if (other === bestLeft || other === bestRight) continue;
      const distance = 0.5 * (
        (matrix.get(key(bestLeft, other)) ?? 0)
        + (matrix.get(key(bestRight, other)) ?? 0)
        - pairDistance
      );
      matrix.set(key(mergedId, other), Math.max(0, distance));
    }
    active = active.filter((id) => id !== bestLeft && id !== bestRight);
    active.push(mergedId);
  }
  const finalDistance = matrix.get(key(active[0], active[1])) ?? 0;
  const root: NeighborJoiningNode = {
    length: 0,
    children: [
      { ...(nodes.get(active[0]) as NeighborJoiningNode), length: Math.max(0, finalDistance / 2) },
      { ...(nodes.get(active[1]) as NeighborJoiningNode), length: Math.max(0, finalDistance / 2) },
    ],
  };
  return { root, newick: nodeNewick(root, true) };
}

export function buildLocalTree(
  alignment: AlignmentData,
  indexes: number[],
  segments: [number, number][],
): NeighborJoiningTree {
  const unique = [...new Set(indexes)].filter((index) => index >= 0 && index < alignment.sequences.length);
  const distances = new Array(unique.length ** 2).fill(0);
  for (let left = 0; left < unique.length; left += 1) {
    for (let right = left + 1; right < unique.length; right += 1) {
      const distance = 1 - pairwiseIdentitySegmentsSampled(
        alignment.sequences[unique[left]].sequence,
        alignment.sequences[unique[right]].sequence,
        segments,
      );
      distances[left * unique.length + right] = distance;
      distances[right * unique.length + left] = distance;
    }
  }
  return neighborJoining(unique.map((index) => alignment.sequences[index].name), distances);
}

export function breakpointHotspotTest(
  events: RdpEvent[],
  length: number,
  bins = 48,
  replicates = 999,
  seed = 1511506142,
): HotspotTest {
  const retained = events.filter((event) => event.decision !== "rejected");
  const counts = new Int32Array(bins);
  for (const event of retained) {
    counts[Math.min(bins - 1, Math.floor(event.start / Math.max(1, length) * bins))] += 1;
    counts[Math.min(bins - 1, Math.floor((event.end % Math.max(1, length)) / Math.max(1, length) * bins))] += 1;
  }
  const observedMaximum = Math.max(0, ...counts);
  const random = xorshift(seed);
  let exceedances = 0;
  for (let replicate = 0; replicate < replicates; replicate += 1) {
    const simulated = new Int32Array(bins);
    for (let breakpoint = 0; breakpoint < retained.length * 2; breakpoint += 1) {
      simulated[Math.min(bins - 1, Math.floor(random() * bins))] += 1;
    }
    if (Math.max(0, ...simulated) >= observedMaximum) exceedances += 1;
  }
  return {
    observedMaximum,
    expectedPerBin: retained.length * 2 / bins,
    empiricalP: (exceedances + 1) / (replicates + 1),
    bins,
    replicates,
  };
}

function windowSegments(center: number, half: number, length: number, circular: boolean): [number, number][] {
  if (!circular) return [[Math.max(0, center - half), Math.min(length, center + half)]];
  if (half * 2 >= length) return [[0, length]];
  const start = center - half;
  const end = center + half;
  if (start < 0) return [[0, end], [length + start, length]];
  if (end > length) return [[start, length], [0, end - length]];
  return [[start, end]];
}

export function evidenceProfile(
  alignment: AlignmentData,
  event: RdpEvent,
  window: number,
  points = 100,
  circular = event.wraps,
): EvidencePoint[] {
  const recombinant = alignment.sequences[event.recombinant]?.sequence ?? "";
  const major = alignment.sequences[event.majorParent]?.sequence ?? "";
  const minor = alignment.sequences[event.minorParent]?.sequence ?? "";
  const half = Math.max(10, Math.floor(window / 2));
  const step = Math.max(1, Math.floor(alignment.length / points));
  const profile: EvidencePoint[] = [];
  for (let position = 0; position < alignment.length; position += step) {
    const segments = windowSegments(position, half, alignment.length, circular);
    profile.push({
      position,
      recombinantMajor: pairwiseIdentitySegments(recombinant, major, segments),
      recombinantMinor: pairwiseIdentitySegments(recombinant, minor, segments),
      parentParent: pairwiseIdentitySegments(major, minor, segments),
    });
  }
  return profile;
}

export function alignmentStats(alignment: AlignmentData): {
  variableSites: number;
  gaps: number;
  ambiguities: number;
  meanIdentity: number;
  sampled: boolean;
} {
  const sequenceSampleCount = Math.min(alignment.sequences.length, 256);
  const siteSampleCount = Math.min(alignment.length, 10_000, Math.max(1, Math.floor(2_000_000 / Math.max(1, sequenceSampleCount))));
  const sequenceIndexes = Array.from({ length: sequenceSampleCount }, (_, index) => (
    sequenceSampleCount === alignment.sequences.length
      ? index
      : Math.floor(index * (alignment.sequences.length - 1) / Math.max(1, sequenceSampleCount - 1))
  ));
  let variableSites = 0;
  let gaps = 0;
  let ambiguities = 0;
  for (let sample = 0; sample < siteSampleCount; sample += 1) {
    const site = siteSampleCount === alignment.length
      ? sample
      : Math.floor(sample * (alignment.length - 1) / Math.max(1, siteSampleCount - 1));
    const bases = new Set<string>();
    for (const sequenceIndex of sequenceIndexes) {
      const record = alignment.sequences[sequenceIndex];
      const base = record.sequence[site];
      if (base === "-") gaps += 1;
      else if (!canonicalBase(base)) ambiguities += 1;
      else bases.add(base);
    }
    if (bases.size > 1) variableSites += 1;
  }
  const siteScale = alignment.length / Math.max(1, siteSampleCount);
  const cellScale = alignment.length * alignment.sequences.length / Math.max(1, siteSampleCount * sequenceSampleCount);
  variableSites = Math.round(variableSites * siteScale);
  gaps = Math.round(gaps * cellScale);
  ambiguities = Math.round(ambiguities * cellScale);
  let identity = 0;
  let pairs = 0;
  const identityCount = Math.min(alignment.sequences.length, 32);
  const identityIndexes = Array.from({ length: identityCount }, (_, index) => (
    identityCount === alignment.sequences.length
      ? index
      : Math.floor(index * (alignment.sequences.length - 1) / Math.max(1, identityCount - 1))
  ));
  for (let a = 0; a < identityIndexes.length; a += 1) {
    for (let b = a + 1; b < identityIndexes.length; b += 1) {
      identity += pairwiseIdentitySampled(
        alignment.sequences[identityIndexes[a]].sequence,
        alignment.sequences[identityIndexes[b]].sequence,
        2048,
      );
      pairs += 1;
    }
  }
  return {
    variableSites,
    gaps,
    ambiguities,
    meanIdentity: pairs ? identity / pairs : 1,
    sampled: sequenceSampleCount < alignment.sequences.length || siteSampleCount < alignment.length || identityIndexes.length < alignment.sequences.length || alignment.length > 2048,
  };
}

export function toFasta(records: SequenceRecord[], wrap = 80): string {
  return records
    .map((record) => {
      const chunks = record.sequence.match(new RegExp(`.{1,${wrap}}`, "g")) ?? [];
      return `>${record.name}\n${chunks.join("\n")}`;
    })
    .join("\n");
}

export function exportRecombinationFree(
  alignment: AlignmentData,
  events: RdpEvent[],
  mode: "remove" | "mask" | "mask-codon" | "split" | "partition",
): { filename: string; content: string }[] {
  const accepted = events.filter((event) => event.decision !== "rejected");
  if (mode === "remove") {
    const recombinantIndexes = new Set(accepted.map((event) => event.recombinant));
    return [{
      filename: "rdp-clean-nonrecombinant.fasta",
      content: toFasta(alignment.sequences.filter((_, index) => !recombinantIndexes.has(index))),
    }];
  }
  if (mode === "mask" || mode === "mask-codon") {
    const records = alignment.sequences.map((record, index) => {
      const sequence = [...record.sequence];
      for (const event of accepted.filter((candidate) => candidate.recombinant === index)) {
        for (const [start, end] of eventSegments(event, alignment.length)) {
          let maskStart = start;
          let maskEnd = end;
          if (mode === "mask-codon") {
            const coding = (alignment.features ?? []).find((feature) => feature.type.toLowerCase() === "cds" && feature.start < end && feature.end > start);
            const origin = coding
              ? (coding.strand === "-" ? coding.end - (coding.phase ?? 0) : coding.start + (coding.phase ?? 0))
              : 0;
            maskStart = Math.max(0, origin + Math.floor((start - origin) / 3) * 3);
            maskEnd = Math.min(alignment.length, origin + Math.ceil((end - origin) / 3) * 3);
          }
          for (let site = maskStart; site < maskEnd; site += 1) sequence[site] = "N";
        }
      }
      return { ...record, sequence: sequence.join("") };
    });
    return [{ filename: mode === "mask-codon" ? "rdp-clean-codon-masked.fasta" : "rdp-clean-masked.fasta", content: toFasta(records) }];
  }
  if (mode === "split") {
    const records: SequenceRecord[] = [];
    for (let index = 0; index < alignment.sequences.length; index += 1) {
      const record = alignment.sequences[index];
      const breakpoints = new Set([0, alignment.length]);
      accepted
        .filter((event) => event.recombinant === index)
        .forEach((event) => { breakpoints.add(event.start); breakpoints.add(event.end); });
      const sorted = [...breakpoints].sort((a, b) => a - b);
      for (let part = 0; part < sorted.length - 1; part += 1) {
        records.push({
          name: `${record.name}_${sorted[part] + 1}-${sorted[part + 1]}`,
          sequence: record.sequence.slice(sorted[part], sorted[part + 1]),
        });
      }
    }
    return [{ filename: "rdp-clean-split.fasta", content: toFasta(records) }];
  }
  const breakpoints = [...new Set([0, alignment.length, ...accepted.flatMap((event) => [event.start, event.end])])]
    .sort((a, b) => a - b);
  const files: { filename: string; content: string }[] = [];
  for (let part = 0; part < breakpoints.length - 1; part += 1) {
    const start = breakpoints[part];
    const end = breakpoints[part + 1];
    files.push({
      filename: `rdp-partition-${String(part + 1).padStart(2, "0")}-${start + 1}-${end}.fasta`,
      content: toFasta(alignment.sequences.map((record) => ({ ...record, sequence: record.sequence.slice(start, end) }))),
    });
  }
  return files;
}

export const DEFAULT_OPTIONS: AnalysisOptions = {
  mode: "exploratory",
  circular: false,
  window: 120,
  step: 5,
  alpha: 0.05,
  correction: "bonferroni",
  minMethods: 3,
  candidateParents: 8,
  methods: [...PRIMARY_METHODS],
  exhaustive: false,
  polishBreakpoints: true,
  checkMisalignment: true,
  bootstrapReplicates: 100,
  randomSeed: 1511506142,
};

export interface RdpProject {
  schema: "rdp-web/0.4";
  alignment: AlignmentData;
  options: AnalysisOptions;
  events: RdpEvent[];
  metrics: {
    elapsedMs: number;
    comparisons: number;
    engine: string;
    matrixMode?: string;
    parentSamples?: number;
    timing?: { distanceMs: number; scanMs: number; statisticsMs: number; diagnosticsMs?: number };
    diagnostics?: AlignmentDiagnostics;
  } | null;
  distance: number[];
}

function finiteNumber(value: unknown, fallback: number): number {
  return typeof value === "number" && Number.isFinite(value) ? value : fallback;
}

export function serializeProject(project: Omit<RdpProject, "schema">): string {
  return JSON.stringify({ schema: "rdp-web/0.4", ...project }, null, 2);
}

export function parseProject(text: string): RdpProject {
  const raw = JSON.parse(text) as Record<string, unknown>;
  if (typeof raw.schema !== "string" || !raw.schema.startsWith("rdp-web/")) {
    throw new Error("This JSON file is not an RDP Web project.");
  }
  const rawAlignment = raw.alignment as Partial<AlignmentData> | undefined;
  if (!rawAlignment || !Array.isArray(rawAlignment.sequences)) {
    throw new Error("The project does not contain an alignment.");
  }
  const format = (["FASTA", "CLUSTAL", "PHYLIP", "NEXUS", "generated"] as const)
    .includes(rawAlignment.format as AlignmentData["format"])
    ? rawAlignment.format as AlignmentData["format"]
    : "generated";
  const alignmentBase = finalizeAlignment(
    typeof rawAlignment.name === "string" ? rawAlignment.name : "Imported RDP Web project",
    format,
    rawAlignment.sequences.map((record, index) => {
      const candidate = record as Partial<SequenceRecord>;
      const role = candidate.role === "query" || candidate.role === "reference" || candidate.role === "both"
        ? candidate.role
        : "both";
      return {
        name: typeof candidate.name === "string" ? candidate.name : `Sequence_${index + 1}`,
        sequence: typeof candidate.sequence === "string" ? candidate.sequence : "",
        role,
      };
    }),
  );
  const features = Array.isArray(rawAlignment.features) ? rawAlignment.features.flatMap((value, index): GenomeFeature[] => {
    const feature = value as Partial<GenomeFeature>;
    const start = Math.max(0, Math.min(alignmentBase.length - 1, Math.trunc(finiteNumber(feature.start, 0))));
    const end = Math.max(start + 1, Math.min(alignmentBase.length, Math.trunc(finiteNumber(feature.end, start + 1))));
    return [{
      id: typeof feature.id === "string" ? feature.id : `imported-feature-${index + 1}`,
      type: typeof feature.type === "string" ? feature.type : "feature",
      start,
      end,
      strand: feature.strand === "+" || feature.strand === "-" ? feature.strand : ".",
      phase: feature.phase === 0 || feature.phase === 1 || feature.phase === 2 ? feature.phase : undefined,
      name: typeof feature.name === "string" ? feature.name : `feature ${index + 1}`,
      source: typeof feature.source === "string" ? feature.source : "imported",
      attributes: feature.attributes && typeof feature.attributes === "object" ? feature.attributes : {},
    }];
  }) : [];
  const alignment: AlignmentData = { ...alignmentBase, features };
  const rawOptions = (raw.options ?? {}) as Partial<AnalysisOptions>;
  const methods = Array.isArray(rawOptions.methods)
    ? rawOptions.methods.filter((method): method is MethodName => PRIMARY_METHODS.includes(method as MethodName))
    : [...PRIMARY_METHODS];
  const options: AnalysisOptions = {
    ...DEFAULT_OPTIONS,
    ...rawOptions,
    methods: methods.length ? methods : [...PRIMARY_METHODS],
    mode: rawOptions.mode === "query-reference" ? "query-reference" : "exploratory",
    correction: rawOptions.correction === "holm" || rawOptions.correction === "none" ? rawOptions.correction : "bonferroni",
    bootstrapReplicates: Math.max(0, Math.min(1000, Math.trunc(finiteNumber(rawOptions.bootstrapReplicates, DEFAULT_OPTIONS.bootstrapReplicates)))),
    randomSeed: Math.trunc(finiteNumber(rawOptions.randomSeed, DEFAULT_OPTIONS.randomSeed)) >>> 0,
  };
  const rawEvents = Array.isArray(raw.events) ? raw.events : [];
  const events = rawEvents.flatMap((value, index): RdpEvent[] => {
    const event = value as Partial<RdpEvent>;
    const recombinant = Math.trunc(finiteNumber(event.recombinant, -1));
    const majorParent = Math.trunc(finiteNumber(event.majorParent, -1));
    const minorParent = Math.trunc(finiteNumber(event.minorParent, -1));
    if ([recombinant, majorParent, minorParent].some((sequence) => sequence < 0 || sequence >= alignment.sequences.length)) return [];
    const start = Math.max(0, Math.min(alignment.length - 1, Math.trunc(finiteNumber(event.start, 0))));
    const requestedWrap = event.wraps === true;
    const rawEnd = Math.trunc(finiteNumber(event.end, alignment.length));
    const circularEnd = Math.max(0, Math.min(alignment.length - 1, rawEnd));
    const linearEnd = Math.max(start + 1, Math.min(alignment.length, rawEnd));
    const wraps = requestedWrap && start > circularEnd;
    const end = wraps ? circularEnd : linearEnd;
    const evidence = Array.isArray(event.evidence) ? event.evidence.flatMap((value): MethodEvidence[] => {
      const item = value as Partial<MethodEvidence>;
      if (!PRIMARY_METHODS.includes(item.method as MethodName)) return [];
      const pValue = finiteNumber(item.pValue, 1);
      const correctedP = finiteNumber(item.correctedP, pValue);
      return [{
        method: item.method as MethodName,
        pValue,
        correctedP,
        score: finiteNumber(item.score, -Math.log10(Math.max(Number.MIN_VALUE, pValue))),
        supported: item.supported === true,
        statistic: finiteNumber(item.statistic, 0),
        statisticLabel: typeof item.statisticLabel === "string" ? item.statisticLabel : "imported statistic",
        calibration: typeof item.calibration === "string" ? item.calibration : "legacy project",
      }];
    }) : [];
    const confidence = Math.max(2, Math.floor(options.window / 12));
    const history = Array.isArray(event.history) ? event.history.flatMap((value, historyIndex): EventAuditEntry[] => {
      const entry = value as Partial<EventAuditEntry>;
      if (typeof entry.action !== "string") return [];
      return [{
        id: typeof entry.id === "string" ? entry.id : `imported-${index + 1}-history-${historyIndex + 1}`,
        timestamp: typeof entry.timestamp === "string" ? entry.timestamp : new Date(0).toISOString(),
        action: entry.action,
        summary: typeof entry.summary === "string" ? entry.summary : "Imported project history entry.",
      }];
    }) : [];
    return [{
      id: typeof event.id === "string" ? event.id : `imported-${index + 1}`,
      recombinant,
      majorParent,
      minorParent,
      start,
      end,
      wraps,
      confidenceStart: Array.isArray(event.confidenceStart) && event.confidenceStart.length === 2
        ? [finiteNumber(event.confidenceStart[0], Math.max(0, start - confidence)), finiteNumber(event.confidenceStart[1], Math.min(alignment.length, start + confidence))]
        : [Math.max(0, start - confidence), Math.min(alignment.length, start + confidence)],
      confidenceEnd: Array.isArray(event.confidenceEnd) && event.confidenceEnd.length === 2
        ? [finiteNumber(event.confidenceEnd[0], Math.max(0, end - confidence)), finiteNumber(event.confidenceEnd[1], Math.min(alignment.length, end + confidence))]
        : [Math.max(0, end - confidence), Math.min(alignment.length, end + confidence)],
      breakpointModel: event.breakpointModel && typeof event.breakpointModel === "object"
        ? {
            method: event.breakpointModel.method === "two-state-hmm" || event.breakpointModel.method === "manual"
              ? event.breakpointModel.method
              : "local-chi-square",
            informativeSites: Math.trunc(finiteNumber(event.breakpointModel.informativeSites, 0)),
            stateSwitches: Math.trunc(finiteNumber(event.breakpointModel.stateSwitches, 0)) || undefined,
            majorFit: finiteNumber(event.breakpointModel.majorFit, 0) || undefined,
            minorFit: finiteNumber(event.breakpointModel.minorFit, 0) || undefined,
          }
        : undefined,
      evidence,
      chiSquare: finiteNumber(event.chiSquare, 0),
      informativeSites: Math.trunc(finiteNumber(event.informativeSites, 0)),
      decision: event.decision === "accepted" || event.decision === "rejected" ? event.decision : "unreviewed",
      warnings: Array.isArray(event.warnings) ? event.warnings.filter((item): item is string => typeof item === "string") : [],
      note: typeof event.note === "string" ? event.note : "",
      source: event.source === "manual" || event.source === "example" ? event.source : "wasm",
      groupId: typeof event.groupId === "string" && event.groupId.trim() ? event.groupId.trim() : null,
      history,
      evidenceStale: event.evidenceStale === true,
      diagnostics: event.diagnostics && typeof event.diagnostics === "object"
        ? {
            tractVariableDensity: finiteNumber(event.diagnostics.tractVariableDensity, 0),
            backgroundVariableDensity: finiteNumber(event.diagnostics.backgroundVariableDensity, 0),
            rateRatio: finiteNumber(event.diagnostics.rateRatio, 1),
            parentConflictRate: finiteNumber(event.diagnostics.parentConflictRate, 0),
            parentDiscriminatingSites: Math.trunc(finiteNumber(event.diagnostics.parentDiscriminatingSites, 0)),
            diffuseIncompatibility: event.diagnostics.diffuseIncompatibility === true,
          }
        : {
            tractVariableDensity: 0,
            backgroundVariableDensity: 0,
            rateRatio: 1,
            parentConflictRate: 0,
            parentDiscriminatingSites: 0,
            diffuseIncompatibility: false,
          },
    }];
  });
  const rawMetrics = raw.metrics as RdpProject["metrics"] | undefined;
  const metrics = rawMetrics && typeof rawMetrics.engine === "string"
    ? {
        elapsedMs: finiteNumber(rawMetrics.elapsedMs, 0),
        comparisons: Math.trunc(finiteNumber(rawMetrics.comparisons, 0)),
        engine: rawMetrics.engine,
        matrixMode: typeof rawMetrics.matrixMode === "string" ? rawMetrics.matrixMode : undefined,
        parentSamples: finiteNumber(rawMetrics.parentSamples, 0) || undefined,
        timing: rawMetrics.timing && typeof rawMetrics.timing === "object"
          ? {
              distanceMs: finiteNumber(rawMetrics.timing.distanceMs, 0),
              scanMs: finiteNumber(rawMetrics.timing.scanMs, 0),
              statisticsMs: finiteNumber(rawMetrics.timing.statisticsMs, 0),
              diagnosticsMs: finiteNumber(rawMetrics.timing.diagnosticsMs, 0) || undefined,
            }
          : undefined,
        diagnostics: rawMetrics.diagnostics && typeof rawMetrics.diagnostics === "object"
          ? {
              sampledSequences: Math.trunc(finiteNumber(rawMetrics.diagnostics.sampledSequences, 0)),
              sampledBiallelicSites: Math.trunc(finiteNumber(rawMetrics.diagnostics.sampledBiallelicSites, 0)),
              testedSitePairs: Math.trunc(finiteNumber(rawMetrics.diagnostics.testedSitePairs, 0)),
              incompatibleSitePairs: Math.trunc(finiteNumber(rawMetrics.diagnostics.incompatibleSitePairs, 0)),
              fourGameteFraction: finiteNumber(rawMetrics.diagnostics.fourGameteFraction, 0),
              nearIncompatibility: finiteNumber(rawMetrics.diagnostics.nearIncompatibility, 0),
              farIncompatibility: finiteNumber(rawMetrics.diagnostics.farIncompatibility, 0),
              proximityRatio: finiteNumber(rawMetrics.diagnostics.proximityRatio, 1),
              ambiguityFraction: finiteNumber(rawMetrics.diagnostics.ambiguityFraction, 0),
            }
          : undefined,
      }
    : null;
  const distance = Array.isArray(raw.distance)
    ? raw.distance.map((value) => finiteNumber(value, 0))
    : [];
  return { schema: "rdp-web/0.4", alignment, options, events, metrics, distance };
}
