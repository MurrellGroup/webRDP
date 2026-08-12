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
// Only these families are allowed in production scans until their complete
// author-source discovery path has been ported. This prevents a simplified
// stand-in from silently participating in an analysis labelled RDP parity.
export const SOURCE_READY_METHODS: MethodName[] = ["RDP", "GENECONV", "BootScan", "MaxChi", "Chimaera", "SiScan", "3Seq"];
export type EventDecision = "unreviewed" | "accepted" | "rejected";

export interface SequenceRecord {
  name: string;
  sequence: string;
  role?: "query" | "reference" | "both";
  referenceGroup?: string;
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
  correctionScope?: string;
}

export interface MethodSignal {
  method: MethodName | "shared-screen";
  start: number;
  end: number;
  wraps: boolean;
  statistic: number;
  locator: string;
  sourceRoutine?: string;
  sourceChi?: {
    track: number;
    targetSlot: number | null;
    informativeSites: number;
    halfWindow: number;
    boundaryStatistics: [number, number];
    boundaryRanks: [number, number];
    growthWidths: [number, number];
    direction: 1 | -1;
  };
  sourceGeneconv?: {
    track: number;
    targetSlot: number;
    minorSlot: number;
    majorSlot: number;
    fragmentScore: number;
    informativeSites: number;
    matchingSites: number;
    mismatchSites: number;
    mismatchPenalty: number;
    rawP: number;
    startRank: number;
    endRank: number;
  };
  sourceBootscan?: {
    topology: number;
    baselineTopology: number;
    bootstrapSupport: number;
    bootstrapReplicates: number;
    runWindows: number;
    tractPairMatches: number;
    backgroundPairMatches: number;
    tractInformativeSites: number;
    informativeSites: number;
    rawP: number;
    window: number;
    step: number;
    relationshipMode: "distance";
  };
  sourceThreeSeq?: {
    target: number;
    majorParent: number;
    minorParent: number;
    direction: 1 | -1;
    upSteps: number;
    downSteps: number;
    descent: number;
    informativeSites: number;
    cycle: number;
    rawStart: number;
    rawEnd: number;
    rawP: number;
    probabilityMode: "exact-table" | "siegmund-discrete" | "scaled-table" | "unavailable";
    sourceWrap?: boolean;
    linearComplement?: boolean;
    splitRefined?: boolean;
    fullDescent?: number;
    splitInformativeSites?: number;
  };
  sourceSiScan?: {
    rawP: number;
    rawStart: number;
    rawEnd: number;
    runWindows: number;
    outgroupSourcePath: string;
    positionMode: "triplet-variable" | "quartet-variable" | "all";
    gapMode: "strip" | "fifth-state";
    window: number;
    step: number;
    topologyTriplet: [number, number, number];
    recombinant: number;
    majorParent: number;
    minorParent: number;
  };
  outgroup?: number | null;
  outgroupMode?: "nearest" | "most-divergent" | "randomized" | "manual";
  outgroupSampled?: boolean;
  permutations?: number;
  scanPermutations?: number;
  pattern?: number;
  scoreFamily?: "pattern" | "sum";
  baselineTopology?: number;
  inferredTopology?: number;
  profile?: Array<{
    position: number;
    z: number;
    topology: number;
    baselineTopology: number;
    pattern?: number;
    scoreFamily?: "pattern" | "sum";
  }>;
}

export interface BreakpointModel {
  method: "burt-hmm" | "two-state-hmm" | "local-chi-square" | "manual";
  informativeSites: number;
  stateSwitches?: number;
  majorFit?: number;
  minorFit?: number;
  states?: number;
  logLikelihood?: number;
  viterbiLogLikelihood?: number;
  bic?: number;
  aic?: number;
  criterion?: string;
  randomStarts?: number;
  iterations?: number;
  winningRestart?: number;
  selectedState?: number;
  posteriorThreshold?: number;
  sourceParity?: boolean;
  sourceCompatibility?: string;
  sourceRoutines?: string[];
  sequenceOrder?: number[];
  stateDominantCategories?: number[];
  circularPadding?: { offset: number; fittedSites: number; croppedSites: number };
  candidateBreakpoints?: [number, number];
  polishedBreakpoints?: [number, number];
  polishDecision?: {
    startAdopted: boolean;
    endAdopted: boolean;
    sameSwitchResolved: boolean;
    startMissingBoundary?: boolean;
    endMissingBoundary?: boolean;
    revertedForInformation: boolean;
    insideVariableSites: number;
    outsideVariableSites: number;
    startWithin99: boolean;
    endWithin99: boolean;
    startVariableSiteDistance?: number;
    endVariableSiteDistance?: number;
  };
  confidence99Start?: [number, number];
  confidence99End?: [number, number];
  emissions?: number[][];
  transitions?: number[][];
  switches?: Array<{ position: number; informativeIndex?: number; fromState: number; toState: number; confidence95: [number, number]; confidence99?: [number, number]; sourceCoordinates?: number[]; matchedStart?: boolean; matchedEnd?: boolean }>;
  posteriorTrace?: Array<{ position: number; informativeIndex?: number; state: number; probabilities: number[] }>;
  modelSelection?: Array<{ states: number; logLikelihood: number; bic: number; aic: number; iterations: number; winningRestart: number }>;
}

export interface AncestralCluster {
  inference: "rdp5-three-set" | "manual";
  representativeId: string;
  memberEventIds: string[];
  sequenceMembers: number[];
  confidence: number;
  evidenceCounts: { phylogenetic: number; distance: number; detectableSignal: number; sourceSimilarity?: number };
  partialOverprint: boolean;
  sourceMerge?: { threshold: number; pairDistances: Array<{ eventIds: [string, string]; distance: number; belowThreshold: boolean }> };
  pairwise: Array<Record<string, unknown>>;
}

export interface CoRecombinantSet {
  presumedRecombinant: number;
  parents: number[];
  sequenceMembers: number[];
  testedSequences: number;
  requiredEvidenceSets: number;
  evidence: Array<{
    sequence: number;
    sets: number;
    phylogenetic: boolean;
    distance: boolean;
    detectableSignal: boolean;
    bestCorrelation?: { r: number; pValue: number; inversion: number };
    topologyMargin?: number;
    treeBootstrap?: {
      replicates: number;
      cutoff: number;
      cohortTaxa: number;
      sourceSequenceCount: number;
      included: boolean;
      exactSiteBootstrap: boolean;
      sourceScore: number;
    };
    regionEvidence?: Array<{
      pair: string;
      phylogenetic: boolean;
      movesTogether: boolean;
      sisterTogether: boolean;
      topologyMargin: number;
      treeSourceScore: number;
      bootstrapSupport: number;
      bootstrapReplicates: number;
      bootstrapCutoff: number;
      treeExcluded: boolean;
      correlationR: number;
      correlationP: number;
      correlationInversion: number;
      correlationPermutations?: number[];
      correlationSdmFiltered: boolean;
    }>;
  }>;
}

export interface AnalysisComponentReference {
  originIndex: number;
  kind: "remainder" | "extracted-tract";
  lineage: string[];
  sourceEventId?: string;
  parentLineage?: string[];
  start?: number;
  end?: number;
  wraps?: boolean;
  erasedEventIds: string[];
}

export interface EventComponentProvenance {
  reconstruction: "rdp5-signal-disassembly";
  appliedEventIds: string[];
  recombinant: AnalysisComponentReference;
  majorParent: AnalysisComponentReference;
  minorParent: AnalysisComponentReference;
}

export interface StructuralBreakpointUncertainty {
  source: "rdp5-erased-signal-boundary";
  originalStart: number;
  originalEnd: number;
  originalWraps: boolean;
  piece: number;
  pieces: number;
  uncertainStart: boolean;
  uncertainEnd: boolean;
  adjacentEventIds: string[];
}

export interface RecombinantIdentificationTest {
  id: string;
  label: string;
  sourceRoutine: string;
  direction: "lower" | "higher";
  values: Array<number | null>;
  points: number[];
  fullWeight: number;
  partialWeight: number;
  winnerIndexes: number[];
  decisive: boolean;
}

export interface RecombinantIdentificationOrientation {
  recombinant: number;
  majorParent: number;
  minorParent: number;
  affinitySwitch: number;
  candidateIndex: number;
  sourcePoints: number;
  sourceScore: number;
  sourceShare: number;
}

export interface RecombinantIdentification {
  inference: "rdp5-source-profile-consensus" | "rdp5-source-distance-consensus";
  candidates: number[];
  recommended: number;
  recommendedMajorParent: number;
  recommendedMinorParent: number;
  confidence: number;
  ambiguous: boolean;
  sourceThreshold: number;
  orientations: RecombinantIdentificationOrientation[];
  tests: RecombinantIdentificationTest[];
  cohortSize: number;
  sourceSequenceCount: number;
  sampled: boolean;
  treeEvidence: boolean;
  bootstrapReplicates: number;
  bootstrapCutoff: number;
  quartetCohortSize?: number;
  quartetCounts?: number[];
  dmaxWasmAccelerated?: boolean;
  sourceTieBreak?: string;
  sourceTieBreakValues?: number[] | null;
  implementedComponents: string[];
  pendingComponents: string[];
}

export interface EventAuditEntry {
  id: string;
  timestamp: string;
  action: string;
  summary: string;
}

export interface ProjectAuditEntry extends EventAuditEntry {
  eventId?: string;
  eventSnapshot?: string;
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
  proximityStatistic: number;
  proximityPermutationP: number;
  proximityPermutationReplicates: number;
  ambiguityFraction: number;
  phiPValue?: number;
  phiStatistic?: number;
  phiMean?: number;
  phiVariance?: number;
  phiZ?: number;
  phiInformativeSites?: number;
  phiTotalInformativeSites?: number;
  phiK?: number;
  phiWindow?: number;
  phiSubsampled?: boolean;
  phiValidNormalApproximation?: boolean;
  phiCompatibility?: string;
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
  methodSignals?: MethodSignal[];
  evidence: MethodEvidence[];
  chiSquare: number;
  informativeSites: number;
  decision: EventDecision;
  warnings: string[];
  note: string;
  source: "wasm" | "example" | "manual";
  groupId: string | null;
  ancestralCluster?: AncestralCluster;
  coRecombinantSets?: CoRecombinantSet[];
  componentProvenance?: EventComponentProvenance;
  structuralUncertainty?: StructuralBreakpointUncertainty;
  recombinantIdentification?: RecombinantIdentification;
  alternativeParents?: number[];
  hypothesisTests?: number;
  recalculationNote?: string;
  history: EventAuditEntry[];
  evidenceStale: boolean;
  diagnostics: EventDiagnostics;
}

export interface AnalysisOptions {
  mode: "exploratory" | "query-reference";
  testReferences: boolean;
  circular: boolean;
  window: number;
  rdpWindow: number;
  rdpSignalsPerTriplet: number;
  chiSignalsPerTriplet: number;
  geneconvSignalsPerTriplet: number;
  geneconvGScale: number;
  bootscanWindow: number;
  bootscanStep: number;
  bootscanCutoff: number;
  bootscanSignals: number;
  threeSeqExactOperations: number;
  siskanOutgroupMode: "nearest" | "most-divergent" | "randomized" | "manual";
  siskanOutgroupSequence: number | null;
  siskanPositionMode: "triplet-variable" | "quartet-variable" | "all";
  siskanGapMode: "strip" | "fifth-state";
  siskanWindow: number;
  siskanStep: number;
  siskanScanPermutations: number;
  siskanPValuePermutations: number;
  step: number;
  alpha: number;
  correction: "bonferroni" | "holm" | "none";
  minMethods: number;
  candidateParents: number;
  methods: MethodName[];
  exhaustive: boolean;
  cyclicDetection: boolean;
  maximumDetectionCycles: number;
  polishBreakpoints: boolean;
  burtMode: "rdp5-source" | "manual-step-up";
  burtRandomStarts: number;
  burtMaxIterations: number;
  burtMaxStates: number;
  burtExhaustiveModels: boolean;
  burtPosteriorThreshold: number;
  ancestralClustering: boolean;
  clusterFlankVnps: number;
  clusterMinimumSets: 1 | 2 | 3;
  clusterCorrelationAlpha: number;
  clusterCorrelationR: number;
  clusterSignalOverlap: number;
  clusterTopologyMargin: number;
  clusterBootstrapReplicates: number;
  clusterBootstrapCutoff: number;
  clusterBootstrapBlocks: number;
  clusterTreeTaxaLimit: number;
  clusterMinimumConfidence: number;
  clusterSourceSimilarity: number;
  clusterReciprocal: boolean;
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
      `Sequences must already be aligned. “${mismatch.name}” has ${mismatch.sequence.length.toLocaleString("en-US")} sites; expected ${length.toLocaleString("en-US")}.`,
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
  const alignment = finalizeAlignment("Mosaic virus tutorial · 12 × 2,400 nt", "generated", [
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
  return { ...alignment, createdAt: 1_767_225_600_000 };
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
      statisticLabel: ["identity shift", "fragment score", "bootstrap topology support", "boundary χ²", "boundary χ²", "Sister-Scanning category/sum Z", "maximum HGRW descent"][index],
      calibration: ["binomial", "source finite-G fragment / KA", "seeded p-distance bootstrap + window sign", "χ²", "binary-triplet χ²", "RDP5 vertical permutation Z", "exact HGRW first-passage DP"][index],
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
    const tail = `${(event.start + 1).toLocaleString("en-US")}–${length.toLocaleString("en-US")}`;
    return event.end > 0
      ? `${tail} ↻ 1–${event.end.toLocaleString("en-US")}`
      : `${tail} ↻ origin`;
  }
  return `${(event.start + 1).toLocaleString("en-US")}–${event.end.toLocaleString("en-US")}`;
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
  testReferences: false,
  circular: false,
  window: 120,
  rdpWindow: 30,
  rdpSignalsPerTriplet: 128,
  chiSignalsPerTriplet: 24,
  geneconvSignalsPerTriplet: 64,
  geneconvGScale: 1,
  // RDP5 automatic BOOTSCAN defaults from the supplied source/options form.
  bootscanWindow: 200,
  bootscanStep: 20,
  bootscanCutoff: 0.7,
  bootscanSignals: 20_000,
  // Maximum work for an on-demand Seq3PVals-equivalent calculation. Larger
  // walks follow RDP5 GetTSPVal's SiegmundDiscrete branch.
  threeSeqExactOperations: 1_000_000,
  siskanOutgroupMode: "nearest",
  siskanOutgroupSequence: null,
  siskanPositionMode: "triplet-variable",
  siskanGapMode: "strip",
  siskanWindow: 200,
  siskanStep: 20,
  siskanScanPermutations: 100,
  siskanPValuePermutations: 1000,
  step: 5,
  alpha: 0.05,
  correction: "bonferroni",
  minMethods: 3,
  candidateParents: 8,
  methods: [...SOURCE_READY_METHODS],
  // RDP5 parity default: enumerate every concrete unordered sequence triplet.
  // Parent shortlisting remains an explicit approximate opt-in for previews.
  exhaustive: true,
  // RDP5 manual §4.1.6: repeatedly apply the strongest event, split its
  // recombinant lineage into remainder/tract components, and redo only the
  // triplets touched by that split until no supported signal remains.
  cyclicDetection: true,
  maximumDetectionCycles: 250,
  polishBreakpoints: true,
  burtMode: "rdp5-source",
  burtRandomStarts: 21,
  burtMaxIterations: 100,
  burtMaxStates: 20,
  burtExhaustiveModels: false,
  burtPosteriorThreshold: 0.995,
  ancestralClustering: true,
  clusterFlankVnps: 60,
  clusterMinimumSets: 2,
  clusterCorrelationAlpha: 0.05,
  clusterCorrelationR: 0.83,
  clusterSignalOverlap: 0.3,
  clusterTopologyMargin: 0.005,
  clusterBootstrapReplicates: 100,
  clusterBootstrapCutoff: 0.5,
  clusterBootstrapBlocks: 128,
  clusterTreeTaxaLimit: 32,
  clusterMinimumConfidence: 0.55,
  clusterSourceSimilarity: 0.1,
  clusterReciprocal: false,
  checkMisalignment: true,
  bootstrapReplicates: 100,
  randomSeed: 1511506142,
};

export interface RdpProject {
  schema: "rdp-web/0.5";
  alignment: AlignmentData;
  options: AnalysisOptions;
  events: RdpEvent[];
  metrics: {
    elapsedMs: number;
    comparisons: number;
    engine: string;
    matrixMode?: string;
    tripletMode?: "all-concrete-triplets" | "approximate-parent-shortlist";
    concreteTripletInputs?: boolean;
    parentSamples?: number;
    timing?: { distanceMs: number; scanMs: number; statisticsMs: number; diagnosticsMs?: number; clusteringMs?: number };
    diagnostics?: AlignmentDiagnostics;
    disassembly?: { appliedEvents: number; components: number; erasedCanonicalBases: number };
    rdpSignalTruncations?: number;
    chiSignalTruncations?: number;
    tripletKernelCalls?: { rdp: number; geneconv: number; sourceChi: number; threeSeq: number; siscan: number };
    bootscanSignalTruncations?: number;
    bootscanBatch?: {
      calls: number;
      triplets: number;
      usedPairs: number;
      windows: number;
      replicates: number;
      workspaceBytes: number;
      relationshipMode: "distance";
    };
    geneconvSignalTruncations?: number;
    detectionCycle?: {
      enabled: boolean;
      eventsApplied: number;
      passes: number;
      initialComparisons: number;
      redoComparisons: number;
      stoppedBecause: "no-detectable-signals" | "cycle-cap";
      maximumCycles: number;
    };
  } | null;
  distance: number[];
  auditLog: ProjectAuditEntry[];
}

function finiteNumber(value: unknown, fallback: number): number {
  return typeof value === "number" && Number.isFinite(value) ? value : fallback;
}

export function serializeProject(project: Omit<RdpProject, "schema" | "auditLog"> & { auditLog?: ProjectAuditEntry[] }): string {
  return JSON.stringify({ schema: "rdp-web/0.5", ...project, auditLog: project.auditLog ?? [] }, null, 2);
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
        referenceGroup: typeof candidate.referenceGroup === "string" && candidate.referenceGroup.trim() ? candidate.referenceGroup.trim() : undefined,
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
    ? rawOptions.methods.filter((method): method is MethodName => SOURCE_READY_METHODS.includes(method as MethodName))
    : [...SOURCE_READY_METHODS];
  const options: AnalysisOptions = {
    ...DEFAULT_OPTIONS,
    ...rawOptions,
    methods: methods.length ? methods : [...SOURCE_READY_METHODS],
    mode: rawOptions.mode === "query-reference" ? "query-reference" : "exploratory",
    testReferences: rawOptions.testReferences === true,
    correction: rawOptions.correction === "holm" || rawOptions.correction === "none" ? rawOptions.correction : "bonferroni",
    window: Math.max(12, Math.min(Math.max(12, alignment.length), Math.trunc(finiteNumber(rawOptions.window, DEFAULT_OPTIONS.window)))),
    step: Math.max(1, Math.min(Math.max(1, alignment.length), Math.trunc(finiteNumber(rawOptions.step, DEFAULT_OPTIONS.step)))),
    alpha: Math.max(1e-12, Math.min(1, finiteNumber(rawOptions.alpha, DEFAULT_OPTIONS.alpha))),
    minMethods: Math.max(1, Math.min(Math.max(1, methods.length), Math.trunc(finiteNumber(rawOptions.minMethods, DEFAULT_OPTIONS.minMethods)))),
    candidateParents: Math.max(2, Math.min(300, Math.trunc(finiteNumber(rawOptions.candidateParents, DEFAULT_OPTIONS.candidateParents)))),
    cyclicDetection: rawOptions.cyclicDetection !== false,
    maximumDetectionCycles: Math.max(1, Math.min(1000, Math.trunc(finiteNumber(rawOptions.maximumDetectionCycles, DEFAULT_OPTIONS.maximumDetectionCycles)))),
    rdpWindow: Math.max(5, Math.min(300, Math.trunc(finiteNumber(rawOptions.rdpWindow, DEFAULT_OPTIONS.rdpWindow)))),
    rdpSignalsPerTriplet: Math.max(1, Math.min(256, Math.trunc(finiteNumber(rawOptions.rdpSignalsPerTriplet, DEFAULT_OPTIONS.rdpSignalsPerTriplet)))),
    chiSignalsPerTriplet: Math.max(1, Math.min(256, Math.trunc(finiteNumber(rawOptions.chiSignalsPerTriplet, DEFAULT_OPTIONS.chiSignalsPerTriplet)))),
    geneconvSignalsPerTriplet: Math.max(1, Math.min(256, Math.trunc(finiteNumber(rawOptions.geneconvSignalsPerTriplet, DEFAULT_OPTIONS.geneconvSignalsPerTriplet)))),
    geneconvGScale: Math.max(0, Math.min(100, finiteNumber(rawOptions.geneconvGScale, DEFAULT_OPTIONS.geneconvGScale))),
    bootscanWindow: Math.max(5, Math.min(Math.max(5, Math.floor(alignment.length / 2)), Math.trunc(finiteNumber(rawOptions.bootscanWindow, DEFAULT_OPTIONS.bootscanWindow)))),
    bootscanStep: Math.max(1, Math.min(Math.max(1, Math.floor(alignment.length / 4)), Math.trunc(finiteNumber(rawOptions.bootscanStep, DEFAULT_OPTIONS.bootscanStep)))),
    bootscanCutoff: Math.max(0.5, Math.min(0.999, finiteNumber(rawOptions.bootscanCutoff, DEFAULT_OPTIONS.bootscanCutoff))),
    bootscanSignals: Math.max(128, Math.min(50_000, Math.trunc(finiteNumber(rawOptions.bootscanSignals, DEFAULT_OPTIONS.bootscanSignals)))),
    threeSeqExactOperations: Math.max(10_000, Math.min(20_000_000, Math.trunc(finiteNumber(rawOptions.threeSeqExactOperations, DEFAULT_OPTIONS.threeSeqExactOperations)))),
    siskanOutgroupMode: rawOptions.siskanOutgroupMode === "most-divergent" || rawOptions.siskanOutgroupMode === "randomized" || rawOptions.siskanOutgroupMode === "manual"
      ? rawOptions.siskanOutgroupMode
      : "nearest",
    siskanOutgroupSequence: Number.isFinite(rawOptions.siskanOutgroupSequence)
      && (rawOptions.siskanOutgroupSequence as number) >= 0
      && (rawOptions.siskanOutgroupSequence as number) < alignment.sequences.length
      ? Math.trunc(rawOptions.siskanOutgroupSequence as number)
      : null,
    siskanPositionMode: rawOptions.siskanPositionMode === "quartet-variable" || rawOptions.siskanPositionMode === "all"
      ? rawOptions.siskanPositionMode
      : "triplet-variable",
    siskanGapMode: rawOptions.siskanGapMode === "fifth-state" ? "fifth-state" : "strip",
    siskanWindow: Math.max(12, Math.min(Math.max(12, alignment.length), Math.trunc(finiteNumber(rawOptions.siskanWindow, DEFAULT_OPTIONS.siskanWindow)))),
    siskanStep: Math.max(1, Math.min(Math.max(1, alignment.length), Math.trunc(finiteNumber(rawOptions.siskanStep, DEFAULT_OPTIONS.siskanStep)))),
    siskanScanPermutations: Math.max(2, Math.min(1000, Math.trunc(finiteNumber(rawOptions.siskanScanPermutations, DEFAULT_OPTIONS.siskanScanPermutations)))),
    siskanPValuePermutations: Math.max(
      Math.max(2, Math.min(1000, Math.trunc(finiteNumber(rawOptions.siskanScanPermutations, DEFAULT_OPTIONS.siskanScanPermutations)))),
      Math.min(10_000, Math.trunc(finiteNumber(rawOptions.siskanPValuePermutations, DEFAULT_OPTIONS.siskanPValuePermutations))),
    ),
    bootstrapReplicates: Math.max(0, Math.min(1000, Math.trunc(finiteNumber(rawOptions.bootstrapReplicates, DEFAULT_OPTIONS.bootstrapReplicates)))),
    randomSeed: Math.trunc(finiteNumber(rawOptions.randomSeed, DEFAULT_OPTIONS.randomSeed)) >>> 0,
    burtMode: rawOptions.burtMode === "manual-step-up" ? "manual-step-up" : "rdp5-source",
    burtRandomStarts: Math.max(1, Math.min(64, Math.trunc(finiteNumber(rawOptions.burtRandomStarts, DEFAULT_OPTIONS.burtRandomStarts)))),
    burtMaxIterations: Math.max(2, Math.min(250, Math.trunc(finiteNumber(rawOptions.burtMaxIterations, DEFAULT_OPTIONS.burtMaxIterations)))),
    burtMaxStates: Math.max(2, Math.min(20, Math.trunc(finiteNumber(rawOptions.burtMaxStates, DEFAULT_OPTIONS.burtMaxStates)))),
    burtExhaustiveModels: rawOptions.burtExhaustiveModels === true,
    burtPosteriorThreshold: Math.max(0.5, Math.min(0.9999, finiteNumber(rawOptions.burtPosteriorThreshold, DEFAULT_OPTIONS.burtPosteriorThreshold))),
    ancestralClustering: rawOptions.ancestralClustering !== false,
    clusterFlankVnps: Math.max(4, Math.min(200, Math.trunc(finiteNumber(rawOptions.clusterFlankVnps, DEFAULT_OPTIONS.clusterFlankVnps)))),
    clusterMinimumSets: ([1, 2, 3].includes(Math.trunc(finiteNumber(rawOptions.clusterMinimumSets, 2))) ? Math.trunc(finiteNumber(rawOptions.clusterMinimumSets, 2)) : 2) as 1 | 2 | 3,
    clusterCorrelationAlpha: Math.max(1e-6, Math.min(0.5, finiteNumber(rawOptions.clusterCorrelationAlpha, DEFAULT_OPTIONS.clusterCorrelationAlpha))),
    clusterCorrelationR: Math.max(0, Math.min(0.999999, finiteNumber(rawOptions.clusterCorrelationR, DEFAULT_OPTIONS.clusterCorrelationR))),
    clusterSignalOverlap: Math.max(0.05, Math.min(1, finiteNumber(rawOptions.clusterSignalOverlap, DEFAULT_OPTIONS.clusterSignalOverlap))),
    clusterTopologyMargin: Math.max(0, Math.min(1, finiteNumber(rawOptions.clusterTopologyMargin, DEFAULT_OPTIONS.clusterTopologyMargin))),
    clusterBootstrapReplicates: Math.max(0, Math.min(1000, Math.trunc(finiteNumber(rawOptions.clusterBootstrapReplicates, DEFAULT_OPTIONS.clusterBootstrapReplicates)))),
    clusterBootstrapCutoff: Math.max(0, Math.min(1, finiteNumber(rawOptions.clusterBootstrapCutoff, DEFAULT_OPTIONS.clusterBootstrapCutoff))),
    clusterBootstrapBlocks: Math.max(8, Math.min(4096, Math.trunc(finiteNumber(rawOptions.clusterBootstrapBlocks, DEFAULT_OPTIONS.clusterBootstrapBlocks)))),
    clusterTreeTaxaLimit: Math.max(4, Math.min(300, Math.trunc(finiteNumber(rawOptions.clusterTreeTaxaLimit, DEFAULT_OPTIONS.clusterTreeTaxaLimit)))),
    clusterMinimumConfidence: Math.max(0, Math.min(1, finiteNumber(rawOptions.clusterMinimumConfidence, DEFAULT_OPTIONS.clusterMinimumConfidence))),
    clusterSourceSimilarity: Math.max(0, Math.min(1, finiteNumber(rawOptions.clusterSourceSimilarity, DEFAULT_OPTIONS.clusterSourceSimilarity))),
    clusterReciprocal: rawOptions.clusterReciprocal === true,
  };
  const rawEvents = Array.isArray(raw.events) ? raw.events : [];
  const parseComponentReference = (value: unknown): AnalysisComponentReference | undefined => {
    if (!value || typeof value !== "object") return undefined;
    const candidate = value as Partial<AnalysisComponentReference>;
    const originIndex = Math.trunc(finiteNumber(candidate.originIndex, -1));
    if (originIndex < 0 || originIndex >= alignment.sequences.length) return undefined;
    const kind = candidate.kind === "extracted-tract" ? "extracted-tract" : candidate.kind === "remainder" ? "remainder" : null;
    if (!kind) return undefined;
    const lineage = Array.isArray(candidate.lineage)
      ? candidate.lineage.filter((item): item is string => typeof item === "string" && item.length > 0)
      : [];
    const parentLineage = Array.isArray(candidate.parentLineage)
      ? candidate.parentLineage.filter((item): item is string => typeof item === "string" && item.length > 0)
      : undefined;
    return {
      originIndex,
      kind,
      lineage,
      sourceEventId: typeof candidate.sourceEventId === "string" ? candidate.sourceEventId : undefined,
      parentLineage,
      start: Number.isFinite(candidate.start) ? Math.max(0, Math.min(alignment.length - 1, Math.trunc(candidate.start as number))) : undefined,
      end: Number.isFinite(candidate.end) ? Math.max(0, Math.min(alignment.length, Math.trunc(candidate.end as number))) : undefined,
      wraps: candidate.wraps === true,
      erasedEventIds: Array.isArray(candidate.erasedEventIds)
        ? candidate.erasedEventIds.filter((item): item is string => typeof item === "string" && item.length > 0)
        : [],
    };
  };
  const parseRecombinantIdentification = (value: unknown): RecombinantIdentification | undefined => {
    if (!value || typeof value !== "object") return undefined;
    const candidate = value as Partial<RecombinantIdentification>;
    const sequenceIndex = (item: unknown): number | null => {
      const index = Math.trunc(finiteNumber(item, -1));
      return index >= 0 && index < alignment.sequences.length ? index : null;
    };
    const candidates = Array.isArray(candidate.candidates)
      ? candidate.candidates.flatMap((item) => {
          const index = sequenceIndex(item);
          return index === null ? [] : [index];
        })
      : [];
    const recommended = sequenceIndex(candidate.recommended);
    const recommendedMajorParent = sequenceIndex(candidate.recommendedMajorParent);
    const recommendedMinorParent = sequenceIndex(candidate.recommendedMinorParent);
    if (candidates.length !== 3 || new Set(candidates).size !== 3 || recommended === null || recommendedMajorParent === null || recommendedMinorParent === null) return undefined;
    const orientations = Array.isArray(candidate.orientations) ? candidate.orientations.flatMap((item): RecombinantIdentificationOrientation[] => {
      if (!item || typeof item !== "object") return [];
      const recombinant = sequenceIndex(item.recombinant);
      const majorParent = sequenceIndex(item.majorParent);
      const minorParent = sequenceIndex(item.minorParent);
      if (recombinant === null || majorParent === null || minorParent === null || new Set([recombinant, majorParent, minorParent]).size !== 3) return [];
      return [{
        recombinant,
        majorParent,
        minorParent,
        affinitySwitch: finiteNumber(item.affinitySwitch, 0),
        candidateIndex: Math.max(0, Math.min(2, Math.trunc(finiteNumber(item.candidateIndex, 0)))),
        sourcePoints: Math.max(0, finiteNumber(item.sourcePoints, 0)),
        sourceScore: Math.max(0, Math.min(100, finiteNumber(item.sourceScore, 0))),
        sourceShare: Math.max(0, Math.min(1, finiteNumber(item.sourceShare, 0))),
      }];
    }) : [];
    if (orientations.length !== 3) return undefined;
    const tests = Array.isArray(candidate.tests) ? candidate.tests.flatMap((item): RecombinantIdentificationTest[] => {
      if (!item || typeof item !== "object" || typeof item.id !== "string") return [];
      const values = Array.isArray(item.values) ? item.values.slice(0, 3).map((entry) => typeof entry === "number" && Number.isFinite(entry) ? entry : null) : [];
      const points = Array.isArray(item.points) ? item.points.slice(0, 3).map((entry) => finiteNumber(entry, 0)) : [];
      if (values.length !== 3 || points.length !== 3) return [];
      return [{
        id: item.id,
        label: typeof item.label === "string" ? item.label : item.id,
        sourceRoutine: typeof item.sourceRoutine === "string" ? item.sourceRoutine : "imported source test",
        direction: item.direction === "higher" ? "higher" : "lower",
        values,
        points,
        fullWeight: Math.max(0, finiteNumber(item.fullWeight, 0)),
        partialWeight: Math.max(0, finiteNumber(item.partialWeight, 0)),
        winnerIndexes: Array.isArray(item.winnerIndexes) ? item.winnerIndexes.map((entry) => Math.trunc(finiteNumber(entry, -1))).filter((entry) => entry >= 0 && entry < 3) : [],
        decisive: item.decisive === true,
      }];
    }) : [];
    return {
      inference: candidate.inference === "rdp5-source-profile-consensus" ? "rdp5-source-profile-consensus" : "rdp5-source-distance-consensus",
      candidates,
      recommended,
      recommendedMajorParent,
      recommendedMinorParent,
      confidence: Math.max(0, Math.min(1, finiteNumber(candidate.confidence, 1 / 3))),
      ambiguous: candidate.ambiguous === true,
      sourceThreshold: Math.max(0, Math.min(1, finiteNumber(candidate.sourceThreshold, 0.6))),
      orientations,
      tests,
      cohortSize: Math.max(3, Math.trunc(finiteNumber(candidate.cohortSize, 3))),
      sourceSequenceCount: Math.max(3, Math.trunc(finiteNumber(candidate.sourceSequenceCount, alignment.sequences.length))),
      sampled: candidate.sampled === true,
      treeEvidence: candidate.treeEvidence === true,
      bootstrapReplicates: Math.max(0, Math.trunc(finiteNumber(candidate.bootstrapReplicates, 0))),
      bootstrapCutoff: Math.max(0, Math.min(1, finiteNumber(candidate.bootstrapCutoff, 0.5))),
      quartetCohortSize: Math.max(0, Math.trunc(finiteNumber(candidate.quartetCohortSize, 0))),
      quartetCounts: Array.isArray(candidate.quartetCounts) ? candidate.quartetCounts.slice(0, 3).map((item) => Math.max(0, Math.trunc(finiteNumber(item, 0)))) : [],
      dmaxWasmAccelerated: candidate.dmaxWasmAccelerated === true,
      sourceTieBreak: typeof candidate.sourceTieBreak === "string" ? candidate.sourceTieBreak : undefined,
      sourceTieBreakValues: candidate.sourceTieBreakValues === null
        ? null
        : Array.isArray(candidate.sourceTieBreakValues)
          ? candidate.sourceTieBreakValues.slice(0, 3).map((item) => finiteNumber(item, 0))
          : undefined,
      implementedComponents: Array.isArray(candidate.implementedComponents) ? candidate.implementedComponents.filter((item): item is string => typeof item === "string") : [],
      pendingComponents: Array.isArray(candidate.pendingComponents) ? candidate.pendingComponents.filter((item): item is string => typeof item === "string") : [],
    };
  };
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
        correctionScope: typeof item.correctionScope === "string" ? item.correctionScope : undefined,
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
            method: event.breakpointModel.method === "burt-hmm" || event.breakpointModel.method === "two-state-hmm" || event.breakpointModel.method === "manual"
              ? event.breakpointModel.method
              : "local-chi-square",
            informativeSites: Math.trunc(finiteNumber(event.breakpointModel.informativeSites, 0)),
            stateSwitches: Math.trunc(finiteNumber(event.breakpointModel.stateSwitches, 0)) || undefined,
            majorFit: finiteNumber(event.breakpointModel.majorFit, 0) || undefined,
            minorFit: finiteNumber(event.breakpointModel.minorFit, 0) || undefined,
            states: Math.trunc(finiteNumber(event.breakpointModel.states, 0)) || undefined,
            logLikelihood: Number.isFinite(event.breakpointModel.logLikelihood) ? event.breakpointModel.logLikelihood : undefined,
            viterbiLogLikelihood: Number.isFinite(event.breakpointModel.viterbiLogLikelihood) ? event.breakpointModel.viterbiLogLikelihood : undefined,
            bic: Number.isFinite(event.breakpointModel.bic) ? event.breakpointModel.bic : undefined,
            aic: Number.isFinite(event.breakpointModel.aic) ? event.breakpointModel.aic : undefined,
            criterion: typeof event.breakpointModel.criterion === "string" ? event.breakpointModel.criterion : undefined,
            randomStarts: Math.trunc(finiteNumber(event.breakpointModel.randomStarts, 0)) || undefined,
            iterations: Math.trunc(finiteNumber(event.breakpointModel.iterations, 0)) || undefined,
            winningRestart: Math.trunc(finiteNumber(event.breakpointModel.winningRestart, 0)) || undefined,
            selectedState: Number.isFinite(event.breakpointModel.selectedState) ? Math.trunc(event.breakpointModel.selectedState as number) : undefined,
            posteriorThreshold: finiteNumber(event.breakpointModel.posteriorThreshold, 0) || undefined,
            sourceParity: event.breakpointModel.sourceParity === true,
            sourceCompatibility: typeof event.breakpointModel.sourceCompatibility === "string" ? event.breakpointModel.sourceCompatibility : undefined,
            sourceRoutines: Array.isArray(event.breakpointModel.sourceRoutines) ? event.breakpointModel.sourceRoutines.filter((entry): entry is string => typeof entry === "string") : undefined,
            sequenceOrder: Array.isArray(event.breakpointModel.sequenceOrder) ? event.breakpointModel.sequenceOrder.map((entry) => Math.trunc(finiteNumber(entry, 0))) : undefined,
            stateDominantCategories: Array.isArray(event.breakpointModel.stateDominantCategories) ? event.breakpointModel.stateDominantCategories.map((entry) => Math.trunc(finiteNumber(entry, 0))) : undefined,
            circularPadding: event.breakpointModel.circularPadding && typeof event.breakpointModel.circularPadding === "object"
              ? {
                  offset: Math.trunc(finiteNumber(event.breakpointModel.circularPadding.offset, 0)),
                  fittedSites: Math.trunc(finiteNumber(event.breakpointModel.circularPadding.fittedSites, 0)),
                  croppedSites: Math.trunc(finiteNumber(event.breakpointModel.circularPadding.croppedSites, 0)),
                }
              : undefined,
            candidateBreakpoints: Array.isArray(event.breakpointModel.candidateBreakpoints) && event.breakpointModel.candidateBreakpoints.length === 2
              ? [finiteNumber(event.breakpointModel.candidateBreakpoints[0], start), finiteNumber(event.breakpointModel.candidateBreakpoints[1], end)]
              : undefined,
            polishedBreakpoints: Array.isArray(event.breakpointModel.polishedBreakpoints) && event.breakpointModel.polishedBreakpoints.length === 2
              ? [finiteNumber(event.breakpointModel.polishedBreakpoints[0], start), finiteNumber(event.breakpointModel.polishedBreakpoints[1], end)]
              : undefined,
            polishDecision: event.breakpointModel.polishDecision && typeof event.breakpointModel.polishDecision === "object"
              ? {
                  startAdopted: event.breakpointModel.polishDecision.startAdopted === true,
                  endAdopted: event.breakpointModel.polishDecision.endAdopted === true,
                  sameSwitchResolved: event.breakpointModel.polishDecision.sameSwitchResolved === true,
                  startMissingBoundary: event.breakpointModel.polishDecision.startMissingBoundary === true,
                  endMissingBoundary: event.breakpointModel.polishDecision.endMissingBoundary === true,
                  revertedForInformation: event.breakpointModel.polishDecision.revertedForInformation === true,
                  insideVariableSites: Math.trunc(finiteNumber(event.breakpointModel.polishDecision.insideVariableSites, 0)),
                  outsideVariableSites: Math.trunc(finiteNumber(event.breakpointModel.polishDecision.outsideVariableSites, 0)),
                  startWithin99: event.breakpointModel.polishDecision.startWithin99 === true,
                  endWithin99: event.breakpointModel.polishDecision.endWithin99 === true,
                  startVariableSiteDistance: Number.isFinite(event.breakpointModel.polishDecision.startVariableSiteDistance) ? event.breakpointModel.polishDecision.startVariableSiteDistance : undefined,
                  endVariableSiteDistance: Number.isFinite(event.breakpointModel.polishDecision.endVariableSiteDistance) ? event.breakpointModel.polishDecision.endVariableSiteDistance : undefined,
                }
              : undefined,
            confidence99Start: Array.isArray(event.breakpointModel.confidence99Start) && event.breakpointModel.confidence99Start.length === 2
              ? [finiteNumber(event.breakpointModel.confidence99Start[0], 0), finiteNumber(event.breakpointModel.confidence99Start[1], alignment.length)]
              : undefined,
            confidence99End: Array.isArray(event.breakpointModel.confidence99End) && event.breakpointModel.confidence99End.length === 2
              ? [finiteNumber(event.breakpointModel.confidence99End[0], 0), finiteNumber(event.breakpointModel.confidence99End[1], alignment.length)]
              : undefined,
            emissions: Array.isArray(event.breakpointModel.emissions)
              ? event.breakpointModel.emissions.map((row) => Array.isArray(row) ? row.map((item) => finiteNumber(item, 0)) : [])
              : undefined,
            transitions: Array.isArray(event.breakpointModel.transitions)
              ? event.breakpointModel.transitions.map((row) => Array.isArray(row) ? row.map((item) => finiteNumber(item, 0)) : [])
              : undefined,
            switches: Array.isArray(event.breakpointModel.switches) ? event.breakpointModel.switches.flatMap((entry) => {
              if (!entry || typeof entry !== "object" || !Array.isArray(entry.confidence95) || entry.confidence95.length !== 2) return [];
              return [{
                position: Math.trunc(finiteNumber(entry.position, 0)),
                informativeIndex: Number.isFinite(entry.informativeIndex) ? Math.trunc(entry.informativeIndex as number) : undefined,
                fromState: Math.trunc(finiteNumber(entry.fromState, 0)),
                toState: Math.trunc(finiteNumber(entry.toState, 0)),
                confidence95: [finiteNumber(entry.confidence95[0], 0), finiteNumber(entry.confidence95[1], alignment.length)] as [number, number],
                confidence99: Array.isArray(entry.confidence99) && entry.confidence99.length === 2
                  ? [finiteNumber(entry.confidence99[0], 0), finiteNumber(entry.confidence99[1], alignment.length)] as [number, number]
                  : undefined,
                sourceCoordinates: Array.isArray(entry.sourceCoordinates) ? entry.sourceCoordinates.map((value) => finiteNumber(value, 0)) : undefined,
                matchedStart: entry.matchedStart === true,
                matchedEnd: entry.matchedEnd === true,
              }];
            }) : undefined,
            posteriorTrace: Array.isArray(event.breakpointModel.posteriorTrace) ? event.breakpointModel.posteriorTrace.flatMap((entry) => {
              if (!entry || typeof entry !== "object" || !Array.isArray(entry.probabilities)) return [];
              return [{
                position: Math.trunc(finiteNumber(entry.position, 0)),
                informativeIndex: Number.isFinite(entry.informativeIndex) ? Math.trunc(entry.informativeIndex as number) : undefined,
                state: Math.trunc(finiteNumber(entry.state, 0)),
                probabilities: entry.probabilities.map((item) => finiteNumber(item, 0)),
              }];
            }) : undefined,
            modelSelection: Array.isArray(event.breakpointModel.modelSelection) ? event.breakpointModel.modelSelection.flatMap((entry) => {
              if (!entry || typeof entry !== "object") return [];
              return [{
                states: Math.trunc(finiteNumber(entry.states, 0)),
                logLikelihood: finiteNumber(entry.logLikelihood, 0),
                bic: finiteNumber(entry.bic, 0),
                aic: finiteNumber(entry.aic, 0),
                iterations: Math.trunc(finiteNumber(entry.iterations, 0)),
                winningRestart: Math.trunc(finiteNumber(entry.winningRestart, 0)),
              }];
            }) : undefined,
          }
        : undefined,
      methodSignals: Array.isArray(event.methodSignals) ? event.methodSignals.flatMap((signal) => {
        if (!signal || typeof signal !== "object") return [];
        const method = signal.method === "shared-screen" || PRIMARY_METHODS.includes(signal.method as MethodName)
          ? signal.method as MethodName | "shared-screen"
          : null;
        if (!method) return [];
        const signalStart = Math.max(0, Math.min(alignment.length - 1, Math.trunc(finiteNumber(signal.start, start))));
        const signalEnd = Math.max(0, Math.min(alignment.length, Math.trunc(finiteNumber(signal.end, end))));
        return [{
          method,
          start: signalStart,
          end: signalEnd,
          wraps: signal.wraps === true && signalStart > signalEnd,
          statistic: finiteNumber(signal.statistic, 0),
          locator: typeof signal.locator === "string" ? signal.locator : "imported locator",
          sourceRoutine: typeof signal.sourceRoutine === "string" ? signal.sourceRoutine : undefined,
          sourceChi: signal.sourceChi && typeof signal.sourceChi === "object" ? {
            track: Math.max(0, Math.min(2, Math.trunc(finiteNumber(signal.sourceChi.track, 0)))),
            targetSlot: signal.sourceChi.targetSlot === null
              ? null
              : Math.max(0, Math.min(2, Math.trunc(finiteNumber(signal.sourceChi.targetSlot, 0)))),
            informativeSites: Math.max(0, Math.trunc(finiteNumber(signal.sourceChi.informativeSites, 0))),
            halfWindow: Math.max(0, Math.trunc(finiteNumber(signal.sourceChi.halfWindow, 0))),
            boundaryStatistics: [
              Math.max(0, finiteNumber(signal.sourceChi.boundaryStatistics?.[0], 0)),
              Math.max(0, finiteNumber(signal.sourceChi.boundaryStatistics?.[1], 0)),
            ] as [number, number],
            boundaryRanks: [
              Math.max(0, Math.trunc(finiteNumber(signal.sourceChi.boundaryRanks?.[0], 0))),
              Math.max(0, Math.trunc(finiteNumber(signal.sourceChi.boundaryRanks?.[1], 0))),
            ] as [number, number],
            growthWidths: [
              Math.max(0, Math.trunc(finiteNumber(signal.sourceChi.growthWidths?.[0], 0))),
              Math.max(0, Math.trunc(finiteNumber(signal.sourceChi.growthWidths?.[1], 0))),
            ] as [number, number],
            direction: signal.sourceChi.direction === -1 ? -1 as const : 1 as const,
          } : undefined,
          sourceGeneconv: signal.sourceGeneconv && typeof signal.sourceGeneconv === "object" ? {
            track: Math.max(0, Math.min(5, Math.trunc(finiteNumber(signal.sourceGeneconv.track, 0)))),
            targetSlot: Math.max(0, Math.min(2, Math.trunc(finiteNumber(signal.sourceGeneconv.targetSlot, 0)))),
            minorSlot: Math.max(0, Math.min(2, Math.trunc(finiteNumber(signal.sourceGeneconv.minorSlot, 1)))),
            majorSlot: Math.max(0, Math.min(2, Math.trunc(finiteNumber(signal.sourceGeneconv.majorSlot, 2)))),
            fragmentScore: Math.max(0, Math.trunc(finiteNumber(signal.sourceGeneconv.fragmentScore, signal.statistic))),
            informativeSites: Math.max(0, Math.trunc(finiteNumber(signal.sourceGeneconv.informativeSites, 0))),
            matchingSites: Math.max(0, Math.trunc(finiteNumber(signal.sourceGeneconv.matchingSites, 0))),
            mismatchSites: Math.max(0, Math.trunc(finiteNumber(signal.sourceGeneconv.mismatchSites, 0))),
            mismatchPenalty: Math.max(1, Math.trunc(finiteNumber(signal.sourceGeneconv.mismatchPenalty, 1))),
            rawP: Math.max(Number.MIN_VALUE, Math.min(1, finiteNumber(signal.sourceGeneconv.rawP, 1))),
            startRank: Math.max(0, Math.trunc(finiteNumber(signal.sourceGeneconv.startRank, 0))),
            endRank: Math.max(0, Math.trunc(finiteNumber(signal.sourceGeneconv.endRank, 0))),
          } : undefined,
          sourceBootscan: signal.sourceBootscan && typeof signal.sourceBootscan === "object" ? {
            topology: Math.max(0, Math.min(2, Math.trunc(finiteNumber(signal.sourceBootscan.topology, 0)))),
            baselineTopology: Math.max(0, Math.min(2, Math.trunc(finiteNumber(signal.sourceBootscan.baselineTopology, 0)))),
            bootstrapSupport: Math.max(0, Math.min(1, finiteNumber(signal.sourceBootscan.bootstrapSupport, 0))),
            bootstrapReplicates: Math.max(2, Math.min(1000, Math.trunc(finiteNumber(signal.sourceBootscan.bootstrapReplicates, 100)))),
            runWindows: Math.max(1, Math.trunc(finiteNumber(signal.sourceBootscan.runWindows, 1))),
            tractPairMatches: Math.max(0, Math.trunc(finiteNumber(signal.sourceBootscan.tractPairMatches, 0))),
            backgroundPairMatches: Math.max(0, Math.trunc(finiteNumber(signal.sourceBootscan.backgroundPairMatches, 0))),
            tractInformativeSites: Math.max(0, Math.trunc(finiteNumber(signal.sourceBootscan.tractInformativeSites, 0))),
            informativeSites: Math.max(0, Math.trunc(finiteNumber(signal.sourceBootscan.informativeSites, 0))),
            rawP: Math.max(Number.MIN_VALUE, Math.min(1, finiteNumber(signal.sourceBootscan.rawP, 1))),
            window: Math.max(5, Math.trunc(finiteNumber(signal.sourceBootscan.window, DEFAULT_OPTIONS.bootscanWindow))),
            step: Math.max(1, Math.trunc(finiteNumber(signal.sourceBootscan.step, DEFAULT_OPTIONS.bootscanStep))),
            relationshipMode: "distance" as const,
          } : undefined,
          sourceThreeSeq: signal.sourceThreeSeq && typeof signal.sourceThreeSeq === "object" ? {
            target: Math.max(0, Math.min(alignment.sequences.length - 1, Math.trunc(finiteNumber(signal.sourceThreeSeq.target, finiteNumber(event.recombinant, 0))))),
            majorParent: Math.max(0, Math.min(alignment.sequences.length - 1, Math.trunc(finiteNumber(signal.sourceThreeSeq.majorParent, finiteNumber(event.majorParent, 0))))),
            minorParent: Math.max(0, Math.min(alignment.sequences.length - 1, Math.trunc(finiteNumber(signal.sourceThreeSeq.minorParent, finiteNumber(event.minorParent, 0))))),
            direction: signal.sourceThreeSeq.direction === -1 ? -1 as const : 1 as const,
            upSteps: Math.max(0, Math.trunc(finiteNumber(signal.sourceThreeSeq.upSteps, 0))),
            downSteps: Math.max(0, Math.trunc(finiteNumber(signal.sourceThreeSeq.downSteps, 0))),
            descent: Math.max(0, Math.trunc(finiteNumber(signal.sourceThreeSeq.descent, signal.statistic))),
            informativeSites: Math.max(0, Math.trunc(finiteNumber(signal.sourceThreeSeq.informativeSites, 0))),
            cycle: Math.max(0, Math.min(2, Math.trunc(finiteNumber(signal.sourceThreeSeq.cycle, 0)))),
            rawStart: Math.max(0, Math.min(alignment.length - 1, Math.trunc(finiteNumber(signal.sourceThreeSeq.rawStart, signalStart)))),
            rawEnd: Math.max(0, Math.min(alignment.length, Math.trunc(finiteNumber(signal.sourceThreeSeq.rawEnd, signalEnd)))),
            rawP: Math.max(Number.MIN_VALUE, Math.min(1, finiteNumber(signal.sourceThreeSeq.rawP, 1))),
            probabilityMode: signal.sourceThreeSeq.probabilityMode === "siegmund-discrete"
              || signal.sourceThreeSeq.probabilityMode === "scaled-table"
              || signal.sourceThreeSeq.probabilityMode === "unavailable"
              ? signal.sourceThreeSeq.probabilityMode
              : "exact-table" as const,
            sourceWrap: signal.sourceThreeSeq.sourceWrap === true,
            linearComplement: signal.sourceThreeSeq.linearComplement === true,
            splitRefined: signal.sourceThreeSeq.splitRefined === true,
            fullDescent: Number.isFinite(signal.sourceThreeSeq.fullDescent)
              ? Math.max(0, Math.trunc(signal.sourceThreeSeq.fullDescent as number))
              : undefined,
            splitInformativeSites: Number.isFinite(signal.sourceThreeSeq.splitInformativeSites)
              ? Math.max(0, Math.trunc(signal.sourceThreeSeq.splitInformativeSites as number))
              : undefined,
          } : undefined,
          sourceSiScan: signal.sourceSiScan && typeof signal.sourceSiScan === "object"
            && Array.isArray(signal.sourceSiScan.topologyTriplet)
            && signal.sourceSiScan.topologyTriplet.length === 3 ? {
              rawP: Math.max(Number.MIN_VALUE, Math.min(1, finiteNumber(signal.sourceSiScan.rawP, 1))),
              rawStart: Math.max(0, Math.min(alignment.length - 1, Math.trunc(finiteNumber(signal.sourceSiScan.rawStart, signalStart)))),
              rawEnd: Math.max(0, Math.min(alignment.length, Math.trunc(finiteNumber(signal.sourceSiScan.rawEnd, signalEnd)))),
              runWindows: Math.max(1, Math.trunc(finiteNumber(signal.sourceSiScan.runWindows, 1))),
              outgroupSourcePath: typeof signal.sourceSiScan.outgroupSourcePath === "string"
                ? signal.sourceSiScan.outgroupSourcePath
                : "imported source outgroup rule",
              positionMode: signal.sourceSiScan.positionMode === "quartet-variable" || signal.sourceSiScan.positionMode === "all"
                ? signal.sourceSiScan.positionMode
                : "triplet-variable" as const,
              gapMode: signal.sourceSiScan.gapMode === "fifth-state" ? "fifth-state" as const : "strip" as const,
              window: Math.max(12, Math.min(alignment.length, Math.trunc(finiteNumber(signal.sourceSiScan.window, DEFAULT_OPTIONS.siskanWindow)))),
              step: Math.max(1, Math.min(alignment.length, Math.trunc(finiteNumber(signal.sourceSiScan.step, DEFAULT_OPTIONS.siskanStep)))),
              topologyTriplet: signal.sourceSiScan.topologyTriplet.map((value) => (
                Math.max(0, Math.min(alignment.sequences.length - 1, Math.trunc(finiteNumber(value, 0))))
              )) as [number, number, number],
              recombinant: Math.max(0, Math.min(alignment.sequences.length - 1, Math.trunc(finiteNumber(signal.sourceSiScan.recombinant, finiteNumber(event.recombinant, 0))))),
              majorParent: Math.max(0, Math.min(alignment.sequences.length - 1, Math.trunc(finiteNumber(signal.sourceSiScan.majorParent, finiteNumber(event.majorParent, 0))))),
              minorParent: Math.max(0, Math.min(alignment.sequences.length - 1, Math.trunc(finiteNumber(signal.sourceSiScan.minorParent, finiteNumber(event.minorParent, 0))))),
            } : undefined,
          outgroup: signal.outgroup === null
            ? null
            : Number.isFinite(signal.outgroup) && (signal.outgroup as number) >= 0 && (signal.outgroup as number) < alignment.sequences.length
              ? Math.trunc(signal.outgroup as number)
              : undefined,
          outgroupMode: signal.outgroupMode === "most-divergent" || signal.outgroupMode === "randomized" || signal.outgroupMode === "manual"
            ? signal.outgroupMode
            : signal.outgroupMode === "nearest" ? "nearest" : undefined,
          outgroupSampled: signal.outgroupSampled === true,
          permutations: Number.isFinite(signal.permutations) ? Math.max(2, Math.min(10_000, Math.trunc(signal.permutations as number))) : undefined,
          scanPermutations: Number.isFinite(signal.scanPermutations) ? Math.max(2, Math.min(1000, Math.trunc(signal.scanPermutations as number))) : undefined,
          pattern: Number.isFinite(signal.pattern) ? Math.max(0, Math.min(15, Math.trunc(signal.pattern as number))) : undefined,
          scoreFamily: signal.scoreFamily === "sum" ? "sum" : signal.scoreFamily === "pattern" ? "pattern" : undefined,
          baselineTopology: Number.isFinite(signal.baselineTopology) ? Math.max(0, Math.min(2, Math.trunc(signal.baselineTopology as number))) : undefined,
          inferredTopology: Number.isFinite(signal.inferredTopology) ? Math.max(0, Math.min(2, Math.trunc(signal.inferredTopology as number))) : undefined,
          profile: Array.isArray(signal.profile) ? signal.profile.slice(0, 512).flatMap((entry) => {
            if (!entry || typeof entry !== "object") return [];
            return [{
              position: Math.max(0, Math.min(alignment.length - 1, Math.trunc(finiteNumber(entry.position, 0)))),
              z: finiteNumber(entry.z, 0),
              topology: Math.max(0, Math.min(2, Math.trunc(finiteNumber(entry.topology, 0)))),
              baselineTopology: Math.max(0, Math.min(2, Math.trunc(finiteNumber(entry.baselineTopology, 0)))),
              pattern: Number.isFinite(entry.pattern) ? Math.max(0, Math.min(15, Math.trunc(entry.pattern as number))) : undefined,
              scoreFamily: entry.scoreFamily === "sum" ? "sum" as const : entry.scoreFamily === "pattern" ? "pattern" as const : undefined,
            }];
          }) : undefined,
        }];
      }) : undefined,
      evidence,
      chiSquare: finiteNumber(event.chiSquare, 0),
      informativeSites: Math.trunc(finiteNumber(event.informativeSites, 0)),
      decision: event.decision === "accepted" || event.decision === "rejected" ? event.decision : "unreviewed",
      warnings: Array.isArray(event.warnings) ? event.warnings.filter((item): item is string => typeof item === "string") : [],
      note: typeof event.note === "string" ? event.note : "",
      source: event.source === "manual" || event.source === "example" ? event.source : "wasm",
      groupId: typeof event.groupId === "string" && event.groupId.trim() ? event.groupId.trim() : null,
      ancestralCluster: event.ancestralCluster && typeof event.ancestralCluster === "object"
        && typeof event.ancestralCluster.representativeId === "string"
        ? {
            inference: event.ancestralCluster.inference === "manual" ? "manual" : "rdp5-three-set",
            representativeId: event.ancestralCluster.representativeId,
            memberEventIds: Array.isArray(event.ancestralCluster.memberEventIds) ? event.ancestralCluster.memberEventIds.filter((item): item is string => typeof item === "string") : [],
            sequenceMembers: Array.isArray(event.ancestralCluster.sequenceMembers) ? event.ancestralCluster.sequenceMembers.map((item) => Math.trunc(finiteNumber(item, -1))).filter((item) => item >= 0 && item < alignment.sequences.length) : [],
            confidence: Math.max(0, Math.min(1, finiteNumber(event.ancestralCluster.confidence, 0))),
            evidenceCounts: event.ancestralCluster.evidenceCounts && typeof event.ancestralCluster.evidenceCounts === "object" ? {
              phylogenetic: Math.trunc(finiteNumber(event.ancestralCluster.evidenceCounts.phylogenetic, 0)),
              distance: Math.trunc(finiteNumber(event.ancestralCluster.evidenceCounts.distance, 0)),
              detectableSignal: Math.trunc(finiteNumber(event.ancestralCluster.evidenceCounts.detectableSignal, 0)),
              sourceSimilarity: Math.trunc(finiteNumber(event.ancestralCluster.evidenceCounts.sourceSimilarity, 0)),
            } : { phylogenetic: 0, distance: 0, detectableSignal: 0 },
            partialOverprint: event.ancestralCluster.partialOverprint === true,
            sourceMerge: event.ancestralCluster.sourceMerge && typeof event.ancestralCluster.sourceMerge === "object"
              ? {
                  threshold: Math.max(0, Math.min(1, finiteNumber(event.ancestralCluster.sourceMerge.threshold, options.clusterSourceSimilarity))),
                  pairDistances: Array.isArray(event.ancestralCluster.sourceMerge.pairDistances) ? event.ancestralCluster.sourceMerge.pairDistances.flatMap((entry) => {
                    if (!entry || typeof entry !== "object" || !Array.isArray(entry.eventIds) || entry.eventIds.length !== 2) return [];
                    return [{
                      eventIds: [String(entry.eventIds[0]), String(entry.eventIds[1])] as [string, string],
                      distance: finiteNumber(entry.distance, 100),
                      belowThreshold: entry.belowThreshold === true,
                    }];
                  }) : [],
                }
              : undefined,
            pairwise: Array.isArray(event.ancestralCluster.pairwise) ? event.ancestralCluster.pairwise.filter((item): item is Record<string, unknown> => Boolean(item) && typeof item === "object") : [],
          }
        : undefined,
      coRecombinantSets: Array.isArray(event.coRecombinantSets) ? event.coRecombinantSets.flatMap((set) => {
        if (!set || typeof set !== "object") return [];
        const presumedRecombinant = Math.trunc(finiteNumber(set.presumedRecombinant, -1));
        if (presumedRecombinant < 0 || presumedRecombinant >= alignment.sequences.length) return [];
        const parents = Array.isArray(set.parents)
          ? set.parents.map((item) => Math.trunc(finiteNumber(item, -1))).filter((item) => item >= 0 && item < alignment.sequences.length)
          : [];
        const sequenceMembers = Array.isArray(set.sequenceMembers)
          ? [...new Set(set.sequenceMembers.map((item) => Math.trunc(finiteNumber(item, -1))).filter((item) => item >= 0 && item < alignment.sequences.length))]
          : [presumedRecombinant];
        const setEvidence = Array.isArray(set.evidence) ? set.evidence.flatMap((entry) => {
          if (!entry || typeof entry !== "object") return [];
          const sequence = Math.trunc(finiteNumber(entry.sequence, -1));
          if (sequence < 0 || sequence >= alignment.sequences.length) return [];
          const bestCorrelation = entry.bestCorrelation && typeof entry.bestCorrelation === "object"
            ? {
                r: finiteNumber(entry.bestCorrelation.r, 0),
                pValue: Math.max(0, Math.min(1, finiteNumber(entry.bestCorrelation.pValue, 1))),
                inversion: Math.trunc(finiteNumber(entry.bestCorrelation.inversion, 0)),
              }
            : undefined;
          const treeBootstrap = entry.treeBootstrap && typeof entry.treeBootstrap === "object"
            ? {
                replicates: Math.max(0, Math.min(1000, Math.trunc(finiteNumber(entry.treeBootstrap.replicates, 0)))),
                cutoff: Math.max(0, Math.min(1, finiteNumber(entry.treeBootstrap.cutoff, 0.5))),
                cohortTaxa: Math.max(0, Math.trunc(finiteNumber(entry.treeBootstrap.cohortTaxa, 0))),
                sourceSequenceCount: Math.max(0, Math.trunc(finiteNumber(entry.treeBootstrap.sourceSequenceCount, alignment.sequences.length))),
                included: entry.treeBootstrap.included === true,
                exactSiteBootstrap: entry.treeBootstrap.exactSiteBootstrap === true,
                sourceScore: finiteNumber(entry.treeBootstrap.sourceScore, 0),
              }
            : undefined;
          const regionEvidence = Array.isArray(entry.regionEvidence) ? entry.regionEvidence.flatMap((region) => {
            if (!region || typeof region !== "object") return [];
            return [{
              pair: typeof region.pair === "string" ? region.pair : "regional comparison",
              phylogenetic: region.phylogenetic === true,
              movesTogether: region.movesTogether === true,
              sisterTogether: region.sisterTogether === true,
              topologyMargin: finiteNumber(region.topologyMargin, 0),
              treeSourceScore: finiteNumber(region.treeSourceScore, 0),
              bootstrapSupport: finiteNumber(region.bootstrapSupport, 0),
              bootstrapReplicates: Math.max(0, Math.trunc(finiteNumber(region.bootstrapReplicates, 0))),
              bootstrapCutoff: Math.max(0, Math.min(1, finiteNumber(region.bootstrapCutoff, 0.5))),
              treeExcluded: region.treeExcluded === true,
              correlationR: finiteNumber(region.correlationR, 0),
              correlationP: Math.max(0, Math.min(1, finiteNumber(region.correlationP, 1))),
              correlationInversion: Math.trunc(finiteNumber(region.correlationInversion, 0)),
              correlationPermutations: Array.isArray(region.correlationPermutations)
                ? region.correlationPermutations.slice(0, 6).map((value) => finiteNumber(value, 0))
                : undefined,
              correlationSdmFiltered: region.correlationSdmFiltered === true,
            }];
          }) : undefined;
          return [{
            sequence,
            sets: Math.max(0, Math.min(3, Math.trunc(finiteNumber(entry.sets, 0)))),
            phylogenetic: entry.phylogenetic === true,
            distance: entry.distance === true,
            detectableSignal: entry.detectableSignal === true,
            bestCorrelation,
            topologyMargin: finiteNumber(entry.topologyMargin, 0),
            treeBootstrap,
            regionEvidence,
          }];
        }) : [];
        return [{
          presumedRecombinant,
          parents,
          sequenceMembers,
          testedSequences: Math.max(0, Math.trunc(finiteNumber(set.testedSequences, 0))),
          requiredEvidenceSets: Math.max(1, Math.min(3, Math.trunc(finiteNumber(set.requiredEvidenceSets, 2)))),
          evidence: setEvidence,
        }];
      }) : undefined,
      componentProvenance: event.componentProvenance && typeof event.componentProvenance === "object"
        ? (() => {
            const recombinantComponent = parseComponentReference(event.componentProvenance.recombinant);
            const majorParentComponent = parseComponentReference(event.componentProvenance.majorParent);
            const minorParentComponent = parseComponentReference(event.componentProvenance.minorParent);
            if (!recombinantComponent || !majorParentComponent || !minorParentComponent) return undefined;
            return {
              reconstruction: "rdp5-signal-disassembly" as const,
              appliedEventIds: Array.isArray(event.componentProvenance.appliedEventIds)
                ? event.componentProvenance.appliedEventIds.filter((item): item is string => typeof item === "string" && item.length > 0)
                : [],
              recombinant: recombinantComponent,
              majorParent: majorParentComponent,
              minorParent: minorParentComponent,
            };
          })()
        : undefined,
      structuralUncertainty: event.structuralUncertainty && typeof event.structuralUncertainty === "object"
        ? {
            source: "rdp5-erased-signal-boundary",
            originalStart: Math.max(0, Math.min(alignment.length - 1, Math.trunc(finiteNumber(event.structuralUncertainty.originalStart, start)))),
            originalEnd: Math.max(0, Math.min(alignment.length, Math.trunc(finiteNumber(event.structuralUncertainty.originalEnd, end)))),
            originalWraps: event.structuralUncertainty.originalWraps === true,
            piece: Math.max(1, Math.trunc(finiteNumber(event.structuralUncertainty.piece, 1))),
            pieces: Math.max(1, Math.trunc(finiteNumber(event.structuralUncertainty.pieces, 1))),
            uncertainStart: event.structuralUncertainty.uncertainStart === true,
            uncertainEnd: event.structuralUncertainty.uncertainEnd === true,
            adjacentEventIds: Array.isArray(event.structuralUncertainty.adjacentEventIds)
              ? event.structuralUncertainty.adjacentEventIds.filter((item): item is string => typeof item === "string" && item.length > 0)
              : [],
          }
        : undefined,
      recombinantIdentification: parseRecombinantIdentification(event.recombinantIdentification),
      alternativeParents: Array.isArray(event.alternativeParents)
        ? [...new Set(event.alternativeParents.map((item) => Math.trunc(finiteNumber(item, -1))).filter((item) => item >= 0 && item < alignment.sequences.length && item !== recombinant && item !== majorParent && item !== minorParent))]
        : [],
      hypothesisTests: Math.max(1, Math.trunc(finiteNumber(event.hypothesisTests, 1))),
      recalculationNote: typeof event.recalculationNote === "string" ? event.recalculationNote : undefined,
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
        tripletMode: (rawMetrics.tripletMode === "approximate-parent-shortlist" || (rawMetrics.tripletMode === undefined && !options.exhaustive)
          ? "approximate-parent-shortlist"
          : "all-concrete-triplets") as "all-concrete-triplets" | "approximate-parent-shortlist",
        concreteTripletInputs: rawMetrics.concreteTripletInputs !== false,
        parentSamples: finiteNumber(rawMetrics.parentSamples, 0) || undefined,
        rdpSignalTruncations: Math.max(0, Math.trunc(finiteNumber(rawMetrics.rdpSignalTruncations, 0))) || undefined,
        geneconvSignalTruncations: Math.max(0, Math.trunc(finiteNumber(rawMetrics.geneconvSignalTruncations, 0))) || undefined,
        chiSignalTruncations: Math.max(0, Math.trunc(finiteNumber(rawMetrics.chiSignalTruncations, 0))) || undefined,
        bootscanSignalTruncations: Math.max(0, Math.trunc(finiteNumber(rawMetrics.bootscanSignalTruncations, 0))) || undefined,
        bootscanBatch: rawMetrics.bootscanBatch && typeof rawMetrics.bootscanBatch === "object"
          ? {
              calls: Math.max(0, Math.trunc(finiteNumber(rawMetrics.bootscanBatch.calls, 0))),
              triplets: Math.max(0, Math.trunc(finiteNumber(rawMetrics.bootscanBatch.triplets, 0))),
              usedPairs: Math.max(0, Math.trunc(finiteNumber(rawMetrics.bootscanBatch.usedPairs, 0))),
              windows: Math.max(0, Math.trunc(finiteNumber(rawMetrics.bootscanBatch.windows, 0))),
              replicates: Math.max(0, Math.trunc(finiteNumber(rawMetrics.bootscanBatch.replicates, 0))),
              workspaceBytes: Math.max(0, Math.trunc(finiteNumber(rawMetrics.bootscanBatch.workspaceBytes, 0))),
              relationshipMode: "distance" as const,
            }
          : undefined,
        tripletKernelCalls: rawMetrics.tripletKernelCalls && typeof rawMetrics.tripletKernelCalls === "object"
            ? {
              rdp: Math.max(0, Math.trunc(finiteNumber(rawMetrics.tripletKernelCalls.rdp, 0))),
              geneconv: Math.max(0, Math.trunc(finiteNumber(rawMetrics.tripletKernelCalls.geneconv, 0))),
              sourceChi: Math.max(0, Math.trunc(finiteNumber(rawMetrics.tripletKernelCalls.sourceChi, 0))),
              threeSeq: Math.max(0, Math.trunc(finiteNumber(rawMetrics.tripletKernelCalls.threeSeq, 0))),
              siscan: Math.max(0, Math.trunc(finiteNumber(rawMetrics.tripletKernelCalls.siscan, 0))),
            }
          : undefined,
        detectionCycle: rawMetrics.detectionCycle && typeof rawMetrics.detectionCycle === "object" && rawMetrics.detectionCycle.enabled === true
          ? {
              enabled: true,
              eventsApplied: Math.max(0, Math.trunc(finiteNumber(rawMetrics.detectionCycle.eventsApplied, 0))),
              passes: Math.max(1, Math.trunc(finiteNumber(rawMetrics.detectionCycle.passes, 1))),
              initialComparisons: Math.max(0, Math.trunc(finiteNumber(rawMetrics.detectionCycle.initialComparisons, 0))),
              redoComparisons: Math.max(0, Math.trunc(finiteNumber(rawMetrics.detectionCycle.redoComparisons, 0))),
              stoppedBecause: rawMetrics.detectionCycle.stoppedBecause === "cycle-cap" ? "cycle-cap" as const : "no-detectable-signals" as const,
              maximumCycles: Math.max(1, Math.min(1000, Math.trunc(finiteNumber(rawMetrics.detectionCycle.maximumCycles, DEFAULT_OPTIONS.maximumDetectionCycles)))),
            }
          : undefined,
        timing: rawMetrics.timing && typeof rawMetrics.timing === "object"
          ? {
              distanceMs: finiteNumber(rawMetrics.timing.distanceMs, 0),
              scanMs: finiteNumber(rawMetrics.timing.scanMs, 0),
              statisticsMs: finiteNumber(rawMetrics.timing.statisticsMs, 0),
              diagnosticsMs: finiteNumber(rawMetrics.timing.diagnosticsMs, 0) || undefined,
              clusteringMs: finiteNumber(rawMetrics.timing.clusteringMs, 0) || undefined,
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
              proximityStatistic: finiteNumber(rawMetrics.diagnostics.proximityStatistic, 0),
              proximityPermutationP: finiteNumber(rawMetrics.diagnostics.proximityPermutationP, 1),
              proximityPermutationReplicates: Math.trunc(finiteNumber(rawMetrics.diagnostics.proximityPermutationReplicates, 0)),
              ambiguityFraction: finiteNumber(rawMetrics.diagnostics.ambiguityFraction, 0),
              phiPValue: finiteNumber(rawMetrics.diagnostics.phiPValue, 1),
              phiStatistic: finiteNumber(rawMetrics.diagnostics.phiStatistic, 0),
              phiMean: finiteNumber(rawMetrics.diagnostics.phiMean, 0),
              phiVariance: Math.max(0, finiteNumber(rawMetrics.diagnostics.phiVariance, 0)),
              phiZ: finiteNumber(rawMetrics.diagnostics.phiZ, 0),
              phiInformativeSites: Math.max(0, Math.trunc(finiteNumber(rawMetrics.diagnostics.phiInformativeSites, 0))),
              phiTotalInformativeSites: Math.max(0, Math.trunc(finiteNumber(rawMetrics.diagnostics.phiTotalInformativeSites, 0))),
              phiK: Math.max(0, Math.trunc(finiteNumber(rawMetrics.diagnostics.phiK, 0))),
              phiWindow: Math.max(0, Math.trunc(finiteNumber(rawMetrics.diagnostics.phiWindow, 0))),
              phiSubsampled: rawMetrics.diagnostics.phiSubsampled === true,
              phiValidNormalApproximation: rawMetrics.diagnostics.phiValidNormalApproximation === true,
              phiCompatibility: typeof rawMetrics.diagnostics.phiCompatibility === "string"
                ? rawMetrics.diagnostics.phiCompatibility.slice(0, 240)
                : undefined,
            }
          : undefined,
        disassembly: rawMetrics.disassembly && typeof rawMetrics.disassembly === "object"
          ? {
              appliedEvents: Math.max(0, Math.trunc(finiteNumber(rawMetrics.disassembly.appliedEvents, 0))),
              components: Math.max(0, Math.trunc(finiteNumber(rawMetrics.disassembly.components, 0))),
              erasedCanonicalBases: Math.max(0, Math.trunc(finiteNumber(rawMetrics.disassembly.erasedCanonicalBases, 0))),
            }
          : undefined,
      }
    : null;
  const distance = Array.isArray(raw.distance)
    ? raw.distance.map((value) => finiteNumber(value, 0))
    : [];
  const auditLog = Array.isArray(raw.auditLog) ? raw.auditLog.flatMap((value, index): ProjectAuditEntry[] => {
    const entry = value as Partial<ProjectAuditEntry>;
    if (typeof entry.action !== "string") return [];
    return [{
      id: typeof entry.id === "string" ? entry.id : `imported-project-audit-${index + 1}`,
      timestamp: typeof entry.timestamp === "string" ? entry.timestamp : new Date(0).toISOString(),
      action: entry.action,
      summary: typeof entry.summary === "string" ? entry.summary : "Imported project audit entry.",
      eventId: typeof entry.eventId === "string" ? entry.eventId : undefined,
      eventSnapshot: typeof entry.eventSnapshot === "string" ? entry.eventSnapshot : undefined,
    }];
  }) : [];
  return { schema: "rdp-web/0.5", alignment, options, events, metrics, distance, auditLog };
}
