import { AlignmentData, AnalysisOptions, makeDemoAlignment } from "./rdp-core";

export interface ExampleTruth {
  recombinant: string;
  donor: string;
  region: string;
  note: string;
}

export interface ExampleDataset {
  id: string;
  title: string;
  organism: string;
  complexity: "Starter" | "Intermediate" | "Advanced" | "Stress test";
  sequenceCount: number;
  length: number;
  description: string;
  challenge: string;
  tags: string[];
  truth: ExampleTruth[];
  recommendedOptions: Partial<AnalysisOptions>;
  generate: () => AlignmentData;
}

function random(seed: number): () => number {
  let state = seed >>> 0;
  return () => {
    state = (Math.imul(state, 1664525) + 1013904223) >>> 0;
    return state / 0x1_0000_0000;
  };
}

function ancestor(length: number, seed: number): string {
  const rng = random(seed);
  const bases = "ACGT";
  const output = new Array<string>(length);
  for (let index = 0; index < length; index += 1) output[index] = bases[Math.floor(rng() * 4)];
  return output.join("");
}

function mutate(sequence: string, rate: number, seed: number): string {
  const rng = random(seed);
  const output = sequence.split("");
  const alternatives: Record<string, string> = { A: "CGT", C: "AGT", G: "ACT", T: "ACG" };
  for (let index = 0; index < output.length; index += 1) {
    if (rng() < rate) output[index] = alternatives[output[index]][Math.floor(rng() * 3)];
  }
  return output.join("");
}

function mosaic(target: string, donor: string, start: number, end: number): string {
  if (start <= end) return `${target.slice(0, start)}${donor.slice(start, end)}${target.slice(end)}`;
  return `${donor.slice(0, end)}${target.slice(end, start)}${donor.slice(start)}`;
}

function generatedAlignment(
  name: string,
  sequences: Array<{ name: string; sequence: string; role?: "query" | "reference" | "both" }>,
): AlignmentData {
  return {
    name,
    format: "generated",
    sequences,
    length: sequences[0]?.sequence.length ?? 0,
    createdAt: Date.now(),
  };
}

function minimalTriplet(): AlignmentData {
  const root = ancestor(1_200, 101);
  const major = mutate(root, 0.035, 102);
  const minor = mutate(root, 0.13, 103);
  const recombinant = mutate(mosaic(major, minor, 360, 760), 0.006, 104);
  return generatedAlignment("Minimal triplet control", [
    { name: "Mosaic_query", sequence: recombinant, role: "query" },
    { name: "Major_parent", sequence: major, role: "reference" },
    { name: "Minor_parent", sequence: minor, role: "reference" },
    ...Array.from({ length: 5 }, (_, index) => ({
      name: `Background_${index + 1}`,
      sequence: mutate(index < 3 ? major : minor, 0.012, 110 + index),
      role: "reference" as const,
    })),
  ]);
}

function virusFamily(): AlignmentData {
  const root = ancestor(12_000, 201);
  const templates = [mutate(root, 0.035, 202), mutate(root, 0.075, 203), mutate(root, 0.115, 204)];
  const sequences = Array.from({ length: 36 }, (_, index) => {
    const clade = Math.floor(index / 12);
    return {
      name: `Virus_${String.fromCharCode(65 + clade)}_${String(index % 12 + 1).padStart(2, "0")}`,
      sequence: mutate(templates[clade], 0.008 + (index % 4) * 0.002, 220 + index),
      role: (index % 12 < 2 ? "reference" : "both") as "reference" | "both",
    };
  });
  sequences[9].sequence = mosaic(sequences[9].sequence, sequences[16].sequence, 2_100, 4_650);
  sequences[22].sequence = mosaic(sequences[22].sequence, sequences[31].sequence, 7_250, 10_300);
  sequences[35].sequence = mosaic(sequences[35].sequence, sequences[3].sequence, 4_800, 6_150);
  return generatedAlignment("Three-lineage RNA virus family", sequences);
}

function circularFamily(): AlignmentData {
  const root = ancestor(3_200, 301);
  const cladeA = mutate(root, 0.045, 302);
  const cladeB = mutate(root, 0.13, 303);
  const sequences = Array.from({ length: 28 }, (_, index) => ({
    name: `Circular_${index < 14 ? "A" : "B"}_${String(index % 14 + 1).padStart(2, "0")}`,
    sequence: mutate(index < 14 ? cladeA : cladeB, 0.012, 320 + index),
    role: "both" as const,
  }));
  sequences[7].sequence = mosaic(sequences[7].sequence, sequences[19].sequence, 2_720, 470);
  sequences[11].sequence = mosaic(sequences[11].sequence, sequences[22].sequence, 1_180, 1_940);
  return generatedAlignment("Circular ssDNA family", sequences);
}

function segmentedFamily(): AlignmentData {
  const segmentLength = 6_000;
  const roots = Array.from({ length: 4 }, (_, index) => ancestor(segmentLength, 400 + index));
  const cladeA = roots.map((root, index) => mutate(root, 0.035, 410 + index));
  const cladeB = roots.map((root, index) => mutate(root, 0.095, 420 + index));
  const sequences = Array.from({ length: 48 }, (_, index) => {
    const clade = index < 24 ? cladeA : cladeB;
    return {
      name: `Segmented_${index < 24 ? "A" : "B"}_${String(index % 24 + 1).padStart(2, "0")}`,
      sequence: clade.map((segment, segmentIndex) => mutate(segment, 0.006, 450 + index * 7 + segmentIndex)).join(""),
      role: "both" as const,
    };
  });
  sequences[8].sequence = mosaic(sequences[8].sequence, sequences[31].sequence, 6_000, 12_000);
  sequences[18].sequence = mosaic(sequences[18].sequence, sequences[39].sequence, 18_000, 24_000);
  sequences[42].sequence = mosaic(sequences[42].sequence, sequences[4].sequence, 12_000, 18_000);
  return generatedAlignment("Four-segment virus family", sequences);
}

function nestedMosaics(): AlignmentData {
  const root = ancestor(18_000, 501);
  const templates = [mutate(root, 0.035, 502), mutate(root, 0.08, 503), mutate(root, 0.14, 504), mutate(root, 0.19, 505)];
  const sequences = Array.from({ length: 60 }, (_, index) => ({
    name: `Deep_clade_${Math.floor(index / 15) + 1}_${String(index % 15 + 1).padStart(2, "0")}`,
    sequence: mutate(templates[Math.floor(index / 15)], 0.008, 520 + index),
    role: "both" as const,
  }));
  sequences[12].sequence = mosaic(sequences[12].sequence, sequences[22].sequence, 2_400, 12_600);
  sequences[12].sequence = mosaic(sequences[12].sequence, sequences[49].sequence, 6_050, 7_900);
  sequences[37].sequence = mosaic(sequences[37].sequence, sequences[4].sequence, 10_100, 15_800);
  sequences[55].sequence = mosaic(sequences[55].sequence, sequences[26].sequence, 16_200, 1_100);
  return generatedAlignment("Nested and origin-spanning mosaics", sequences);
}

function bacterialCore(): AlignmentData {
  const root = ancestor(80_000, 601);
  const templates = [mutate(root, 0.025, 602), mutate(root, 0.05, 603), mutate(root, 0.082, 604)];
  const sequences = Array.from({ length: 72 }, (_, index) => ({
    name: `Bacterium_ST${Math.floor(index / 24) + 1}_${String(index % 24 + 1).padStart(2, "0")}`,
    sequence: mutate(templates[Math.floor(index / 24)], 0.0035, 620 + index),
    role: "both" as const,
  }));
  sequences[14].sequence = mosaic(sequences[14].sequence, sequences[39].sequence, 7_500, 14_200);
  sequences[21].sequence = mosaic(sequences[21].sequence, sequences[61].sequence, 31_000, 44_500);
  sequences[46].sequence = mosaic(sequences[46].sequence, sequences[5].sequence, 58_200, 65_300);
  sequences[68].sequence = mosaic(sequences[68].sequence, sequences[28].sequence, 18_000, 24_400);
  return generatedAlignment("Bacterial core-genome HGT panel", sequences);
}

function outbreakStress(): AlignmentData {
  const root = ancestor(12_000, 701);
  const cladeA = mutate(root, 0.012, 702);
  const cladeB = mutate(root, 0.032, 703);
  const sequences = Array.from({ length: 520 }, (_, index) => ({
    name: `Outbreak_${String(index + 1).padStart(4, "0")}`,
    sequence: mutate(index < 360 ? cladeA : cladeB, 0.0015 + (index % 5) * 0.0002, 720 + index),
    role: (index < 40 ? "reference" : "query") as "reference" | "query",
  }));
  for (let index = 0; index < 24; index += 1) {
    const recombinant = 40 + index * 11;
    const donor = 380 + index * 5;
    const start = 900 + (index * 431) % 8_000;
    sequences[recombinant].sequence = mosaic(sequences[recombinant].sequence, sequences[donor].sequence, start, Math.min(11_700, start + 1_200 + (index % 4) * 350));
  }
  return generatedAlignment("Large outbreak performance panel", sequences);
}

export const EXAMPLE_DATASETS: ExampleDataset[] = [
  {
    id: "minimal-triplet",
    title: "Minimal triplet control",
    organism: "Generic nucleotide alignment",
    complexity: "Starter",
    sequenceCount: 8,
    length: 1_200,
    description: "One clean mosaic, two known parents, and five nearby controls.",
    challenge: "Use this to learn event selection, method evidence, breakpoint editing, and local-tree switching.",
    tags: ["known parents", "single event", "fast"],
    truth: [{ recombinant: "Mosaic_query", donor: "Minor_parent", region: "361–760", note: "Single internal donor tract." }],
    recommendedOptions: { mode: "exploratory", circular: false, window: 120, step: 20, candidateParents: 7 },
    generate: minimalTriplet,
  },
  {
    id: "tutorial-virus",
    title: "RDP Web tutorial control",
    organism: "Synthetic virus-like genomes",
    complexity: "Starter",
    sequenceCount: 12,
    length: 2_400,
    description: "The original deterministic positive control, now loaded without pre-authored method results.",
    challenge: "A moderate-divergence tract tests basic localization and parent assignment.",
    tags: ["tutorial", "known truth", "small"],
    truth: [{ recombinant: "Mosaic-X (and Mosaic-Y)", donor: "Beta-01 lineage", region: "783–1,538", note: "Known internal donor tract; Alpha-01 is the major-parent lineage." }],
    recommendedOptions: { mode: "exploratory", circular: false, window: 180, step: 30 },
    generate: makeDemoAlignment,
  },
  {
    id: "virus-family",
    title: "Three-lineage virus family",
    organism: "Synthetic RNA virus family",
    complexity: "Intermediate",
    sequenceCount: 36,
    length: 12_000,
    description: "Three divergent lineages with three independently introduced mosaics.",
    challenge: "Tests multiple events, incomplete parent sampling, and whether tree placement changes by region.",
    tags: ["virus family", "multiple events", "deep lineages"],
    truth: [
      { recombinant: "Virus_A_10", donor: "Virus_B_05", region: "2,101–4,650", note: "A→B internal tract." },
      { recombinant: "Virus_B_11", donor: "Virus_C_08", region: "7,251–10,300", note: "B→C internal tract." },
      { recombinant: "Virus_C_12", donor: "Virus_A_04", region: "4,801–6,150", note: "C→A short tract." },
    ],
    recommendedOptions: { mode: "exploratory", circular: false, window: 420, step: 70, candidateParents: 12 },
    generate: virusFamily,
  },
  {
    id: "circular-family",
    title: "Circular ssDNA family",
    organism: "Synthetic circular virus family",
    complexity: "Intermediate",
    sequenceCount: 28,
    length: 3_200,
    description: "Two clades with an origin-spanning recombinant and one internal event.",
    challenge: "Exercises circular coordinates, complementary arcs, and three-region tree comparison.",
    tags: ["circular", "origin spanning", "virus"],
    truth: [
      { recombinant: "Circular_A_08", donor: "Circular_B_06", region: "2,721–3,200 + 1–470", note: "Origin-spanning tract." },
      { recombinant: "Circular_A_12", donor: "Circular_B_09", region: "1,181–1,940", note: "Internal tract." },
    ],
    recommendedOptions: { mode: "exploratory", circular: true, window: 180, step: 30, candidateParents: 12 },
    generate: circularFamily,
  },
  {
    id: "segmented-family",
    title: "Four-segment virus family",
    organism: "Synthetic segmented virus",
    complexity: "Advanced",
    sequenceCount: 48,
    length: 24_000,
    description: "Four concatenated 6 kb segments with three complete-segment ancestry switches.",
    challenge: "A reassortment-like challenge: interpret segment-boundary signals manually; automatic reassortment labeling is not claimed.",
    tags: ["segments", "reassortment challenge", "large"],
    truth: [
      { recombinant: "Segmented_A_09", donor: "Segmented_B_08", region: "6,001–12,000", note: "Complete segment 2 replacement." },
      { recombinant: "Segmented_A_19", donor: "Segmented_B_16", region: "18,001–24,000", note: "Complete segment 4 replacement." },
      { recombinant: "Segmented_B_19", donor: "Segmented_A_05", region: "12,001–18,000", note: "Complete segment 3 replacement." },
    ],
    recommendedOptions: { mode: "exploratory", circular: false, window: 600, step: 100, candidateParents: 14 },
    generate: segmentedFamily,
  },
  {
    id: "nested-mosaics",
    title: "Nested and overlapping mosaics",
    organism: "Synthetic deeply divergent family",
    complexity: "Advanced",
    sequenceCount: 60,
    length: 18_000,
    description: "Four clades with long, nested, overlapping, and origin-spanning donor tracts.",
    challenge: "Exercises sequential erase/extract disassembly, nested-event ordering, co-recombinant grouping, and selective affected-triplet rescanning.",
    tags: ["nested", "overlap", "unknown ancestry"],
    truth: [
      { recombinant: "Deep_clade_1_13", donor: "Deep_clade_2_08", region: "2,401–12,600", note: "Long tract containing a nested second donor." },
      { recombinant: "Deep_clade_1_13", donor: "Deep_clade_4_05", region: "6,051–7,900", note: "Nested donor tract." },
      { recombinant: "Deep_clade_3_08", donor: "Deep_clade_1_05", region: "10,101–15,800", note: "Overlapping-lineage challenge." },
      { recombinant: "Deep_clade_4_11", donor: "Deep_clade_2_12", region: "16,201–18,000 + 1–1,100", note: "Origin-spanning tract." },
    ],
    recommendedOptions: { mode: "exploratory", circular: true, window: 540, step: 90, candidateParents: 16 },
    generate: nestedMosaics,
  },
  {
    id: "bacterial-core",
    title: "Bacterial core-genome HGT panel",
    organism: "Synthetic bacterial core genome",
    complexity: "Advanced",
    sequenceCount: 72,
    length: 80_000,
    description: "Three sequence types with four horizontal-transfer tracts across an 80 kb core alignment.",
    challenge: "Tests long alignments, sparse within-lineage variation, broad donor tracts, and sampled tree inference.",
    tags: ["bacteria", "HGT", "80 kb"],
    truth: [
      { recombinant: "Bacterium_ST1_15", donor: "Bacterium_ST2_16", region: "7,501–14,200", note: "ST2 donor tract." },
      { recombinant: "Bacterium_ST1_22", donor: "Bacterium_ST3_14", region: "31,001–44,500", note: "Long ST3 donor tract." },
      { recombinant: "Bacterium_ST2_23", donor: "Bacterium_ST1_06", region: "58,201–65,300", note: "ST1 donor tract." },
      { recombinant: "Bacterium_ST3_21", donor: "Bacterium_ST2_05", region: "18,001–24,400", note: "ST2 donor tract." },
    ],
    recommendedOptions: { mode: "exploratory", circular: false, window: 1_500, step: 250, candidateParents: 16, bootstrapReplicates: 50 },
    generate: bacterialCore,
  },
  {
    id: "outbreak-stress",
    title: "Large outbreak performance panel",
    organism: "Synthetic low-diversity outbreak",
    complexity: "Stress test",
    sequenceCount: 520,
    length: 12_000,
    description: "More than 500 low-diversity genomes with 24 introduced donor tracts.",
    challenge: "For benchmarking the sampled parent screen, worker responsiveness, large-table behavior, and memory use—not numerical validation.",
    tags: ["520 genomes", "performance", "low diversity"],
    truth: [{ recombinant: "24 known mosaics", donor: "Outbreak clade B", region: "Variable 1.2–2.25 kb tracts", note: "Truth formula is documented in the example generator." }],
    recommendedOptions: { mode: "exploratory", circular: false, window: 480, step: 120, candidateParents: 10, bootstrapReplicates: 25 },
    generate: outbreakStress,
  },
];
