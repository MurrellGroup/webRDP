"use client";

import {
  ChangeEvent,
  DragEvent,
  ReactNode,
  useCallback,
  useDeferredValue,
  useEffect,
  useMemo,
  useRef,
  useState,
} from "react";
import {
  AlignmentData,
  AlignmentDiagnostics,
  AnalysisOptions,
  DEFAULT_OPTIONS,
  EventDecision,
  MethodName,
  NeighborJoiningNode,
  PRIMARY_METHODS,
  ProjectAuditEntry,
  RdpProject,
  RdpEvent,
  alignmentStats,
  breakpointHotspotTest,
  buildLocalTree,
  evidenceProfile,
  eventLength,
  eventSegments,
  exportRecombinationFree,
  formatEventRegion,
  makeDemoAlignment,
  pairwiseIdentitySampled,
  pairwiseIdentitySegments,
  parseAlignment,
  parseGenomeAnnotations,
  parseProject,
  serializeProject,
  toFasta,
  toGff3,
} from "./rdp-core";
import { EXAMPLE_DATASETS, ExampleDataset } from "./example-datasets";
import { AutosaveRecord, clearAutosave, loadAutosave, saveAutosave } from "./local-project-store";
import {
  affinityDescription,
  classifyParentAffinity,
  parentInformativeSites,
} from "./alignment-highlighter";
import { formatClockTime, formatDateTime, formatInteger } from "./format";
import { buildReconstructionModel, type ReconstructionRelationship } from "./reconstruction";
import { layoutNeighborJoiningTree } from "./tree-layout";
import {
  computeBreakpointPairDensity,
  computeLocalDiscordanceMatrices,
  computeRegionSeparationMatrices,
  eventContainsPosition,
} from "./pattern-matrices";
import {
  AUTO_RESOLVE_PRESETS,
  applyAutoResolutionPlan,
  filterResolvedEventDuplicates,
  planAutoResolution,
  rescanTargetsForBarrier,
  type AutoResolvePresetName,
  type AutoResolveSettings,
} from "./auto-resolve";

type Tab = "explore" | "reconstruction" | "trees" | "alignment" | "patterns" | "export" | "methods";
type RunState = "idle" | "running" | "complete" | "error";
type ExportScope = "accepted-fresh" | "all-fresh" | "all-retained";

interface HistoryFrame {
  label: string;
  events: RdpEvent[];
  selectedId: string | null;
}

interface RunMetrics {
  elapsedMs: number;
  comparisons: number;
  engine: string;
  matrixMode?: string;
  parentSamples?: number;
  timing?: { distanceMs: number; scanMs: number; statisticsMs: number; diagnosticsMs?: number };
  diagnostics?: AlignmentDiagnostics;
}

interface AutoResolveStatus {
  state: "idle" | "resolving" | "rescanning" | "complete" | "paused" | "error";
  message: string;
  round: number;
  processed: number;
  accepted: number;
  rejected: number;
  reviewed: number;
  rescans: number;
}

interface AutoResolveSession {
  settings: AutoResolveSettings;
  profileLabel: string;
  round: number;
  processedEventIds: string[];
  processed: number;
  accepted: number;
  rejected: number;
  reviewed: number;
  rescans: number;
}

interface AnalysisLaunchConfig {
  excludedTargets?: number[];
  excludedParents?: number[];
  retainedEvents?: RdpEvent[];
  filterResolvedAgainst?: RdpEvent[];
  resetEditHistory?: boolean;
  auditAction?: string;
  auditContext?: string;
  onComplete?: (events: RdpEvent[]) => void;
}

const METHOD_META: Record<MethodName, { family: string; detail: string; limit: string; citation: string }> = {
  RDP: {
    family: "Identity / triplet",
    detail: "Informative-site triplet scan with binomial assessment of identity runs.",
    limit: "Current candidate generation is shared with the triplet CUSUM screen; numerical parity with desktop RDP is not established.",
    citation: "https://doi.org/10.1093/bioinformatics/16.6.562",
  },
  GENECONV: {
    family: "Substitution distribution",
    detail: "Unusually long fragments of pairwise identity relative to the alignment background.",
    limit: "Uses a G-scale-0 run bound, not the complete GENECONV mismatch-penalty and permutation model.",
    citation: "https://doi.org/10.1006/viro.1999.0058",
  },
  BootScan: {
    family: "Phylogenetic / windowed",
    detail: "Windowed support for changes in the recombinant’s closest relative.",
    limit: "Uses seeded p-distance triplet resampling, not full multi-taxon neighbor-joining bootstrap parity.",
    citation: "https://doi.org/10.1089/aid.2005.21.98",
  },
  MaxChi: {
    family: "Substitution distribution",
    detail: "Maximum chi-square contrasts across variable sites around candidate breakpoints.",
    limit: "Evaluates shared candidate boundaries; it is not yet an independent MaxChi candidate scan.",
    citation: "https://doi.org/10.1007/BF00182389",
  },
  Chimaera: {
    family: "Substitution distribution",
    detail: "Two-state refinement of MaxChi-style breakpoint evidence in sequence triplets.",
    limit: "Evaluates shared candidate boundaries; full Chimaera scan and calibration parity remain validation work.",
    citation: "https://doi.org/10.1073/pnas.241370698",
  },
  SiScan: {
    family: "Site-category / permutation",
    detail: "Fast oriented category-Z evidence; outgroup permutations remain on the parity track.",
    limit: "This is a category-Z surrogate without the complete SiScan permutation and outgroup procedure.",
    citation: "https://doi.org/10.1093/bioinformatics/16.7.573",
  },
  "3Seq": {
    family: "Non-parametric triplet",
    detail: "Maximum-descent hypergeometric random walk with exact bounded dynamic-program calibration.",
    limit: "Exact DP is operation-bounded and falls back to a conservative bound; candidates still come from the shared screen.",
    citation: "https://doi.org/10.1093/molbev/msx263",
  },
};

const REFERENCES = [
  {
    year: "2021",
    tag: "RDP5",
    title: "RDP5: a computer program for analyzing recombination in, and removing signals of recombination from, nucleotide sequence datasets",
    authors: "Martin et al.",
    href: "https://doi.org/10.1093/ve/veaa087",
    note: "Automation, false-positive flags, query/reference mode, recombination-free exports, and performance.",
  },
  {
    year: "2025",
    tag: "Tutorial",
    title: "Recombination Analysis of Geminiviruses Using RDP",
    authors: "Sattar et al.",
    href: "https://pubmed.ncbi.nlm.nih.gov/40064777/",
    note: "Current practical protocol spanning dataset preparation, event characterization, and recombination-free outputs.",
  },
  {
    year: "2017",
    tag: "Tutorial",
    title: "Detecting and Analyzing Genetic Recombination Using RDP4",
    authors: "Martin, Murrell, Khoosal & Muhire",
    href: "https://doi.org/10.1007/978-1-4939-6622-6_17",
    note: "Practical workflow for formulating, testing, and refining recombination hypotheses.",
  },
  {
    year: "2015",
    tag: "RDP4",
    title: "RDP4: Detection and analysis of recombination patterns in virus genomes",
    authors: "Martin et al.",
    href: "https://doi.org/10.1093/ve/vev003",
    note: "Interactive verification, breakpoint matrices, recombination-aware phylogenetics, and pattern analyses.",
  },
  {
    year: "2011",
    tag: "Review",
    title: "Analysing recombination in nucleotide sequences",
    authors: "Martin, Lemey & Posada",
    href: "https://doi.org/10.1111/j.1755-0998.2011.03026.x",
    note: "Method selection, failure modes, and best-practice interpretation.",
  },
  {
    year: "2010",
    tag: "RDP3",
    title: "RDP3: a flexible and fast computer program for analyzing recombination",
    authors: "Martin et al.",
    href: "https://pubmed.ncbi.nlm.nih.gov/20798170/",
    note: "Integrated exploratory analysis and event characterization.",
  },
  {
    year: "2005",
    tag: "RDP2",
    title: "RDP2: recombination detection and analysis from sequence alignments",
    authors: "Martin, Williamson & Posada",
    href: "https://pubmed.ncbi.nlm.nih.gov/15377507/",
    note: "Early multi-method graphical analysis environment.",
  },
  {
    year: "2000",
    tag: "RDP",
    title: "RDP: detection of recombination amongst aligned sequences",
    authors: "Martin & Rybicki",
    href: "https://doi.org/10.1093/bioinformatics/16.6.562",
    note: "Original informative-site triplet method.",
  },
  {
    year: "2024",
    tag: "Benchmark",
    title: "Evaluation of recombination detection methods for viral sequencing",
    authors: "Jaya et al.",
    href: "https://pmc.ncbi.nlm.nih.gov/articles/PMC10734630/",
    note: "Simulation-based comparison across diversity, frequency, and dataset scale.",
  },
  {
    year: "Manual",
    tag: "Guide",
    title: "RDP5 Instruction Manual",
    authors: "Darren P. Martin",
    href: "https://web.cbio.uct.ac.za/~darren/RDP5Manual.pdf",
    note: "Complete method settings, displays, review workflow, export modes, and step-by-step tutorial.",
  },
];

function formatP(value: number): string {
  if (!Number.isFinite(value)) return "—";
  if (value === 0 || value < 1e-99) return "<1×10⁻⁹⁹";
  if (value < 0.001) {
    const [coefficient, exponent] = value.toExponential(1).split("e");
    return `${coefficient}×10${Number(exponent).toString().replace("-", "⁻")}`;
  }
  return value.toPrecision(2);
}

function tutorialTruthEvent(): RdpEvent {
  return {
    id: "tutorial-known-truth",
    recombinant: 0,
    majorParent: 2,
    minorParent: 5,
    start: 782,
    end: 1538,
    wraps: false,
    confidenceStart: [782, 782],
    confidenceEnd: [1538, 1538],
    breakpointModel: { method: "manual", informativeSites: 0 },
    evidence: [],
    chiSquare: 0,
    informativeSites: 0,
    decision: "unreviewed",
    warnings: ["Known synthetic truth for orientation only; this is not a scan result. Recalculate or run a full scan before interpretation."],
    note: "Tutorial truth: a Beta-derived tract was inserted at sites 783–1,538.",
    source: "example",
    groupId: "known-truth",
    alternativeParents: [6, 7, 11],
    hypothesisTests: 1,
    history: [{
      id: "tutorial-known-truth-history",
      timestamp: "2026-01-01T00:00:00.000Z",
      action: "Loaded known synthetic truth",
      summary: "No method p-values were pre-authored; this event exists only to demonstrate review and tree views.",
    }],
    evidenceStale: true,
    diagnostics: { tractVariableDensity: 0, backgroundVariableDensity: 0, rateRatio: 1, parentConflictRate: 0, parentDiscriminatingSites: 0, diffuseIncompatibility: false },
  };
}

function download(filename: string, content: string, type = "text/plain"): void {
  const blob = new Blob([content], { type });
  const url = URL.createObjectURL(blob);
  const anchor = document.createElement("a");
  anchor.href = url;
  anchor.download = filename;
  anchor.click();
  URL.revokeObjectURL(url);
}

function Icon({ name, size = 18 }: { name: string; size?: number }) {
  const paths: Record<string, ReactNode> = {
    run: <path d="M8 5v14l11-7z" />,
    stop: <path d="M7 7h10v10H7z" />,
    upload: <><path d="M12 16V4m0 0L7.5 8.5M12 4l4.5 4.5" /><path d="M5 15v4h14v-4" /></>,
    chevron: <path d="m9 6 6 6-6 6" />,
    check: <path d="m5 12 4 4L19 6" />,
    close: <path d="m6 6 12 12M18 6 6 18" />,
    help: <><circle cx="12" cy="12" r="9" /><path d="M9.8 9a2.4 2.4 0 1 1 3.3 2.2c-.8.4-1.1.9-1.1 1.8m0 3h.01" /></>,
    download: <><path d="M12 4v11m0 0 4-4m-4 4-4-4" /><path d="M5 19h14" /></>,
    settings: <><circle cx="12" cy="12" r="3" /><path d="M19 13.5V10.5l-2-.7-.5-1.2.9-1.9-2.1-2.1-1.9.9-1.2-.5-.7-2h-3l-.7 2-1.2.5-1.9-.9-2.1 2.1.9 1.9-.5 1.2-2 .7v3l2 .7.5 1.2-.9 1.9 2.1 2.1 1.9-.9 1.2.5.7 2h3l.7-2 1.2-.5 1.9.9 2.1-2.1-.9-1.9.5-1.2z" /></>,
    undo: <><path d="M9 7 4 12l5 5"/><path d="M5 12h8a6 6 0 0 1 6 6"/></>,
    redo: <><path d="m15 7 5 5-5 5"/><path d="M19 12h-8a6 6 0 0 0-6 6"/></>,
    expand: <><path d="M9 4H4v5M15 4h5v5M20 15v5h-5M4 15v5h5"/><path d="m4 9 5-5m6 0 5 5m0 6-5 5M9 20l-5-5"/></>,
    collapse: <><path d="M9 9H4V4M15 9h5V4M20 15h-5v5M4 15h5v5"/><path d="m4 4 5 5m6 0 5-5m0 16-5-5M9 15l-5 5"/></>,
    star: <path d="m12 3 2.8 5.7 6.2.9-4.5 4.4 1.1 6.2-5.6-2.9-5.6 2.9 1.1-6.2L3 9.6l6.2-.9z" />,
  };
  return (
    <svg className="icon" viewBox="0 0 24 24" width={size} height={size} aria-hidden="true">
      {paths[name] ?? paths.chevron}
    </svg>
  );
}

function Segmented<T extends string>({ value, options, onChange }: {
  value: T;
  options: { value: T; label: string }[];
  onChange: (value: T) => void;
}) {
  return (
    <div className="segmented">
      {options.map((option) => (
        <button
          type="button"
          key={option.value}
          className={value === option.value ? "active" : ""}
          onClick={() => onChange(option.value)}
        >
          {option.label}
        </button>
      ))}
    </div>
  );
}

function Panel({ title, action, children, className = "" }: {
  title: string;
  action?: ReactNode;
  children: ReactNode;
  className?: string;
}) {
  const [expanded, setExpanded] = useState(false);
  const panelRef = useRef<HTMLElement>(null);
  useEffect(() => {
    if (!expanded) return;
    const previousFocus = document.activeElement as HTMLElement | null;
    const close = (event: globalThis.KeyboardEvent) => {
      if (event.key === "Escape") setExpanded(false);
    };
    document.body.classList.add("panel-is-expanded");
    window.addEventListener("keydown", close);
    panelRef.current?.focus();
    return () => {
      document.body.classList.remove("panel-is-expanded");
      window.removeEventListener("keydown", close);
      previousFocus?.focus();
    };
  }, [expanded]);
  return <>
    {expanded && <button type="button" className="panel-expanded-backdrop" aria-label={`Close expanded ${title}`} onClick={() => setExpanded(false)}/>}
    <section ref={panelRef} tabIndex={expanded ? -1 : undefined} className={`panel ${expanded ? "panel-expanded" : ""} ${className}`}>
      <div className="panel-title"><h2>{title}</h2><div className="panel-actions">{action}<button type="button" className="panel-expand-button" aria-label={`${expanded ? "Exit full screen for" : "Expand"} ${title}`} aria-expanded={expanded} title={expanded ? "Exit full screen (Esc)" : "Expand panel to full screen"} onClick={() => setExpanded((current) => !current)}><Icon name={expanded ? "collapse" : "expand"} size={15}/></button></div></div>
      <div className="panel-body" tabIndex={0} role="region" aria-label={`${title} content`}>
        {children}
      </div>
    </section>
  </>;
}

function EvidencePlot({ alignment, event, window, circular, onUpdate }: {
  alignment: AlignmentData;
  event: RdpEvent;
  window: number;
  circular: boolean;
  onUpdate?: (patch: Partial<RdpEvent>, action?: string) => void;
}) {
  const [draggingBreakpoint, setDraggingBreakpoint] = useState<"start" | "end" | null>(null);
  const [preview, setPreview] = useState<{ start: number; end: number } | null>(null);
  const displayEvent = preview ? { ...event, ...preview } : event;
  const profile = useMemo(
    () => evidenceProfile(alignment, event, window, 130, circular),
    [alignment, circular, event, window],
  );
  const width = 900;
  const height = 222;
  const left = 42;
  const right = 16;
  const top = 18;
  const bottom = 32;
  const innerWidth = width - left - right;
  const innerHeight = height - top - bottom;
  const yMin = Math.max(0, Math.min(...profile.flatMap((point) => [point.recombinantMajor, point.recombinantMinor, point.parentParent])) - 0.04);
  const x = (position: number) => left + (position / alignment.length) * innerWidth;
  const y = (identity: number) => top + (1 - (identity - yMin) / Math.max(0.001, 1 - yMin)) * innerHeight;
  const path = (key: "recombinantMajor" | "recombinantMinor" | "parentParent") =>
    profile.map((point, index) => `${index ? "L" : "M"}${x(point.position).toFixed(1)},${y(point[key]).toFixed(1)}`).join(" ");
  return (
    <div className="plot-wrap">
      <svg
        viewBox={`0 0 ${width} ${height}`}
        role="img"
        aria-label="Sliding-window pairwise identity plot"
        onPointerMove={(pointer) => {
          if (!draggingBreakpoint || !preview) return;
          const rectangle = pointer.currentTarget.getBoundingClientRect();
          const viewPosition = (pointer.clientX - rectangle.left) / Math.max(1, rectangle.width) * width;
          const requested = Math.round((viewPosition - left) / innerWidth * alignment.length);
          if (draggingBreakpoint === "start") {
            const start = event.wraps
              ? Math.max(preview.end + 1, Math.min(alignment.length - 1, requested))
              : Math.max(0, Math.min(preview.end - 1, requested));
            setPreview({ ...preview, start });
          } else {
            const end = event.wraps
              ? Math.max(0, Math.min(preview.start - 1, requested))
              : Math.max(preview.start + 1, Math.min(alignment.length, requested));
            setPreview({ ...preview, end });
          }
        }}
        onPointerUp={() => {
          if (draggingBreakpoint && preview && (preview.start !== event.start || preview.end !== event.end)) {
            onUpdate?.({
              ...preview,
              confidenceStart: [preview.start, preview.start],
              confidenceEnd: [preview.end, preview.end],
              breakpointModel: { method: "manual", informativeSites: event.informativeSites },
            }, "Dragged breakpoint handle");
          }
          setDraggingBreakpoint(null);
          setPreview(null);
        }}
        onPointerCancel={() => { setDraggingBreakpoint(null); setPreview(null); }}
      >
        <defs>
          <linearGradient id="tract" x1="0" y1="0" x2="0" y2="1"><stop stopColor="#ff7a55" stopOpacity=".16"/><stop offset="1" stopColor="#ff7a55" stopOpacity=".03"/></linearGradient>
        </defs>
        {[0, 0.25, 0.5, 0.75, 1].map((fraction) => (
          <line key={fraction} x1={left} x2={width - right} y1={top + fraction * innerHeight} y2={top + fraction * innerHeight} className="grid-line" />
        ))}
        {eventSegments(displayEvent, alignment.length).map(([start, end], index) => (
          <rect key={`${start}-${end}-${index}`} x={x(start)} y={top} width={Math.max(1, x(end) - x(start))} height={innerHeight} fill="url(#tract)" />
        ))}
        <line x1={x(displayEvent.start)} x2={x(displayEvent.start)} y1={top} y2={top + innerHeight} className="breakpoint-line" />
        <line x1={x(displayEvent.end)} x2={x(displayEvent.end)} y1={top} y2={top + innerHeight} className="breakpoint-line" />
        {onUpdate && <>
          <line x1={x(displayEvent.start)} x2={x(displayEvent.start)} y1={top} y2={top + innerHeight} className="breakpoint-handle" onPointerDown={(pointer) => { pointer.currentTarget.setPointerCapture(pointer.pointerId); setDraggingBreakpoint("start"); setPreview({ start: event.start, end: event.end }); }}/>
          <line x1={x(displayEvent.end)} x2={x(displayEvent.end)} y1={top} y2={top + innerHeight} className="breakpoint-handle" onPointerDown={(pointer) => { pointer.currentTarget.setPointerCapture(pointer.pointerId); setDraggingBreakpoint("end"); setPreview({ start: event.start, end: event.end }); }}/>
        </>}
        <path d={path("parentParent")} className="plot-line parent" />
        <path d={path("recombinantMajor")} className="plot-line major" />
        <path d={path("recombinantMinor")} className="plot-line minor" />
        <text x="5" y={top + 5} className="axis-label">100%</text>
        <text x="5" y={top + innerHeight} className="axis-label">{Math.round(yMin * 100)}%</text>
        {[0, 0.25, 0.5, 0.75, 1].map((fraction) => (
          <text key={fraction} x={left + fraction * innerWidth} y={height - 8} textAnchor={fraction === 0 ? "start" : fraction === 1 ? "end" : "middle"} className="axis-label">
            {Math.round(fraction * alignment.length).toLocaleString("en-US")}
          </text>
        ))}
      </svg>
      <div className="plot-legend">
        <span><i className="legend-dot major"/>Recombinant ↔ major</span>
        <span><i className="legend-dot minor"/>Recombinant ↔ minor</span>
        <span><i className="legend-dot parent"/>Major ↔ minor</span>
        <span className="plot-hint">Drag the dashed breakpoint handles for live manual refinement</span>
      </div>
    </div>
  );
}

function Overview({ alignment, events, selectedId, onSelect }: {
  alignment: AlignmentData;
  events: RdpEvent[];
  selectedId: string | null;
  onSelect: (id: string) => void;
}) {
  const eventsByRecombinant = useMemo(() => {
    const index = new Map<number, RdpEvent[]>();
    events.forEach((event) => {
      if (event.decision === "rejected") return;
      const existing = index.get(event.recombinant);
      if (existing) existing.push(event);
      else index.set(event.recombinant, [event]);
    });
    return index;
  }, [events]);
  const displayedIndexes = useMemo(() => [...new Set([
    ...eventsByRecombinant.keys(),
    ...Array.from({ length: Math.min(400, alignment.sequences.length) }, (_, index) => index),
  ])].sort((left, right) => left - right), [alignment.sequences.length, eventsByRecombinant]);
  return (
    <div className="overview">
      <div className="overview-axis">
        <span>1</span><span>{Math.round(alignment.length * 0.25).toLocaleString("en-US")}</span><span>{Math.round(alignment.length * 0.5).toLocaleString("en-US")}</span><span>{Math.round(alignment.length * 0.75).toLocaleString("en-US")}</span><span>{alignment.length.toLocaleString("en-US")} nt</span>
      </div>
      <div className="overview-scroll">
        {displayedIndexes.map((index) => {
          const sequence = alignment.sequences[index];
          const sequenceEvents = eventsByRecombinant.get(index) ?? [];
          return (
            <div className="sequence-track" key={`${sequence.name}-${index}`}>
              <button className="sequence-name" type="button" title={sequence.name}>{sequence.name}</button>
              <div className="track-line">
                <div className="track-base" />
                {sequenceEvents.flatMap((event) => eventSegments(event, alignment.length).map(([start, end], segment) => (
                  <button
                    type="button"
                    key={`${event.id}-${segment}`}
                    className={`event-block ${event.id === selectedId ? "selected" : ""} ${event.decision} ${event.wraps ? "wraps" : ""}`}
                    style={{ left: `${(start / alignment.length) * 100}%`, width: `${Math.max(0.8, ((end - start) / alignment.length) * 100)}%` }}
                    onClick={() => onSelect(event.id)}
                    title={`Event ${formatEventRegion(event, alignment.length)}`}
                  ><span /></button>
                )))}
              </div>
              <span className={`role-dot ${sequence.role ?? "both"}`}>{sequence.role === "query" ? "Q" : sequence.role === "reference" ? "R" : "B"}</span>
            </div>
          );
        })}
        {displayedIndexes.length < alignment.sequences.length && <div className="overview-cap">Showing {displayedIndexes.length.toLocaleString("en-US")} sequences, including every sequence with a retained event. Use Sequence roles search for the rest.</div>}
      </div>
    </div>
  );
}

function EventTable({ alignment, events, selectedId, onSelect }: {
  alignment: AlignmentData;
  events: RdpEvent[];
  selectedId: string | null;
  onSelect: (id: string) => void;
}) {
  const [query, setQuery] = useState("");
  const [method, setMethod] = useState<MethodName | "all">("all");
  const [decision, setDecision] = useState<EventDecision | "all">("all");
  const filtered = events.filter((event) => {
    const names = [event.recombinant, event.majorParent, event.minorParent]
      .map((index) => alignment.sequences[index]?.name ?? "")
      .join(" ")
      .toLowerCase();
    const matchesQuery = !query || names.includes(query.toLowerCase()) || formatEventRegion(event, alignment.length).includes(query);
    const matchesMethod = method === "all" || event.evidence.some((item) => item.method === method && item.supported);
    const matchesDecision = decision === "all" || event.decision === decision;
    return matchesQuery && matchesMethod && matchesDecision;
  });
  return (
    <div>
      <div className="event-filters">
        <label><span>Find event</span><input value={query} onChange={(value) => setQuery(value.target.value)} placeholder="Sequence or coordinate…" /></label>
        <label><span>Supported by</span><select value={method} onChange={(value) => setMethod(value.target.value as MethodName | "all")}><option value="all">Any enabled method</option>{PRIMARY_METHODS.map((item) => <option value={item} key={item}>{item}</option>)}</select></label>
        <label><span>Review</span><select value={decision} onChange={(value) => setDecision(value.target.value as EventDecision | "all")}><option value="all">Any decision</option><option value="unreviewed">Unreviewed</option><option value="accepted">Accepted</option><option value="rejected">Rejected</option></select></label>
        <span>{filtered.length} of {events.length} hypotheses</span>
      </div>
      <div className="table-scroll">
      <table className="event-table">
        <thead><tr><th>Event</th><th>Recombinant</th><th>Parent assignment</th><th>Region</th><th>Supporting methods</th><th>Best adjusted p</th><th>Review</th></tr></thead>
        <tbody>
          {filtered.length === 0 ? <tr><td colSpan={7} className="empty-cell">{events.length ? "No hypotheses match these filters." : "No events pass the current consensus threshold."}</td></tr> : filtered.map((event) => {
            const index = events.findIndex((candidate) => candidate.id === event.id);
            const supported = event.evidence.filter((item) => item.supported);
            const bestP = event.evidence.length ? Math.min(...event.evidence.map((item) => item.correctedP)) : Number.NaN;
            return (
              <tr key={event.id} className={event.id === selectedId ? "selected" : ""} onClick={() => onSelect(event.id)} onKeyDown={(key) => { if (key.key === "Enter" || key.key === " ") { key.preventDefault(); onSelect(event.id); } }} tabIndex={0} aria-selected={event.id === selectedId}>
                <td><span className="event-number" title={event.groupId ? `Grouped as ${event.groupId}` : "Ungrouped event"}>{index + 1}{event.groupId && <i/>}</span><small className={`source-label ${event.source}`}>{event.source === "wasm" ? "scan" : event.source === "example" ? "known truth" : "manual"}</small></td>
                <td><b className="primary-name">{alignment.sequences[event.recombinant]?.name ?? "—"}</b></td>
                <td><div className="parent-pair"><span><i className="parent-role major"/>Major · {alignment.sequences[event.majorParent]?.name ?? "—"}</span><span><i className="parent-role minor"/>Minor · {alignment.sequences[event.minorParent]?.name ?? "—"}</span></div></td>
                <td className="mono">{formatEventRegion(event, alignment.length)}</td>
                <td><div className="method-chips">{supported.length ? supported.map((item) => <span key={item.method} title={`${item.method}: adjusted p ${formatP(item.correctedP)}`}>{item.method}</span>) : <em>{event.evidenceStale ? "needs calculation" : "none"}</em>}</div></td>
                <td className="mono">{formatP(bestP)}</td>
                <td><span className={`decision-label ${event.decision}`}>{event.decision}</span>{event.evidenceStale && <small className="stale-label">stale</small>}</td>
              </tr>
            );
          })}
        </tbody>
      </table>
      </div>
    </div>
  );
}

function MethodSpecificPlot({ alignment, event, method, window }: { alignment: AlignmentData; event: RdpEvent; method: MethodName; window: number }) {
  const points = useMemo(() => {
    const recombinant = alignment.sequences[event.recombinant].sequence;
    const major = alignment.sequences[event.majorParent].sequence;
    const minor = alignment.sequences[event.minorParent].sequence;
    const canonical = (base: string) => base === "A" || base === "C" || base === "G" || base === "T";
    const stride = Math.max(1, Math.ceil(alignment.length / 180));
    if (method === "3Seq") {
      const output: Array<{ position: number; value: number }> = [];
      let walk = 0;
      let maximum = 0;
      for (let site = 0; site < alignment.length; site += 1) {
        const r = recombinant[site];
        const a = major[site];
        const b = minor[site];
        if (canonical(r) && canonical(a) && canonical(b) && a !== b) {
          if (r === a) walk += 1;
          else if (r === b) walk -= 1;
          maximum = Math.max(maximum, walk);
        }
        if (site % stride === 0 || site === alignment.length - 1) output.push({ position: site, value: maximum - walk });
      }
      return { values: output, label: "HGRW descent from previous maximum", baseline: 0 };
    }
    const output: Array<{ position: number; value: number }> = [];
    const half = Math.max(10, Math.floor(window / 2));
    for (let center = 0; center < alignment.length; center += stride) {
      const start = Math.max(0, center - half);
      const end = Math.min(alignment.length, center + half);
      let majorMatches = 0;
      let minorMatches = 0;
      let discriminating = 0;
      let validSites = 0;
      let identityMajor = 0;
      let identityMinor = 0;
      let run = 0;
      let bestRun = 0;
      const sides = [[0, 0], [0, 0]];
      for (let site = start; site < end; site += 1) {
        const r = recombinant[site];
        const a = major[site];
        const b = minor[site];
        if (!canonical(r) || !canonical(a) || !canonical(b)) continue;
        validSites += 1;
        if (r === a) identityMajor += 1;
        if (r === b) identityMinor += 1;
        if (a === b) { run = 0; continue; }
        discriminating += 1;
        const side = site < center ? 0 : 1;
        if (r === a) {
          majorMatches += 1;
          sides[side][0] += 1;
          run = 0;
        } else if (r === b) {
          minorMatches += 1;
          sides[side][1] += 1;
          run += 1;
          bestRun = Math.max(bestRun, run);
        } else run = 0;
      }
      const [leftMajor, leftMinor] = sides[0];
      const [rightMajor, rightMinor] = sides[1];
      const total = leftMajor + leftMinor + rightMajor + rightMinor;
      const determinant = leftMajor * rightMinor - leftMinor * rightMajor;
      const denominator = Math.max(1, (leftMajor + leftMinor) * (rightMajor + rightMinor) * (leftMajor + rightMajor) * (leftMinor + rightMinor));
      const chi = total * determinant * determinant / denominator;
      const minorLeft = leftMinor / Math.max(1, leftMajor + leftMinor);
      const minorRight = rightMinor / Math.max(1, rightMajor + rightMinor);
      const value = method === "RDP"
        ? (minorMatches - majorMatches) / Math.max(1, discriminating)
        : method === "GENECONV"
          ? bestRun
          : method === "BootScan"
            ? (identityMinor - identityMajor) / Math.max(1, validSites)
            : method === "MaxChi"
              ? chi
              : method === "Chimaera"
                ? Math.abs(minorLeft - minorRight)
                : (minorMatches - majorMatches) / Math.sqrt(Math.max(1, discriminating));
      output.push({ position: center, value });
    }
    const labels: Record<Exclude<MethodName, "3Seq">, string> = {
      RDP: "Minor-minus-major match fraction at parent-discriminating sites",
      GENECONV: "Longest local minor-parent concordant run",
      BootScan: "Minor-minus-major window identity",
      MaxChi: "Local 2×2 boundary χ² profile",
      Chimaera: "Binary parent-state contrast across boundary",
      SiScan: "Oriented parent-category Z profile",
    };
    return { values: output, label: labels[method], baseline: method === "RDP" || method === "BootScan" || method === "SiScan" ? 0 : undefined };
  }, [alignment, event, method, window]);
  const width = 900;
  const height = 150;
  const left = 42;
  const right = 12;
  const top = 12;
  const bottom = 25;
  const minimum = Math.min(...points.values.map((point) => point.value), points.baseline ?? Number.POSITIVE_INFINITY);
  const maximum = Math.max(...points.values.map((point) => point.value), points.baseline ?? Number.NEGATIVE_INFINITY);
  const padding = Math.max(0.001, (maximum - minimum) * 0.08);
  const low = minimum - padding;
  const high = maximum + padding;
  const x = (position: number) => left + position / Math.max(1, alignment.length) * (width - left - right);
  const y = (value: number) => top + (high - value) / Math.max(0.001, high - low) * (height - top - bottom);
  const path = points.values.map((point, index) => `${index ? "L" : "M"}${x(point.position).toFixed(1)},${y(point.value).toFixed(1)}`).join(" ");
  return <div className="method-profile">
    <div><b>{method}-specific exploratory profile</b><span>{points.label}</span></div>
    <svg viewBox={`0 0 ${width} ${height}`} role="img" aria-label={`${method} profile across the alignment`}>
      {eventSegments(event, alignment.length).map(([start, end], index) => <rect key={index} x={x(start)} y={top} width={Math.max(1, x(end) - x(start))} height={height - top - bottom} className="method-profile-tract"/>)}
      {points.baseline !== undefined && <line x1={left} x2={width - right} y1={y(points.baseline)} y2={y(points.baseline)} className="method-profile-zero"/>}
      <path d={path} className="method-profile-line"/>
      <line x1={x(event.start)} x2={x(event.start)} y1={top} y2={height - bottom} className="breakpoint-line"/>
      <line x1={x(event.end)} x2={x(event.end)} y1={top} y2={height - bottom} className="breakpoint-line"/>
      <text x={4} y={top + 5} className="axis-label">{maximum.toPrecision(3)}</text>
      <text x={4} y={height - bottom} className="axis-label">{minimum.toPrecision(3)}</text>
      <text x={left} y={height - 7} className="axis-label">1</text>
      <text x={width - right} y={height - 7} textAnchor="end" className="axis-label">{alignment.length.toLocaleString("en-US")}</text>
    </svg>
    <p>This bounded profile is calculated only for the selected triplet and is intended for visual review. The saved p-value above comes from the worker statistic and its documented calibration.</p>
  </div>;
}

function MethodEvidencePanel({ alignment, event, window, selectedMethod, onSelectMethod }: {
  alignment: AlignmentData;
  event: RdpEvent;
  window: number;
  selectedMethod: MethodName;
  onSelectMethod: (method: MethodName) => void;
}) {
  const evidence = event.evidence.find((item) => item.method === selectedMethod);
  const supported = event.evidence.filter((item) => item.supported).length;
  return <div className="method-result-explorer">
    <div className="method-result-tabs" role="tablist" aria-label="Method results">
      {PRIMARY_METHODS.map((method) => {
        const item = event.evidence.find((candidate) => candidate.method === method);
        return <button type="button" role="tab" aria-selected={selectedMethod === method} key={method} className={`${selectedMethod === method ? "active" : ""} ${item?.supported ? "supported" : item ? "not-supported" : "not-run"}`} onClick={() => onSelectMethod(method)}><span>{item?.supported ? "✓" : item ? "·" : "—"}</span><b>{method}</b><small>{item ? formatP(item.correctedP) : "not run"}</small></button>;
      })}
    </div>
    <div className="method-result-detail">
      <div className="method-summary">
        <span className="eyebrow">{supported}/{event.evidence.length || 0} enabled methods support this hypothesis</span>
        <h3>{selectedMethod} result</h3>
        <p>{METHOD_META[selectedMethod].detail}</p>
      </div>
      {evidence ? <>
        <div className="method-numbers">
          <article><span>Decision at α</span><b className={evidence.supported ? "positive" : "negative"}>{evidence.supported ? "Supports" : "Does not support"}</b><small>{event.evidenceStale ? "Saved result is stale after edits" : "Current hypothesis"}</small></article>
          <article><span>Adjusted p</span><b>{formatP(evidence.correctedP)}</b><small>{evidence.correctionScope ?? "Saved scan correction"}</small></article>
          <article><span>Raw method p</span><b>{formatP(evidence.pValue)}</b><small>Before experiment-wide correction</small></article>
          <article><span>{evidence.statisticLabel}</span><b>{Number.isFinite(evidence.statistic) ? evidence.statistic.toPrecision(4) : "—"}</b><small>{evidence.calibration}</small></article>
        </div>
        {event.recalculationNote && <div className="method-scope-note">{event.recalculationNote}</div>}
      </> : <div className="method-not-run"><b>{selectedMethod} has no result for this hypothesis.</b><span>Enable the method and run a full scan, or recalculate the selected manual hypothesis.</span></div>}
      <div className="method-limit"><strong>Interpretation limit</strong><span>{METHOD_META[selectedMethod].limit}</span><a href={METHOD_META[selectedMethod].citation} target="_blank" rel="noreferrer">Primary paper ↗</a></div>
    </div>
    <MethodSpecificPlot alignment={alignment} event={event} method={selectedMethod} window={window}/>
  </div>;
}

function AlignmentViewer({ alignment, event }: { alignment: AlignmentData; event: RdpEvent | null }) {
  const initial = event ? Math.max(0, event.start - 45) : 0;
  const [cursor, setCursor] = useState(initial);
  const [mode, setMode] = useState<"bases" | "affinity">(event ? "affinity" : "bases");
  const [columnCount, setColumnCount] = useState<48 | 96 | 144>(96);
  const [informativeOnly, setInformativeOnly] = useState(Boolean(event));
  const [parentIndexes, setParentIndexes] = useState<number[]>(() => {
    const suggested = event ? [event.majorParent, event.minorParent] : [0, 1];
    return [...new Set(suggested)].filter((index) => index >= 0 && index < alignment.sequences.length).slice(0, 6);
  });
  const parentColors = ["#39b99b", "#ff7655", "#a778d2", "#4c8ed9", "#d69b35", "#65a857"];
  const informativeSupported = alignment.length <= 2_000_000;
  const informativeSites = useMemo(() => mode === "affinity" && informativeOnly && informativeSupported
    ? parentInformativeSites(alignment.sequences.map((record) => record.sequence), parentIndexes, alignment.length)
    : null, [alignment.length, alignment.sequences, informativeOnly, informativeSupported, mode, parentIndexes]);
  const lowerBound = (values: number[], target: number) => {
    let left = 0;
    let right = values.length;
    while (left < right) {
      const middle = (left + right) >>> 1;
      if (values[middle] < target) left = middle + 1;
      else right = middle;
    }
    return left;
  };
  const informativeStart = informativeSites ? Math.min(lowerBound(informativeSites, cursor), Math.max(0, informativeSites.length - 1)) : 0;
  const shownSites = informativeSites
    ? informativeSites.slice(informativeStart, informativeStart + columnCount)
    : Array.from({ length: Math.min(columnCount, alignment.length - cursor) }, (_, offset) => cursor + offset);
  const visibleLimit = Math.min(250, alignment.sequences.length);
  const suggestedIndexes = event ? [event.recombinant, event.majorParent, event.minorParent] : [];
  const indexes = [...new Set([...parentIndexes, ...suggestedIndexes, ...Array.from({ length: visibleLimit }, (_, index) => index)])];
  const firstSite = shownSites[0] ?? cursor;
  const lastSite = shownSites.at(-1) ?? cursor;
  const parentSlotBySequence = new Map(parentIndexes.map((sequenceIndex, parentSlot) => [sequenceIndex, parentSlot]));
  const baseClass = (base: string) => `base base-${base === "-" ? "gap" : "ACGT".includes(base) ? base : "amb"}`;
  const move = (direction: -1 | 1) => {
    if (informativeSites?.length) {
      const next = Math.max(0, Math.min(informativeSites.length - 1, informativeStart + direction * columnCount));
      setCursor(informativeSites[next]);
      return;
    }
    setCursor((current) => Math.max(0, Math.min(Math.max(0, alignment.length - columnCount), current + direction * columnCount)));
  };
  const setParent = (sequenceIndex: number) => {
    setParentIndexes((current) => {
      if (current.includes(sequenceIndex)) return current.length > 2 ? current.filter((index) => index !== sequenceIndex) : current;
      return current.length < 6 ? [...current, sequenceIndex] : [...current.slice(1), sequenceIndex];
    });
  };
  const resetSuggestedParents = () => {
    if (!event) return;
    setParentIndexes([...new Set([event.majorParent, event.minorParent])]);
    setMode("affinity");
    setInformativeOnly(true);
  };
  return (
    <Panel title="Nucleotide alignment" action={<span className="panel-caption">{informativeSites ? `${shownSites.length} parent-informative columns · ` : "Sites "}{formatInteger(firstSite + 1)}–{formatInteger(lastSite + 1)}</span>}>
      <div className="alignment-modebar">
        <div><span className="control-label">Display</span><Segmented value={mode} options={[{ value: "bases", label: "Base colors" }, { value: "affinity", label: "Parent highlighter" }]} onChange={setMode}/></div>
        <label className="alignment-check"><input type="checkbox" checked={informativeOnly} disabled={mode !== "affinity" || !informativeSupported} onChange={(value) => setInformativeOnly(value.target.checked)}/><span>Parent-informative sites only</span></label>
        <label className="alignment-zoom"><span>Columns</span><select value={columnCount} onChange={(value) => setColumnCount(Number(value.target.value) as 48 | 96 | 144)}><option value={48}>48 · large</option><option value={96}>96 · medium</option><option value={144}>144 · compact</option></select></label>
        {event && <button type="button" className="small-button" onClick={resetSuggestedParents}>Use inferred parents</button>}
      </div>
      <div className="parent-picker">
        <div className="parent-picker-copy"><b>Selected parents</b><span>Click the star beside any sequence to add or remove it. Two to six parents can be compared.</span></div>
        <div className="parent-chips">{parentIndexes.map((sequenceIndex, parentSlot) => <button type="button" key={sequenceIndex} style={{ "--parent-color": parentColors[parentSlot] } as React.CSSProperties} onClick={() => setParent(sequenceIndex)} disabled={parentIndexes.length <= 2} title={parentIndexes.length <= 2 ? "At least two parents are required" : "Remove parent"}><i/>{alignment.sequences[sequenceIndex]?.name ?? `Sequence ${sequenceIndex + 1}`}<span>×</span></button>)}</div>
        <select aria-label="Add a parent sequence" value="" disabled={parentIndexes.length >= 6} onChange={(value) => { if (value.target.value) setParent(Number(value.target.value)); }}><option value="">＋ Add parent…</option>{alignment.sequences.map((record, index) => parentIndexes.includes(index) ? null : <option key={index} value={index}>{record.name}</option>)}</select>
      </div>
      {mode === "affinity" && <div className="affinity-legend"><span className="legend-intro">Target cell matches:</span>{parentIndexes.map((sequenceIndex, parentSlot) => <span key={sequenceIndex}><i style={{ background: parentColors[parentSlot] }}/>{alignment.sequences[sequenceIndex].name}</span>)}<span><i className="shared"/>multiple parents</span><span><i className="novel"/>none / novel</span><span><i className="missing"/>gap / ambiguous</span></div>}
      {!informativeSupported && mode === "affinity" && <div className="alignment-large-note">Informative-only indexing is disabled above 2,000,000 sites to keep memory bounded; highlighter colors remain available in coordinate windows.</div>}
      <div className="alignment-toolbar">
        <button type="button" className="small-button" onClick={() => move(-1)}>← Previous</button>
        <input aria-label="Alignment position" type="range" min={0} max={Math.max(0, alignment.length - 1)} value={Math.min(cursor, Math.max(0, alignment.length - 1))} onChange={(eventValue) => setCursor(Number(eventValue.target.value))} />
        <button type="button" className="small-button" onClick={() => move(1)}>Next →</button>
      </div>
      <div className={`alignment-grid ${mode === "affinity" ? "affinity-mode" : ""}`} role="region" aria-label="Scrollable sequence alignment">
        <div className="alignment-ruler" style={{ gridTemplateColumns: `190px repeat(${shownSites.length}, 10px)`, minWidth: `${190 + shownSites.length * 10}px` }}><span />{shownSites.map((site) => <i key={site} title={`Site ${formatInteger(site + 1)}`}>{(site + 1) % 10 === 0 ? "·" : ""}</i>)}</div>
        {informativeSites?.length === 0 && <div className="alignment-no-sites"><b>No parent-informative columns</b><span>The selected parents have no callable A/C/G/T differences. Add a more divergent parent or turn off informative-only filtering.</span></div>}
        {indexes.map((sequenceIndex) => {
          const record = alignment.sequences[sequenceIndex];
          const relation = event?.recombinant === sequenceIndex ? "Recombinant" : event?.majorParent === sequenceIndex ? "Major parent" : event?.minorParent === sequenceIndex ? "Minor parent" : "";
          const parentSlot = parentSlotBySequence.get(sequenceIndex);
          return (
            <div className={`alignment-row ${relation ? "relevant" : ""} ${parentSlot !== undefined ? "parent-source" : ""}`} style={{ gridTemplateColumns: `190px repeat(${shownSites.length}, 10px)`, minWidth: `${190 + shownSites.length * 10}px`, borderLeftColor: parentSlot !== undefined ? parentColors[parentSlot] : "transparent" }} key={`${record.name}-${sequenceIndex}`}>
              <div className="alignment-label" title={`${record.name}${relation ? ` · ${relation}` : ""}`}><button type="button" className={parentSlot !== undefined ? "parent-toggle active" : "parent-toggle"} style={parentSlot !== undefined ? { color: parentColors[parentSlot] } : undefined} aria-pressed={parentSlot !== undefined} title={parentSlot !== undefined && parentIndexes.length <= 2 ? "At least two parents are required" : parentSlot !== undefined ? "Remove from parent set" : "Add to parent set"} aria-label={`${parentSlot !== undefined ? "Remove" : "Use"} ${record.name} ${parentSlot !== undefined ? "from" : "as"} parent set`} onClick={() => setParent(sequenceIndex)}><Icon name="star" size={12}/></button><span><b>{record.name}</b><small>{parentSlot !== undefined ? `Parent ${parentSlot + 1}${relation ? ` · ${relation}` : ""}` : relation}</small></span></div>
              {shownSites.map((site) => {
                const base = record.sequence[site];
                if (mode !== "affinity" || parentSlot !== undefined) return <i key={site} className={`${baseClass(base)} ${parentSlot !== undefined ? `parent-base parent-${parentSlot}` : ""}`} style={parentSlot !== undefined ? { background: parentColors[parentSlot] } : undefined} title={`${record.name} · site ${formatInteger(site + 1)} · ${base}`}>{base}</i>;
                const affinity = classifyParentAffinity(base, parentIndexes.map((index) => alignment.sequences[index]?.sequence[site]));
                const uniqueSlot = affinity.kind === "unique" ? affinity.parentSlots[0] : null;
                return <i key={site} className={`base affinity-${affinity.kind}${uniqueSlot !== null ? ` affinity-parent-${uniqueSlot}` : ""}`} style={uniqueSlot !== null ? { background: parentColors[uniqueSlot] } : undefined} title={`${record.name} · site ${formatInteger(site + 1)} · ${base} · ${affinityDescription(affinity, parentIndexes.map((index) => alignment.sequences[index].name))}`}>{base}</i>;
              })}
            </div>
          );
        })}
        {indexes.length < alignment.sequences.length && <div className="alignment-cap">Showing {formatInteger(indexes.length)} of {formatInteger(alignment.sequences.length)} sequences. Parent selections and event sequences are always retained; use the role filter below for the full collection.</div>}
      </div>
    </Panel>
  );
}

function AnnotationPanel({ alignment, onOpen }: { alignment: AlignmentData; onOpen: () => void }) {
  const features = alignment.features ?? [];
  return <Panel title="Genome annotations" action={<div className="annotation-actions"><button type="button" className="small-button" onClick={onOpen}>Open GFF / GenBank / BED</button>{features.length > 0 && <button type="button" className="small-button" onClick={() => download("rdp-web-annotations.gff3", toGff3(features, alignment.name))}>GFF3 ↓</button>}</div>}>
    {features.length === 0 ? <div className="empty-state compact">Import annotations to map genes, CDS features, and named regions onto event coordinates.</div> : <>
      <div className="annotation-track">{features.slice(0, 2_000).map((feature, index) => <button type="button" key={feature.id} title={`${feature.name} · ${feature.type} · ${feature.start + 1}–${feature.end} (${feature.strand})`} style={{ left: `${feature.start / alignment.length * 100}%`, width: `${Math.max(0.5, (feature.end - feature.start) / alignment.length * 100)}%`, top: `${8 + (index % 4) * 20}px` }}><span>{feature.name}</span></button>)}</div>
      <div className="annotation-list">{features.slice(0, 200).map((feature) => <div key={`row-${feature.id}`}><b>{feature.name}</b><span>{feature.type}</span><code>{feature.start + 1}–{feature.end}{feature.strand !== "." ? ` · ${feature.strand}` : ""}</code></div>)}</div>
      {features.length > 200 && <p className="annotation-overflow">Showing the first 200 of {features.length.toLocaleString("en-US")} features.</p>}
    </>}
  </Panel>;
}

interface HeatmapLayer {
  values: Float32Array;
  label: string;
  minimum: number;
  maximum: number;
  colors: string[];
  logarithmic?: boolean;
  format: (value: number) => string;
}

const COUNT_COLORS = ["#f7faf8", "#dcebe7", "#9bcfc3", "#4e9d97", "#356a82", "#26345d"];
const DISCORDANCE_COLORS = ["#f7f9f8", "#e2e2ec", "#b9b3d5", "#8179b4", "#594789", "#382757"];
const RESIDUAL_COLORS = ["#b64d72", "#e59aaa", "#f1d8d8", "#eef1ef", "#b8ded8", "#5ea99f", "#176d78"];

function parseHexColor(color: string): [number, number, number] {
  return [Number.parseInt(color.slice(1, 3), 16), Number.parseInt(color.slice(3, 5), 16), Number.parseInt(color.slice(5, 7), 16)];
}

function interpolateColors(colors: string[], fraction: number): string {
  const bounded = Math.max(0, Math.min(1, fraction));
  const position = bounded * (colors.length - 1);
  const leftIndex = Math.min(colors.length - 1, Math.floor(position));
  const rightIndex = Math.min(colors.length - 1, leftIndex + 1);
  const local = position - leftIndex;
  const left = parseHexColor(colors[leftIndex]);
  const right = parseHexColor(colors[rightIndex]);
  const channels = left.map((channel, index) => Math.round(channel + (right[index] - channel) * local));
  return `rgb(${channels[0]} ${channels[1]} ${channels[2]})`;
}

function layerColor(layer: HeatmapLayer, value: number): string {
  const span = Math.max(1e-12, layer.maximum - layer.minimum);
  const linear = Math.max(0, Math.min(1, (value - layer.minimum) / span));
  const normalized = layer.logarithmic && layer.minimum === 0
    ? Math.log1p(Math.max(0, value)) / Math.log1p(Math.max(1e-12, layer.maximum))
    : linear;
  return interpolateColors(layer.colors, normalized);
}

function matrixBinLabel(bin: number, resolution: number, length: number): string {
  const start = Math.floor((bin / resolution) * length) + 1;
  const end = Math.max(start, Math.floor(((bin + 1) / resolution) * length));
  return `${formatInteger(start)}–${formatInteger(end)}`;
}

function GenomePositionHeatmap({
  resolution,
  length,
  upperLayer,
  lowerLayer,
  ariaLabel,
  onActivate,
}: {
  resolution: number;
  length: number;
  upperLayer: HeatmapLayer;
  lowerLayer?: HeatmapLayer;
  ariaLabel: string;
  onActivate?: (row: number, column: number) => void;
}) {
  const canvasRef = useRef<HTMLCanvasElement>(null);
  const [activeCell, setActiveCell] = useState<{ row: number; column: number } | null>(null);
  useEffect(() => {
    const canvas = canvasRef.current;
    if (!canvas) return;
    const draw = () => {
      const logicalSize = Math.max(420, Math.round(canvas.getBoundingClientRect().width || 640));
      const pixelRatio = Math.max(1, Math.min(2, window.devicePixelRatio || 1));
      canvas.width = Math.round(logicalSize * pixelRatio);
      canvas.height = Math.round(logicalSize * pixelRatio);
      const context = canvas.getContext("2d", { alpha: false });
      if (!context) return;
      context.setTransform(pixelRatio, 0, 0, pixelRatio, 0, 0);
      context.fillStyle = "#ffffff";
      context.fillRect(0, 0, logicalSize, logicalSize);
      const cellSize = logicalSize / resolution;
      for (let row = 0; row < resolution; row += 1) {
        for (let column = 0; column < resolution; column += 1) {
          if (lowerLayer && row === column) {
            context.fillStyle = "#ffffff";
          } else {
            const layer = lowerLayer && row > column ? lowerLayer : upperLayer;
            context.fillStyle = layerColor(layer, layer.values[row * resolution + column] ?? 0);
          }
          context.fillRect(column * cellSize, row * cellSize, Math.ceil(cellSize + 0.15), Math.ceil(cellSize + 0.15));
        }
      }
      if (lowerLayer) {
        context.strokeStyle = "rgb(255 255 255 / 92%)";
        context.lineWidth = Math.max(2, logicalSize / 180);
        context.beginPath();
        context.moveTo(0, 0);
        context.lineTo(logicalSize, logicalSize);
        context.stroke();
      }
    };
    draw();
    const observer = typeof ResizeObserver === "undefined" ? null : new ResizeObserver(draw);
    observer?.observe(canvas);
    window.addEventListener("resize", draw);
    return () => {
      observer?.disconnect();
      window.removeEventListener("resize", draw);
    };
  }, [lowerLayer, resolution, upperLayer]);

  const updatePointer = (clientX: number, clientY: number) => {
    const canvas = canvasRef.current;
    if (!canvas) return;
    const bounds = canvas.getBoundingClientRect();
    const column = Math.max(0, Math.min(resolution - 1, Math.floor(((clientX - bounds.left) / bounds.width) * resolution)));
    const row = Math.max(0, Math.min(resolution - 1, Math.floor(((clientY - bounds.top) / bounds.height) * resolution)));
    setActiveCell({ row, column });
  };
  const activeLayer = activeCell
    ? lowerLayer && activeCell.row > activeCell.column ? lowerLayer : upperLayer
    : null;
  const activeValue = activeCell && activeLayer ? activeLayer.values[activeCell.row * resolution + activeCell.column] ?? 0 : null;
  const ticks = [0, 0.25, 0.5, 0.75, 1].map((fraction) => formatInteger(Math.max(1, Math.round(fraction * length))));
  return <div className="genome-heatmap">
    <div className="genome-heatmap-axis top" aria-hidden="true">{ticks.map((tick, index) => <span key={`${tick}-${index}`}>{tick}</span>)}</div>
    <div className="genome-heatmap-y-axis" aria-hidden="true">{ticks.map((tick, index) => <span key={`${tick}-${index}`}>{tick}</span>)}</div>
    <div className="genome-heatmap-canvas-wrap">
      <canvas
        ref={canvasRef}
        className="genome-heatmap-canvas"
        role="img"
        tabIndex={0}
        aria-label={ariaLabel}
        onPointerMove={(event) => updatePointer(event.clientX, event.clientY)}
        onPointerLeave={() => setActiveCell(null)}
        onClick={() => activeCell && onActivate?.(activeCell.row, activeCell.column)}
        onFocus={() => setActiveCell((current) => current ?? { row: Math.floor(resolution / 2), column: Math.floor(resolution / 2) })}
        onKeyDown={(event) => {
          const current = activeCell ?? { row: Math.floor(resolution / 2), column: Math.floor(resolution / 2) };
          if (event.key === "Enter" || event.key === " ") {
            event.preventDefault();
            onActivate?.(current.row, current.column);
            return;
          }
          const movement: Record<string, [number, number]> = { ArrowUp: [-1, 0], ArrowDown: [1, 0], ArrowLeft: [0, -1], ArrowRight: [0, 1] };
          const delta = movement[event.key];
          if (!delta) return;
          event.preventDefault();
          setActiveCell({ row: Math.max(0, Math.min(resolution - 1, current.row + delta[0])), column: Math.max(0, Math.min(resolution - 1, current.column + delta[1])) });
        }}
      >Genome-position matrix. Use a browser with canvas support.</canvas>
      {activeCell && <><i className="matrix-crosshair horizontal" style={{ top: `${((activeCell.row + 0.5) / resolution) * 100}%` }}/><i className="matrix-crosshair vertical" style={{ left: `${((activeCell.column + 0.5) / resolution) * 100}%` }}/></>}
    </div>
    <div className="genome-heatmap-readout" aria-live="polite">{activeCell && activeLayer && activeValue !== null ? <><b>{activeLayer.label}</b><span>x {matrixBinLabel(activeCell.column, resolution, length)}</span><span>y {matrixBinLabel(activeCell.row, resolution, length)}</span><code>{activeLayer.format(activeValue)}</code>{onActivate && <em>Click / Enter to select a matching event</em>}</> : <><b>Inspect the matrix</b><span>Hover, focus, or use arrow keys for exact window coordinates.</span></>}</div>
    <div className={`genome-heatmap-legends ${lowerLayer ? "split" : ""}`}>
      <div><span>{upperLayer.label}</span><i style={{ background: `linear-gradient(90deg, ${upperLayer.colors.join(", ")})` }}/><small>{upperLayer.format(upperLayer.minimum)}</small><small>{upperLayer.format(upperLayer.maximum)}</small></div>
      {lowerLayer && <div><span>{lowerLayer.label}</span><i style={{ background: `linear-gradient(90deg, ${lowerLayer.colors.join(", ")})` }}/><small>{lowerLayer.format(lowerLayer.minimum)}</small><small>{lowerLayer.format(lowerLayer.maximum)}</small></div>}
    </div>
  </div>;
}

function RdpPatternMatrices({ alignment, events, onSelect }: { alignment: AlignmentData; events: RdpEvent[]; onSelect: (id: string) => void }) {
  const [resolution, setResolution] = useState<48 | 64 | 96>(96);
  const [scope, setScope] = useState<"retained" | "accepted">("retained");
  const [cohortSize, setCohortSize] = useState<12 | 18 | 24>(18);
  const scopedEvents = useMemo(() => events.filter((event) => event.decision !== "rejected" && (scope === "retained" || (event.decision === "accepted" && !event.evidenceStale))), [events, scope]);
  const effectiveResolution = scopedEvents.length > 2_000 ? Math.min(48, resolution) : scopedEvents.length > 500 ? Math.min(64, resolution) : resolution;
  const breakpointPairs = useMemo(() => computeBreakpointPairDensity(scopedEvents, alignment.length, effectiveResolution), [alignment.length, effectiveResolution, scopedEvents]);
  const regionSeparation = useMemo(() => computeRegionSeparationMatrices(scopedEvents, alignment.length, effectiveResolution), [alignment.length, effectiveResolution, scopedEvents]);
  const localDiscordance = useMemo(() => computeLocalDiscordanceMatrices(alignment.sequences.map((record) => record.sequence), alignment.length, effectiveResolution, cohortSize), [alignment.length, alignment.sequences, cohortSize, effectiveResolution]);
  const boundaryBin = (position: number) => Math.max(0, Math.min(effectiveResolution - 1, Math.floor((Math.max(0, Math.min(alignment.length - 1, position)) / Math.max(1, alignment.length)) * effectiveResolution)));
  const selectBreakpointPair = (row: number, column: number) => {
    const match = scopedEvents.find((event) => {
      const start = boundaryBin(event.start);
      const end = boundaryBin(event.end >= alignment.length ? alignment.length - 1 : event.end);
      return (start === row && end === column) || (start === column && end === row);
    });
    if (match) onSelect(match.id);
  };
  const selectSeparatingEvent = (row: number, column: number) => {
    const rowPosition = ((row + 0.5) / effectiveResolution) * alignment.length;
    const columnPosition = ((column + 0.5) / effectiveResolution) * alignment.length;
    const match = scopedEvents.find((event) => eventContainsPosition(event, rowPosition) !== eventContainsPosition(event, columnPosition));
    if (match) onSelect(match.id);
  };
  const integer = (value: number) => `${formatInteger(Math.round(value))} event${Math.round(value) === 1 ? "" : "s"}`;
  const signed = (value: number) => `${value >= 0 ? "+" : ""}${value.toFixed(2)} z`;
  return <>
    <div className="matrix-suite-header">
      <div><span className="eyebrow">RDP-style genome pattern matrices</span><h1>Read recombination across genomic coordinates—not event IDs.</h1><p>The matrix semantics and triangular layouts follow the dense views in RDP4 Figure 2. The palettes are perceptually ordered and avoid the original rainbow scale.</p></div>
      <div className="matrix-suite-controls">
        <label><span>Event scope</span><select value={scope} onChange={(event) => setScope(event.target.value as "retained" | "accepted")}><option value="retained">All retained hypotheses</option><option value="accepted">Accepted + fresh only</option></select></label>
        <label><span>Matrix bins</span><select value={resolution} onChange={(event) => setResolution(Number(event.target.value) as 48 | 64 | 96)}><option value={48}>48 · fast overview</option><option value={64}>64 · balanced</option><option value={96}>96 · high density</option></select></label>
        <label><span>Distance cohort</span><select value={cohortSize} onChange={(event) => setCohortSize(Number(event.target.value) as 12 | 18 | 24)}><option value={12}>12 sequences</option><option value={18}>18 sequences</option><option value={24}>24 sequences</option></select></label>
        <a href="https://academic.oup.com/ve/article/1/1/vev003/2568683" target="_blank" rel="noreferrer">RDP4 Figure 2 ↗</a>
      </div>
      <div className="matrix-suite-status"><span><b>{formatInteger(scopedEvents.length)}</b> events aggregated</span><span><b>{effectiveResolution} × {effectiveResolution}</b> genomic bins</span><span><b>{formatInteger(localDiscordance.pairCount)}</b> sequence pairs / window</span>{effectiveResolution !== resolution && <span className="limited"><b>Adaptive cap</b> resolution reduced for interactive speed</span>}</div>
    </div>
    <div className="matrix-suite-grid">
      <Panel title="Breakpoint-pair density" action={<span className="panel-caption">RDP4 Fig. 2c semantics · symmetric endpoint counts</span>}>
        <div className="matrix-card-body">
          <GenomePositionHeatmap resolution={effectiveResolution} length={alignment.length} ariaLabel="Breakpoint-pair counts by genomic position on both axes" upperLayer={{ values: breakpointPairs.values, label: "Breakpoint-pair count", minimum: 0, maximum: Math.max(1, breakpointPairs.maximum), colors: COUNT_COLORS, logarithmic: true, format: integer }} onActivate={selectBreakpointPair}/>
          <div className="matrix-interpretation"><b>Paired hotspots become symmetric islands.</b><p>Every event contributes its inferred start/end bin and the mirrored end/start bin. Empty space stays quiet; single events remain visible through logarithmic colour scaling.</p></div>
        </div>
      </Panel>
      <Panel title="Recombination region separation" action={<span className="panel-caption">Observed counts above · tract-placement residuals below</span>}>
        <div className="matrix-card-body">
          <GenomePositionHeatmap resolution={effectiveResolution} length={alignment.length} ariaLabel="Split triangular recombination region matrix by genomic position" upperLayer={{ values: regionSeparation.observed, label: "Upper · events separating windows", minimum: 0, maximum: Math.max(1, regionSeparation.maximumObserved), colors: COUNT_COLORS, logarithmic: true, format: integer }} lowerLayer={{ values: regionSeparation.standardizedResidual, label: "Lower · excess / deficit vs random tracts", minimum: -Math.max(2.5, Math.min(6, regionSeparation.maximumAbsoluteResidual)), maximum: Math.max(2.5, Math.min(6, regionSeparation.maximumAbsoluteResidual)), colors: RESIDUAL_COLORS, format: signed }} onActivate={selectSeparatingEvent}/>
          <div className="matrix-interpretation"><b>Upper triangle is observed; lower triangle is comparative.</b><p>The lower half uses a circular random-placement null that preserves every tract length. Teal means more separation than expected; rose means less. It is a fast analytical analogue, not a claim of numerical RDP4 permutation parity.</p></div>
        </div>
      </Panel>
      <Panel title="Local distance-profile discordance" className="matrix-wide" action={<span className="panel-caption">Fast PDDM-style phylogenetic signal proxy</span>}>
        <div className="matrix-card-body wide">
          <GenomePositionHeatmap resolution={effectiveResolution} length={alignment.length} ariaLabel="Split triangular local distance-profile discordance matrix by genomic position" upperLayer={{ values: localDiscordance.rmsDeviation, label: "Upper · RMS p-distance deviation", minimum: 0, maximum: Math.max(0.001, localDiscordance.maximumRmsDeviation), colors: DISCORDANCE_COLORS, format: (value) => value.toFixed(4) }} lowerLayer={{ values: localDiscordance.correlationLoss, label: "Lower · 1 − profile correlation", minimum: 0, maximum: Math.max(0.05, localDiscordance.maximumCorrelationLoss), colors: DISCORDANCE_COLORS, format: (value) => value.toFixed(3) }}/>
          <div className="matrix-interpretation"><b>Blocks expose regions with internally consistent but mutually discordant ancestry signals.</b><p>Each genomic window is represented by {formatInteger(localDiscordance.pairCount)} sampled pairwise distance values from {localDiscordance.sequenceIndexes.length} evenly distributed sequences and at most {localDiscordance.sampledSitesPerWindow} sites. RDP4 Figure 2e used SH and Robinson–Foulds matrices; this browser-fast view deliberately labels its distance-profile proxies.</p></div>
        </div>
      </Panel>
    </div>
  </>;
}

function DistanceMatrix({ alignment, matrix }: { alignment: AlignmentData; matrix: number[] }) {
  const sourceDimension = Number.isInteger(Math.sqrt(matrix.length)) ? Math.sqrt(matrix.length) : 0;
  const count = Math.min(64, alignment.sequences.length, sourceDimension || 32);
  const sourceStride = matrix.length === alignment.sequences.length ** 2 ? alignment.sequences.length : sourceDimension;
  const valueAt = (row: number, column: number) => {
    if (sourceStride >= count) return matrix[row * sourceStride + column] ?? 0;
    return 1 - pairwiseIdentitySampled(alignment.sequences[row].sequence, alignment.sequences[column].sequence);
  };
  const values = Array.from({ length: count ** 2 }, (_, flatIndex) => valueAt(Math.floor(flatIndex / count), flatIndex % count));
  const max = Math.max(0.001, ...values);
  return (
    <div className="matrix-layout">
      <div className="heatmap" style={{ gridTemplateColumns: `110px repeat(${count}, 11px)` }}>
        <span />
        {alignment.sequences.slice(0, count).map((record, index) => <span className="column-label" key={record.name} title={record.name}>{index % Math.max(1, Math.ceil(count / 32)) === 0 ? index + 1 : ""}</span>)}
        {alignment.sequences.slice(0, count).flatMap((record, row) => [
          <span className="row-label" key={`label-${record.name}`} title={record.name}>{record.name}</span>,
          ...alignment.sequences.slice(0, count).map((_, column) => {
            const value = values[row * count + column] ?? 0;
            const intensity = value / max;
            return <span key={`${row}-${column}`} className="heat-cell" title={`${record.name} × ${alignment.sequences[column].name}: ${(value * 100).toFixed(2)}% distance`} style={{ background: `color-mix(in srgb, #ff7050 ${Math.round(intensity * 88)}%, #f2f4f3)` }} />;
          }),
        ])}
      </div>
      <div className="matrix-key"><span>{formatInteger(count)} sequences · identical</span><i /><span>{(max * 100).toFixed(1)}% divergent</span></div>
    </div>
  );
}

function TopologyCards({ alignment, event }: { alignment: AlignmentData; event: RdpEvent | null }) {
  if (!event) return <div className="empty-state compact">Select an event to compare local topologies.</div>;
  const indexes = [event.recombinant, event.majorParent, event.minorParent];
  const labels = indexes.map((index) => alignment.sequences[index].name);
  const regions = event.wraps
    ? [
        { label: "Circular flank", segments: [[event.end, event.start] as [number, number]] },
        { label: "Origin-spanning tract", segments: eventSegments(event, alignment.length) },
      ]
    : [
        { label: "Left flank", segments: [[0, event.start] as [number, number]] },
        { label: "Recombinant tract", segments: [[event.start, event.end] as [number, number]] },
        { label: "Right flank", segments: [[event.end, alignment.length] as [number, number]] },
      ];
  return (
    <div className="topology-cards">
      {regions.map((region) => {
        const pairs = [
          { a: 0, b: 1, identity: pairwiseIdentitySegments(alignment.sequences[indexes[0]].sequence, alignment.sequences[indexes[1]].sequence, region.segments) },
          { a: 0, b: 2, identity: pairwiseIdentitySegments(alignment.sequences[indexes[0]].sequence, alignment.sequences[indexes[2]].sequence, region.segments) },
          { a: 1, b: 2, identity: pairwiseIdentitySegments(alignment.sequences[indexes[1]].sequence, alignment.sequences[indexes[2]].sequence, region.segments) },
        ].sort((left, right) => right.identity - left.identity);
        return (
          <article key={region.label} className={region.label.includes("tract") ? "highlight" : ""}>
            <span>{region.label}</span>
            <div className="mini-tree"><i /><b>{labels[pairs[0].a]}</b><b>{labels[pairs[0].b]}</b><em>{labels[3 - pairs[0].a - pairs[0].b]}</em></div>
            <small>{(pairs[0].identity * 100).toFixed(1)}% closest-pair identity</small>
          </article>
        );
      })}
    </div>
  );
}

function ChallengeDiagnostics({ diagnostics, event }: { diagnostics?: AlignmentDiagnostics; event: RdpEvent | null }) {
  if (!diagnostics) return <div className="empty-state compact">Run a scan to calculate bounded four-gamete and rate-variation diagnostics.</div>;
  const values = [
    ["Four-gamete incompatible", `${(diagnostics.fourGameteFraction * 100).toFixed(1)}%`, `${diagnostics.incompatibleSitePairs.toLocaleString("en-US")} / ${diagnostics.testedSitePairs.toLocaleString("en-US")} sampled site pairs`],
    ["Proximity permutation", formatP(diagnostics.proximityPermutationP), `One-sided four-gamete distance-clustering test · Δ ${(diagnostics.proximityStatistic * 100).toFixed(2)} points · ${diagnostics.proximityPermutationReplicates} seeded permutations`],
    ["Near / far contrast", diagnostics.proximityRatio.toFixed(2), "Incompatibility ratio; related to PHI’s question but not a PHI implementation"],
    ["Canonical coverage", `${((1 - diagnostics.ambiguityFraction) * 100).toFixed(2)}%`, `${diagnostics.sampledSequences} sequences × ${diagnostics.sampledBiallelicSites} biallelic sites sampled`],
    ["Tract rate ratio", event ? event.diagnostics.rateRatio.toFixed(2) : "—", event ? `${(event.diagnostics.tractVariableDensity * 100).toFixed(1)}% tract vs ${(event.diagnostics.backgroundVariableDensity * 100).toFixed(1)}% background variable sites` : "Select an event"],
    ["Parent-conflict rate", event ? `${(event.diagnostics.parentConflictRate * 100).toFixed(1)}%` : "—", event ? `${event.diagnostics.parentDiscriminatingSites.toLocaleString("en-US")} parent-discriminating tract sites` : "Select an event"],
  ];
  return <div className="diagnostic-grid">{values.map(([label, value, note]) => <article key={label}><span>{label}</span><b>{value}</b><small>{note}</small></article>)}</div>;
}

function TreeSvg({ root, roles = {}, markedNames = [], onToggleMark }: { root: NeighborJoiningNode; roles?: Record<string, "recombinant" | "major" | "minor">; markedNames?: string[]; onToggleMark?: (name: string) => void }) {
  const layout = useMemo(() => layoutNeighborJoiningTree(root), [root]);
  return <div className="tree-svg-wrap">
    <svg viewBox={`0 0 ${layout.width} ${layout.height}`} role="img" aria-label="Connected neighbor-joining tree">
      <path className="tree-branches" d={layout.path}/>
      {layout.joints.map((joint, index) => <circle className="tree-joint" key={`${joint.x}-${joint.y}-${index}`} cx={joint.x} cy={joint.y} r="1.8"/>)}
      {layout.labels.map((label) => <g key={`${label.name}-${label.y}`} className={markedNames.includes(label.name) ? "tree-leaf-marked" : undefined}>{markedNames.includes(label.name) && <circle className="tree-mark" cx={label.x - 4} cy={label.y} r="4"/>}<text role={onToggleMark ? "button" : undefined} tabIndex={onToggleMark ? 0 : undefined} aria-label={onToggleMark ? `${markedNames.includes(label.name) ? "Unmark" : "Mark"} ${label.name} across every tree` : undefined} className={`${roles[label.name] ? `tree-label-${roles[label.name]}` : ""} ${markedNames.includes(label.name) ? "tree-label-marked" : ""}`} x={label.x} y={label.y + 3} onClick={() => onToggleMark?.(label.name)} onKeyDown={(keyEvent) => { if (onToggleMark && (keyEvent.key === "Enter" || keyEvent.key === " ")) { keyEvent.preventDefault(); onToggleMark(label.name); } }}>{label.name}</text></g>)}
    </svg>
    {layout.zeroLengthBranches > 0 && <div className="tree-zero-note">{layout.zeroLengthBranches} zero-length branch{layout.zeroLengthBranches === 1 ? "" : "es"} collapse at a drawn node; all descendants remain connected.</div>}
  </div>;
}

function LocalTrees({ alignment, events, event, onSelectEvent }: { alignment: AlignmentData; events: RdpEvent[]; event: RdpEvent | null; onSelectEvent: (id: string) => void }) {
  const [cohortSize, setCohortSize] = useState(14);
  const [cohortMode, setCohortMode] = useState<"nearest" | "first">("nearest");
  const [markedNames, setMarkedNames] = useState<string[]>([]);
  const trees = useMemo(() => {
    if (!event) return null;
    const triad = [event.recombinant, event.majorParent, event.minorParent];
    const remaining = alignment.sequences.map((record, index) => ({
      index,
      score: cohortMode === "nearest"
        ? triad.reduce((sum, triadIndex) => sum + pairwiseIdentitySampled(record.sequence, alignment.sequences[triadIndex].sequence), 0) / triad.length
        : -index,
    })).filter(({ index }) => !triad.includes(index)).sort((left, right) => right.score - left.score);
    const selected = [...triad, ...remaining.slice(0, Math.max(0, cohortSize - triad.length)).map(({ index }) => index)];
    const regions = event.wraps
      ? [
          { label: "Origin-spanning tract", slug: "tract", segments: eventSegments(event, alignment.length), highlight: true },
          { label: "Complementary background", slug: "background", segments: [[event.end, event.start] as [number, number]], highlight: false },
        ]
      : [
          ...(event.start > 1 ? [{ label: "Left flank", slug: "left-flank", segments: [[0, event.start] as [number, number]], highlight: false }] : []),
          { label: "Recombinant tract", slug: "tract", segments: [[event.start, event.end] as [number, number]], highlight: true },
          ...(event.end < alignment.length - 1 ? [{ label: "Right flank", slug: "right-flank", segments: [[event.end, alignment.length] as [number, number]], highlight: false }] : []),
        ];
    return {
      selected,
      regions: regions.map((region) => ({ ...region, tree: buildLocalTree(alignment, selected, region.segments) })),
    };
  }, [alignment, cohortMode, cohortSize, event]);
  if (!event || !trees) return <div className="empty-state tree-empty"><span className="empty-mark">⑂</span><h2>Select an event to compare local trees</h2><p>Tree views are hypothesis-centered: choose an event in Events, or use the selector above after a scan.</p></div>;
  const recombinant = alignment.sequences[event.recombinant].name;
  const major = alignment.sequences[event.majorParent].name;
  const minor = alignment.sequences[event.minorParent].name;
  const roles: Record<string, "recombinant" | "major" | "minor"> = { [recombinant]: "recombinant", [major]: "major", [minor]: "minor" };
  const allNewick = trees.regions.map((region) => `[${region.label}]\n${region.tree.newick}`).join("\n\n");
  const toggleMark = (name: string) => setMarkedNames((current) => current.includes(name) ? current.filter((candidate) => candidate !== name) : [...current, name]);
  return <div className="tree-explorer">
    <div className="tree-hero">
      <div><span className="eyebrow">Selected recombination hypothesis</span><h1>{recombinant}</h1><p>Does the recombinant move from the major-parent neighborhood in its flanks to the minor-parent neighborhood inside the proposed tract?</p></div>
      <div className="tree-hypothesis">
        <span><i className="parent-role major"/>Major parent<b>{major}</b></span>
        <span><i className="parent-role minor"/>Minor parent<b>{minor}</b></span>
        <span><i className="parent-role recombinant"/>Region<b>{formatEventRegion(event, alignment.length)}</b></span>
      </div>
    </div>
    <div className="tree-toolbar">
      <label><span>Event hypothesis</span><select value={event.id} onChange={(value) => onSelectEvent(value.target.value)}>{events.map((candidate, index) => <option value={candidate.id} key={candidate.id}>E{index + 1} · {alignment.sequences[candidate.recombinant].name} · {formatEventRegion(candidate, alignment.length)}</option>)}</select></label>
      <label><span>Context sequences</span><select value={cohortSize} onChange={(value) => setCohortSize(Number(value.target.value))}>{[8, 14, 20, 24].filter((count) => count <= Math.max(8, alignment.sequences.length)).map((count) => <option value={count} key={count}>Up to {count}</option>)}</select></label>
      <label><span>Cohort strategy</span><select value={cohortMode} onChange={(value) => setCohortMode(value.target.value as "nearest" | "first")}><option value="nearest">Nearest to event triad</option><option value="first">First sequences in alignment</option></select></label>
      <div className="tree-toolbar-actions"><button type="button" className="small-button" onClick={() => download(`rdp-event-${event.id}-local-trees.nwk`, allNewick)}>All Newick ↓</button><button type="button" className="small-button" disabled={!markedNames.length} onClick={() => setMarkedNames([])}>Clear marks</button></div>
    </div>
    <div className="tree-mark-help"><b>Linked marking</b><span>Click any leaf name—or a cohort chip below—to mark that sequence in every regional tree. This makes topology switches much easier to track.</span><div>{markedNames.map((name) => <button type="button" key={name} onClick={() => toggleMark(name)}>{name} ×</button>)}</div></div>
    <TopologyCards alignment={alignment} event={event}/>
    <div className={`tree-comparison regions-${trees.regions.length}`}>{trees.regions.map((item) => <article key={item.label} className={item.highlight ? "tract-tree" : ""}><div><div><b>{item.label}</b><span>{item.segments.map(([start, end]) => `${start + 1}–${end}`).join(" + ")}</span></div><button type="button" onClick={() => download(`rdp-${item.slug}.nwk`, item.tree.newick)}>Newick ↓</button></div><TreeSvg root={item.tree.root} roles={roles} markedNames={markedNames} onToggleMark={toggleMark}/><code title={item.tree.newick}>{item.tree.newick}</code></article>)}</div>
    <div className="tree-cohort"><strong>Tree cohort · {trees.selected.length} sequences</strong><div>{trees.selected.map((index) => { const name = alignment.sequences[index].name; return <button type="button" aria-pressed={markedNames.includes(name)} onClick={() => toggleMark(name)} className={`${index === event.recombinant ? "recombinant" : index === event.majorParent ? "major" : index === event.minorParent ? "minor" : ""} ${markedNames.includes(name) ? "marked" : ""}`} key={index}>{name}</button>; })}</div></div>
    <div className="tree-caveat"><strong>Exploratory local phylogenies</strong><span>Neighbor-joining on uncorrected canonical-site p-distances, with at most 8,192 sampled sites per comparison. These trees are interactive verification aids—not ML trees, topology tests, or independent proof of recombination.</span></div>
  </div>;
}

function AutoRange({ label, value, minimum, maximum, step = 1, display, onChange }: {
  label: string;
  value: number;
  minimum: number;
  maximum: number;
  step?: number;
  display?: (value: number) => string;
  onChange: (value: number) => void;
}) {
  return <label className="auto-range"><span>{label}<output>{display ? display(value) : value}</output></span><input type="range" min={minimum} max={maximum} step={step} value={value} onChange={(event) => onChange(Number(event.target.value))}/></label>;
}

function AutoSwitch({ label, detail, checked, onChange }: { label: string; detail: string; checked: boolean; onChange: (checked: boolean) => void }) {
  return <label className="auto-switch"><input type="checkbox" checked={checked} onChange={(event) => onChange(event.target.checked)}/><i/><span><b>{label}</b><small>{detail}</small></span></label>;
}

function AutoResolvePanel({ alignment, events, status, rescanning, onRun }: {
  alignment: AlignmentData;
  events: RdpEvent[];
  status: AutoResolveStatus;
  rescanning: boolean;
  onRun: (settings: AutoResolveSettings, profileLabel: string) => void;
}) {
  const [profile, setProfile] = useState<AutoResolvePresetName | "custom">("conservative");
  const [settings, setSettings] = useState<AutoResolveSettings>(() => ({ ...AUTO_RESOLVE_PRESETS.conservative }));
  const deferredSettings = useDeferredValue(settings);
  const plan = useMemo(() => planAutoResolution(events, alignment.length, deferredSettings), [alignment.length, deferredSettings, events]);
  const visibleEntries = plan.entries.filter((entry) => entry.recommendation !== "keep").slice(0, 8);
  const firstBarrier = plan.barriers[0];
  const running = status.state === "resolving" || status.state === "rescanning" || rescanning;
  const setValue = <K extends keyof AutoResolveSettings>(key: K, value: AutoResolveSettings[K]) => {
    setProfile("custom");
    setSettings((current) => ({ ...current, [key]: value }));
  };
  const selectPreset = (next: AutoResolvePresetName) => {
    setProfile(next);
    setSettings({ ...AUTO_RESOLVE_PRESETS[next] });
  };
  return <Panel title="Heuristic auto-resolver" action={<span className="panel-caption">Dry-run first · reversible as one workflow action</span>}>
    <div className="auto-resolve-shell">
      <div className="auto-resolve-intro">
        <div><span className="eyebrow">Explainable ordered decisions</span><h2>Resolve until the dependency graph says “rescan”.</h2><p>Each event is scored in inference order. Automatic acceptance requires every configured scientific gate; ambiguous or stale evidence stays with the analyst. A rescan is inserted only when a changed decision can alter a later recombinant, parent proxy, overlapping tract, or co-recombinant group.</p></div>
        <div className="auto-preset-picker" role="radiogroup" aria-label="Auto-resolve profile">
          {(["conservative", "balanced", "aggressive"] as AutoResolvePresetName[]).map((name) => <button type="button" role="radio" aria-checked={profile === name} className={profile === name ? "active" : ""} key={name} onClick={() => selectPreset(name)}><b>{name}</b><small>{name === "conservative" ? "High-confidence automation; earlier targeted rescans" : name === "balanced" ? "Moderate gates; adaptive scan scope" : "Resolve more; batch more; fewer rescans"}</small></button>)}
          {profile === "custom" && <span className="custom-profile">Custom advanced profile</span>}
        </div>
      </div>
      <div className="auto-plan-summary" aria-label="Auto-resolution dry-run summary">
        <article className="accept"><b>{plan.acceptCount}</b><span>auto-accept</span></article>
        <article className="reject"><b>{plan.rejectCount}</b><span>auto-reject</span></article>
        <article className="review"><b>{plan.reviewCount}</b><span>analyst review</span></article>
        <article><b>{plan.keepCount}</b><span>locked decisions</span></article>
        <article className={firstBarrier ? "rescan" : "quiet"}><b>{plan.barriers.length}</b><span>next rescan barrier</span></article>
      </div>
      <div className="auto-run-row">
        <div className={firstBarrier ? "auto-rescan-forecast active" : "auto-rescan-forecast"}>
          <span>{firstBarrier ? `First barrier after E${firstBarrier.afterEventIndex + 1}` : settings.rescanStrategy === "off" ? "Automatic rescans disabled" : "No dependency-triggered rescan predicted"}</span>
          <b>{firstBarrier ? `${firstBarrier.risk}/100 risk · ${firstBarrier.impactedTargetIndexes.length} target${firstBarrier.impactedTargetIndexes.length === 1 ? "" : "s"}` : "Queue can run continuously"}</b>
          {firstBarrier && <small>{firstBarrier.reasons.slice(0, 2).join("; ")}</small>}
        </div>
        <button type="button" className="run-button auto-run-button" disabled={events.length === 0 || running} onClick={() => onRun(settings, profile === "custom" ? "custom" : profile)}>{running ? status.state === "rescanning" ? "Rescanning impacted targets…" : "Resolving queue…" : "Run auto-resolve"}<span>→</span></button>
      </div>
      {status.state !== "idle" && <div className={`auto-status ${status.state}`} role="status"><span>{status.state === "complete" ? "✓" : status.state === "error" ? "!" : status.state === "paused" ? "Ⅱ" : "↻"}</span><div><b>{status.message}</b><small>Round {status.round + 1} · {status.processed} processed · {status.accepted} accepted · {status.rejected} rejected · {status.reviewed} held · {status.rescans} rescan{status.rescans === 1 ? "" : "s"}</small></div></div>}
      <div className="auto-plan-list" aria-label="First auto-resolution recommendations">
        {visibleEntries.length === 0 ? <div className="empty-state compact">Nothing is eligible for automatic resolution under this profile.</div> : visibleEntries.map((entry) => {
          const event = events[entry.eventIndex];
          return <article className={entry.recommendation} key={entry.eventId}><span>E{entry.eventIndex + 1}</span><div><b>{alignment.sequences[event.recombinant]?.name ?? "Unknown recombinant"}</b><small>{entry.reasons[0]}</small>{entry.impactRisk > 0 && <em>downstream risk {entry.impactRisk}/100 · {entry.impactedTargetIndexes.length} target{entry.impactedTargetIndexes.length === 1 ? "" : "s"}</em>}</div><strong>{entry.score}</strong><i>{entry.recommendation}</i></article>;
        })}
        {plan.entries.filter((entry) => entry.recommendation !== "keep").length > visibleEntries.length && <div className="auto-plan-cap">Showing the first {visibleEntries.length} actionable recommendations; the run replans after every rescan.</div>}
      </div>
      <details className="auto-advanced"><summary><span>Advanced decision model</span><small>Precisely tune evidence gates, score weights and rescan sensitivity</small></summary>
        <div className="auto-advanced-grid">
          <section><h3>Hard decision gates</h3><p>Every gate must pass before auto-acceptance, irrespective of the composite score.</p>
            <AutoRange label="Accept score" value={settings.acceptScore} minimum={50} maximum={99} onChange={(value) => setValue("acceptScore", value)}/>
            <AutoRange label="Reject score" value={settings.rejectScore} minimum={0} maximum={60} onChange={(value) => setValue("rejectScore", value)}/>
            <AutoRange label="Supporting methods" value={settings.minimumSupportingMethods} minimum={1} maximum={7} onChange={(value) => setValue("minimumSupportingMethods", value)}/>
            <label className="auto-select"><span>Maximum adjusted P</span><select value={settings.maximumCorrectedP} onChange={(event) => setValue("maximumCorrectedP", Number(event.target.value))}><option value={0.1}>0.10</option><option value={0.05}>0.05</option><option value={0.01}>0.01</option><option value={0.001}>0.001</option><option value={0.0001}>0.0001</option></select></label>
            <AutoRange label="Informative sites" value={settings.minimumInformativeSites} minimum={0} maximum={200} step={4} onChange={(value) => setValue("minimumInformativeSites", value)}/>
            <AutoRange label="Breakpoint uncertainty" value={settings.maximumBreakpointUncertainty} minimum={0.05} maximum={1} step={0.05} display={(value) => `${Math.round(value * 100)}% tract`} onChange={(value) => setValue("maximumBreakpointUncertainty", value)}/>
            <AutoRange label="Parent-conflict ceiling" value={settings.maximumParentConflict} minimum={0.02} maximum={0.8} step={0.02} display={(value) => `${Math.round(value * 100)}%`} onChange={(value) => setValue("maximumParentConflict", value)}/>
            <AutoRange label="Rate-density fold ceiling" value={settings.maximumRateFold} minimum={1.2} maximum={10} step={0.2} display={(value) => `${value.toFixed(1)}×`} onChange={(value) => setValue("maximumRateFold", value)}/>
            <AutoSwitch label="Block high-risk warnings" detail="Misalignment, homoplasy, diffuse incompatibility, rate-shift and incomplete-scan flags prevent acceptance." checked={settings.blockSevereWarnings} onChange={(value) => setValue("blockSevereWarnings", value)}/>
            <AutoSwitch label="Revisit analyst decisions" detail="Off by default: accepted and rejected events remain locked unless their evidence is stale." checked={settings.revisitReviewed} onChange={(value) => setValue("revisitReviewed", value)}/>
          </section>
          <section><h3>Composite score weights</h3><p>Weights are relative; setting a weight to zero removes that evidence family from the score.</p>
            <AutoRange label="Method consensus" value={settings.consensusWeight} minimum={0} maximum={50} onChange={(value) => setValue("consensusWeight", value)}/>
            <AutoRange label="Adjusted significance" value={settings.significanceWeight} minimum={0} maximum={50} onChange={(value) => setValue("significanceWeight", value)}/>
            <AutoRange label="Informative-site depth" value={settings.informationWeight} minimum={0} maximum={50} onChange={(value) => setValue("informationWeight", value)}/>
            <AutoRange label="Breakpoint precision" value={settings.breakpointWeight} minimum={0} maximum={50} onChange={(value) => setValue("breakpointWeight", value)}/>
            <AutoRange label="False-positive diagnostics" value={settings.diagnosticsWeight} minimum={0} maximum={50} onChange={(value) => setValue("diagnosticsWeight", value)}/>
            <div className="auto-score-note"><b>Safety invariant</b><span>Stale or uncalibrated evidence is never auto-accepted. The review band is intentionally preserved between the accept and reject thresholds.</span></div>
          </section>
          <section><h3>Dependency-aware rescanning</h3><p>The planner accumulates risk and pauses at the first barrier, withholds accepted mosaic sequences from parent search, rescans, then rebuilds the remaining plan from new signals.</p>
            <label className="auto-select"><span>Rescan strategy</span><select value={settings.rescanStrategy} onChange={(event) => setValue("rescanStrategy", event.target.value as AutoResolveSettings["rescanStrategy"])}><option value="off">Off · decisions only</option><option value="targeted">Impacted recombinant targets only</option><option value="adaptive">Targeted unless impact is broad</option><option value="full">All unresolved targets</option></select></label>
            <AutoRange label="Trigger risk" value={settings.rescanRiskThreshold} minimum={0} maximum={100} onChange={(value) => setValue("rescanRiskThreshold", value)}/>
            <AutoRange label="Same-sequence overlap trigger" value={settings.overlapTriggerFraction} minimum={0} maximum={1} step={0.05} display={(value) => `${Math.round(value * 100)}%`} onChange={(value) => setValue("overlapTriggerFraction", value)}/>
            <AutoRange label="Minimum queue gap" value={settings.minimumEventsBetweenRescans} minimum={1} maximum={25} display={(value) => `${value} event${value === 1 ? "" : "s"}`} onChange={(value) => setValue("minimumEventsBetweenRescans", value)}/>
            <AutoRange label="Maximum rescan rounds" value={settings.maximumRescanRounds} minimum={0} maximum={8} onChange={(value) => setValue("maximumRescanRounds", value)}/>
            <AutoRange label="Adaptive full-scan threshold" value={settings.adaptiveFullTargetFraction} minimum={0.05} maximum={1} step={0.05} display={(value) => `${Math.round(value * 100)}% sequences`} onChange={(value) => setValue("adaptiveFullTargetFraction", value)}/>
            <AutoRange label="Same recombinant risk" value={settings.sameRecombinantRisk} minimum={0} maximum={100} onChange={(value) => setValue("sameRecombinantRisk", value)}/>
            <AutoRange label="Recombinant-parent risk" value={settings.recombinantParentRisk} minimum={0} maximum={100} onChange={(value) => setValue("recombinantParentRisk", value)}/>
            <AutoRange label="Shared-group risk" value={settings.groupedEventRisk} minimum={0} maximum={100} onChange={(value) => setValue("groupedEventRisk", value)}/>
            <AutoRange label="Withdraw accepted-event risk" value={settings.acceptedWithdrawalRisk} minimum={0} maximum={100} onChange={(value) => setValue("acceptedWithdrawalRisk", value)}/>
          </section>
        </div>
      </details>
    </div>
  </Panel>;
}

function ReconstructionWorkspace({ alignment, events, selectedId, onSelect, onDecision, onRescan, rescanning, autoResolveStatus, onAutoResolve }: {
  alignment: AlignmentData;
  events: RdpEvent[];
  selectedId: string | null;
  onSelect: (id: string, view?: "explore" | "trees" | "alignment") => void;
  onDecision: (id: string, decision: EventDecision) => void;
  onRescan: () => void;
  rescanning: boolean;
  autoResolveStatus: AutoResolveStatus;
  onAutoResolve: (settings: AutoResolveSettings, profileLabel: string) => void;
}) {
  const [queueLimit, setQueueLimit] = useState(250);
  const model = useMemo(() => buildReconstructionModel(events, alignment.length), [alignment.length, events]);
  const relationsByEvent = useMemo(() => {
    const map = new Map<number, ReconstructionRelationship[]>();
    model.relationships.forEach((relationship) => {
      map.set(relationship.fromIndex, [...(map.get(relationship.fromIndex) ?? []), relationship]);
      map.set(relationship.toIndex, [...(map.get(relationship.toIndex) ?? []), relationship]);
    });
    return map;
  }, [model.relationships]);
  const accepted = events.filter((event) => event.decision === "accepted" && !event.evidenceStale).length;
  const reviewed = events.filter((event) => event.decision !== "unreviewed").length;
  const grouped = events.filter((event) => event.groupId).length;
  const canRescan = accepted > 0 && !rescanning;
  const relationLabel = (kind: "possible-overprint" | "recombinant-parent" | "event-group", otherIndex: number, overlapBases: number) => {
    if (kind === "possible-overprint") return `possible overprint with E${otherIndex + 1} · ${formatInteger(overlapBases)} nt`;
    if (kind === "recombinant-parent") return `recombinant-parent dependency · E${otherIndex + 1}`;
    return `same co-recombinant group as E${otherIndex + 1}`;
  };
  return <div className="reconstruction-workspace">
    <div className="reconstruction-hero">
      <div><span className="eyebrow">Global event reconstruction</span><h1>Build the collection-level mosaic in review order.</h1><p>Desktop RDP separates raw signals from unique events, reviews the strongest events first, groups co-recombinant descendants, then rescans after edits. This workspace makes that dependency chain explicit.</p></div>
      <div className="reconstruction-hero-actions">
        <button type="button" className="run-button" disabled={model.nextReviewIndex === null} onClick={() => model.nextReviewIndex !== null && onSelect(events[model.nextReviewIndex].id, "explore")}>{model.nextReviewIndex === null ? "Review complete" : `Review next · E${model.nextReviewIndex + 1}`} <span>→</span></button>
        <button type="button" className="open-button" disabled={!canRescan} onClick={onRescan}>↻ Rescan unresolved</button>
      </div>
    </div>
    <div className="reconstruction-stages" aria-label="Event reconstruction workflow">
      {[
        ["1", "Scan", events.length ? `${events.length} signals retained` : "Run the multi-method scan", events.length > 0],
        ["2", "Reconcile", grouped ? `${grouped} grouped hypotheses` : "Group shared ancestry", grouped > 0 || events.length <= 1],
        ["3", "Review in order", `${reviewed}/${events.length} decided`, events.length > 0 && reviewed === events.length],
        ["4", "Rescan after edits", model.staleFromIndex === null ? "No stale edit chain" : `Required from E${model.staleFromIndex + 1}`, model.staleFromIndex === null],
        ["5", "Global mosaic", `${accepted} accepted + fresh`, accepted > 0 && model.nextReviewIndex === null],
      ].map(([number, label, note, complete]) => <article className={complete ? "complete" : ""} key={String(label)}><span>{complete ? "✓" : number}</span><div><b>{label}</b><small>{note}</small></div></article>)}
    </div>
    <AutoResolvePanel alignment={alignment} events={events} status={autoResolveStatus} rescanning={rescanning} onRun={onAutoResolve}/>
    {model.staleFromIndex !== null && <div className="reconstruction-warning"><strong>Downstream characterization may now be stale.</strong><span>E{model.staleFromIndex + 1} was edited after its evidence was calculated. RDP’s ordered workflow requires the remaining signals to be rescanned; {model.downstreamIndexes.length} later retained event{model.downstreamIndexes.length === 1 ? "" : "s"} are marked downstream, not silently treated as independent.{!canRescan && " Recalculate and accept the edited event before rescanning."}</span><button type="button" disabled={!canRescan} title={canRescan ? "Keep accepted fresh events and rescan unresolved targets" : "Recalculate and accept at least one event first"} onClick={onRescan}>Rescan unresolved sequences</button></div>}
    <Panel title="Ordered event queue" action={<span className="panel-caption">Inference/review order · not literal historical time</span>}>
      <div className="event-queue">
        {events.length === 0 ? <div className="empty-state compact">Run a scan or create a manual hypothesis to begin event reconstruction.</div> : events.slice(0, queueLimit).map((event, eventIndex) => {
          const relationships = relationsByEvent.get(eventIndex) ?? [];
          const downstream = model.downstreamIndexes.includes(eventIndex);
          return <article key={event.id} className={`${event.id === selectedId ? "selected" : ""} ${event.decision} ${event.evidenceStale ? "stale" : ""}`}>
            <button type="button" className="queue-number" onClick={() => onSelect(event.id, "explore")}><span>E{eventIndex + 1}</span><small>{event.source === "wasm" ? "scan" : event.source === "manual" ? "manual" : "truth"}</small></button>
            <div className="queue-summary"><b>{alignment.sequences[event.recombinant]?.name ?? "Unknown recombinant"}</b><span>{formatEventRegion(event, alignment.length)} · {alignment.sequences[event.majorParent]?.name ?? "?"} ↔ {alignment.sequences[event.minorParent]?.name ?? "?"}</span><div>{event.groupId && <em>{event.groupId}</em>}{event.evidenceStale && <em className="stale">recalculate</em>}{downstream && <em className="downstream">downstream of edit</em>}{relationships.slice(0, 3).map((relationship, relationIndex) => { const otherIndex = relationship.fromIndex === eventIndex ? relationship.toIndex : relationship.fromIndex; return <em className={relationship.kind} key={`${relationship.kind}-${otherIndex}-${relationIndex}`}>{relationLabel(relationship.kind, otherIndex, relationship.overlapBases)}</em>; })}</div></div>
            <div className="queue-decision"><span className={`decision-label ${event.decision}`}>{event.decision}</span><button type="button" onClick={() => onSelect(event.id, "alignment")}>Alignment</button><button type="button" onClick={() => onSelect(event.id, "trees")}>Trees</button></div>
            <div className="queue-review"><button type="button" disabled={rescanning} className={event.decision === "rejected" ? "active reject" : "reject"} onClick={() => onDecision(event.id, event.decision === "rejected" ? "unreviewed" : "rejected")}>Reject</button><button type="button" disabled={rescanning} className={event.decision === "accepted" ? "active accept" : "accept"} onClick={() => onDecision(event.id, event.decision === "accepted" ? "unreviewed" : "accepted")}>Accept</button></div>
          </article>;
        })}{events.length > queueLimit && <button type="button" className="queue-load-more" onClick={() => setQueueLimit((current) => Math.min(events.length, current + 250))}>Show {formatInteger(Math.min(250, events.length - queueLimit))} more events · {formatInteger(events.length - queueLimit)} hidden</button>}
      </div>
    </Panel>
    <Panel title="Global mosaic map" action={<span className="panel-caption">All retained events · nested tracts use separate lanes</span>}>
      <div className="global-mosaic">
        <div className="mosaic-axis"><span>1</span><span>{formatInteger(Math.round(alignment.length / 2))}</span><span>{formatInteger(alignment.length)} nt</span></div>
        {model.sequenceRows.length === 0 ? <div className="empty-state compact">Retained events will be assembled into sequence-level mosaics here.</div> : model.sequenceRows.slice(0, 300).map((row) => <div className="mosaic-row" key={row.sequenceIndex}><button type="button" title={alignment.sequences[row.sequenceIndex]?.name}>{alignment.sequences[row.sequenceIndex]?.name ?? `Sequence ${row.sequenceIndex + 1}`}</button><div className="mosaic-track">{row.eventIndexes.flatMap((eventIndex, laneIndex) => {
          const event = events[eventIndex];
          return eventSegments(event, alignment.length).map(([start, end], segmentIndex) => <button type="button" key={`${event.id}-${segmentIndex}`} className={`${event.decision} ${event.evidenceStale ? "stale" : ""} ${event.id === selectedId ? "selected" : ""}`} style={{ left: `${start / alignment.length * 100}%`, width: `${Math.max(.45, (end - start) / alignment.length * 100)}%`, top: `${4 + (laneIndex % 4) * 9}px` }} title={`E${eventIndex + 1} · ${formatEventRegion(event, alignment.length)} · ${event.decision}`} onClick={() => onSelect(event.id, "explore")}><span>E{eventIndex + 1}</span></button>);
        })}</div></div>)}{model.sequenceRows.length > 300 && <div className="mosaic-cap">Showing the first 300 recombinant sequences; project export retains all {formatInteger(model.sequenceRows.length)}.</div>}
      </div>
    </Panel>
    <div className="reconstruction-bottom-grid">
      <Panel title="Nested events & parent dependencies" action={<span className="panel-caption">Hypotheses to verify in trees</span>}>
        <div className="dependency-list">{model.relationships.length === 0 ? <div className="empty-state compact">No overlapping, grouped, or recombinant-parent dependencies are currently inferred.</div> : model.relationships.slice(0, 500).map((relationship, index) => {
          const from = events[relationship.fromIndex];
          const to = events[relationship.toIndex];
          return <button type="button" key={`${relationship.kind}-${relationship.fromIndex}-${relationship.toIndex}-${index}`} onClick={() => onSelect(to.id, "trees")}><span className={relationship.kind}>{relationship.kind === "possible-overprint" ? "Overlap" : relationship.kind === "recombinant-parent" ? "Nested parent" : "Event group"}</span><b>E{relationship.fromIndex + 1} → E{relationship.toIndex + 1}</b><small>{relationship.kind === "possible-overprint" ? `${alignment.sequences[to.recombinant]?.name} shares ${formatInteger(relationship.overlapBases)} tract bases; test possible overprinting.` : relationship.kind === "recombinant-parent" ? `${alignment.sequences[from.recombinant]?.name} is used as a parent proxy in another event.` : `${from.groupId} links co-recombinant descendants.`}</small><i>Compare trees →</i></button>;
        })}{(model.relationships.length > 500 || model.relationshipsTruncated) && <div className="dependency-cap">Showing 500 of {model.relationshipsTruncated ? "at least " : ""}{formatInteger(model.relationships.length)} derived relationships to keep the interface responsive{model.relationshipsTruncated ? "; derivation is capped at 20,000 relationships" : ""}.</div>}</div>
      </Panel>
      <Panel title="How to reconstruct defensibly" action={<a className="panel-help-link" href="https://web.cbio.uct.ac.za/~darren/RDP5Manual.pdf" target="_blank" rel="noreferrer">RDP5 manual ↗</a>}>
        <ol className="reconstruction-guide"><li><b>Start with the strongest/earliest characterized signal.</b><span>Later assignments can depend on how earlier recombinant descendants were grouped.</span></li><li><b>Verify the alignment and all three local trees.</b><span>Track the recombinant, major proxy, and minor proxy across tract and background trees.</span></li><li><b>Group co-recombinant descendants and annotate nesting.</b><span>A sampled “parent” is a proxy relative, not necessarily the historical donor.</span></li><li><b>After breakpoint, parent, recombinant, or grouping edits, rescan.</b><span>Do not trust the untouched tail of the queue as if it were independent of the edit.</span></li><li><b>Accept a global mosaic only after the queue is resolved.</b><span>Exported provenance keeps rejected alternatives, edits, and the project ledger.</span></li></ol>
      </Panel>
    </div>
  </div>;
}

function ExampleLibrary({ onClose, onLoad }: { onClose: () => void; onLoad: (example: ExampleDataset) => void }) {
  const [filter, setFilter] = useState<"All" | ExampleDataset["complexity"]>("All");
  const visible = EXAMPLE_DATASETS.filter((example) => filter === "All" || example.complexity === filter);
  return <div className="modal-backdrop" role="presentation" onMouseDown={onClose}><section className="modal example-modal" role="dialog" aria-modal="true" aria-labelledby="examples-title" onMouseDown={(event) => event.stopPropagation()}>
    <button className="modal-close" type="button" onClick={onClose} aria-label="Close"><Icon name="close"/></button>
    <div className="example-heading"><div><span className="eyebrow">Synthetic, deterministic, truth annotated</span><h2 id="examples-title">Example dataset library</h2><p>Learn the workflow on a clean triplet, explore circular and nested events, or stress the sampled parent screen with bacterial and 500-genome panels. Examples contain no empirical sequences.</p></div><div className="example-filter">{(["All", "Starter", "Intermediate", "Advanced", "Stress test"] as const).map((item) => <button type="button" className={filter === item ? "active" : ""} onClick={() => setFilter(item)} key={item}>{item}</button>)}</div></div>
    <div className="example-grid">{visible.map((example) => <article key={example.id} className={`complexity-${example.complexity.toLowerCase().replace(" ", "-")}`}>
      <div className="example-card-head"><span>{example.complexity}</span><small>{example.organism}</small></div>
      <h3>{example.title}</h3>
      <div className="example-size"><b>{example.sequenceCount.toLocaleString("en-US")}</b><span>sequences</span><b>{example.length.toLocaleString("en-US")}</b><span>sites</span></div>
      <p>{example.description}</p>
      <div className="example-tags">{example.tags.map((tag) => <span key={tag}>{tag}</span>)}</div>
      <details><summary>Known truth · {example.truth.length} item{example.truth.length === 1 ? "" : "s"}</summary><div>{example.truth.map((truth, index) => <p key={`${truth.recombinant}-${index}`}><b>{truth.recombinant}</b><span>{truth.donor} · {truth.region}</span><small>{truth.note}</small></p>)}</div></details>
      <div className="example-challenge"><strong>Why use it</strong><span>{example.challenge}</span></div>
      <button type="button" className="run-inline" onClick={() => onLoad(example)}>Generate and load locally <span>→</span></button>
    </article>)}</div>
  </section></div>;
}

function Inspector({ alignment, event, onDecision, onUpdate, onNavigate, onDuplicate, onDelete, onRecalculate, onOpenTrees, onOpenMethod, recalculating, canPrevious, canNext, groupIds }: {
  alignment: AlignmentData;
  event: RdpEvent | null;
  onDecision: (decision: EventDecision) => void;
  onUpdate: (patch: Partial<RdpEvent>, action?: string) => void;
  onNavigate: (direction: -1 | 1) => void;
  onDuplicate: () => void;
  onDelete: () => void;
  onRecalculate: () => void;
  onOpenTrees: () => void;
  onOpenMethod: (method: MethodName) => void;
  recalculating: boolean;
  canPrevious: boolean;
  canNext: boolean;
  groupIds: string[];
}) {
  if (!event) return <aside className="inspector"><div className="empty-inspector"><span>↗</span><h2>No event selected</h2><p>Run a scan or choose an event in the overview.</p></div></aside>;
  const recombinant = alignment.sequences[event.recombinant]?.name;
  let newGroupId = `group-${groupIds.length + 1}`;
  let groupSuffix = groupIds.length + 1;
  while (groupIds.includes(newGroupId)) {
    groupSuffix += 1;
    newGroupId = `group-${groupSuffix}`;
  }
  return (
    <aside className="inspector">
      <div className="inspector-header">
        <div><span className="eyebrow">Selected hypothesis</span><h2>{recombinant}</h2></div>
        <div className="event-nav"><button type="button" title="Duplicate event" onClick={onDuplicate}>⧉</button><button type="button" title="Delete event (undoable)" onClick={onDelete}>⌫</button><button type="button" disabled={!canPrevious} onClick={() => onNavigate(-1)}>←</button><button type="button" disabled={!canNext} onClick={() => onNavigate(1)}>→</button></div>
      </div>
      <button type="button" className="inspect-trees-button" onClick={onOpenTrees}><span>⑂</span><div><b>Compare local trees</b><small>Flanks vs recombinant tract</small></div><i>→</i></button>
      <div className="parent-map assignment-map">
        <label><span>Recombinant</span><select value={event.recombinant} onChange={(value) => onUpdate({ recombinant: Number(value.target.value) }, "Reassigned recombinant")}>
          {alignment.sequences.map((record, index) => <option key={`${record.name}-${index}`} value={index} disabled={index === event.majorParent || index === event.minorParent}>{record.name}</option>)}
        </select></label>
        <label><span>Major parent</span><select value={event.majorParent} onChange={(value) => onUpdate({ majorParent: Number(value.target.value) }, "Reassigned major parent")}>
          {alignment.sequences.map((record, index) => <option key={`${record.name}-${index}`} value={index} disabled={index === event.recombinant || index === event.minorParent}>{record.name}</option>)}
        </select></label>
        <label><span>Minor parent</span><select value={event.minorParent} onChange={(value) => onUpdate({ minorParent: Number(value.target.value) }, "Reassigned minor parent")}>
          {alignment.sequences.map((record, index) => <option key={`${record.name}-${index}`} value={index} disabled={index === event.recombinant || index === event.majorParent}>{record.name}</option>)}
        </select></label>
      </div>
      {(event.alternativeParents?.length ?? 0) > 0 && <div className="alternative-parents"><span>Alternative minor parents retained by deduplication</span><div>{event.alternativeParents?.map((index) => <button type="button" key={index} onClick={() => onUpdate({ minorParent: index }, `Selected alternative minor parent ${alignment.sequences[index]?.name}`)}>{alignment.sequences[index]?.name ?? `Sequence ${index + 1}`}</button>)}</div></div>}
      {event.evidenceStale && <div className="stale-box"><strong>Evidence needs recalculation</strong><span>Assignments or breakpoints changed after the scan. The saved method statistics remain visible for comparison but no longer describe this edited hypothesis.</span><button type="button" onClick={onRecalculate} disabled={recalculating}>{recalculating ? "Recalculating…" : "Recalculate this exact hypothesis"}</button></div>}
      <div className="inspector-section">
        <div className="section-heading"><h3>Breakpoint interval</h3><span className="mono">{eventLength(event, alignment.length).toLocaleString("en-US")} nt{event.wraps ? " · ↻" : ""}</span></div>
        <div className="breakpoint-inputs">
          <label>Start<input type="number" min={event.wraps ? event.end + 2 : 1} max={event.wraps ? alignment.length : event.end} value={event.start + 1} onChange={(value) => {
            const requested = Number(value.target.value) - 1;
            const start = event.wraps
              ? Math.max(event.end + 1, Math.min(alignment.length - 1, requested))
              : Math.max(0, Math.min(event.end - 1, requested));
            onUpdate({ start, confidenceStart: [start, start], breakpointModel: { method: "manual", informativeSites: event.informativeSites } }, "Edited left breakpoint");
          }}/><small>{event.confidenceStart[0] + 1}–{event.confidenceStart[1]}</small></label>
          <label>End<input type="number" min={1} max={event.wraps ? Math.max(1, event.start - 1) : alignment.length} value={event.end} onChange={(value) => {
            const requested = Number(value.target.value);
            const end = event.wraps
              ? Math.max(1, Math.min(Math.max(1, event.start - 1), requested))
              : Math.max(event.start + 1, Math.min(alignment.length, requested));
            onUpdate({ end, confidenceEnd: [end, end], breakpointModel: { method: "manual", informativeSites: event.informativeSites } }, "Edited right breakpoint");
          }}/><small>{event.confidenceEnd[0] + 1}–{event.confidenceEnd[1]}</small></label>
        </div>
        {event.breakpointModel && <div className="breakpoint-model"><b>{event.breakpointModel.method === "two-state-hmm" ? "Windowless two-state HMM" : event.breakpointModel.method === "manual" ? "Manually edited" : "Local χ² refinement"}</b><span>{event.breakpointModel.informativeSites.toLocaleString("en-US")} informative sites{event.breakpointModel.stateSwitches !== undefined ? ` · ${event.breakpointModel.stateSwitches} state switches` : ""}</span></div>}
        <button className="text-button" type="button" onClick={() => {
          const nextStart = event.end === alignment.length ? 0 : event.end;
          const nextEnd = event.start;
          onUpdate({
            start: nextStart,
            end: nextEnd,
            wraps: nextStart > nextEnd,
            confidenceStart: [nextStart, nextStart],
            confidenceEnd: [nextEnd, nextEnd],
            majorParent: event.minorParent,
            minorParent: event.majorParent,
            breakpointModel: { method: "manual", informativeSites: event.informativeSites },
          }, "Selected complementary circular arc");
        }}>↻ Use complementary circular arc</button>
        <button className="text-button" type="button" onClick={() => onUpdate({ majorParent: event.minorParent, minorParent: event.majorParent }, "Swapped parental assignments")}>⇄ Swap parental assignments</button>
        <div className="group-controls">
          <label><span>Event group</span><select value={event.groupId ?? ""} onChange={(value) => onUpdate({ groupId: value.target.value || null }, value.target.value ? `Assigned ${value.target.value}` : "Ungrouped event")}><option value="">Ungrouped</option>{groupIds.map((groupId) => <option value={groupId} key={groupId}>{groupId}</option>)}</select></label>
          <button type="button" onClick={() => onUpdate({ groupId: newGroupId }, `Created ${newGroupId}`)}>New group</button>
        </div>
      </div>
      <div className="inspector-section">
        <div className="section-heading"><h3>Method concordance</h3><span>{event.evidence.filter((item) => item.supported).length}/{event.evidence.length}</span></div>
        <div className="evidence-list">
          {event.evidence.length ? event.evidence.map((item) => <button type="button" key={item.method} onClick={() => onOpenMethod(item.method)} title={`Open ${item.method} result · ${item.statisticLabel}: ${item.statistic.toPrecision(4)} · ${item.calibration}`}><span className={item.supported ? "support-dot" : "support-dot muted"}/><b>{item.method}</b><code>{formatP(item.correctedP)}</code><i>→</i></button>) : <p>No method evidence yet. This is a manual or known-truth hypothesis.</p>}
        </div>
        {!event.evidenceStale && <button className="text-button" type="button" onClick={onRecalculate} disabled={recalculating}>{recalculating ? "Recalculating evidence…" : "↻ Recalculate exact edited hypothesis"}</button>}
        <div className="event-diagnostics"><span><b>{event.diagnostics.rateRatio.toFixed(2)}×</b> tract/background variable-site density</span><span><b>{(event.diagnostics.parentConflictRate * 100).toFixed(1)}%</b> parent-conflict sites</span></div>
      </div>
      {event.warnings.length > 0 && <div className="warning-box"><strong>Review flags</strong>{event.warnings.map((warning) => <span key={warning}>{warning}</span>)}</div>}
      <div className="inspector-section">
        <label className="note-label">Analyst note<textarea key={`${event.id}-${event.note}`} rows={3} defaultValue={event.note} onBlur={(value) => { if (value.target.value !== event.note) onUpdate({ note: value.target.value }, "Edited analyst note"); }} placeholder="Record why this hypothesis was accepted, changed, or rejected…" /></label>
        <details className="event-history"><summary>Audit trail · {event.history.length} entr{event.history.length === 1 ? "y" : "ies"}</summary><div>{[...event.history].reverse().map((entry) => <article key={entry.id}><b>{entry.action}</b><span>{entry.summary}</span><time dateTime={entry.timestamp}>{formatDateTime(entry.timestamp)}</time></article>)}</div></details>
      </div>
      <div className="decision-actions">
        <button type="button" className={event.decision === "rejected" ? "reject active" : "reject"} onClick={() => onDecision(event.decision === "rejected" ? "unreviewed" : "rejected")}><Icon name="close" size={16}/> Reject <kbd>X</kbd></button>
        <button type="button" className={event.decision === "accepted" ? "accept active" : "accept"} onClick={() => onDecision(event.decision === "accepted" ? "unreviewed" : "accepted")}><Icon name="check" size={16}/> Accept <kbd>A</kbd></button>
      </div>
      <p className="alpha-note">Scientific alpha. Treat inferred events as hypotheses until method-parity validation and independent review are complete.</p>
    </aside>
  );
}

export default function Home() {
  const [alignment, setAlignment] = useState<AlignmentData>(() => makeDemoAlignment());
  const [events, setEvents] = useState<RdpEvent[]>(() => [tutorialTruthEvent()]);
  const [selectedId, setSelectedId] = useState<string | null>("tutorial-known-truth");
  const [tab, setTab] = useState<Tab>("explore");
  const [options, setOptions] = useState<AnalysisOptions>(() => ({ ...DEFAULT_OPTIONS, mode: "query-reference" }));
  const [runState, setRunState] = useState<RunState>("idle");
  const [progress, setProgress] = useState(0);
  const [phase, setPhase] = useState("Ready");
  const [metrics, setMetrics] = useState<RunMetrics | null>(null);
  const [distanceMatrix, setDistanceMatrix] = useState<number[]>([]);
  const [pasteOpen, setPasteOpen] = useState(false);
  const [pasteText, setPasteText] = useState("");
  const [tutorialOpen, setTutorialOpen] = useState(false);
  const [examplesOpen, setExamplesOpen] = useState(false);
  const [tutorialStep, setTutorialStep] = useState(0);
  const [toast, setToast] = useState<string | null>(null);
  const [methodFilter, setMethodFilter] = useState("");
  const [selectedMethod, setSelectedMethod] = useState<MethodName>("RDP");
  const [exportScope, setExportScope] = useState<ExportScope>("accepted-fresh");
  const [loadedExampleId, setLoadedExampleId] = useState<string | null>("tutorial-virus");
  const [autoResolveStatus, setAutoResolveStatus] = useState<AutoResolveStatus>({ state: "idle", message: "Ready to preview", round: 0, processed: 0, accepted: 0, rejected: 0, reviewed: 0, rescans: 0 });
  const [auditLog, setAuditLog] = useState<ProjectAuditEntry[]>([{
    id: "project-audit-initial",
    timestamp: "2026-01-01T00:00:00.000Z",
    action: "Loaded tutorial dataset",
    summary: "Loaded deterministic synthetic alignment and one explicitly labeled known-truth hypothesis.",
    eventId: "tutorial-known-truth",
  }]);
  const [autosaveCandidate, setAutosaveCandidate] = useState<AutosaveRecord | null>(null);
  const [autosaveReady, setAutosaveReady] = useState(false);
  const [lastAutosavedAt, setLastAutosavedAt] = useState<number | null>(null);
  const [roleFilter, setRoleFilter] = useState("");
  const [dragging, setDragging] = useState(false);
  const [undoStack, setUndoStack] = useState<HistoryFrame[]>([]);
  const [redoStack, setRedoStack] = useState<HistoryFrame[]>([]);
  const workerRef = useRef<Worker | null>(null);
  const fileRef = useRef<HTMLInputElement>(null);
  const annotationRef = useRef<HTMLInputElement>(null);
  const jobRef = useRef(0);
  const partialEventsRef = useRef<RdpEvent[]>([]);
  const autoResolveSessionRef = useRef<AutoResolveSession | null>(null);
  const autoResolveRunnerRef = useRef<((inputEvents: RdpEvent[], session: AutoResolveSession) => void) | null>(null);

  const stats = useMemo(() => alignmentStats(alignment), [alignment]);
  const selectedIndex = events.findIndex((event) => event.id === selectedId);
  const selectedEvent = selectedIndex >= 0 ? events[selectedIndex] : null;
  const reviewed = events.filter((event) => event.decision !== "unreviewed").length;
  const groupIds = useMemo(() => [...new Set(events.map((event) => event.groupId).filter((value): value is string => Boolean(value)))].sort(), [events]);
  const hotspot = useMemo(() => breakpointHotspotTest(events, alignment.length, 48, 499, options.randomSeed), [alignment.length, events, options.randomSeed]);
  const matchingRoleRows = useMemo(() => alignment.sequences
    .map((record, index) => ({ record, index }))
    .filter(({ record }) => record.name.toLowerCase().includes(roleFilter.toLowerCase())), [alignment.sequences, roleFilter]);
  const exportableEvents = useMemo(() => events.filter((event) => {
    if (event.decision === "rejected") return false;
    if (exportScope === "accepted-fresh") return event.decision === "accepted" && !event.evidenceStale;
    if (exportScope === "all-fresh") return !event.evidenceStale;
    return true;
  }), [events, exportScope]);

  const showToast = useCallback((message: string) => {
    setToast(message);
    window.setTimeout(() => setToast(null), 3200);
  }, []);

  const appendAudit = useCallback((action: string, summary: string, eventId?: string, eventSnapshot?: string) => {
    setAuditLog((current) => [...current, {
      id: `project-audit-${Date.now()}-${current.length + 1}`,
      timestamp: new Date().toISOString(),
      action,
      summary,
      eventId,
      eventSnapshot,
    }]);
  }, []);

  const loadAlignment = useCallback((next: AlignmentData) => {
    workerRef.current?.terminate();
    autoResolveSessionRef.current = null;
    setAutoResolveStatus({ state: "idle", message: "Ready to preview", round: 0, processed: 0, accepted: 0, rejected: 0, reviewed: 0, rescans: 0 });
    setAlignment(next);
    setEvents([]);
    setUndoStack([]);
    setRedoStack([]);
    setSelectedId(null);
    setDistanceMatrix([]);
    setMetrics(null);
    setRunState("idle");
    setPhase("Ready to scan");
    setTab("explore");
    setLoadedExampleId(null);
    setAuditLog([{
      id: `project-audit-load-${Date.now()}`,
      timestamp: new Date().toISOString(),
      action: "Loaded alignment",
      summary: `${next.sequences.length.toLocaleString("en-US")} sequences × ${next.length.toLocaleString("en-US")} sites loaded locally.`,
    }]);
    showToast(`${next.sequences.length} aligned sequences loaded locally`);
  }, [showToast]);

  const loadProject = useCallback((project: ReturnType<typeof parseProject>) => {
    workerRef.current?.terminate();
    autoResolveSessionRef.current = null;
    setAutoResolveStatus({ state: "idle", message: "Ready to preview", round: 0, processed: 0, accepted: 0, rejected: 0, reviewed: 0, rescans: 0 });
    setAlignment(project.alignment);
    setOptions(project.options);
    setEvents(project.events);
    setUndoStack([]);
    setRedoStack([]);
    setSelectedId(project.events[0]?.id ?? null);
    setDistanceMatrix(project.distance);
    setMetrics(project.metrics);
    setAuditLog(project.auditLog);
    setLoadedExampleId(null);
    setRunState(project.events.length ? "complete" : "idle");
    setPhase(project.events.length ? `Restored · ${project.events.length} event${project.events.length === 1 ? "" : "s"}` : "Ready to scan");
    setTab("explore");
    showToast(`Project restored: ${project.alignment.sequences.length} sequences and ${project.events.length} events`);
  }, [showToast]);

  const loadExample = useCallback((example: ExampleDataset) => {
    setExamplesOpen(false);
    setPhase(`Generating ${example.title}`);
    setRunState("idle");
    window.setTimeout(() => {
      const next = example.generate();
      loadAlignment(next);
      setLoadedExampleId(example.id);
      setOptions((current) => ({ ...DEFAULT_OPTIONS, ...example.recommendedOptions, methods: [...current.methods] }));
      if (example.id === "tutorial-virus") {
        const truth = tutorialTruthEvent();
        setEvents([truth]);
        setSelectedId(truth.id);
      }
      setAuditLog([{
        id: `project-audit-example-${Date.now()}`,
        timestamp: new Date().toISOString(),
        action: "Generated example dataset",
        summary: `${example.title}: ${example.sequenceCount.toLocaleString("en-US")} synthetic sequences × ${example.length.toLocaleString("en-US")} sites.`,
      }]);
      showToast(`${example.title} generated locally · ready to scan`);
    }, 0);
  }, [loadAlignment, showToast]);

  const restoreAutosave = useCallback(() => {
    if (!autosaveCandidate) return;
    try {
      loadProject(parseProject(JSON.stringify(autosaveCandidate.project)));
      setAutosaveCandidate(null);
      setLastAutosavedAt(autosaveCandidate.savedAt);
      appendAudit("Restored browser autosave", `Recovered the local project saved ${formatDateTime(autosaveCandidate.savedAt)}.`);
    } catch (error) {
      showToast(error instanceof Error ? error.message : "The local autosave could not be restored.");
    }
  }, [appendAudit, autosaveCandidate, loadProject, showToast]);

  const dismissAutosave = useCallback(() => {
    setAutosaveCandidate(null);
    void clearAutosave();
    showToast("Previous autosave dismissed; the current project will become the new checkpoint");
  }, [showToast]);

  const readFile = useCallback(async (file: File) => {
    try {
      const text = await file.text();
      if (/\.(json|rdpweb)$/i.test(file.name) || text.trimStart().startsWith("{")) {
        loadProject(parseProject(text));
      } else {
        loadAlignment(parseAlignment(text, file.name));
      }
    } catch (error) {
      showToast(error instanceof Error ? error.message : "Could not parse the file.");
    }
  }, [loadAlignment, loadProject, showToast]);

  const readAnnotationFile = useCallback(async (file: File) => {
    try {
      const features = parseGenomeAnnotations(await file.text(), file.name, alignment.length);
      setAlignment((current) => ({ ...current, features }));
      showToast(`${features.length.toLocaleString("en-US")} annotations mapped to alignment coordinates`);
    } catch (error) {
      showToast(error instanceof Error ? error.message : "Could not parse annotations.");
    }
  }, [alignment.length, showToast]);

  const launchAnalysis = useCallback((config: AnalysisLaunchConfig = {}) => {
    const excludedTargets = config.excludedTargets ?? [];
    const excludedParents = config.excludedParents ?? [];
    const retainedEvents = config.retainedEvents ?? [];
    const filterResults = (candidateEvents: RdpEvent[]) => config.filterResolvedAgainst?.length
      ? filterResolvedEventDuplicates(candidateEvents, config.filterResolvedAgainst, alignment.length)
      : candidateEvents;
    workerRef.current?.terminate();
    const worker = new Worker(new URL("rdp-worker.js", document.baseURI), {
      type: "module",
      name: "rdp-analysis",
    });
    workerRef.current = worker;
    const jobId = Date.now();
    jobRef.current = jobId;
    setRunState("running");
    partialEventsRef.current = [];
    setProgress(0);
    setPhase("Preparing WebAssembly engine");
    setMetrics(null);
    worker.onmessage = (message) => {
      const payload = message.data;
      if (payload.jobId !== jobRef.current) return;
      if (payload.type === "progress") {
        setProgress(payload.progress);
        setPhase(payload.phase);
      } else if (payload.type === "partial") {
        const partialEvents = filterResults(payload.events);
        partialEventsRef.current = retainedEvents.length ? [...retainedEvents, ...partialEvents] : partialEvents;
      } else if (payload.type === "result") {
        partialEventsRef.current = [];
        const resultEvents = filterResults(payload.events);
        const nextEvents = retainedEvents.length
          ? [...retainedEvents, ...resultEvents]
          : resultEvents;
        setEvents(nextEvents);
        if (config.resetEditHistory !== false) {
          setUndoStack([]);
          setRedoStack([]);
        }
        setSelectedId(nextEvents[0]?.id ?? null);
        setDistanceMatrix(payload.distance);
        setMetrics({
          elapsedMs: payload.elapsedMs,
          comparisons: payload.comparisons,
          engine: payload.engine,
          matrixMode: payload.matrixMode,
          parentSamples: payload.parentSamples,
          timing: payload.timing,
          diagnostics: payload.diagnostics,
        });
        setProgress(1);
        setPhase(`Complete · ${nextEvents.length} event${nextEvents.length === 1 ? "" : "s"}`);
        setRunState("complete");
        appendAudit(config.auditAction ?? (retainedEvents.length ? "Completed unresolved-sequence rescan" : "Completed full scan"), `${config.auditContext ? `${config.auditContext} ` : ""}${payload.comparisons.toLocaleString("en-US")} triplets tested; ${resultEvents.length} new consensus hypotheses retained${payload.events.length !== resultEvents.length ? ` after suppressing ${payload.events.length - resultEvents.length} duplicate of an accepted event` : ""}${retainedEvents.length ? ` alongside ${retainedEvents.length} preserved hypotheses` : ""}.`);
        worker.terminate();
        showToast(resultEvents.length ? `Scan complete: ${resultEvents.length} new consensus event${resultEvents.length === 1 ? "" : "s"}` : retainedEvents.length ? "Rescan complete: no additional events passed the filters" : "Scan complete: no events passed the current filters");
        if (config.onComplete) window.setTimeout(() => config.onComplete?.(nextEvents), 0);
      } else if (payload.type === "error") {
        setRunState("error");
        setPhase(payload.message);
        const autoSession = autoResolveSessionRef.current;
        if (autoSession) {
          autoResolveSessionRef.current = null;
          setAutoResolveStatus({ state: "error", message: payload.message, round: autoSession.round, processed: autoSession.processed, accepted: autoSession.accepted, rejected: autoSession.rejected, reviewed: autoSession.reviewed, rescans: autoSession.rescans });
        }
        worker.terminate();
        showToast(payload.message);
      }
    };
    worker.onerror = () => {
      setRunState("error");
      setPhase("Analysis worker failed");
      const autoSession = autoResolveSessionRef.current;
      if (autoSession) {
        autoResolveSessionRef.current = null;
        setAutoResolveStatus({ state: "error", message: "The impacted-target rescan failed.", round: autoSession.round, processed: autoSession.processed, accepted: autoSession.accepted, rejected: autoSession.rejected, reviewed: autoSession.reviewed, rescans: autoSession.rescans });
      }
      showToast("The local analysis worker stopped unexpectedly.");
    };
    worker.postMessage({ type: "analyze", jobId, alignment, options, excludedTargets, excludedParents });
  }, [alignment, appendAudit, options, showToast]);

  const runAnalysis = useCallback(() => {
    autoResolveSessionRef.current = null;
    setAutoResolveStatus({ state: "idle", message: "Ready to preview", round: 0, processed: 0, accepted: 0, rejected: 0, reviewed: 0, rescans: 0 });
    launchAnalysis();
  }, [launchAnalysis]);

  const runIterativeAnalysis = useCallback(() => {
    autoResolveSessionRef.current = null;
    setAutoResolveStatus({ state: "idle", message: "Ready to preview", round: 0, processed: 0, accepted: 0, rejected: 0, reviewed: 0, rescans: 0 });
    const retained = events.filter((event) => event.decision === "accepted" && !event.evidenceStale);
    if (!retained.length) {
      showToast("Accept and recalculate at least one event before rescanning unresolved sequences");
      return;
    }
    launchAnalysis({ excludedTargets: [...new Set(retained.map((event) => event.recombinant))], retainedEvents: retained });
  }, [events, launchAnalysis, showToast]);

  const continueAutoResolve = useCallback((inputEvents: RdpEvent[], session: AutoResolveSession) => {
    if (inputEvents.length === 0) {
      autoResolveSessionRef.current = null;
      setAutoResolveStatus({ state: "complete", message: "No unresolved hypotheses remain.", round: session.round, processed: session.processed, accepted: session.accepted, rejected: session.rejected, reviewed: session.reviewed, rescans: session.rescans });
      return;
    }
    const plan = planAutoResolution(inputEvents, alignment.length, session.settings);
    const barrier = plan.barriers[0];
    const throughEventIndex = barrier?.afterEventIndex ?? inputEvents.length - 1;
    const timestamp = new Date().toISOString();
    const applied = applyAutoResolutionPlan(inputEvents, plan, throughEventIndex, session.profileLabel, timestamp);
    const processedNow = plan.entries.filter((entry) => entry.eventIndex <= throughEventIndex && entry.recommendation !== "keep").map((entry) => entry.eventId);
    const newlyReviewed = plan.entries.filter((entry) => entry.eventIndex <= throughEventIndex && entry.recommendation === "review" && !session.processedEventIds.includes(entry.eventId)).length;
    const processedEventIds = [...new Set([...session.processedEventIds, ...processedNow])];
    const nextSession: AutoResolveSession = {
      ...session,
      processedEventIds,
      processed: processedEventIds.length,
      accepted: session.accepted + applied.accepted,
      rejected: session.rejected + applied.rejected,
      reviewed: session.reviewed + newlyReviewed,
    };
    setEvents(applied.events);
    if (applied.changedIndexes.length) setSelectedId(applied.events[applied.changedIndexes[0]]?.id ?? selectedId);
    appendAudit("Auto-resolved reconstruction queue", `${session.profileLabel} profile processed ${processedNow.length} eligible hypotheses through E${throughEventIndex + 1}: ${applied.accepted} accepted, ${applied.rejected} rejected, and ${newlyReviewed} newly held for analyst review.${barrier ? ` A ${barrier.risk}/100 dependency-risk barrier paused the queue.` : " No dependency-triggered rescan barrier remained."}`);

    if (!barrier) {
      autoResolveSessionRef.current = null;
      setAutoResolveStatus({ state: "complete", message: applied.changedIndexes.length ? "Queue resolved as far as the configured evidence gates allow." : "No further automatic decisions pass the configured gates.", round: nextSession.round, processed: nextSession.processed, accepted: nextSession.accepted, rejected: nextSession.rejected, reviewed: nextSession.reviewed, rescans: nextSession.rescans });
      showToast(`Auto-resolve complete · ${nextSession.accepted} accepted, ${nextSession.rejected} rejected, ${nextSession.reviewed} held`);
      return;
    }

    if (session.rescans >= session.settings.maximumRescanRounds) {
      autoResolveSessionRef.current = null;
      setAutoResolveStatus({ state: "paused", message: `Paused at E${barrier.afterEventIndex + 1}: the configured rescan-round cap was reached.`, round: nextSession.round, processed: nextSession.processed, accepted: nextSession.accepted, rejected: nextSession.rejected, reviewed: nextSession.reviewed, rescans: nextSession.rescans });
      appendAudit("Paused auto-resolution", `Stopped at dependency barrier E${barrier.afterEventIndex + 1} because the maximum of ${session.settings.maximumRescanRounds} rescan rounds was reached.`);
      showToast("Auto-resolve paused at the configured rescan cap");
      return;
    }

    const targetPlan = rescanTargetsForBarrier(barrier, alignment.sequences.length, session.settings);
    const acceptedEvents = applied.events.filter((event) => event.decision === "accepted" && !event.evidenceStale);
    const acceptedTargets = new Set(acceptedEvents.map((event) => event.recombinant));
    const targetIndexes = targetPlan.targetIndexes;
    if (targetIndexes.length === 0) {
      autoResolveSessionRef.current = null;
      setAutoResolveStatus({ state: "complete", message: "The predicted downstream targets are already resolved; no rescan was needed.", round: nextSession.round, processed: nextSession.processed, accepted: nextSession.accepted, rejected: nextSession.rejected, reviewed: nextSession.reviewed, rescans: nextSession.rescans });
      return;
    }
    const targetSet = new Set(targetIndexes);
    const preservedEvents = applied.events.filter((event) => {
      if (event.decision === "accepted" && !event.evidenceStale) return true;
      if (event.decision === "rejected") return true;
      return !targetSet.has(event.recombinant);
    });
    const excludedTargets = Array.from({ length: alignment.sequences.length }, (_, index) => index)
      .filter((index) => !targetSet.has(index));
    const rescanningSession = { ...nextSession, rescans: nextSession.rescans + 1 };
    autoResolveSessionRef.current = rescanningSession;
    setAutoResolveStatus({ state: "rescanning", message: `${targetPlan.scope === "full" ? "Broad" : "Targeted"} rescan after E${barrier.afterEventIndex + 1}: ${targetIndexes.length} downstream target${targetIndexes.length === 1 ? "" : "s"}.`, round: rescanningSession.round, processed: rescanningSession.processed, accepted: rescanningSession.accepted, rejected: rescanningSession.rejected, reviewed: rescanningSession.reviewed, rescans: rescanningSession.rescans });
    launchAnalysis({
      excludedTargets,
      excludedParents: [...acceptedTargets],
      retainedEvents: preservedEvents,
      filterResolvedAgainst: acceptedEvents,
      resetEditHistory: false,
      auditAction: "Completed auto-resolve dependency rescan",
      auditContext: `${targetPlan.scope === "full" ? "Full unresolved-target" : "Impacted-target"} rescan after E${barrier.afterEventIndex + 1}; ${targetIndexes.length} targets selected by ${barrier.risk}/100 accumulated dependency risk, with ${acceptedTargets.size} accepted mosaic parent prox${acceptedTargets.size === 1 ? "y" : "ies"} excluded from parent search.`,
      onComplete: (rescannedEvents) => {
        const activeSession = autoResolveSessionRef.current;
        if (!activeSession) return;
        const nextRound = { ...activeSession, round: activeSession.round + 1 };
        autoResolveSessionRef.current = nextRound;
        setAutoResolveStatus({ state: "resolving", message: "Rescan complete; rebuilding the remaining queue plan.", round: nextRound.round, processed: nextRound.processed, accepted: nextRound.accepted, rejected: nextRound.rejected, reviewed: nextRound.reviewed, rescans: nextRound.rescans });
        autoResolveRunnerRef.current?.(rescannedEvents, nextRound);
      },
    });
  }, [alignment.length, alignment.sequences.length, appendAudit, launchAnalysis, selectedId, showToast]);
  useEffect(() => {
    autoResolveRunnerRef.current = continueAutoResolve;
    return () => { autoResolveRunnerRef.current = null; };
  }, [continueAutoResolve]);

  const runAutoResolve = useCallback((settings: AutoResolveSettings, profileLabel: string) => {
    if (runState === "running" || events.length === 0) return;
    const session: AutoResolveSession = { settings: { ...settings }, profileLabel, round: 0, processedEventIds: [], processed: 0, accepted: 0, rejected: 0, reviewed: 0, rescans: 0 };
    setUndoStack((current) => [...current.slice(-99), { label: `Run ${profileLabel} auto-resolve`, events, selectedId }]);
    setRedoStack([]);
    autoResolveSessionRef.current = session;
    setAutoResolveStatus({ state: "resolving", message: `Processing the queue with the ${profileLabel} profile.`, round: 0, processed: 0, accepted: 0, rejected: 0, reviewed: 0, rescans: 0 });
    continueAutoResolve(events, session);
  }, [continueAutoResolve, events, runState, selectedId]);

  const cancelAnalysis = useCallback(() => {
    const partial = partialEventsRef.current;
    const autoSession = autoResolveSessionRef.current;
    if (autoSession) {
      autoResolveSessionRef.current = null;
      setAutoResolveStatus({ state: "paused", message: "Stopped by the analyst; completed decisions and recoverable partial results were preserved.", round: autoSession.round, processed: autoSession.processed, accepted: autoSession.accepted, rejected: autoSession.rejected, reviewed: autoSession.reviewed, rescans: autoSession.rescans });
    }
    workerRef.current?.terminate();
    workerRef.current = null;
    jobRef.current += 1;
    setRunState("idle");
    setProgress(0);
    if (partial.length) {
      setEvents(partial);
      setSelectedId(partial[0].id);
      setPhase(`Cancelled · recovered ${partial.length} partial candidate${partial.length === 1 ? "" : "s"}`);
      appendAudit("Cancelled scan with partial recovery", `${partial.length} candidate hypotheses were checkpointed without complete method calibration.`);
      showToast(`Scan stopped; ${partial.length} partial candidate${partial.length === 1 ? "" : "s"} recovered for recalculation`);
    } else {
      setPhase("Cancelled before any candidate checkpoint");
    }
    partialEventsRef.current = [];
  }, [appendAudit, showToast]);

  const updateSelected = useCallback((patch: Partial<RdpEvent>, action = "Edited event") => {
    if (!selectedId) return;
    if (autoResolveSessionRef.current) {
      showToast("Stop auto-resolve before editing an event");
      return;
    }
    const currentEvent = events.find((event) => event.id === selectedId);
    if (!currentEvent) return;
    const unchanged = Object.entries(patch).every(([key, value]) => Object.is(currentEvent[key as keyof RdpEvent], value));
    if (unchanged) return;
    const scientificFields = new Set(["recombinant", "majorParent", "minorParent", "start", "end", "wraps", "groupId"]);
    const makesEvidenceStale = Object.keys(patch).some((key) => scientificFields.has(key));
    const normalizedPatch: Partial<RdpEvent> = makesEvidenceStale ? {
      ...patch,
      confidenceStart: patch.confidenceStart ?? [patch.start ?? currentEvent.start, patch.start ?? currentEvent.start],
      confidenceEnd: patch.confidenceEnd ?? [patch.end ?? currentEvent.end, patch.end ?? currentEvent.end],
    } : patch;
    const changedFields = Object.keys(patch).filter((key) => key !== "breakpointModel").join(", ");
    const selectedOrder = events.findIndex((event) => event.id === selectedId);
    const downstreamWarning = `Review-order dependency: an earlier event (E${selectedOrder + 1}) changed; rescan unresolved signals before relying on this characterization.`;
    const nextEvents = events.map((event, eventIndex) => {
      if (event.id === selectedId) return {
        ...event,
        ...normalizedPatch,
        evidenceStale: normalizedPatch.evidenceStale ?? (event.evidenceStale || makesEvidenceStale),
        history: [...(event.history ?? []), {
          id: `audit-${Date.now()}-${event.history?.length ?? 0}`,
          timestamp: new Date().toISOString(),
          action,
          summary: changedFields ? `Changed ${changedFields}.` : "Updated the event hypothesis.",
        }],
      };
      if (makesEvidenceStale && eventIndex > selectedOrder && event.decision !== "rejected") return {
        ...event,
        evidenceStale: true,
        warnings: event.warnings.includes(downstreamWarning) ? event.warnings : [...event.warnings, downstreamWarning],
      };
      return event;
    });
    setUndoStack((current) => [...current.slice(-99), { label: action, events, selectedId }]);
    setRedoStack([]);
    setEvents(nextEvents);
    appendAudit(action, changedFields ? `Changed ${changedFields}.` : "Updated the event hypothesis.", selectedId);
  }, [appendAudit, events, selectedId, showToast]);

  const setEventDecision = useCallback((eventId: string, decision: EventDecision) => {
    if (autoResolveSessionRef.current) {
      showToast("Stop auto-resolve before changing a decision");
      return;
    }
    const currentEvent = events.find((event) => event.id === eventId);
    if (!currentEvent || currentEvent.decision === decision) return;
    const action = decision === "accepted" ? "Accepted event" : decision === "rejected" ? "Rejected event" : "Reset review decision";
    setUndoStack((current) => [...current.slice(-99), { label: action, events, selectedId }]);
    setRedoStack([]);
    setEvents((current) => current.map((event) => event.id === eventId ? {
      ...event,
      decision,
      history: [...event.history, { id: `audit-${Date.now()}-${event.history.length}`, timestamp: new Date().toISOString(), action, summary: `Decision changed to ${decision}.` }],
    } : event));
    setSelectedId(eventId);
    appendAudit(action, `Decision changed to ${decision}.`, eventId);
  }, [appendAudit, events, selectedId, showToast]);

  const recalculateSelected = useCallback(() => {
    if (!selectedEvent) return;
    if (autoResolveSessionRef.current) {
      showToast("Stop auto-resolve before recalculating an individual event");
      return;
    }
    workerRef.current?.terminate();
    const worker = new Worker(new URL("rdp-worker.js", document.baseURI), { type: "module", name: "rdp-event-recalculation" });
    workerRef.current = worker;
    const jobId = Date.now();
    jobRef.current = jobId;
    setRunState("running");
    setProgress(0.25);
    setPhase("Recalculating edited hypothesis");
    worker.onmessage = (message) => {
      const payload = message.data;
      if (payload.jobId !== jobRef.current) return;
      if (payload.type === "recalculated") {
        updateSelected(payload.patch, "Recalculated exact hypothesis");
        setMetrics((current) => current
          ? { ...current, diagnostics: payload.diagnostics }
          : { elapsedMs: payload.elapsedMs, comparisons: payload.patch.hypothesisTests ?? 1, engine: "WebAssembly manual-event recalculation", diagnostics: payload.diagnostics });
        setProgress(1);
        setRunState("complete");
        setPhase(`Recalculated · ${payload.elapsedMs.toFixed(1)} ms`);
        worker.terminate();
        showToast("Method evidence recalculated for the edited assignment and breakpoints");
      } else if (payload.type === "error") {
        setRunState("error");
        setPhase(payload.message);
        worker.terminate();
        showToast(payload.message);
      }
    };
    worker.onerror = () => {
      setRunState("error");
      setPhase("Event recalculation failed");
      showToast("The local recalculation worker stopped unexpectedly.");
    };
    worker.postMessage({ type: "recalculate", jobId, alignment, options, event: selectedEvent, comparisons: metrics?.comparisons ?? selectedEvent.hypothesisTests ?? 1 });
  }, [alignment, metrics?.comparisons, options, selectedEvent, showToast, updateSelected]);

  const undo = useCallback(() => {
    if (autoResolveSessionRef.current) {
      showToast("Stop auto-resolve before undoing the workflow");
      return;
    }
    const frame = undoStack.at(-1);
    if (!frame) return;
    setRedoStack((current) => [...current.slice(-99), { label: frame.label, events, selectedId }]);
    setUndoStack((current) => current.slice(0, -1));
    setEvents(frame.events);
    setSelectedId(frame.selectedId);
    appendAudit("Undo", `Restored state before: ${frame.label}.`, frame.selectedId ?? undefined);
  }, [appendAudit, events, selectedId, showToast, undoStack]);

  const redo = useCallback(() => {
    if (autoResolveSessionRef.current) {
      showToast("Stop auto-resolve before redoing the workflow");
      return;
    }
    const frame = redoStack.at(-1);
    if (!frame) return;
    setUndoStack((current) => [...current.slice(-99), { label: frame.label, events, selectedId }]);
    setRedoStack((current) => current.slice(0, -1));
    setEvents(frame.events);
    setSelectedId(frame.selectedId);
    appendAudit("Redo", `Reapplied: ${frame.label}.`, frame.selectedId ?? undefined);
  }, [appendAudit, events, redoStack, selectedId, showToast]);

  const createManualEvent = useCallback(() => {
    if (autoResolveSessionRef.current) {
      showToast("Stop auto-resolve before creating a manual event");
      return;
    }
    const start = Math.max(0, Math.floor(alignment.length * 0.25));
    const end = Math.max(start + 1, Math.floor(alignment.length * 0.5));
    const id = `manual-${Date.now()}`;
    const event: RdpEvent = {
      id,
      recombinant: 0,
      majorParent: 1,
      minorParent: 2,
      start,
      end,
      wraps: false,
      confidenceStart: [start, start],
      confidenceEnd: [end, end],
      breakpointModel: { method: "manual", informativeSites: 0 },
      evidence: [],
      chiSquare: 0,
      informativeSites: 0,
      decision: "unreviewed",
      warnings: ["Manual hypothesis; method evidence has not been calculated"],
      note: "",
      source: "manual",
      groupId: null,
      history: [{ id: `${id}-history-1`, timestamp: new Date().toISOString(), action: "Created manual event", summary: "Initialized an analyst-defined hypothesis." }],
      evidenceStale: true,
      diagnostics: { tractVariableDensity: 0, backgroundVariableDensity: 0, rateRatio: 1, parentConflictRate: 0, parentDiscriminatingSites: 0, diffuseIncompatibility: false },
    };
    setUndoStack((current) => [...current.slice(-99), { label: "Create manual event", events, selectedId }]);
    setRedoStack([]);
    setEvents([...events, event]);
    setSelectedId(id);
    setTab("explore");
    appendAudit("Created manual event", `Initialized manual hypothesis ${id}.`, id);
  }, [alignment.length, appendAudit, events, selectedId, showToast]);

  const duplicateSelected = useCallback(() => {
    if (!selectedEvent) return;
    if (autoResolveSessionRef.current) {
      showToast("Stop auto-resolve before duplicating an event");
      return;
    }
    const id = `manual-${Date.now()}`;
    const duplicate: RdpEvent = {
      ...selectedEvent,
      id,
      source: "manual",
      decision: "unreviewed",
      note: selectedEvent.note ? `${selectedEvent.note}\nDuplicated for an alternative interpretation.` : "Duplicated for an alternative interpretation.",
      history: [...selectedEvent.history, { id: `${id}-history-${selectedEvent.history.length + 1}`, timestamp: new Date().toISOString(), action: "Duplicated event", summary: `Copied from ${selectedEvent.id}.` }],
    };
    setUndoStack((current) => [...current.slice(-99), { label: "Duplicate event", events, selectedId }]);
    setRedoStack([]);
    setEvents([...events, duplicate]);
    setSelectedId(id);
    appendAudit("Duplicated event", `Copied ${selectedEvent.id} to ${id}.`, id);
  }, [appendAudit, events, selectedEvent, selectedId, showToast]);

  const deleteSelected = useCallback(() => {
    if (!selectedEvent) return;
    if (autoResolveSessionRef.current) {
      showToast("Stop auto-resolve before deleting an event");
      return;
    }
    const index = events.findIndex((event) => event.id === selectedEvent.id);
    const remaining = events.filter((event) => event.id !== selectedEvent.id);
    setUndoStack((current) => [...current.slice(-99), { label: "Delete event", events, selectedId }]);
    setRedoStack([]);
    setEvents(remaining);
    setSelectedId(remaining[Math.min(index, remaining.length - 1)]?.id ?? null);
    appendAudit("Deleted event", `Deleted ${selectedEvent.id}; complete event provenance remains in this project audit ledger.`, selectedEvent.id, JSON.stringify(selectedEvent));
  }, [appendAudit, events, selectedEvent, selectedId, showToast]);

  const navigateEvent = useCallback((direction: -1 | 1) => {
    if (!events.length) return;
    const next = Math.max(0, Math.min(events.length - 1, selectedIndex + direction));
    setSelectedId(events[next].id);
  }, [events, selectedIndex]);

  useEffect(() => {
    const listener = (keyEvent: globalThis.KeyboardEvent) => {
      const target = keyEvent.target as HTMLElement | null;
      if (target?.matches("input, textarea, select")) return;
      if (keyEvent.key.toLowerCase() === "j") navigateEvent(1);
      if (keyEvent.key.toLowerCase() === "k") navigateEvent(-1);
      if (keyEvent.key.toLowerCase() === "a" && selectedEvent) updateSelected({ decision: selectedEvent.decision === "accepted" ? "unreviewed" : "accepted" });
      if (keyEvent.key.toLowerCase() === "x" && selectedEvent) updateSelected({ decision: selectedEvent.decision === "rejected" ? "unreviewed" : "rejected" });
      if ((keyEvent.metaKey || keyEvent.ctrlKey) && keyEvent.key.toLowerCase() === "z") {
        keyEvent.preventDefault();
        if (keyEvent.shiftKey) redo(); else undo();
      }
      if ((keyEvent.metaKey || keyEvent.ctrlKey) && keyEvent.key === "Enter") runAnalysis();
    };
    window.addEventListener("keydown", listener);
    return () => window.removeEventListener("keydown", listener);
  }, [navigateEvent, redo, runAnalysis, selectedEvent, undo, updateSelected]);

  useEffect(() => () => workerRef.current?.terminate(), []);

  useEffect(() => {
    let cancelled = false;
    loadAutosave()
      .then((record) => {
        if (!cancelled && record?.project?.schema?.startsWith("rdp-web/")) setAutosaveCandidate(record);
      })
      .catch(() => undefined)
      .finally(() => { if (!cancelled) setAutosaveReady(true); });
    return () => { cancelled = true; };
  }, []);

  useEffect(() => {
    if (!autosaveReady || autosaveCandidate) return;
    const timeout = window.setTimeout(() => {
      const project: RdpProject = {
        schema: "rdp-web/0.5",
        alignment,
        options,
        events,
        metrics,
        distance: distanceMatrix,
        auditLog,
      };
      saveAutosave(project).then(setLastAutosavedAt).catch(() => undefined);
    }, 1_200);
    return () => window.clearTimeout(timeout);
  }, [alignment, auditLog, autosaveCandidate, autosaveReady, distanceMatrix, events, metrics, options]);

  const setOption = <K extends keyof AnalysisOptions>(key: K, value: AnalysisOptions[K]) => {
    setOptions((current) => ({ ...current, [key]: value }));
  };

  const toggleMethod = (method: MethodName) => {
    setOptions((current) => ({
      ...current,
      methods: current.methods.includes(method)
        ? current.methods.filter((item) => item !== method)
        : [...current.methods, method],
    }));
  };

  const autoGroupReferences = () => {
    setAlignment((current) => ({
      ...current,
      sequences: current.sequences.map((record) => {
        if (record.role === "query") return record;
        const inferred = record.name.replace(/(?:[-_.\s]?\d+)+$/u, "").replace(/[-_.\s]+$/u, "").trim();
        return { ...record, referenceGroup: inferred || "Reference" };
      }),
    }));
    appendAudit("Auto-grouped reference sequences", "Inferred reference-group labels from sequence-name prefixes.");
    showToast("Reference groups inferred from sequence names; edit any label inline");
  };

  const onDrop = (dragEvent: DragEvent<HTMLElement>) => {
    dragEvent.preventDefault();
    setDragging(false);
    const file = dragEvent.dataTransfer.files[0];
    if (file) readFile(file);
  };

  const exportResults = (kind: "json" | "csv") => {
    if (kind === "json") {
      download("rdp-web-project.rdpweb", serializeProject({ alignment, options, events, metrics, distance: distanceMatrix, auditLog }), "application/json");
      return;
    }
    const rows = [["event", "group", "recombinant", "major_parent", "minor_parent", "start", "end", "wraps_origin", "breakpoint_model", "evidence_stale", "methods", "best_corrected_p", "decision", "warnings", "audit_entries"]];
    events.forEach((event, index) => rows.push([
      String(index + 1),
      event.groupId ?? "",
      alignment.sequences[event.recombinant].name,
      alignment.sequences[event.majorParent].name,
      alignment.sequences[event.minorParent].name,
      String(event.start + 1),
      String(event.end),
      String(event.wraps),
      event.breakpointModel?.method ?? "raw",
      String(event.evidenceStale),
      event.evidence.filter((item) => item.supported).map((item) => item.method).join(";"),
      event.evidence.length ? String(Math.min(...event.evidence.map((item) => item.correctedP))) : "",
      event.decision,
      event.warnings.join(";"),
      String(event.history.length),
    ]));
    download("rdp-web-events.csv", rows.map((row) => row.map((cell) => `"${cell.replaceAll('"', '""')}"`).join(",")).join("\n"), "text/csv");
  };

  const exportClean = (mode: "remove" | "mask" | "mask-codon" | "split" | "partition") => {
    if (!exportableEvents.length) {
      showToast("No fresh hypotheses match the selected export safety scope");
      return;
    }
    const files = exportRecombinationFree(alignment, exportableEvents, mode);
    files.forEach((file, index) => window.setTimeout(() => download(file.filename, file.content), index * 120));
    appendAudit("Exported recombination-aware alignment", `${mode} export used ${exportableEvents.length} fresh ${exportScope} hypotheses.`);
  };

  return (
    <main
      className={`app-shell ${dragging ? "dragging" : ""}`}
      onDragOver={(event) => { event.preventDefault(); setDragging(true); }}
      onDragLeave={(event) => { if (!event.currentTarget.contains(event.relatedTarget as Node)) setDragging(false); }}
      onDrop={onDrop}
    >
      <header className="topbar">
        <div className="brand"><span className="brand-mark"><i/><i/><i/></span><div><strong>RDP <em>Web</em></strong><small>Recombination workbench</small></div></div>
        <div className="engine-status"><span className={runState === "running" ? "pulse" : ""}/><b>{runState === "running" ? phase : "Local WebAssembly"}</b><small>Sequences never leave this tab</small></div>
        <div className="top-actions">
          <div className="history-actions" aria-label="Event edit history">
            <button type="button" className="ghost-button" disabled={!undoStack.length} onClick={undo} title={undoStack.length ? `Undo ${undoStack.at(-1)?.label}` : "Nothing to undo"} aria-label="Undo event edit"><Icon name="undo"/></button>
            <button type="button" className="ghost-button" disabled={!redoStack.length} onClick={redo} title={redoStack.length ? `Redo ${redoStack.at(-1)?.label}` : "Nothing to redo"} aria-label="Redo event edit"><Icon name="redo"/></button>
          </div>
          <button type="button" className="ghost-button examples-button" onClick={() => setExamplesOpen(true)}>▦ Examples</button>
          <button type="button" className="ghost-button" onClick={() => setTutorialOpen(true)}><Icon name="help"/> Guide</button>
          <input ref={fileRef} hidden type="file" accept=".fa,.fasta,.fas,.fna,.aln,.phy,.phylip,.nex,.nexus,.txt,.json,.rdpweb" onChange={(event: ChangeEvent<HTMLInputElement>) => { const file = event.target.files?.[0]; if (file) readFile(file); event.target.value = ""; }}/>
          <input ref={annotationRef} hidden type="file" accept=".gff,.gff3,.gb,.gbk,.genbank,.bed,.txt" onChange={(event: ChangeEvent<HTMLInputElement>) => { const file = event.target.files?.[0]; if (file) readAnnotationFile(file); event.target.value = ""; }}/>
          <button type="button" className="open-button" onClick={() => fileRef.current?.click()}><Icon name="upload"/> Open data / project</button>
          {runState === "running" ? (
            <button type="button" className="run-button stop" onClick={cancelAnalysis}><Icon name="stop"/> Stop <span>{Math.round(progress * 100)}%</span></button>
          ) : (
            <button type="button" className="run-button" onClick={runAnalysis} disabled={options.methods.length === 0}><Icon name="run"/> Run scan <kbd>⌘↵</kbd></button>
          )}
        </div>
      </header>

      {runState === "running" && <div className="global-progress"><i style={{ width: `${progress * 100}%` }}/></div>}

      <div className="workbench">
        <aside className="sidebar">
          <div className="sidebar-scroll">
            <section className="dataset-card">
              <div className="dataset-heading"><span className="eyebrow">Dataset</span><div><button type="button" onClick={() => setExamplesOpen(true)}>Examples</button><button type="button" onClick={() => setPasteOpen(true)}>Paste</button></div></div>
              <h2 title={alignment.name}>{alignment.name}</h2>
              <div className="dataset-stats"><span><b>{alignment.sequences.length}</b> sequences</span><span><b>{alignment.length.toLocaleString("en-US")}</b> sites</span><span title={stats.sampled ? "Stratified estimate for responsive large-data loading" : "Exact count"}><b>{stats.sampled ? "≈" : ""}{stats.variableSites.toLocaleString("en-US")}</b> variable</span><span title={stats.sampled ? "Stratified sequence/site estimate" : "Exact mean"}><b>{stats.sampled ? "≈" : ""}{(stats.meanIdentity * 100).toFixed(1)}%</b> mean ID</span></div>
              <div className="quality-line"><i style={{ width: `${Math.max(2, 100 - ((stats.gaps + stats.ambiguities) / (alignment.length * alignment.sequences.length)) * 100)}%` }}/><span>{stats.gaps + stats.ambiguities === 0 ? "No gaps or ambiguities" : `${(stats.gaps + stats.ambiguities).toLocaleString("en-US")} gaps / ambiguities`}</span></div>
              <button type="button" className="example-link" onClick={() => setExamplesOpen(true)}>{loadedExampleId ? `Loaded example · ${EXAMPLE_DATASETS.find((example) => example.id === loadedExampleId)?.complexity ?? "synthetic"} · browse library` : "Browse synthetic examples and stress tests"}</button>
            </section>

            <section className="settings-section">
              <div className="sidebar-title"><h3>Scan design</h3><Icon name="settings" size={16}/></div>
              <Segmented value={options.mode} options={[{ value: "exploratory", label: "Exploratory" }, { value: "query-reference", label: "Query ↔ Ref" }]} onChange={(value) => setOption("mode", value)} />
              <p className="setting-help">{options.mode === "exploratory" ? "Every sequence is tested against locally plausible parents." : "Only Q sequences are tested against sequences marked R or B."}</p>
            </section>

            <section className="settings-section">
              <div className="sidebar-title"><h3>Primary methods</h3><span>{options.methods.length}/7</span></div>
              <div className="method-toggles">
                {PRIMARY_METHODS.map((method) => <label key={method}><input type="checkbox" checked={options.methods.includes(method)} onChange={() => toggleMethod(method)}/><i/><span><b>{method}</b><small>{METHOD_META[method].family}</small></span><a href={METHOD_META[method].citation} target="_blank" rel="noreferrer" onClick={(event) => event.stopPropagation()}>↗</a></label>)}
              </div>
            </section>

            <section className="settings-section controls">
              <div className="sidebar-title"><h3>Detection</h3><button type="button" onClick={() => setOptions(DEFAULT_OPTIONS)}>Defaults</button></div>
              <label><span>Window size <b>{options.window} nt</b></span><input type="range" min={30} max={600} step={10} value={options.window} onChange={(event) => setOption("window", Number(event.target.value))}/></label>
              <label><span>Candidate parents <b>{options.exhaustive ? "All" : options.candidateParents}</b></span><input type="range" min={3} max={20} step={1} disabled={options.exhaustive} value={options.candidateParents} onChange={(event) => setOption("candidateParents", Number(event.target.value))}/></label>
              <label className="select-label"><span>Multiple testing</span><select value={options.correction} onChange={(event) => setOption("correction", event.target.value as AnalysisOptions["correction"])}><option value="bonferroni">Bonferroni</option><option value="holm">Holm step-down</option><option value="none">None (local p)</option></select></label>
              <label><span>Minimum support <b>{options.minMethods} methods</b></span><input type="range" min={1} max={Math.max(1, options.methods.length)} value={Math.min(options.minMethods, Math.max(1, options.methods.length))} onChange={(event) => setOption("minMethods", Number(event.target.value))}/></label>
              <details><summary>Advanced controls</summary><div className="advanced-controls">
                <label className="switch"><input type="checkbox" checked={options.circular} onChange={(event) => setOption("circular", event.target.checked)}/><i/><span>Circular genome</span></label>
                <label className="switch"><input type="checkbox" checked={options.polishBreakpoints} onChange={(event) => setOption("polishBreakpoints", event.target.checked)}/><i/><span>HMM-polish breakpoints</span></label>
                <label className="switch"><input type="checkbox" checked={options.checkMisalignment} onChange={(event) => setOption("checkMisalignment", event.target.checked)}/><i/><span>Flag alignment artefacts</span></label>
                <label className="switch"><input type="checkbox" checked={options.exhaustive} onChange={(event) => setOption("exhaustive", event.target.checked)}/><i/><span>Exhaustive parent search</span></label>
                {options.mode === "query-reference" && <label className="switch"><input type="checkbox" checked={options.testReferences} onChange={(event) => setOption("testReferences", event.target.checked)}/><i/><span>Also test references as recombinants</span></label>}
                <label><span>Bootstrap replicates <b>{options.bootstrapReplicates}</b></span><input type="range" min={0} max={500} step={25} value={options.bootstrapReplicates} onChange={(event) => setOption("bootstrapReplicates", Number(event.target.value))}/></label>
                <label className="select-label"><span>Random seed</span><input type="number" min={0} max={4294967295} value={options.randomSeed} onChange={(event) => setOption("randomSeed", Number(event.target.value) >>> 0)}/></label>
                <label className="select-label"><span>Global α</span><select value={options.alpha} onChange={(event) => setOption("alpha", Number(event.target.value))}><option value={0.1}>0.10</option><option value={0.05}>0.05</option><option value={0.01}>0.01</option><option value={0.001}>0.001</option></select></label>
              </div></details>
            </section>
          </div>
          <div className="sidebar-footer"><span>{runState === "running" ? phase : metrics ? `${metrics.comparisons.toLocaleString("en-US")} triplets · ${(metrics.elapsedMs / 1000).toFixed(2)} s` : "Ready for local analysis"}</span><small>{lastAutosavedAt ? `Browser autosaved ${formatClockTime(lastAutosavedAt)}` : metrics?.engine ?? "Worker-isolated · autosave initializing"}</small></div>
        </aside>

        <section className="workspace">
          {autosaveCandidate && <div className="autosave-banner"><div><strong>Unsaved local project found</strong><span>{autosaveCandidate.project.alignment.name} · {autosaveCandidate.project.events.length} hypotheses · saved {formatDateTime(autosaveCandidate.savedAt)}</span></div><button type="button" onClick={dismissAutosave}>Start fresh</button><button type="button" className="restore" onClick={restoreAutosave}>Restore autosave</button></div>}
          <nav className="tabs" aria-label="Analysis views">
            {(["explore", "reconstruction", "trees", "alignment", "patterns", "export", "methods"] as Tab[]).map((item) => <button type="button" key={item} className={tab === item ? "active" : ""} onClick={() => setTab(item)}>{item === "explore" ? "Events & evidence" : item === "reconstruction" ? "Global reconstruction" : item === "trees" ? "⑂ Local trees" : item === "alignment" ? "Alignment" : item === "patterns" ? "Patterns & matrices" : item === "export" ? "Export" : "Methods & papers"}{item === "explore" && events.length > 0 && <span>{events.length}</span>}{item === "reconstruction" && events.some((event) => event.evidenceStale) && <span>!</span>}</button>)}
          </nav>

          <div className="workspace-content">
            {tab === "explore" && <>
              {selectedEvent && <div className="selected-context"><div><span>Selected hypothesis · E{selectedIndex + 1}</span><b>{alignment.sequences[selectedEvent.recombinant].name}</b></div><div><span>Major parent</span><b>{alignment.sequences[selectedEvent.majorParent].name}</b></div><div><span>Minor parent</span><b>{alignment.sequences[selectedEvent.minorParent].name}</b></div><div><span>Region</span><b>{formatEventRegion(selectedEvent, alignment.length)}</b></div><button type="button" onClick={() => setTab("trees")}>Compare trees <span>→</span></button></div>}
              <Panel title="Recombination map" action={<div className="map-summary"><span className="status-chip">{events.filter((event) => event.decision !== "rejected").length} visible events</span><span>{reviewed}/{events.length} reviewed</span>{events.some((event) => event.decision === "accepted" && !event.evidenceStale) && <button type="button" className="small-button" onClick={runIterativeAnalysis} title="Keep accepted hypotheses and exclude their recombinant sequences from the next target pass">↻ Rescan unresolved</button>}</div>}>
                <Overview alignment={alignment} events={events} selectedId={selectedId} onSelect={setSelectedId}/>
              </Panel>
              {selectedEvent ? <>
                <Panel title="Method-by-method result" action={<span className={`decision-label ${selectedEvent.evidenceStale ? "unreviewed" : selectedEvent.decision}`}>{selectedEvent.evidenceStale ? "evidence stale" : selectedEvent.decision}</span>}><MethodEvidencePanel alignment={alignment} event={selectedEvent} window={options.window} selectedMethod={selectedMethod} onSelectMethod={setSelectedMethod}/></Panel>
                <Panel title="Identity context & breakpoint editor" action={<div className="region-chip"><span>Context plot—not a method-specific test</span><b>{formatEventRegion(selectedEvent, alignment.length)}</b></div>}><EvidencePlot alignment={alignment} event={selectedEvent} window={options.window} circular={options.circular} onUpdate={updateSelected}/></Panel>
              </> : <div className="empty-state"><span className="empty-mark">∿</span><h2>No hypothesis selected</h2><p>Run the enabled methods, create a manual hypothesis, or load a truth-annotated example.</p><div className="empty-actions"><button type="button" className="run-inline" onClick={runAnalysis}><Icon name="run"/> Run local scan</button><button type="button" className="small-button" onClick={() => setExamplesOpen(true)}>Browse examples</button></div></div>}
              <Panel title="Event hypotheses" action={<div className="table-actions"><button className="small-button" type="button" onClick={createManualEvent}>＋ Manual event</button><button className="small-button" type="button" onClick={() => exportResults("csv")}><Icon name="download" size={14}/> CSV</button></div>}><EventTable alignment={alignment} events={events} selectedId={selectedId} onSelect={setSelectedId}/></Panel>
            </>}

            {tab === "reconstruction" && <ReconstructionWorkspace alignment={alignment} events={events} selectedId={selectedId} onSelect={(id, view = "explore") => { setSelectedId(id); setTab(view); }} onDecision={setEventDecision} onRescan={runIterativeAnalysis} rescanning={runState === "running"} autoResolveStatus={autoResolveStatus} onAutoResolve={runAutoResolve}/>}

            {tab === "trees" && <LocalTrees alignment={alignment} events={events} event={selectedEvent} onSelectEvent={setSelectedId}/>}

            {tab === "alignment" && <>
              <AlignmentViewer key={`${alignment.createdAt}-${selectedEvent?.id ?? "none"}`} alignment={alignment} event={selectedEvent}/>
              <AnnotationPanel alignment={alignment} onOpen={() => annotationRef.current?.click()}/>
              <Panel title="Sequence roles & reference groups" action={<div className="role-actions"><input className="role-filter" aria-label="Filter sequence roles" placeholder="Find a sequence…" value={roleFilter} onChange={(event) => setRoleFilter(event.target.value)}/><button type="button" className="small-button" onClick={autoGroupReferences}>Auto-group by name</button></div>}>
                <div className="role-help">Groups diversify the fast parent shortlist. References can also be tested as recombinants from Advanced controls.</div>
                <div className="role-list">{matchingRoleRows.slice(0, 500).map(({ record, index }) => <div key={`${record.name}-${index}`}><span className={`role-dot ${record.role ?? "both"}`}>{record.role === "query" ? "Q" : record.role === "reference" ? "R" : "B"}</span><b>{record.name}</b><input className="reference-group-input" aria-label={`Reference group for ${record.name}`} placeholder="Unassigned group" value={record.referenceGroup ?? ""} disabled={record.role === "query"} onChange={(value) => setAlignment((current) => ({ ...current, sequences: current.sequences.map((candidate, candidateIndex) => candidateIndex === index ? { ...candidate, referenceGroup: value.target.value || undefined } : candidate) }))}/><code>{record.sequence.length.toLocaleString("en-US")} nt</code><button type="button" onClick={() => setAlignment((current) => ({ ...current, sequences: current.sequences.map((candidate, candidateIndex) => candidateIndex === index ? { ...candidate, role: candidate.role === "query" ? "reference" : candidate.role === "reference" ? "both" : "query" } : candidate) }))}>{record.role === "query" ? "Query" : record.role === "reference" ? "Reference" : "Both"}</button></div>)}</div>
                {matchingRoleRows.length > 500 && <p className="role-overflow">Showing the first 500 matches. Refine the sequence-name filter to edit another role.</p>}
              </Panel>
            </>}

            {tab === "patterns" && <>
              <RdpPatternMatrices alignment={alignment} events={events} onSelect={setSelectedId}/>
              <Panel title="Breakpoint density" action={<span className="panel-caption">Uniform-location null · seed {options.randomSeed}</span>}>
                <div className="density-chart">{Array.from({ length: 48 }, (_, bin) => {
                  const start = (bin / 48) * alignment.length;
                  const end = ((bin + 1) / 48) * alignment.length;
                  const count = events.filter((event) => event.decision !== "rejected" && ((event.start >= start && event.start < end) || (event.end >= start && event.end < end))).length;
                  const max = Math.max(1, ...Array.from({ length: 48 }, (_, other) => events.filter((event) => event.decision !== "rejected" && ((event.start >= (other / 48) * alignment.length && event.start < ((other + 1) / 48) * alignment.length) || (event.end >= (other / 48) * alignment.length && event.end < ((other + 1) / 48) * alignment.length))).length));
                  return <i key={bin} title={`${Math.round(start) + 1}–${Math.round(end)}: ${count} breakpoints`} style={{ height: `${Math.max(2, (count / max) * 100)}%` }}/>} )}</div>
                <div className="density-axis"><span>1</span><span>{Math.round(alignment.length / 2).toLocaleString("en-US")}</span><span>{alignment.length.toLocaleString("en-US")} nt</span></div>
                <div className="hotspot-result"><span><b>{hotspot.observedMaximum}</b> maximum breakpoints / bin</span><span><b>{hotspot.expectedPerBin.toFixed(2)}</b> expected under uniform null</span><span><b>{formatP(hotspot.empiricalP)}</b> empirical hotspot p · {hotspot.replicates} replicates</span></div>
              </Panel>
              <Panel title="Pairwise sequence p-distance" action={<span className="panel-caption">Sequence × sequence · first 24 · canonical sites</span>}><DistanceMatrix alignment={alignment} matrix={distanceMatrix}/></Panel>
              <Panel title="False-positive challenge diagnostics" action={<span className="panel-caption">Review evidence · never a hidden veto</span>}><ChallengeDiagnostics diagnostics={metrics?.diagnostics} event={selectedEvent}/></Panel>
            </>}

            {tab === "export" && <>
              <div className="export-hero"><span className="eyebrow">Recombination-aware outputs</span><h1>Take a defensible dataset downstream.</h1><p>The safe default uses only accepted hypotheses whose method evidence is fresh. Project and CSV exports always preserve every hypothesis and audit entry.</p><div className="review-meter"><i style={{ width: `${events.length ? (reviewed / events.length) * 100 : 0}%` }}/><span>{reviewed} of {events.length} events reviewed · {events.filter((event) => event.evidenceStale).length} stale</span></div></div>
              <div className="export-scope"><div><span className="eyebrow">Safety scope for sequence-changing exports</span><h2>{exportableEvents.length} hypothesis{exportableEvents.length === 1 ? "" : "es"} will be applied</h2></div><Segmented value={exportScope} options={[{ value: "accepted-fresh", label: "Accepted + fresh" }, { value: "all-fresh", label: "All fresh" }, { value: "all-retained", label: "All retained" }]} onChange={setExportScope}/><p>{exportScope === "accepted-fresh" ? "Publication-safe default: ignores unreviewed, rejected, and edited-but-unrecalculated events." : exportScope === "all-fresh" ? "Exploratory: includes fresh accepted and unreviewed hypotheses, while excluding rejected and stale events." : "Unsafe override: includes unreviewed and stale hypotheses; rejected events remain excluded. The project export preserves this choice in its ledger."}</p></div>
              <div className="export-grid">
                <article><span className="export-icon">—</span><h3>Remove recombinants</h3><p>Exclude every sequence carrying an in-scope event.</p><button disabled={!exportableEvents.length} type="button" onClick={() => exportClean("remove")}><Icon name="download"/> FASTA</button></article>
                <article><span className="export-icon">N</span><h3>Mask recombinant tracts</h3><p>Replace in-scope recombinant fragments with N characters.</p><button disabled={!exportableEvents.length} type="button" onClick={() => exportClean("mask")}><Icon name="download"/> FASTA</button></article>
                <article><span className="export-icon">3</span><h3>Codon-aware masking</h3><p>Expand in-scope tracts to CDS phase-aligned codon boundaries before masking.</p><button disabled={!exportableEvents.length} type="button" onClick={() => exportClean("mask-codon")}><Icon name="download"/> FASTA</button></article>
                <article><span className="export-icon">÷</span><h3>Split mosaic sequences</h3><p>Disassemble recombinants at in-scope breakpoints.</p><button disabled={!exportableEvents.length} type="button" onClick={() => exportClean("split")}><Icon name="download"/> FASTA</button></article>
                <article><span className="export-icon">▤</span><h3>Breakpoint partitions</h3><p>Create non-overlapping sub-alignments bounded by in-scope breakpoints.</p><button disabled={!exportableEvents.length} type="button" onClick={() => exportClean("partition")}><Icon name="download"/> FASTA set</button></article>
              </div>
              <Panel title="Results & provenance"><><div className="provenance-actions"><button type="button" onClick={() => exportResults("json")}><Icon name="download"/> Restorable project <small>.rdpweb schema 0.5 · all events + immutable project ledger</small></button><button type="button" onClick={() => exportResults("csv")}><Icon name="download"/> Event table CSV <small>all hypotheses, one event per row</small></button><button type="button" onClick={() => download("rdp-web-input.fasta", toFasta(alignment.sequences))}><Icon name="download"/> Input FASTA <small>normalized alignment</small></button></div><div className="run-provenance"><span><b>{auditLog.length}</b> project audit entries</span>{metrics && <><span><b>{metrics.comparisons.toLocaleString("en-US")}</b> triplets</span><span><b>{metrics.elapsedMs.toFixed(1)} ms</b> wall time</span><span><b>{metrics.matrixMode ?? "exact"}</b> parent screen</span>{metrics.timing && <span><b>{metrics.timing.distanceMs.toFixed(1)} / {metrics.timing.scanMs.toFixed(1)} / {metrics.timing.statisticsMs.toFixed(1)} / {(metrics.timing.diagnosticsMs ?? 0).toFixed(1)} ms</b> distance / scan / evidence / diagnostics</span>}</>}</div><details className="project-ledger"><summary>Project audit ledger · {auditLog.length} entries</summary><div>{[...auditLog].reverse().slice(0, 100).map((entry) => <article key={entry.id}><div><b>{entry.action}</b>{entry.eventSnapshot && <em>event tombstone saved</em>}</div><span>{entry.summary}</span><time dateTime={entry.timestamp}>{formatDateTime(entry.timestamp)}</time></article>)}</div></details></></Panel>
            </>}

            {tab === "methods" && <>
              <div className="methods-heading"><div><span className="eyebrow">Scientific basis</span><h1>Methods, limits, and primary sources.</h1><p>This MIT-licensed implementation is clean-room: papers and public documentation define behavior; no GPL OpenRDP or proprietary RDP source is incorporated.</p></div><input aria-label="Filter methods and papers" placeholder="Filter methods or papers…" value={methodFilter} onChange={(event) => setMethodFilter(event.target.value)}/></div>
              <div className="fidelity-banner"><strong>Validation state · scientific alpha</strong><p>Method-specific WebAssembly kernels, seeded triplet bootstraps, exact bounded 3SEQ dynamic programming, circular events, a windowless two-state HMM, review workflow, and restorable projects are operational. Full BURT fitting, method-specific permutation parity, and RDP5 numerical validation are not yet established; results must be independently verified.</p><a href="https://github.com/PoonLab/OpenRDP" target="_blank" rel="noreferrer">Why OpenRDP code is not bundled ↗</a></div>
              <Panel title="Primary exploratory methods" action={<span className="panel-caption">Seven-method consensus</span>}><div className="method-catalog">{PRIMARY_METHODS.filter((method) => `${method} ${METHOD_META[method].family} ${METHOD_META[method].detail}`.toLowerCase().includes(methodFilter.toLowerCase())).map((method) => <article key={method}><div><span className={options.methods.includes(method) ? "support-dot" : "support-dot muted"}/><h3>{method}</h3><em>{METHOD_META[method].family}</em></div><p>{METHOD_META[method].detail}</p><a href={METHOD_META[method].citation} target="_blank" rel="noreferrer">Primary method paper ↗</a></article>)}</div></Panel>
              <Panel title="Secondary verification views"><div className="secondary-methods">{[
                ["BURT precursor", "Operational windowless two-state Viterbi polishing; full 2–20-state fitting remains on the parity track", "RDP5 Manual §8.13", "operational precursor"],
                ["LARD", "Likelihood-ratio breakpoint comparison", "Holmes et al. 1999"],
                ["PhylPro", "Local phylogenetic profiles", "Weiller 1998"],
                ["VisRD", "Quartet scanning", "Lemey et al. 2009"],
                ["TOPAL / DSS", "Difference-of-sums-of-squares scan", "McGuire & Wright 2000"],
                ["Distance plots", "Windowed evolutionary-distance contrast", "RDP5 Manual §8.11"],
                ["PHI / 4-gamete", "False-positive and homoplasy flags", "Bruen et al. 2006; McVean et al. 2002"],
              ].filter((item) => item.join(" ").toLowerCase().includes(methodFilter.toLowerCase())).map((item) => <div key={item[0]}><b>{item[0]}</b><span>{item[1]}</span><small>{item[2]}</small><em>{item[3] ?? "view / validation track"}</em></div>)}</div></Panel>
              <Panel title="RDP lineage & practical guidance"><div className="reference-list">{REFERENCES.filter((reference) => `${reference.title} ${reference.authors} ${reference.tag}`.toLowerCase().includes(methodFilter.toLowerCase())).map((reference) => <a href={reference.href} target="_blank" rel="noreferrer" key={reference.title}><span>{reference.year}</span><div><em>{reference.tag}</em><h3>{reference.title}</h3><p>{reference.authors} · {reference.note}</p></div><b>↗</b></a>)}</div></Panel>
            </>}
          </div>
        </section>

        <Inspector alignment={alignment} event={selectedEvent} onDecision={(decision) => updateSelected({ decision }, decision === "accepted" ? "Accepted event" : decision === "rejected" ? "Rejected event" : "Reset review decision")} onUpdate={updateSelected} onNavigate={navigateEvent} onDuplicate={duplicateSelected} onDelete={deleteSelected} onRecalculate={recalculateSelected} onOpenTrees={() => setTab("trees")} onOpenMethod={(method) => { setSelectedMethod(method); setTab("explore"); }} recalculating={runState === "running"} canPrevious={selectedIndex > 0} canNext={selectedIndex >= 0 && selectedIndex < events.length - 1} groupIds={groupIds}/>
      </div>

      {examplesOpen && <ExampleLibrary onClose={() => setExamplesOpen(false)} onLoad={loadExample}/>}
      {pasteOpen && <div className="modal-backdrop" role="presentation" onMouseDown={() => setPasteOpen(false)}><section className="modal paste-modal" role="dialog" aria-modal="true" aria-labelledby="paste-title" onMouseDown={(event) => event.stopPropagation()}><button className="modal-close" type="button" onClick={() => setPasteOpen(false)} aria-label="Close"><Icon name="close"/></button><span className="eyebrow">Local import</span><h2 id="paste-title">Paste an aligned dataset</h2><p>FASTA, CLUSTAL, sequential/interleaved PHYLIP, and NEXUS MATRIX are recognized. Sequences are parsed in your browser.</p><textarea autoFocus spellCheck={false} value={pasteText} onChange={(event) => setPasteText(event.target.value)} placeholder={">sequence_A\nACGTACGT…\n>sequence_B\nACGTTCGT…"}/><div className="modal-actions"><button type="button" className="ghost-button" onClick={() => setPasteOpen(false)}>Cancel</button><button type="button" className="run-button" onClick={() => { try { loadAlignment(parseAlignment(pasteText)); setPasteOpen(false); setPasteText(""); } catch (error) { showToast(error instanceof Error ? error.message : "Could not parse alignment"); } }}>Load alignment</button></div></section></div>}

      {tutorialOpen && <div className="modal-backdrop" role="presentation" onMouseDown={() => setTutorialOpen(false)}><section className="modal tutorial-modal" role="dialog" aria-modal="true" aria-labelledby="tutorial-title" onMouseDown={(event) => event.stopPropagation()}><button className="modal-close" type="button" onClick={() => setTutorialOpen(false)} aria-label="Close"><Icon name="close"/></button><div className="tutorial-rail">{["Load", "Scan", "Verify", "Reconstruct", "Export"].map((step, index) => <button type="button" key={step} className={index === tutorialStep ? "active" : index < tutorialStep ? "done" : ""} onClick={() => setTutorialStep(index)}><span>{index < tutorialStep ? "✓" : index + 1}</span>{step}</button>)}</div><div className="tutorial-body"><span className="eyebrow">RDP5-style workflow · {tutorialStep + 1}/5</span><h2 id="tutorial-title">{["Start from a defensible alignment.", "Screen with multiple signals.", "Treat every event as a hypothesis.", "Reconstruct the global mosaic in order.", "Remove only the signals you retained."][tutorialStep]}</h2><p>{[
        "Use homologous, pre-aligned nucleotide sequences. Inspect divergent and gap-rich regions: misalignment is a major source of false positives. Mark query/reference roles only when that design is biologically justified.",
        "The fast mode prunes parent candidates with a WebAssembly distance pass before triplet scans. Use exhaustive mode for small definitive analyses; require concordance across methods and retain multiple-testing correction.",
        "For each event, compare the parent-affinity alignment, breakpoint localization, parental assignments, all local trees, method p-values, and alignment quality. Accept, reject, or edit—then record the rationale.",
        "Review the strongest/earliest characterized event first. Group co-recombinant descendants, inspect possible overprinting and recombinant-parent dependencies, and rescan unresolved signals after any scientific edit.",
        "Choose whether to remove recombinant sequences, mask fragments, split mosaics, or partition the alignment. Export the restorable project with parameters, rejected alternatives, and review decisions for provenance.",
      ][tutorialStep]}</p><div className="tutorial-tip"><b>{["Current dataset", "Performance", "Keyboard", "Review order", "Reproducibility"][tutorialStep]}</b><span>{[
        `${alignment.sequences.length} sequences × ${alignment.length.toLocaleString("en-US")} sites are loaded; the included tutorial has a known tract at 783–1,538.`,
        "Heavy computation runs in a worker, so plots and controls stay responsive. Stop cancels immediately.",
        "J/K: next/previous event · A: accept · X: reject · ⌘/Ctrl+Enter: run.",
        "An earlier edit can invalidate later characterizations. The Global reconstruction tab marks that dependency and keeps stale events out of safe exports.",
        "The project JSON records the normalized alignment, options, evidence, warnings, notes, and decisions.",
      ][tutorialStep]}</span></div><div className="modal-actions"><button type="button" className="ghost-button" disabled={tutorialStep === 0} onClick={() => setTutorialStep((current) => Math.max(0, current - 1))}>Back</button>{tutorialStep < 4 ? <button type="button" className="run-button" onClick={() => setTutorialStep((current) => current + 1)}>Continue <span>→</span></button> : <button type="button" className="run-button" onClick={() => { setTutorialOpen(false); setTab("reconstruction"); }}>Open reconstruction</button>}</div></div></section></div>}

      {dragging && <div className="drop-overlay"><Icon name="upload" size={32}/><strong>Drop data or a saved project to open it locally</strong><span>FASTA · CLUSTAL · PHYLIP · NEXUS · RDPWEB</span></div>}
      {toast && <div className="toast" role="status">{toast}</div>}
    </main>
  );
}
