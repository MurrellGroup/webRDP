"use client";

import {
  ChangeEvent,
  DragEvent,
  ReactNode,
  useCallback,
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
  RdpEvent,
  alignmentStats,
  breakpointHotspotTest,
  buildLocalTree,
  demoEvent,
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

type Tab = "explore" | "alignment" | "matrices" | "export" | "methods";
type RunState = "idle" | "running" | "complete" | "error";

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

const METHOD_META: Record<MethodName, { family: string; detail: string; citation: string }> = {
  RDP: {
    family: "Identity / triplet",
    detail: "Informative-site triplet scan with binomial assessment of identity runs.",
    citation: "https://doi.org/10.1093/bioinformatics/16.6.562",
  },
  GENECONV: {
    family: "Substitution distribution",
    detail: "Unusually long fragments of pairwise identity relative to the alignment background.",
    citation: "https://doi.org/10.1006/viro.1999.0058",
  },
  BootScan: {
    family: "Phylogenetic / windowed",
    detail: "Windowed support for changes in the recombinant’s closest relative.",
    citation: "https://doi.org/10.1089/aid.2005.21.98",
  },
  MaxChi: {
    family: "Substitution distribution",
    detail: "Maximum chi-square contrasts across variable sites around candidate breakpoints.",
    citation: "https://doi.org/10.1007/BF00182389",
  },
  Chimaera: {
    family: "Substitution distribution",
    detail: "Two-state refinement of MaxChi-style breakpoint evidence in sequence triplets.",
    citation: "https://doi.org/10.1073/pnas.241370698",
  },
  SiScan: {
    family: "Site-category / permutation",
    detail: "Fast oriented category-Z evidence; outgroup permutations remain on the parity track.",
    citation: "https://doi.org/10.1093/bioinformatics/16.7.573",
  },
  "3Seq": {
    family: "Non-parametric triplet",
    detail: "Maximum-descent hypergeometric random walk with exact bounded dynamic-program calibration.",
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
  return (
    <section className={`panel ${className}`}>
      <div className="panel-title"><h2>{title}</h2>{action}</div>
      {children}
    </section>
  );
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
            onUpdate?.({ ...preview, breakpointModel: { method: "manual", informativeSites: event.informativeSites } }, "Dragged breakpoint handle");
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
            {Math.round(fraction * alignment.length).toLocaleString()}
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
  const displayedIndexes = [...new Set([
    ...events.map((event) => event.recombinant),
    ...Array.from({ length: Math.min(400, alignment.sequences.length) }, (_, index) => index),
  ])].sort((left, right) => left - right);
  return (
    <div className="overview">
      <div className="overview-axis">
        <span>1</span><span>{Math.round(alignment.length * 0.25).toLocaleString()}</span><span>{Math.round(alignment.length * 0.5).toLocaleString()}</span><span>{Math.round(alignment.length * 0.75).toLocaleString()}</span><span>{alignment.length.toLocaleString()} nt</span>
      </div>
      <div className="overview-scroll">
        {displayedIndexes.map((index) => {
          const sequence = alignment.sequences[index];
          const sequenceEvents = events.filter((event) => event.recombinant === index && event.decision !== "rejected");
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
        {displayedIndexes.length < alignment.sequences.length && <div className="overview-cap">Showing {displayedIndexes.length.toLocaleString()} sequences, including every sequence with a retained event. Use Sequence roles search for the rest.</div>}
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
  return (
    <div className="table-scroll">
      <table className="event-table">
        <thead><tr><th>Event</th><th>Recombinant</th><th>Minor parent</th><th>Region</th><th>Methods</th><th>Best corrected p</th><th>Review</th></tr></thead>
        <tbody>
          {events.length === 0 ? <tr><td colSpan={7} className="empty-cell">No events pass the current consensus threshold.</td></tr> : events.map((event, index) => {
            const supported = event.evidence.filter((item) => item.supported).length;
            const bestP = Math.min(...event.evidence.map((item) => item.correctedP));
            return (
              <tr key={event.id} className={event.id === selectedId ? "selected" : ""} onClick={() => onSelect(event.id)}>
                <td><span className="event-number" title={event.groupId ? `Grouped as ${event.groupId}` : "Ungrouped event"}>{index + 1}{event.groupId && <i/>}</span></td>
                <td>{alignment.sequences[event.recombinant]?.name ?? "—"}</td>
                <td>{alignment.sequences[event.minorParent]?.name ?? "—"}</td>
                <td className="mono">{formatEventRegion(event, alignment.length)}</td>
                <td><span className="method-count">{supported}/{event.evidence.length}</span></td>
                <td className="mono">{formatP(bestP)}</td>
                <td><span className={`decision-label ${event.decision}`}>{event.decision}</span></td>
              </tr>
            );
          })}
        </tbody>
      </table>
    </div>
  );
}

function AlignmentViewer({ alignment, event }: { alignment: AlignmentData; event: RdpEvent | null }) {
  const initial = event ? Math.max(0, event.start - 45) : 0;
  const [cursor, setCursor] = useState(initial);
  const width = 96;
  const end = Math.min(alignment.length, cursor + width);
  const visibleLimit = Math.min(200, alignment.sequences.length);
  const indexes = event
    ? [...new Set([event.recombinant, event.majorParent, event.minorParent, ...Array.from({ length: visibleLimit }, (_, index) => index)])]
    : Array.from({ length: visibleLimit }, (_, index) => index);
  const baseClass = (base: string) => `base base-${base === "-" ? "gap" : "ACGT".includes(base) ? base : "amb"}`;
  return (
    <Panel title="Nucleotide alignment" action={<span className="panel-caption">Sites {cursor + 1}–{end}</span>}>
      <div className="alignment-toolbar">
        <button type="button" className="small-button" onClick={() => setCursor(Math.max(0, cursor - width))}>← Previous</button>
        <input aria-label="Alignment position" type="range" min={0} max={Math.max(0, alignment.length - width)} value={cursor} onChange={(eventValue) => setCursor(Number(eventValue.target.value))} />
        <button type="button" className="small-button" onClick={() => setCursor(Math.min(alignment.length - width, cursor + width))}>Next →</button>
      </div>
      <div className="alignment-grid" role="region" aria-label="Scrollable sequence alignment">
        <div className="alignment-ruler"><span />{Array.from({ length: end - cursor }, (_, offset) => <i key={offset}>{(cursor + offset + 1) % 10 === 0 ? "·" : ""}</i>)}</div>
        {indexes.map((sequenceIndex) => {
          const record = alignment.sequences[sequenceIndex];
          const relation = event?.recombinant === sequenceIndex ? "Recombinant" : event?.majorParent === sequenceIndex ? "Major parent" : event?.minorParent === sequenceIndex ? "Minor parent" : "";
          return (
            <div className={`alignment-row ${relation ? "relevant" : ""}`} key={`${record.name}-${sequenceIndex}`}>
              <span title={`${record.name}${relation ? ` · ${relation}` : ""}`}><b>{record.name}</b><small>{relation}</small></span>
              {[...record.sequence.slice(cursor, end)].map((base, offset) => <i key={offset} className={baseClass(base)}>{base}</i>)}
            </div>
          );
        })}
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
      {features.length > 200 && <p className="annotation-overflow">Showing the first 200 of {features.length.toLocaleString()} features.</p>}
    </>}
  </Panel>;
}

function DistanceMatrix({ alignment, matrix }: { alignment: AlignmentData; matrix: number[] }) {
  const count = Math.min(24, alignment.sequences.length);
  const fullMatrix = matrix.length === alignment.sequences.length ** 2;
  const compactMatrix = matrix.length === count ** 2;
  const valueAt = (row: number, column: number) => {
    if (fullMatrix) return matrix[row * alignment.sequences.length + column] ?? 0;
    if (compactMatrix) return matrix[row * count + column] ?? 0;
    return 1 - pairwiseIdentitySampled(alignment.sequences[row].sequence, alignment.sequences[column].sequence);
  };
  const values = Array.from({ length: count ** 2 }, (_, flatIndex) => valueAt(Math.floor(flatIndex / count), flatIndex % count));
  const max = Math.max(0.001, ...values);
  return (
    <div className="matrix-layout">
      <div className="heatmap" style={{ gridTemplateColumns: `104px repeat(${count}, minmax(18px, 1fr))` }}>
        <span />
        {alignment.sequences.slice(0, count).map((record) => <span className="column-label" key={record.name} title={record.name}>{record.name.slice(0, 3)}</span>)}
        {alignment.sequences.slice(0, count).flatMap((record, row) => [
          <span className="row-label" key={`label-${record.name}`} title={record.name}>{record.name}</span>,
          ...alignment.sequences.slice(0, count).map((_, column) => {
            const value = values[row * count + column] ?? 0;
            const intensity = value / max;
            return <span key={`${row}-${column}`} className="heat-cell" title={`${record.name} × ${alignment.sequences[column].name}: ${(value * 100).toFixed(2)}% distance`} style={{ background: `color-mix(in srgb, #ff7050 ${Math.round(intensity * 88)}%, #f2f4f3)` }} />;
          }),
        ])}
      </div>
      <div className="matrix-key"><span>Identical</span><i /><span>{(max * 100).toFixed(1)}% divergent</span></div>
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
    ["Four-gamete incompatible", `${(diagnostics.fourGameteFraction * 100).toFixed(1)}%`, `${diagnostics.incompatibleSitePairs.toLocaleString()} / ${diagnostics.testedSitePairs.toLocaleString()} sampled site pairs`],
    ["Near / far contrast", diagnostics.proximityRatio.toFixed(2), "PHI-style proximity diagnostic; descriptive, not a PHI p-value"],
    ["Canonical coverage", `${((1 - diagnostics.ambiguityFraction) * 100).toFixed(2)}%`, `${diagnostics.sampledSequences} sequences × ${diagnostics.sampledBiallelicSites} biallelic sites sampled`],
    ["Tract rate ratio", event ? event.diagnostics.rateRatio.toFixed(2) : "—", event ? `${(event.diagnostics.tractVariableDensity * 100).toFixed(1)}% tract vs ${(event.diagnostics.backgroundVariableDensity * 100).toFixed(1)}% background variable sites` : "Select an event"],
    ["Parent-conflict rate", event ? `${(event.diagnostics.parentConflictRate * 100).toFixed(1)}%` : "—", event ? `${event.diagnostics.parentDiscriminatingSites.toLocaleString()} parent-discriminating tract sites` : "Select an event"],
  ];
  return <div className="diagnostic-grid">{values.map(([label, value, note]) => <article key={label}><span>{label}</span><b>{value}</b><small>{note}</small></article>)}</div>;
}

function TreeSvg({ root }: { root: NeighborJoiningNode }) {
  const leafCount = (node: NeighborJoiningNode): number => node.children?.length
    ? node.children.reduce((total, child) => total + leafCount(child), 0)
    : 1;
  const farthest = (node: NeighborJoiningNode, distance = 0): number => {
    const here = distance + node.length;
    return node.children?.length ? Math.max(...node.children.map((child) => farthest(child, here))) : here;
  };
  const width = 720;
  const leaves = leafCount(root);
  const height = Math.max(150, leaves * 22 + 24);
  const left = 18;
  const labelWidth = 180;
  const scale = (width - left - labelWidth) / Math.max(0.000001, farthest(root));
  const edges: Array<{ x1: number; y1: number; x2: number; y2: number }> = [];
  const labels: Array<{ x: number; y: number; name: string }> = [];
  let leaf = 0;
  function place(node: NeighborJoiningNode, parentDistance: number): { x: number; y: number } {
    const distance = parentDistance + node.length;
    const x = left + distance * scale;
    if (!node.children?.length) {
      const y = 14 + (leaf + 0.5) * ((height - 24) / Math.max(1, leaves));
      leaf += 1;
      labels.push({ x: x + 5, y, name: node.name ?? "unnamed" });
      return { x, y };
    }
    const children = node.children.map((child) => place(child, distance));
    const y = children.reduce((sum, child) => sum + child.y, 0) / children.length;
    const childYs = children.map((child) => child.y);
    edges.push({ x1: x, y1: Math.min(...childYs), x2: x, y2: Math.max(...childYs) });
    children.forEach((child) => edges.push({ x1: x, y1: child.y, x2: child.x, y2: child.y }));
    return { x, y };
  }
  place(root, 0);
  return <div className="tree-svg-wrap"><svg viewBox={`0 0 ${width} ${height}`} role="img" aria-label="Neighbor-joining tree">{edges.map((edge, index) => <line key={index} {...edge}/>) }{labels.map((label) => <text key={`${label.name}-${label.y}`} x={label.x} y={label.y + 3}>{label.name}</text>)}</svg></div>;
}

function LocalTrees({ alignment, event }: { alignment: AlignmentData; event: RdpEvent | null }) {
  const trees = useMemo(() => {
    if (!event) return null;
    const selected = [event.recombinant, event.majorParent, event.minorParent];
    for (let index = 0; index < alignment.sequences.length && selected.length < 14; index += 1) {
      if (!selected.includes(index)) selected.push(index);
    }
    const tract = eventSegments(event, alignment.length);
    const background: [number, number][] = event.wraps
      ? [[event.end, event.start]]
      : [[0, event.start], [event.end, alignment.length]].filter(([start, end]) => end > start) as [number, number][];
    return [
      { label: "Recombinant tract", filename: "rdp-local-tract.nwk", tree: buildLocalTree(alignment, selected, tract) },
      { label: "Background / flanks", filename: "rdp-local-background.nwk", tree: buildLocalTree(alignment, selected, background) },
    ];
  }, [alignment, event]);
  if (!trees) return <div className="empty-state compact">Select an event to infer local trees.</div>;
  return <div className="tree-comparison">{trees.map((item) => <article key={item.label}><div><b>{item.label}</b><button type="button" onClick={() => download(item.filename, item.tree.newick)}>Newick ↓</button></div><TreeSvg root={item.tree.root}/><code title={item.tree.newick}>{item.tree.newick}</code></article>)}</div>;
}

function BreakpointMatrix({ alignment, events, onSelect }: { alignment: AlignmentData; events: RdpEvent[]; onSelect: (id: string) => void }) {
  const visible = events.filter((event) => event.decision !== "rejected").slice(0, 32);
  if (!visible.length) return <div className="empty-state compact">Retained events will appear as a breakpoint-pair similarity matrix.</div>;
  const circularDistance = (left: number, right: number) => {
    const distance = Math.abs(left - right);
    return Math.min(distance, alignment.length - distance);
  };
  const similarity = (left: RdpEvent, right: RdpEvent) => {
    const direct = circularDistance(left.start, right.start) + circularDistance(left.end % alignment.length, right.end % alignment.length);
    const reversed = circularDistance(left.start, right.end % alignment.length) + circularDistance(left.end % alignment.length, right.start);
    return 1 - Math.min(1, Math.min(direct, reversed) / Math.max(1, alignment.length));
  };
  return <div className="breakpoint-matrix" style={{ gridTemplateColumns: `88px repeat(${visible.length}, 18px)` }}><span/>{visible.map((_, index) => <b key={`column-${index}`}>{index + 1}</b>)}{visible.flatMap((event, row) => [<span key={`row-${event.id}`} title={alignment.sequences[event.recombinant].name}>E{row + 1} · {alignment.sequences[event.recombinant].name}</span>, ...visible.map((other, column) => {
    const value = similarity(event, other);
    return <button type="button" key={`${event.id}-${other.id}`} title={`Events ${row + 1} × ${column + 1}: ${(value * 100).toFixed(1)}% breakpoint similarity`} onClick={() => onSelect(event.id)} style={{ background: `color-mix(in srgb, #168d79 ${Math.round(value * 88)}%, #f1f4f2)` }}/>
  })])}</div>;
}

function Inspector({ alignment, event, onDecision, onUpdate, onNavigate, onDuplicate, onDelete, onRecalculate, recalculating, canPrevious, canNext, groupIds }: {
  alignment: AlignmentData;
  event: RdpEvent | null;
  onDecision: (decision: EventDecision) => void;
  onUpdate: (patch: Partial<RdpEvent>, action?: string) => void;
  onNavigate: (direction: -1 | 1) => void;
  onDuplicate: () => void;
  onDelete: () => void;
  onRecalculate: () => void;
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
      {event.evidenceStale && <div className="stale-box"><strong>Evidence needs recalculation</strong><span>Assignments or breakpoints changed after the scan. The saved method statistics remain visible for comparison but no longer describe this edited hypothesis.</span><button type="button" onClick={onRecalculate} disabled={recalculating}>{recalculating ? "Recalculating…" : "Recalculate this exact hypothesis"}</button></div>}
      <div className="inspector-section">
        <div className="section-heading"><h3>Breakpoint interval</h3><span className="mono">{eventLength(event, alignment.length).toLocaleString()} nt{event.wraps ? " · ↻" : ""}</span></div>
        <div className="breakpoint-inputs">
          <label>Start<input type="number" min={event.wraps ? event.end + 2 : 1} max={event.wraps ? alignment.length : event.end} value={event.start + 1} onChange={(value) => {
            const requested = Number(value.target.value) - 1;
            const start = event.wraps
              ? Math.max(event.end + 1, Math.min(alignment.length - 1, requested))
              : Math.max(0, Math.min(event.end - 1, requested));
            onUpdate({ start, breakpointModel: { method: "manual", informativeSites: event.informativeSites } }, "Edited left breakpoint");
          }}/><small>{event.confidenceStart[0] + 1}–{event.confidenceStart[1]}</small></label>
          <label>End<input type="number" min={1} max={event.wraps ? Math.max(1, event.start - 1) : alignment.length} value={event.end} onChange={(value) => {
            const requested = Number(value.target.value);
            const end = event.wraps
              ? Math.max(1, Math.min(Math.max(1, event.start - 1), requested))
              : Math.max(event.start + 1, Math.min(alignment.length, requested));
            onUpdate({ end, breakpointModel: { method: "manual", informativeSites: event.informativeSites } }, "Edited right breakpoint");
          }}/><small>{event.confidenceEnd[0] + 1}–{event.confidenceEnd[1]}</small></label>
        </div>
        {event.breakpointModel && <div className="breakpoint-model"><b>{event.breakpointModel.method === "two-state-hmm" ? "Windowless two-state HMM" : event.breakpointModel.method === "manual" ? "Manually edited" : "Local χ² refinement"}</b><span>{event.breakpointModel.informativeSites.toLocaleString()} informative sites{event.breakpointModel.stateSwitches !== undefined ? ` · ${event.breakpointModel.stateSwitches} state switches` : ""}</span></div>}
        <button className="text-button" type="button" onClick={() => {
          const nextStart = event.end === alignment.length ? 0 : event.end;
          const nextEnd = event.start;
          onUpdate({
            start: nextStart,
            end: nextEnd,
            wraps: nextStart > nextEnd,
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
          {event.evidence.map((item) => <a key={item.method} href={METHOD_META[item.method].citation} target="_blank" rel="noreferrer" title={`${item.statisticLabel}: ${item.statistic.toPrecision(4)} · ${item.calibration}`}><span className={item.supported ? "support-dot" : "support-dot muted"}/><b>{item.method}</b><code>{formatP(item.correctedP)}</code></a>)}
        </div>
        {!event.evidenceStale && <button className="text-button" type="button" onClick={onRecalculate} disabled={recalculating}>{recalculating ? "Recalculating evidence…" : "↻ Recalculate exact edited hypothesis"}</button>}
        <div className="event-diagnostics"><span><b>{event.diagnostics.rateRatio.toFixed(2)}×</b> tract/background variable-site density</span><span><b>{(event.diagnostics.parentConflictRate * 100).toFixed(1)}%</b> parent-conflict sites</span></div>
      </div>
      {event.warnings.length > 0 && <div className="warning-box"><strong>Review flags</strong>{event.warnings.map((warning) => <span key={warning}>{warning}</span>)}</div>}
      <div className="inspector-section">
        <label className="note-label">Analyst note<textarea key={`${event.id}-${event.note}`} rows={3} defaultValue={event.note} onBlur={(value) => { if (value.target.value !== event.note) onUpdate({ note: value.target.value }, "Edited analyst note"); }} placeholder="Record why this hypothesis was accepted, changed, or rejected…" /></label>
        <details className="event-history"><summary>Audit trail · {event.history.length} entr{event.history.length === 1 ? "y" : "ies"}</summary><div>{[...event.history].reverse().map((entry) => <article key={entry.id}><b>{entry.action}</b><span>{entry.summary}</span><time dateTime={entry.timestamp}>{new Date(entry.timestamp).toLocaleString()}</time></article>)}</div></details>
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
  const [events, setEvents] = useState<RdpEvent[]>(() => [demoEvent()]);
  const [selectedId, setSelectedId] = useState<string | null>("example-1");
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
  const [tutorialStep, setTutorialStep] = useState(0);
  const [toast, setToast] = useState<string | null>(null);
  const [methodFilter, setMethodFilter] = useState("");
  const [roleFilter, setRoleFilter] = useState("");
  const [dragging, setDragging] = useState(false);
  const [undoStack, setUndoStack] = useState<HistoryFrame[]>([]);
  const [redoStack, setRedoStack] = useState<HistoryFrame[]>([]);
  const workerRef = useRef<Worker | null>(null);
  const fileRef = useRef<HTMLInputElement>(null);
  const annotationRef = useRef<HTMLInputElement>(null);
  const jobRef = useRef(0);

  const stats = useMemo(() => alignmentStats(alignment), [alignment]);
  const selectedIndex = events.findIndex((event) => event.id === selectedId);
  const selectedEvent = selectedIndex >= 0 ? events[selectedIndex] : null;
  const reviewed = events.filter((event) => event.decision !== "unreviewed").length;
  const groupIds = useMemo(() => [...new Set(events.map((event) => event.groupId).filter((value): value is string => Boolean(value)))].sort(), [events]);
  const hotspot = useMemo(() => breakpointHotspotTest(events, alignment.length, 48, 499, options.randomSeed), [alignment.length, events, options.randomSeed]);
  const matchingRoleRows = useMemo(() => alignment.sequences
    .map((record, index) => ({ record, index }))
    .filter(({ record }) => record.name.toLowerCase().includes(roleFilter.toLowerCase())), [alignment.sequences, roleFilter]);

  const showToast = useCallback((message: string) => {
    setToast(message);
    window.setTimeout(() => setToast(null), 3200);
  }, []);

  const loadAlignment = useCallback((next: AlignmentData) => {
    workerRef.current?.terminate();
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
    showToast(`${next.sequences.length} aligned sequences loaded locally`);
  }, [showToast]);

  const loadProject = useCallback((project: ReturnType<typeof parseProject>) => {
    workerRef.current?.terminate();
    setAlignment(project.alignment);
    setOptions(project.options);
    setEvents(project.events);
    setUndoStack([]);
    setRedoStack([]);
    setSelectedId(project.events[0]?.id ?? null);
    setDistanceMatrix(project.distance);
    setMetrics(project.metrics);
    setRunState(project.events.length ? "complete" : "idle");
    setPhase(project.events.length ? `Restored · ${project.events.length} event${project.events.length === 1 ? "" : "s"}` : "Ready to scan");
    setTab("explore");
    showToast(`Project restored: ${project.alignment.sequences.length} sequences and ${project.events.length} events`);
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
      showToast(`${features.length.toLocaleString()} annotations mapped to alignment coordinates`);
    } catch (error) {
      showToast(error instanceof Error ? error.message : "Could not parse annotations.");
    }
  }, [alignment.length, showToast]);

  const runAnalysis = useCallback(() => {
    workerRef.current?.terminate();
    const worker = new Worker(new URL("rdp-worker.js", document.baseURI), {
      type: "module",
      name: "rdp-analysis",
    });
    workerRef.current = worker;
    const jobId = Date.now();
    jobRef.current = jobId;
    setRunState("running");
    setProgress(0);
    setPhase("Preparing WebAssembly engine");
    setMetrics(null);
    worker.onmessage = (message) => {
      const payload = message.data;
      if (payload.jobId !== jobRef.current) return;
      if (payload.type === "progress") {
        setProgress(payload.progress);
        setPhase(payload.phase);
      } else if (payload.type === "result") {
        setEvents(payload.events);
        setUndoStack([]);
        setRedoStack([]);
        setSelectedId(payload.events[0]?.id ?? null);
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
        setPhase(`Complete · ${payload.events.length} event${payload.events.length === 1 ? "" : "s"}`);
        setRunState("complete");
        worker.terminate();
        showToast(payload.events.length ? `Scan complete: ${payload.events.length} consensus event${payload.events.length === 1 ? "" : "s"}` : "Scan complete: no events passed the current filters");
      } else if (payload.type === "error") {
        setRunState("error");
        setPhase(payload.message);
        worker.terminate();
        showToast(payload.message);
      }
    };
    worker.onerror = () => {
      setRunState("error");
      setPhase("Analysis worker failed");
      showToast("The local analysis worker stopped unexpectedly.");
    };
    worker.postMessage({ type: "analyze", jobId, alignment, options });
  }, [alignment, options, showToast]);

  const cancelAnalysis = useCallback(() => {
    workerRef.current?.terminate();
    workerRef.current = null;
    jobRef.current += 1;
    setRunState("idle");
    setProgress(0);
    setPhase("Cancelled");
  }, []);

  const updateSelected = useCallback((patch: Partial<RdpEvent>, action = "Edited event") => {
    if (!selectedId) return;
    const currentEvent = events.find((event) => event.id === selectedId);
    if (!currentEvent) return;
    const unchanged = Object.entries(patch).every(([key, value]) => Object.is(currentEvent[key as keyof RdpEvent], value));
    if (unchanged) return;
    const scientificFields = new Set(["recombinant", "majorParent", "minorParent", "start", "end", "wraps"]);
    const makesEvidenceStale = Object.keys(patch).some((key) => scientificFields.has(key));
    const changedFields = Object.keys(patch).filter((key) => key !== "breakpointModel").join(", ");
    const nextEvents = events.map((event) => event.id === selectedId ? {
      ...event,
      ...patch,
      evidenceStale: patch.evidenceStale ?? (event.evidenceStale || makesEvidenceStale),
      history: [...(event.history ?? []), {
        id: `audit-${Date.now()}-${event.history?.length ?? 0}`,
        timestamp: new Date().toISOString(),
        action,
        summary: changedFields ? `Changed ${changedFields}.` : "Updated the event hypothesis.",
      }],
    } : event);
    setUndoStack((current) => [...current.slice(-99), { label: action, events, selectedId }]);
    setRedoStack([]);
    setEvents(nextEvents);
  }, [events, selectedId]);

  const recalculateSelected = useCallback(() => {
    if (!selectedEvent) return;
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
          : { elapsedMs: payload.elapsedMs, comparisons: 1, engine: "WebAssembly manual-event recalculation", diagnostics: payload.diagnostics });
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
    worker.postMessage({ type: "recalculate", jobId, alignment, options, event: selectedEvent });
  }, [alignment, options, selectedEvent, showToast, updateSelected]);

  const undo = useCallback(() => {
    const frame = undoStack.at(-1);
    if (!frame) return;
    setRedoStack((current) => [...current.slice(-99), { label: frame.label, events, selectedId }]);
    setUndoStack((current) => current.slice(0, -1));
    setEvents(frame.events);
    setSelectedId(frame.selectedId);
  }, [events, selectedId, undoStack]);

  const redo = useCallback(() => {
    const frame = redoStack.at(-1);
    if (!frame) return;
    setUndoStack((current) => [...current.slice(-99), { label: frame.label, events, selectedId }]);
    setRedoStack((current) => current.slice(0, -1));
    setEvents(frame.events);
    setSelectedId(frame.selectedId);
  }, [events, redoStack, selectedId]);

  const createManualEvent = useCallback(() => {
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
  }, [alignment.length, events, selectedId]);

  const duplicateSelected = useCallback(() => {
    if (!selectedEvent) return;
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
  }, [events, selectedEvent, selectedId]);

  const deleteSelected = useCallback(() => {
    if (!selectedEvent) return;
    const index = events.findIndex((event) => event.id === selectedEvent.id);
    const remaining = events.filter((event) => event.id !== selectedEvent.id);
    setUndoStack((current) => [...current.slice(-99), { label: "Delete event", events, selectedId }]);
    setRedoStack([]);
    setEvents(remaining);
    setSelectedId(remaining[Math.min(index, remaining.length - 1)]?.id ?? null);
  }, [events, selectedEvent, selectedId]);

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

  const onDrop = (dragEvent: DragEvent<HTMLElement>) => {
    dragEvent.preventDefault();
    setDragging(false);
    const file = dragEvent.dataTransfer.files[0];
    if (file) readFile(file);
  };

  const exportResults = (kind: "json" | "csv") => {
    if (kind === "json") {
      download("rdp-web-project.rdpweb", serializeProject({ alignment, options, events, metrics, distance: distanceMatrix }), "application/json");
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
    const files = exportRecombinationFree(alignment, events, mode);
    files.forEach((file, index) => window.setTimeout(() => download(file.filename, file.content), index * 120));
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
              <div className="dataset-heading"><span className="eyebrow">Dataset</span><button type="button" onClick={() => setPasteOpen(true)}>Paste</button></div>
              <h2 title={alignment.name}>{alignment.name}</h2>
              <div className="dataset-stats"><span><b>{alignment.sequences.length}</b> sequences</span><span><b>{alignment.length.toLocaleString()}</b> sites</span><span title={stats.sampled ? "Stratified estimate for responsive large-data loading" : "Exact count"}><b>{stats.sampled ? "≈" : ""}{stats.variableSites.toLocaleString()}</b> variable</span><span title={stats.sampled ? "Stratified sequence/site estimate" : "Exact mean"}><b>{stats.sampled ? "≈" : ""}{(stats.meanIdentity * 100).toFixed(1)}%</b> mean ID</span></div>
              <div className="quality-line"><i style={{ width: `${Math.max(2, 100 - ((stats.gaps + stats.ambiguities) / (alignment.length * alignment.sequences.length)) * 100)}%` }}/><span>{stats.gaps + stats.ambiguities === 0 ? "No gaps or ambiguities" : `${(stats.gaps + stats.ambiguities).toLocaleString()} gaps / ambiguities`}</span></div>
              <button type="button" className="example-link" onClick={() => { setAlignment(makeDemoAlignment()); setEvents([demoEvent()]); setUndoStack([]); setRedoStack([]); setSelectedId("example-1"); setMetrics(null); setRunState("idle"); showToast("Synthetic positive-control dataset restored"); }}>Reset tutorial dataset</button>
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
                <label><span>Bootstrap replicates <b>{options.bootstrapReplicates}</b></span><input type="range" min={0} max={500} step={25} value={options.bootstrapReplicates} onChange={(event) => setOption("bootstrapReplicates", Number(event.target.value))}/></label>
                <label className="select-label"><span>Random seed</span><input type="number" min={0} max={4294967295} value={options.randomSeed} onChange={(event) => setOption("randomSeed", Number(event.target.value) >>> 0)}/></label>
                <label className="select-label"><span>Global α</span><select value={options.alpha} onChange={(event) => setOption("alpha", Number(event.target.value))}><option value={0.1}>0.10</option><option value={0.05}>0.05</option><option value={0.01}>0.01</option><option value={0.001}>0.001</option></select></label>
              </div></details>
            </section>
          </div>
          <div className="sidebar-footer"><span>{runState === "running" ? phase : metrics ? `${metrics.comparisons.toLocaleString()} triplets · ${(metrics.elapsedMs / 1000).toFixed(2)} s` : "Ready for local analysis"}</span><small>{metrics?.engine ?? "Worker-isolated"}</small></div>
        </aside>

        <section className="workspace">
          <nav className="tabs" aria-label="Analysis views">
            {(["explore", "alignment", "matrices", "export", "methods"] as Tab[]).map((item) => <button type="button" key={item} className={tab === item ? "active" : ""} onClick={() => setTab(item)}>{item === "explore" ? "Explore" : item === "alignment" ? "Alignment" : item === "matrices" ? "Trees & matrices" : item === "export" ? "Export" : "Methods & papers"}{item === "explore" && events.length > 0 && <span>{events.length}</span>}</button>)}
          </nav>

          <div className="workspace-content">
            {tab === "explore" && <>
              <Panel title="Recombination map" action={<div className="map-summary"><span className="status-chip">{events.filter((event) => event.decision !== "rejected").length} visible events</span><span>{reviewed}/{events.length} reviewed</span></div>}>
                <Overview alignment={alignment} events={events} selectedId={selectedId} onSelect={setSelectedId}/>
              </Panel>
              {selectedEvent ? <Panel title="Evidence across the alignment" action={<div className="region-chip"><span>{alignment.sequences[selectedEvent.recombinant].name}</span><b>{formatEventRegion(selectedEvent, alignment.length)}</b></div>}><EvidencePlot alignment={alignment} event={selectedEvent} window={options.window} circular={options.circular} onUpdate={updateSelected}/></Panel> : <div className="empty-state"><span className="empty-mark">∿</span><h2>No signal selected</h2><p>Run the enabled methods, or lower the consensus threshold if a scan returns no events.</p><button type="button" className="run-inline" onClick={runAnalysis}><Icon name="run"/> Run local scan</button></div>}
              <Panel title="Event hypotheses" action={<div className="table-actions"><button className="small-button" type="button" onClick={createManualEvent}>＋ Manual event</button><button className="small-button" type="button" onClick={() => exportResults("csv")}><Icon name="download" size={14}/> CSV</button></div>}><EventTable alignment={alignment} events={events} selectedId={selectedId} onSelect={setSelectedId}/></Panel>
            </>}

            {tab === "alignment" && <>
              <AlignmentViewer key={`${alignment.createdAt}-${selectedEvent?.id ?? "none"}`} alignment={alignment} event={selectedEvent}/>
              <AnnotationPanel alignment={alignment} onOpen={() => annotationRef.current?.click()}/>
              <Panel title="Sequence roles" action={<input className="role-filter" aria-label="Filter sequence roles" placeholder="Find a sequence…" value={roleFilter} onChange={(event) => setRoleFilter(event.target.value)}/>}>
                <div className="role-list">{matchingRoleRows.slice(0, 500).map(({ record, index }) => <div key={`${record.name}-${index}`}><span className={`role-dot ${record.role ?? "both"}`}>{record.role === "query" ? "Q" : record.role === "reference" ? "R" : "B"}</span><b>{record.name}</b><code>{record.sequence.length.toLocaleString()} nt</code><button type="button" onClick={() => setAlignment((current) => ({ ...current, sequences: current.sequences.map((candidate, candidateIndex) => candidateIndex === index ? { ...candidate, role: candidate.role === "query" ? "reference" : candidate.role === "reference" ? "both" : "query" } : candidate) }))}>{record.role === "query" ? "Query" : record.role === "reference" ? "Reference" : "Both"}</button></div>)}</div>
                {matchingRoleRows.length > 500 && <p className="role-overflow">Showing the first 500 matches. Refine the sequence-name filter to edit another role.</p>}
              </Panel>
            </>}

            {tab === "matrices" && <>
              <Panel title="Pairwise p-distance matrix" action={<span className="panel-caption">Uncorrected · canonical sites</span>}><DistanceMatrix alignment={alignment} matrix={distanceMatrix}/></Panel>
              <Panel title="Local topology contrast" action={<span className="panel-caption">Closest pair by region</span>}><TopologyCards alignment={alignment} event={selectedEvent}/></Panel>
              <Panel title="Local neighbor-joining trees" action={<span className="panel-caption">Up to 14 sequences · p-distance · Newick export</span>}><LocalTrees alignment={alignment} event={selectedEvent}/></Panel>
              <Panel title="Breakpoint density" action={<span className="panel-caption">Uniform-location null · seed {options.randomSeed}</span>}>
                <div className="density-chart">{Array.from({ length: 48 }, (_, bin) => {
                  const start = (bin / 48) * alignment.length;
                  const end = ((bin + 1) / 48) * alignment.length;
                  const count = events.filter((event) => event.decision !== "rejected" && ((event.start >= start && event.start < end) || (event.end >= start && event.end < end))).length;
                  const max = Math.max(1, ...Array.from({ length: 48 }, (_, other) => events.filter((event) => event.decision !== "rejected" && ((event.start >= (other / 48) * alignment.length && event.start < ((other + 1) / 48) * alignment.length) || (event.end >= (other / 48) * alignment.length && event.end < ((other + 1) / 48) * alignment.length))).length));
                  return <i key={bin} title={`${Math.round(start) + 1}–${Math.round(end)}: ${count} breakpoints`} style={{ height: `${Math.max(2, (count / max) * 100)}%` }}/>} )}</div>
                <div className="density-axis"><span>1</span><span>{Math.round(alignment.length / 2).toLocaleString()}</span><span>{alignment.length.toLocaleString()} nt</span></div>
                <div className="hotspot-result"><span><b>{hotspot.observedMaximum}</b> maximum breakpoints / bin</span><span><b>{hotspot.expectedPerBin.toFixed(2)}</b> expected under uniform null</span><span><b>{formatP(hotspot.empiricalP)}</b> empirical hotspot p · {hotspot.replicates} replicates</span></div>
              </Panel>
              <Panel title="Breakpoint pair matrix" action={<span className="panel-caption">Circular-aware · first 32 retained events</span>}><BreakpointMatrix alignment={alignment} events={events} onSelect={setSelectedId}/></Panel>
              <Panel title="False-positive challenge diagnostics" action={<span className="panel-caption">Review evidence · never a hidden veto</span>}><ChallengeDiagnostics diagnostics={metrics?.diagnostics} event={selectedEvent}/></Panel>
            </>}

            {tab === "export" && <>
              <div className="export-hero"><span className="eyebrow">Recombination-aware outputs</span><h1>Take a defensible dataset downstream.</h1><p>Exports use accepted and unreviewed hypotheses; rejected events are excluded. Review all events before publication-grade use.</p><div className="review-meter"><i style={{ width: `${events.length ? (reviewed / events.length) * 100 : 0}%` }}/><span>{reviewed} of {events.length} events reviewed</span></div></div>
              <div className="export-grid">
                <article><span className="export-icon">—</span><h3>Remove recombinants</h3><p>Exclude every sequence carrying a retained event.</p><button type="button" onClick={() => exportClean("remove")}><Icon name="download"/> FASTA</button></article>
                <article><span className="export-icon">N</span><h3>Mask recombinant tracts</h3><p>Replace retained recombinant fragments with N characters.</p><button type="button" onClick={() => exportClean("mask")}><Icon name="download"/> FASTA</button></article>
                <article><span className="export-icon">3</span><h3>Codon-aware masking</h3><p>Expand retained tracts to CDS phase-aligned codon boundaries before masking.</p><button type="button" onClick={() => exportClean("mask-codon")}><Icon name="download"/> FASTA</button></article>
                <article><span className="export-icon">÷</span><h3>Split mosaic sequences</h3><p>Disassemble recombinants at retained breakpoints.</p><button type="button" onClick={() => exportClean("split")}><Icon name="download"/> FASTA</button></article>
                <article><span className="export-icon">▤</span><h3>Breakpoint partitions</h3><p>Create non-overlapping sub-alignments bounded by breakpoints.</p><button type="button" onClick={() => exportClean("partition")}><Icon name="download"/> FASTA set</button></article>
              </div>
              <Panel title="Results & provenance"><><div className="provenance-actions"><button type="button" onClick={() => exportResults("json")}><Icon name="download"/> Restorable project <small>.rdpweb · alignment, options, evidence, decisions</small></button><button type="button" onClick={() => exportResults("csv")}><Icon name="download"/> Event table CSV <small>tidy, one event per row</small></button><button type="button" onClick={() => download("rdp-web-input.fasta", toFasta(alignment.sequences))}><Icon name="download"/> Input FASTA <small>normalized alignment</small></button></div>{metrics && <div className="run-provenance"><span><b>{metrics.comparisons.toLocaleString()}</b> triplets</span><span><b>{metrics.elapsedMs.toFixed(1)} ms</b> wall time</span><span><b>{metrics.matrixMode ?? "exact"}</b> parent screen</span>{metrics.timing && <span><b>{metrics.timing.distanceMs.toFixed(1)} / {metrics.timing.scanMs.toFixed(1)} / {metrics.timing.statisticsMs.toFixed(1)} / {(metrics.timing.diagnosticsMs ?? 0).toFixed(1)} ms</b> distance / scan / evidence / diagnostics</span>}</div>}</></Panel>
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

        <Inspector alignment={alignment} event={selectedEvent} onDecision={(decision) => updateSelected({ decision }, decision === "accepted" ? "Accepted event" : decision === "rejected" ? "Rejected event" : "Reset review decision")} onUpdate={updateSelected} onNavigate={navigateEvent} onDuplicate={duplicateSelected} onDelete={deleteSelected} onRecalculate={recalculateSelected} recalculating={runState === "running"} canPrevious={selectedIndex > 0} canNext={selectedIndex >= 0 && selectedIndex < events.length - 1} groupIds={groupIds}/>
      </div>

      {pasteOpen && <div className="modal-backdrop" role="presentation" onMouseDown={() => setPasteOpen(false)}><section className="modal paste-modal" role="dialog" aria-modal="true" aria-labelledby="paste-title" onMouseDown={(event) => event.stopPropagation()}><button className="modal-close" type="button" onClick={() => setPasteOpen(false)} aria-label="Close"><Icon name="close"/></button><span className="eyebrow">Local import</span><h2 id="paste-title">Paste an aligned dataset</h2><p>FASTA, CLUSTAL, sequential/interleaved PHYLIP, and NEXUS MATRIX are recognized. Sequences are parsed in your browser.</p><textarea autoFocus spellCheck={false} value={pasteText} onChange={(event) => setPasteText(event.target.value)} placeholder={">sequence_A\nACGTACGT…\n>sequence_B\nACGTTCGT…"}/><div className="modal-actions"><button type="button" className="ghost-button" onClick={() => setPasteOpen(false)}>Cancel</button><button type="button" className="run-button" onClick={() => { try { loadAlignment(parseAlignment(pasteText)); setPasteOpen(false); setPasteText(""); } catch (error) { showToast(error instanceof Error ? error.message : "Could not parse alignment"); } }}>Load alignment</button></div></section></div>}

      {tutorialOpen && <div className="modal-backdrop" role="presentation" onMouseDown={() => setTutorialOpen(false)}><section className="modal tutorial-modal" role="dialog" aria-modal="true" aria-labelledby="tutorial-title" onMouseDown={(event) => event.stopPropagation()}><button className="modal-close" type="button" onClick={() => setTutorialOpen(false)} aria-label="Close"><Icon name="close"/></button><div className="tutorial-rail">{["Load", "Scan", "Verify", "Export"].map((step, index) => <button type="button" key={step} className={index === tutorialStep ? "active" : index < tutorialStep ? "done" : ""} onClick={() => setTutorialStep(index)}><span>{index < tutorialStep ? "✓" : index + 1}</span>{step}</button>)}</div><div className="tutorial-body"><span className="eyebrow">RDP5-style workflow · {tutorialStep + 1}/4</span><h2 id="tutorial-title">{["Start from a defensible alignment.", "Screen with multiple signals.", "Treat every event as a hypothesis.", "Remove the signal you reviewed."][tutorialStep]}</h2><p>{[
        "Use homologous, pre-aligned nucleotide sequences. Inspect divergent and gap-rich regions: misalignment is a major source of false positives. Mark query/reference roles only when that design is biologically justified.",
        "The fast mode prunes parent candidates with a WebAssembly distance pass before triplet scans. Use exhaustive mode for small definitive analyses; require concordance across methods and retain multiple-testing correction.",
        "For each event, compare identity curves, breakpoint localization, parental assignments, local topologies, method p-values, and alignment quality. Accept, reject, or edit—then record the rationale.",
        "Choose whether to remove recombinant sequences, mask fragments, split mosaics, or partition the alignment. Export the project JSON with parameters and review decisions for provenance.",
      ][tutorialStep]}</p><div className="tutorial-tip"><b>{["Current dataset", "Performance", "Keyboard", "Reproducibility"][tutorialStep]}</b><span>{[
        `${alignment.sequences.length} sequences × ${alignment.length.toLocaleString()} sites are loaded; the included tutorial has a known tract at 783–1,538.`,
        "Heavy computation runs in a worker, so plots and controls stay responsive. Stop cancels immediately.",
        "J/K: next/previous event · A: accept · X: reject · ⌘/Ctrl+Enter: run.",
        "The project JSON records the normalized alignment, options, evidence, warnings, notes, and decisions.",
      ][tutorialStep]}</span></div><div className="modal-actions"><button type="button" className="ghost-button" disabled={tutorialStep === 0} onClick={() => setTutorialStep((current) => Math.max(0, current - 1))}>Back</button>{tutorialStep < 3 ? <button type="button" className="run-button" onClick={() => setTutorialStep((current) => current + 1)}>Continue <span>→</span></button> : <button type="button" className="run-button" onClick={() => { setTutorialOpen(false); setTab("explore"); }}>Open workbench</button>}</div></div></section></div>}

      {dragging && <div className="drop-overlay"><Icon name="upload" size={32}/><strong>Drop data or a saved project to open it locally</strong><span>FASTA · CLUSTAL · PHYLIP · NEXUS · RDPWEB</span></div>}
      {toast && <div className="toast" role="status">{toast}</div>}
    </main>
  );
}
