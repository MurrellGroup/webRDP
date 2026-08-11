import assert from "node:assert/strict";
import { readFile } from "node:fs/promises";
import test from "node:test";

const [page, styles, worker] = await Promise.all([
  readFile(new URL("../app/page.tsx", import.meta.url), "utf8"),
  readFile(new URL("../app/globals.css", import.meta.url), "utf8"),
  readFile(new URL("../public/rdp-worker.js", import.meta.url), "utf8"),
]);

test("workspace owns vertical scrolling and generic panels expand without leaking paint", () => {
  assert.match(styles, /\.workspace-content\s*\{[^}]*overflow-y:\s*auto/s);
  assert.match(styles, /\.workspace-content\s*\{[^}]*height:\s*0/s);
  const panelRule = styles.match(/\.panel\s*\{([^}]*)\}/)?.[1] ?? "";
  assert.doesNotMatch(panelRule, /max-height/);
  assert.match(panelRule, /overflow:\s*hidden/);
  assert.match(panelRule, /isolation:\s*isolate/);
  const panelBodyRule = styles.match(/\.panel-body\s*\{([^}]*)\}/)?.[1] ?? "";
  assert.match(panelBodyRule, /overflow:\s*visible/);
  const overviewRule = styles.match(/\.overview-scroll\s*\{([^}]*)\}/)?.[1] ?? "";
  assert.match(overviewRule, /overflow:\s*auto/);
  assert.match(overviewRule, /contain:\s*layout paint/);
});

test("every Panel exposes an accessible full-screen control", () => {
  assert.match(page, /className="panel-expand-button"/);
  assert.match(page, /aria-expanded=\{expanded\}/);
  assert.match(page, /event\.key === "Escape"/);
  assert.match(styles, /\.panel\.panel-expanded\s*\{[^}]*position:\s*fixed/s);
  assert.match(styles, /\.panel-expanded \.panel-body\s*\{[^}]*overflow:\s*auto/s);
});

test("alignment highlighter, reconstruction, and connected tree renderer are wired", () => {
  assert.match(page, /Parent highlighter/);
  assert.match(page, /classifyParentAffinity/);
  assert.match(page, /Global reconstruction/);
  assert.match(page, /buildReconstructionModel/);
  assert.match(page, /className="tree-branches"/);
  assert.match(page, /layoutNeighborJoiningTree/);
});

test("global reconstruction exposes tunable ordered auto-resolution and dependency rescans", () => {
  assert.match(page, /Heuristic auto-resolver/);
  assert.match(page, /conservative/);
  assert.match(page, /balanced/);
  assert.match(page, /aggressive/);
  assert.match(page, /Advanced decision model/);
  assert.match(page, /planAutoResolution/);
  assert.match(page, /applyAutoResolutionPlan/);
  assert.match(page, /rescanTargetsForBarrier/);
  assert.match(page, /filterResolvedEventDuplicates/);
  assert.match(page, /Impacted recombinant targets only/);
  assert.match(page, /Revisit analyst decisions/);
  assert.match(worker, /excludedParents/);
  assert.match(styles, /\.auto-resolve-shell\s*\{[^}]*display:\s*grid/s);
  assert.match(styles, /\.auto-advanced-grid\s*\{[^}]*grid-template-columns:\s*repeat\(3/s);
});

test("RDP-style genome-position matrices replace the event-similarity grid", () => {
  assert.match(page, /Breakpoint-pair density/);
  assert.match(page, /Recombination region separation/);
  assert.match(page, /Local distance-profile discordance/);
  assert.match(page, /RDP4 Figure 2/);
  assert.match(page, /GenomePositionHeatmap/);
  assert.match(page, /useState<48 \| 64 \| 96>\(96\)/);
  assert.doesNotMatch(page, /Events \$\{row \+ 1\} × \$\{column \+ 1\}/);
  assert.match(styles, /\.genome-heatmap-canvas-wrap\s*\{[^}]*overflow:\s*hidden/s);
});

test("server-rendered text does not use ambient locale formatting", () => {
  assert.doesNotMatch(page, /\.toLocaleString\(\)/);
  assert.doesNotMatch(page, /\.toLocaleTimeString/);
  assert.doesNotMatch(page, /new Date\([^)]*\)\.toLocale/);
});
