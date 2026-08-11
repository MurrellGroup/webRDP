import assert from "node:assert/strict";
import { readFile } from "node:fs/promises";
import test from "node:test";

const [page, styles, worker, recombinantIdentification, sisterScan, phi, burt] = await Promise.all([
  readFile(new URL("../app/page.tsx", import.meta.url), "utf8"),
  readFile(new URL("../app/globals.css", import.meta.url), "utf8"),
  readFile(new URL("../public/rdp-worker.js", import.meta.url), "utf8"),
  readFile(new URL("../public/rdp-recombinant-identification.js", import.meta.url), "utf8"),
  readFile(new URL("../public/rdp-siscan.js", import.meta.url), "utf8"),
  readFile(new URL("../public/rdp-phi.js", import.meta.url), "utf8"),
  readFile(new URL("../public/rdp-burt.js", import.meta.url), "utf8"),
]);

test("the false-positive studio uses the source PHI statistic rather than a PHI-labelled surrogate", () => {
  assert.match(page, /RDP5 PHI analytic p/);
  assert.match(page, /source PHITest2 moments/);
  assert.match(worker, /sourcePhiTest/);
  assert.match(phi, /export function sourcePairIncompatibility/);
  assert.match(phi, /export function sourcePhiAnalyticMeanVariance/);
  assert.match(phi, /RDP5 PHITest2\/PHI\/pair_score/);
});

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

test("source SiScan controls, provenance, and interactive trace are wired", () => {
  assert.match(page, /SiScan source controls/);
  assert.match(page, /Nearest outlier · RDP5 default/);
  assert.match(page, /One analyst-selected sequence/);
  assert.match(page, /Final-region permutations/);
  assert.match(page, /RDP5 source topology-run Z/);
  assert.match(page, /selectedSignal\?\.profile/);
  assert.match(worker, /runSourceSiScan/);
  assert.match(worker, /sourceSiScanProfile/);
  assert.match(worker, /sourceResult\?\.regions/);
  assert.match(sisterScan, /export function sourceSiScanPattern/);
  assert.match(sisterScan, /buildPermutationPrefix/);
  assert.match(sisterScan, /GetSSOL \+ Get3Score\/GetPScores2 \+ DoPerms3P/);
});

test("source BURT exposes the desktop polish path and an interactive posterior workbench", () => {
  assert.match(burt, /export function buildSourceBurtWorkingSet/);
  assert.match(burt, /export function sourceBurtSwitches/);
  assert.match(burt, /export function matchSourceBreakpoint/);
  assert.match(burt, /export function polishSourceBreakpointPair/);
  assert.match(burt, /DoHMMCyclesSerial.*GetLaticePathP.*ForwardCP.*ReverseCP.*MatchBPtoCI.*PolishBP/s);
  assert.match(page, /BURT posterior evidence/);
  assert.match(page, /Use as start/);
  assert.match(page, /Use as end/);
  assert.match(page, /Source circular padding applied/);
  assert.match(page, /0\.995 \/ 0\.999 fixed/);
  assert.match(styles, /\.burt-evidence-plot/);
  assert.match(styles, /\.burt-switch-actions\s*\{[^}]*overflow:\s*auto/s);
});

test("full scans default to every concrete sequence triplet with approximate pruning clearly isolated", () => {
  assert.match(page, /Every unordered set of three actual sequences is screened/);
  assert.match(page, /All concrete triplets/);
  assert.match(page, /Use approximate parent shortlist/);
  assert.match(page, /not RDP5 triplet parity/);
  assert.match(page, /three explicit sequences each/);
  assert.match(worker, /all-concrete-triplets/);
  assert.match(worker, /No alignment consensus or rest-of-alignment proxy/);
  assert.match(worker, /concreteTripletInputs:\s*true/);
  assert.match(styles, /\.triplet-coverage-setting/);
  assert.match(styles, /\.approximate-setting/);
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

test("source-guided review studio exposes the ordered RDP5 refinement workflow", () => {
  assert.match(page, /RDP5-style hypothesis refinement/);
  assert.match(page, /Ordered reconstruction queue/);
  assert.match(page, /Best unresolved/);
  assert.match(page, /Skip accepted/);
  assert.match(page, /Method-by-method confirmation/);
  assert.match(page, /Role-assignment challenge/);
  assert.match(page, /RDP5 recombinant identification/);
  assert.match(page, /Source test ledger/);
  assert.match(page, /Role-consensus confidence/);
  assert.match(page, /Apply these breakpoints to all/);
  assert.match(page, /Accept group/);
  assert.match(page, /Tract vs combined background/);
  assert.match(page, /midpoint-oriented for readability but has no inferred root/);
  assert.match(page, /navigateReviewEvent/);
  assert.match(page, /roleAssignmentTrials/);
  assert.match(page, /GENECONV G-scale/);
  assert.match(page, /RDP signals retained\/triplet/);
  assert.match(page, /Tree bootstrap/);
  assert.match(page, /Collapse branches below/);
  assert.match(page, /Bootstrap site blocks/);
  assert.match(page, /Tree cohort cap/);
  assert.match(page, /Signal-disassembly lineage/);
  assert.match(styles, /\.review-queue-list\s*\{[^}]*overflow:\s*auto/s);
  assert.match(styles, /\.review-workspace\s*\{[^}]*display:\s*grid/s);
  assert.match(styles, /\.role-source-tests\s*\{[^}]*display:\s*grid/s);
  assert.match(worker, /identifyRecombinantRoles/);
  assert.match(recombinantIdentification, /export function sourcePhPrScores/);
  assert.match(recombinantIdentification, /export function sourceTrpScores/);
  assert.match(recombinantIdentification, /export function sourceOuCheckScores/);
  assert.match(recombinantIdentification, /export function sourceHistoricalSetMembers/);
  assert.match(recombinantIdentification, /export function sourceParsimonyScores/);
  assert.match(recombinantIdentification, /export function sourceConflictScores/);
  assert.match(recombinantIdentification, /export function sourceSetDistanceScores/);
  assert.match(recombinantIdentification, /export function sourceDmaxScores/);
  assert.match(recombinantIdentification, /MakeTrpGroups \+ MakeTrpScore/);
  assert.match(recombinantIdentification, /FinalTrim \+ RCompatC\/RCompatD \+ MakeConsensusC/);
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
