import assert from "node:assert/strict";
import { readFile } from "node:fs/promises";
import test from "node:test";

const [page, styles] = await Promise.all([
  readFile(new URL("../app/page.tsx", import.meta.url), "utf8"),
  readFile(new URL("../app/globals.css", import.meta.url), "utf8"),
]);

test("workspace owns vertical scrolling and generic panels expand in document flow", () => {
  assert.match(styles, /\.workspace-content\s*\{[^}]*overflow-y:\s*auto/s);
  assert.match(styles, /\.workspace-content\s*\{[^}]*height:\s*0/s);
  const panelRule = styles.match(/\.panel\s*\{([^}]*)\}/)?.[1] ?? "";
  assert.doesNotMatch(panelRule, /max-height/);
  assert.match(panelRule, /overflow:\s*visible/);
  const panelBodyRule = styles.match(/\.panel-body\s*\{([^}]*)\}/)?.[1] ?? "";
  assert.match(panelBodyRule, /overflow:\s*visible/);
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

test("server-rendered text does not use ambient locale formatting", () => {
  assert.doesNotMatch(page, /\.toLocaleString\(\)/);
  assert.doesNotMatch(page, /\.toLocaleTimeString/);
  assert.doesNotMatch(page, /new Date\([^)]*\)\.toLocale/);
});
