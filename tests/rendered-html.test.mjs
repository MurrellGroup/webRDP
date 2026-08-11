import assert from "node:assert/strict";
import fs from "node:fs";
import test from "node:test";

const developmentPreviewMeta =
  /<meta(?=[^>]*\bname=["']codex-preview["'])(?=[^>]*\bcontent=["']development["'])[^>]*>/i;

test("renders development preview metadata", async () => {
  const workerUrl = new URL("../dist/server/index.js", import.meta.url);
  workerUrl.searchParams.set("test", `${process.pid}-${Date.now()}`);
  const { default: worker } = await import(workerUrl.href);

  const response = await worker.fetch(
    new Request("http://localhost/", {
      headers: { accept: "text/html" },
    }),
    {
      ASSETS: {
        fetch: async () => new Response("Not found", { status: 404 }),
      },
    },
    {
      waitUntil() {},
      passThroughOnException() {},
    },
  );

  assert.equal(response.status, 200);
  assert.match(
    response.headers.get("content-type") ?? "",
    /^text\/html\b/i,
  );
  const html = await response.text();
  assert.match(html, developmentPreviewMeta);
  assert.match(html, /Events &amp; evidence/);
  assert.match(html, /Local trees/);
  assert.match(html, /Method-by-method result/);
  assert.match(html, /Known synthetic truth for orientation only/);
  assert.match(
    html,
    /<div(?=[^>]*class="panel-body")(?=[^>]*tabindex="0")(?=[^>]*role="region")(?=[^>]*aria-label="Method-by-method result content")[^>]*>/,
  );
});

test("every constrained application surface has a reachable scroll owner", () => {
  const css = fs.readFileSync(new URL("../app/globals.css", import.meta.url), "utf8");
  assert.match(css, /\.app-shell\s*\{[^}]*height:\s*100dvh[^}]*min-height:\s*0[^}]*overflow:\s*hidden/s);
  assert.match(css, /\.workbench\s*\{[^}]*min-height:\s*0[^}]*overflow:\s*hidden/s);
  assert.match(css, /\.workspace-content\s*\{[^}]*flex:\s*1 1 0[^}]*height:\s*0[^}]*min-height:\s*0[^}]*overflow:\s*auto/s);
  assert.match(css, /\.panel\s*\{[^}]*max-height:\s*calc\(100dvh - 150px\)[^}]*min-height:\s*0[^}]*display:\s*flex[^}]*flex-direction:\s*column[^}]*overflow:\s*hidden/s);
  assert.match(css, /\.panel-body\s*\{[^}]*flex:\s*1 1 auto[^}]*min-height:\s*0[^}]*overflow:\s*auto[^}]*overscroll-behavior:\s*auto/s);
  assert.match(css, /\.method-result-tabs\s*\{[^}]*position:\s*sticky[^}]*overflow-y:\s*auto/s);
  assert.match(css, /\.method-profile\s*\{[^}]*grid-column:\s*2/s);
  assert.match(css, /\.sidebar-scroll\s*\{[^}]*min-height:\s*0[^}]*overflow-y:\s*auto/s);
  assert.match(css, /\.inspector\s*\{[^}]*overflow-y:\s*auto[^}]*min-height:\s*0/s);
  assert.match(css, /\.modal\s*\{[^}]*max-height:\s*calc\(100dvh - 40px\)[^}]*overflow:\s*auto/s);
  assert.match(css, /\.tutorial-modal\s*\{[^}]*overflow:\s*auto/s);
  assert.doesNotMatch(css, /\.tutorial-modal\s*\{[^}]*overflow:\s*hidden/s);
  assert.match(css, /@media \(max-width:\s*980px\)[\s\S]*?body\s*\{\s*overflow:\s*auto/s);
  assert.match(css, /@media \(max-width:\s*980px\)[\s\S]*?\.workspace-content\s*\{[^}]*height:\s*auto[^}]*overflow:\s*visible/s);
  assert.match(css, /@media \(max-width:\s*980px\)[\s\S]*?\.panel\s*\{[^}]*max-height:\s*none/s);
  assert.match(css, /@media \(max-width:\s*980px\)[\s\S]*?\.panel-body\s*\{[^}]*overflow:\s*visible/s);
});
