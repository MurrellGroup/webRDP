# RDP Web

RDP Web is a local-first, browser-based recombination analysis workbench. It
combines optimized WebAssembly kernels with an interactive workflow for
alignment/project import, method-specific screening, event verification,
breakpoint editing, matrices/topology checks, and recombination-aware exports.
Version 0.4 adds exact bounded 3SEQ calibration, seeded p-distance
bootstrapping, challenge diagnostics, local NJ trees, breakpoint matrices,
reproducible hotspot tests, manual-event recalculation, undo/redo with an audit
trail, and GFF3/GenBank/BED annotation workflows.

The project is MIT licensed and maintains a strict clean-room boundary from
OpenRDP (GPL-3.0) and the proprietary Windows RDP package. See
`SCIENTIFIC_BASIS.md` for the implementation/validation ledger and
`THIRD_PARTY_NOTICES.md` for licensing details.

## Local development

```sh
npm install
npm run wasm:build
npm run dev
```

The checked-in `public/wasm/rdp.wasm` is built from `assembly/index.ts`. The
module worker in `public/rdp-worker.js` owns the WebAssembly instance and sends
only compact progress/results messages to the interface. `public/rdp-statistics.js`
contains dependency-free probability calibrations.

## GitHub Pages

The repository exports to static files and is safe at both a domain root and a
project subpath. To verify a project-site build locally:

```sh
RDP_BASE_PATH=/your-repository npm run github:build
```

The artifact is written to `out/`. The included
`.github/workflows/github-pages.yml` rebuilds the WebAssembly binary, runs the
scientific regression suite, creates the static export, and deploys it when the
repository's Pages source is set to **GitHub Actions**. The workflow derives the
project path from `GITHUB_REPOSITORY`; repositories named `*.github.io` build at
the domain root.

Worker and WASM URLs are resolved relative to the document/worker locations,
so a repository name or later custom-domain change does not require source
edits.

## Tests and benchmarks

```sh
npm run test:wasm
npm run bench:wasm -- 100 10000
npm run lint
```

The tests cover scalar/packed distance equivalence, exact HGRW tails, seeded
bootstraps, sampled large-dataset parent selection, windowless two-state HMM
breakpoint localization, origin-spanning circular events, rate-variation and
gap-block challenges, exact manual-event recalculation, NJ/hotspot routines, a
513-sequence large-path integration fixture, annotation-aware masking, and
lossless `.rdpweb` project round-trips.

## Current status

This is a functional scientific alpha, not yet a validated drop-in replacement
for RDP5. The application itself surfaces that distinction prominently. See
`ROADMAP.md` for the replacement gates rather than treating the current method
names as a parity claim.
