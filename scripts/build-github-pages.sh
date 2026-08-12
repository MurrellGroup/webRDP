#!/usr/bin/env bash
set -euo pipefail

project_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${project_root}"

if [[ -n "${RDP_BASE_PATH:-}" ]]; then
  pages_base="${RDP_BASE_PATH}"
elif [[ -n "${GITHUB_REPOSITORY:-}" ]]; then
  repository_name="${GITHUB_REPOSITORY##*/}"
  if [[ "${repository_name}" == *.github.io ]]; then
    pages_base=""
  else
    pages_base="/${repository_name}"
  fi
else
  pages_base=""
fi

export RDP_GITHUB_PAGES=1
export RDP_BASE_PATH="${pages_base}"
memory_shim="${project_root}/scripts/node-memory-usage-shim.cjs"
export NODE_OPTIONS="--require=${memory_shim}${NODE_OPTIONS:+ ${NODE_OPTIONS}}"

"${project_root}/node_modules/.bin/next" build --webpack
touch "${project_root}/out/.nojekyll"

[[ -f "${project_root}/out/index.html" ]]
[[ -f "${project_root}/out/rdp-worker.js" ]]
[[ -f "${project_root}/out/rdp-statistics.js" ]]
[[ -f "${project_root}/out/rdp-siscan.js" ]]
[[ -f "${project_root}/out/rdp-phi.js" ]]
[[ -f "${project_root}/out/rdp-bootstrap-tree.js" ]]
[[ -f "${project_root}/out/rdp-burt.js" ]]
[[ -f "${project_root}/out/rdp-clustering.js" ]]
[[ -f "${project_root}/out/rdp-disassembly.js" ]]
[[ -f "${project_root}/out/rdp-recombinant-identification.js" ]]
[[ -f "${project_root}/out/wasm/rdp.wasm" ]]

if [[ -n "${pages_base}" ]]; then
  grep -Fq "${pages_base}/_next/" "${project_root}/out/index.html"
fi

node --input-type=module - "${project_root}/out/wasm/rdp.wasm" <<'NODE'
import { readFile } from "node:fs/promises";
const [wasmPath] = process.argv.slice(2);
const { instance } = await WebAssembly.instantiate(await readFile(wasmPath));
if (
  typeof instance.exports.distance_matrix_packed !== "function"
  || typeof instance.exports.nearest_candidates_sampled !== "function"
  || typeof instance.exports.triplet_counts !== "function"
  || typeof instance.exports.scan_rdp5_triplet_all !== "function"
  || typeof instance.exports.scan_rdp5_triplet_all_packed !== "function"
  || typeof instance.exports.scan_source_chi_all !== "function"
  || typeof instance.exports.scan_source_chi_all_packed !== "function"
  || typeof instance.exports.scan_source_three_seq_triplet !== "function"
  || typeof instance.exports.scan_source_three_seq_triplet_packed !== "function"
  || typeof instance.exports.source_three_seq_workspace_bytes !== "function"
  || typeof instance.exports.scan_source_three_seq_triplet_mode !== "function"
  || typeof instance.exports.scan_source_three_seq_triplet_packed_mode !== "function"
) {
  throw new Error("The GitHub Pages artifact is missing the optimized RDP Web kernels.");
}
NODE

echo "GitHub Pages artifact ready in out${pages_base:+ for ${pages_base}}."
