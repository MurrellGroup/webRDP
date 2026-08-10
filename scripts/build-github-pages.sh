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
  || typeof instance.exports.method_stats !== "function"
  || typeof instance.exports.triplet_counts !== "function"
  || typeof instance.exports.hmm_polish !== "function"
) {
  throw new Error("The GitHub Pages artifact is missing the optimized RDP Web kernels.");
}
NODE

echo "GitHub Pages artifact ready in out${pages_base:+ for ${pages_base}}."
