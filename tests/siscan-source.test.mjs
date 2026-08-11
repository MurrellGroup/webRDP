import assert from "node:assert/strict";
import test from "node:test";
import {
  buildSourceSiScanRandomization,
  runSourceSiScan,
  selectSourceSiScanOutgroup,
  sourceNormalZ,
  sourceSiScanPattern,
} from "../public/rdp-siscan.js";

function exactMosaic(length = 240, start = 80, end = 160) {
  const encoded = new Uint8Array(4 * length);
  for (let site = 0; site < length; site += 1) {
    encoded[site] = site >= start && site < end ? 1 : 0; // recombinant
    encoded[length + site] = 0; // major parent
    encoded[2 * length + site] = 1; // minor parent
    encoded[3 * length + site] = 2; // fourth sequence
  }
  return encoded;
}

test("SiScan uses the supplied source's 15-category quartet table", () => {
  assert.equal(sourceSiScanPattern(0, 0, 1, 0), 12);
  assert.equal(sourceSiScanPattern(0, 0, 1, 1), 8);
  assert.equal(sourceSiScanPattern(0, 1, 2, 3), 1);
  assert.equal(sourceSiScanPattern(0, 1, 2, 4), 0, "ambiguous outgroups are stripped by default");
  assert.equal(sourceSiScanPattern(5, 5, 0, 5, { gapsAsState: true }), 12);
});

test("GetSSOL direct fallback chooses the closest eligible fourth sequence", () => {
  const length = 20;
  const encoded = new Uint8Array(5 * length);
  encoded.fill(0);
  encoded.fill(1, 2 * length, 3 * length);
  for (let site = 0; site < length; site += 1) {
    encoded[3 * length + site] = site < 2 ? 1 : 0;
    encoded[4 * length + site] = site < 10 ? 1 : 0;
  }
  const selected = selectSourceSiScanOutgroup(encoded, length, 5, [0, 1, 2], {
    outgroupMode: "nearest",
    candidatePool: [3, 4],
  });
  assert.equal(selected.index, 3);
  assert.match(selected.sourcePath, /direct-distance fallback/);
});

test("RDP5 Sister-Scanning recovers a known topology tract and every calibrated run", () => {
  const encoded = exactMosaic();
  const options = {
    window: 40,
    step: 10,
    scanPermutations: 100,
    pValuePermutations: 1000,
    seed: 91,
    outgroupMode: "nearest",
    positionMode: "triplet-variable",
    candidatePool: [3],
  };
  const result = runSourceSiScan(encoded, 240, 4, [0, 1, 2], options);
  assert.ok(result);
  assert.equal(result.start, 80);
  assert.equal(result.end, 160);
  assert.equal(result.outgroupIndex, 3);
  assert.equal(result.baselineTopology, 0);
  assert.equal(result.inferredTopology, 1);
  assert.ok(result.z > 10);
  assert.ok(result.rawP < 1e-10);
  assert.ok(result.regions.some((region) => region.start === 80 && region.end === 160));
  assert.match(result.sourceRoutine, /GetSSOL.*DoPerms3P.*ShrinkRegionC/);
});

test("streamed large-genome randomization is identical to the cached source table", () => {
  const encoded = exactMosaic();
  const cached = buildSourceSiScanRandomization(240, 1000, 19);
  const streamed = buildSourceSiScanRandomization(10_000, 1000, 19);
  assert.ok(cached.values);
  assert.equal(streamed.values, null);
  const common = {
    window: 40,
    step: 10,
    scanPermutations: 100,
    pValuePermutations: 1000,
    seed: 19,
    candidatePool: [3],
  };
  const left = runSourceSiScan(encoded, 240, 4, [0, 1, 2], { ...common, randomization: cached });
  const right = runSourceSiScan(encoded, 240, 4, [0, 1, 2], { ...common, randomization: streamed });
  const reference = runSourceSiScan(encoded, 240, 4, [0, 1, 2], {
    ...common,
    randomization: cached,
    referencePermutationPath: true,
  });
  assert.equal(left.rawP, right.rawP);
  assert.equal(left.z, right.z);
  assert.equal(left.pattern, right.pattern);
  assert.deepEqual(left.windows, reference.windows, "prefix-range optimization must exactly match direct DoPerms3P enumeration");
  assert.equal(left.rawP, reference.rawP);
});

test("NormalZ retains the desktop two-sided tail convention", () => {
  assert.ok(Math.abs(sourceNormalZ(1.96) - 0.0499958) < 1e-6);
  assert.equal(sourceNormalZ(0), 1);
  assert.ok(sourceNormalZ(8) < 1e-12);
});
