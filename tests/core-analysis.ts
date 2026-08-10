import assert from "node:assert/strict";
import {
  alignmentStats,
  breakpointHotspotTest,
  demoEvent,
  exportRecombinationFree,
  makeDemoAlignment,
  neighborJoining,
  parseGenomeAnnotations,
  toGff3,
} from "../app/rdp-core";

const names = ["A", "B", "C", "D"];
const distances = [
  0, 0.1, 0.8, 0.8,
  0.1, 0, 0.8, 0.8,
  0.8, 0.8, 0, 0.1,
  0.8, 0.8, 0.1, 0,
];
const tree = neighborJoining(names, distances);
assert.match(tree.newick, /\(A:0\.05000000,B:0\.05000000\)|\(B:0\.05000000,A:0\.05000000\)/);
assert.match(tree.newick, /\(C:0\.05000000,D:0\.05000000\)|\(D:0\.05000000,C:0\.05000000\)/);
assert.ok(tree.newick.endsWith(";"));

const clustered = Array.from({ length: 12 }, (_, index) => ({
  ...demoEvent(),
  id: `clustered-${index}`,
  start: 100 + index,
  end: 500 + index,
}));
const first = breakpointHotspotTest(clustered, 2_400, 48, 499, 1234);
const second = breakpointHotspotTest(clustered, 2_400, 48, 499, 1234);
assert.deepEqual(first, second, "hotspot permutations must be seed-reproducible");
assert.ok(first.observedMaximum >= 12);
assert.ok(first.empiricalP <= 0.01);

const features = parseGenomeAnnotations([
  "##gff-version 3",
  "genome\ttest\tCDS\t101\t2000\t.\t+\t0\tID=cds-1;Name=polyprotein",
].join("\n"), "fixture.gff3", 2_400);
assert.equal(features[0].start, 100);
assert.equal(features[0].end, 2_000);
assert.equal(features[0].name, "polyprotein");
assert.match(toGff3(features), /CDS\t101\t2000/);
const genbankFeatures = parseGenomeAnnotations([
  "FEATURES             Location/Qualifiers",
  "     CDS             complement(300..900)",
  "                     /gene=\"rep\"",
].join("\n"), "fixture.gbk", 2_400);
assert.equal(genbankFeatures[0].strand, "-");
assert.equal(genbankFeatures[0].name, "rep");
const bedFeatures = parseGenomeAnnotations("genome\t10\t40\tleader\t.\t+", "fixture.bed", 2_400);
assert.equal(bedFeatures[0].name, "leader");
const annotated = { ...makeDemoAlignment(), features };
const codonMasked = exportRecombinationFree(annotated, [{ ...demoEvent(), decision: "accepted" }], "mask-codon")[0];
const maskedSites = codonMasked.content.match(/N/g)?.length ?? 0;
assert.ok(maskedSites >= 756);
assert.equal(maskedSites % 3, 0);

const largeSummary = alignmentStats({
  name: "bounded summary fixture",
  format: "generated",
  length: 1_000,
  createdAt: 1,
  sequences: Array.from({ length: 1_000 }, (_, index) => ({
    name: `S${index}`,
    sequence: index % 10 === 0 ? `${"A".repeat(999)}C` : "A".repeat(1_000),
  })),
});
assert.equal(largeSummary.sampled, true);
assert.ok(largeSummary.meanIdentity > 0.99);
