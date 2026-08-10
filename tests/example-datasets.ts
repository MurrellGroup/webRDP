import assert from "node:assert/strict";
import { EXAMPLE_DATASETS } from "../app/example-datasets";

assert.equal(EXAMPLE_DATASETS.length, 8);
assert.ok(EXAMPLE_DATASETS.some((example) => example.organism.toLowerCase().includes("bacterial")));
assert.ok(EXAMPLE_DATASETS.some((example) => example.sequenceCount >= 500));
assert.ok(EXAMPLE_DATASETS.some((example) => example.recommendedOptions.circular));

for (const example of EXAMPLE_DATASETS) {
  const alignment = example.generate();
  assert.equal(alignment.sequences.length, example.sequenceCount, `${example.id}: sequence count`);
  assert.equal(alignment.length, example.length, `${example.id}: declared alignment length`);
  assert.ok(alignment.sequences.every((record) => record.sequence.length === example.length), `${example.id}: aligned records`);
  assert.ok(example.truth.length > 0, `${example.id}: truth annotations`);
}
