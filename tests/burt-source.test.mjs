import assert from "node:assert/strict";
import test from "node:test";
import { collectBurtObservations, fitBurtTriplet, fitCategoricalHMM } from "../public/rdp-burt.js";

function categoricalTriplet(categories) {
  const length = categories.length;
  const encoded = new Uint8Array(3 * length);
  for (let site = 0; site < length; site += 1) {
    const bases = categories[site] === 0
      ? [0, 0, 1]
      : categories[site] === 1
        ? [1, 0, 0]
        : [0, 1, 0];
    for (let sequence = 0; sequence < 3; sequence += 1) encoded[sequence * length + site] = bases[sequence];
  }
  return encoded;
}

test("BURT source categories follow the supplied BenHMM recoding", () => {
  const encoded = categoricalTriplet([0, 1, 2]);
  const result = collectBurtObservations(encoded, 3, 0, 1, 2);
  assert.deepEqual([...result.observations], [0, 1, 2]);
  assert.deepEqual([...result.positions], [0, 1, 2]);
});

test("RDP5 source-mode BURT deterministically recovers a two-switch tract", () => {
  const categories = Uint8Array.from({ length: 900 }, (_, site) => site < 300 || site >= 600 ? 0 : 2);
  const encoded = categoricalTriplet(categories);
  const options = { sourceParity: true, randomStarts: 21, maxIterations: 100, posteriorThreshold: 0.995, seed: 1234 };
  const first = fitBurtTriplet(encoded, 900, 0, 1, 2, 280, 620, options);
  const second = fitBurtTriplet(encoded, 900, 0, 1, 2, 280, 620, options);
  assert.ok(first);
  assert.equal(first.model.states, 3);
  assert.equal(first.model.sourceCompatibility, "RDP5 BenHMM + DoHMMCyclesSerial");
  assert.equal(first.model.randomStarts, 21);
  assert.equal(first.model.stateSwitches, 2);
  assert.ok(Math.abs(first.start - 300) <= 1);
  assert.ok(Math.abs(first.end - 600) <= 1);
  assert.ok(first.confidenceStart[0] <= 300 && first.confidenceStart[1] >= 300);
  assert.ok(first.confidenceEnd[0] <= 600 && first.confidenceEnd[1] >= 600);
  assert.deepEqual(first, second, "seeded source mode must be byte-for-byte deterministic");
});

test("manual 2–20-state mode selects three categorical regimes by BIC", () => {
  const observations = Uint8Array.from({ length: 1_200 }, (_, site) => site < 400 ? 0 : site < 800 ? 1 : 2);
  const fitted = fitCategoricalHMM(observations, {
    sourceParity: false,
    randomStarts: 10,
    maxIterations: 100,
    maxStates: 6,
    exhaustiveModels: true,
    seed: 19,
  });
  assert.equal(fitted.selected.stateCount, 3);
  assert.equal(fitted.ledger.length, 5);
  assert.ok(fitted.ledger.find((entry) => entry.states === 3).bic < fitted.ledger.find((entry) => entry.states === 2).bic);
});
