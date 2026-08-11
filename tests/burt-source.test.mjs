import assert from "node:assert/strict";
import test from "node:test";
import {
  buildSourceBurtWorkingSet,
  collectBurtObservations,
  fitBurtTriplet,
  fitCategoricalHMM,
  matchSourceBreakpoint,
  polishSourceBreakpointPair,
  vbLong,
} from "../public/rdp-burt.js";

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

test("BURT preserves VB6 midpoint rounding and the source circular sentinel layout", () => {
  assert.equal(vbLong(2.5), 2);
  assert.equal(vbLong(3.5), 4);
  const positions = Int32Array.from({ length: 10 }, (_, index) => index);
  const observations = Uint8Array.from([1, 1, 1, 1, 1, 2, 2, 2, 2, 2]);
  const working = buildSourceBurtWorkingSet(positions, observations, 10, true);
  assert.equal(working.circularOffset, 5);
  assert.equal(working.sourceLength, 20);
  assert.deepEqual([...working.observations.slice(0, 6)], [2, 2, 2, 2, 2, 0], "RDP5 leaves a zero sentinel after the rotated leading half-copy");
  assert.deepEqual([...working.observations.slice(6, 16)], [...observations]);
  assert.deepEqual([...working.observations.slice(16, 21)], [1, 1, 1, 1, 1]);
  assert.deepEqual([...working.xDiffPos.slice(1, 11)], [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]);
});

test("MatchBPtoCI gives an in-99%-CI switch priority and keeps strict source ties", () => {
  const xPosDiff = Int32Array.from({ length: 101 }, (_, index) => index);
  const switches = [
    { position: 20, confidence95: [15, 25], confidence99: [10, 30], sourceCoordinates: [10, 30, 20, 15, 25] },
    { position: 26, confidence95: [24, 28], confidence99: [40, 50], sourceCoordinates: [40, 50, 26, 24, 28] },
  ];
  const inside = matchSourceBreakpoint(switches, 25, xPosDiff, 100, false);
  assert.equal(inside.switchIndex, 0, "a switch whose 99% CI contains the proposal beats a closer outside switch");
  assert.equal(inside.within99, true);
  const outside = matchSourceBreakpoint(switches, 35, xPosDiff, 100, false);
  assert.equal(outside.switchIndex, 1);
  assert.equal(outside.within99, false);
  assert.ok(outside.signedCoordinates.every((value) => value <= 0), "outside matches carry the desktop negative-CI signal");
  const tied = matchSourceBreakpoint([
    { position: 45, confidence95: [40, 50], confidence99: [35, 65], sourceCoordinates: [35, 65, 45, 40, 50] },
    { position: 55, confidence95: [50, 60], confidence99: [35, 65], sourceCoordinates: [35, 65, 55, 50, 60] },
  ], 50, xPosDiff, 100, false);
  assert.equal(tied.switchIndex, 0, "strict comparisons preserve the first switch on an exact VNP-distance tie");
  const wrapped = matchSourceBreakpoint([
    { position: 98, confidence95: [94, 5], confidence99: [90, 10], sourceCoordinates: [90, 10, 98, 94, 5] },
  ], 5, xPosDiff, 100, true);
  assert.equal(wrapped.within99, true);
});

test("RDP5 source-mode BURT deterministically recovers a two-switch tract", () => {
  const categories = Uint8Array.from({ length: 900 }, (_, site) => site < 300 || site >= 600 ? 0 : 2);
  const encoded = categoricalTriplet(categories);
  const options = { sourceParity: true, randomStarts: 21, maxIterations: 100, posteriorThreshold: 0.995, seed: 1234 };
  const first = fitBurtTriplet(encoded, 900, 0, 1, 2, 280, 620, options);
  const second = fitBurtTriplet(encoded, 900, 0, 1, 2, 280, 620, options);
  assert.ok(first);
  assert.equal(first.model.states, 3);
  assert.equal(first.model.sourceCompatibility, "RDP5 BenHMM + DoHMMCyclesSerial + MatchBPtoCI + PolishBP");
  assert.equal(first.model.randomStarts, 21);
  assert.equal(first.model.stateSwitches, 2);
  assert.ok(Math.abs(first.start - 300) <= 1);
  assert.ok(Math.abs(first.end - 600) <= 1);
  assert.ok(first.confidenceStart[0] <= 300 && first.confidenceStart[1] >= 300);
  assert.ok(first.confidenceEnd[0] <= 600 && first.confidenceEnd[1] >= 600);
  assert.equal(first.model.polishDecision.startAdopted, true);
  assert.equal(first.model.polishDecision.endAdopted, true);
  assert.equal(first.model.polishDecision.startWithin99, false, "PolishBP may adopt an outside-CI switch within half a tract length");
  assert.deepEqual(first.model.candidateBreakpoints, [280, 620]);
  assert.deepEqual(first.model.polishedBreakpoints, [300, 600]);
  assert.ok(first.model.posteriorTrace.every((entry) => entry.probabilities.length === 3));
  assert.deepEqual(first, second, "seeded source mode must be byte-for-byte deterministic");
});

test("source BURT circular padding recovers an origin-spanning event", () => {
  const length = 900;
  const categories = Uint8Array.from({ length }, (_, site) => site < 150 || site >= 750 ? 2 : 0);
  const result = fitBurtTriplet(categoricalTriplet(categories), length, 0, 1, 2, 750, 150, {
    sourceParity: true,
    circular: true,
    randomStarts: 21,
    maxIterations: 100,
    seed: 1234,
  });
  assert.ok(result);
  assert.equal(result.start, 750);
  assert.equal(result.end, 150);
  assert.deepEqual(result.model.switches.map((entry) => entry.position), [150, 750]);
  assert.deepEqual(result.model.circularPadding, { offset: 450, fittedSites: 1801, croppedSites: 901 });
  assert.equal(result.model.polishDecision.startWithin99, true);
  assert.equal(result.model.polishDecision.endWithin99, true);
});

test("PolishBP resolves two candidates mapped to one switch without moving both", () => {
  const categories = Uint8Array.from({ length: 900 }, (_, site) => site < 300 || site >= 600 ? 0 : 2);
  const result = fitBurtTriplet(categoricalTriplet(categories), 900, 0, 1, 2, 282, 325, {
    sourceParity: true,
    randomStarts: 21,
    maxIterations: 100,
    seed: 1234,
  });
  assert.ok(result);
  assert.equal(result.model.polishDecision.sameSwitchResolved, true);
  assert.equal(result.start, 300);
  assert.equal(result.end, 325);
  assert.equal(result.model.polishDecision.startAdopted, true);
  assert.equal(result.model.polishDecision.endAdopted, false);
});

test("PolishBP snaps breakpoints to a missing-data edge", () => {
  const categories = Uint8Array.from({ length: 900 }, (_, site) => site < 300 || site >= 600 ? 0 : 2);
  const encoded = categoricalTriplet(categories);
  for (let site = 295; site <= 299; site += 1) {
    for (let sequence = 0; sequence < 3; sequence += 1) encoded[sequence * 900 + site] = 4;
  }
  const result = fitBurtTriplet(encoded, 900, 0, 1, 2, 280, 620, {
    sourceParity: true,
    randomStarts: 21,
    maxIterations: 100,
    seed: 1234,
  });
  assert.ok(result);
  assert.equal(result.start, 301);
  assert.equal(result.model.polishDecision.startMissingBoundary, true);
});

test("PolishBP restores the proposal when a polished partition has fewer than three variable sites on one side", () => {
  const encoded = categoricalTriplet(Uint8Array.from({ length: 100 }, (_, site) => site % 2 ? 0 : 2));
  const xPosDiff = Int32Array.from({ length: 101 }, (_, index) => index);
  const switches = [
    { position: 2, confidence95: [1, 5], confidence99: [1, 5], sourceCoordinates: [1, 5, 2, 1, 5] },
    { position: 99, confidence95: [95, 100], confidence99: [95, 100], sourceCoordinates: [95, 100, 99, 95, 100] },
  ];
  const result = polishSourceBreakpointPair(encoded, 100, [0, 1, 2], 3, 98, switches, xPosDiff, false);
  assert.equal(result.revertedForInformation, true);
  assert.equal(result.information.outside, 2);
  assert.deepEqual([result.start, result.end], [3, 98]);
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
