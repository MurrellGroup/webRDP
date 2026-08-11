// BURT implementation for RDP Web.
//
// The default compatibility path is translated from the author-supplied
// RDP5 BenHMM/DoHMMCyclesSerial source with permission: a three-state model,
// 21 starts, Viterbi training, and 0.995/0.999 posterior intervals. The manual
// 2..20-state step-up interpretation remains available as an explicit mode.

const LOG_FLOOR = Math.log(1e-300);

function clampProbability(value) {
  return Math.max(1e-12, Math.min(1 - 1e-12, value));
}

function makeRng(seed) {
  let state = (seed >>> 0) || 0x6d2b79f5;
  return () => {
    state ^= state << 13;
    state ^= state >>> 17;
    state ^= state << 5;
    return (state >>> 0) / 4294967296;
  };
}

// The Windows build's C runtime rand() is the MSVC 15-bit LCG. Matching it
// makes random-start selection reproducible against source-derived fixtures.
function makeMsvcRng(seed) {
  let state = seed >>> 0;
  return () => {
    state = (Math.imul(state, 214013) + 2531011) >>> 0;
    return ((state >>> 16) & 0x7fff) / 0x7fff;
  };
}

function logSumExp(values, count = values.length) {
  let maximum = -Infinity;
  for (let index = 0; index < count; index += 1) maximum = Math.max(maximum, values[index]);
  if (!Number.isFinite(maximum)) return LOG_FLOOR;
  let total = 0;
  for (let index = 0; index < count; index += 1) total += Math.exp(values[index] - maximum);
  return maximum + Math.log(Math.max(1e-300, total));
}

function normalise(values, offset, count) {
  let total = 0;
  for (let index = 0; index < count; index += 1) total += values[offset + index];
  if (!(total > 0)) {
    for (let index = 0; index < count; index += 1) values[offset + index] = 1 / count;
    return;
  }
  for (let index = 0; index < count; index += 1) values[offset + index] /= total;
}

export function collectBurtObservations(encoded, nSites, sequence1, sequence2, sequence3) {
  const positions = [];
  const observations = [];
  const firstOffset = sequence1 * nSites;
  const secondOffset = sequence2 * nSites;
  const thirdOffset = sequence3 * nSites;
  for (let site = 0; site < nSites; site += 1) {
    const first = encoded[firstOffset + site];
    const second = encoded[secondOffset + site];
    const third = encoded[thirdOffset + site];
    if (first >= 4 || second >= 4 || third >= 4) continue;
    if (first === second && second !== third) {
      positions.push(site);
      observations.push(0); // A: sequence 1 = sequence 2
    } else if (first === third && first !== second) {
      positions.push(site);
      observations.push(2); // source category C: sequence 1 = sequence 3
    } else if (second === third && second !== first) {
      positions.push(site);
      observations.push(1); // source category B: sequence 2 = sequence 3
    }
  }
  return {
    positions: Int32Array.from(positions),
    observations: Uint8Array.from(observations),
  };
}

function initialiseModel(stateCount, random, options) {
  const initial = new Float64Array(stateCount);
  const transition = new Float64Array(stateCount * stateCount);
  const emission = new Float64Array(stateCount * 3);
  if (options.sourceParity === true) {
    const switchProbability = Math.max(1e-9, Math.min(0.999999, 5 / Math.max(6, options.alignmentLength ?? 1)));
    const imbalance = (Math.floor(random() * 3) + 1) / 10;
    const anchors = [];
    const used = new Uint8Array(3);
    // DoHMMCyclesSerial does not shuffle [0,1,2]. It repeatedly draws
    // floor(6*rand)-2 until an in-range, unused emission anchor appears. The
    // rejected draws materially change every later random start, so preserve
    // this odd-looking source loop exactly.
    while (anchors.length < stateCount) {
      const candidate = Math.floor(6 * random()) - 2;
      if (candidate < 0 || candidate >= 3 || used[candidate]) continue;
      used[candidate] = 1;
      anchors.push(candidate);
    }
    for (let state = 0; state < stateCount; state += 1) {
      initial[state] = 1 / stateCount;
      for (let next = 0; next < stateCount; next += 1) {
        transition[state * stateCount + next] = state === next
          ? 1 - switchProbability
          : switchProbability / Math.max(1, stateCount - 1);
      }
      for (let category = 0; category < 3; category += 1) {
        emission[state * 3 + category] = category === anchors[state % anchors.length]
          ? 1 / 3 + imbalance * 2
          : 1 / 3 - imbalance;
      }
    }
    return { initial, transition, emission };
  }
  for (let state = 0; state < stateCount; state += 1) {
    initial[state] = 0.25 + random();
    for (let next = 0; next < stateCount; next += 1) {
      // Long ancestry tracts are expected, while the random component lets
      // distinct starts explore non-symmetric transition matrices.
      transition[state * stateCount + next] = (state === next ? stateCount * 5 : 0.08) + random();
    }
    const anchor = state % 3;
    for (let category = 0; category < 3; category += 1) {
      emission[state * 3 + category] = (category === anchor ? 2.5 : 0.25) + random() * 1.5;
    }
    normalise(transition, state * stateCount, stateCount);
    normalise(emission, state * 3, 3);
  }
  normalise(initial, 0, stateCount);
  return { initial, transition, emission };
}

function viterbi(observations, model, sourceBacktrace = false) {
  const { initial, transition, emission } = model;
  const stateCount = initial.length;
  const siteCount = observations.length;
  const predecessors = new Uint8Array(siteCount * stateCount);
  let previous = new Float64Array(stateCount);
  let current = new Float64Array(stateCount);
  for (let state = 0; state < stateCount; state += 1) {
    previous[state] = Math.log(clampProbability(initial[state]))
      + Math.log(clampProbability(emission[state * 3 + observations[0]]));
  }
  for (let site = 1; site < siteCount; site += 1) {
    const observation = observations[site];
    for (let state = 0; state < stateCount; state += 1) {
      let bestState = 0;
      let best = -Infinity;
      for (let prior = 0; prior < stateCount; prior += 1) {
        const score = previous[prior] + Math.log(clampProbability(transition[prior * stateCount + state]));
        if (score > best) {
          best = score;
          bestState = prior;
        }
      }
      current[state] = best + Math.log(clampProbability(emission[state * 3 + observation]));
      predecessors[site * stateCount + state] = bestState;
    }
    [previous, current] = [current, previous];
  }
  let finalState = 0;
  let score = previous[0];
  for (let state = 1; state < stateCount; state += 1) {
    if (previous[state] > score) {
      score = previous[state];
      finalState = state;
    }
  }
  const path = new Uint8Array(siteCount);
  if (sourceBacktrace && siteCount > 1) {
    // GetLaticePathP stores the predecessor of the best terminal state at the
    // terminal index. BenHMM's one-based XDiffPos mapping compensates for this
    // historical offset, so the compatibility path preserves it.
    path[siteCount - 1] = predecessors[(siteCount - 1) * stateCount + finalState];
    for (let site = siteCount - 2; site >= 0; site -= 1) {
      path[site] = site === 0 ? 0 : predecessors[site * stateCount + path[site + 1]];
    }
  } else {
    path[siteCount - 1] = finalState;
    for (let site = siteCount - 1; site > 0; site -= 1) {
      path[site - 1] = predecessors[site * stateCount + path[site]];
    }
  }
  return { path, score };
}

function updateModel(observations, path, stateCount, options) {
  const initial = new Float64Array(stateCount);
  const transition = new Float64Array(stateCount * stateCount);
  const emission = new Float64Array(stateCount * 3);
  const occupancy = new Uint32Array(stateCount);
  // Weak symmetric pseudocounts keep every fitted model proper. A modest
  // diagonal prior prevents single-site states without forcing a switch rate.
  const sourceParity = options.sourceParity === true;
  const pseudocount = sourceParity ? 0.01 : 0.25;
  initial.fill(sourceParity ? 1 / stateCount : 0.25);
  transition.fill(sourceParity ? 0.01 : 0.05);
  emission.fill(pseudocount);
  if (!sourceParity) {
    initial[path[0]] += 1;
    for (let state = 0; state < stateCount; state += 1) transition[state * stateCount + state] += 0.5;
  }
  for (let site = 0; site < observations.length; site += 1) {
    const state = path[site];
    occupancy[state] += 1;
    emission[state * 3 + observations[site]] += 1;
    if (site > 0) transition[path[site - 1] * stateCount + state] += 1;
  }
  normalise(initial, 0, stateCount);
  for (let state = 0; state < stateCount; state += 1) {
    normalise(transition, state * stateCount, stateCount);
    normalise(emission, state * 3, 3);
  }
  return { initial, transition, emission, occupancy };
}

function samePath(left, right) {
  if (!left || left.length !== right.length) return false;
  for (let index = 0; index < left.length; index += 1) if (left[index] !== right[index]) return false;
  return true;
}

function forwardLogLikelihood(observations, model) {
  const { initial, transition, emission } = model;
  const stateCount = initial.length;
  let previous = new Float64Array(stateCount);
  let current = new Float64Array(stateCount);
  const terms = new Float64Array(stateCount);
  for (let state = 0; state < stateCount; state += 1) {
    previous[state] = Math.log(clampProbability(initial[state]))
      + Math.log(clampProbability(emission[state * 3 + observations[0]]));
  }
  for (let site = 1; site < observations.length; site += 1) {
    for (let state = 0; state < stateCount; state += 1) {
      for (let prior = 0; prior < stateCount; prior += 1) {
        terms[prior] = previous[prior] + Math.log(clampProbability(transition[prior * stateCount + state]));
      }
      current[state] = logSumExp(terms, stateCount)
        + Math.log(clampProbability(emission[state * 3 + observations[site]]));
    }
    [previous, current] = [current, previous];
  }
  return logSumExp(previous, stateCount);
}

function fitStateCount(observations, stateCount, options, seed) {
  let best = null;
  const randomStarts = Math.max(1, Math.min(64, Math.trunc(options.randomStarts ?? (options.sourceParity === true ? 21 : 10))));
  const maxIterations = Math.max(2, Math.min(250, Math.trunc(options.maxIterations ?? (options.sourceParity === true ? 100 : 60))));
  const sourceParity = options.sourceParity === true;
  const sourceRandom = sourceParity ? makeMsvcRng(seed) : null;
  let sourcePathMaximum = 0;
  for (let restart = 0; restart < randomStarts; restart += 1) {
    const random = sourceRandom ?? makeRng((seed ^ Math.imul(stateCount + 1, 0x9e3779b1) ^ Math.imul(restart + 1, 0x85ebca6b)) >>> 0);
    let model = initialiseModel(stateCount, random, options);
    let previousPath = null;
    let decoded = null;
    let iterations = 0;
    for (; iterations < maxIterations; iterations += 1) {
      decoded = viterbi(observations, model, sourceParity);
      if (sourceParity ? decoded.score === sourcePathMaximum : samePath(previousPath, decoded.path)) break;
      if (sourceParity) sourcePathMaximum = decoded.score;
      previousPath = decoded.path;
      model = updateModel(observations, decoded.path, stateCount, options);
    }
    decoded = decoded ?? viterbi(observations, model, sourceParity);
    let occupancy = new Uint32Array(stateCount);
    if (!sourceParity) {
      decoded = viterbi(observations, model);
      const updated = updateModel(observations, decoded.path, stateCount, options);
      model = { initial: updated.initial, transition: updated.transition, emission: updated.emission };
      occupancy = updated.occupancy;
    } else {
      for (const state of decoded.path) occupancy[state] += 1;
    }
    const logLikelihood = forwardLogLikelihood(observations, model);
    const selectionScore = sourceParity ? decoded.score : logLikelihood;
    if (!best || selectionScore > best.selectionScore) {
      best = { ...model, path: decoded.path, occupancy, logLikelihood, viterbiLogLikelihood: decoded.score, selectionScore, iterations: Math.min(maxIterations, iterations + 1), restart };
    }
  }
  const parameterCount = stateCount - 1 + stateCount * (stateCount - 1) + stateCount * 2;
  const bic = -2 * best.logLikelihood + parameterCount * Math.log(Math.max(2, observations.length));
  const aic = -2 * best.logLikelihood + 2 * parameterCount;
  return { ...best, stateCount, parameterCount, bic, aic };
}

function forwardBackward(observations, model) {
  const { initial, transition, emission } = model;
  const stateCount = initial.length;
  const siteCount = observations.length;
  const alpha = new Float64Array(siteCount * stateCount);
  const beta = new Float64Array(siteCount * stateCount);
  const terms = new Float64Array(stateCount);
  for (let state = 0; state < stateCount; state += 1) {
    alpha[state] = Math.log(clampProbability(initial[state]))
      + Math.log(clampProbability(emission[state * 3 + observations[0]]));
  }
  for (let site = 1; site < siteCount; site += 1) {
    const offset = site * stateCount;
    const priorOffset = offset - stateCount;
    for (let state = 0; state < stateCount; state += 1) {
      for (let prior = 0; prior < stateCount; prior += 1) {
        terms[prior] = alpha[priorOffset + prior]
          + Math.log(clampProbability(transition[prior * stateCount + state]));
      }
      alpha[offset + state] = logSumExp(terms, stateCount)
        + Math.log(clampProbability(emission[state * 3 + observations[site]]));
    }
  }
  const finalOffset = (siteCount - 1) * stateCount;
  const logLikelihood = logSumExp(alpha.subarray(finalOffset, finalOffset + stateCount));
  for (let state = 0; state < stateCount; state += 1) beta[finalOffset + state] = 0;
  for (let site = siteCount - 2; site >= 0; site -= 1) {
    const offset = site * stateCount;
    const nextOffset = offset + stateCount;
    for (let state = 0; state < stateCount; state += 1) {
      for (let next = 0; next < stateCount; next += 1) {
        terms[next] = Math.log(clampProbability(transition[state * stateCount + next]))
          + Math.log(clampProbability(emission[next * 3 + observations[site + 1]]))
          + beta[nextOffset + next];
      }
      beta[offset + state] = logSumExp(terms, stateCount);
    }
  }
  const posterior = new Float32Array(siteCount * stateCount);
  for (let site = 0; site < siteCount; site += 1) {
    const offset = site * stateCount;
    let total = 0;
    for (let state = 0; state < stateCount; state += 1) {
      const probability = Math.exp(alpha[offset + state] + beta[offset + state] - logLikelihood);
      posterior[offset + state] = Number.isFinite(probability) ? probability : 0;
      total += posterior[offset + state];
    }
    if (total > 0) for (let state = 0; state < stateCount; state += 1) posterior[offset + state] /= total;
  }
  return { posterior, logLikelihood };
}

function intervalOverlap(start, end, targetStart, targetEnd) {
  return Math.max(0, Math.min(end, targetEnd) - Math.max(start, targetStart));
}

function runBounds(positions, firstIndex, lastIndex, nSites) {
  const previous = firstIndex > 0 ? positions[firstIndex - 1] : -1;
  const first = positions[firstIndex];
  const last = positions[lastIndex];
  const next = lastIndex + 1 < positions.length ? positions[lastIndex + 1] : nSites;
  return {
    start: previous >= 0 ? Math.floor((previous + first + 1) / 2) : 0,
    end: next < nSites ? Math.floor((last + next + 1) / 2) : nSites,
  };
}

function chooseRun(path, positions, nSites, targetStart, targetEnd) {
  let best = null;
  let first = 0;
  while (first < path.length) {
    let last = first;
    while (last + 1 < path.length && path[last + 1] === path[first]) last += 1;
    const bounds = runBounds(positions, first, last, nSites);
    const overlap = intervalOverlap(bounds.start, bounds.end, targetStart, targetEnd);
    const targetLength = Math.max(1, targetEnd - targetStart);
    const union = Math.max(bounds.end, targetEnd) - Math.min(bounds.start, targetStart);
    const score = overlap / targetLength + overlap / Math.max(1, union);
    const candidate = { first, last, state: path[first], ...bounds, overlap, score };
    if (!best || candidate.score > best.score || (candidate.score === best.score && overlap > best.overlap)) best = candidate;
    first = last + 1;
  }
  return best;
}

function switchConfidence(positions, posterior, stateCount, boundaryIndex, fromState, toState, threshold, nSites) {
  let left = Math.max(0, boundaryIndex - 1);
  while (left > 0 && posterior[left * stateCount + fromState] < threshold) left -= 1;
  let right = Math.min(positions.length - 1, boundaryIndex);
  while (right + 1 < positions.length && posterior[right * stateCount + toState] < threshold) right += 1;
  const lowSite = posterior[left * stateCount + fromState] >= threshold
    ? positions[left]
    : (left > 0 ? positions[left - 1] : 0);
  const highSite = posterior[right * stateCount + toState] >= threshold
    ? positions[right]
    : (right + 1 < positions.length ? positions[right + 1] : nSites);
  return [Math.max(0, lowSite), Math.min(nSites, Math.max(lowSite, highSite))];
}

function switchConfidenceAnyState(positions, posterior, stateCount, boundaryIndex, threshold, nSites) {
  let left = Math.max(0, boundaryIndex - 1);
  const confident = (site) => {
    for (let state = 0; state < stateCount; state += 1) {
      if (posterior[site * stateCount + state] > threshold) return true;
    }
    return false;
  };
  while (left > 0 && !confident(left)) left -= 1;
  let right = Math.min(positions.length - 1, boundaryIndex);
  while (right + 1 < positions.length && !confident(right)) right += 1;
  const low = left > 0 ? positions[left - 1] + 1 : 0;
  const high = right + 1 < positions.length ? positions[right + 1] - 1 : nSites;
  return [Math.max(0, low), Math.min(nSites, Math.max(low, high))];
}

function confidenceContains(confidence, position) {
  return position >= confidence[0] && position <= confidence[1];
}

function nearestInformativeIndex(positions, target) {
  let low = 0;
  let high = positions.length - 1;
  while (low < high) {
    const middle = Math.floor((low + high) / 2);
    if (positions[middle] < target) low = middle + 1;
    else high = middle;
  }
  if (low > 0 && Math.abs(positions[low - 1] - target) <= Math.abs(positions[low] - target)) return low - 1;
  return low;
}

function selectSwitch(switches, target, positions) {
  if (!switches.length) return null;
  const inside = switches.filter((entry) => confidenceContains(entry.confidence95, target));
  const pool = inside.length ? inside : switches;
  const targetIndex = nearestInformativeIndex(positions, target);
  return [...pool].sort((left, right) =>
    Math.abs(left.informativeIndex - targetIndex) - Math.abs(right.informativeIndex - targetIndex)
      || (left.confidence95[1] - left.confidence95[0]) - (right.confidence95[1] - right.confidence95[0]),
  )[0];
}

function partitionHasEnoughInformation(positions, start, end, minimum = 3) {
  let inside = 0;
  for (const position of positions) if (position >= start && position < end) inside += 1;
  return inside >= minimum && positions.length - inside >= minimum;
}

function samplePosteriorTrace(positions, path, posterior, stateCount, maximumPoints = 240) {
  const stride = Math.max(1, Math.ceil(positions.length / maximumPoints));
  const trace = [];
  for (let site = 0; site < positions.length; site += stride) {
    const probabilities = [];
    for (let state = 0; state < stateCount; state += 1) probabilities.push(Number(posterior[site * stateCount + state].toFixed(5)));
    trace.push({ position: positions[site], state: path[site], probabilities });
  }
  if (trace.at(-1)?.position !== positions.at(-1)) {
    const site = positions.length - 1;
    const probabilities = [];
    for (let state = 0; state < stateCount; state += 1) probabilities.push(Number(posterior[site * stateCount + state].toFixed(5)));
    trace.push({ position: positions[site], state: path[site], probabilities });
  }
  return trace;
}

function pairCategory(sortedTriplet, leftSequence, rightSequence) {
  const left = sortedTriplet.indexOf(leftSequence);
  const right = sortedTriplet.indexOf(rightSequence);
  const pair = [left, right].sort((a, b) => a - b).join(":");
  if (pair === "0:1") return 0;
  if (pair === "0:2") return 2;
  return 1;
}

export function fitCategoricalHMM(observations, options = {}) {
  if (!(observations instanceof Uint8Array)) observations = Uint8Array.from(observations);
  if (observations.length < 8) return null;
  const sourceParity = options.sourceParity === true;
  const maximumStates = sourceParity
    ? 3
    : Math.max(2, Math.min(20, Math.trunc(options.maxStates ?? 20), Math.floor(observations.length / 3)));
  const criterionName = options.criterion === "aic" ? "aic" : "bic";
  const patience = Math.max(1, Math.min(19, Math.trunc(options.stepUpPatience ?? 2)));
  const exhaustive = options.exhaustiveModels === true;
  const seed = (options.seed ?? 0x5a17c0de) >>> 0;
  const ledger = [];
  let selected = null;
  let failures = 0;
  const firstStateCount = sourceParity ? 3 : 2;
  for (let stateCount = firstStateCount; stateCount <= maximumStates; stateCount += 1) {
    const fit = fitStateCount(observations, stateCount, options, seed);
    ledger.push({
      states: stateCount,
      logLikelihood: fit.logLikelihood,
      bic: fit.bic,
      aic: fit.aic,
      iterations: fit.iterations,
      winningRestart: fit.restart + 1,
    });
    if (!selected || fit[criterionName] < selected[criterionName] - 1e-8) {
      selected = fit;
      failures = 0;
    } else {
      failures += 1;
      if (!exhaustive && failures >= patience) break;
    }
  }
  return { selected, ledger, criterion: criterionName };
}

export function fitBurtTriplet(encoded, nSites, sequence1, sequence2, sequence3, candidateStart, candidateEnd, options = {}) {
  const sortedTriplet = [sequence1, sequence2, sequence3].sort((left, right) => left - right);
  const collected = collectBurtObservations(encoded, nSites, sortedTriplet[0], sortedTriplet[1], sortedTriplet[2]);
  const realPositions = collected.positions;
  let { positions, observations } = collected;
  if (observations.length < 8 || candidateEnd <= candidateStart) return null;
  const sourceParity = options.sourceParity !== false;
  if (sourceParity) {
    // BenHMM uses one-based XDiffPos alongside a zero-based RecodeB array and
    // passes SLen as the count (therefore processing its zero-filled sentinel).
    // Preserve that observable offset only in source-compatibility mode.
    const sourcePositions = new Int32Array(positions.length + 1);
    const sourceObservations = new Uint8Array(observations.length + 1);
    sourceObservations.set(observations);
    for (let index = 0; index < positions.length; index += 1) sourcePositions[index + 1] = positions[index];
    positions = sourcePositions;
    observations = sourceObservations;
  }
  const fitOptions = { ...options, sourceParity, alignmentLength: nSites };
  const fitted = fitCategoricalHMM(observations, fitOptions);
  if (!fitted?.selected) return null;
  const model = fitted.selected;
  const posteriorResult = forwardBackward(observations, model);
  let selectedRun = chooseRun(model.path, positions, nSites, candidateStart, candidateEnd);
  if (!selectedRun || selectedRun.overlap <= 0 || selectedRun.last - selectedRun.first + 1 < 2) {
    selectedRun = { first: 0, last: model.path.length - 1, state: model.path[Math.floor(model.path.length / 2)], start: candidateStart, end: candidateEnd, overlap: candidateEnd - candidateStart };
  }
  const threshold = Math.max(0.5, Math.min(0.9999, options.posteriorThreshold ?? (sourceParity ? 0.995 : 0.95)));
  const startFrom = selectedRun.first > 0 ? model.path[selectedRun.first - 1] : selectedRun.state;
  const endTo = selectedRun.last + 1 < model.path.length ? model.path[selectedRun.last + 1] : selectedRun.state;
  let confidenceStart = selectedRun.first > 0
    ? (sourceParity
      ? switchConfidenceAnyState(positions, posteriorResult.posterior, model.stateCount, selectedRun.first, threshold, nSites)
      : switchConfidence(positions, posteriorResult.posterior, model.stateCount, selectedRun.first, startFrom, selectedRun.state, threshold, nSites))
    : [0, positions[0]];
  let confidenceEnd = selectedRun.last + 1 < model.path.length
    ? (sourceParity
      ? switchConfidenceAnyState(positions, posteriorResult.posterior, model.stateCount, selectedRun.last + 1, threshold, nSites)
      : switchConfidence(positions, posteriorResult.posterior, model.stateCount, selectedRun.last + 1, selectedRun.state, endTo, threshold, nSites))
    : [positions.at(-1), nSites];
  const switches = [];
  const switchFirst = sourceParity ? 3 : 1;
  const switchLast = sourceParity ? model.path.length - 1 : model.path.length;
  for (let site = switchFirst; site < switchLast; site += 1) {
    if (model.path[site] === model.path[site - 1]) continue;
    const ci = sourceParity
      ? switchConfidenceAnyState(positions, posteriorResult.posterior, model.stateCount, site, threshold, nSites)
      : switchConfidence(positions, posteriorResult.posterior, model.stateCount, site, model.path[site - 1], model.path[site], threshold, nSites);
    const ci99 = sourceParity
      ? switchConfidenceAnyState(positions, posteriorResult.posterior, model.stateCount, site, 0.999, nSites)
      : switchConfidence(positions, posteriorResult.posterior, model.stateCount, site, model.path[site - 1], model.path[site], 0.99, nSites);
    switches.push({
      position: Math.floor((positions[site - 1] + positions[site] + 1) / 2),
      informativeIndex: site,
      fromState: model.path[site - 1],
      toState: model.path[site],
      confidence95: ci,
      confidence99: ci99,
    });
  }
  let selectedStart = selectedRun.start;
  let selectedEnd = selectedRun.end;
  if (sourceParity && switches.length >= 2) {
    const leftSwitch = selectSwitch(switches, candidateStart, positions);
    const rightSwitch = selectSwitch(switches, candidateEnd, positions);
    if (leftSwitch && rightSwitch && rightSwitch.position > leftSwitch.position) {
      selectedStart = leftSwitch.position;
      selectedEnd = rightSwitch.position;
      confidenceStart = leftSwitch.confidence95;
      confidenceEnd = rightSwitch.confidence95;
    }
  }
  // PolishBP rejects a proposed two-breakpoint partition when either side has
  // too little variable-site evidence. Preserve the original candidate rather
  // than returning an overfit HMM switch pair in that case.
  if (!partitionHasEnoughInformation(realPositions, selectedStart, selectedEnd)) {
    selectedStart = candidateStart;
    selectedEnd = candidateEnd;
  }
  const selectedLeftSwitch = selectSwitch(switches, selectedStart, positions);
  const selectedRightSwitch = selectSwitch(switches, selectedEnd, positions);
  if (selectedLeftSwitch) confidenceStart = selectedLeftSwitch.confidence95;
  if (selectedRightSwitch) confidenceEnd = selectedRightSwitch.confidence95;
  const emissions = [];
  const transitions = [];
  for (let state = 0; state < model.stateCount; state += 1) {
    emissions.push(Array.from(model.emission.subarray(state * 3, state * 3 + 3), (value) => Number(value.toFixed(6))));
    transitions.push(Array.from(model.transition.subarray(state * model.stateCount, (state + 1) * model.stateCount), (value) => Number(value.toFixed(6))));
  }
  const majorCategory = pairCategory(sortedTriplet, sequence1, sequence2);
  const minorCategory = pairCategory(sortedTriplet, sequence1, sequence3);
  let majorState = 0;
  for (let state = 1; state < model.stateCount; state += 1) {
    if (model.emission[state * 3 + majorCategory] > model.emission[majorState * 3 + majorCategory]) majorState = state;
  }
  return {
    start: selectedStart,
    end: selectedEnd,
    confidenceStart,
    confidenceEnd,
    model: {
      method: "burt-hmm",
      informativeSites: collected.observations.length,
      states: model.stateCount,
      stateSwitches: switches.length,
      majorFit: model.emission[majorState * 3 + majorCategory],
      minorFit: model.emission[selectedRun.state * 3 + minorCategory],
      logLikelihood: posteriorResult.logLikelihood,
      bic: model.bic,
      aic: model.aic,
      criterion: sourceParity ? "RDP5 source (maximum Viterbi likelihood)" : fitted.criterion.toUpperCase(),
      randomStarts: Math.max(1, Math.trunc(options.randomStarts ?? (sourceParity ? 21 : 10))),
      iterations: model.iterations,
      selectedState: selectedRun.state,
      posteriorThreshold: threshold,
      sourceParity,
      sourceCompatibility: sourceParity ? "RDP5 BenHMM + DoHMMCyclesSerial" : undefined,
      confidence99Start: selectedLeftSwitch?.confidence99,
      confidence99End: selectedRightSwitch?.confidence99,
      emissions,
      transitions,
      switches,
      posteriorTrace: samplePosteriorTrace(positions, model.path, posteriorResult.posterior, model.stateCount),
      modelSelection: fitted.ledger,
    },
  };
}
