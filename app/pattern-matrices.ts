import type { RdpEvent } from "./rdp-core";

export interface BreakpointPairDensity {
  values: Float32Array;
  maximum: number;
}

export interface RegionSeparationMatrices {
  observed: Float32Array;
  standardizedResidual: Float32Array;
  maximumObserved: number;
  maximumAbsoluteResidual: number;
}

export interface LocalDiscordanceMatrices {
  rmsDeviation: Float32Array;
  correlationLoss: Float32Array;
  maximumRmsDeviation: number;
  maximumCorrelationLoss: number;
  sequenceIndexes: number[];
  pairCount: number;
  sampledSitesPerWindow: number;
}

function canonical(base: string | undefined): boolean {
  return base === "A" || base === "C" || base === "G" || base === "T";
}

function eventTractLength(event: RdpEvent, length: number): number {
  if (length <= 0) return 0;
  return event.wraps
    ? Math.max(0, length - event.start + event.end)
    : Math.max(0, event.end - event.start);
}

export function eventContainsPosition(event: RdpEvent, position: number): boolean {
  return event.wraps
    ? position >= event.start || position < event.end
    : position >= event.start && position < event.end;
}

function binForBoundary(position: number, length: number, resolution: number): number {
  if (length <= 0) return 0;
  const bounded = Math.max(0, Math.min(length - Number.EPSILON, position));
  return Math.max(0, Math.min(resolution - 1, Math.floor((bounded / length) * resolution)));
}

export function computeBreakpointPairDensity(
  events: RdpEvent[],
  length: number,
  resolution: number,
): BreakpointPairDensity {
  const values = new Float32Array(resolution * resolution);
  let maximum = 0;
  for (const event of events) {
    const startBin = binForBoundary(event.start, length, resolution);
    const endBoundary = event.end >= length ? length - Number.EPSILON : event.end;
    const endBin = binForBoundary(endBoundary, length, resolution);
    const forward = startBin * resolution + endBin;
    values[forward] += 1;
    maximum = Math.max(maximum, values[forward]);
    if (startBin !== endBin) {
      const reverse = endBin * resolution + startBin;
      values[reverse] += 1;
      maximum = Math.max(maximum, values[reverse]);
    }
  }
  return { values, maximum };
}

export function computeRegionSeparationMatrices(
  events: RdpEvent[],
  length: number,
  resolution: number,
): RegionSeparationMatrices {
  const observed = new Float32Array(resolution * resolution);
  const standardizedResidual = new Float32Array(resolution * resolution);
  if (length <= 0 || resolution <= 1 || events.length === 0) {
    return { observed, standardizedResidual, maximumObserved: 0, maximumAbsoluteResidual: 0 };
  }

  const membership = events.map((event) => {
    const bins = new Uint8Array(resolution);
    for (let bin = 0; bin < resolution; bin += 1) {
      const position = ((bin + 0.5) / resolution) * length;
      bins[bin] = eventContainsPosition(event, position) ? 1 : 0;
    }
    return bins;
  });
  const expectedByDistance = new Float64Array(Math.floor(resolution / 2) + 1);
  const varianceByDistance = new Float64Array(expectedByDistance.length);
  for (const event of events) {
    const tract = Math.min(length, eventTractLength(event, length));
    for (let distanceBins = 1; distanceBins < expectedByDistance.length; distanceBins += 1) {
      const distance = (distanceBins / resolution) * length;
      const probability = Math.max(0, Math.min(1, (2 * Math.min(distance, tract, length - tract)) / length));
      expectedByDistance[distanceBins] += probability;
      varianceByDistance[distanceBins] += probability * (1 - probability);
    }
  }

  let maximumObserved = 0;
  let maximumAbsoluteResidual = 0;
  for (let row = 0; row < resolution; row += 1) {
    for (let column = row + 1; column < resolution; column += 1) {
      let count = 0;
      for (const bins of membership) count += bins[row] === bins[column] ? 0 : 1;
      const directDistance = column - row;
      const distanceBins = Math.min(directDistance, resolution - directDistance);
      const expected = expectedByDistance[distanceBins];
      const variance = varianceByDistance[distanceBins];
      const residual = variance > 1e-9 ? (count - expected) / Math.sqrt(variance) : 0;
      const forward = row * resolution + column;
      const reverse = column * resolution + row;
      observed[forward] = count;
      observed[reverse] = count;
      standardizedResidual[forward] = residual;
      standardizedResidual[reverse] = residual;
      maximumObserved = Math.max(maximumObserved, count);
      maximumAbsoluteResidual = Math.max(maximumAbsoluteResidual, Math.abs(residual));
    }
  }
  return { observed, standardizedResidual, maximumObserved, maximumAbsoluteResidual };
}

function evenlySpacedIndexes(total: number, requested: number): number[] {
  if (total <= requested) return Array.from({ length: total }, (_, index) => index);
  if (requested <= 1) return [0];
  return [...new Set(Array.from({ length: requested }, (_, index) => Math.round((index / (requested - 1)) * (total - 1))))];
}

function correlationLoss(left: Float32Array, right: Float32Array): number {
  if (left.length === 0 || left.length !== right.length) return 0;
  let leftMean = 0;
  let rightMean = 0;
  for (let index = 0; index < left.length; index += 1) {
    leftMean += left[index];
    rightMean += right[index];
  }
  leftMean /= left.length;
  rightMean /= right.length;
  let covariance = 0;
  let leftVariance = 0;
  let rightVariance = 0;
  for (let index = 0; index < left.length; index += 1) {
    const leftDelta = left[index] - leftMean;
    const rightDelta = right[index] - rightMean;
    covariance += leftDelta * rightDelta;
    leftVariance += leftDelta * leftDelta;
    rightVariance += rightDelta * rightDelta;
  }
  if (leftVariance <= 1e-12 || rightVariance <= 1e-12) {
    return leftVariance <= 1e-12 && rightVariance <= 1e-12 && Math.abs(leftMean - rightMean) <= 1e-9 ? 0 : 1;
  }
  const correlation = covariance / Math.sqrt(leftVariance * rightVariance);
  return Math.max(0, Math.min(2, 1 - correlation));
}

export function computeLocalDiscordanceMatrices(
  sequences: string[],
  length: number,
  resolution: number,
  maximumSequences = 18,
  maximumSitesPerWindow = 72,
): LocalDiscordanceMatrices {
  const sequenceIndexes = evenlySpacedIndexes(sequences.length, Math.max(3, maximumSequences));
  const pairs: Array<[number, number]> = [];
  for (let left = 0; left < sequenceIndexes.length; left += 1) {
    for (let right = left + 1; right < sequenceIndexes.length; right += 1) pairs.push([sequenceIndexes[left], sequenceIndexes[right]]);
  }
  const profiles: Float32Array[] = [];
  let sampledSitesPerWindow = 0;
  for (let windowIndex = 0; windowIndex < resolution; windowIndex += 1) {
    const start = Math.floor((windowIndex / resolution) * length);
    const end = Math.max(start + 1, Math.floor(((windowIndex + 1) / resolution) * length));
    const windowLength = Math.max(1, end - start);
    const sampleCount = Math.min(maximumSitesPerWindow, windowLength);
    sampledSitesPerWindow = Math.max(sampledSitesPerWindow, sampleCount);
    const sites = Array.from({ length: sampleCount }, (_, sampleIndex) => Math.min(end - 1, Math.floor(start + ((sampleIndex + 0.5) / sampleCount) * windowLength)));
    const profile = new Float32Array(pairs.length);
    pairs.forEach(([leftIndex, rightIndex], pairIndex) => {
      let callable = 0;
      let different = 0;
      for (const site of sites) {
        const leftBase = sequences[leftIndex]?.[site];
        const rightBase = sequences[rightIndex]?.[site];
        if (!canonical(leftBase) || !canonical(rightBase)) continue;
        callable += 1;
        if (leftBase !== rightBase) different += 1;
      }
      profile[pairIndex] = callable > 0 ? different / callable : 0;
    });
    profiles.push(profile);
  }

  const rmsDeviation = new Float32Array(resolution * resolution);
  const correlationLossValues = new Float32Array(resolution * resolution);
  let maximumRmsDeviation = 0;
  let maximumCorrelationLoss = 0;
  for (let row = 0; row < resolution; row += 1) {
    for (let column = row + 1; column < resolution; column += 1) {
      let squaredDifference = 0;
      for (let pairIndex = 0; pairIndex < pairs.length; pairIndex += 1) {
        const difference = profiles[row][pairIndex] - profiles[column][pairIndex];
        squaredDifference += difference * difference;
      }
      const rms = pairs.length > 0 ? Math.sqrt(squaredDifference / pairs.length) : 0;
      const loss = correlationLoss(profiles[row], profiles[column]);
      const forward = row * resolution + column;
      const reverse = column * resolution + row;
      rmsDeviation[forward] = rms;
      rmsDeviation[reverse] = rms;
      correlationLossValues[forward] = loss;
      correlationLossValues[reverse] = loss;
      maximumRmsDeviation = Math.max(maximumRmsDeviation, rms);
      maximumCorrelationLoss = Math.max(maximumCorrelationLoss, loss);
    }
  }
  return {
    rmsDeviation,
    correlationLoss: correlationLossValues,
    maximumRmsDeviation,
    maximumCorrelationLoss,
    sequenceIndexes,
    pairCount: pairs.length,
    sampledSitesPerWindow,
  };
}

