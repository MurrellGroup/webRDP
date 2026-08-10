// Clean-room probability calibrations for the RDP Web method-family kernels.
// This module is deliberately dependency-free so the worker remains a single
// static-site download and can run unchanged from a GitHub Pages subpath.

export const METHODS = ["RDP", "GENECONV", "BootScan", "MaxChi", "Chimaera", "SiScan", "3Seq"];

function clampProbability(value) {
  if (!Number.isFinite(value)) return 1;
  return Math.max(Number.MIN_VALUE, Math.min(1, value));
}

export function erfc(value) {
  const z = Math.abs(value);
  const t = 1 / (1 + z / 2);
  const answer = t * Math.exp(
    -z * z - 1.26551223
      + t * (1.00002368
        + t * (0.37409196
          + t * (0.09678418
            + t * (-0.18628806
              + t * (0.27886807
                + t * (-1.13520398
                  + t * (1.48851587
                    + t * (-0.82215223 + t * 0.17087277)))))))),
  );
  return value >= 0 ? answer : 2 - answer;
}

export function chiSquareP(statistic) {
  return clampProbability(erfc(Math.sqrt(Math.max(0, statistic) / 2)));
}

export function normalTwoSided(z) {
  return clampProbability(erfc(Math.abs(z) / Math.SQRT2));
}

// Lanczos log-gamma and the regularized incomplete beta give a stable exact
// binomial tail without summing thousands of tiny probability masses.
function logGamma(value) {
  const coefficients = [
    676.5203681218851,
    -1259.1392167224028,
    771.3234287776531,
    -176.6150291621406,
    12.507343278686905,
    -0.13857109526572012,
    9.984369578019572e-6,
    1.5056327351493116e-7,
  ];
  if (value < 0.5) {
    return Math.log(Math.PI) - Math.log(Math.sin(Math.PI * value)) - logGamma(1 - value);
  }
  const z = value - 1;
  let series = 0.9999999999998099;
  for (let index = 0; index < coefficients.length; index += 1) {
    series += coefficients[index] / (z + index + 1);
  }
  const shifted = z + coefficients.length - 0.5;
  return 0.5 * Math.log(2 * Math.PI) + (z + 0.5) * Math.log(shifted) - shifted + Math.log(series);
}

function betaFraction(a, b, x) {
  const maxIterations = 220;
  const epsilon = 3e-14;
  const tiny = 1e-300;
  const qab = a + b;
  const qap = a + 1;
  const qam = a - 1;
  let c = 1;
  let d = 1 - (qab * x) / qap;
  if (Math.abs(d) < tiny) d = tiny;
  d = 1 / d;
  let result = d;
  for (let iteration = 1; iteration <= maxIterations; iteration += 1) {
    const doubled = 2 * iteration;
    let numerator = (iteration * (b - iteration) * x)
      / ((qam + doubled) * (a + doubled));
    d = 1 + numerator * d;
    if (Math.abs(d) < tiny) d = tiny;
    c = 1 + numerator / c;
    if (Math.abs(c) < tiny) c = tiny;
    d = 1 / d;
    result *= d * c;
    numerator = -((a + iteration) * (qab + iteration) * x)
      / ((a + doubled) * (qap + doubled));
    d = 1 + numerator * d;
    if (Math.abs(d) < tiny) d = tiny;
    c = 1 + numerator / c;
    if (Math.abs(c) < tiny) c = tiny;
    d = 1 / d;
    const delta = d * c;
    result *= delta;
    if (Math.abs(delta - 1) < epsilon) break;
  }
  return result;
}

function regularizedBeta(x, a, b) {
  if (x <= 0) return 0;
  if (x >= 1) return 1;
  const front = Math.exp(
    logGamma(a + b) - logGamma(a) - logGamma(b)
      + a * Math.log(x) + b * Math.log1p(-x),
  );
  if (x < (a + 1) / (a + b + 2)) {
    return clampProbability((front * betaFraction(a, b, x)) / a);
  }
  return clampProbability(1 - (front * betaFraction(b, a, 1 - x)) / b);
}

export function binomialUpper(successes, trials, probability) {
  if (successes <= 0) return 1;
  if (successes > trials) return Number.MIN_VALUE;
  if (trials <= 0) return 1;
  const p = Math.max(1e-12, Math.min(1 - 1e-12, probability));
  return regularizedBeta(p, successes, trials - successes + 1);
}

const threeSeqCache = new Map();

// Exact tail probability for the maximum descent of a hypergeometric random
// walk with a fixed number of up/down steps. The dynamic program follows the
// drawdown from the running maximum and accumulates first-passage probability
// directly, avoiding catastrophic cancellation for very small tails. The
// bounded work guard keeps browser analyses responsive; larger cases retain a
// conservative finite-sample bound until lookup-table generation is added.
export function threeSeqExactP(upSteps, downSteps, observedDescent, maxOperations = 4_000_000) {
  const up = Math.max(0, Math.trunc(upSteps));
  const down = Math.max(0, Math.trunc(downSteps));
  const threshold = Math.max(0, Math.trunc(observedDescent));
  if (threshold <= 0) return { p: 1, exact: true };
  if (threshold > down || down === 0) return { p: Number.MIN_VALUE, exact: true };
  const operations = (up + 1) * (down + 1) * threshold;
  if (operations > maxOperations) return { p: null, exact: false };
  const key = `${up}:${down}:${threshold}`;
  const cached = threeSeqCache.get(key);
  if (cached !== undefined) return { p: cached, exact: true };

  const stride = threshold;
  let current = new Float64Array((up + 1) * stride);
  let next = new Float64Array(current.length);
  current[0] = 1;
  let tail = 0;
  const total = up + down;
  for (let step = 0; step < total; step += 1) {
    next.fill(0);
    const firstUp = Math.max(0, step - down);
    const lastUp = Math.min(up, step);
    for (let usedUp = firstUp; usedUp <= lastUp; usedUp += 1) {
      const usedDown = step - usedUp;
      const remaining = total - step;
      const upProbability = (up - usedUp) / remaining;
      const downProbability = (down - usedDown) / remaining;
      const offset = usedUp * stride;
      for (let drawdown = 0; drawdown < stride; drawdown += 1) {
        const probability = current[offset + drawdown];
        if (probability === 0) continue;
        if (upProbability > 0) {
          const nextDrawdown = drawdown > 0 ? drawdown - 1 : 0;
          next[(usedUp + 1) * stride + nextDrawdown] += probability * upProbability;
        }
        if (downProbability > 0) {
          if (drawdown + 1 >= threshold) tail += probability * downProbability;
          else next[offset + drawdown + 1] += probability * downProbability;
        }
      }
    }
    [current, next] = [next, current];
  }
  const p = clampProbability(tail);
  threeSeqCache.set(key, p);
  return { p, exact: true };
}

function identityRunUpperBound(run, eligible, backgroundMatches) {
  if (run <= 0 || eligible <= 0) return 1;
  const background = Math.max(1e-9, Math.min(1 - 1e-9, backgroundMatches / eligible));
  return clampProbability(Math.max(1, eligible - run + 1) * Math.pow(background, run));
}

function corrected(raw, correction, comparisons) {
  return correction === "none"
    ? clampProbability(raw)
    : clampProbability(raw * Math.max(1, comparisons));
}

export function methodEvidence(candidate, stats, options, comparisons, nSites) {
  const insideTotal = candidate.insideMinor + candidate.insideMajor;
  const outsideTotal = candidate.outsideMinor + candidate.outsideMajor;
  const backgroundMinor = candidate.outsideMinor / Math.max(1, outsideTotal);
  const windowTests = Math.max(1, Math.floor(Math.max(0, nSites - options.window) / Math.max(1, options.step)) + 1);
  const partitionTests = Math.max(1, nSites - 1);
  const sisterZ = stats.siskanScore / Math.sqrt(Math.max(1, stats.siskanSites));
  const descent = stats.threeSeqDescent ?? stats.threeSeqBridge / 1000;
  const exactThreeSeq = threeSeqExactP(
    stats.threeSeqMajorSites ?? Math.floor(stats.threeSeqSites / 2),
    stats.threeSeqMinorSites ?? Math.ceil(stats.threeSeqSites / 2),
    descent,
    options.threeSeqMaxOperations ?? 4_000_000,
  );
  const threeSeqBound = Math.min(1, 2 * Math.exp(
    (-2 * descent * descent) / Math.max(1, stats.threeSeqSites),
  ));
  const bootstrapAvailable = stats.bootscanBootstrapReplicates > 0;
  const bootstrapP = bootstrapAvailable
    ? binomialUpper(stats.bootscanBootstrapConsistent, stats.bootscanBootstrapReplicates, 0.5)
    : 1;
  const windowSignP = binomialUpper(stats.bootscanConsistent, stats.bootscanWindows, 0.5);

  const calculations = {
    RDP: {
      p: clampProbability(binomialUpper(candidate.insideMinor, insideTotal, backgroundMinor) * windowTests),
      statistic: candidate.insideMinor / Math.max(1, insideTotal) - backgroundMinor,
      statisticLabel: "identity shift",
      calibration: "binomial · window-corrected",
    },
    GENECONV: {
      p: identityRunUpperBound(stats.genconvRun, stats.genconvEligible, stats.genconvMatches),
      statistic: stats.genconvRun,
      statisticLabel: "concordant run",
      calibration: "G-scale 0 run bound",
    },
    BootScan: {
      p: bootstrapAvailable ? Math.max(windowSignP, bootstrapP) : windowSignP,
      statistic: bootstrapAvailable
        ? stats.bootscanBootstrapConsistent / Math.max(1, stats.bootscanBootstrapReplicates)
        : stats.bootscanConsistent / Math.max(1, stats.bootscanWindows),
      statisticLabel: bootstrapAvailable ? "bootstrap topology support" : "topology agreement",
      calibration: bootstrapAvailable
        ? "seeded p-distance bootstrap + window sign"
        : "distance-window sign test",
    },
    MaxChi: {
      p: clampProbability(chiSquareP(stats.maxChi) * partitionTests),
      statistic: stats.maxChi,
      statisticLabel: "minimum boundary χ²",
      calibration: "χ² · partition-corrected",
    },
    Chimaera: {
      p: clampProbability(chiSquareP(stats.chimaera) * partitionTests),
      statistic: stats.chimaera,
      statisticLabel: "minimum boundary χ²",
      calibration: "binary-triplet χ²",
    },
    SiScan: {
      p: clampProbability(normalTwoSided(sisterZ) * windowTests),
      statistic: sisterZ,
      statisticLabel: "oriented category Z",
      calibration: "fast category-Z surrogate",
    },
    "3Seq": {
      p: exactThreeSeq.p ?? clampProbability(threeSeqBound),
      statistic: descent,
      statisticLabel: "maximum HGRW descent",
      calibration: exactThreeSeq.exact ? "exact HGRW first-passage DP" : "finite-sample exponential bound",
    },
  };

  return METHODS
    .filter((method) => options.methods.includes(method))
    .map((method) => {
      const calculation = calculations[method];
      const pValue = clampProbability(calculation.p);
      const correctedP = corrected(pValue, options.correction, comparisons);
      return {
        method,
        pValue,
        correctedP,
        score: -Math.log10(pValue),
        supported: correctedP <= options.alpha,
        statistic: calculation.statistic,
        statisticLabel: calculation.statisticLabel,
        calibration: calculation.calibration,
      };
    });
}
