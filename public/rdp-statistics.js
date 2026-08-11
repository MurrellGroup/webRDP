// Probability calibrations for the RDP Web method-family kernels. Source-mode
// calculations port the RDP5 implementation supplied by the original authors;
// this dependency-free module still runs unchanged from a GitHub Pages subpath.

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

export function studentTTwoSided(tStatistic, degreesOfFreedom) {
  const t = Math.abs(Number(tStatistic));
  const df = Number(degreesOfFreedom);
  if (!Number.isFinite(t) || !(df > 0)) return 1;
  if (t === 0) return 1;
  return clampProbability(regularizedBeta(df / (df + t * t), df / 2, 0.5));
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

// Direct G-scale 0 specialization of RDP5 CalcKMaxP + GCCalcPValP.
// In the source, lambda=-log(Q), K=P and the Karlin-Altschul-like
// probability is 1-exp[-exp{-(lambda*S-log(K*L))}]. With P=mismatches/L
// and Q=matches/L this reduces to 1-exp[-mismatches*Q^S].
export function geneconvSourceG0Probability(run, eligible, matchingSites) {
  const length = Math.trunc(eligible);
  const matches = Math.trunc(matchingSites);
  const mismatches = length - matches;
  if (run <= 3 || length <= 0 || matches <= 0 || mismatches <= 0) return 1;
  const q = matches / length;
  const poissonMean = mismatches * Math.pow(q, run);
  return clampProbability(-Math.expm1(-poissonMean));
}

export function geneconvSourceProbability(score, eligible, matchingSites, gScale = 1) {
  if (!(gScale > 0)) return geneconvSourceG0Probability(score, eligible, matchingSites);
  const length = Math.trunc(eligible);
  const matches = Math.trunc(matchingSites);
  const mismatches = length - matches;
  if (score <= 3 || length <= 0 || matches <= 0 || mismatches <= 0) return 1;
  const mismatchProbability = mismatches / length;
  const matchProbability = matches / length;
  const mismatchPenalty = Math.floor(length * gScale / mismatches) + 1;
  const weightedMismatch = mismatchPenalty * mismatchProbability;
  if (!(weightedMismatch > matchProbability)) return 1;
  let z = Math.exp((2 * Math.log(weightedMismatch / matchProbability)) / (mismatchPenalty + 1));
  for (let iteration = 0; iteration < 10_000; iteration += 1) {
    const inversePower = Math.pow(z, -mismatchPenalty);
    const residual = matchProbability * z + mismatchProbability * inversePower - 1;
    const derivative = matchProbability - weightedMismatch * inversePower / z;
    if (!(Math.abs(derivative) > 1e-15)) return 1;
    const delta = residual / derivative;
    z -= delta;
    if (!(z > 1) || !Number.isFinite(z)) return 1;
    if (Math.abs(delta) <= 1e-6 && Math.abs(residual) <= 1e-6) break;
  }
  const lambda = Math.log(z);
  const k = (z - 1) * (
    matchProbability
      - weightedMismatch * Math.exp(-(mismatchPenalty + 1) * lambda)
  );
  if (!(lambda > 0) || !(k > 0)) return 1;
  const kaScore = lambda * score - Math.log(k * length);
  if (!(kaScore > 0)) return 1;
  const poissonMean = Math.exp(-Math.min(kaScore, 745));
  return clampProbability(kaScore < 32 ? -Math.expm1(-poissonMean) : poissonMean);
}

// FastRecCheckMC/FastRecCheckChim apply the one-degree-of-freedom chi-square
// tail to the best peak, multiply by the number of half-window placements in
// compressed informative-site space, then by the three pairwise triplet
// orientations. Experiment-wide correction is applied separately below.
export function sourceChiWindowProbability(statistic, informativeSites, fullWindow, exactHalfWindow = null) {
  const halfWindow = exactHalfWindow === null
    ? Math.max(8, Math.floor(fullWindow / 2))
    : Math.max(1, Math.floor(exactHalfWindow));
  const placements = Math.max(1, informativeSites / halfWindow);
  return clampProbability(chiSquareP(statistic) * placements * 3);
}

export function rdp5SourceProbability(source) {
  if (!source || source.tractSites <= 0 || source.common <= 0 || source.mediumSites <= 0) return null;
  const originalLength = Math.trunc(source.tractSites);
  let tractLength = originalLength;
  let common = Math.trunc(source.common);
  let exponent = 1;
  if (tractLength >= 170) {
    exponent = tractLength / 169;
    const different = Math.round((tractLength - common) * 169 / tractLength);
    tractLength = 169;
    common = tractLength - different;
  }
  const identity = Math.max(1e-12, Math.min(1 - 1e-12, source.mediumSites / Math.max(1, source.informativeSites ?? 1)));
  const probabilitySites = Math.max(1, source.probabilitySites ?? source.informativeSites ?? 1);
  let probability = binomialUpper(common, tractLength, identity) * (probabilitySites / tractLength);
  if (exponent > 1.000001) probability = probability > 0 ? Math.pow(probability, exponent) : 0.05;
  return clampProbability(Math.max(1e-300, probability));
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
  const sourceSiScanAvailable = Number.isFinite(stats.siskanSourceZ)
    && Number.isFinite(stats.siskanSourceP)
    && typeof stats.siskanSourceRoutine === "string";
  const sisterZ = sourceSiScanAvailable
    ? stats.siskanSourceZ
    : stats.siskanScore / Math.sqrt(Math.max(1, stats.siskanSites));
  const descent = stats.threeSeqDescent ?? stats.threeSeqBridge / 1000;
  const exactThreeSeq = options.methods.includes("3Seq")
    ? threeSeqExactP(
        stats.threeSeqMajorSites ?? Math.floor(stats.threeSeqSites / 2),
        stats.threeSeqMinorSites ?? Math.ceil(stats.threeSeqSites / 2),
        descent,
        options.threeSeqMaxOperations ?? 4_000_000,
      )
    : { p: null, exact: false };
  const threeSeqBound = Math.min(1, 2 * Math.exp(
    (-2 * descent * descent) / Math.max(1, stats.threeSeqSites),
  ));
  const bootstrapAvailable = stats.bootscanBootstrapReplicates > 0;
  const bootstrapP = bootstrapAvailable
    ? binomialUpper(stats.bootscanBootstrapConsistent, stats.bootscanBootstrapReplicates, 0.5)
    : 1;
  const windowSignP = binomialUpper(stats.bootscanConsistent, stats.bootscanWindows, 0.5);
  const sourceRdpP = rdp5SourceProbability(stats.rdpSource
    ? { ...stats.rdpSource, informativeSites: candidate.informative }
    : null);
  const hasSignalLedger = Array.isArray(candidate.methodSignals);
  const locatedMethods = new Set((candidate.methodSignals ?? []).map((signal) => signal.method));
  const maxChiSignal = (candidate.methodSignals ?? []).find((signal) => signal.method === "MaxChi" && signal.sourceChi);
  const chimaeraSignal = (candidate.methodSignals ?? []).find((signal) => signal.method === "Chimaera" && signal.sourceChi);

  const calculations = {
    RDP: {
      p: sourceRdpP ?? clampProbability(binomialUpper(candidate.insideMinor, insideTotal, backgroundMinor) * windowTests),
      statistic: candidate.insideMinor / Math.max(1, insideTotal) - backgroundMinor,
      statisticLabel: "identity shift",
      calibration: sourceRdpP === null ? "binomial · window-corrected fallback" : "RDP5 ProbCalcP/P2-equivalent binomial tail",
    },
    GENECONV: {
      p: geneconvSourceProbability(stats.genconvRun, stats.genconvEligible, stats.genconvMatches, options.geneconvGScale ?? 1),
      statistic: stats.genconvRun,
      statisticLabel: "fragment score",
      calibration: `RDP5 CalcKMaxP/GCCalcPValP · G-scale ${options.geneconvGScale ?? 1}`,
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
      p: sourceChiWindowProbability(
        stats.maxChi,
        maxChiSignal?.sourceChi?.informativeSites ?? stats.genconvEligible,
        options.window,
        maxChiSignal?.sourceChi?.halfWindow ?? null,
      ),
      statistic: stats.maxChi,
      statisticLabel: "minimum boundary χ²",
      calibration: maxChiSignal?.sourceChi
        ? `RDP5 ChiPVal2P · ${maxChiSignal.sourceChi.informativeSites} compressed sites · half-window ${maxChiSignal.sourceChi.halfWindow} × 3`
        : "RDP5 ChiPVal2P · informative half-window × 3",
    },
    Chimaera: {
      p: sourceChiWindowProbability(
        stats.chimaera,
        chimaeraSignal?.sourceChi?.informativeSites ?? stats.threeSeqSites,
        options.window,
        chimaeraSignal?.sourceChi?.halfWindow ?? null,
      ),
      statistic: stats.chimaera,
      statisticLabel: "minimum boundary χ²",
      calibration: chimaeraSignal?.sourceChi
        ? `RDP5 ChiPVal2P · ${chimaeraSignal.sourceChi.informativeSites} binary sites · half-window ${chimaeraSignal.sourceChi.halfWindow} × 3`
        : "RDP5 ChiPVal2P · binary half-window × 3",
    },
    SiScan: {
      p: sourceSiScanAvailable
        ? clampProbability(stats.siskanSourceP)
        : clampProbability(normalTwoSided(sisterZ) * windowTests),
      statistic: sisterZ,
      statisticLabel: sourceSiScanAvailable ? "Sister-Scanning category/sum Z" : "oriented category Z",
      calibration: sourceSiScanAvailable
        ? `RDP5 vertical permutation Z (${stats.siskanPValuePermutations} replicates) · ${stats.siskanOutgroupMode} outgroup · ${stats.siskanPositionMode}`
        : "fast category-Z locator fallback",
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
      const coLocated = !hasSignalLedger || locatedMethods.has(method);
      if (!coLocated) {
        return {
          method,
          pValue: 1,
          correctedP: 1,
          score: 0,
          supported: false,
          statistic: calculation.statistic,
          statisticLabel: calculation.statisticLabel,
          calibration: "no co-located discovery signal",
          correctionScope: "Not entered into the retained signal family",
        };
      }
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
        correctionScope: options.correction === "none"
          ? "Unadjusted"
          : options.correction === "holm"
            ? `Holm family across ${Math.max(1, comparisons).toLocaleString()} scanned hypotheses`
            : `Bonferroni across ${Math.max(1, comparisons).toLocaleString()} scanned triplets`,
      };
    });
}
