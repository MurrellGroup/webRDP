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
// walk with a fixed number of up/down steps. This is the same probability
// represented by RDP5 Seq3PVals/Get3SeqPvalC, evaluated as a first-passage DP
// rather than materialising the desktop program's four-dimensional YTable.
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

function threeSeqApproxNormCdf(value) {
  const b1 = 0.31938153;
  const b2 = -0.356563782;
  const b3 = 1.781477937;
  const b4 = -1.821255978;
  const b5 = 1.330274429;
  const p = 0.2316419;
  const c = 0.39894228;
  const x = Number(value);
  if (x > 8) return 1;
  if (x < -8) return 0;
  if (x >= 0) {
    const t = 1 / (1 + p * x);
    return 1 - c * Math.exp(-x * x / 2) * t
      * (t * (t * (t * (t * b5 + b4) + b3) + b2) + b1);
  }
  const t = 1 / (1 - p * x);
  return c * Math.exp(-x * x / 2) * t
    * (t * (t * (t * (t * b5 + b4) + b3) + b2) + b1);
}

// Direct port of RDP5 SiegmundDiscrete and ApproxNu. GetTSPVal uses this when
// a tuple is larger than the installed exact probability table.
export function threeSeqSiegmundP(upSteps, downSteps, observedDescent) {
  const up = Math.max(0, Math.trunc(upSteps));
  const down = Math.max(0, Math.trunc(downSteps));
  const threshold = Math.max(0, Math.trunc(observedDescent));
  const total = up + down;
  if (threshold <= 0 || up === 0 || down === 0 || total === 0) return 1;
  const boundary = threshold - 0.5;
  const displacement = down - up;
  const exponent = Math.min(700, -2 * boundary * (boundary - displacement) / total);
  const p1 = Math.min(1e200, Math.exp(exponent));
  const p2 = p1 * (2 * (2 * boundary - displacement) * (boundary - displacement) / total + 1);
  const a = 2 * (2 * boundary - displacement) / total;
  if (!(a > 0) || !(p2 > 0) || !Number.isFinite(p2)) return null;
  const half = a / 2;
  const cdf = threeSeqApproxNormCdf(half);
  const density = Math.exp(-0.5 * half * half) / Math.sqrt(2 * Math.PI);
  const denominator = a * (density + a * cdf / 2);
  if (!(denominator > 0)) return null;
  const nu = ((cdf - 0.5) * 2) / denominator;
  const poisson = nu * nu * p2;
  if (!(poisson > 0) || !Number.isFinite(poisson)) return null;
  const probability = -Math.expm1(-poisson);
  return probability > 0 && probability <= 1 ? probability : null;
}

// Source dispatch used by production 3Seq. Small walks use the exact
// Seq3PVals-equivalent DP. Large walks follow GetTSPVal's SiegmundDiscrete
// branch; its final scaled-table fallback is retained for approximation edge
// cases where the discrete expression is outside its numerical domain.
export function threeSeqSourceP(upSteps, downSteps, observedDescent, maxOperations = 4_000_000) {
  const up = Math.max(0, Math.trunc(upSteps));
  const down = Math.max(0, Math.trunc(downSteps));
  const threshold = Math.max(0, Math.trunc(observedDescent));
  if (threshold <= 0 || up === 0 || down === 0) return { p: 1, mode: "exact-table" };
  const exact = threeSeqExactP(up, down, threshold, maxOperations);
  if (exact.exact) return { p: exact.p, mode: "exact-table" };
  const siegmund = threeSeqSiegmundP(up, down, threshold);
  if (siegmund !== null) return { p: clampProbability(siegmund), mode: "siegmund-discrete" };

  const tableLimit = Math.max(8, Math.floor(Math.cbrt(Math.max(512, maxOperations))) - 2);
  const factor = Math.max(up, down, threshold) / tableLimit;
  if (!(factor > 1)) return { p: 1, mode: "unavailable" };
  const scaledUp = Math.max(0, Math.floor(up / factor));
  const scaledDown = Math.max(0, Math.floor(down / factor));
  const scaledThreshold = Math.max(1, Math.floor(threshold / factor));
  const scaled = threeSeqExactP(scaledUp, scaledDown, scaledThreshold, maxOperations);
  if (!scaled.exact || scaled.p === null) return { p: 1, mode: "unavailable" };
  return { p: clampProbability(Math.pow(scaled.p, factor)), mode: "scaled-table" };
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
  const sourceRdpP = rdp5SourceProbability(stats.rdpSource
    ? { ...stats.rdpSource, informativeSites: candidate.informative }
    : null);
  const hasSignalLedger = Array.isArray(candidate.methodSignals);
  const locatedMethods = new Set((candidate.methodSignals ?? []).map((signal) => signal.method));
  const geneconvSignal = (candidate.methodSignals ?? []).find((signal) => signal.method === "GENECONV" && signal.sourceGeneconv);
  const bootscanSignal = (candidate.methodSignals ?? []).find((signal) => signal.method === "BootScan" && signal.sourceBootscan);
  const maxChiSignal = (candidate.methodSignals ?? []).find((signal) => signal.method === "MaxChi" && signal.sourceChi);
  const chimaeraSignal = (candidate.methodSignals ?? []).find((signal) => signal.method === "Chimaera" && signal.sourceChi);
  const threeSeqSignal = (candidate.methodSignals ?? []).find((signal) => signal.method === "3Seq" && signal.sourceThreeSeq);

  const calculations = {
    RDP: {
      p: sourceRdpP ?? clampProbability(binomialUpper(candidate.insideMinor, insideTotal, backgroundMinor) * windowTests),
      statistic: candidate.insideMinor / Math.max(1, insideTotal) - backgroundMinor,
      statisticLabel: "identity shift",
      calibration: sourceRdpP === null ? "binomial · window-corrected fallback" : "RDP5 ProbCalcP/P2-equivalent binomial tail",
    },
    GENECONV: {
      p: geneconvSignal?.sourceGeneconv?.rawP
        ?? geneconvSourceProbability(stats.genconvRun, stats.genconvEligible, stats.genconvMatches, options.geneconvGScale ?? 1),
      statistic: geneconvSignal?.sourceGeneconv?.fragmentScore ?? stats.genconvRun,
      statisticLabel: "fragment score",
      calibration: geneconvSignal?.sourceGeneconv
        ? `RDP5 six-track fragment queue · ${geneconvSignal.sourceGeneconv.informativeSites} compressed sites · G-scale ${Math.round(options.geneconvGScale ?? 1)}`
        : `RDP5 CalcKMaxP/GCCalcPValP · G-scale ${Math.round(options.geneconvGScale ?? 1)}`,
    },
    BootScan: {
      p: bootscanSignal?.sourceBootscan?.rawP ?? 1,
      statistic: bootscanSignal?.sourceBootscan?.bootstrapSupport ?? 0,
      statisticLabel: bootscanSignal?.sourceBootscan
        ? "maximum topology bootstrap"
        : "source batch signal unavailable",
      calibration: bootscanSignal?.sourceBootscan
        ? `RDP5 RecScan ${bootscanSignal.sourceBootscan.relationshipMode} batch · ${bootscanSignal.sourceBootscan.bootstrapReplicates} SEQBOOT2 replicates · ${bootscanSignal.sourceBootscan.window}/${bootscanSignal.sourceBootscan.step} nt window/step · MakeScoresBS/ProbCalc`
        : "RDP5 source BootScan batch signal unavailable · recalculate this hypothesis",
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
      p: threeSeqSignal?.sourceThreeSeq?.rawP ?? 1,
      statistic: threeSeqSignal?.sourceThreeSeq?.descent ?? descent,
      statisticLabel: threeSeqSignal?.sourceThreeSeq ? "maximum HGRW excursion" : "source signal unavailable",
      calibration: threeSeqSignal?.sourceThreeSeq
        ? `RDP5 FindSubSeqTS/Seq3PVals · ${threeSeqSignal.sourceThreeSeq.informativeSites} compressed sites · ${threeSeqSignal.sourceThreeSeq.probabilityMode}`
        : "RDP5 source 3Seq signal unavailable · recalculate this hypothesis",
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
