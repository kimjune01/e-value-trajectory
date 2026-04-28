# Experiment design

## Hindsight advantage

We know the ground truth for all conditions. This is not a prospective study. We generate data from known systems, compute e-values, and ask: does the trajectory shape match the system dynamics? Hindsight lets us validate the diagnostic before applying it to real data where ground truth is unknown.

## Prereg audit

Audited against [june.kim/prereg-audit](https://june.kim/prereg-audit). Gaps identified and addressed below.

## Assumptions that would invalidate results (Descartes Q3)

- All DGPs are Gaussian. If the trajectory diagnostic only works for normal distributions, the thesis is narrow. **Mitigation**: add one heavy-tailed condition (t-distribution, df=3) and one discrete condition (Bernoulli with time-varying p) to each claim.
- E-values are computed with λ targeting the true alternative μ=0.3. Misspecified λ could cause oscillation artifacts. **Mitigation**: run sensitivity analysis with λ at 0.5× and 2× the optimal value. If oscillation appears/disappears with λ, the diagnostic is an artifact of the betting strategy, not the system.

## Mechanism (Hume Q4)

The e-value at each step is e_t = exp(λX_t - λ²/2). When X_t is drawn from the alternative (μ > 0), the exponent is positive in expectation and E_t grows. When X_t is drawn from the null or reversed alternative (μ ≤ 0), the exponent is negative in expectation and E_t shrinks. In a cyclic system where μ(t) alternates sign, the growth and shrinkage alternate, producing oscillation in E_t with the same period as μ(t). The mechanism is multiplicative compounding of signed log-evidence.

## Competing explanations (Chamberlin Q8)

What else could cause oscillation in E_t besides system cyclicity?
1. **Misspecified λ**: too large λ amplifies noise into apparent oscillation. Tested via sensitivity analysis.
2. **Finite-sample noise**: random fluctuations in E_t that look periodic by coincidence. Tested by comparing autocorrelation structure of cyclic condition to null condition — null should show no peak.
3. **Non-stationarity in the null**: if the null DGP drifts, E_t oscillates even without a cycle. Tested by verifying condition (c) shows flat autocorrelation.

## Self-deception risk (Feynman Q14)

The shape classifier in Claim 2 uses autocorrelation at lag T/2. In real data, T is unknown. This bakes in oracle knowledge. **Mitigation**: add a "blind" variant of Claim 2 where the classifier uses a bank of lags (T/4, T/2, T, 2T) and must identify the correct lag as well as the condition. If the blind classifier fails, the diagnostic requires prior knowledge of the cycle period and is not self-sufficient.

## Severity (Mayo Q17)

"Indistinguishable trajectory shapes" requires a metric. **Metric**: two-sample Kolmogorov-Smirnov test on the distribution of rolling-window slopes between conditions (a) and (b). If KS p > 0.05, the shapes are indistinguishable and the thesis fails. For "smooth divergence," the coefficient of variation of rolling slopes must be < 0.3 (consistently positive). For "oscillation," coefficient of variation must be > 1.0 (sign-alternating).

## Power (Ioannidis Q16)

100 replications per condition. Effect size for distinguishing conditions: the KS test between rolling-slope distributions of (a) vs (b) should detect a difference with power > 0.95 at α=0.05. If preliminary runs show power < 0.80, increase to 500 replications.

## Trail (Gwern Q18)

All replications, all conditions, all parameter sensitivity runs will be published in `data/` (gitignored but reproducible via `generate_synthetic.py` with fixed seeds). The `report.md` will include every figure, every failed condition, and every sensitivity result. No curation.

Seeds: `np.random.seed(42)` for primary runs. Seeds 1–100 for replications.

---

## Pre-registration

### Claim 1: Cyclic systems produce oscillating e-value trajectories

**Setup:**
- Generate N=10,000 observations from five systems:
  - **(a) Stationary effect**: i.i.d. draws from N(μ=0.3, σ=1). True effect, no feedback.
  - **(b) Cyclic effect**: effect size oscillates sinusoidally, μ(t) = 0.3 × sin(2πt/T), period T=500.
  - **(c) Null**: i.i.d. draws from N(0, 1). No effect.
  - **(g) Heavy-tailed stationary**: i.i.d. draws from t(df=3, ncp=0.3). True effect, heavy tails.
  - **(h) Heavy-tailed cyclic**: t(df=3) with ncp(t) = 0.3 × sin(2πt/500). Cyclic, heavy tails.

**E-value computation:**
- Simple normal-mean e-value: e_t = exp(λ × X_t - λ²/2) for each observation, with λ chosen to target the alternative μ=0.3.
- Running product E_t = ∏ e_i gives the cumulative e-value trajectory.

**Prediction:**
- (a, g) E_t diverges upward. Rolling-window slope CV < 0.3.
- (b, h) E_t oscillates. Rolling-window slope CV > 1.0. Autocorrelation of log(E_t) differences peaks at lag T/2 ± 50.
- (c) E_t trends toward zero. Rolling-window slope consistently negative.

**Falsification:**
- If KS test on rolling-slope distributions of (a) vs (b) gives p > 0.05, the thesis is wrong.
- If (b) oscillates but autocorrelation peak is not at lag T/2 ± 50, the oscillation is noise, not the cycle.
- If (g) and (h) fail to replicate the pattern seen in (a) and (b), the diagnostic is Gaussian-specific.

**Sensitivity:**
- Repeat Claim 1 with λ at 0.5× and 2× optimal. Report whether oscillation in (b) appears/disappears.

**Measurement:**
- Autocorrelation of log(E_t) differences at lags 1 to 1000.
- Rolling-window slope (window = T/2 = 250) and its coefficient of variation.
- KS test between slope distributions of (a) vs (b).
- Visual: plot all five trajectories on the same axes (log scale).

### Claim 2: Trajectory shape classifies faster than threshold crossing

**Setup:**
- Same five conditions as above.
- Three classifiers race:
  - **Threshold classifier**: standard e-value test, reject null when E_t > 1/α (α=0.05, threshold = 20).
  - **Oracle shape classifier**: logistic regression on trajectory features (slope, variance, autocorrelation at lag T/2) computed over a rolling window of size W=200. Knows T.
  - **Blind shape classifier**: same features but uses a bank of lags (T/4, T/2, T, 2T) and must identify the correct lag. Does not know T.

**Prediction:**
- Oracle classifier identifies condition (b) as "cyclic" before threshold classifier reaches a verdict.
- Blind classifier identifies condition (b) as "cyclic" before threshold classifier, but later than oracle.
- Both classifiers identify condition (c) as "null" before E_t has drifted far enough to be conclusive.

**Falsification:**
- If threshold crossing is consistently faster than both shape classifiers, trajectory analysis adds no value.
- If blind classifier fails entirely, the diagnostic requires oracle knowledge of T and is not self-sufficient.

**Measurement:**
- For each condition, record the observation index at which each classifier first commits to the correct label.
- Compare median stopping times across 100 replications per condition.

### Claim 3: Regime changes in the system appear as regime changes in the evidence

**Setup:**
- Generate N=10,000 observations with a regime switch at t=5,000:
  - **(d) Effect → null**: μ=0.3 for t<5,000, μ=0 for t≥5,000.
  - **(e) Null → effect**: μ=0 for t<5,000, μ=0.3 for t≥5,000.
  - **(f) Effect → reversed**: μ=0.3 for t<5,000, μ=-0.3 for t≥5,000.

**Prediction:**
- (d) Slope of log(E_t) in rolling window flips from positive to near-zero at t≈5,000.
- (e) Slope flips from near-zero to positive.
- (f) Slope flips from positive to negative. Clearest case.

**Falsification:**
- If the slope change is not detectable within 500 observations of the true switch point (|t_detected - 5,000| > 500), the trajectory is too noisy to be diagnostic.

**Measurement:**
- CUSUM and Bayesian online changepoint detection on rolling slope of log(E_t).
- Detection delay: |t_detected - 5,000|. Target: under 500 observations.
- False alarm rate: how often does the detector fire in stationary conditions (a) and (c)?

## Implementation order

1. `generate_synthetic.py` — produce all conditions (a)–(h), with fixed seeds
2. `compute_evalues.py` — sequential e-value computation
3. `plot_trajectories.py` — visual comparison (the figure that makes or breaks the thesis)
4. `classify_shape.py` — autocorrelation, rolling slope, shape classifier (oracle + blind)
5. `detect_regime.py` — changepoint detection on evidence trajectories
6. `sensitivity.py` — λ sensitivity analysis
7. `report.md` — all results, all figures, all failures, verdict
