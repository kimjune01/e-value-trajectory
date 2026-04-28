# Experiment design

## Hindsight advantage

We know the ground truth for all conditions. This is not a prospective study. We generate data from known systems, compute e-values, and ask: does the trajectory shape match the system dynamics? Hindsight lets us validate the diagnostic before applying it to real data where ground truth is unknown.

## Pre-registration

### Claim 1: Cyclic systems produce oscillating e-value trajectories

**Setup:**
- Generate N=10,000 observations from three systems:
  - **(a) Stationary effect**: i.i.d. draws from N(μ=0.3, σ=1). True effect, no feedback.
  - **(b) Cyclic effect**: effect size oscillates sinusoidally, μ(t) = 0.3 × sin(2πt/T), period T=500. Feedback system where the effect reverses periodically.
  - **(c) Null**: i.i.d. draws from N(0, 1). No effect.

**E-value computation:**
- Use the simple normal-mean e-value: e_t = exp(λ × X_t - λ²/2) for each observation, with λ chosen to target the alternative μ=0.3.
- Running product E_t = ∏ e_i gives the cumulative e-value trajectory.

**Prediction:**
- (a) E_t diverges upward (smooth exponential growth on log scale)
- (b) E_t oscillates: grows during positive-effect phases, shrinks during negative-effect phases. The period of the oscillation in E_t should match the period T of the underlying cycle.
- (c) E_t trends toward zero (supermartingale convergence under null)

**Falsification:**
- If (a) and (b) produce indistinguishable trajectory shapes, the thesis is wrong. The e-value is not transparent to system dynamics.
- If (b) oscillates but with a period unrelated to T, the e-value reflects noise, not the cycle.

**Measurement:**
- Autocorrelation of log(E_t) differences. Cyclic system should show a peak at lag T/2. Stationary system should show no autocorrelation structure.
- Slope of log(E_t) in rolling windows of size T/2. Cyclic system should show sign-alternating slopes. Stationary system should show consistently positive slope.
- Visual: plot all three trajectories on the same axes (log scale).

### Claim 2: Trajectory shape classifies faster than threshold crossing

**Setup:**
- Same three conditions as above.
- Two classifiers race:
  - **Threshold classifier**: standard e-value test, reject null when E_t > 1/α (α=0.05, so threshold = 20).
  - **Shape classifier**: logistic regression on trajectory features (slope, variance, autocorrelation at lag T/2) computed over a rolling window of size W=200.

**Prediction:**
- Shape classifier identifies condition (b) as "cyclic" before threshold classifier rejects or fails to reject.
- Shape classifier identifies condition (c) as "null" before E_t has drifted far enough to be conclusive.

**Falsification:**
- If threshold crossing is consistently faster than shape classification, trajectory analysis adds no value.

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
- (d) E_t grows then plateaus. Slope of log(E_t) in rolling window flips from positive to near-zero at t≈5,000.
- (e) E_t stalls then grows. Slope flips from near-zero to positive.
- (f) E_t grows then shrinks. Slope flips from positive to negative. This is the clearest case: the evidence trajectory reverses, signaling the regime change.

**Falsification:**
- If the slope change is not detectable within T/4 observations of the true switch point, the trajectory is too noisy to be diagnostic.

**Measurement:**
- CUSUM or Bayesian changepoint detection on log(E_t) slope. Compare detected changepoint to true changepoint at t=5,000.
- Detection delay: |t_detected - 5,000|. Target: under 500 observations.

## Implementation order

1. `generate_synthetic.py` — produce all conditions (a)–(f)
2. `compute_evalues.py` — sequential e-value computation
3. `plot_trajectories.py` — visual comparison (the figure that makes or breaks the thesis)
4. `classify_shape.py` — autocorrelation, rolling slope, shape classifier
5. `detect_regime.py` — changepoint detection on evidence trajectories
6. `report.md` — results, figures, verdict
