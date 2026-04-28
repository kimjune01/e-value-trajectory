# Experiment design

## Hindsight advantage

We know the ground truth for all conditions. This is not a prospective study. We generate data from known systems, compute e-values, and ask: does the trajectory shape match the system dynamics? Hindsight lets us validate the diagnostic before applying it to real data where ground truth is unknown.

## Prereg audit

Audited against [june.kim/prereg-audit](https://june.kim/prereg-audit). Gaps identified and addressed inline.

## The mechanism (Hume Q4)

The e-value at each step is e_t = exp(λX_t - λ²/2). The expected log-evidence per step is:

    E[log e_t] = λμ(t) - λ²/2

When μ(t) is constant (stationary effect), E[log e_t] is constant and positive (if λ matches the alternative). E_t grows exponentially. When μ(t) oscillates (cyclic system), E[log e_t] oscillates with the same period. E_t grows during phases where λμ(t) > λ²/2 and shrinks otherwise. The oscillation in the evidence is a direct algebraic consequence of the oscillation in the data-generating process, not a statistical artifact.

When λ is misspecified relative to the current μ(t) (which happens for half of each cycle), the e-value is still valid (supermartingale under the null) but the trajectory shape reflects the misspecification. This is a feature: the shape tells you the betting strategy doesn't match the system's current regime.

## Assumptions that would invalidate results (Descartes Q3)

- All DGPs are univariate. Multivariate or high-dimensional systems are out of scope for this experiment.
- E-values are computed with fixed λ. Adaptive λ (which would track the system) is a different experiment.

## Competing explanations (Chamberlin Q8)

What else could cause oscillation in E_t besides system cyclicity?
1. **Misspecified λ**: too large λ amplifies noise into apparent oscillation. Tested via sensitivity analysis.
2. **Finite-sample noise**: random fluctuations that look periodic by coincidence. Tested by comparing spectral analysis of cyclic condition to null condition.
3. **Autocorrelation in the DGP**: an AR process can produce quasi-periodic trajectories without explicit cycles. This is why Claim 1 includes an AR condition.

## Self-deception risk (Feynman Q14)

Previous draft baked oracle knowledge of the cycle period T into every analysis component. Fixed: all periodicity detection now uses spectral methods (periodogram) that don't require knowing T. The hypothesis is that the e-value trajectory's periodogram shows a peak at the system's true period, not that a pre-tuned detector finds what it was tuned to find.

## Severity (Mayo Q17)

"Indistinguishable" is tested via TOST (two one-sided tests) equivalence testing on rolling-slope distributions, not KS accept-the-null. Equivalence margin: Cohen's d = 0.2 (small effect).

## Trail (Gwern Q18)

All replications, all conditions, all sensitivity runs published via reproducible scripts with fixed seeds. `report.md` will include every figure and every failure. No curation. Seeds: `np.random.seed(42)` for primary runs. Seeds 1–100 for replications.

---

## Pre-registration

### Claim 1: Cyclic systems produce oscillating e-value trajectories

**Setup — five DGPs, three with feedback:**
- **(a) Stationary effect**: i.i.d. draws from N(μ=0.3, σ=1). True effect, no dynamics.
- **(b) Feedback cycle (Lotka-Volterra)**: discretized predator-prey system. X_t is noisy observation of prey population. The oscillation emerges from the dynamics (prey grows → predator grows → prey shrinks → predator shrinks → repeat). Period is not injected; it emerges from the parameters.
- **(c) Null**: i.i.d. draws from N(0, 1). No effect.
- **(d) AR(1) with drift**: X_t = 0.9·X_{t-1} + 0.3 + ε_t, ε ~ N(0,1). Autocorrelated, mean-reverting, no explicit cycle. Tests whether autocorrelation alone produces oscillation artifacts.
- **(e) Regime-switching**: hidden Markov model with two states (μ=0.3 and μ=-0.3), transition probability 0.005 per step. Random regime changes, not periodic.

**E-value computation:**
- e_t = exp(λ × X_t - λ²/2), λ = 0.3 (targeting N(0.3, 1) alternative).
- Running product E_t = ∏ e_i.

**Prediction:**
- (a) log(E_t) grows linearly. Periodogram of Δlog(E_t) shows no peak.
- (b) log(E_t) oscillates. Periodogram of Δlog(E_t) shows a peak at the emergent Lotka-Volterra period.
- (c) log(E_t) trends negative. Periodogram flat.
- (d) log(E_t) grows (positive mean) but periodogram shows no peak (autocorrelation ≠ periodicity).
- (e) log(E_t) shows regime-dependent growth/shrinkage. Periodogram shows no clean peak (random switching, not periodic).

**Falsification:**
- If (a) and (b) produce indistinguishable periodograms (TOST equivalence, d < 0.2 on peak height), the thesis is wrong.
- If (d) produces a periodogram peak comparable to (b), autocorrelation confounds the diagnostic.
- If (e) produces a clean peak, the diagnostic can't distinguish periodic from random switching.

**Sensitivity:**
- Repeat with λ at 0.5× and 2× optimal. Report whether oscillation pattern changes.

**Measurement:**
- Periodogram of Δlog(E_t) for all five conditions (no oracle parameters needed).
- Same periodogram on raw X_t, for comparison. The thesis holds if the e-value periodogram reveals structure the raw-data periodogram also shows (validating transparency) or reveals it more cleanly (e-value as amplifier).
- TOST equivalence test between periodogram peak heights of (a) vs (b).
- Visual: plot all five trajectories on the same axes (log scale).

### Claim 2: Trajectory shape classifies faster than threshold crossing

**Setup:**
- Same five conditions.
- Two classifiers race:
  - **Threshold classifier**: standard e-value test, reject null when E_t > 20 (α=0.05).
  - **Spectral classifier**: sliding-window periodogram on Δlog(E_t) (window = 500, no knowledge of T). Classifies as "periodic" when peak power exceeds 3× median power; "stationary" when slope of log(E_t) is consistently positive; "null" when slope is consistently negative.

**Prediction:**
- Spectral classifier identifies (b) as periodic before threshold classifier reaches a verdict.
- Spectral classifier identifies (c) as null before E_t has drifted far enough to be conclusive.

**Falsification:**
- If threshold crossing is consistently faster than spectral classification, trajectory analysis adds no value over standard e-value testing.

**Measurement:**
- For each condition, record the observation index at which each classifier first commits to the correct label.
- Compare median stopping times across 100 replications. Report interquartile range.

### Claim 3: Regime changes in the system appear as regime changes in the evidence

**Setup:**
- Generate N=10,000 observations with a regime switch at t=5,000:
  - **(f) Effect → null**: μ=0.3 for t<5,000, μ=0 for t≥5,000.
  - **(g) Null → effect**: μ=0 for t<5,000, μ=0.3 for t≥5,000.
  - **(h) Effect → reversed**: μ=0.3 for t<5,000, μ=-0.3 for t≥5,000.

**Prediction:**
- (f) Slope of log(E_t) in sliding window flips from positive to near-zero at t≈5,000.
- (g) Slope flips from near-zero to positive.
- (h) Slope flips from positive to negative.

**Falsification:**
- If changepoint detection delay |t_detected - 5,000| > 500 observations, the trajectory is too noisy to be diagnostic.

**Measurement:**
- CUSUM and Bayesian online changepoint detection on sliding-window slope of log(E_t).
- Detection delay and false alarm rate (tested against stationary conditions a and c).
- Head-to-head: same changepoint detectors on raw X_t vs on log(E_t). Does the e-value trajectory make the regime change easier or harder to detect than the raw data?

## Implementation order

1. `generate_synthetic.py` — produce all conditions (a)–(h), with fixed seeds
2. `compute_evalues.py` — sequential e-value computation
3. `plot_trajectories.py` — visual comparison
4. `spectral_analysis.py` — periodogram on Δlog(E_t) and raw X_t
5. `classify_shape.py` — spectral classifier vs threshold classifier race
6. `detect_regime.py` — changepoint detection on evidence trajectories vs raw data
7. `sensitivity.py` — λ sensitivity analysis
8. `report.md` — all results, all figures, all failures, verdict
