# Experiment design

## Hindsight advantage

We know the ground truth for all conditions. This is not a prospective study. We generate data from known systems, compute e-values, and ask: does the trajectory shape match the system dynamics? Hindsight lets us validate the diagnostic before applying it to real data where ground truth is unknown.

## Prereg audit

Audited against [june.kim/prereg-audit](https://june.kim/prereg-audit). Gaps identified and addressed inline.

## Escape hatches

**Programming errors.** If a sanity check fails (e.g., E_t under condition (c) trends positive instead of negative, or Lotka-Volterra diverges to infinity), the bug is in the implementation, not the thesis. Fix the code, re-run from scratch, log the fix in WORKLOG.md. A code fix is not a protocol change. The distinction: fixing a bug that prevents the experiment from running is maintenance; changing the DGP, the classifier, or the success criterion after seeing results is p-hacking.

**Futility.** If after implementing conditions (a)–(c), the periodogram of condition (b) shows no detectable peak above the noise floor of condition (c), stop. The mechanism doesn't work empirically even though the math says it should. Diagnose: is N too small? Is λ too far from the current μ(t) for too long? Is the Lotka-Volterra period too long relative to N? Log the diagnosis in WORKLOG.md and either (1) adjust a single parameter with justification and re-run (logged as a new experiment, not a revision of the original), or (2) report the negative result.

**Scope reduction.** If Claim 1 fails, Claims 2 and 3 are moot. Don't run them. Report Claim 1's failure and the diagnosis. Partial results are results.

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

**Question:** Does the periodogram of an e-value trajectory show a spectral peak at the system's cycle period, distinguishable from noise, autocorrelation, and random regime switching?
**In plain English:** If the thing you're studying goes in circles, can you see the circles in your evidence?

**N = 10,000 observations for all conditions.**

**Setup — five DGPs:**
- **(a) Stationary effect**: i.i.d. draws from N(μ=0.3, σ=1). True effect, no dynamics.
- **(b) Deterministic feedback cycle (Lotka-Volterra)**: discretized predator-prey system (Euler, dt=0.01). Parameters: α=1.0, β=0.1, δ=0.075, γ=1.5, initial (x₀, y₀) = (10, 5). X_t is noisy observation of prey population: X_t = x_t + ε, ε ~ N(0, 1). Expected emergent period ≈ 500 steps. If actual period < 100 or > 2000, adjust parameters per escape hatch and log as new experiment.
- **(b2) Stochastic feedback cycle (stochastic Lotka-Volterra)**: same parameters as (b), but with process noise added to the dynamics: dx = (αx - βxy)dt + σ_x · x · dW₁, dy = (δxy - γy)dt + σ_y · y · dW₂. Process noise σ_x = σ_y = 0.05. The period, amplitude, and phase all wander. X_t = x_t + ε, ε ~ N(0, 1). This is the fair test: can the diagnostic detect a noisy cycle where the period itself is stochastic?
- **(c) Null**: i.i.d. draws from N(0, 1). No effect.
- **(d) AR(1) with drift**: X_t = 0.9·X_{t-1} + 0.3 + ε_t, ε ~ N(0,1). Autocorrelated, mean-reverting, no explicit cycle. Tests whether autocorrelation alone produces oscillation artifacts.
- **(e) Regime-switching**: hidden Markov model with two states (μ=0.3 and μ=-0.3), transition probability 0.005 per step. Random regime changes, not periodic.

**E-value computation:**
- e_t = exp(λ × X_t - λ²/2), λ = 0.3 (targeting N(0.3, 1) alternative).
- Running product E_t = ∏ e_i.

**Prediction:**
- (a) log(E_t) grows linearly. Periodogram of Δlog(E_t) shows no peak.
- (b) log(E_t) oscillates. Periodogram of Δlog(E_t) shows a sharp peak at the emergent Lotka-Volterra period.
- (b2) log(E_t) oscillates with a broader, noisier peak in the periodogram. Peak should still be detectable but wider than (b)'s. If the peak vanishes entirely, the diagnostic fails on stochastic cycles.
- (c) log(E_t) trends negative. Periodogram flat.
- (d) log(E_t) grows (positive mean) but periodogram shows no peak (autocorrelation ≠ periodicity).
- (e) log(E_t) shows regime-dependent growth/shrinkage. Periodogram shows no clean peak (random switching, not periodic).

**Falsification:**
- If (b) and (c) produce indistinguishable periodograms (TOST equivalence on peak height, margin d = 0.2), the thesis is wrong even for deterministic cycles.
- If (b) passes but (b2) is indistinguishable from (c), the diagnostic only works on deterministic systems and is not useful for real feedback systems where the period wanders.
- If (d) produces a periodogram peak comparable to (b), autocorrelation confounds the diagnostic.
- If (e) produces a clean peak, the diagnostic can't distinguish periodic from random switching.

**Sensitivity:**
- Repeat with λ at 0.5× and 2× optimal. Report whether oscillation pattern changes.

**Measurement:**
- Periodogram of Δlog(E_t) for all five conditions (no oracle parameters needed).
- Same periodogram on raw X_t, for comparison. The thesis holds if the e-value periodogram reveals structure the raw-data periodogram also shows (validating transparency) or reveals it more cleanly (e-value as amplifier).
- TOST equivalence test between periodogram peak heights of (b) vs (c) noise floor.
- Visual: plot all five trajectories on the same axes (log scale).

### Claim 2: Trajectory shape classifies faster than threshold crossing

**Question:** Does a spectral classifier on the e-value trajectory identify cyclic systems earlier than a standard e-value threshold test, and does e-value's variable-method flexibility outperform a fixed-model Bayesian sequential test?
**In plain English:** Does the shape of your bet tell you something that the size of your bet doesn't?

**Setup:**
- Same five conditions, N=10,000.
- All three classifiers output the same three labels: **reject null** (effect detected), **periodic** (cyclic dynamics detected), **null** (no effect).
  - **Threshold classifier**: reject null when E_t > 20 (α=0.05). Labels "null" when log(E_t) slope is negative over a 500-observation window. Labels "periodic" never (no mechanism to detect periodicity). Baseline: can't diagnose cyclicity.
  - **Bayesian sequential classifier**: online Bayesian posterior under a fixed N(0.3, 1) likelihood. Labels "reject null" when posterior P(μ > 0) > 0.95. Labels "null" when P(μ > 0) < 0.05. Labels "periodic" via sliding-window periodogram on the posterior mean trajectory (same spectral method as below, but on the posterior, not the e-value). Baseline: fixed model, can detect periodicity through the posterior.
  - **Spectral e-value classifier**: sliding-window periodogram on Δlog(E_t) (window = 500, no knowledge of T). Labels "periodic" when peak power exceeds 3× median power. Labels "reject null" when slope of log(E_t) is consistently positive over a 500-observation window with no spectral peak. Labels "null" when slope is consistently negative.

**Prediction:**
- For condition (b), spectral e-value classifier labels "periodic" before threshold classifier labels "reject null," and at least as early as the Bayesian classifier labels "periodic."
- For condition (c), all three classifiers reach "null" at similar times.
- For conditions (a) and (d), threshold classifier labels "reject null" first (no periodicity to detect).

**Falsification:**
- If Bayesian classifier detects periodicity consistently earlier than the spectral e-value classifier, e-values add no advantage over standard Bayesian sequential testing for this task.
- If spectral e-value classifier never labels (b) as "periodic" across 100 replications, it can't detect the cyclicity the thesis claims is visible.

**Measurement:**
- For each condition, record the observation index at which each classifier first commits to any label.
- Compare median stopping times across 100 replications. Report interquartile range and variance.

### Claim 3: Regime changes in the system appear as regime changes in the evidence

**Question:** When the data-generating process switches regime mid-stream, does the slope of log(E_t) change detectably within 500 observations of the true switch point?
**In plain English:** If the rules change halfway through, does your evidence notice?

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
- CUSUM on sliding-window slope of log(E_t). Reference value k = 0.5 × (mean slope under effect). Threshold h = 5 (standard). These are pinned; no post-hoc tuning.
- Bayesian online changepoint detection (Adams & MacKay 2007) with hazard rate = 1/2000 (expecting ~1 change per 2000 observations). Pinned.
- Detection delay and false alarm rate (tested against stationary conditions a and c).
- Head-to-head: same detectors, same parameters, on raw X_t vs on log(E_t).

**Multiple comparisons:** Across all claims, 5 DGP conditions × 3 regime variants × 2 detection methods × 3 λ values = ~90 tests. Apply Bonferroni correction (α = 0.05/90 ≈ 0.0006) to all significance tests. Report both corrected and uncorrected p-values.

## Implementation order

1. `generate_synthetic.py` — produce all conditions (a)–(h), with fixed seeds
2. `compute_evalues.py` — sequential e-value computation
3. `plot_trajectories.py` — visual comparison
4. `spectral_analysis.py` — periodogram on Δlog(E_t) and raw X_t
5. `classify_shape.py` — spectral classifier vs threshold classifier race
6. `detect_regime.py` — changepoint detection on evidence trajectories vs raw data
7. `sensitivity.py` — λ sensitivity analysis
8. `report.md` — all results, all figures, all failures, verdict
