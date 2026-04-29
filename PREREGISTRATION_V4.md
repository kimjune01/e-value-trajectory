# Pre-registration: Robustness and sensitivity analysis

Fourth pre-registration. V3 achieved F1 = 1.000 on clean synthetic data. V4 tests whether V3's result survives adversarial conditions.

## Primary endpoint

- **Primary metric:** macro-F1 over five clean-label classes (null, convergent, divergent, oscillatory, aperiodic).
- **Primary robustness rule:** pass if lower 95% bootstrap CI of macro-F1 > 0.8 at the mildest severity of each degradation.
- **Stress-test rule:** report macro-F1 at all severity levels. Harsher severities are expected to degrade and are descriptive, not pass/fail.
- **Failure rule:** fail if macro-F1 point estimate < 0.6 at the mildest severity of any degradation.
- **Secondary metrics:** per-class recall, null false positive rate, full confusion matrix, label entropy for mixed cases.
- **Unit of resampling:** replicate trajectory, not time point.

## Seed policy

- Calibration seeds: 80000–80999 (1000 reps, same as V3).
- Evaluation seeds: 99999+ per condition, disjoint from calibration.
- Ablations use identical generated trajectories (same seeds) across original and ablated classifiers.
- No threshold retuning after seeing V4 results.

## Analysis 1: Amplitude sensitivity curves

For each forcing pattern, sweep amplitude from 0 to 2× the V3 value in 20 steps. At each amplitude, run 100 reps on the composed signal. Report detection rate (correct label) as a function of amplitude.

**Second stage: equalized difficulty.** From the sweep, estimate the amplitude at which each bin reaches 80% detection rate on composed signal. Rerun the classifier at those matched amplitudes (100 reps per bin + 100 null reps). Report five-class macro-F1 at equalized difficulty.

**Prediction:** Detection thresholds differ across bins. Divergence needs more signal than oscillation. At equalized difficulty, macro-F1 will be lower than V3's 1.000 but above 0.8.

## Analysis 2: Decision-tree ablation

Remove each component one at a time. Use identical seeds and trajectories across all variants. Report paired deltas in macro-F1 and per-class recall.

1. Remove monotone test (skip to curvature)
2. Remove curvature test (monotone → periodicity)
3. Remove period > 10 filter
4. Remove 0-1 test (PE alone for aperiodic)
5. Remove PE (0-1 test alone for aperiodic)
6. Remove Mann-Kendall threshold (use only Theil-Sen slope)
7. Perturb curvature R_curve threshold by ±20%

**Prediction:** Curvature removal breaks convergence/divergence. Period filter removal breaks aperiodic/oscillatory. Each component is load-bearing. Threshold perturbation ±20% degrades but doesn't collapse.

## Analysis 3: Mixed dynamics (descriptive only)

Generate conditions with superimposed forcing. No correct label exists. Do not fold into global F1. Report label distributions, hierarchy-trigger frequencies, and dominance ratios.

| Condition | Forcing |
|---|---|
| trend + oscillation | A_div · (t/N − 0.5) + A_osc · sin(2πt/T) |
| oscillation + aperiodic | A_osc · sin(2πt/T) + A_aper · z̃_t |
| convergence + oscillation | A_conv · exp(−t/τ) + A_osc · sin(2πt/T) |
| weak divergence + heavy-tailed noise | A_div · (t/N − 0.5) + t-distributed noise (df=3) |

Dominance ratio: vary the relative amplitudes of the two forcings at ratios 0.2:0.8, 0.4:0.6, 0.5:0.5, 0.6:0.4, 0.8:0.2. Report which label the classifier picks at each ratio.

100 reps per condition per ratio. Report label distribution and entropy.

**Prediction:** The classifier picks whichever pattern's test fires first in the hierarchy. At 0.5:0.5, the earlier test in the hierarchy wins.

## Analysis 4: Degraded conditions (severity sweeps)

Test the classifier under degradations at multiple severity levels. Each degradation is applied to all five clean-label conditions (null, convergent, divergent, oscillatory, aperiodic) so macro-F1 is computable. 500 reps per degradation level × 5 conditions.

| Degradation | Severity levels |
|---|---|
| Autocorrelated null (AR(1)) | φ = {0.2, 0.5, 0.8} |
| Correlated streams (shared noise) | ρ = {0.1, 0.3, 0.6} |
| Missing data (random dropout) | {5%, 10%, 25%, 50%} |
| Misspecified distribution (Normal → t) | df = {3, 5, 10} |
| Nonstationary baseline (slow drift) | drift = {0.005, 0.01, 0.02} · t/N |

Report per-condition macro-F1 and null false positive rate at each severity level. Include paired comparisons against standardized sum baseline at every severity level.

**Predictions:**
- Autocorrelated nulls: false positive rate for monotone detection increases with φ. At φ = 0.8, null FPR > 10%.
- Correlated streams: composition advantage shrinks with ρ. At ρ = 0.6, composed F1 approaches standardized sum F1.
- Missing data: graceful degradation up to 25%. At 50%, classification degrades substantially (stress-test level, not primary pass/fail).
- Misspecified distribution: minimal impact at df = 10 and df = 5. At df = 3, heavy tails inflate feature estimates.
- Nonstationary baseline: mimics divergence. At drift = 0.02, null FPR for divergence > 20%.

## Analysis 5: Period filter sensitivity

Sweep the period filter threshold from 2 to 50 in steps of 2. At each threshold, run 100 reps of aperiodic and 100 reps of oscillatory. Report classification rates for both. Identify the crossover point.

**Prediction:** Below period ≈ 8, aperiodic misclassified as oscillatory. Above period ≈ 20, short-period oscillatory patterns missed.

## E-value validity under misspecification

For the misspecified distribution condition (Normal → t), check e-value validity operationally:
- Empirical type-I rate: under null with t-distributed Normal stream, does max_t E_t exceed 20 (anytime exceedance) in ≤5% of reps?
- Supermartingale diagnostic: compute mean of e_t at each time point across 1000 reps. Report max over time of this mean. If max mean > 1.05, the e-value is inflated.

This is separate from classifier F1. The e-value can be invalid (inflated type-I) while the classifier still works, or vice versa.

## Experiment parameters

- N = 10,000 observations per stream
- K = 5 streams (same as V3)
- 1000 null calibration reps (seeds 80000–80999)
- Analysis 1: 100 reps × 20 amplitude steps × 4 bins + 100 reps at equalized amplitudes
- Analysis 2: 100 reps × 7 ablation variants, paired seeds
- Analysis 3: 100 reps × 4 conditions × 4 dominance ratios
- Analysis 4: 500 reps × 5 conditions × 3–4 severity levels per degradation
- Analysis 5: 100 reps × 25 threshold values × 2 conditions
- Bootstrap CIs: 1000 bootstrap samples for all F1 estimates
- All thresholds from V3's 1000-rep null calibration, frozen

## Outputs

- `results/tables/v4_amplitude.csv`
- `results/tables/v4_ablation.csv`
- `results/tables/v4_mixed.csv`
- `results/tables/v4_degraded.csv`
- `results/tables/v4_period_filter.csv`
- `results/tables/v4_validity.csv`
- `results/plots/v4_sensitivity_curves.png`
- `results/plots/v4_ablation.png`
- `results/plots/v4_degraded_heatmap.png`
- `results/plots/v4_mixed_dominance.png`

## Falsification

- If lower 95% CI of macro-F1 < 0.8 under any degradation at its mildest severity, the classifier is fragile.
- If amplitude sensitivity shows >90% detection at 50% of V3 amplitude for all bins, V3 was too easy.
- If removing the period > 10 filter has no effect on aperiodic classification, it was unnecessary.
- If correlated streams at ρ = 0.3 reduce composed F1 to standardized sum level, Fisher information advantage is an independence artifact.
- If e-value type-I rate exceeds 10% under t(df=3) misspecification, the Normal-assumption e-value is invalid for heavy-tailed data.
