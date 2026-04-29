# Pre-registration: Robustness and sensitivity analysis

Fourth pre-registration. V1 tested single-stream diagnostics. V2 tested heterogeneous composition for periodicity. V3 tested four-bin classification on composed trajectories (F1 = 1.000 on clean synthetic data). V4 tests whether V3's result survives adversarial conditions.

## What V3 established

The kill-condition classifier (monotone → curvature → periodicity → aperiodic → null) achieves perfect classification on composed e-value trajectories from five heterogeneous streams sharing a forcing signal, across five labels (null, convergent, divergent, oscillatory, aperiodic). Standardized sums fail on divergence and aperiodic. Individual streams detect only oscillation.

## What V3 left open (per codex adversarial review)

1. **Manual amplitudes.** Divergence amplitude (1.5) is 15× oscillation amplitude (0.1). Difficulty is not equalized across bins.
2. **Period > 10 filter is post-hoc.** Added after observing aperiodic misclassified as oscillatory. No sensitivity analysis.
3. **Synthetic archetypes are too clean.** Pure forcing, independent streams, stationary nulls. Real data has mixed dynamics, correlated streams, and nonstationary baselines.
4. **No confidence intervals.** 100 eval reps per condition, no bootstrap CIs on F1.

## Claim 6: Robustness under adversarial conditions

**Question:** Does V3's classifier maintain F1 > 0.8 under equalized difficulty, threshold perturbation, mixed dynamics, and degraded conditions?

### Analysis 1: Amplitude sensitivity curves

For each forcing pattern, sweep amplitude from 0 to 2× the V3 value in 20 steps. At each amplitude, run 100 reps on the composed signal. Report detection rate (correct label) as a function of amplitude. Plot per-bin sensitivity curves.

**Prediction:** Each bin has a detection threshold below which classification fails. The thresholds differ across bins (divergence needs more signal than oscillation). The curves should be sigmoid-shaped with a sharp transition.

**What this reveals:** Whether V3's amplitudes are in the easy regime (plateau) or near the cliff (transition). If all four bins are on the plateau, the amplitudes were too generous.

### Analysis 2: Decision-tree ablation

Remove each decision-tree component one at a time and measure the impact on F1:

1. Remove monotone test (skip straight to curvature)
2. Remove curvature test (monotone → periodicity)
3. Remove period > 10 filter
4. Remove 0-1 test (rely on PE alone for aperiodic)
5. Remove PE (rely on 0-1 test alone for aperiodic)

**Prediction:** Curvature removal breaks convergence/divergence separation. Period filter removal breaks aperiodic/oscillatory separation. Each component is load-bearing.

### Analysis 3: Mixed dynamics

Generate conditions where two forcing patterns are superimposed:

| Condition | Forcing | Expected label |
|---|---|---|
| trend + oscillation | A_div · (t/N − 0.5) + A_osc · sin(2πt/T) | ambiguous |
| oscillation + aperiodic | A_osc · sin(2πt/T) + A_aper · z̃_t | ambiguous |
| convergence + oscillation | A_conv · exp(−t/τ) + A_osc · sin(2πt/T) | ambiguous |
| weak divergence + noise | A_div · (t/N − 0.5) + heavy-tailed noise | divergent or null |

100 reps per condition. Report the classifier's label distribution. There is no "correct" answer for mixed dynamics; the point is to see how the classifier degrades.

**Prediction:** The classifier picks the dominant pattern. When both are near-threshold, it picks whichever test fires first in the hierarchy.

### Analysis 4: Degraded conditions

Test the classifier under five realistic degradations:

1. **Autocorrelated null.** Replace i.i.d. null with AR(1) φ = 0.5 null. Does the null false positive rate increase?
2. **Correlated streams.** Add shared noise (ρ = 0.3) across all five streams. Does composition still outperform standardized sum?
3. **Missing data.** Randomly drop 10% of observations per stream. Does classification degrade?
4. **Misspecified distribution.** Replace the Normal stream with a t-distribution (df=5). Does the e-value remain valid?
5. **Nonstationary baseline.** Add a slow drift (0.01 · t/N) to the null mean. Does null false positive rate increase?

100 reps per degradation × 5 conditions. Report per-condition F1 under each degradation.

**Prediction:** Autocorrelated nulls increase false positive rate for monotone detection. Correlated streams reduce the composition advantage. Missing data degrades gracefully. Misspecified distribution has minimal impact (e-values are robust to mild misspecification). Nonstationary baseline is the hardest — it mimics divergence.

### Analysis 5: Period filter sensitivity

Sweep the period filter threshold from 2 to 50 in steps of 2. At each threshold, run 100 reps of the aperiodic condition and 100 reps of the oscillatory condition. Report:
- Aperiodic correct classification rate
- Oscillatory correct classification rate
- The threshold at which aperiodic and oscillatory cross (if they do)

**Prediction:** Below period ≈ 8, aperiodic gets misclassified as oscillatory (logistic map period-2 structure). Above period ≈ 20, some short-period oscillatory patterns get missed.

### Bootstrap confidence intervals

All F1 scores and detection rates reported with 95% bootstrap CIs (1000 bootstrap samples). Per-condition accuracy also reported with Wilson intervals. This addresses V3's missing uncertainty quantification.

### Prereg checklist audit

Audited against [the prereg checklist](https://june.kim/the-prereg-checklist). 18/20 questions answered. Two skipped:
- Q12 (Kuhn): alternative bin structures not tested — that's a different experiment.
- Q15 (Pearl): V4 makes no causal claims; it's a robustness test.

### Experiment parameters

- N = 10,000 observations per stream
- K = 5 streams (same as V3)
- 1000 null calibration reps (same seeds as V3: 80000–80999)
- 100 eval reps per condition per analysis
- V3 amplitudes as baseline: null=0, divergence=1.5, convergence=1.5, oscillation=0.1, aperiodic=0.5
- All thresholds from V3's 1000-rep null calibration

### Outputs

- `results/tables/sensitivity_amplitude.csv` — per-bin detection rate × amplitude
- `results/tables/ablation.csv` — F1 per ablation condition
- `results/tables/mixed_dynamics.csv` — label distribution per mixed condition
- `results/tables/degraded.csv` — per-condition F1 under each degradation
- `results/tables/period_filter.csv` — classification rates × period threshold
- `results/plots/sensitivity_curves.png` — amplitude sensitivity per bin
- `results/plots/ablation.png` — F1 bar chart per ablation
- `results/plots/degraded_heatmap.png` — F1 heatmap: degradation × condition

### Falsification

- If F1 drops below 0.6 under any single degradation, the classifier is fragile.
- If amplitude sensitivity shows all bins are on the plateau (>90% detection at 50% of V3 amplitude), V3's amplitudes were too generous and the result is weaker than it appears.
- If removing the period > 10 filter has no effect, it was unnecessary.
- If correlated streams reduce composed F1 to the level of standardized sum, the Fisher information advantage is an independence artifact.

### References

All V1–V3 references apply. No new citations needed — this is a stress test of existing methods.
