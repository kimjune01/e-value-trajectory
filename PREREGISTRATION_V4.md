# Pre-registration: Robustness analysis

V3 achieved F1 = 1.000 on clean synthetic data. V4 asks: how fragile is it?

## Primary endpoint

- **Metric:** macro-F1 over five classes (null, convergent, divergent, oscillatory, aperiodic).
- **Pass:** lower 95% bootstrap CI > 0.8 at mildest severity of each perturbation.
- **Fail:** point estimate < 0.6 at mildest severity.
- **Harsher severities are descriptive, not pass/fail.**
- Bootstrap: 1000 samples. Unit of resampling: trajectory.

## Seed policy

Calibration: 80000–80999 (1000 reps). Evaluation: 99999+, disjoint. Ablations use identical trajectories (paired seeds). No threshold retuning after V4 starts.

## Design

One pipeline, many perturbations. Each perturbation modifies one aspect of V3's setup, runs the same generate → compose → classify → score loop, and reports macro-F1 + confusion matrix + paired comparison against standardized sum.

### Perturbation grid

| Category | Perturbation | Severity levels | Reps |
|---|---|---|---|
| **Amplitude** | Sweep each bin 0 to 2× V3 | 20 steps | 100 |
| **Amplitude** | Equalized difficulty (80% detection per bin + null) | 1 | 100 |
| **Ablation** | Remove curvature test | 1 | 100 |
| **Ablation** | Remove period > 10 filter | 1 | 100 |
| **Ablation** | Remove 0-1 test (PE only) | 1 | 100 |
| **Ablation** | Remove PE (0-1 only) | 1 | 100 |
| **Ablation** | Remove Mann-Kendall (slope only) | 1 | 100 |
| **Ablation** | Curvature threshold ±20% | 2 | 100 |
| **Degradation** | Autocorrelated null (AR(1)) | φ = {0.2, 0.5, 0.8} | 500 |
| **Degradation** | Correlated streams | ρ = {0.1, 0.3, 0.6} | 500 |
| **Degradation** | Missing data | {5%, 10%, 25%} | 500 |
| **Degradation** | Misspecified distribution (Normal → t) | df = {3, 5, 10} | 500 |
| **Degradation** | Nonstationary baseline | drift = {0.005, 0.01, 0.02}·t/N | 500 |
| **Filter** | Period threshold sweep | {2, 4, 6, ..., 50} | 100 |
| **Mixed** | Two forcings superimposed | ratio = {0.2:0.8, 0.5:0.5, 0.8:0.2} | 100 |

Degradations applied to all five classes so macro-F1 is computable. Mixed dynamics are descriptive only (no correct label, report label distribution and entropy).

## Predictions

- Amplitude: detection thresholds differ across bins. At equalized difficulty, macro-F1 > 0.8 but < 1.0.
- Ablation: curvature removal breaks convergence/divergence. Period filter removal breaks aperiodic/oscillatory.
- Autocorrelated null: FPR increases with φ. At φ = 0.8, null FPR > 10%.
- Correlated streams: composition advantage shrinks. At ρ = 0.6, composed F1 approaches standardized sum.
- Missing data: graceful up to 25%.
- Misspecified distribution: minimal at df ≥ 5. At df = 3, feature estimates inflate.
- Nonstationary baseline: mimics divergence. At drift = 0.02, null FPR > 20%.
- Period filter: below ~8, aperiodic → oscillatory. Above ~20, short-period oscillatory missed.
- Mixed: hierarchy determines label. At 0.5:0.5, earlier test wins.

## E-value validity

Separate from classifier F1. For the t-distribution misspecification under null:
- Anytime exceedance: does max_t E_t exceed 20 in ≤5% of reps?
- Supermartingale: max over time of mean(e_t) across 1000 reps. Inflated if > 1.05.

## Falsification

- Lower CI < 0.8 at mildest severity of any degradation → fragile.
- All bins > 90% detection at 50% of V3 amplitude → V3 was too easy.
- Removing period filter has no effect → filter unnecessary.
- ρ = 0.3 reduces composed F1 to standardized sum → Fisher advantage is independence artifact.
- Type-I rate > 10% under t(df=3) → Normal-assumption e-value invalid for heavy tails.

## Outputs

- `results/tables/v4_grid.csv` — all perturbations, all metrics, one table
- `results/plots/v4_sensitivity.png` — amplitude curves per bin
- `results/plots/v4_degraded.png` — F1 heatmap: degradation × severity
- `results/plots/v4_ablation.png` — paired delta bar chart
- `results/plots/v4_mixed.png` — label distribution at each dominance ratio
