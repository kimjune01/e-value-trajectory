# Pre-registration: Robustness analysis

V3 achieved F1 = 1.000 on clean synthetic data. V4 asks: how fragile is it?

## Primary endpoint

- **Metric:** macro-F1 over five classes (null, convergent, divergent, oscillatory, aperiodic).
- **Pass:** lower 95% bootstrap CI > 0.8 at mildest severity of each degradation.
- **Inconclusive:** lower CI ≤ 0.8 but point estimate ≥ 0.6.
- **Fail:** point estimate < 0.6 at mildest severity.
- **Harsher severities are descriptive, not pass/fail.**
- Bootstrap: 1000 samples. Unit of resampling: trajectory.

## Seed policy

Calibration: 80000–80999 (1000 reps). Evaluation: 99999+, disjoint. Ablations use identical trajectories (paired seeds). No threshold retuning after V4 starts.

## Design

One pipeline, many perturbations. Each perturbation modifies one aspect of V3's setup, runs the same generate → compose → classify → score loop. Degradations and ablations report macro-F1 + confusion matrix + paired comparison against standardized sum. Mixed dynamics report label distribution and entropy (no correct label).

### Perturbation grid

| Category | Perturbation | Severity levels | Reps |
|---|---|---|---|
| **Amplitude** | Sweep each bin 0 to 2× V3 | 21 steps (inclusive endpoints) | 100 |
| **Amplitude** | Equalized difficulty (80% detection per bin + null) | 1 | 100 |
| **Ablation** | Remove curvature test | 1 | 100 |
| **Ablation** | Remove period > 10 filter | 1 | 100 |
| **Ablation** | Remove 0-1 test (PE only) | 1 | 100 |
| **Ablation** | Remove PE (0-1 only) | 1 | 100 |
| **Ablation** | Remove monotone test (skip to curvature) | 1 | 100 |
| **Ablation** | Remove Mann-Kendall (slope only) | 1 | 100 |
| **Ablation** | Curvature threshold ±20% | 2 | 100 |
| **Degradation** | Autocorrelated null (AR(1)) | φ = {0.2, 0.5, 0.8} | 500 |
| **Degradation** | Correlated streams | ρ = {0.1, 0.3, 0.6} | 500 |
| **Degradation** | Missing data | {5%, 10%, 25%} | 500 |
| **Degradation** | Misspecified distribution (Normal → t) | df = {3, 5, 10} | 500 |
| **Degradation** | Nonstationary baseline | drift_coef = {0.005, 0.01, 0.02} | 500 |
| **Filter** | Period threshold sweep | {2, 4, 6, ..., 50} | 100 |
| **Mixed** | trend+oscillation | ratio = {0.2:0.8, 0.4:0.6, 0.5:0.5, 0.6:0.4, 0.8:0.2} | 100 |
| **Mixed** | oscillation+aperiodic | ratio = {0.2:0.8, 0.4:0.6, 0.5:0.5, 0.6:0.4, 0.8:0.2} | 100 |
| **Mixed** | convergence+oscillation | ratio = {0.2:0.8, 0.4:0.6, 0.5:0.5, 0.6:0.4, 0.8:0.2} | 100 |
| **Mixed** | weak divergence+heavy-tailed noise (df=3) | ratio = {0.2:0.8, 0.4:0.6, 0.5:0.5, 0.6:0.4, 0.8:0.2} | 100 |

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

## Operational appendix

Pins every implementation decision. Two independent implementers reading this appendix should produce identical results given the same seeds.

### V3 normative reference

Pinned commit: `63fda12` (`~/e-value-trajectory`). Authoritative files: `src/fourbin.py`, `configs/conditions.yaml`. N=10,000 per stream, K=5 streams (temperature/Normal, fish_count/Poisson, inter_event/Exponential, turbidity/Bernoulli, dissolved_oxygen/Lognormal). Classifier thresholds from 1000-rep null calibration (seeds 80000–80999). Base class amplitudes from `src/fourbin.py`: null=0.0, divergence=1.5, convergence=1.5, oscillation=0.1, aperiodic=0.5. These are distinct from per-stream amplitudes in `conditions.yaml`.

### Classifier

Do NOT delegate to V3 `classify()`. Reimplement the decision tree with two tunable knobs:

- `period_min` (default 10): minimum period for oscillatory classification.
- `curve_mult` (default 1.0): multiplier on `R_curve_q05` threshold.

Decision tree (identical logic to V3 but with knobs):
```
if |S| > S_q99 and |Z_MK| > 2.58 and |med_last - med_first| > 3*MAD_dx:
    if R_curve < R_curve_q05 * curve_mult:
        → convergent
    else:
        → divergent
elif G_spec > G_spec_q99 and Q_spec > 5 and period_min < peak_period < N/8:
    → oscillatory
elif K_01 > 0.8 and 0.55 <= PE <= 0.95:
    → aperiodic
else:
    → null
```

### Ablation operations

Each ablation neutralizes one feature before classification. The tree runs unchanged; only the input differs.

| Ablation | Operation | Effect |
|---|---|---|
| Remove monotone | Set S=inf, Z_MK=inf | Monotone check always passes; classifier skips to curvature |
| Remove Mann-Kendall | Set Z_MK=inf | MK always passes; monotone depends only on S and median diff |
| Remove curvature | Set R_curve=0.0 | Curvature always below threshold; monotone always → convergent |
| Remove period filter | Set period_min=0 | All spectral peaks accepted regardless of frequency |
| Remove 0-1 test | Set K_01=1.0 | 0-1 test always passes; aperiodic depends only on PE |
| Remove PE | Set PE=0.75 | PE always in range; aperiodic depends only on K_01 |
| Curvature ±20% | Set curve_mult=0.8 or 1.2 | Stricter or looser curvature separation |

Logic: to "remove" a test in a boolean AND chain, set the neutralized variable to a value that always satisfies the check (not one that always fails it).

### Seed schedule

- Null calibration: seeds 80000–80999 (1000 reps).
- Evaluation: seed = 99999 + rep + CONDITIONS.index(condition) * 10000. CONDITIONS order is from `src/fourbin.py`: ["null", "divergence", "convergence", "oscillation", "aperiodic"]. Reps are per condition per severity. Disjoint from calibration.
- Random draw order within each rep: (1) AR(1) eps if autocorrelated, (2) shared latent z if correlated, (3) stream generation in stream_names order, (4) MCAR masks per stream, (5) t-noise replacement if misspecified. This order is normative for seed reproducibility.
- Ablations: same seeds as baseline evaluation (paired comparison).
- Bootstrap: per-sample seed = `int(hashlib.sha256(f"{category},{perturbation},{severity},{signal},{b}".encode()).hexdigest()[:8], 16)` where `b` is the bootstrap iteration index (0–999). This ensures each of the 1000 samples draws a different resample. Stable across Python processes.

### Amplitude sweep

- 21 steps from 0 to 2× V3 amplitude, inclusive (np.linspace(0, 2*base, 21)).
- One bin swept at a time; other bins stay at V3 amplitude.
- Report per-bin detection rate (fraction of reps where the swept bin's ground truth matches prediction) AND five-class macro-F1.

### Equalized difficulty

- After the sweep, for each bin find the amplitude closest to 80% per-bin detection rate (from the sweep data).
- Run the classifier at those four amplitudes + null (amplitude 0). Report five-class macro-F1.

### Degradation operators

Each degradation modifies the data generation for ALL five conditions, not just null. The forcing and stream generation proceed as in V3, then the degradation is applied.

**Autocorrelated (AR(1)):**
- Generate AR(1) noise: `ar[0] = 0; ar[t] = phi * ar[t-1] + eps[t] * sqrt(1 - phi^2)` where `eps ~ N(0,1)`.
- Marginal variance of `ar` is 1 (by construction from the `sqrt(1 - phi^2)` scaling).
- Add `0.1 * ar` to the forcing signal before stream generation. This adds autocorrelated structure to all conditions, not just null.

**Correlated streams:**
- Generate one shared latent `z ~ N(0, 1)` of length N per rep.
- Before generating each stream, add correlated noise to forcing: `forcing_corr = forcing + rho * 0.3 * z`.
- Do NOT scale down the deterministic forcing. The signal amplitude stays unchanged; only correlated noise is added. This isolates the effect of correlation from signal attenuation.
- Generate streams from `forcing_corr`. Marginal distributions are preserved because corruption enters through the forcing, not the observations.

**Missing data (MCAR):**
- For each stream independently, draw a mask: `mask = rng.random(N) < frac`.
- Set `log_e[mask] = 0` (neutral evidence: the missing observation contributes nothing to the composed trajectory). Do NOT interpolate or drop.
- The observation `x` is also masked for the standardized sum: `x[mask] = NaN`. Per-stream z-scores use `nanmean` and `nanstd`. Missing z-scores (where x was NaN) are filled with 0 before summing across streams (neutral contribution). Sum uses `np.sum`, not `np.nansum`.

**Misspecified distribution (Normal → t):**
- Replace only the temperature/Normal stream's noise with t-distributed noise.
- Scale: `x = forcing + rng.standard_t(df, N) * sigma * sqrt((df-2) / df)` for df > 2, preserving marginal variance equal to V3's Normal(0, sigma^2). (Standard t has variance df/(df-2); multiply by sqrt((df-2)/df) to normalize to variance 1, then by sigma.)
- E-value formula unchanged (still uses Normal lambda). This is the misspecification.
- All other streams (Poisson, Exponential, Bernoulli, Lognormal) unchanged.
- Classifier still uses V3 thresholds (no recalibration).

**Nonstationary baseline:**
- Let `drift_coef` ∈ {0.005, 0.01, 0.02}. Add `drift_coef * t / N` to the forcing for all conditions, where `t = np.arange(N)`.
- At t=N, the total drift added equals `drift_coef`. Additive, positive. Applied before stream generation.

### Mixed dynamics

- Two forcings are summed: `forcing = w1 * A1 * f1(t) + w2 * A2 * f2(t)` where A1, A2 are V3 amplitudes, f1/f2 are the forcing functions, and w1+w2=1.
- Ratios: (w1, w2) in {(0.2, 0.8), (0.4, 0.6), (0.5, 0.5), (0.6, 0.4), (0.8, 0.2)}.
- "weak divergence + heavy-tailed noise": `forcing = w1 * A_div * f_div(t) + (1-w1) * 0.3 * rng.standard_t(3, N)`. The t-noise is the second component, not a null amplitude.
- No correct label. Report label distribution and Shannon entropy (log base 2).

### E-value validity

- Run 1000 reps under null with the temperature stream using t(df=3) noise (same as misspecified degradation). Other streams use their V3 null distributions.
- For each rep: compute `cum_log_e = cumsum(composed_log_e)` and `max_E = exp(max(cum_log_e))`.
- Anytime exceedance rate: fraction of reps where `max_E > 20`.
- Supermartingale diagnostic: accumulate `mean_Et[t] = mean over reps of exp(cum_log_e[t])` at each time point. This checks the **cumulative** e-value process E_t, not the single-step e_t. Report `max over t of mean_Et`. Inflated if > 1.05.
- Separate seeds from degradation reps.

### Macro-F1

- Standard macro-averaged F1: compute precision and recall per class, F1 per class, average across all 5 classes.
- Zero-division: if TP+FP=0 or TP+FN=0, F1 for that class is 0.
- Computed over pooled trajectories (all reps within a cell).

### Bootstrap CI

- Stratified percentile bootstrap: resample trajectories within each true class (preserving class balance), recompute macro-F1 on each bootstrap sample.
- 1000 bootstrap samples. CI = [2.5th percentile, 97.5th percentile].
- Per-sample seed: `int(hashlib.sha256(f"{category},{perturbation},{severity},{signal},{b}".encode()).hexdigest()[:8], 16)` where `b` = bootstrap iteration index. Each sample gets a unique seed.

### Standardized sum baseline

- Per stream: `z = (x - mean(x)) / std(x)`. For missing data: use `nanmean` and `nanstd`.
- Sum across streams: `z_sum = sum of z_k`.
- **Separate null calibration for z_sum:** run 1000 null reps, compute features on z_sum, extract thresholds (S_q99, R_curve_q05, G_spec_q99) from the z_sum null distribution. The standardized sum operates on a different scale than composed log_e; applying composed thresholds would produce ~100% misclassification.
- Run the V4 classifier with z_sum-calibrated thresholds on `z_sum` features.
- Report macro-F1 and paired delta (composed F1 minus standardized F1) per cell.

### Output schema

One CSV: `results/tables/v4_grid.csv`. Required columns:

```
category, perturbation, severity, signal, macro_f1, ci_lo, ci_hi,
null_fpr, label_entropy, label_counts, detection_rate, paired_delta
```

- `signal`: "composed" or "standardized_sum".
- `macro_f1`: empty for mixed dynamics.
- `label_counts`: JSON dict of label → count.
- `detection_rate`: per-bin detection for amplitude sweep; empty otherwise.
- `paired_delta`: composed F1 minus standardized F1; empty for standardized rows and mixed.

## Falsification

- Lower CI < 0.8 at mildest severity of any degradation → fragile.
- All bins > 90% detection at 50% of V3 amplitude → V3 was too easy.
- Removing period filter has no effect → filter unnecessary.
- ρ = 0.3 reduces composed F1 to standardized sum → Fisher advantage is independence artifact.
- Type-I rate > 10% under t(df=3) → Normal-assumption e-value invalid for heavy tails.

## Outputs

- `results/tables/v4_grid.csv` — all perturbations, all metrics, one table
- `results/plots/v4_sensitivity.png` — 4 subplots (one per bin), x=amplitude, y=detection rate, with 95% CI shaded. Horizontal line at 80% and 90%.
- `results/plots/v4_degraded.png` — heatmap, rows=degradation type, cols=severity, cell color=composed macro-F1, text=F1 value. Colorscale: red(<0.6) → yellow(0.8) → green(>0.9).
- `results/plots/v4_ablation.png` — horizontal bar chart, one bar per ablation, x=paired delta in macro-F1 vs baseline, sorted by magnitude. Error bars: bootstrap CI on the delta, computed by resampling paired trajectories by shared (condition, rep) and computing delta-F1 on each bootstrap sample.
- `results/plots/v4_mixed.png` — stacked bar chart per mixed condition, faceted by ratio, bars=label counts, colored by label.
