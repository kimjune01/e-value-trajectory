# E-value trajectory diagnostics: experiment report

All results, all figures, all failures. Nothing curated.

## Setup

Seven conditions (a–e, i, b2) for Claims 1–2, three regime-switch conditions (f–h) for Claim 3. N = 10,000 observations per condition, 100 replications each. E-value: e_t = exp(λX_t − λ²/2), λ = 0.3. All parameters in `configs/conditions.yaml`.

## The affine transform result

The per-step e-value increment is Δlog(E_t) = λX_t − λ²/2. This is an affine transform of the raw data. Consequences:

1. The periodogram of Δlog(E_t) is λ² × the periodogram of X_t, shifted by a constant that vanishes after centering.
2. The peak/median ratio is invariant to λ (confirmed: identical across λ = 0.15, 0.3, 0.6 for all conditions).
3. TOST equivalence test: mean difference in peak/median ratio between e-value and raw = 0.0000 for all seven conditions (p ≈ 0).

**For a single univariate stream, the e-value periodogram carries exactly the same spectral information as the raw-data periodogram.** This is not a failure — it's a theorem. The e-value's value proposition is composability across heterogeneous experiments (different DGPs, scales, test statistics), not spectral amplification within a single stream.

## Claim 1: Cyclic systems produce oscillating e-value trajectories

**Verdict: Supported, with qualifications.**

### Periodogram classification works

All cyclic conditions produce spectral peaks significantly above the null noise floor (Mann-Whitney U = 10,000, p = 1.28 × 10⁻³⁴, surviving Bonferroni correction at α/90 = 0.00056):

| Condition | Median peak/median | IQR | Detected period |
|---|---|---|---|
| (c) Null | 13.1 | 11.8–14.1 | random |
| (a) Stationary | 12.7 | 11.9–13.8 | random |
| (i) Sinusoidal | 318.9 | 300.5–343.2 | 500 ✓ |
| (b) LV | 2217.9 | 2173.6–2277.0 | 12.3 |
| (b2) Stochastic LV | 249.6 | 206.2–288.7 | 12.3 |
| (d) AR(1) | 920.1 | 787.8–1056.9 | ~196 |
| (e) Regime switching | varies | — | none clean |

Full-length periodogram (10,000 points). E-value and raw periodograms are identical (affine transform).

### Qualifications

1. **LV period is ~12, not ~500.** The chosen parameters (α=1.0, β=0.1, δ=0.075, γ=1.5) produce fast oscillations. This is an emergent property, not a bug. The diagnostic detects the actual period.

2. **AR(1) confound triggers.** AR(1) with φ = 0.9 produces a strong spectral peak (median peak/median = 920, broad but above any reasonable threshold). The spectral classifier cannot distinguish autocorrelation from true periodicity using peak height alone. Peak *width* could differentiate them (AR(1) is broad, LV is sharp), but the pre-registered classifier doesn't check width.

3. **Sensitivity to λ: none.** Periodogram peak/median ratios are numerically identical at λ = 0.15, 0.3, and 0.6. Direct consequence of the affine transform.

![Periodograms](results/plots/periodograms.png)

### Pre-registration prediction scorecard

| Prediction | Result |
|---|---|
| (a) No spectral peak | ✓ peak/median = 12.7 (noise floor) |
| (i) Sharp peak at period 500 | ✓ period 500 detected perfectly |
| (b) Sharp peak at LV period | ✓ period 12.3 (not 500 as predicted, but detected) |
| (b2) Broader, noisier peak | ✓ peak/median 250 vs 2218 for deterministic |
| (c) Flat periodogram | ✓ peak/median = 13.1 (noise floor) |
| (d) No peak (autocorrelation ≠ periodicity) | ✗ peak/median = 920 (strong peak) |
| (e) No clean peak | ✓ mixed, no single dominant period |

Prediction (d) failed: AR(1) does produce a spectral peak. Lemma 4 says the spectral density is "continuous with no discrete spike," which is true — but a broad continuous peak still produces a high peak/median ratio in a finite-sample periodogram. The lemma is correct; the prediction confused "no discrete spike" with "no detectable peak."

## Claim 2: Trajectory shape classifies faster than threshold crossing

**Verdict: The comparison is categorical, not competitive.**

### Stopping times (calibrated peak_threshold = 15)

| Condition | Threshold | Bayesian | Spectral EV |
|---|---|---|---|
| (c) Null | 500 | 52 | 499 |
| (a) Stationary | 53 | 12 | 499 |
| (d) AR(1) | 22 | 5 | 499 |
| (i) Sinusoidal | 35 | 10 | 499 |
| (b) LV | 48 | 7 | 499 |
| (b2) Stochastic LV | 42 | 8 | 499 |
| (e) Regime switch | 66 | 12 | 499 |

Median stopping time across 100 reps. Spectral EV always commits at t = 499 (first full window).

### Label accuracy

| Condition | Ground truth | Threshold | Bayesian | Spectral EV |
|---|---|---|---|---|
| (c) Null | null | 92% null ✓ | 41% reject, 27% periodic, 32% null | 100% null ✓ |
| (a) Stationary | reject_null | 100% reject ✓ | 96% reject ✓ | 100% reject ✓ |
| (d) AR(1) | reject_null | 90% reject ✓ | 62% reject ✓ | 100% periodic ✗ |
| (i) Sinusoidal | periodic | 100% reject ✗ | 93% reject ✗ | 53% periodic ✓ |
| (b) LV | periodic | 100% reject ✗ | 94% reject ✗ | 100% periodic ✓ |
| (b2) Stoch LV | periodic | 100% reject ✗ | 91% reject ✗ | 88% periodic ✓ |
| (e) Regime | ambiguous | 81% reject | 91% reject | 46% periodic |

![Classifier labels](results/plots/classifier_labels.png)

### Interpretation

The three classifiers answer different questions:

- **Threshold** and **Bayesian** answer "is there an effect?" They commit fast (10–50 steps) but have no mechanism to detect periodicity. They are correct about effect presence but blind to dynamics.
- **Spectral EV** answers "is the evidence oscillating?" It's the only classifier that labels periodic correctly on LV conditions. But it requires a minimum of 500 observations, and it can't distinguish AR(1) from true periodicity.

The pre-registered prediction — "spectral e-value classifier labels 'periodic' before threshold classifier labels 'reject null'" — is **false**. The spectral classifier always fires at t = 499 (the earliest possible), while the threshold classifier fires at t ≈ 50. The spectral classifier is categorically slower because it needs a full window of data to compute a periodogram.

The useful comparison isn't speed — it's what each classifier can see. Only the spectral classifier detects periodicity. Only the threshold/Bayesian classifiers give fast effect-present/absent decisions.

### Pre-registration calibration failure

The pre-registered spectral peak threshold of 3.0 produced 100% false positive rate on all conditions (including null). A 500-point periodogram with ~250 bins has an expected max/median ratio of ~9 under white noise. The threshold was set without accounting for multiple comparisons across frequency bins. Calibrated to 15.0 using the null condition's empirical distribution (p99 = 12.6).

Similarly, the Bayesian classifier's periodogram on the posterior mean produces 27% false "periodic" labels on the null condition. The posterior mean is autocorrelated by construction (each value shares most observations with its predecessor), inflating spectral peaks. This is an inherent limitation of periodogram-on-posterior approaches.

## Claim 3: Regime changes appear as regime changes in the evidence

**Verdict: Supported for CUSUM on e-value increments. Bayesian CPD fails.**

### CUSUM detection (calibrated h = 20)

| Condition | Signal | Median delay | IQR | Within 500 | FA rate |
|---|---|---|---|---|---|
| (f) Effect → null | e-value | 260 | 208–306 | 99/100 | — |
| (f) Effect → null | raw | 4751 | 4623–4829 | 0/100 | — |
| (g) Null → effect | e-value | 256 | 200–319 | 99/100 | — |
| (g) Null → effect | raw | 4757 | 4629–4842 | 0/100 | — |
| (h) Effect → reversed | e-value | 118 | 100–132 | 100/100 | — |
| (h) Effect → reversed | raw | 4798 | 4655–4860 | 0/100 | — |
| (a) Stationary | e-value | — | — | — | 0% |
| (a) Stationary | raw | — | — | — | 100% |
| (c) Null | e-value | — | — | — | 2% |
| (c) Null | raw | — | — | — | 100% |

Detection delay = |t_detected − 5000|. All e-value delays within the pre-registered 500-step budget (p < 10⁻⁴⁸, surviving Bonferroni).

![Changepoint detection](results/plots/changepoint_delays.png)

### Why CUSUM on raw data fails here

CUSUM on raw X_t with the same parameters (k = 0.0225, h = 20) has 100% false alarm rate because these parameters are calibrated for e-value variance (σ² = 0.09), not raw variance (σ² = 1.0). The e-value transformation compresses noise by factor λ = 0.3. With equivalently scaled parameters (h ≈ 220 for raw data), detection performance would be the same — the affine transform guarantees it. The head-to-head with identical parameters is not a fair comparison; it tests parameter sensitivity, not signal quality.

### Bayesian CPD

Poor across all conditions. False alarm rate 22–54%. Detection delay IQR spans thousands of steps. Hazard rate 1/2000 expects changepoints too frequently for a 10,000-step sequence. The pre-registered parameters don't work. Unlike CUSUM, the Bayesian CPD wasn't recalibrated because the poor performance is more fundamental — the model assumes stationary segments with a known observation model, but the e-value increment distribution changes in ways the Normal-Normal conjugate model doesn't capture well.

### Pre-registration calibration failure

CUSUM h = 5 (pre-registered as "standard") produced 100% false alarm rate on all conditions. Root cause: h = 5 gives ARL₀ ≈ 1,083 for signals with σ = 0.3. Over 10,000 observations, at least one false alarm is near-certain. Calibrated to h = 20 (null FA rate 2%).

## Figures

| Figure | Location |
|---|---|
| E-value trajectories | `results/plots/trajectories.png` |
| Periodograms | `results/plots/periodograms.png` |
| Classifier stopping times | `results/plots/classifier_stopping_times.png` |
| Classifier label distributions | `results/plots/classifier_labels.png` |
| Changepoint detection delay | `results/plots/changepoint_delays.png` |
| Regime-switch trajectories | `results/plots/changepoint_trajectories.png` |
| Sensitivity (λ variants) | `results/plots/sensitivity.png` |

## Tables

| Table | Location |
|---|---|
| Spectral summary (full-length) | `results/tables/spectral_summary.csv` |
| Classifier race | `results/tables/classifier_race.csv` |
| Changepoint detection | `results/tables/changepoint_detection.csv` |
| Sensitivity analysis | `results/tables/sensitivity.csv` |

## Failures

1. **Prediction (d) wrong.** AR(1) produces a detectable spectral peak. The prediction confused "no discrete spike in the spectral density" (Lemma 4, correct) with "no detectable peak in a finite-sample periodogram" (wrong). A broad continuous peak with enough power exceeds any fixed peak/median threshold.

2. **Spectral threshold = 3 miscalibrated.** 100% false positive rate. Required post-hoc calibration to 15.

3. **CUSUM h = 5 miscalibrated.** 100% false alarm rate. Required post-hoc calibration to 20.

4. **Bayesian CPD fails.** 22–54% false alarm rate, huge variance in detection delay. Pre-registered parameters (hazard = 1/2000) are inappropriate.

5. **Bayesian classifier posterior mean periodogram fails.** 27% false "periodic" on null. The posterior mean is too autocorrelated for naive periodogram analysis.

6. **Claim 2 prediction wrong.** "Spectral e-value labels periodic before threshold labels reject_null" is false. Spectral always fires at t = 499; threshold fires at t ≈ 50.

7. **Sinusoidal borderline.** 53% periodic detection in 500-point windows (1 cycle per window is barely enough for spectral resolution).

## Verdict

### What works

E-value trajectories are transparent to system dynamics. The periodogram of Δlog(E_t) detects the cycle period of Lotka-Volterra systems (both deterministic and stochastic), and CUSUM on e-value increments detects regime switches within 260 steps (median), well within the 500-step budget.

### What doesn't add value

For a single univariate normal stream, the e-value adds nothing the raw data doesn't already carry. Δlog(e_t) = λX_t − λ²/2 is an affine transform. Same periodogram, same spectral peaks, same classification power. This is mathematically necessary and was confirmed empirically (TOST equivalence p ≈ 0, sensitivity invariant to λ).

### What matters

The e-value's value isn't spectral amplification — it's **composability**. E-values from heterogeneous experiments (different DGPs, scales, test statistics) can be multiplied. Raw data from different experiments can't be concatenated. The trajectory of a composed e-value should still carry dynamics the individual streams don't. This is Candidate Claim 4 in EXPERIMENT.md — a new experiment, not tested here.

### The honest summary

The experiment confirms that e-value trajectory shape reflects system dynamics (Claims 1 and 3 supported). It also shows that for a single stream, this reflection is an affine transform of the raw data — no amplification, no advantage. The spectral classifier detects periodicity that threshold/Bayesian classifiers can't see, but it's slower and can't distinguish autocorrelation from true cycles. Two of three pre-registered detection parameters were miscalibrated (spectral threshold, CUSUM h), requiring post-hoc correction. Bayesian CPD failed entirely. The experiment answered the pre-registered questions, but the interesting question — whether composing e-values across heterogeneous experiments preserves dynamics that no single stream reveals — remains untested.
