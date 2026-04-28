# Pre-registration: Heterogeneous e-value composition

Secondary pre-registration. Claims 1–3 tested in V1 (see `PREREGISTRATION_V1.md` and `report.md`). This document covers Claim 4 only.

## What V1 established

For a single univariate stream, the e-value periodogram is an affine transform of the raw-data periodogram: Δlog(e_t) = λX_t − λ²/2. No spectral advantage. The e-value's value proposition is composability, not amplification.

## Claim 4: Heterogeneous e-value composition recovers shared periodic structure

Heterogeneous e-value composition provides a common scale for aggregating weak periodic evidence across experiments with different observation models. In a simulated five-stream setting, individual e-value periodograms remain below the calibrated detection threshold, while the periodogram of the summed log e-values reliably recovers the shared period. This does not give a spectral advantage within any single stream; the gain comes from valid cross-experiment composition.

### Hypotheses

- **H0:** All five streams follow their null distributions, independent across streams and time, with no tidal component.
- **H1:** Each stream has a weak shared sinusoidal component with known period T=500 and known phase (φ=0).

### Scenario: Water quality monitoring

Five sensors at a coastal station, each measuring a different quantity, all sharing a tidal cycle (period T=500, matching V1's timescale). Each stream is individually too noisy to detect the period; composed e-values reveal it.

**Independence assumption:** This demo assumes independent streams. Real co-located sensors would likely be correlated under both null and alternative. The validity argument depends on independence; the claim is narrow to this setting.

### Streams with time-indexed alternatives

Each e-value uses the time-indexed likelihood ratio e_t = p(X_t | θ_alt,t) / p(X_t | θ_0), where θ_alt,t is fixed before observing X_t. Under H0, E[e_t | past] = 1. The alternative encodes a pre-specified hypothesis about the tidal cycle, not knowledge of the data.

| # | Stream | Distribution | H0 | Link | Alternative |
|---|---|---|---|---|---|
| 1 | Temperature anomaly | Normal | N(0, 1) | identity | μ_alt,t = A₁·sin(2πt/T) |
| 2 | Fish count per trawl | Poisson | Pois(μ₀) | log | μ_alt,t = μ₀·exp(A₂·sin(2πt/T)) |
| 3 | Inter-event time | Exponential | Exp(r₀) | log | r_alt,t = r₀·exp(A₃·sin(2πt/T)) |
| 4 | Turbidity exceedance | Bernoulli | Bern(p₀) | logit | logit(p_alt,t) = logit(p₀) + A₄·sin(2πt/T) |
| 5 | Dissolved oxygen | Lognormal | LN(μ₀, σ) | identity on log | δ_t = A₅·sin(2πt/T) added to log-mean |

Log/logit links ensure Poisson rates stay positive, Bernoulli probabilities stay in (0,1), Exponential rates stay positive.

### E-value formulas

| # | log(e_t) |
|---|---|
| 1 | μ_alt,t · X_t − μ_alt,t² / 2 |
| 2 | X_t · log(μ_alt,t / μ₀) − (μ_alt,t − μ₀) |
| 3 | log(r_alt,t / r₀) − (r_alt,t − r₀) · X_t |
| 4 | X_t · log(p_alt,t / p₀) + (1−X_t) · log((1−p_alt,t) / (1−p₀)) |
| 5 | δ_t · (log X_t − μ₀) / σ² − δ_t² / (2σ²) |

### Composition

log(E_composed,t) = Σ_k log(e_t^(k)). Valid because streams are independent and each e_t^(k) satisfies E[e_t] ≤ 1 under H0. The product of valid e-values is a valid e-value.

### Amplitude calibration

Set per-stream amplitudes so each individual peak/median ratio ≈ 8–10 (below the calibrated detection threshold of 15 from V1). With K=5, the composed peak/median should be well above threshold.

### Baseline comparison

**Standardized raw-data sum:** z-score each stream (subtract mean, divide by std), sum, take periodogram. Expected result: comparable spectral detection power to e-value composition in this toy setting. The standardized sum also works for spectral detection but does not give a valid sequential test — it is not an e-value.

Both methods reported in the results table.

### Null controls

Run all five streams under global null (no tidal component, A_k = 0 for all k), 100 replications. Report:
- Per-stream peak/median distribution under null
- Composed peak/median distribution under null (95th and 99th percentile)
- Empirical null threshold for composed signal (do not borrow the single-stream threshold of 15)

### Predictions

1. Each individual stream's peak/median at T=500 is below the single-stream detection threshold (15) in ≥ 80% of replications.
2. The composed e-value's peak/median at T=500 exceeds the empirical null threshold (99th percentile of composed null distribution) in ≥ 90% of replications.
3. The true period T=500 is the dominant peak (rank 1) in the composed periodogram in ≥ 80% of replications.
4. The standardized raw-data sum shows comparable spectral detection to the composed e-value (TOST equivalence of peak/median ratios, margin d=0.2).
5. Per-stream e-value means under null are within 5% of 1.0.

### Falsification

- If the composed peak/median does not exceed the null threshold in ≥ 90% of reps, composition does not reliably recover the shared period.
- If standardized-sum detection is strictly superior (not merely equivalent), the e-value adds no value even for composition.
- If per-stream e-value null means deviate substantially from 1.0, the e-value construction is invalid.

### Measurement

Replication summaries (100 reps, both signal and null):
- Median and IQR of peak/median at period 500, per stream and composed
- P(individual stream exceeds threshold)
- P(composed exceeds threshold)
- P(standardized-sum exceeds threshold)
- Rank of true period among periodogram peaks
- Empirical null mean of e-values per stream

### Experiment parameters

- N = 10,000 observations per stream per replication
- K = 5 streams
- 100 replications (signal) + 100 replications (null)
- Period T = 500
- Phase φ = 0 (shared, known)
- Base seed: 42 (same as V1)
- Specific per-stream parameters (μ₀, r₀, p₀, σ, amplitudes A_k) to be set in `configs/conditions.yaml` before running

### Outputs

- `results/tables/composition.csv` — per-replication summary statistics
- `results/plots/composition_periodograms.png` — 2×3 grid: 5 individual + 1 composed (rep 0, illustrative)
- `results/plots/composition_trajectories.png` — composed cumulative log(E_t) showing oscillation structure

### What this does not test

- Correlated streams (real-world sensors sharing environmental forcing would be correlated)
- Unknown or misspecified period (the alternative encodes the true T=500)
- Adaptive composition (λ or alternative updated based on observations)
- Multi-round experimental design (iterating perturb → classify → select next experiment)

These are separate experiments.
