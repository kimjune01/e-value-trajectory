# Pre-registration: Forcing-pattern classification on composed e-value trajectories

Third pre-registration. V1 tested single-stream diagnostics (Claims 1–3). V2 tested heterogeneous composition for periodicity detection (Claim 4). This document covers Claim 5: classifying composed trajectories into forcing patterns using a hierarchical decision rule.

## Framing

This experiment classifies **shared forcing patterns**, not intrinsic dynamical regimes. The composed e-value trajectory is a noisy scalar time series produced by aggregating evidence across heterogeneous experiments. We classify its morphology — trending, decaying, oscillating, aperiodically forced, or null — using a fixed decision tree with null-calibrated thresholds.

The detection methods for each pattern are individually established (see citations below). No single paper assembles them into a unified four-bin pipeline. The pipeline itself is the contribution.

## What prior work establishes (no experiment needed)

### Single-stream classification is invariant

For a single univariate stream with fixed alternative, log(e_t) = cX_t + d (affine transform). All classification methods are invariant under invertible affine transforms: Lyapunov exponents unchanged (Kantz & Schreiber, 2004, §3.1–3.3), 0-1 test invariant under scaling/translation (Gottwald & Melbourne, 2009), periodogram peak ratios cancel the scale factor, trend tests sign-preserved for c > 0. No experiment needed.

### E-value composition is valid

The product of independent e-values is a valid e-value (Vovk & Wang, 2021; Grünwald et al., 2024). V2 demonstrated 3.4× power advantage of composed e-values over standardized sums for periodicity detection, via Fisher information weighting.

## What this experiment tests

Composed e-value trajectories are not an affine transform of any single raw stream. V2 showed composition works for oscillation. This experiment tests whether it works for all four forcing patterns, using a single hierarchical classifier that assigns exactly one label per trajectory.

## Claim 5: Hierarchical forcing-pattern classification

**Question:** When K heterogeneous experiments share a common forcing pattern, does a hierarchical classifier on the composed e-value trajectory correctly identify the pattern — and does composition improve classification over individual streams?

### Hypotheses

- **H0:** All K streams follow their null distributions (no shared forcing).
- **H1-diverge:** Shared linear drift in all streams' means.
- **H1-converge:** Shared exponential decay (perturbation returns to baseline).
- **H1-oscillate:** Shared periodic forcing (reuse V2 setup, period T=500).
- **H1-chaos:** Shared chaotic forcing from a logistic map.

### Streams

K=5 streams, same distributions and fixed alternatives as V2 (Normal, Poisson, Exponential, Bernoulli, Lognormal). N=10,000 observations.

### Forcing generators

| Pattern | Forcing on stream k's mean | Parameters |
|---|---|---|
| Null | None (constant at null) | — |
| Divergence | A_k · (t/N − 0.5) | Centered linear ramp |
| Convergence | A_k · exp(−t/τ) | τ = 2000, decays toward null |
| Oscillation | A_k · sin(2πt/T) | T = 500, reuse V2 |
| Chaos | A_k · z̃_t (standardized logistic map) | r = 3.9, z̃ = (z − mean)/std after 1000-step burn-in |

Chaos forcing: discard first 1000 iterates. Standardize z̃_t to zero mean and unit variance. Vary z_0 per replication (z_0 = 0.1 + 0.008·rep) to produce distinct chaotic trajectories, not just noise replications on a single trajectory.

Data generation uses log/logit links (same as V2) to keep parameters valid.

### Amplitude calibration (locked)

Calibration uses a separate seed range (seeds 90000–90099, 100 reps) from evaluation (seeds 99999+, 100 reps).

Procedure:
1. For each forcing pattern and each stream, sweep amplitudes.
2. Select the smallest A_k where the median individual-stream detection statistic is between the null 75th and 90th percentile (below reliable detection but above floor).
3. Freeze all amplitudes before running evaluation reps.
4. Do not recalibrate on evaluation data.

### Preprocessing

Applied to each trajectory before computing features:

1. Winsorize at 0.5% and 99.5%.
2. Robust standardize: z_t = (x_t − median(x)) / (1.4826 · MAD(x)).
3. First differences: dx_t = z_t − z_{t−1}.

### Feature vector

Computed on the preprocessed trajectory z_t:

| Feature | Definition |
|---|---|
| S | Theil-Sen slope of z_t vs t, scaled by MAD(dx) |
| Z_MK | Mann-Kendall trend z-score |
| ρ_env | Theil-Sen slope of log(rolling RMS envelope), window w=500 |
| E_ratio | RMS(z, last 20%) / RMS(z, first 20%) |
| G_spec | Max multitaper spectral peak / total spectral power |
| Q_spec | Peak frequency / spectral bandwidth (quality factor) |
| K_01 | Median 0-1 chaos statistic over 100 random c ∈ (π/5, 4π/5) |
| LLE | Rosenstein largest Lyapunov exponent estimate |
| PE | Normalized permutation entropy, order m=5 |

### Hierarchical decision rule

Labels are not symmetric. Divergence and oscillation can fool chaos statistics. Check simpler patterns first, assign the most conservative chaos label last.

```
1. DIVERGENT
   if |S| > q99(null) AND |Z_MK| > 2.58
   AND |median(last 10%) − median(first 10%)| > 3·MAD(dx)
   → divergent

2. CONVERGENT
   if ρ_env < q01(null)
   AND E_ratio < q05(null)
   AND |median(last 20%) − baseline| < 1·MAD(z)
   → convergent

3. OSCILLATORY
   if G_spec > q99(null) AND Q_spec > 5
   AND dominant period has ≥ 8 observed cycles in N
   → oscillatory

4. CHAOTIC
   if K_01 > 0.8
   AND LLE > q99(null)
   AND PE ∈ [0.55, 0.95]
   AND not already classified above
   → chaotic

5. OTHERWISE → null
```

All q_XX(null) thresholds computed from 1000 null-calibration replications (seeds 80000–80999). Evaluation replications (seeds 99999+) are disjoint.

### Baselines

For each condition, run the same hierarchical classifier on:
1. Each individual stream's log(e_t) — 5 individual classifications per rep
2. Composed log-e-value — the method under test
3. Standardized raw-data sum — z-score each stream, sum, apply same classifier
4. Best individual stream — the stream with highest detection statistic, per rep

### Predictions

1. **Divergence:** Individual Mann-Kendall detection <50% of reps. Composed >90%.
2. **Convergence:** Individual envelope-decay detection <50%. Composed >85%.
3. **Oscillation:** Replicates V2. Individual <35%. Composed >90%.
4. **Chaos:** Individual 0-1 test ambiguous (K ∈ 0.3–0.7). Composed K > 0.8 in >70% of reps. This is the weakest prediction — chaotic forcing may not survive noise compression.
5. **Null:** Per-method false positive rate <5%. Overall misclassification rate <10%.
6. **Confusion matrix:** Composed classifier achieves macro-F1 > 0.8 across all five labels. Standardized sum achieves lower macro-F1 (replicating V2's Fisher information advantage).

### Falsification

- If composed macro-F1 < 0.6, the classifier does not work on composed trajectories.
- If standardized-sum macro-F1 equals or exceeds composed macro-F1, the Fisher information advantage from V2 does not generalize.
- If chaos detection fails entirely (K_01 indeterminate on composed, <50% correct label), report as a negative result for that bin. The other three bins can still succeed independently.
- If the classifier assigns the wrong label to >20% of reps for any non-null condition, that bin's generator or detector needs revision.

### Null calibration

1000 replications under global null (seeds 80000–80999). Compute per-feature null distributions. Set thresholds at the quantiles specified in the decision rule. Report null false positive rates per branch of the decision tree.

### What this does not test

- **Multi-round loop:** Classify → generate hypothesis → test again. Deferred.
- **Correlated streams:** Independence assumed, same as V2.
- **Mixed regimes / transitions:** Systems that change bin mid-stream (e.g., oscillation → chaos via bifurcation).
- **Intrinsic dynamics:** The generators inject external forcing, not autonomous dynamical systems. The classifier detects forcing patterns, not the generating mechanism.

### Experiment parameters

- N = 10,000 observations per stream per replication
- K = 5 streams
- 1000 null-calibration reps (seeds 80000–80999)
- 100 amplitude-calibration reps (seeds 90000–90099)
- 100 evaluation reps per condition (seeds 99999+, 5 conditions)
- Logistic map: r = 3.9, z_0 varied per rep, 1000-step burn-in, standardized
- Decay time constant: τ = 2000
- Oscillation period: T = 500
- Per-stream amplitudes: frozen after calibration phase

### Outputs

- `results/tables/fourbin.csv` — per-rep: condition, signal type, feature values, assigned label
- `results/tables/fourbin_confusion.csv` — 5×5 confusion matrices for each signal type
- `results/plots/fourbin_confusion.png` — heatmap confusion matrices (composed vs individual vs standardized)
- `results/plots/fourbin_features.png` — feature distributions by condition (violin plots)
- `results/plots/fourbin_01test.png` — K_01 distributions by condition

### References

Abarbanel, H. D. I. (1996). *Analysis of Observed Chaotic Data*. Springer.

Falconer, I., Gottwald, G. A., Melbourne, I., & Wormnes, K. (2007). "Application of the 0-1 test for chaos to experimental data." *SIAM J. Appl. Dyn. Syst.*, 6(2), 395–402.

Gottwald, G. A., & Melbourne, I. (2004). "A new test for chaos in deterministic systems." *Proc. R. Soc. Lond. A*, 460, 603–611.

Gottwald, G. A., & Melbourne, I. (2009). "On the implementation of the 0-1 test for chaos." *SIAM J. Appl. Dyn. Syst.*, 8(1), 129–145.

Grünwald, P., de Heide, R., & Koolen, W. M. (2024). "Safe testing." *J. R. Stat. Soc. B*, 86(5), 1091–1128.

Hegger, R., Kantz, H., & Schreiber, T. (1999). "Practical implementation of nonlinear time series methods: The TISEAN package." *Chaos*, 9(2), 413–435.

Kantz, H. (1994). "A robust method to estimate the maximal Lyapunov exponent of a time series." *Physics Letters A*, 185(1), 77–87.

Kantz, H., & Schreiber, T. (2004). *Nonlinear Time Series Analysis*, 2nd ed. Cambridge University Press.

Mann, H. B. (1945). "Nonparametric tests against trend." *Econometrica*, 13(3), 245–259.

Ramdas, A., Grünwald, P., Vovk, V., & Shafer, G. (2023). "Game-theoretic statistics and safe anytime-valid inference." *Statist. Sci.*, 38(4), 576–601.

Rosenstein, M. T., Collins, J. J., & De Luca, C. J. (1993). "A practical method for calculating largest Lyapunov exponents from small data sets." *Physica D*, 65(1–2), 117–134.

Strogatz, S. H. (2014). *Nonlinear Dynamics and Chaos*, 2nd ed. Westview Press.

Vovk, V., & Wang, R. (2021). "E-values: Calibration, combination, and applications." *Ann. Statist.*, 49(3), 1736–1754.

Wolf, A., Swift, J. B., Swinney, H. L., & Vastano, J. A. (1985). "Determining Lyapunov exponents from a time series." *Physica D*, 16(3), 285–317.

Zanin, M., & Olivares, F. (2021). "Ordinal patterns-based methodologies for distinguishing chaos from noise in discrete time series." *Commun. Phys.*, 4, 190.
