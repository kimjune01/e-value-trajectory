# Pre-registration: Four-bin classification on composed e-value trajectories

Third pre-registration. V1 tested single-stream diagnostics (Claims 1–3). V2 tested heterogeneous composition for periodicity detection (Claim 4). This document covers Claim 5: full four-bin dynamical classification on composed trajectories.

## What prior work establishes (no experiment needed)

### Single-stream e-value trajectories inherit dynamical classification from raw data

For a single univariate stream with fixed alternative, log(e_t) = cX_t + d (affine transform). Dynamical classification is invariant under invertible affine transforms: Lyapunov exponents are unchanged (log-distance growth adds only a constant from |c|), the 0-1 test is invariant under scaling and translation, periodogram peak ratios cancel the scale factor, and trend tests are sign-preserved for c > 0 (Kantz & Schreiber, 2004, §3.1–3.3). No experiment needed — classifying log(e_t) gives the same result as classifying X_t.

### The four bins and their detection methods are established

Each regime has canonical detection methods from the nonlinear time series analysis literature:

**Convergence (stable fixed point).** Negative largest Lyapunov exponent. Exponentially decaying autocorrelation. Perturbations shrink.
- Rosenstein, M. T., Collins, J. J., & De Luca, C. J. (1993). "A practical method for calculating largest Lyapunov exponents from small data sets." *Physica D*, 65(1–2), 117–134.
- Kantz, H. (1994). "A robust method to estimate the maximal Lyapunov exponent of a time series." *Physics Letters A*, 185(1), 77–87.

**Divergence (unstable/trending).** Positive trend, non-stationary. Unit root or explosive root.
- Mann, H. B. (1945). "Nonparametric tests against trend." *Econometrica*, 13(3), 245–259.
- Dickey, D. A., & Fuller, W. A. (1979). "Distribution of the estimators for autoregressive time series with a unit root." *JASA*, 74(366), 427–431.

**Oscillation (limit cycle).** Spectral peaks, periodic autocorrelation. Tested in V1 (Claim 1) and V2 (Claim 4).
- Periodogram analysis. See V1 report for implementation.

**Chaos (strange attractor).** Positive largest Lyapunov exponent. Bounded but aperiodic. Broadband spectrum.
- Gottwald, G. A., & Melbourne, I. (2004). "A new test for chaos in deterministic systems." *Proc. R. Soc. Lond. A*, 460, 603–611.
- Gottwald, G. A., & Melbourne, I. (2009). "On the implementation of the 0-1 test for chaos." *SIAM J. Appl. Dyn. Syst.*, 8(1), 129–145.
- Wolf, A., Swift, J. B., Swinney, H. L., & Vastano, J. A. (1985). "Determining Lyapunov exponents from a time series." *Physica D*, 16(3), 285–317.

**General references:**
- Kantz, H., & Schreiber, T. (2004). *Nonlinear Time Series Analysis*, 2nd ed. Cambridge University Press.
- Strogatz, S. H. (2014). *Nonlinear Dynamics and Chaos*, 2nd ed. Westview Press.
- Abarbanel, H. D. I. (1996). *Analysis of Observed Chaotic Data*. Springer.
- Hegger, R., Kantz, H., & Schreiber, T. (1999). "Practical implementation of nonlinear time series methods: The TISEAN package." *Chaos*, 9(2), 413–435.

No single paper or textbook presents the four-bin classification as a unified diagnostic pipeline. The pipeline is assembled from the above components.

### E-value composition is valid across heterogeneous experiments

The product of independent e-values is a valid e-value (Vovk & Wang, 2021; Grünwald et al., 2024). V2 demonstrated that composed e-value trajectories recover shared periodic dynamics with 3.4× the power of standardized-sum aggregation, via Fisher information weighting.

- Vovk, V., & Wang, R. (2021). "E-values: Calibration, combination, and applications." *Ann. Statist.*, 49(3), 1736–1754.
- Grünwald, P., de Heide, R., & Koolen, W. M. (2024). "Safe testing." *J. R. Stat. Soc. B*, 86(5), 1091–1128.
- Ramdas, A., Grünwald, P., Vovk, V., & Shafer, G. (2023). "Game-theoretic statistics and safe anytime-valid inference." *Statist. Sci.*, 38(4), 576–601.

## What this experiment tests

The affine transform argument guarantees that single-stream classification is invariant. But composed e-value trajectories are **not** an affine transform of any single raw stream — they are sums of nonlinearly-transformed heterogeneous signals. V2 showed composition works for oscillation detection (periodograms). This experiment tests whether it works for all four bins.

## Claim 5: Four-bin classification on composed e-value trajectories

**Question:** When K heterogeneous experiments share a common dynamical regime, does the composed e-value trajectory classify correctly into convergent, divergent, oscillatory, or chaotic — and does composition improve detection over individual streams?

### Hypotheses

- **H0:** All K streams follow their null distributions (no shared dynamics).
- **H1-converge:** Each stream has a shared perturbation that decays exponentially (mean reverts to null).
- **H1-diverge:** Each stream has a shared trend (mean drifts linearly).
- **H1-oscillate:** Each stream has a shared periodic forcing (reuse V2 setup, period T=500).
- **H1-chaos:** Each stream is driven by a shared chaotic forcing (logistic map, Lorenz x-component, or Rössler).

### Design

K=5 streams, same distributions as V2 (Normal, Poisson, Exponential, Bernoulli, Lognormal). Same fixed alternatives, same null parameters. N=10,000 observations. 100 replications per condition.

The shared forcing differs by regime:

| Regime | Forcing on stream k's mean | Parameters |
|---|---|---|
| Null | None (constant at null) | — |
| Convergence | A_k · exp(-t/τ) | τ = 2000 (slow decay) |
| Divergence | A_k · t/N | Linear ramp, amplitude A_k |
| Oscillation | A_k · sin(2πt/T) | T=500, reuse V2 |
| Chaos | A_k · z_t (shared logistic map) | z_{t+1} = r·z_t·(1−z_t), r=3.9 |

Amplitudes A_k calibrated per stream so individual detection is unreliable (below threshold), same methodology as V2.

For chaos: the logistic map at r=3.9 produces bounded aperiodic dynamics with positive Lyapunov exponent (ln 2 ≈ 0.693 at r=4; slightly lower at r=3.9). The shared z_t sequence drives all five streams' means, creating correlated chaotic modulation. Each stream adds its own distribution-specific noise.

### Classifiers

Per-stream and composed, apply all four detection methods:

| Method | Target bin | Implementation |
|---|---|---|
| Largest Lyapunov exponent (Rosenstein) | Chaos (+) vs convergence (−) | Embedding dimension via false nearest neighbors, lag via mutual information |
| 0-1 test (Gottwald & Melbourne) | Chaos (K≈1) vs regular (K≈0) | Implementation per Gottwald & Melbourne (2009) |
| Periodogram peak/median ratio | Oscillation (high) vs noise (low) | Same as V1/V2, threshold from null calibration |
| Mann-Kendall trend test | Divergence (significant trend) vs stationary (no trend) | Two-sided, α=0.05 |

### Composition

Same as V2: log(E_composed,t) = Σ_k log(e_t^(k)). Also compute standardized-sum baseline for comparison.

### Predictions

1. **Convergence:** Individual streams' Lyapunov exponents are indistinguishable from null (weak signal). Composed trajectory has a more negative Lyapunov exponent (or lower variance in the estimate) than individual streams.
2. **Divergence:** Individual Mann-Kendall tests detect the trend in <50% of reps. Composed Mann-Kendall detects in >90%.
3. **Oscillation:** Replicates V2 result. Individual periodogram detection <35%. Composed >90%.
4. **Chaos:** Individual 0-1 test scores are ambiguous (K between 0.3–0.7). Composed 0-1 test gives K > 0.9 in >80% of reps. Composed Lyapunov exponent is positive and distinguishable from null.
5. **Null:** All classifiers on composed trajectory return null/no-signal in >95% of reps (false positive control).

### Falsification

- If composed trajectories fail to classify correctly in ≥2 of the 4 bins, composition does not generalize beyond oscillation detection.
- If standardized sums match or exceed composed e-values across all bins, the Fisher information advantage from V2 was specific to periodicity, not general.
- If chaos detection fails entirely on composed trajectories (K near 0.5, Lyapunov exponent indeterminate), the chaotic forcing may not survive the e-value transformation's noise compression.

### Null controls

Run all five streams under global null (no shared forcing), 100 replications. Report per-method false positive rates for individual and composed signals. Compute empirical null thresholds for each classifier.

### What this does not test

- **Multi-round loop:** Classify → generate hypothesis → run next experiment → classify again. Whether this converges on causal structure is a separate experiment.
- **Correlated streams:** Independence assumed, same as V2.
- **Unknown regime:** The experiment assigns a known regime to each condition. A blind classifier that determines the regime from the trajectory without knowing the ground truth is deferred.
- **Mixed regimes:** Systems that transition between bins (e.g., oscillation → chaos via bifurcation) are not tested.

### Experiment parameters

- N = 10,000 observations per stream per replication
- K = 5 streams
- 100 replications per condition (5 conditions: null, converge, diverge, oscillate, chaos)
- 100 null replications for threshold calibration
- Base seed: 42
- Logistic map: r = 3.9, z_0 = 0.4
- Decay time constant: τ = 2000
- Oscillation period: T = 500
- Per-stream amplitudes: to be calibrated before running (target: individual detection <35%)

### Outputs

- `results/tables/fourbin.csv` — per-replication summary: condition, rep, signal (individual/composed/standardized), method, test statistic, classification label
- `results/plots/fourbin_detection.png` — detection rates by condition × method × signal type
- `results/plots/fourbin_lyapunov.png` — Lyapunov exponent distributions by condition (individual vs composed)
- `results/plots/fourbin_01test.png` — 0-1 test K distributions by condition

### References

Abarbanel, H. D. I. (1996). *Analysis of Observed Chaotic Data*. Springer.

Dickey, D. A., & Fuller, W. A. (1979). "Distribution of the estimators for autoregressive time series with a unit root." *JASA*, 74(366), 427–431.

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
