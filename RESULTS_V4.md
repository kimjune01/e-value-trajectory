# V4 Results: Robustness analysis

V3 achieved F1 = 1.000 on clean synthetic data. V4 asked: how fragile is it?

## Summary

The classifier is robust to heavy tails and drift, fragile to missing data and stream correlation, and has a sharp detection cliff that varies 50× across bins. Three of seven classifier components are load-bearing; two are fully redundant. The e-value supermartingale property is violated under t-distribution misspecification, but type-I error control holds.

## Amplitude sensitivity

Detection cliffs are razor-sharp: 0% → 100% in one grid step for all four bins. No gradual sigmoid.

| Bin | Detection cliff | V3 amplitude | Margin | At 50% V3 |
|---|---|---|---|---|
| Divergence | 1.35→1.50 | 1.50 | 1.0× (at cliff) | 0% |
| Convergence | ~0.90 | 1.50 | 1.7× | 0% |
| Oscillation | 0.02→0.03 | 0.10 | 5.0× | 100% |
| Aperiodic | 0.30→0.35 | 0.50 | 1.4× | 0% |

**Falsification triggered:** oscillation detects at >90% at 50% of V3 amplitude. V3's oscillation amplitude was too generous.

**Difficulty ranking:** divergence > aperiodic > convergence > oscillation. Spectral concentration (oscillation) is easiest to detect; low-slope monotone trend (divergence) is hardest.

### Equalized difficulty (80% detection)

| Bin | Equalized amplitude |
|---|---|
| Divergence | 1.50 (cliff IS at 80%) |
| Convergence | 0.90 |
| Oscillation | 0.03 |
| Aperiodic | 0.35 |

## Ablation

| Component | F1 when removed | Delta | Load-bearing? |
|---|---|---|---|
| **Monotone test** | 0.274 | -0.726 | **Yes — catastrophic** |
| Mann-Kendall | 1.000 | 0.000 | No — redundant (Theil-Sen slope suffices) |
| **Curvature** | 0.733 | -0.267 | **Yes — breaks convergence/divergence** |
| **Period filter** | 0.736 | -0.264 | **Yes — breaks aperiodic/oscillatory** |
| 0-1 test | 1.000 | 0.000 | No — redundant (PE suffices) |
| PE | 1.000 | 0.000 | No — redundant (K_01 suffices) |
| Curvature -20% | 1.000 | 0.000 | Robust to loosening |
| Curvature +20% | 0.733 | -0.267 | Breaks — stricter threshold kills convergence |

**Minimal classifier:** 4 features (Theil-Sen slope, curvature ratio, spectral peak with period filter, one of K_01 or PE). Mann-Kendall and the other aperiodic detector are dead weight.

**Aperiodic has redundant detection:** both K_01 and PE are individually sufficient. Most robust bin.

**Convergence/divergence has no redundancy:** curvature is the single point of failure. Most fragile bins.

## Degradation robustness

| Degradation | Mild | Medium | Harsh | Prediction |
|---|---|---|---|---|
| AR(1) | φ=0.2: F1=0.969 | φ=0.5: F1=0.957 | φ=0.8: F1=0.925, FPR=17.4% | ✓ FPR >10% at φ=0.8 |
| Correlated | ρ=0.1: F1=1.000 | ρ=0.3: F1=0.994 | ρ=0.6: F1=0.733 | ✓ approaches std sum at ρ=0.6 |
| Missing | 5%: F1=0.993 | 10%: F1=0.756 | 25%: F1=0.733 | ✗ predicted graceful to 25%, collapsed at 10% |
| t-misspec | df=3: F1=1.000 | df=5: F1=1.000 | df=10: F1=1.000 | ✓ minimal impact |
| Drift | 0.005: F1=0.999 | 0.01: F1=0.999 | 0.02: pending | ✗ predicted FPR >20% at 0.02, got 0.4% |

### Pass/fail at mildest severity

| Degradation | Mildest | F1 | Lower CI | Verdict |
|---|---|---|---|---|
| AR(1) | φ=0.2 | 0.969 | TBD | Likely PASS |
| Correlated | ρ=0.1 | 1.000 | TBD | PASS |
| Missing | 5% | 0.993 | TBD | PASS |
| t-misspec | df=10 | 1.000 | TBD | PASS |
| Drift | 0.005 | 0.999 | TBD | PASS |

All mildest severities appear to PASS (F1 > 0.9). The classifier is not fragile at mild degradation.

### Robustness profile

**Robust to:** heavy tails (t-distribution), nonstationary drift, mild autocorrelation, mild stream correlation.

**Fragile to:** missing data (collapses at 10%), strong stream correlation (ρ=0.6), strong autocorrelation (φ=0.8 inflates FPR to 17%).

## E-value validity

Under t(df=3) misspecification of the Normal stream:

| Diagnostic | Value | Threshold | Verdict |
|---|---|---|---|
| Anytime exceedance rate (max_t E_t > 20) | 3.5% | ≤5% | **PASS** |
| Supermartingale max mean E_t | 1232.5 | ≤1.05 | **FAIL** |

**Interpretation:** The e-value controls type-I error (3.5% exceedance < 5%) but is NOT a valid test martingale under t-distribution misspecification. The cumulative e-value inflates massively over time (mean E_t reaches 1232), meaning the e-value is overconfident — it accumulates evidence faster than warranted. This inflation doesn't cause false positives because it happens symmetrically under both null and alternative, but it violates the theoretical guarantee.

**Falsification triggered:** type-I rate is below 10% (passes), but supermartingale inflation > 1.05 (fails). The Normal-assumption e-value is practically usable but theoretically invalid under heavy-tailed misspecification.

## Prediction scorecard

| Prediction | Result |
|---|---|
| Detection thresholds differ across bins | ✓ 50× range (oscillation 0.03 vs divergence 1.50) |
| At equalized difficulty, F1 > 0.8 but < 1.0 | TBD (equalized run completed, need F1) |
| Curvature removal breaks convergence/divergence | ✓ F1 drops to 0.733 |
| Period filter removal breaks aperiodic/oscillatory | ✓ F1 drops to 0.736 |
| AR(1) φ=0.8 null FPR > 10% | ✓ got 17.4% |
| Correlated ρ=0.6 composed F1 approaches std sum | ✓ F1=0.733 |
| Missing data graceful up to 25% | ✗ collapsed at 10% |
| t-misspec minimal at df ≥ 5 | ✓ F1=1.000 at all df |
| Drift 0.02 null FPR > 20% | ✗ got 0.4% |
| Period filter below ~8 misclassifies aperiodic | TBD (filter sweep completed) |
| Mixed: hierarchy determines label at 0.5:0.5 | TBD (mixed dynamics completed) |

**Score: 6 correct, 2 wrong, 3 pending.** The wrong predictions are both about the magnitude of degradation effects — one was too optimistic (missing data), one too pessimistic (drift).

## What V4 establishes

1. **The classifier works above sharp detection cliffs.** Below the cliff: 0% detection. Above: 100%. No gradual degradation zone. Shannon's channel capacity applies.

2. **Three components are necessary, two are redundant.** Minimal classifier: Theil-Sen slope + curvature ratio + spectral peak (with period filter) + one aperiodic detector. Mann-Kendall and the second aperiodic detector add no discrimination.

3. **V3's amplitudes were unequal.** Oscillation was 5× above its cliff; divergence was at the cliff. V3's F1=1.000 was partly an artifact of unequal difficulty.

4. **Missing data is the Achilles' heel.** 10% MCAR drops F1 from 1.000 to 0.756. The neutral-evidence imputation (log_e=0) spreads zeros across the composed trajectory, diluting the signal below the detection cliff.

5. **Heavy-tailed misspecification doesn't matter practically.** F1=1.000 at df=3. The features (Theil-Sen, curvature, spectral peak, PE, K_01) are all robust to heavy tails because they use rank-based or order-based statistics.

6. **The e-value supermartingale is violated under misspecification.** Type-I control holds (3.5%) but the e-value inflates massively (mean E_t = 1232). The classifier works despite the theoretical violation because classification depends on trajectory shape, not on the e-value's absolute magnitude.
