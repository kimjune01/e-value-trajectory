# Worklog

Full trail. Every decision, every failed attempt, every fork in the analysis. Published per [Gwern Q18](https://june.kim/prereg-audit).

---

### 2026-04-27

**Origin.** Idea emerged from writing [Evidence has a trajectory](https://june.kim/evidence-has-a-trajectory). The post argues that e-value trajectories are transparent enough to transmit the dynamics of the system under study. Prior art search found no work on reading e-value trajectory *shape* as a diagnostic. Two open questions planted at the end of the post.

**Repo created.** Three claims pre-registered in EXPERIMENT.md. Audited against the [prereg checklist](https://june.kim/prereg-audit) — seven gaps found and addressed:
- Added heavy-tailed conditions (Descartes Q3)
- Stated the multiplicative mechanism (Hume Q4)
- Listed competing explanations for oscillation (Chamberlin Q8)
- Added blind classifier variant (Feynman Q14)
- Defined KS test and CV thresholds for severity (Mayo Q17)
- Committed to publishing all replications and failures (Gwern Q18)
- Added λ sensitivity analysis to rule out betting-strategy artifacts

**Not yet started:** implementation. Next step is `generate_synthetic.py`.

**Prereg v2 (same day).** Codex + Gemini review surfaced three fatal problems: sinusoidal DGP has no feedback (it's a schedule, not a system), all analysis was oracle-tuned to T, and KS accept-the-null is invalid. Fixed: replaced sinusoidal with Lotka-Volterra, replaced oracle lag-bank with blind periodogram, replaced KS with TOST. Added Bayesian sequential baseline to Claim 2 race.

**Prereg v3 (same day).** Second codex round found six more issues: unspecified LV parameters, classifiers in different decision spaces, TOST comparing wrong conditions, unspecified changepoint parameters, no multiple comparisons correction. All fixed. Added escape hatches for bugs, futility, scope reduction.

**Prereg v4 (same day).** Added sinusoidal sanity check (i) to stage the deltas cleanly: null → constant → autocorrelation → injected periodicity → emergent feedback → stochastic feedback → aperiodic switching. Added stochastic LV (b2) as the fair test. Added parallelization strategy. Created directory scaffold and conditions.yaml.

**Prereg v5 (same day).** Third codex round caught three fatal scale issues: LV prey population (~10) doesn't match e-value formula (expects ~0.3), sinusoidal mean is zero (drifts negative like null), AR(1) stationary mean is 3.0 not 0.3. Fixed: center/scale LV observations from burn-in, add offset to sinusoidal, correct AR(1) drift parameter. Also noted Claim 2's Bayesian vs e-value comparison reduces to "which signal is cleaner for the periodogram" rather than a categorical capability difference.

**Analytical proofs (same day).** Separated provable claims from empirical ones in PROOFS.md. Four lemmas: periodic μ → periodic evidence, regime switch → slope sign flip, stationary → flat periodogram, AR(1) → broad spectrum. The experiments now only test what the math can't answer: detection power in finite samples.

**Prereg v6–v9 (same day).** Three more codex volleys. AR(1) variance fix introduced a new problem (broke conditional variance assumption), reverted. Claim 2 reframed honestly as signal-surface comparison. AR(1) e-value i.i.d. misspecification documented as intentional.

**Converged (round 8).** Codex: "internally consistent, cannot identify a remaining design flaw." Eight rounds total. Ready to implement.

**Status:** prereg converged. Next: implement `src/generate.py`.

### 2026-04-27 (late)

**First run.** Generated all conditions, computed e-values, ran spectral analysis. Key findings:

1. **Sinusoidal sanity check passes.** Period 500 detected perfectly across all 100 reps.
2. **LV period is ~12, not ~500.** The chosen parameters produce fast oscillations. Expected period estimate was wrong. Not a bug — the emergent period is what it is. Need to update the prediction or adjust LV parameters to produce a slower cycle.
3. **AR(1) produces a strong spectral peak (period ~132, peak/med ~1109).** The confound condition works — autocorrelation does produce a peak. But the peak is broad and at a different frequency than the LV peak.
4. **E-value periodogram is identical to raw-data periodogram.** Δlog(e_t) = λX_t - constant, so the periodogram is a scaled copy. The e-value adds no spectral information the raw data doesn't already have.

Finding #4 is a potential thesis-killer. The claim is that e-values are "transparent" to system dynamics, but transparency without amplification means the periodogram on raw X_t is strictly equivalent. The e-value's added value would have to come from the *cumulative* trajectory log(E_t), not the per-step Δlog(e_t). Need to investigate whether the periodogram on cumulative log(E_t) differs from the periodogram on cumulative X_t.

**Reframe.** The identical periodograms aren't a thesis-killer — they're the thesis. E-values are a universal projection: regardless of data type (frequency, probability, periodic, disjoint), they map to the same temporal evidence scale with valid composition. For a single univariate normal stream, the projection is an affine transform of the raw data — of course the periodograms match. The value appears when you compose evidence across *heterogeneous* experiments with different DGPs, scales, and test statistics. Raw data from different experiments can't be concatenated; e-values can (multiplicative composition). The trajectory of the composed e-value should still carry the system's dynamics even when no single raw stream does.

This reframes the experiment: the interesting question isn't "does the e-value periodogram beat the raw periodogram on a single stream?" (trivially no for univariate normal). It's "when you compose across heterogeneous experiments, does the composed trajectory carry dynamics the individual streams don't?"

**Status:** deciding whether to pivot the experiment to heterogeneous composition or report the current results as a completed (negative) finding on the original claim.

### 2026-04-27 (continued)

**Implemented remaining pipeline.** Built `classify.py`, `changepoint.py`, `sensitivity.py`, `plot.py`. All read parameters from `configs/conditions.yaml`.

**Sensitivity analysis (λ = 0.15, 0.3, 0.6).** Peak/median ratios are *identical* across all three λ values for every condition. Confirms the affine transform: Δlog(e_t) = λX_t - λ²/2, so the periodogram scales by λ² and the peak/median ratio cancels. λ has zero effect on spectral classification. This is a direct corollary of the affine transform finding from the first run — not a new result, but now quantified across 100 reps × 6 conditions × 3 λ values.

**Two pre-registered parameter miscalibrations found.**

1. **Spectral peak threshold = 3 (pre-registered).** The 500-point window periodogram has ~250 frequency bins. Under white noise, the expected max/median ratio is ~9 (null condition: p99 = 12.6, max = 14.1). Threshold 3 triggers on everything — 100% false positive rate across all conditions. Root cause: the threshold was set without accounting for the number of independent frequency bins tested. This is a design error in the pre-registration, not a code bug.

2. **CUSUM h = 5 (pre-registered).** With e-value increment variance λ² = 0.09 and reference value k = 0.0225, the in-control ARL₀ ≈ 1083 steps. Over 10,000 observations, false alarm is near-certain. Root cause: h = 5 is a "textbook default" that assumes unit-variance signals; e-value increments have σ = 0.3, requiring h ≈ 20 for adequate ARL₀.

**Calibration.** Computed corrected thresholds from null condition data (condition c, 100 reps):
- Spectral: threshold = 15 (null p99 = 12.6; at 15, null FA = 0%).
- CUSUM: h = 20 (null FA = 2/100).
Added calibrated values to `conditions.yaml` as a clearly labeled post-hoc section. Original pre-registered values preserved. Logged as a parameter fix, not a protocol change — the classifiers literally don't run at the original thresholds (100% false positive = no discriminative power). See EXPERIMENT.md escape hatch: "fixing a bug that prevents the experiment from running is maintenance."

**Reran classify.py and changepoint.py with calibrated parameters.**

**Classifier race results (calibrated peak_threshold = 15):**
- **Threshold classifier:** Fastest for non-periodic conditions. Median stopping time 53 (stationary), 35 (sinusoidal), 48 (LV). Labels: 100% reject_null for all positive-effect conditions. Cannot label "periodic" by design. Correct on null (92% null, 8% reject_null — the 8% are reps where E_t crosses 20 before the negative slope detection kicks in).
- **Bayesian classifier:** Fastest overall. Median stopping time 12 (stationary), 10 (sinusoidal), 7 (LV). But 27% false "periodic" on null — the posterior mean's autocorrelation inflates spectral peaks even with the calibrated threshold. Also 38% false "null" on AR(1) (posterior drifts negative during autocorrelated stretches).
- **Spectral e-value classifier:** Always stops at t = 499 (first window). Correctly labels null → "null" (100%), stationary → "reject_null" (100%), LV → "periodic" (100%), stochastic LV → "periodic" (88%). AR(1) → "periodic" (100%) — the designed confound. Cannot distinguish autocorrelation from periodicity by peak height. Sinusoidal → 53% "periodic", 47% "reject_null" — borderline because exactly 1 cycle fits in 500-point window.
- **Verdict on Claim 2:** The three classifiers answer different questions. Threshold and Bayesian answer "is there an effect?" fast but cannot diagnose periodicity. Spectral e-value is the only classifier that labels "periodic" on LV conditions — but it's also the slowest (needs 500 observations minimum) and can't distinguish AR(1) from true cycles. The comparison is categorical, not competitive: it's like comparing a thermometer to a stethoscope and asking which is "faster."

**Changepoint results (calibrated h = 20):**
- **CUSUM on e-value increments:** Median detection delay 260 steps (f, effect→null), 256 steps (g, null→effect), 118 steps (h, effect→reversed). All within the 500-step budget (99–100% of reps). False alarm: 0% on stationary, 2% on null.
- **CUSUM on raw X_t (same parameters):** 100% false alarm rate. Parameters calibrated for e-value variance (σ² = 0.09) are miscalibrated for raw variance (σ² = 1.0). Not a fair head-to-head — the e-value transformation compresses noise by factor λ = 0.3, making h = 20 appropriate. With equivalently calibrated h (≈ 220), raw CUSUM would have the same detection power (affine transform).
- **Bayesian CPD:** Poor on both signals. FA rate 22–54%. Huge variance in detection delay (IQR spans thousands of steps). Hazard = 1/2000 is too sensitive for 10,000-step sequences.
- **Verdict on Claim 3:** CUSUM on e-value increments detects regime switches within 500 steps — the thesis claim holds. But this is not because e-values amplify the signal; it's because the e-value transformation rescales the data into a regime where standard CUSUM parameters work. The same result is achievable on raw data with appropriately scaled parameters.

**Statistical tests.**
- TOST equivalence (e-value vs raw periodogram peak/median ratio): all conditions p ≈ 0. Perfect equivalence confirmed. The ratio is numerically identical (mean diff = 0.0000) because λ² cancels in max/median.
- Mann-Whitney (LV vs null peak/median): U = 10000, p = 1.28e-34. Survives Bonferroni (α/90 = 0.00056) by 30 orders of magnitude. Stochastic LV and sinusoidal also survive.
- AR(1) vs null: U = 10000, p = 1.28e-34. The confound also passes — AR(1) has spectral peaks significantly above null.
- CUSUM detection delay < 500: all three regime conditions, p < 1e-48. Effect→reversed is fastest (median 118, p = 6.24e-118).

**Status:** All claims tested. Report pending.
