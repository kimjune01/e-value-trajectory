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
