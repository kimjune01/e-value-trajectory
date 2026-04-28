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

**Status:** prereg is converging. No implementation yet. Next: one more codex volley, then implement `src/generate.py`.
