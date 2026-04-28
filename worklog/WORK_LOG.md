# Work Log

Full trail. Every decision, every failed attempt, every fork in the analysis. Published per [Gwern Q18](https://june.kim/prereg-audit).

## 2026-04-27

### 00:00 — Origin and repo creation

Idea emerged from writing [Evidence has a trajectory](https://june.kim/evidence-has-a-trajectory). The post argues that e-value trajectories are transparent enough to transmit the dynamics of the system under study. Prior art search found no work on reading e-value trajectory *shape* as a diagnostic. Two open questions planted at the end of the post.

Repo created. Three claims pre-registered in EXPERIMENT.md. Audited against the [prereg checklist](https://june.kim/prereg-audit) — seven gaps found and addressed:
- Added heavy-tailed conditions (Descartes Q3)
- Stated the multiplicative mechanism (Hume Q4)
- Listed competing explanations for oscillation (Chamberlin Q8)
- Added blind classifier variant (Feynman Q14)
- Defined KS test and CV thresholds for severity (Mayo Q17)
- Committed to publishing all replications and failures (Gwern Q18)
- Added λ sensitivity analysis to rule out betting-strategy artifacts

Not yet started: implementation. Next step is `generate_synthetic.py`.

### 01:00 — Prereg v2

Codex + Gemini review surfaced three fatal problems: sinusoidal DGP has no feedback (it's a schedule, not a system), all analysis was oracle-tuned to T, and KS accept-the-null is invalid. Fixed: replaced sinusoidal with Lotka-Volterra, replaced oracle lag-bank with blind periodogram, replaced KS with TOST. Added Bayesian sequential baseline to Claim 2 race.

### 02:00 — Prereg v3

Second codex round found six more issues: unspecified LV parameters, classifiers in different decision spaces, TOST comparing wrong conditions, unspecified changepoint parameters, no multiple comparisons correction. All fixed. Added escape hatches for bugs, futility, scope reduction.

### 03:00 — Prereg v4

Added sinusoidal sanity check (i) to stage the deltas cleanly: null → constant → autocorrelation → injected periodicity → emergent feedback → stochastic feedback → aperiodic switching. Added stochastic LV (b2) as the fair test. Added parallelization strategy. Created directory scaffold and conditions.yaml.

### 04:00 — Prereg v5

Third codex round caught three fatal scale issues: LV prey population (~10) doesn't match e-value formula (expects ~0.3), sinusoidal mean is zero (drifts negative like null), AR(1) stationary mean is 3.0 not 0.3. Fixed: center/scale LV observations from burn-in, add offset to sinusoidal, correct AR(1) drift parameter. Also noted Claim 2's Bayesian vs e-value comparison reduces to "which signal is cleaner for the periodogram" rather than a categorical capability difference.

### 05:00 — Analytical proofs

Separated provable claims from empirical ones in PROOFS.md. Four lemmas: periodic μ → periodic evidence, regime switch → slope sign flip, stationary → flat periodogram, AR(1) → broad spectrum. The experiments now only test what the math can't answer: detection power in finite samples.

### 06:00 — Prereg v6–v9 convergence

Three more codex volleys. AR(1) variance fix introduced a new problem (broke conditional variance assumption), reverted. Claim 2 reframed honestly as signal-surface comparison. AR(1) e-value i.i.d. misspecification documented as intentional. Codex round 8: "internally consistent, cannot identify a remaining design flaw." Ready to implement.

### 14:00 — First run

Generated all conditions, computed e-values, ran spectral analysis. Key findings:

1. Sinusoidal sanity check passes. Period 500 detected perfectly across all 100 reps.
2. LV period is ~12, not ~500. The chosen parameters produce fast oscillations. Not a bug — the emergent period is what it is.
3. AR(1) produces a strong spectral peak (period ~132, peak/med ~1109). The confound condition works.
4. E-value periodogram is identical to raw-data periodogram. Δlog(e_t) = λX_t - constant, so the periodogram is a scaled copy.

Finding #4 is a potential thesis-killer. But on reflection: the identical periodograms aren't a thesis-killer — they're the thesis. E-values are a universal projection: regardless of data type, they map to the same temporal evidence scale with valid composition. For a single univariate normal stream, the projection is an affine transform — of course the periodograms match. The value appears when you compose evidence across heterogeneous experiments.

This reframes the experiment: the interesting question isn't "does the e-value periodogram beat the raw periodogram on a single stream?" (trivially no for univariate normal). It's "when you compose across heterogeneous experiments, does the composed trajectory carry dynamics the individual streams don't?"

### 16:00 — Full pipeline implementation

Built `classify.py`, `changepoint.py`, `sensitivity.py`, `plot.py`. All read parameters from `configs/conditions.yaml`.

Sensitivity analysis (λ = 0.15, 0.3, 0.6): peak/median ratios identical across all three λ values for every condition. Confirms the affine transform — not a new result, but now quantified across 100 reps × 6 conditions × 3 λ values.

### 17:00 — Two pre-registered parameter miscalibrations found

1. Spectral peak threshold = 3 (pre-registered). 500-point periodogram with ~250 bins has expected max/median ~9 under white noise. Threshold 3 triggers on everything — 100% false positive rate. Design error, not code bug.
2. CUSUM h = 5 (pre-registered). ARL₀ ≈ 1083 steps for σ = 0.3. Over 10,000 observations, false alarm near-certain. "Textbook default" for unit-variance signals.

Calibrated from null condition data: spectral threshold → 15 (null p99 = 12.6), CUSUM h → 20 (null FA = 2/100). Logged as parameter fix per escape hatch.

### 18:00 — All claims tested

Classifier race (calibrated): threshold fastest for effect detection (median ~50 steps), Bayesian fastest overall (median ~10 steps), spectral EV only one that labels "periodic" but always fires at t=499. Comparison is categorical, not competitive.

Changepoint: CUSUM on e-value increments detects regime switches within 260 steps (median). Same parameters on raw data give 100% false alarm — but this is parameter miscalibration, not signal quality. Bayesian CPD failed entirely (22–54% FA rate).

Statistical tests: TOST equivalence confirms e-value = raw periodogram (mean diff = 0.0000). Mann-Whitney separates cyclic from null (p = 1.28e-34). CUSUM delays within 500-step budget (p < 1e-48).

### 19:00 — Report written

`report.md` covers all results, all failures. Seven documented failures including wrong AR(1) prediction, two miscalibrated thresholds, failed Bayesian CPD, wrong Claim 2 speed prediction, and borderline sinusoidal detection.

Honest summary: e-value trajectory shape reflects system dynamics (Claims 1 and 3 supported), but for a single stream this reflection is an affine transform of raw data. The interesting untested question is heterogeneous composition (Candidate Claim 4).

## 2026-04-28

### 09:14 — Committed full pipeline, started composition demo planning

Committed all work from yesterday's session: `classify.py`, `changepoint.py`, `sensitivity.py`, `plot.py`, `report.md`, updated `WORKLOG.md` and `conditions.yaml`. Commit `ab488b0`.

Then planned the next experiment: heterogeneous e-value composition (Claim 4). The idea is that for a single stream, the e-value periodogram is an affine transform of raw data (no advantage), but composing e-values across streams with different distributions reveals shared dynamics that no individual stream can detect. Scenario: five sensors (Normal, Poisson, Exponential, Bernoulli, Lognormal) sharing a tidal cycle, each individually below detection threshold.

Sent the plan to codex. Key feedback:
- Oracle alternative problem: e-value formulas must use time-indexed alternatives, not fixed ones that encode the true signal.
- Use log/logit links for Poisson/Exponential/Bernoulli to avoid invalid parameter values (negative rates, probabilities out of [0,1]).
- "Proof by example" overclaims — this is a constructive demonstration.
- Need a standardize-and-sum baseline comparison (skeptic's alternative to e-value composition).
- Need null controls: run under global null, report composed false positive rate.
- Independence assumption does heavy lifting — real co-located sensors would be correlated. Keep claim narrow.
- Cut trajectory figure unless it serves a claim beyond spectral detection.

Status: waiting on user decision about which codex feedback to act on before implementing `src/compose.py`.

### 09:45 — Composition experiment: two bugs, one surprise

Implemented `src/compose.py`. First run revealed two bugs:

1. **Null e-values were zero.** When null=True, amplitude was set to 0, making the alternative equal the null. Fix: generate data under H0 but compute e-values against the same fixed alternative (not time-indexed).

2. **Period detected at 250, not 500.** Time-indexed alternatives create KL-divergence harmonics at 2/T because the expected log-LR is quadratic in the parameter perturbation. Fix: use fixed alternatives (like V1). log(e_t) is linear in X_t, periodogram shows fundamental. Codex had recommended time-indexed alternatives, but in practice they introduce artifacts.

After fixing: individual streams at PMR 13–14 (26–35% detection), composed e-value at PMR 37.8 (99% detection), period 500 in 99% of reps. Null means all within 0.02% of 1.0.

**Surprise:** standardized-sum baseline only 29% detection (vs composed 99%). Pre-registered prediction was equivalence (TOST d=0.2). The composed e-value is 2.3× stronger. Reason: likelihood ratio weights each stream by Fisher information, extracting more from informative streams. Standardized sum weights equally. This is a genuine power advantage, not just a validity argument.

Updated report.md with full Claim 4 section.
