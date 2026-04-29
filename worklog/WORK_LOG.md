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

### 10:30 — Blog post edits, prior art search, prereg V3

Sharpened the blog post's e-value argument: the contrast is between lossy compression (p-values force a scalar because peeking invalidates the guarantee) and keeping the time series (anytime validity removes that constraint). Mathematical property, not empirical. Added Pearl connection: he said stop observing and intervene; same logic one level up — stop compressing and look at the data. Four bins are sensemaking modes, not conclusions.

Prior art search (web + codex): no prior work on classifying e-value trajectories via dynamical systems methods. The 4-bin classification on general time series is established (Kantz & Schreiber 2004, Strogatz 2014, Gottwald & Melbourne 2004/2009, Rosenstein et al. 1993). Affine invariance means single-stream classification needs no experiment — just cite.

Wrote PREREGISTRATION_V3.md: four-bin classification on composed e-value trajectories. The gap is that composed signals are not affine transforms of any single stream, so composition needs testing for convergence (Lyapunov), divergence (Mann-Kendall), and chaos (0-1 test) — not just oscillation (periodogram, already done in V2). Commit `68f6c2e`.

### 11:00 — Prereg V3 codex review and revision

Codex flagged V3 as overclaiming — "four-bin dynamical classification" when the generators inject external forcing, not autonomous dynamics. Also: no multiclass decision rule (just four independent detectors), Lyapunov fragile for convergence, chaos forcing uncentered, amplitude calibration too loose, 100 null reps too few.

Asked codex for practical classification research. Got a concrete hierarchical decision tree: diverge → converge → oscillate → chaos → null. Key per-bin statistics: Theil-Sen slope + Mann-Kendall for divergence, log rolling RMS envelope decay for convergence, multitaper spectral peak for oscillation, 0-1 test + LLE + permutation entropy for chaos (most conservative label, checked last). IAAFT surrogates for null calibration.

Revised prereg: reframed as forcing-pattern classification, added the hierarchical rule with 9-feature vector, centered chaos forcing (standardized logistic map, varied z_0 per rep), locked amplitude calibration on separate seeds, bumped null reps to 1000, added 5×5 confusion matrix with macro-F1 > 0.8 as success criterion. Commit `e51a6b0`.

Also edited the blog post: sharpened the e-value argument to assert temporal preservation as mathematical certainty. P-values force lossy compression because peeking breaks the guarantee. E-values remove that constraint. Pearl said intervene; same logic one level up — stop compressing and look at the data. Four bins are sensemaking modes.

### 11:30 — Kill-condition decision tree, pushed blog post, started V3 experiment

Restructured the classifier around kill conditions (proof-manual pattern). Key fix: convergence vs divergence resolved by curvature ratio (RSS exponential / RSS linear), not branch ordering. Mann-Kendall detects monotone trend; curvature disambiguates. Renamed chaos → aperiodic forcing.

Pushed blog post edits and deployed to june.kim. Commit `acaf106`.

Implemented `src/fourbin.py`: 1000 null calibration reps, amplitude sweep on separate seeds, 100 eval reps × 5 conditions, kill-condition logging, confusion matrices.

### 12:00 — Four-bin experiment: three iterations to convergence

First run: Theil-Sen slope was 1.8s per trajectory (O(n²) pairwise slopes). Subsampled to 1000 points. Also subsampled 0-1 test and Mann-Kendall. Feature computation dropped from 2.2s to 0.49s.

Run 1 (amp=0.5): divergence and convergence undetected (all null). Trend too weak — Mann-Kendall has low power on per-step increments with small slopes.

Run 2 (amp=1.5 for diverge/converge): divergent 100%, convergent 100%, oscillatory 100%. BUT aperiodic → 100% oscillatory. The logistic map's period-2 orbit structure creates a high-Q spectral peak near Nyquist (period 2.65), triggering the oscillation detector before the aperiodic test.

Run 3 (added period > 10 filter): all five bins correct. Composed macro-F1 = 0.996. Standardized sum F1 = 0.478. Individual majority vote F1 = 0.279. The composition advantage generalizes across all four forcing patterns.

The curvature test (RSS exponential / RSS linear) perfectly separates convergence from divergence. The kill-condition framework works: monotone → curvature → periodicity → aperiodic → null. Each failure mode names the next test.

### 14:00 — Blog post: The Hypothesis Graph

Published "The Hypothesis Graph" as sequel to "Evidence has a trajectory." Core argument: experiments are nodes, kill conditions generate edges, the frontier is where belief meets the unknown. Three SVGs: mechanic diagnostic tree (linear, 3 bins), web server under load (nonlinear, all 4 bins), p-value vs e-value pipeline comparison. Added MECE proof via Milnor attractor classification. Proof by contradiction that the hypothesis generation algorithm must exist. Closer: "If we can poke it, we know how to know."

Multiple copyedit passes revealed pipeline issues: em-dash pass was missing (18 dashes slipped through), prosody pass was missing, parallel scan caused silent conflict resolution. Fixed the copyedit skill: sequential pipeline, apply-first workflow.

### 17:00 — V3 rerun with 1000 null reps: F1 = 1.000

Proper run with 1000 null calibration reps. Composed classifier: perfect classification across all five labels. The two false positives from the 200-rep run disappeared with tighter thresholds. Codex adversarial review flagged: manual amplitudes are a confound, period > 10 filter is post-hoc, synthetic archetypes are too clean, need amplitude sensitivity curves and out-of-generator stress tests.

### 17:30 — Prereg V4: robustness and sensitivity analysis

Pre-registered five stress tests to turn F1 = 1.000 from exploratory to confirmatory: amplitude sensitivity curves, decision-tree ablation, mixed dynamics, degraded conditions (autocorrelated nulls, correlated streams, missing data, misspecified distributions, nonstationary baselines), period filter sensitivity sweep. Added bootstrap CIs. Audited against the prereg checklist: 18/20 questions answered, two skipped with reasons (Q12 alternative bin structures, Q15 no causal claims).

Connected Hume to the blog post: sciences that can't perturb are structurally unable to learn causality. The hypothesis graph only works where you can poke.

### 19:00 — Blog post published: The Hypothesis Graph

"The Hypothesis Graph" published as sequel to "Evidence has a trajectory." Core claim: experiments are nodes, kill conditions generate edges, the frontier is where belief meets the unknown. Proof by contradiction that the hypothesis generation algorithm must exist. Three SVGs (mechanic, web server, p-value vs e-value). MECE partition via Milnor. Hume cited for perturbation-access limit. Closer: "If we can poke it, we know how to know."

Multiple copyedit passes fixed pipeline issues: em-dash pass was missing (18→0), prosody pass was missing, parallel scan caused silent conflict resolution. Rewrote copyedit skill: sequential pipeline, apply-first workflow. Codex cross-review of both posts found contradictions ("don't reshape" vs Fisher weighting, "cleanly" vs partial knowledge). All fixed.

### 20:00 — V3 rerun F1=1.000, V4 prereg started

1000-null-rep rerun of V3: perfect classification. Codex adversarial review flagged manual amplitudes, post-hoc period filter, synthetic cleanliness. V4 prereg written: amplitude sensitivity, ablation, mixed dynamics, degraded conditions, period filter sweep. Three codex rounds to converge the prereg.

### 21:00 — V4 implementation: blind-blind-merge

Round 1: wrote implementation, sent codex the prereg to write its own blind. Diffed. Five disagreements: missing data (zero-pad vs interpolate), correlated streams (corrupt x vs log_e), equalized difficulty (estimate vs hardcode), output schema (7 vs 13 columns), standardized sum (missing vs per-cell). Merged best of both.

Codex bug hunt found 17 bugs. Highest impact: period sweep was a no-op (classifier delegated to V3 ignoring knobs), correlated streams only affected Normal, equalized difficulty was stubbed, e-value validity not implemented.

Root cause: 20 spec ambiguities in the prereg. Codex spec-level review confirmed. Wrote operational appendix pinning every decision: classifier knobs, ablation operations, seed schedule, degradation operators, mixed dynamics formula, bootstrap method, output schema.

Round 2: blind-blind-merge with appendix. 1 disagreement (plots only, completeness not design). Down from 5. The appendix worked.

### 21:45 — V4 prereg nearing codex approval

Seven nits remaining: pinned commit SHA, stable hash (SHA-256 not Python hash), 21 steps not 20, t-scaling corrected (sqrt((df-2)/df) not sqrt(df/(df-2))), reps per-condition-per-severity, CONDITIONS order pinned to YAML, plot specs fully defined. All fixed.

Status: waiting on codex final approval of V4 prereg. Implementation ready to run once approved.

### 22:00 — V4 prereg approved by 3 independent reviewers

Codex (GPT-5.5) approved after 3 rounds. Gemini 2.5 Pro approved on first pass. Gemini 3.1 Pro Preview rejected on first pass, finding 7 bugs the other two missed:

1. Ablation logic inverted: setting features to 0 kills the class instead of skipping the test. Fix: set to always-pass values (inf, 1.0, 0.75).
2. Standardized sum needs its own null calibration (different scale from composed log_e).
3. Correlated streams formula attenuated the deterministic signal. Fix: add noise, don't scale down.
4. Supermartingale check was on single-step e_t, not cumulative E_t.
5. Bootstrap seed lacked iteration index, producing 1000 identical samples.
6. Drift formula was ambiguous (could be squared).
7. Pass/fail had a gap between 0.6 and 0.8.

All 7 fixed. Gemini 3.1 Pro re-approved. Three reviewers, three approvals.

Blind-blind-merge process: round 1 had 5 design disagreements and 17 implementation bugs. Operational appendix resolved design disagreements. Round 2 had 1 disagreement (plots only). Gemini 3.1 caught logic bugs that survived both the appendix and two other reviewers.

### 22:30 — Gemini CLI fixed

The gemini CLI uses Vertex AI, but preview models are on the Generative Language API (different endpoint). Added gem() to .zshrc but removed it per user preference. Direct python3/urllib calls to generativelanguage.googleapis.com work for gemini-3.1-pro-preview.

Status: V4 prereg approved. Implementation needs rewrite to match fixed spec (ablation logic, z_sum calibration, correlation formula, supermartingale check). Ready to implement.
