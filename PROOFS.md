# Analytical results

Casual proofs that justify the experimental design. These establish what must be true in expectation. The experiments test whether these results survive finite-sample noise.

## Lemma 1: Periodic μ produces periodic evidence

**Setup.** Observations X_t ~ N(μ(t), 1) where μ(t) is periodic with period T. E-value at each step: e_t = exp(λX_t - λ²/2). Cumulative log-evidence: L_t = Σ log(e_i) = Σ (λX_i - λ²/2).

**Claim.** The expected growth rate of L_t is periodic with the same period T.

**Proof.** The expected log-evidence at step t is:

    E[log e_t] = E[λX_t - λ²/2]
               = λE[X_t] - λ²/2
               = λμ(t) - λ²/2

Since μ(t) has period T, λμ(t) has period T. Subtracting the constant λ²/2 doesn't change the period. So E[log e_t] is periodic with period T.

The expected slope of L_t in any window is the average of E[log e_t] over that window. A window that spans a full cycle averages to λμ̄ - λ²/2 where μ̄ is the mean of μ(t) over one period. A window that spans a half-cycle at the peak averages higher; at the trough, lower. The rolling slope oscillates with period T.

**What this doesn't prove.** That the oscillation is detectable above sampling noise in finite N. That's the experiment.

## Lemma 2: Regime switch produces slope sign flip

**Setup.** μ(t) = μ₁ for t < t*, μ(t) = μ₂ for t ≥ t*. Same e-value as above.

**Claim.** The expected slope of log(E_t) changes from λμ₁ - λ²/2 to λμ₂ - λ²/2 at t*.

**Proof.**

    Before t*: E[log e_t] = λμ₁ - λ²/2
    After t*:  E[log e_t] = λμ₂ - λ²/2

With λ = 0.3:

    Effect → null (μ₁=0.3, μ₂=0):   slope flips from +0.045 to -0.045
    Null → effect (μ₁=0, μ₂=0.3):   slope flips from -0.045 to +0.045
    Effect → reversed (μ₁=0.3, μ₂=-0.3): slope flips from +0.045 to -0.135

The sign flip is sharpest in the reversed case (slope magnitude triples). The effect→null case is the hardest to detect because the post-switch slope is small in magnitude.

**What this doesn't prove.** How quickly a changepoint detector notices the flip through noise. A rolling window of size W averages W observations, so the detected slope lags the true switch by up to W/2 steps in the best case. The detection delay is a function of W, the noise level, and the detector's sensitivity. That's the experiment.

## Lemma 3: Stationary effect has no spectral peak

**Setup.** μ(t) = μ (constant). Same e-value.

**Claim.** The periodogram of Δlog(E_t) has no peak.

**Proof.** Δlog(E_t) = log(e_t) = λX_t - λ²/2. Since X_t ~ N(μ, 1) i.i.d., the sequence Δlog(E_t) is i.i.d. with mean λμ - λ²/2 and variance λ². The periodogram of an i.i.d. sequence is flat (white noise). No peak.

**Corollary.** Any spectral peak in the e-value trajectory of a stationary DGP is a finite-sample artifact. The expected periodogram is flat. This is the null baseline for Claim 1: condition (a)'s periodogram should be indistinguishable from white noise.

## Lemma 4: AR(1) doesn't produce periodicity in the evidence

**Setup.** X_t = φX_{t-1} + c + ε_t, |φ| < 1. Same e-value.

**Claim.** The periodogram of Δlog(E_t) shows autocorrelation structure but no discrete spectral peak.

**Proof.** Δlog(E_t) = λX_t - λ²/2. The autocorrelation of X_t under AR(1) is ρ(k) = φᵏ, which decays geometrically. The spectral density of an AR(1) process is continuous:

    S(f) = σ² / |1 - φ·exp(-2πif)|²

This is a smooth function of frequency with a broad peak near f=0 (for positive φ) but no discrete spike. The periodogram of Δlog(E_t) inherits this shape (scaled by λ²). It has a bump at low frequencies but no sharp peak at any specific period.

**What this means.** An AR(1) produces a "warm-colored" periodogram, not a peaked one. A Lotka-Volterra cycle produces a sharp peak. The two are distinguishable spectrally. This is why condition (d) is a useful confound: it tests whether the classifier mistakes broad autocorrelation for a cycle.

## What the experiments add

The lemmas prove the expected behavior. The experiments test:
1. Is the Lotka-Volterra peak detectable above finite-sample noise? (Claim 1)
2. Is spectral detection faster than threshold detection? (Claim 2 — no analytical answer; depends on SNR)
3. How quickly can a detector find the slope flip? (Claim 3 — depends on window size, noise, detector parameters)

If the experiments contradict the lemmas, the implementation has a bug. If the lemmas hold but the experiments show the signal is undetectable, the diagnostic is theoretically correct but practically useless. Both outcomes are informative.
