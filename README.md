# E-value trajectory diagnostics

Can we read the *shape* of an e-value trajectory to learn about the system, not just test the hypothesis?

E-values are used as tests: did the evidence cross the threshold? This project asks what the trajectory shape tells you *before* it crosses. Oscillation, convergence rate, regime changes in the evidence may diagnose properties of the underlying system that a final verdict misses.

## Thesis

From [Evidence has a trajectory](https://june.kim/evidence-has-a-trajectory):

> The four-bin classification applies to e-values when the experiment tracks the system's dynamics. E-values are granular enough to transmit what the system is doing without compressing it to a single number.

The four bins: convergence (effect absent or stable), divergence (effect real and growing), oscillation (effect is cyclic or context-dependent), chaos (experimental design is measuring noise).

## Two claims to test

1. **E-value trajectories from cyclic systems oscillate.** If the data-generating process has feedback, the e-value trajectory should inherit the periodicity. If it doesn't (smooth convergence/divergence regardless of system dynamics), the thesis is wrong.

2. **Trajectory shape classifies faster than threshold crossing.** A simple classifier on the trajectory window (slope, variance, autocorrelation) should distinguish real effects from nulls earlier than waiting for the e-value to cross a rejection boundary.

## Data sources

Data is gitignored. See [SOURCES.md](SOURCES.md) for download instructions.

## Experiment design

See [EXPERIMENT.md](EXPERIMENT.md).
