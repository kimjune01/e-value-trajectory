# Sources: Hypothesis loop (sensemaking formalization)

Prior art search for the mathematical formalization of the perturb → classify → decide → perturb loop. The gap: the pieces exist separately but nobody has unified them.

## Active sequential hypothesis testing (the loop skeleton)

- Chernoff, H. (1959). "Sequential Design of Experiments." *Ann. Math. Statist.*, 30(3), 755–770.
  - Foundational. Choose experiments adaptively as evidence accumulates, stop when confident.
  - High for loop structure. Medium for regime classification. No sensemaking framing.

- Nitinawarat, S., Atia, G., & Veeravalli, V. (2013). "Controlled Sensing for Multihypothesis Testing." *IEEE Trans. Inform. Theory*, 59(9), 5773–5789.
  - Adaptive observation/control action selection for multihypothesis testing. Asymptotic optimality.
  - Very high for perturb-observe-classify-decide with convergence. Hypotheses are fixed statistical models, not dynamical regimes.

- Naghshvar, M., & Javidi, T. (2013). "Active Sequential Hypothesis Testing." *Ann. Statist.*, 41(6), 2703–2738.
  - Policies trading information gain, stopping time, and error risk. Asymptotic optimality.
  - Very high as decision-theoretic formalization of the loop.

## E-value composition (the evidence layer)

- Grünwald, P., de Heide, R., & Koolen, W. M. (2024). "Safe Testing." *J. R. Stat. Soc. B*, 86(5), 1091–1128.
  - E-values compose across optionally continued studies, even when the next study depends on previous outcomes.
  - Supplies the evidential accounting. Does not specify what to perturb next.

- Vovk, V., & Wang, R. (2021). "E-values: Calibration, combination, and applications." *Ann. Statist.*, 49(3), 1736–1754.

- Ramdas, A., Grünwald, P., Vovk, V., & Shafer, G. (2023). "Game-theoretic statistics and safe anytime-valid inference." *Statist. Sci.*, 38(4), 576–601.

## Adaptive perturbation for regime classification

- Kevrekidis et al. (2003). "Adaptive Detection of Instabilities." *Physica D*, various.
  - Active online detection of bifurcations by coupling experiments to feedback/control laws.
  - Very high for "perturb a dynamical system until a qualitative regime boundary is found." Narrow.

- Cheong, S., & Krstic, M. (2015). "Input Design for Discrimination Between Classes of LTI Models." *Automatica*.
  - Designs input signals to discriminate among classes (normal/fault modes).
  - High for perturbation design aimed at class discrimination. Finite LTI models.

- Silk, D., et al. (2011). "Designing Attractive Models via Automated Identification of Chaotic and Oscillatory Dynamical Regimes." *Nature Communications*, 2, 489.
  - Qualitative inference for finding parameter regions with attractor properties (oscillation, chaos, hyperchaos).
  - High for regime classification. Not adaptive perturbation. No convergence theorem for the loop.

## Active causal structure learning (the convergence guarantee)

- He, Y.-B., & Geng, Z. (2008). "Active Learning of Causal Networks with Intervention Experiments and Optimal Designs." *JMLR*, 9, 2523–2547.
  - Start from observational equivalence class, choose interventions to orient edges.
  - Very high for adaptive intervention and causal-structure recovery. Assumes DAGs, not dynamical systems.

- Hauser, A., & Bühlmann, P. (2014). "Two Optimal Strategies for Active Learning of Causal Models from Interventional Data." *Int. J. Approx. Reason.*, 55(4), 926–939.
  - Minimal interventions for full DAG identifiability.

- Mania, H., Jordan, M., & Recht, B. (2022). "Active Learning for Nonlinear System Identification with Guarantees." *JMLR*, 23, 1–30.
  - Active learning for nonlinear dynamical systems. Parametric-rate estimation guarantees.
  - High for adaptive perturbation with guarantees. Estimates dynamics, not qualitative classes.

## Sensemaking / OODA (prose, no theorems)

- Klein, G., et al. (2023). "The Plausibility Transition Model for Sensemaking." *Frontiers in Psychology*, 14, 1160132.
  - Sensemaking as state transitions, stopping when plausibility is sufficient.
  - Psychological/narrative model, not a convergence theorem.

- Brehmer, B. (2005). "The Dynamic OODA Loop."
  - Imports cybernetic feedback/control into command-and-control.
  - Functional architecture, not information-theoretic theorem.

- Weick, K. (1995). *Sensemaking in Organizations*. Sage.
  - No mathematical formalization.

## The gap

Nobody has combined:
1. E-value composition across adaptive experiments (Grünwald)
2. Dynamical regime classification (Silk, Kevrekidis, our four-bin classifier)
3. A perturbation-selection policy (Chernoff, Nitinawarat)
4. A convergence guarantee (He & Geng, Hauser & Bühlmann)

into a unified theory where iterating the loop provably recovers the system's dynamical structure.
