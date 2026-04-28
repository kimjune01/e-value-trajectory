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
