"""CUSUM and Bayesian CPD on regime-switch conditions.

Runs on log(E_t) increments and raw X_t for head-to-head comparison.
Measures detection delay |t_detected - 5000| and false alarm rate.
"""

import numpy as np
import yaml
import csv
from pathlib import Path
from scipy import stats

ROOT = Path(__file__).parent.parent
CONFIG = ROOT / "configs" / "conditions.yaml"
DATA = ROOT / "data"
RESULTS = ROOT / "results" / "tables"


def load_config():
    with open(CONFIG) as f:
        return yaml.safe_load(f)


def load_data(cond_key, rep):
    x = np.loadtxt(DATA / cond_key / f"{rep:04d}.csv")
    ev = np.loadtxt(DATA / cond_key / f"{rep:04d}_evalues.csv", skiprows=1)
    return x, ev[:, 0], ev[:, 1]


def cusum_detect(signal, mu0, k, h):
    """Two-sided CUSUM. Returns first alarm time or len(signal) if none."""
    n = len(signal)
    s_pos = 0.0
    s_neg = 0.0
    for t in range(n):
        s_pos = max(0, s_pos + (signal[t] - mu0 - k))
        s_neg = max(0, s_neg + (mu0 - k - signal[t]))
        if s_pos > h or s_neg > h:
            return t
    return n


def bayesian_cpd_detect(signal, hazard, prior_mu=0.0, prior_var=100.0, obs_var=None, threshold=0.3):
    """Adams & MacKay (2007) online Bayesian changepoint detection.
    Returns first time P(run_length < 10) > threshold after a burn-in of 100 steps."""
    if obs_var is None:
        obs_var = float(np.var(signal[:min(500, len(signal))]))

    n = len(signal)
    max_r = min(n, 500)

    R = np.zeros(max_r + 1)
    R[0] = 1.0
    run_sums = np.zeros(max_r + 1)

    for t in range(n):
        active = min(t + 1, max_r)
        rs = np.arange(active)

        post_prec = 1.0 / prior_var + rs / obs_var
        post_mean = np.where(
            rs > 0,
            (prior_mu / prior_var + run_sums[:active] / obs_var) / post_prec,
            prior_mu,
        )
        post_var = 1.0 / post_prec

        pred_var = obs_var + post_var
        pred_prob = np.exp(-0.5 * (signal[t] - post_mean) ** 2 / pred_var) / np.sqrt(
            2 * np.pi * pred_var
        )

        joint = R[:active] * pred_prob

        new_R = np.zeros(min(active + 1, max_r + 1))
        new_R[0] = hazard * joint.sum()
        growth_limit = min(active, max_r)
        new_R[1 : growth_limit + 1] = (1 - hazard) * joint[:growth_limit]

        evidence = new_R.sum()
        if evidence > 0:
            new_R /= evidence

        new_sums = np.zeros(max_r + 1)
        new_sums[1 : growth_limit + 1] = run_sums[:growth_limit] + signal[t]

        R[:] = 0
        R[: len(new_R)] = new_R
        run_sums = new_sums

        if t > 100:
            short_run_prob = R[:10].sum()
            if short_run_prob > threshold:
                return t

    return n


def main():
    cfg = load_config()
    n_reps = cfg["experiment"]["n_replications"]
    lam = cfg["evalue"]["lambda"]
    cusum_k_factor = cfg["detection"]["cusum"]["reference_k_factor"]
    cusum_h_orig = cfg["detection"]["cusum"]["threshold_h"]
    cusum_h = cfg.get("calibrated", {}).get("cusum", {}).get("threshold_h", cusum_h_orig)
    print(f"  cusum_h: {cusum_h} (pre-registered: {cusum_h_orig})")
    hazard = cfg["detection"]["bayesian_cpd"]["hazard_rate"]
    switch_at = 5000

    mean_slope_effect = lam * 0.3 - lam ** 2 / 2  # 0.045
    cusum_k = cusum_k_factor * mean_slope_effect    # 0.0225

    regime_conditions = cfg["regime_switch_conditions"]
    false_alarm_conditions = ["a_stationary", "c_null"]

    RESULTS.mkdir(parents=True, exist_ok=True)
    rows = []

    for cond_key, cond_cfg in regime_conditions.items():
        print(f"  {cond_cfg['label']}...")
        mu_before = cond_cfg["mu_before"]
        mu_after = cond_cfg["mu_after"]

        ev_mu0 = lam * mu_before - lam ** 2 / 2
        raw_mu0 = mu_before

        for rep in range(n_reps):
            x, log_e, cum_log_E = load_data(cond_key, rep)

            t_cusum_ev = cusum_detect(log_e, ev_mu0, cusum_k, cusum_h)
            t_cusum_raw = cusum_detect(x, raw_mu0, cusum_k, cusum_h)
            t_bayes_ev = bayesian_cpd_detect(log_e, hazard, prior_mu=ev_mu0, obs_var=lam ** 2)
            t_bayes_raw = bayesian_cpd_detect(x, hazard, prior_mu=raw_mu0, obs_var=1.0)

            for method, signal_type, t_det in [
                ("cusum", "evalue", t_cusum_ev),
                ("cusum", "raw", t_cusum_raw),
                ("bayesian_cpd", "evalue", t_bayes_ev),
                ("bayesian_cpd", "raw", t_bayes_raw),
            ]:
                delay = abs(t_det - switch_at) if t_det < len(x) else np.nan
                rows.append({
                    "condition": cond_key,
                    "rep": rep,
                    "method": method,
                    "signal": signal_type,
                    "detection_time": t_det,
                    "delay": delay,
                    "false_alarm": 0,
                })

    print("\n  False alarm analysis...")
    for cond_key in false_alarm_conditions:
        for rep in range(n_reps):
            x, log_e, cum_log_E = load_data(cond_key, rep)

            ev_mu0 = lam * (0.3 if cond_key == "a_stationary" else 0.0) - lam ** 2 / 2
            raw_mu0 = 0.3 if cond_key == "a_stationary" else 0.0

            t_cusum_ev = cusum_detect(log_e, ev_mu0, cusum_k, cusum_h)
            t_cusum_raw = cusum_detect(x, raw_mu0, cusum_k, cusum_h)
            t_bayes_ev = bayesian_cpd_detect(log_e, hazard, prior_mu=ev_mu0, obs_var=lam ** 2)
            t_bayes_raw = bayesian_cpd_detect(x, hazard, prior_mu=raw_mu0, obs_var=1.0)

            n = len(x)
            for method, signal_type, t_det in [
                ("cusum", "evalue", t_cusum_ev),
                ("cusum", "raw", t_cusum_raw),
                ("bayesian_cpd", "evalue", t_bayes_ev),
                ("bayesian_cpd", "raw", t_bayes_raw),
            ]:
                rows.append({
                    "condition": cond_key,
                    "rep": rep,
                    "method": method,
                    "signal": signal_type,
                    "detection_time": t_det,
                    "delay": np.nan,
                    "false_alarm": 1 if t_det < n else 0,
                })

    out_path = RESULTS / "changepoint_detection.csv"
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(
            f, fieldnames=["condition", "rep", "method", "signal", "detection_time", "delay", "false_alarm"]
        )
        writer.writeheader()
        writer.writerows(rows)

    print(f"\nResults written to {out_path}")

    print(f"\n{'Condition':20s} | {'Method':14s} | {'Signal':7s} | {'Med delay':>9s} {'IQR':>12s} | {'FA rate':>7s}")
    print("-" * 85)

    for cond_key in list(regime_conditions.keys()) + false_alarm_conditions:
        for method in ["cusum", "bayesian_cpd"]:
            for signal in ["evalue", "raw"]:
                cond_rows = [
                    r for r in rows
                    if r["condition"] == cond_key and r["method"] == method and r["signal"] == signal
                ]
                if not cond_rows:
                    continue

                if cond_key in regime_conditions:
                    delays = [r["delay"] for r in cond_rows if not np.isnan(r["delay"])]
                    if delays:
                        med = np.median(delays)
                        q25, q75 = np.percentile(delays, [25, 75])
                        delay_str = f"{med:9.0f} {q25:5.0f}-{q75:5.0f}"
                    else:
                        delay_str = f"{'N/A':>9s} {'':>12s}"
                    fa_str = ""
                else:
                    delay_str = f"{'':>9s} {'':>12s}"
                    fa_count = sum(r["false_alarm"] for r in cond_rows)
                    fa_str = f"{fa_count}/{len(cond_rows)}"

                print(f"{cond_key:20s} | {method:14s} | {signal:7s} | {delay_str} | {fa_str:>7s}")


if __name__ == "__main__":
    main()
