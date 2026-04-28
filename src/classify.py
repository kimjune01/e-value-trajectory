"""Three-classifier race: threshold, Bayesian sequential, spectral e-value.

All three output the same labels: reject_null, periodic, null.
Records stopping time (first observation at which a label is committed).
Runs across 100 replications per condition. Reports median + IQR.
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


def periodogram(x):
    n = len(x)
    x_centered = x - x.mean()
    fft_vals = np.fft.rfft(x_centered)
    power = (np.abs(fft_vals) ** 2) / n
    freqs = np.fft.rfftfreq(n)
    return freqs[1:], power[1:]


def threshold_classifier(cum_log_E, window=500):
    """Reject null when E_t > 20. Null when slope negative over window. Never periodic."""
    log_20 = np.log(20)
    n = len(cum_log_E)

    for t in range(n):
        if cum_log_E[t] > log_20:
            return t, "reject_null"
        if t >= window:
            slope = (cum_log_E[t] - cum_log_E[t - window]) / window
            if slope < 0:
                return t, "null"

    return n, "undecided"


def bayesian_classifier(x, window=500, peak_threshold=3.0):
    """Online Bayesian posterior on μ. Periodic via periodogram on posterior mean."""
    n = len(x)
    prior_var = 100.0
    prior_mean = 0.0

    posterior_means = np.zeros(n)
    sum_x = 0.0

    for t in range(n):
        sum_x += x[t]
        count = t + 1
        post_prec = 1.0 / prior_var + count
        post_mean = (prior_mean / prior_var + sum_x) / post_prec
        post_std = 1.0 / np.sqrt(post_prec)
        posterior_means[t] = post_mean

        p_positive = 1.0 - stats.norm.cdf(0, post_mean, post_std)

        if p_positive > 0.95:
            if t >= window:
                freqs, power = periodogram(posterior_means[t - window + 1 : t + 1])
                if np.median(power) > 0:
                    peak_ratio = np.max(power) / np.median(power)
                    if peak_ratio > peak_threshold:
                        return t, "periodic"
            return t, "reject_null"

        if p_positive < 0.05:
            return t, "null"

        if t >= window:
            freqs, power = periodogram(posterior_means[t - window + 1 : t + 1])
            if np.median(power) > 0:
                peak_ratio = np.max(power) / np.median(power)
                if peak_ratio > peak_threshold:
                    return t, "periodic"

    return n, "undecided"


def spectral_evalue_classifier(log_e, cum_log_E, window=500, peak_threshold=3.0):
    """Sliding-window periodogram on Δlog(E_t). Periodic > reject_null > null."""
    n = len(log_e)

    for t in range(window - 1, n):
        segment = log_e[t - window + 1 : t + 1]
        freqs, power = periodogram(segment)
        median_power = np.median(power)
        peak_ratio = np.max(power) / median_power if median_power > 0 else 0

        slope = (cum_log_E[t] - cum_log_E[t - window + 1]) / (window - 1)

        if peak_ratio > peak_threshold:
            return t, "periodic"
        elif slope > 0:
            return t, "reject_null"
        elif slope < 0:
            return t, "null"

    return n, "undecided"


def main():
    cfg = load_config()
    n_reps = cfg["experiment"]["n_replications"]
    window = cfg["detection"]["spectral"]["window"]
    peak_threshold_orig = cfg["detection"]["spectral"]["peak_threshold"]
    peak_threshold = cfg.get("calibrated", {}).get("spectral", {}).get("peak_threshold", peak_threshold_orig)
    print(f"  peak_threshold: {peak_threshold} (pre-registered: {peak_threshold_orig})")

    conditions = [
        "c_null", "a_stationary", "d_ar1", "i_sinusoidal",
        "b_lotka_volterra", "b2_stochastic_lv", "e_regime_switching",
    ]

    classifiers = {
        "threshold": lambda x, log_e, cum: threshold_classifier(cum, window),
        "bayesian": lambda x, log_e, cum: bayesian_classifier(x, window, peak_threshold),
        "spectral_ev": lambda x, log_e, cum: spectral_evalue_classifier(log_e, cum, window, peak_threshold),
    }

    RESULTS.mkdir(parents=True, exist_ok=True)
    rows = []

    for cond in conditions:
        print(f"  {cond}...")
        for rep in range(n_reps):
            x, log_e, cum_log_E = load_data(cond, rep)
            for clf_name, clf_fn in classifiers.items():
                t_stop, label = clf_fn(x, log_e, cum_log_E)
                rows.append({
                    "condition": cond,
                    "rep": rep,
                    "classifier": clf_name,
                    "stopping_time": t_stop,
                    "label": label,
                })

    out_path = RESULTS / "classifier_race.csv"
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["condition", "rep", "classifier", "stopping_time", "label"])
        writer.writeheader()
        writer.writerows(rows)

    print(f"\nResults written to {out_path}")

    print(f"\n{'Condition':25s} | {'Classifier':12s} | {'Med stop':>8s} {'IQR':>12s} | Labels")
    print("-" * 90)
    for cond in conditions:
        for clf_name in classifiers:
            cond_rows = [r for r in rows if r["condition"] == cond and r["classifier"] == clf_name]
            stops = [r["stopping_time"] for r in cond_rows]
            labels = [r["label"] for r in cond_rows]
            label_counts = {}
            for lb in labels:
                label_counts[lb] = label_counts.get(lb, 0) + 1
            label_str = ", ".join(f"{k}:{v}" for k, v in sorted(label_counts.items()))

            med = np.median(stops)
            q25, q75 = np.percentile(stops, [25, 75])
            print(f"{cond:25s} | {clf_name:12s} | {med:8.0f} {q25:5.0f}-{q75:5.0f} | {label_str}")


if __name__ == "__main__":
    main()
