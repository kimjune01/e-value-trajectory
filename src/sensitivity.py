"""Sensitivity analysis: rerun spectral analysis with λ at 0.15, 0.3, 0.6."""

import numpy as np
import yaml
import csv
from pathlib import Path

ROOT = Path(__file__).parent.parent
CONFIG = ROOT / "configs" / "conditions.yaml"
DATA = ROOT / "data"
RESULTS = ROOT / "results" / "tables"


def load_config():
    with open(CONFIG) as f:
        return yaml.safe_load(f)


def periodogram(x):
    n = len(x)
    x_centered = x - x.mean()
    fft_vals = np.fft.rfft(x_centered)
    power = (np.abs(fft_vals) ** 2) / n
    freqs = np.fft.rfftfreq(n)
    return freqs[1:], power[1:]


def peak_stats(freqs, power):
    idx = np.argmax(power)
    peak_freq = freqs[idx]
    peak_period = 1.0 / peak_freq if peak_freq > 0 else np.inf
    peak_power = power[idx]
    median_power = np.median(power)
    ratio = peak_power / median_power if median_power > 0 else np.inf
    return peak_freq, peak_period, peak_power, ratio


def main():
    cfg = load_config()
    n_reps = cfg["experiment"]["n_replications"]
    lambda_base = cfg["evalue"]["lambda"]
    multipliers = cfg["sensitivity"]["lambda_multipliers"]

    conditions = [
        "i_sinusoidal", "b_lotka_volterra", "b2_stochastic_lv",
        "c_null", "a_stationary", "d_ar1",
    ]

    RESULTS.mkdir(parents=True, exist_ok=True)
    rows = []

    for mult in multipliers:
        lam = lambda_base * mult
        print(f"  λ = {lam:.2f} (×{mult})...")

        for cond in conditions:
            periods = []
            ratios = []
            for rep in range(n_reps):
                x = np.loadtxt(DATA / cond / f"{rep:04d}.csv")
                log_e = lam * x - lam ** 2 / 2
                freqs, power = periodogram(log_e)
                _, period, _, ratio = peak_stats(freqs, power)
                periods.append(period)
                ratios.append(ratio)

            med_period = np.median(periods)
            iqr_period_lo, iqr_period_hi = np.percentile(periods, [25, 75])
            med_ratio = np.median(ratios)
            iqr_ratio_lo, iqr_ratio_hi = np.percentile(ratios, [25, 75])

            rows.append({
                "lambda": lam,
                "multiplier": mult,
                "condition": cond,
                "median_period": med_period,
                "iqr_period_lo": iqr_period_lo,
                "iqr_period_hi": iqr_period_hi,
                "median_peak_ratio": med_ratio,
                "iqr_ratio_lo": iqr_ratio_lo,
                "iqr_ratio_hi": iqr_ratio_hi,
            })

    out_path = RESULTS / "sensitivity.csv"
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "lambda", "multiplier", "condition",
                "median_period", "iqr_period_lo", "iqr_period_hi",
                "median_peak_ratio", "iqr_ratio_lo", "iqr_ratio_hi",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)

    print(f"\nResults written to {out_path}")

    print(f"\n{'λ':>5s} | {'Condition':25s} | {'Med period':>10s} {'IQR':>14s} | {'Med peak/med':>12s} {'IQR':>14s}")
    print("-" * 95)
    for row in rows:
        print(
            f"{row['lambda']:5.2f} | {row['condition']:25s} | "
            f"{row['median_period']:10.1f} {row['iqr_period_lo']:6.0f}-{row['iqr_period_hi']:6.0f} | "
            f"{row['median_peak_ratio']:12.1f} {row['iqr_ratio_lo']:6.1f}-{row['iqr_ratio_hi']:6.1f}"
        )


if __name__ == "__main__":
    main()
