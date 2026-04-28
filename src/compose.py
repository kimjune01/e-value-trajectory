"""Claim 4: Heterogeneous e-value composition.

Five streams with different distributions sharing a tidal cycle (period T=500).
Each individually too weak to detect; composed e-values reveal the period.

E-values use fixed alternatives (not time-indexed). This keeps log(e_t) linear
in X_t, so the periodogram shows the fundamental frequency — same affine-transform
property as V1. Time-indexed alternatives introduce KL-divergence harmonics at 2/T.

See PREREGISTRATION_V2.md for the full pre-registration.
"""

import numpy as np
import yaml
import csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.special import expit

ROOT = Path(__file__).parent.parent
CONFIG = ROOT / "configs" / "conditions.yaml"
RESULTS = ROOT / "results"
PLOTS = RESULTS / "plots"
TABLES = RESULTS / "tables"


def load_config():
    with open(CONFIG) as f:
        return yaml.safe_load(f)


def sinusoid(n, period, phase, amplitude):
    t = np.arange(n)
    return amplitude * np.sin(2 * np.pi * t / period + phase)


def generate_and_score(name, cfg, n, period, phase, rng, null=False):
    """Generate observations under H0 or H1, compute log e-values against fixed alternative."""
    dist = cfg["distribution"]
    null_p = cfg["null_params"]
    alt_p = cfg["alt_params"]
    A = cfg["amplitude"]
    s = sinusoid(n, period, phase, A)

    if dist == "normal":
        mu_t = 0.0 if null else s
        x = rng.normal(mu_t, null_p["sigma"], n)
        lam = alt_p["lambda"]
        log_e = lam * x - lam ** 2 / 2

    elif dist == "poisson":
        mu0 = null_p["mu"]
        mu_alt = alt_p["mu"]
        if null:
            x = rng.poisson(mu0, n).astype(float)
        else:
            mu_t = mu0 * np.exp(s)
            x = np.array([rng.poisson(mu_t[i]) for i in range(n)], dtype=float)
        log_e = x * np.log(mu_alt / mu0) - (mu_alt - mu0)

    elif dist == "exponential":
        r0 = null_p["rate"]
        r_alt = alt_p["rate"]
        if null:
            x = rng.exponential(1.0 / r0, n)
        else:
            r_t = r0 * np.exp(s)
            x = np.array([rng.exponential(1.0 / r_t[i]) for i in range(n)])
        log_e = np.log(r_alt / r0) - (r_alt - r0) * x

    elif dist == "bernoulli":
        p0 = null_p["p"]
        p_alt = alt_p["p"]
        if null:
            x = rng.binomial(1, p0, n).astype(float)
        else:
            p_t = expit(np.log(p0 / (1 - p0)) + s)
            x = np.array([rng.binomial(1, p_t[i]) for i in range(n)], dtype=float)
        log_e = x * np.log(p_alt / p0) + (1 - x) * np.log((1 - p_alt) / (1 - p0))

    elif dist == "lognormal":
        mu0 = null_p["mu"]
        sigma = null_p["sigma"]
        delta = alt_p["delta"]
        if null:
            x = rng.lognormal(mu0, sigma, n)
        else:
            mu_t = mu0 + s
            x = np.array([rng.lognormal(mu_t[i], sigma) for i in range(n)])
        log_e = delta * (np.log(x) - mu0) / sigma ** 2 - delta ** 2 / (2 * sigma ** 2)

    else:
        raise ValueError(f"Unknown distribution: {dist}")

    return x, log_e


def periodogram(x):
    n = len(x)
    x_centered = x - x.mean()
    fft_vals = np.fft.rfft(x_centered)
    power = (np.abs(fft_vals) ** 2) / n
    freqs = np.fft.rfftfreq(n)
    return freqs[1:], power[1:]


def peak_median_ratio(power):
    med = np.median(power)
    if med == 0:
        return float("inf")
    return float(np.max(power) / med)


def peak_period(freqs, power):
    return float(1.0 / freqs[np.argmax(power)])


def run_experiment(cfg):
    comp = cfg["composition"]
    n = comp["n_observations"]
    n_reps = comp["n_replications"]
    period = comp["period"]
    phase = comp["phase"]
    base_seed = comp["base_seed"]
    streams = comp["streams"]
    stream_names = list(streams.keys())
    K = len(stream_names)

    print(f"Running composition experiment: K={K} streams, N={n}, {n_reps} reps")

    rows = []

    for condition in ["signal", "null"]:
        is_null = condition == "null"
        for rep in range(n_reps):
            seed = base_seed + 99999 + rep
            rng = np.random.default_rng(seed)

            all_log_e = []
            all_x = []

            for sname in stream_names:
                x, log_e = generate_and_score(
                    sname, streams[sname], n, period, phase, rng, null=is_null
                )
                all_log_e.append(log_e)
                all_x.append(x)

                freqs, power = periodogram(log_e)
                pmr = peak_median_ratio(power)
                pp = peak_period(freqs, power)
                e_mean = float(np.exp(log_e).mean())

                rows.append({
                    "condition": condition,
                    "rep": rep,
                    "signal": f"individual_{sname}",
                    "peak_median_ratio": f"{pmr:.2f}",
                    "peak_period": f"{pp:.1f}",
                    "evalue_null_mean": f"{e_mean:.4f}",
                })

            composed_log_e = np.sum(all_log_e, axis=0)
            freqs, power = periodogram(composed_log_e)
            pmr = peak_median_ratio(power)
            pp = peak_period(freqs, power)
            e_mean = float(np.exp(composed_log_e).mean())

            rows.append({
                "condition": condition,
                "rep": rep,
                "signal": "composed_evalue",
                "peak_median_ratio": f"{pmr:.2f}",
                "peak_period": f"{pp:.1f}",
                "evalue_null_mean": f"{e_mean:.4f}",
            })

            z_scores = [(xi - xi.mean()) / xi.std() for xi in all_x]
            z_sum = np.sum(z_scores, axis=0)
            freqs, power = periodogram(z_sum)
            pmr = peak_median_ratio(power)
            pp = peak_period(freqs, power)

            rows.append({
                "condition": condition,
                "rep": rep,
                "signal": "standardized_sum",
                "peak_median_ratio": f"{pmr:.2f}",
                "peak_period": f"{pp:.1f}",
                "evalue_null_mean": "",
            })

    return rows, stream_names


def summarize(rows, stream_names):
    print("\n" + "=" * 70)
    print("RESULTS")
    print("=" * 70)

    signals = [f"individual_{s}" for s in stream_names] + ["composed_evalue", "standardized_sum"]

    for condition in ["signal", "null"]:
        print(f"\n--- {condition.upper()} ---")
        print(f"{'Signal':<30} {'Median PMR':>12} {'IQR':>20} {'Med Period':>12}")
        print("-" * 76)

        for sig in signals:
            pmrs = [float(r["peak_median_ratio"]) for r in rows
                    if r["condition"] == condition and r["signal"] == sig]
            periods = [float(r["peak_period"]) for r in rows
                       if r["condition"] == condition and r["signal"] == sig]
            if not pmrs:
                continue
            med = np.median(pmrs)
            q25, q75 = np.percentile(pmrs, [25, 75])
            med_period = np.median(periods)
            print(f"{sig:<30} {med:>12.1f} {q25:>8.1f}-{q75:.1f} {med_period:>12.1f}")

    null_composed = [float(r["peak_median_ratio"]) for r in rows
                     if r["condition"] == "null" and r["signal"] == "composed_evalue"]
    p95 = np.percentile(null_composed, 95)
    p99 = np.percentile(null_composed, 99)
    print(f"\nComposed null threshold: p95={p95:.1f}, p99={p99:.1f}")

    sig_composed = [float(r["peak_median_ratio"]) for r in rows
                    if r["condition"] == "signal" and r["signal"] == "composed_evalue"]
    detect_rate = np.mean([p > p99 for p in sig_composed])
    print(f"Composed detection rate (vs p99 null): {detect_rate:.0%}")

    for sname in stream_names:
        individual = [float(r["peak_median_ratio"]) for r in rows
                      if r["condition"] == "signal" and r["signal"] == f"individual_{sname}"]
        rate = np.mean([p > 15.0 for p in individual])
        print(f"  {sname} individual detection rate (vs 15): {rate:.0%}")

    std_sum_sig = [float(r["peak_median_ratio"]) for r in rows
                   if r["condition"] == "signal" and r["signal"] == "standardized_sum"]
    null_std = [float(r["peak_median_ratio"]) for r in rows
                if r["condition"] == "null" and r["signal"] == "standardized_sum"]
    std_p99 = np.percentile(null_std, 99)
    std_detect = np.mean([p > std_p99 for p in std_sum_sig])
    print(f"Standardized sum detection rate (vs p99 null): {std_detect:.0%}")
    print(f"Standardized sum null threshold: p95={np.percentile(null_std, 95):.1f}, p99={std_p99:.1f}")

    print("\nPeriod accuracy (composed, signal condition):")
    sig_periods = [float(r["peak_period"]) for r in rows
                   if r["condition"] == "signal" and r["signal"] == "composed_evalue"]
    correct = np.mean([abs(p - 500) < 5 for p in sig_periods])
    print(f"  Period within +/-5 of 500: {correct:.0%}")

    print("\nE-value null means (null condition):")
    for sname in stream_names:
        means = [float(r["evalue_null_mean"]) for r in rows
                 if r["condition"] == "null" and r["signal"] == f"individual_{sname}"
                 and r["evalue_null_mean"]]
        if means:
            print(f"  {sname}: {np.mean(means):.4f}")


def plot_periodograms(cfg, stream_names):
    comp = cfg["composition"]
    n = comp["n_observations"]
    period = comp["period"]
    phase = comp["phase"]
    rng = np.random.default_rng(comp["base_seed"] + 99999)

    all_log_e = []
    fig, axes = plt.subplots(2, 3, figsize=(16, 10))
    axes_flat = axes.flatten()

    for i, sname in enumerate(stream_names):
        x, log_e = generate_and_score(
            sname, comp["streams"][sname], n, period, phase, rng, null=False
        )
        all_log_e.append(log_e)
        freqs, power = periodogram(log_e)
        periods = 1.0 / freqs
        pmr = peak_median_ratio(power)

        ax = axes_flat[i]
        ax.loglog(periods, power, linewidth=0.5, color="C0")
        ax.axvline(x=500, color="red", linestyle="--", linewidth=0.5, alpha=0.5)
        ax.set_title(f"{comp['streams'][sname]['label']}\npeak/med = {pmr:.1f}", fontsize=10)
        ax.set_xlabel("Period")
        ax.set_ylabel("Power")
        ax.set_xlim(5, 5000)
        ax.grid(True, alpha=0.3)

    composed = np.sum(all_log_e, axis=0)
    freqs, power = periodogram(composed)
    periods = 1.0 / freqs
    pmr = peak_median_ratio(power)

    ax = axes_flat[5]
    ax.loglog(periods, power, linewidth=0.8, color="C3")
    ax.axvline(x=500, color="red", linestyle="--", linewidth=0.8)
    ax.set_title(f"Composed (K=5)\npeak/med = {pmr:.1f}", fontsize=10, fontweight="bold")
    ax.set_xlabel("Period")
    ax.set_ylabel("Power")
    ax.set_xlim(5, 5000)
    ax.grid(True, alpha=0.3)

    fig.suptitle("Individual vs composed periodograms (rep 0, illustrative)", fontsize=12)
    fig.tight_layout()
    fig.savefig(PLOTS / "composition_periodograms.png", dpi=150)
    plt.close(fig)
    print("  composition_periodograms.png")


def plot_trajectories(cfg, stream_names):
    comp = cfg["composition"]
    n = comp["n_observations"]
    period = comp["period"]
    phase = comp["phase"]
    rng = np.random.default_rng(comp["base_seed"] + 99999)

    fig, axes = plt.subplots(len(stream_names) + 1, 1, figsize=(14, 12), sharex=True)

    all_log_e = []
    for i, sname in enumerate(stream_names):
        x, log_e = generate_and_score(
            sname, comp["streams"][sname], n, period, phase, rng, null=False
        )
        all_log_e.append(log_e)
        cum = np.cumsum(log_e)

        ax = axes[i]
        ax.plot(cum, linewidth=0.5, color="C0", alpha=0.7)
        ax.set_ylabel(comp["streams"][sname]["label"], fontsize=8)
        ax.grid(True, alpha=0.3)

    composed = np.sum(all_log_e, axis=0)
    cum_composed = np.cumsum(composed)

    ax = axes[-1]
    ax.plot(cum_composed, linewidth=0.8, color="C3")
    ax.set_ylabel("Composed", fontsize=8, fontweight="bold")
    ax.set_xlabel("Observation")
    ax.grid(True, alpha=0.3)

    fig.suptitle("Cumulative log(E_t) trajectories (rep 0)", fontsize=12)
    fig.tight_layout()
    fig.savefig(PLOTS / "composition_trajectories.png", dpi=150)
    plt.close(fig)
    print("  composition_trajectories.png")


def main():
    cfg = load_config()
    PLOTS.mkdir(parents=True, exist_ok=True)
    TABLES.mkdir(parents=True, exist_ok=True)

    rows, stream_names = run_experiment(cfg)

    csv_path = TABLES / "composition.csv"
    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "condition", "rep", "signal",
            "peak_median_ratio", "peak_period", "evalue_null_mean",
        ])
        writer.writeheader()
        writer.writerows(rows)
    print(f"  {csv_path}")

    summarize(rows, stream_names)

    print("\nGenerating plots...")
    plot_periodograms(cfg, stream_names)
    plot_trajectories(cfg, stream_names)
    print("\nDone.")


if __name__ == "__main__":
    main()
