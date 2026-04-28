"""Periodogram analysis on Δlog(E_t) and raw X_t for all conditions."""

import numpy as np
from pathlib import Path

ROOT = Path(__file__).parent.parent
DATA = ROOT / "data"
RESULTS = ROOT / "results"


def periodogram(x):
    """Compute one-sided periodogram. Returns (frequencies, power)."""
    n = len(x)
    x_centered = x - x.mean()
    fft = np.fft.rfft(x_centered)
    power = (np.abs(fft) ** 2) / n
    freqs = np.fft.rfftfreq(n)
    return freqs[1:], power[1:]  # drop DC


def peak_stats(freqs, power):
    """Return peak frequency, peak period, peak power, and peak/median ratio."""
    idx = np.argmax(power)
    peak_freq = freqs[idx]
    peak_period = 1.0 / peak_freq if peak_freq > 0 else np.inf
    peak_power = power[idx]
    median_power = np.median(power)
    ratio = peak_power / median_power if median_power > 0 else np.inf
    return peak_freq, peak_period, peak_power, ratio


def analyze_one(cond_key, rep=0):
    """Analyze one replication: periodogram on Δlog(E_t) and raw X_t."""
    x_raw = np.loadtxt(DATA / cond_key / f"{rep:04d}.csv")
    ev = np.loadtxt(DATA / cond_key / f"{rep:04d}_evalues.csv", skiprows=1)
    log_e = ev[:, 0]  # Δlog(E_t) = log(e_t)

    freqs_ev, power_ev = periodogram(log_e)
    freqs_raw, power_raw = periodogram(x_raw)

    ev_stats = peak_stats(freqs_ev, power_ev)
    raw_stats = peak_stats(freqs_raw, power_raw)

    return {
        "evalue": {"freqs": freqs_ev, "power": power_ev, "peak": ev_stats},
        "raw": {"freqs": freqs_raw, "power": power_raw, "peak": raw_stats},
    }


def main():
    conditions = [
        "c_null", "a_stationary", "d_ar1", "i_sinusoidal",
        "b_lotka_volterra", "b2_stochastic_lv", "e_regime_switching",
    ]

    print(f"{'Condition':25s} | {'EV peak period':>14s} {'EV peak/med':>11s} | {'Raw peak period':>14s} {'Raw peak/med':>11s}")
    print("-" * 85)

    for cond in conditions:
        results = analyze_one(cond, rep=0)
        ev_p = results["evalue"]["peak"]
        raw_p = results["raw"]["peak"]
        print(f"{cond:25s} | {ev_p[1]:14.1f} {ev_p[3]:11.1f} | {raw_p[1]:14.1f} {raw_p[3]:11.1f}")

    # Aggregate across replications for cyclic conditions
    print("\n--- Claim 1: peak period distribution across 100 reps ---")
    for cond in ["i_sinusoidal", "b_lotka_volterra", "b2_stochastic_lv", "c_null"]:
        periods = []
        ratios = []
        for rep in range(100):
            r = analyze_one(cond, rep)
            periods.append(r["evalue"]["peak"][1])
            ratios.append(r["evalue"]["peak"][3])
        print(f"{cond:25s}  period: {np.median(periods):7.1f} (IQR {np.percentile(periods,25):.0f}-{np.percentile(periods,75):.0f})  "
              f"peak/med: {np.median(ratios):6.1f} (IQR {np.percentile(ratios,25):.1f}-{np.percentile(ratios,75):.1f})")


if __name__ == "__main__":
    main()
