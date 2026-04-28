"""Generate all figures for the experiment."""

import numpy as np
import yaml
import csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

ROOT = Path(__file__).parent.parent
CONFIG = ROOT / "configs" / "conditions.yaml"
DATA = ROOT / "data"
RESULTS = ROOT / "results"
PLOTS = RESULTS / "plots"
TABLES = RESULTS / "tables"


def load_config():
    with open(CONFIG) as f:
        return yaml.safe_load(f)


def load_data(cond_key, rep=0):
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


def load_csv(path):
    with open(path) as f:
        return list(csv.DictReader(f))


def plot_trajectories():
    """Five e-value trajectories on one log-scale plot."""
    conditions = [
        ("a_stationary", "(a) Stationary"),
        ("i_sinusoidal", "(i) Sinusoidal"),
        ("b_lotka_volterra", "(b) Lotka-Volterra"),
        ("b2_stochastic_lv", "(b2) Stochastic LV"),
        ("c_null", "(c) Null"),
        ("d_ar1", "(d) AR(1)"),
        ("e_regime_switching", "(e) Regime switch"),
    ]

    fig, ax = plt.subplots(figsize=(12, 6))
    for cond_key, label in conditions:
        _, _, cum_log_E = load_data(cond_key, rep=0)
        ax.plot(cum_log_E, label=label, linewidth=0.8)

    ax.set_xlabel("Observation")
    ax.set_ylabel("log(E_t)")
    ax.set_title("E-value trajectories across conditions (rep 0)")
    ax.legend(fontsize=8)
    ax.axhline(y=np.log(20), color="gray", linestyle="--", linewidth=0.5, label="E_t = 20")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(PLOTS / "trajectories.png", dpi=150)
    plt.close(fig)
    print("  trajectories.png")


def plot_periodograms():
    """Periodograms side by side for all conditions."""
    conditions = [
        ("a_stationary", "(a) Stationary"),
        ("i_sinusoidal", "(i) Sinusoidal"),
        ("b_lotka_volterra", "(b) LV"),
        ("b2_stochastic_lv", "(b2) Stoch LV"),
        ("c_null", "(c) Null"),
        ("d_ar1", "(d) AR(1)"),
        ("e_regime_switching", "(e) Regime"),
    ]

    fig, axes = plt.subplots(2, 4, figsize=(16, 8))
    axes = axes.flatten()

    for i, (cond_key, label) in enumerate(conditions):
        _, log_e, _ = load_data(cond_key, rep=0)
        freqs, power = periodogram(log_e)
        periods = 1.0 / freqs

        ax = axes[i]
        ax.loglog(periods, power, linewidth=0.5)
        ax.set_title(label, fontsize=10)
        ax.set_xlabel("Period")
        ax.set_ylabel("Power")
        ax.set_xlim(5, 5000)
        ax.grid(True, alpha=0.3)

    axes[-1].set_visible(False)
    fig.suptitle("Periodograms of Δlog(E_t) (rep 0)", fontsize=12)
    fig.tight_layout()
    fig.savefig(PLOTS / "periodograms.png", dpi=150)
    plt.close(fig)
    print("  periodograms.png")


def plot_classifier_stopping_times():
    """Box plots of stopping times by condition and classifier."""
    csv_path = TABLES / "classifier_race.csv"
    if not csv_path.exists():
        print("  SKIP classifier plot (no data)")
        return

    rows = load_csv(csv_path)

    conditions = [
        "c_null", "a_stationary", "d_ar1", "i_sinusoidal",
        "b_lotka_volterra", "b2_stochastic_lv", "e_regime_switching",
    ]
    classifiers = ["threshold", "bayesian", "spectral_ev"]
    clf_labels = {"threshold": "Threshold", "bayesian": "Bayesian", "spectral_ev": "Spectral EV"}

    fig, axes = plt.subplots(1, len(conditions), figsize=(20, 5), sharey=True)

    for i, cond in enumerate(conditions):
        data_by_clf = {}
        for clf in classifiers:
            stops = [
                int(r["stopping_time"])
                for r in rows
                if r["condition"] == cond and r["classifier"] == clf
            ]
            data_by_clf[clf] = stops

        ax = axes[i]
        bp = ax.boxplot(
            [data_by_clf[c] for c in classifiers],
            labels=[clf_labels[c] for c in classifiers],
            widths=0.6,
        )
        ax.set_title(cond.replace("_", "\n"), fontsize=8)
        ax.tick_params(axis="x", rotation=45, labelsize=7)
        if i == 0:
            ax.set_ylabel("Stopping time")

    fig.suptitle("Classifier stopping times (100 reps)", fontsize=12)
    fig.tight_layout()
    fig.savefig(PLOTS / "classifier_stopping_times.png", dpi=150)
    plt.close(fig)
    print("  classifier_stopping_times.png")


def plot_classifier_labels():
    """Stacked bar chart of label distributions per condition × classifier."""
    csv_path = TABLES / "classifier_race.csv"
    if not csv_path.exists():
        return

    rows = load_csv(csv_path)

    conditions = [
        "c_null", "a_stationary", "d_ar1", "i_sinusoidal",
        "b_lotka_volterra", "b2_stochastic_lv", "e_regime_switching",
    ]
    classifiers = ["threshold", "bayesian", "spectral_ev"]
    all_labels = ["reject_null", "periodic", "null", "undecided"]
    colors = {"reject_null": "#2ca02c", "periodic": "#1f77b4", "null": "#d62728", "undecided": "#999999"}

    fig, axes = plt.subplots(1, len(classifiers), figsize=(15, 5), sharey=True)

    for j, clf in enumerate(classifiers):
        ax = axes[j]
        x_pos = np.arange(len(conditions))
        bottoms = np.zeros(len(conditions))

        for lab in all_labels:
            counts = []
            for cond in conditions:
                c = sum(1 for r in rows if r["condition"] == cond and r["classifier"] == clf and r["label"] == lab)
                counts.append(c)
            counts = np.array(counts)
            if counts.sum() > 0:
                ax.bar(x_pos, counts, bottom=bottoms, label=lab, color=colors[lab], width=0.7)
                bottoms += counts

        ax.set_title(clf)
        ax.set_xticks(x_pos)
        ax.set_xticklabels([c.split("_")[0] for c in conditions], rotation=45, fontsize=7)
        if j == 0:
            ax.set_ylabel("Count (out of 100)")
        ax.legend(fontsize=7)

    fig.suptitle("Label distributions by classifier", fontsize=12)
    fig.tight_layout()
    fig.savefig(PLOTS / "classifier_labels.png", dpi=150)
    plt.close(fig)
    print("  classifier_labels.png")


def plot_changepoint_detection():
    """Changepoint detection delay box plots + example trajectories."""
    csv_path = TABLES / "changepoint_detection.csv"
    if not csv_path.exists():
        print("  SKIP changepoint plot (no data)")
        return

    rows = load_csv(csv_path)

    regime_conds = ["f_effect_to_null", "g_null_to_effect", "h_effect_reversed"]
    methods = ["cusum", "bayesian_cpd"]
    signals = ["evalue", "raw"]

    fig, axes = plt.subplots(1, 3, figsize=(15, 5), sharey=True)

    for i, cond in enumerate(regime_conds):
        ax = axes[i]
        data = []
        labels = []
        for method in methods:
            for signal in signals:
                delays = [
                    float(r["delay"])
                    for r in rows
                    if r["condition"] == cond and r["method"] == method and r["signal"] == signal
                    and r["delay"] != "" and r["delay"] != "nan"
                ]
                delays = [d for d in delays if not np.isnan(d)]
                data.append(delays if delays else [0])
                labels.append(f"{method[:5]}\n{signal[:3]}")

        ax.boxplot(data, labels=labels)
        ax.set_title(cond.replace("_", " "), fontsize=9)
        ax.axhline(y=500, color="red", linestyle="--", linewidth=0.5, label="500-step budget")
        if i == 0:
            ax.set_ylabel("Detection delay |t_det - 5000|")
        ax.tick_params(axis="x", labelsize=7)

    fig.suptitle("Changepoint detection delay (100 reps)", fontsize=12)
    fig.tight_layout()
    fig.savefig(PLOTS / "changepoint_delays.png", dpi=150)
    plt.close(fig)
    print("  changepoint_delays.png")

    fig, axes = plt.subplots(3, 2, figsize=(14, 10))
    for i, cond in enumerate(regime_conds):
        x, log_e, cum_log_E = load_data(cond, rep=0)

        axes[i, 0].plot(cum_log_E, linewidth=0.5)
        axes[i, 0].axvline(x=5000, color="red", linestyle="--", linewidth=0.8)
        axes[i, 0].set_title(f"{cond} — log(E_t)")
        axes[i, 0].set_ylabel("log(E_t)")

        axes[i, 1].plot(x, linewidth=0.3, alpha=0.5)
        window = 200
        rolling_mean = np.convolve(x, np.ones(window) / window, mode="valid")
        axes[i, 1].plot(np.arange(window - 1, len(x)), rolling_mean, linewidth=1, color="red")
        axes[i, 1].axvline(x=5000, color="red", linestyle="--", linewidth=0.8)
        axes[i, 1].set_title(f"{cond} — raw X_t + rolling mean")
        axes[i, 1].set_ylabel("X_t")

    for ax in axes[-1]:
        ax.set_xlabel("Observation")
    fig.suptitle("Regime-switch trajectories (rep 0)", fontsize=12)
    fig.tight_layout()
    fig.savefig(PLOTS / "changepoint_trajectories.png", dpi=150)
    plt.close(fig)
    print("  changepoint_trajectories.png")


def plot_sensitivity():
    """Sensitivity analysis: peak/median ratio across λ values."""
    csv_path = TABLES / "sensitivity.csv"
    if not csv_path.exists():
        print("  SKIP sensitivity plot (no data)")
        return

    rows = load_csv(csv_path)

    conditions = ["i_sinusoidal", "b_lotka_volterra", "b2_stochastic_lv"]
    lambdas = sorted(set(float(r["lambda"]) for r in rows))

    fig, axes = plt.subplots(1, len(conditions), figsize=(14, 5), sharey=True)

    for i, cond in enumerate(conditions):
        ax = axes[i]
        ratios = []
        ratio_errs_lo = []
        ratio_errs_hi = []
        for lam in lambdas:
            r = [row for row in rows if row["condition"] == cond and float(row["lambda"]) == lam]
            if r:
                med = float(r[0]["median_peak_ratio"])
                lo = float(r[0]["iqr_ratio_lo"])
                hi = float(r[0]["iqr_ratio_hi"])
                ratios.append(med)
                ratio_errs_lo.append(med - lo)
                ratio_errs_hi.append(hi - med)

        ax.bar(
            range(len(lambdas)),
            ratios,
            yerr=[ratio_errs_lo, ratio_errs_hi],
            capsize=4,
        )
        ax.set_xticks(range(len(lambdas)))
        ax.set_xticklabels([f"λ={l:.2f}" for l in lambdas])
        ax.set_title(cond.replace("_", "\n"), fontsize=9)
        if i == 0:
            ax.set_ylabel("Peak / median power ratio")
        ax.axhline(y=3, color="red", linestyle="--", linewidth=0.5, label="Detection threshold")

    fig.suptitle("Sensitivity: spectral peak strength vs λ", fontsize=12)
    fig.tight_layout()
    fig.savefig(PLOTS / "sensitivity.png", dpi=150)
    plt.close(fig)
    print("  sensitivity.png")


def main():
    PLOTS.mkdir(parents=True, exist_ok=True)

    print("Generating plots...")
    plot_trajectories()
    plot_periodograms()
    plot_classifier_stopping_times()
    plot_classifier_labels()
    plot_changepoint_detection()
    plot_sensitivity()
    print("Done.")


if __name__ == "__main__":
    main()
