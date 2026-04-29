"""Claim 5: Four-bin forcing-pattern classification on composed e-value trajectories.

Kill-condition decision tree: monotone → curvature → periodicity → aperiodic → null.
See PREREGISTRATION_V3.md.
"""

import numpy as np
import yaml
import csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.special import expit
from scipy.stats import theilslopes
from scipy.signal import periodogram as sp_periodogram
from itertools import product as iterproduct

ROOT = Path(__file__).parent.parent
CONFIG = ROOT / "configs" / "conditions.yaml"
RESULTS = ROOT / "results"
PLOTS = RESULTS / "plots"
TABLES = RESULTS / "tables"

LABELS = ["null", "divergent", "convergent", "oscillatory", "aperiodic"]
CONDITIONS = ["null", "divergence", "convergence", "oscillation", "aperiodic"]


def load_config():
    with open(CONFIG) as f:
        return yaml.safe_load(f)


# ── Forcing generators ──────────────────────────────────────────────

def make_forcing(condition, n, rep):
    t = np.arange(n)
    if condition == "null":
        return np.zeros(n)
    elif condition == "divergence":
        return t / n - 0.5
    elif condition == "convergence":
        return np.exp(-t / 2000.0)
    elif condition == "oscillation":
        return np.sin(2 * np.pi * t / 500)
    elif condition == "aperiodic":
        z = np.zeros(n + 1000)
        z[0] = 0.1 + 0.008 * (rep % 100)
        for i in range(1, len(z)):
            z[i] = 3.9 * z[i - 1] * (1 - z[i - 1])
        z = z[1000:]  # burn-in
        return (z - z.mean()) / z.std()
    else:
        raise ValueError(condition)


# ── Stream generation (reuse V2 logic) ──────────────────────────────

def generate_stream(cfg, forcing, rng):
    dist = cfg["distribution"]
    null_p = cfg["null_params"]
    alt_p = cfg["alt_params"]
    n = len(forcing)

    if dist == "normal":
        x = rng.normal(forcing, null_p["sigma"], n)
        lam = alt_p["lambda"]
        log_e = lam * x - lam ** 2 / 2

    elif dist == "poisson":
        mu0 = null_p["mu"]
        mu_alt = alt_p["mu"]
        mu_t = mu0 * np.exp(forcing)
        x = np.array([rng.poisson(mu_t[i]) for i in range(n)], dtype=float)
        log_e = x * np.log(mu_alt / mu0) - (mu_alt - mu0)

    elif dist == "exponential":
        r0 = null_p["rate"]
        r_alt = alt_p["rate"]
        r_t = r0 * np.exp(forcing)
        x = np.array([rng.exponential(1.0 / r_t[i]) for i in range(n)])
        log_e = np.log(r_alt / r0) - (r_alt - r0) * x

    elif dist == "bernoulli":
        p0 = null_p["p"]
        p_alt = alt_p["p"]
        p_t = expit(np.log(p0 / (1 - p0)) + forcing)
        x = np.array([rng.binomial(1, p_t[i]) for i in range(n)], dtype=float)
        log_e = x * np.log(p_alt / p0) + (1 - x) * np.log((1 - p_alt) / (1 - p0))

    elif dist == "lognormal":
        mu0 = null_p["mu"]
        sigma = null_p["sigma"]
        delta = alt_p["delta"]
        x = np.array([rng.lognormal(mu0 + forcing[i], sigma) for i in range(n)])
        log_e = delta * (np.log(x) - mu0) / sigma ** 2 - delta ** 2 / (2 * sigma ** 2)

    else:
        raise ValueError(dist)

    return x, log_e


# ── Feature computation ─────────────────────────────────────────────

def preprocess(x):
    lo, hi = np.percentile(x, [0.5, 99.5])
    z = np.clip(x, lo, hi)
    med = np.median(z)
    mad = np.median(np.abs(z - med))
    if mad == 0:
        mad = 1.0
    z = (z - med) / (1.4826 * mad)
    return z


def mann_kendall_z(x):
    n = len(x)
    if n > 1000:
        idx = np.linspace(0, n - 1, 1000, dtype=int)
        x = x[idx]
        n = len(x)
    s = 0
    for i in range(n - 1):
        diff = x[i + 1:] - x[i]
        s += np.sum(np.sign(diff))
    var_s = n * (n - 1) * (2 * n + 5) / 18.0
    if s > 0:
        return (s - 1) / np.sqrt(var_s)
    elif s < 0:
        return (s + 1) / np.sqrt(var_s)
    return 0.0


def curvature_ratio(z):
    n = len(z)
    t = np.arange(n, dtype=float)
    # Linear fit
    coeffs_lin = np.polyfit(t, z, 1)
    rss_lin = np.sum((z - np.polyval(coeffs_lin, t)) ** 2)
    # Exponential fit: z ≈ a * exp(-t/tau) + c
    # Use log-linear fit on |z| for robustness
    try:
        z_shifted = z - z[-500:].mean()
        sign = np.sign(z_shifted[:500].mean()) if abs(z_shifted[:500].mean()) > 0.01 else 1.0
        z_pos = sign * z_shifted
        z_pos_clipped = np.maximum(z_pos, 1e-10)
        valid = z_pos > 0.01
        if valid.sum() > 100:
            log_z = np.log(z_pos_clipped[valid])
            t_valid = t[valid]
            coeffs_exp = np.polyfit(t_valid, log_z, 1)
            z_exp_fit = np.exp(coeffs_exp[1]) * np.exp(coeffs_exp[0] * t)
            z_exp_recon = sign * z_exp_fit + z[-500:].mean()
            rss_exp = np.sum((z - z_exp_recon) ** 2)
        else:
            rss_exp = rss_lin
    except Exception:
        rss_exp = rss_lin
    if rss_lin == 0:
        return 1.0
    return rss_exp / rss_lin


def rolling_rms_slope(z, w=500):
    n = len(z)
    n_windows = n // w
    if n_windows < 3:
        return 0.0
    rms = np.array([np.sqrt(np.mean(z[i * w:(i + 1) * w] ** 2)) for i in range(n_windows)])
    log_rms = np.log(rms + 1e-10)
    t_w = np.arange(n_windows, dtype=float)
    slope, _, _, _ = theilslopes(log_rms, t_w)
    return float(slope)


def energy_ratio(z):
    n = len(z)
    fifth = n // 5
    head = np.sqrt(np.mean(z[:fifth] ** 2))
    tail = np.sqrt(np.mean(z[-fifth:] ** 2))
    if head == 0:
        return 1.0
    return tail / head


def spectral_features(z):
    freqs, power = sp_periodogram(z - z.mean(), scaling="spectrum")
    freqs = freqs[1:]
    power = power[1:]
    if len(power) == 0 or power.sum() == 0:
        return 0.0, 0.0, 0.0
    g_spec = float(np.max(power) / power.sum())
    peak_idx = np.argmax(power)
    peak_freq = freqs[peak_idx]
    # Q factor: peak freq / half-power bandwidth
    half_power = power[peak_idx] / 2
    above = power >= half_power
    bandwidth_bins = max(1, above.sum())
    df = freqs[1] - freqs[0] if len(freqs) > 1 else 1.0
    bandwidth = bandwidth_bins * df
    q_spec = float(peak_freq / bandwidth) if bandwidth > 0 else 0.0
    peak_period = float(1.0 / peak_freq) if peak_freq > 0 else float("inf")
    return g_spec, q_spec, peak_period


def zero_one_test(z, n_c=30):
    n = len(z)
    if n > 2000:
        step = n // 2000
        z = z[::step]
        n = len(z)
    rng_c = np.random.default_rng(0)
    c_vals = rng_c.uniform(np.pi / 5, 4 * np.pi / 5, n_c)
    ks = []
    t = np.arange(1, n + 1)
    for c in c_vals:
        p = np.cumsum(z * np.cos(c * t))
        q = np.cumsum(z * np.sin(c * t))
        n_cut = min(n // 10, 200)
        if n_cut < 10:
            n_cut = n // 2
        lags = np.arange(1, n_cut + 1)
        ms = np.array([np.mean((p[j:] - p[:-j]) ** 2 + (q[j:] - q[:-j]) ** 2)
                        for j in lags])
        if np.std(ms) > 0:
            k = np.corrcoef(lags.astype(float), ms)[0, 1]
        else:
            k = 0.0
        ks.append(k)
    return float(np.median(ks))


def permutation_entropy(z, m=5):
    n = len(z)
    from math import factorial
    max_perms = factorial(m)
    counts = {}
    for i in range(n - m + 1):
        pattern = tuple(np.argsort(z[i:i + m]))
        counts[pattern] = counts.get(pattern, 0) + 1
    total = sum(counts.values())
    h = 0.0
    for c in counts.values():
        p = c / total
        if p > 0:
            h -= p * np.log2(p)
    return h / np.log2(max_perms)


def compute_features(x):
    z = preprocess(x)
    dx = np.diff(z)
    mad_dx = np.median(np.abs(dx - np.median(dx)))
    if mad_dx == 0:
        mad_dx = 1.0

    n_z = len(z)
    if n_z > 1000:
        idx = np.linspace(0, n_z - 1, 1000, dtype=int)
        s_slope, _, _, _ = theilslopes(z[idx], idx.astype(float))
    else:
        s_slope, _, _, _ = theilslopes(z, np.arange(n_z, dtype=float))
    s = s_slope / mad_dx

    z_mk = mann_kendall_z(z)
    r_curve = curvature_ratio(z)
    rho_env = rolling_rms_slope(z)
    e_rat = energy_ratio(z)
    g_spec, q_spec, peak_period = spectral_features(z)
    k_01 = zero_one_test(z)
    pe = permutation_entropy(z)

    med_first = np.median(z[:len(z) // 10])
    med_last = np.median(z[-len(z) // 10:])

    return {
        "S": s,
        "Z_MK": z_mk,
        "R_curve": r_curve,
        "rho_env": rho_env,
        "E_ratio": e_rat,
        "G_spec": g_spec,
        "Q_spec": q_spec,
        "K_01": k_01,
        "PE": pe,
        "med_first": med_first,
        "med_last": med_last,
        "peak_period": peak_period,
        "MAD_dx": mad_dx,
    }


# ── Kill-condition classifier ───────────────────────────────────────

def classify(features, thresholds):
    f = features
    th = thresholds

    # TEST 1: Monotone trend
    if abs(f["S"]) > th["S_q99"] and abs(f["Z_MK"]) > 2.58 \
            and abs(f["med_last"] - f["med_first"]) > 3 * f["MAD_dx"]:
        # TEST 2: Curvature
        if f["R_curve"] < th["R_curve_q05"]:
            return "convergent", "monotone→curvature→convergent"
        else:
            return "divergent", "monotone→curvature→divergent"

    # TEST 3: Periodicity
    n = 10000
    if f["G_spec"] > th["G_spec_q99"] and f["Q_spec"] > 5:
        period = f["peak_period"]
        if 10 < period < n / 8:
            return "oscillatory", "periodicity→oscillatory"

    # TEST 4: Aperiodic structure
    if f["K_01"] > 0.8 and f["PE"] >= 0.55 and f["PE"] <= 0.95:
        return "aperiodic", "aperiodic→K01+PE"

    return "null", "null (no test triggered)"


# ── Null calibration ────────────────────────────────────────────────

def run_null_calibration(cfg, n_reps=1000):
    comp = cfg["composition"]
    streams = comp["streams"]
    stream_names = list(streams.keys())
    n = comp["n_observations"]

    print(f"Null calibration: {n_reps} reps...")
    null_features = []

    for rep in range(n_reps):
        seed = 80000 + rep
        rng = np.random.default_rng(seed)
        forcing = np.zeros(n)

        all_log_e = []
        for sname in stream_names:
            x, log_e = generate_stream(streams[sname], forcing, rng)
            all_log_e.append(log_e)

        composed = np.sum(all_log_e, axis=0)
        feats = compute_features(composed)
        null_features.append(feats)

        if (rep + 1) % 200 == 0:
            print(f"  null calibration: {rep + 1}/{n_reps}")

    abs_S = [abs(f["S"]) for f in null_features]
    r_curves = [f["R_curve"] for f in null_features]
    g_specs = [f["G_spec"] for f in null_features]

    thresholds = {
        "S_q99": float(np.percentile(abs_S, 99)),
        "R_curve_q05": float(np.percentile(r_curves, 5)),
        "R_curve_q95": float(np.percentile(r_curves, 95)),
        "G_spec_q99": float(np.percentile(g_specs, 99)),
    }

    print(f"  Thresholds: {thresholds}")
    return thresholds, null_features


# ── Amplitude calibration ───────────────────────────────────────────

def calibrate_amplitudes(cfg, thresholds):
    comp = cfg["composition"]
    streams = comp["streams"]
    stream_names = list(streams.keys())
    n = comp["n_observations"]
    n_reps = 100

    target_stat = {
        "divergence": "Z_MK",
        "convergence": "rho_env",
        "oscillation": "G_spec",
        "aperiodic": "K_01",
    }

    amp_grid = np.array([0.005, 0.01, 0.015, 0.02, 0.03, 0.04, 0.05,
                          0.07, 0.1, 0.15, 0.2, 0.3, 0.5])

    calibrated = {}

    for cond in ["divergence", "convergence", "oscillation", "aperiodic"]:
        stat_name = target_stat[cond]
        print(f"  Calibrating {cond} (stat={stat_name})...")

        best_amp = amp_grid[0]
        for amp in amp_grid:
            stats = []
            for rep in range(min(n_reps, 30)):
                seed = 90000 + rep
                rng = np.random.default_rng(seed)
                forcing = make_forcing(cond, n, rep) * amp

                # Use just the normal stream as proxy for calibration speed
                x, log_e = generate_stream(streams[stream_names[0]], forcing, rng)
                feats = compute_features(log_e)
                stats.append(abs(feats[stat_name]))

            med_stat = np.median(stats)
            # Target: between null 75th and 90th percentile
            # Use composed threshold as rough guide
            if stat_name == "Z_MK":
                target_lo, target_hi = 1.5, 2.2
            elif stat_name == "rho_env":
                target_lo, target_hi = -0.15, -0.05
                med_stat = np.median([f[stat_name] for f in [compute_features(
                    generate_stream(streams[stream_names[0]],
                                    make_forcing(cond, n, r) * amp,
                                    np.random.default_rng(90000 + r))[1])
                    for r in range(20)]])
            elif stat_name == "G_spec":
                target_lo, target_hi = 0.003, 0.008
            elif stat_name == "K_01":
                target_lo, target_hi = 0.3, 0.55

            if stat_name == "rho_env":
                if target_lo <= med_stat <= target_hi:
                    best_amp = amp
                    break
                if med_stat < target_lo:
                    best_amp = amp
            else:
                if target_lo <= med_stat <= target_hi:
                    best_amp = amp
                    break
                if med_stat > target_hi:
                    best_amp = amp

        calibrated[cond] = float(best_amp)
        print(f"    {cond}: amplitude = {best_amp}")

    return calibrated


# ── Evaluation ──────────────────────────────────────────────────────

def run_evaluation(cfg, thresholds, amplitudes, n_reps=100):
    comp = cfg["composition"]
    streams = comp["streams"]
    stream_names = list(streams.keys())
    n = comp["n_observations"]

    rows = []

    for cond in CONDITIONS:
        amp = amplitudes.get(cond, 0.0)
        print(f"  Evaluating {cond} (amp={amp:.4f})...")

        for rep in range(n_reps):
            seed = 99999 + rep + CONDITIONS.index(cond) * 1000
            rng = np.random.default_rng(seed)

            forcing = make_forcing(cond, n, rep) * amp if cond != "null" else np.zeros(n)

            all_log_e = []
            all_x = []

            for sname in stream_names:
                x, log_e = generate_stream(streams[sname], forcing, rng)
                all_log_e.append(log_e)
                all_x.append(x)

            # Individual streams
            for k, sname in enumerate(stream_names):
                feats = compute_features(all_log_e[k])
                label, kill_log = classify(feats, thresholds)
                rows.append({
                    "condition": cond, "rep": rep,
                    "signal": f"individual_{sname}",
                    "label": label, "kill_log": kill_log,
                    **{f"f_{k}": f"{v:.6f}" if isinstance(v, float) else str(v)
                       for k, v in feats.items()},
                })

            # Composed e-value
            composed = np.sum(all_log_e, axis=0)
            feats = compute_features(composed)
            label, kill_log = classify(feats, thresholds)
            rows.append({
                "condition": cond, "rep": rep,
                "signal": "composed",
                "label": label, "kill_log": kill_log,
                **{f"f_{k}": f"{v:.6f}" if isinstance(v, float) else str(v)
                   for k, v in feats.items()},
            })

            # Standardized sum
            z_scores = [(xi - xi.mean()) / max(xi.std(), 1e-10) for xi in all_x]
            z_sum = np.sum(z_scores, axis=0)
            feats = compute_features(z_sum)
            label, kill_log = classify(feats, thresholds)
            rows.append({
                "condition": cond, "rep": rep,
                "signal": "standardized_sum",
                "label": label, "kill_log": kill_log,
                **{f"f_{k}": f"{v:.6f}" if isinstance(v, float) else str(v)
                   for k, v in feats.items()},
            })

        print(f"    {cond}: done")

    return rows


# ── Confusion matrix ────────────────────────────────────────────────

GROUND_TRUTH = {
    "null": "null",
    "divergence": "divergent",
    "convergence": "convergent",
    "oscillation": "oscillatory",
    "aperiodic": "aperiodic",
}


def confusion_matrix(rows, signal_type):
    matrix = {gt: {pred: 0 for pred in LABELS} for gt in LABELS}
    for r in rows:
        if r["signal"] != signal_type:
            continue
        gt = GROUND_TRUTH[r["condition"]]
        pred = r["label"]
        matrix[gt][pred] += 1
    return matrix


def macro_f1(matrix):
    f1s = []
    for label in LABELS:
        tp = matrix[label][label]
        fp = sum(matrix[gt][label] for gt in LABELS if gt != label)
        fn = sum(matrix[label][pred] for pred in LABELS if pred != label)
        prec = tp / (tp + fp) if (tp + fp) > 0 else 0
        rec = tp / (tp + fn) if (tp + fn) > 0 else 0
        f1 = 2 * prec * rec / (prec + rec) if (prec + rec) > 0 else 0
        f1s.append(f1)
    return float(np.mean(f1s))


def print_confusion(matrix, title):
    print(f"\n{title}")
    print(f"{'':>14}", end="")
    for pred in LABELS:
        print(f"{pred:>12}", end="")
    print()
    for gt in LABELS:
        print(f"{gt:>14}", end="")
        for pred in LABELS:
            print(f"{matrix[gt][pred]:>12}", end="")
        print()
    print(f"  Macro-F1: {macro_f1(matrix):.3f}")


# ── Plotting ────────────────────────────────────────────────────────

def plot_confusion_matrices(rows):
    signal_types = ["composed", "standardized_sum"]
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    for i, sig in enumerate(signal_types):
        mat = confusion_matrix(rows, sig)
        data = np.array([[mat[gt][pred] for pred in LABELS] for gt in LABELS])
        totals = data.sum(axis=1, keepdims=True)
        totals[totals == 0] = 1
        pct = data / totals * 100

        ax = axes[i]
        im = ax.imshow(pct, cmap="Blues", vmin=0, vmax=100)
        ax.set_xticks(range(len(LABELS)))
        ax.set_yticks(range(len(LABELS)))
        ax.set_xticklabels(LABELS, rotation=45, fontsize=8, ha="right")
        ax.set_yticklabels(LABELS, fontsize=8)
        ax.set_xlabel("Predicted")
        ax.set_ylabel("True")
        f1 = macro_f1(mat)
        ax.set_title(f"{sig}\nMacro-F1 = {f1:.3f}", fontsize=10)

        for r in range(len(LABELS)):
            for c in range(len(LABELS)):
                color = "white" if pct[r, c] > 50 else "black"
                ax.text(c, r, f"{pct[r, c]:.0f}%", ha="center", va="center",
                        fontsize=9, color=color)

    fig.colorbar(im, ax=axes, shrink=0.8, label="% of true class")
    fig.suptitle("Four-bin classification: confusion matrices", fontsize=12)
    fig.tight_layout()
    fig.savefig(PLOTS / "fourbin_confusion.png", dpi=150)
    plt.close(fig)
    print("  fourbin_confusion.png")


def plot_kill_flow(rows):
    signal_type = "composed"
    fig, axes = plt.subplots(1, len(CONDITIONS), figsize=(20, 5), sharey=True)

    for i, cond in enumerate(CONDITIONS):
        kills = {}
        for r in rows:
            if r["signal"] == signal_type and r["condition"] == cond:
                kl = r["kill_log"]
                kills[kl] = kills.get(kl, 0) + 1

        ax = axes[i]
        if kills:
            sorted_kills = sorted(kills.items(), key=lambda x: -x[1])
            labels_k = [k[:30] for k, _ in sorted_kills]
            counts = [v for _, v in sorted_kills]
            ax.barh(range(len(labels_k)), counts)
            ax.set_yticks(range(len(labels_k)))
            ax.set_yticklabels(labels_k, fontsize=7)
        ax.set_title(cond, fontsize=9)
        ax.set_xlabel("Count")
        if i == 0:
            ax.set_ylabel("Kill path")

    fig.suptitle("Kill-condition flow (composed signal)", fontsize=12)
    fig.tight_layout()
    fig.savefig(PLOTS / "fourbin_kills.png", dpi=150)
    plt.close(fig)
    print("  fourbin_kills.png")


# ── Main ────────────────────────────────────────────────────────────

def main():
    cfg = load_config()
    PLOTS.mkdir(parents=True, exist_ok=True)
    TABLES.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print("PHASE 1: Null calibration (1000 reps)")
    print("=" * 60)
    thresholds, null_features = run_null_calibration(cfg, n_reps=1000)

    print("\n" + "=" * 60)
    print("PHASE 2: Amplitudes (manual, iterate to find weak-signal regime)")
    print("=" * 60)
    amplitudes = {
        "null": 0.0,
        "divergence": 1.5,
        "convergence": 1.5,
        "oscillation": 0.1,
        "aperiodic": 0.5,
    }
    print(f"  Amplitudes: {amplitudes}")

    print("\n" + "=" * 60)
    print("PHASE 3: Evaluation (100 reps × 5 conditions)")
    print("=" * 60)
    rows = run_evaluation(cfg, thresholds, amplitudes, n_reps=100)

    # Save results
    feature_keys = [k for k in rows[0].keys() if k.startswith("f_")]
    fieldnames = ["condition", "rep", "signal", "label", "kill_log"] + feature_keys
    csv_path = TABLES / "fourbin.csv"
    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    print(f"\n  {csv_path}")

    # Confusion matrices
    for sig in ["composed", "standardized_sum"]:
        mat = confusion_matrix(rows, sig)
        print_confusion(mat, f"Confusion matrix: {sig}")

    # Individual aggregate
    individual_labels = [r for r in rows if r["signal"].startswith("individual_")]
    # Majority vote per (condition, rep)
    from collections import Counter
    maj_rows = []
    for cond in CONDITIONS:
        for rep in range(100):
            labels = [r["label"] for r in individual_labels
                      if r["condition"] == cond and r["rep"] == rep]
            if labels:
                majority = Counter(labels).most_common(1)[0][0]
                maj_rows.append({"condition": cond, "rep": rep,
                                 "signal": "individual_majority", "label": majority})

    mat_maj = {gt: {pred: 0 for pred in LABELS} for gt in LABELS}
    for r in maj_rows:
        gt = GROUND_TRUTH[r["condition"]]
        mat_maj[gt][r["label"]] += 1
    print_confusion(mat_maj, "Confusion matrix: individual majority vote")

    # Plots
    print("\nGenerating plots...")
    plot_confusion_matrices(rows)
    plot_kill_flow(rows)

    print("\nDone.")


if __name__ == "__main__":
    main()
