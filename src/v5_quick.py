"""V5 quick tests: H24, H25, H6, H7 — the four easiest hypotheses."""

import sys
import numpy as np
from pathlib import Path

sys.path.insert(0, str(Path.home() / "e-value-trajectory" / "src"))
from fourbin import (
    load_config, run_null_calibration, generate_stream, make_forcing,
    compute_features, LABELS, CONDITIONS, GROUND_TRUTH
)

V3_AMPS = {"null": 0.0, "divergence": 1.5, "convergence": 1.5,
           "oscillation": 0.1, "aperiodic": 0.5}


def classify_full(f, th):
    """V3 classifier — 7 features, original ordering."""
    if (abs(f["S"]) > th["S_q99"] and abs(f["Z_MK"]) > 2.58
            and abs(f["med_last"] - f["med_first"]) > 3 * f["MAD_dx"]):
        if f["R_curve"] < th["R_curve_q05"]:
            return "convergent"
        return "divergent"
    if f["G_spec"] > th["G_spec_q99"] and f["Q_spec"] > 5:
        if 10 < f["peak_period"] < 10000 / 8:
            return "oscillatory"
    if f["K_01"] > 0.8 and 0.55 <= f["PE"] <= 0.95:
        return "aperiodic"
    return "null"


def classify_minimal(f, th):
    """H25: 3-feature classifier — no MK, no K_01, no PE. Just slope, curvature, spectral."""
    if (abs(f["S"]) > th["S_q99"]
            and abs(f["med_last"] - f["med_first"]) > 3 * f["MAD_dx"]):
        if f["R_curve"] < th["R_curve_q05"]:
            return "convergent"
        return "divergent"
    if f["G_spec"] > th["G_spec_q99"] and f["Q_spec"] > 5:
        if 10 < f["peak_period"] < 10000 / 8:
            return "oscillatory"
    if f["PE"] >= 0.55 and f["PE"] <= 0.95:
        return "aperiodic"
    return "null"


def classify_damage_order(f, th):
    """H24: reorder by ablation damage. Monotone first (0.726), then spectral (0.264), then curvature (0.267), then aperiodic."""
    # Monotone first (highest damage when removed)
    is_monotone = (abs(f["S"]) > th["S_q99"]
                   and abs(f["med_last"] - f["med_first"]) > 3 * f["MAD_dx"])
    if is_monotone:
        # Curvature to separate convergent/divergent
        if f["R_curve"] < th["R_curve_q05"]:
            return "convergent"
        return "divergent"
    # Spectral next (second highest damage)
    if f["G_spec"] > th["G_spec_q99"] and f["Q_spec"] > 5:
        if 10 < f["peak_period"] < 10000 / 8:
            return "oscillatory"
    # Aperiodic last
    if f["K_01"] > 0.8 and 0.55 <= f["PE"] <= 0.95:
        return "aperiodic"
    return "null"


def macro_f1(true_labels, pred_labels):
    f1s = []
    for lab in LABELS:
        tp = sum(1 for t, p in zip(true_labels, pred_labels) if t == lab and p == lab)
        fp = sum(1 for t, p in zip(true_labels, pred_labels) if t != lab and p == lab)
        fn = sum(1 for t, p in zip(true_labels, pred_labels) if t == lab and p != lab)
        prec = tp / (tp + fp) if (tp + fp) else 0
        rec = tp / (tp + fn) if (tp + fn) else 0
        f1 = 2 * prec * rec / (prec + rec) if (prec + rec) else 0
        f1s.append(f1)
    return float(np.mean(f1s))


def run_test(cfg, thresholds, classifier_fn, composition="product",
             missing_impute="zero", rho=0.0, n_reps=100):
    """Run one test configuration."""
    comp = cfg["composition"]
    streams = comp["streams"]
    stream_names = list(streams.keys())
    n = comp["n_observations"]

    true_labels = []
    pred_labels = []

    for cond in CONDITIONS:
        amp = V3_AMPS.get(cond, 0.0)
        for rep in range(n_reps):
            seed = 99999 + rep + CONDITIONS.index(cond) * 10000
            rng = np.random.default_rng(seed)

            forcing = make_forcing(cond, n, rep) * amp if cond != "null" else np.zeros(n)

            # Correlated streams
            z_shared = None
            if rho > 0:
                z_shared = rng.normal(0, 1, n)

            all_log_e = []
            all_x = []
            for sname in stream_names:
                if rho > 0:
                    forcing_corr = forcing + rho * 0.3 * z_shared
                else:
                    forcing_corr = forcing
                x, log_e = generate_stream(streams[sname], forcing_corr, rng)
                all_x.append(x)
                all_log_e.append(log_e)

            log_e_arr = np.array(all_log_e)
            x_arr = np.array(all_x)

            # Missing data
            if missing_impute != "none":
                frac = 0.10  # test at 10% missing
                for k in range(len(stream_names)):
                    mask = rng.random(n) < frac
                    if missing_impute == "zero":
                        log_e_arr[k][mask] = 0.0
                    elif missing_impute == "weighted":
                        # Weight surviving observations by 1/(1-frac)
                        log_e_arr[k][~mask] *= 1.0 / (1.0 - frac)
                        log_e_arr[k][mask] = 0.0

            # Composition
            if composition == "product":
                composed = np.sum(log_e_arr, axis=0)  # sum of log = log of product
            elif composition == "average":
                # Arithmetic mean of e-values, then take log
                # mean(e_i) = mean(exp(log_e_i))
                # For numerical stability, work in log space:
                # log(mean(exp(x))) = logsumexp(x) - log(K)
                from scipy.special import logsumexp
                composed = np.array([
                    logsumexp(log_e_arr[:, t]) - np.log(len(stream_names))
                    for t in range(n)
                ])

            feats = compute_features(composed)
            label = classifier_fn(feats, thresholds)
            true_labels.append(GROUND_TRUTH[cond])
            pred_labels.append(label)

    return macro_f1(true_labels, pred_labels)


def main():
    cfg = load_config()
    print("Null calibration (1000 reps)...")
    thresholds, _ = run_null_calibration(cfg, n_reps=1000)

    print("\n" + "=" * 60)
    print("V5 QUICK TESTS")
    print("=" * 60)

    # Baseline: V3 classifier, no degradation
    print("\n--- Baseline (V3 conditions, no degradation) ---")
    f1_base = run_test(cfg, thresholds, classify_full, missing_impute="none")
    print(f"  V3 full classifier:     F1={f1_base:.4f}")

    # H25: 3-feature classifier
    f1_min = run_test(cfg, thresholds, classify_minimal, missing_impute="none")
    print(f"  H25 minimal (3-feat):   F1={f1_min:.4f}")

    # H24: damage-ordered classifier
    f1_dmg = run_test(cfg, thresholds, classify_damage_order, missing_impute="none")
    print(f"  H24 damage-ordered:     F1={f1_dmg:.4f}")

    # H6: weighted imputation at 10% missing
    print("\n--- Missing data (10%) ---")
    f1_zero = run_test(cfg, thresholds, classify_full, missing_impute="zero")
    print(f"  Zero imputation:        F1={f1_zero:.4f}")
    f1_weighted = run_test(cfg, thresholds, classify_full, missing_impute="weighted")
    print(f"  H6 weighted imputation: F1={f1_weighted:.4f}")

    # H7: averaging vs product under correlation
    print("\n--- Correlation (rho=0.6) ---")
    f1_prod_corr = run_test(cfg, thresholds, classify_full, rho=0.6, missing_impute="none")
    print(f"  Product composition:    F1={f1_prod_corr:.4f}")
    f1_avg_corr = run_test(cfg, thresholds, classify_full, composition="average",
                           rho=0.6, missing_impute="none")
    print(f"  H7 average composition: F1={f1_avg_corr:.4f}")

    # H7 bonus: averaging under independence (cost of safety)
    print("\n--- Independence (rho=0) ---")
    f1_avg_ind = run_test(cfg, thresholds, classify_full, composition="average",
                          missing_impute="none")
    print(f"  Average (independent):  F1={f1_avg_ind:.4f}")
    print(f"  Product (independent):  F1={f1_base:.4f}")

    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"  H24 (damage order):     {'SAME' if f1_dmg == f1_base else f'DIFFERENT ({f1_dmg:.4f} vs {f1_base:.4f})'}")
    print(f"  H25 (3-feature):        {'SAME' if f1_min == f1_base else f'DIFFERENT ({f1_min:.4f} vs {f1_base:.4f})'}")
    print(f"  H6  (weighted @ 10%):   {'BETTER' if f1_weighted > f1_zero else 'SAME/WORSE'} ({f1_weighted:.4f} vs {f1_zero:.4f})")
    print(f"  H7  (avg @ rho=0.6):    {'BETTER' if f1_avg_corr > f1_prod_corr else 'SAME/WORSE'} ({f1_avg_corr:.4f} vs {f1_prod_corr:.4f})")
    print(f"  H7  cost (avg @ rho=0): {'SAME' if f1_avg_ind == f1_base else f'LOST ({f1_avg_ind:.4f} vs {f1_base:.4f})'}")


if __name__ == "__main__":
    main()
