"""V4 robustness analysis. Second attempt after bug hunt."""

import sys
import csv
import math
import json
import numpy as np
from pathlib import Path
from collections import Counter

sys.path.insert(0, str(Path.home() / "e-value-trajectory" / "src"))
from fourbin import (
    load_config, run_null_calibration, generate_stream, make_forcing,
    compute_features, LABELS, CONDITIONS, GROUND_TRUTH
)

ROOT = Path.home() / "e-value-trajectory"
OUT_DIR = ROOT / "results" / "tables"
PLOT_DIR = ROOT / "results" / "plots"
V3_AMPS = {"null": 0.0, "divergence": 1.5, "convergence": 1.5,
           "oscillation": 0.1, "aperiodic": 0.5}


def classify_v4(features, thresholds, period_min=10, curve_mult=1.0):
    """V4 classifier with tunable knobs. Never delegates to V3."""
    f = features
    th = thresholds
    s_thresh = th["S_q99"]
    r_lo = th["R_curve_q05"] * curve_mult
    g_thresh = th["G_spec_q99"]

    if abs(f["S"]) > s_thresh and abs(f["Z_MK"]) > 2.58 \
            and abs(f["med_last"] - f["med_first"]) > 3 * f["MAD_dx"]:
        if f["R_curve"] < r_lo:
            return "convergent", "monotone->curvature->convergent"
        return "divergent", "monotone->curvature->divergent"

    if f["G_spec"] > g_thresh and f["Q_spec"] > 5:
        if period_min < f["peak_period"] < 10000 / 8:
            return "oscillatory", "periodicity->oscillatory"

    if f["K_01"] > 0.8 and 0.55 <= f["PE"] <= 0.95:
        return "aperiodic", "aperiodic->K01+PE"

    return "null", "null"


def ablated_classify(features, thresholds, ablation=None, period_min=10, curve_mult=1.0):
    """Classify with one component removed."""
    f = dict(features)
    if ablation == "no_monotone":
        f["S"] = 0.0
        f["Z_MK"] = 0.0
    elif ablation == "no_mk":
        f["Z_MK"] = 0.0
    elif ablation == "no_curvature":
        f["R_curve"] = 1.0  # neutral
    elif ablation == "no_period_filter":
        period_min = 0
    elif ablation == "no_01":
        f["K_01"] = 0.0
    elif ablation == "no_pe":
        f["PE"] = 0.0
    return classify_v4(f, thresholds, period_min=period_min, curve_mult=curve_mult)


def macro_f1(true_labels, pred_labels):
    f1s = []
    for lab in LABELS:
        tp = sum(1 for t, p in zip(true_labels, pred_labels) if t == lab and p == lab)
        fp = sum(1 for t, p in zip(true_labels, pred_labels) if t != lab and p == lab)
        fn = sum(1 for t, p in zip(true_labels, pred_labels) if t == lab and p != lab)
        prec = tp / (tp + fp) if (tp + fp) else 0
        rec = tp / (tp + fn) if (tp + fn) else 0
        f1s.append(2 * prec * rec / (prec + rec) if (prec + rec) else 0)
    return float(np.mean(f1s))


def stratified_bootstrap_ci(true_labels, pred_labels, n_boot=1000, seed=None):
    """Bootstrap stratified by true class."""
    if seed is None:
        seed = hash(str(true_labels[:5])) % 2**31
    rng = np.random.default_rng(seed)
    classes = sorted(set(true_labels))
    indices_by_class = {c: [i for i, t in enumerate(true_labels) if t == c] for c in classes}
    f1s = []
    for _ in range(n_boot):
        boot_true, boot_pred = [], []
        for c in classes:
            idx = indices_by_class[c]
            if not idx:
                continue
            sample = rng.choice(idx, len(idx), replace=True)
            boot_true.extend(true_labels[i] for i in sample)
            boot_pred.extend(pred_labels[i] for i in sample)
        f1s.append(macro_f1(boot_true, boot_pred))
    f1s.sort()
    return float(np.median(f1s)), f1s[int(0.025 * n_boot)], f1s[int(0.975 * n_boot)]


def generate_corrupted(cfg, condition, amp, rng, n, degradation=None, deg_params=None):
    """Generate streams with degradation applied at the forcing/noise level."""
    streams = cfg["composition"]["streams"]
    stream_names = list(streams.keys())

    forcing = make_forcing(condition, n, 0) * amp if condition != "null" else np.zeros(n)

    # Autocorrelated: add AR(1) noise to forcing for ALL conditions
    if degradation == "autocorrelated":
        phi = deg_params["phi"]
        ar = np.zeros(n)
        ar[0] = rng.normal()
        for t in range(1, n):
            ar[t] = phi * ar[t-1] + rng.normal() * math.sqrt(1 - phi**2)
        forcing = forcing + 0.1 * ar  # small autocorrelated perturbation

    # Nonstationary: add drift to forcing for ALL conditions
    if degradation == "nonstationary":
        drift = deg_params["drift"]
        forcing = forcing + drift * np.arange(n) / n

    # Correlated: generate shared latent before streams
    shared_latent = None
    if degradation == "correlated":
        shared_latent = rng.normal(0, 1, n)

    all_x, all_log_e = [], []
    for sname in stream_names:
        s_cfg = streams[sname]

        if degradation == "misspecified" and s_cfg["distribution"] == "normal":
            df = deg_params["df"]
            sigma = s_cfg["null_params"]["sigma"]
            scale = sigma * math.sqrt(df / (df - 2)) if df > 2 else sigma
            x = forcing + rng.standard_t(df, n) * scale
            lam = s_cfg["alt_params"]["lambda"]
            log_e = lam * x - lam**2 / 2
        else:
            x, log_e = generate_stream(s_cfg, forcing, rng)

        # Correlated: inject shared latent into forcing before generation
        if degradation == "correlated" and shared_latent is not None:
            rho = deg_params["rho"]
            # Re-generate with correlated forcing
            corr_forcing = forcing * math.sqrt(1 - rho) + shared_latent * math.sqrt(rho) * 0.3
            if degradation == "misspecified" and s_cfg["distribution"] == "normal":
                pass  # already generated above
            else:
                x, log_e = generate_stream(s_cfg, corr_forcing, rng)

        # Missing: set to neutral evidence (log_e = 0), not interpolation
        if degradation == "missing":
            frac = deg_params["frac"]
            mask = rng.random(n) < frac
            log_e = log_e.copy()
            log_e[mask] = 0.0  # neutral evidence for missing observations

        all_x.append(x)
        all_log_e.append(log_e)

    return all_x, all_log_e


def run_cell(cfg, thresholds, category, perturbation, severity_str,
             conditions, amps, n_reps, seed_base,
             degradation=None, deg_params=None,
             ablation=None, period_min=10, curve_mult=1.0,
             mixed_parts=None):
    """Run one cell: generate, compose, classify, score. Returns rows for both composed and standardized."""
    n = cfg["composition"]["n_observations"]
    rows = []

    for cond in conditions:
        amp = amps.get(cond, 0.0)
        for rep in range(n_reps):
            seed = seed_base + rep + CONDITIONS.index(cond) * 10000 if cond in CONDITIONS else seed_base + rep
            rng = np.random.default_rng(seed)

            if mixed_parts:
                forcing = np.zeros(n)
                for mc, mw in mixed_parts:
                    forcing += mw * V3_AMPS[mc] * make_forcing(mc, n, rep)
                all_x, all_log_e = [], []
                for sname in cfg["composition"]["streams"]:
                    x, log_e = generate_stream(cfg["composition"]["streams"][sname], forcing, rng)
                    all_x.append(x)
                    all_log_e.append(log_e)
            else:
                all_x, all_log_e = generate_corrupted(
                    cfg, cond, amp, rng, n, degradation, deg_params)

            # Composed
            composed = np.sum(all_log_e, axis=0)
            feats = compute_features(composed)
            label, kill = ablated_classify(feats, thresholds, ablation, period_min, curve_mult)
            rows.append({"condition": cond, "rep": rep, "signal": "composed",
                         "label": label, "gt": GROUND_TRUTH.get(cond, "mixed")})

            # Standardized sum
            z_scores = [(x - np.nanmean(x)) / max(np.nanstd(x), 1e-10) for x in all_x]
            z_sum = np.sum(z_scores, axis=0)
            feats = compute_features(z_sum)
            label, kill = ablated_classify(feats, thresholds, ablation, period_min, curve_mult)
            rows.append({"condition": cond, "rep": rep, "signal": "standardized_sum",
                         "label": label, "gt": GROUND_TRUTH.get(cond, "mixed")})

    return rows


def summarize_cell(rows, category, perturbation, severity_str, is_mixed=False):
    """Summarize a cell into output rows."""
    results = []
    for signal in ["composed", "standardized_sum"]:
        sig_rows = [r for r in rows if r["signal"] == signal]
        if not sig_rows:
            continue

        if is_mixed:
            counts = Counter(r["label"] for r in sig_rows)
            entropy = -sum((v/len(sig_rows)) * math.log2(v/len(sig_rows))
                           for v in counts.values() if v > 0)
            results.append({
                "category": category, "perturbation": perturbation,
                "severity": severity_str, "signal": signal,
                "macro_f1": "", "ci_lo": "", "ci_hi": "",
                "null_fpr": "", "label_entropy": f"{entropy:.3f}",
                "label_counts": json.dumps(dict(sorted(counts.items()))),
            })
        else:
            true_l = [r["gt"] for r in sig_rows]
            pred_l = [r["label"] for r in sig_rows]
            f1 = macro_f1(true_l, pred_l)
            _, lo, hi = stratified_bootstrap_ci(true_l, pred_l)
            null_pred = [r["label"] for r in sig_rows if r["gt"] == "null"]
            null_fpr = 1 - sum(1 for p in null_pred if p == "null") / max(len(null_pred), 1)
            results.append({
                "category": category, "perturbation": perturbation,
                "severity": severity_str, "signal": signal,
                "macro_f1": f"{f1:.3f}", "ci_lo": f"{lo:.3f}", "ci_hi": f"{hi:.3f}",
                "null_fpr": f"{null_fpr:.3f}", "label_entropy": "",
                "label_counts": "",
            })

    # Paired delta
    comp = [r for r in rows if r["signal"] == "composed" and r["gt"] != "mixed"]
    std = [r for r in rows if r["signal"] == "standardized_sum" and r["gt"] != "mixed"]
    if comp and std and not is_mixed:
        comp_f1 = macro_f1([r["gt"] for r in comp], [r["label"] for r in comp])
        std_f1 = macro_f1([r["gt"] for r in std], [r["label"] for r in std])
        for r in results:
            r["paired_delta"] = f"{comp_f1 - std_f1:+.3f}" if r["signal"] == "composed" else ""
    else:
        for r in results:
            r["paired_delta"] = ""

    return results


def main():
    cfg = load_config()
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print("Null calibration (1000 reps)...")
    thresholds, _ = run_null_calibration(cfg, n_reps=1000)

    all_results = []
    seed = 99999

    # AMPLITUDE SWEEP (per-bin detection rate)
    print("\nAmplitude sweep...")
    detection_by_bin = {}
    for cond in ["divergence", "convergence", "oscillation", "aperiodic"]:
        detection_by_bin[cond] = []
        base = V3_AMPS[cond]
        for step in range(21):
            amp = base * 2.0 * step / 20
            amps = dict(V3_AMPS)
            amps[cond] = amp
            rows = run_cell(cfg, thresholds, "amplitude", f"{cond}_sweep",
                            f"{amp:.4f}", CONDITIONS, amps, 100, seed)
            seed += 50000
            # Per-bin detection for THIS condition
            comp_rows = [r for r in rows if r["signal"] == "composed" and r["gt"] == GROUND_TRUTH[cond]]
            detect = sum(1 for r in comp_rows if r["label"] == r["gt"]) / max(len(comp_rows), 1)
            detection_by_bin[cond].append((amp, detect))
            summary = summarize_cell(rows, "amplitude", f"{cond}_sweep", f"{amp:.4f}")
            for s in summary:
                s["detection_rate"] = f"{detect:.3f}"
            all_results.extend(summary)
            print(f"  {cond} amp={amp:.3f}: detect={detect:.0%}")

    # EQUALIZED DIFFICULTY
    print("\nEqualized difficulty...")
    eq_amps = dict(V3_AMPS)
    for cond in ["divergence", "convergence", "oscillation", "aperiodic"]:
        pairs = detection_by_bin[cond]
        # Find amplitude closest to 80% detection
        best = min(pairs, key=lambda p: abs(p[1] - 0.8))
        eq_amps[cond] = best[0]
        print(f"  {cond}: amp={best[0]:.4f} (detect={best[1]:.0%})")
    rows = run_cell(cfg, thresholds, "amplitude", "equalized_difficulty",
                    "80pct_target", CONDITIONS, eq_amps, 100, seed)
    seed += 50000
    all_results.extend(summarize_cell(rows, "amplitude", "equalized_difficulty", "80pct_target"))

    # ABLATIONS (paired seeds)
    print("\nAblations...")
    abl_seed = 200000
    base_rows = run_cell(cfg, thresholds, "ablation", "baseline", "none",
                         CONDITIONS, V3_AMPS, 100, abl_seed)
    base_f1 = macro_f1([r["gt"] for r in base_rows if r["signal"] == "composed"],
                       [r["label"] for r in base_rows if r["signal"] == "composed"])

    for abl_key, abl_name in [
        ("no_curvature", "Remove curvature"),
        ("no_period_filter", "Remove period>10"),
        ("no_01", "Remove 0-1 test"),
        ("no_pe", "Remove PE"),
        ("no_monotone", "Remove monotone"),
        ("no_mk", "Remove Mann-Kendall"),
    ]:
        rows = run_cell(cfg, thresholds, "ablation", abl_name, abl_key,
                        CONDITIONS, V3_AMPS, 100, abl_seed, ablation=abl_key)
        summary = summarize_cell(rows, "ablation", abl_name, abl_key)
        abl_f1 = macro_f1([r["gt"] for r in rows if r["signal"] == "composed"],
                          [r["label"] for r in rows if r["signal"] == "composed"])
        for s in summary:
            s["paired_delta"] = f"{abl_f1 - base_f1:+.3f}" if s["signal"] == "composed" else ""
        all_results.extend(summary)
        print(f"  {abl_name}: F1={abl_f1:.3f} (delta={abl_f1-base_f1:+.3f})")

    for mult, name in [(0.8, "Curvature -20%"), (1.2, "Curvature +20%")]:
        rows = run_cell(cfg, thresholds, "ablation", name, f"curve_mult={mult}",
                        CONDITIONS, V3_AMPS, 100, abl_seed, curve_mult=mult)
        all_results.extend(summarize_cell(rows, "ablation", name, f"curve_mult={mult}"))

    # DEGRADATIONS (severity sweeps, all 5 conditions)
    print("\nDegradations...")
    degradations = [
        ("autocorrelated", "phi", [0.2, 0.5, 0.8]),
        ("correlated", "rho", [0.1, 0.3, 0.6]),
        ("missing", "frac", [0.05, 0.10, 0.25]),
        ("misspecified", "df", [10, 5, 3]),
        ("nonstationary", "drift", [0.005, 0.01, 0.02]),
    ]
    for deg_name, param_name, severities in degradations:
        for sev in severities:
            rows = run_cell(cfg, thresholds, "degradation", deg_name,
                            f"{param_name}={sev}", CONDITIONS, V3_AMPS, 500, seed,
                            degradation=deg_name, deg_params={param_name: sev})
            seed += 250000
            summary = summarize_cell(rows, "degradation", deg_name, f"{param_name}={sev}")
            all_results.extend(summary)
            comp_f1 = [s["macro_f1"] for s in summary if s["signal"] == "composed"]
            print(f"  {deg_name} {param_name}={sev}: F1={comp_f1[0] if comp_f1 else '?'}")

    # PERIOD FILTER SWEEP
    print("\nPeriod filter sweep...")
    for min_p in range(2, 52, 2):
        rows = run_cell(cfg, thresholds, "filter", "period_threshold",
                        str(min_p), ["oscillation", "aperiodic"], V3_AMPS, 100, seed,
                        period_min=min_p)
        seed += 20000
        comp = [r for r in rows if r["signal"] == "composed"]
        aper_ok = sum(1 for r in comp if r["gt"] == "aperiodic" and r["label"] == "aperiodic") / 100
        osc_ok = sum(1 for r in comp if r["gt"] == "oscillatory" and r["label"] == "oscillatory") / 100
        all_results.append({
            "category": "filter", "perturbation": "period_threshold",
            "severity": str(min_p), "signal": "composed",
            "macro_f1": "", "ci_lo": "", "ci_hi": "",
            "null_fpr": "", "label_entropy": "",
            "label_counts": f"aper={aper_ok:.2f},osc={osc_ok:.2f}",
            "detection_rate": "", "paired_delta": "",
        })
        print(f"  period>{min_p}: aper={aper_ok:.0%}, osc={osc_ok:.0%}")

    # MIXED DYNAMICS (descriptive)
    print("\nMixed dynamics...")
    mixed_specs = [
        ("trend+oscillation", [("divergence", None), ("oscillation", None)]),
        ("oscillation+aperiodic", [("oscillation", None), ("aperiodic", None)]),
        ("convergence+oscillation", [("convergence", None), ("oscillation", None)]),
        ("divergence+heavy_tail", [("divergence", None)]),
    ]
    ratios = [0.2, 0.4, 0.5, 0.6, 0.8]
    for pair_name, parts_spec in mixed_specs:
        for r in ratios:
            if pair_name == "divergence+heavy_tail":
                # Heavy tail: use t-distributed noise as the second component
                mixed_parts = [("divergence", r)]
                # Generate with t-noise degradation for the non-divergence portion
                rows_list = []
                for rep in range(100):
                    rng = np.random.default_rng(seed + rep)
                    n = cfg["composition"]["n_observations"]
                    forcing = make_forcing("divergence", n, rep) * V3_AMPS["divergence"] * r
                    forcing += rng.standard_t(3, n) * 0.3 * (1 - r)
                    all_x, all_log_e = [], []
                    for sname in cfg["composition"]["streams"]:
                        x, log_e = generate_stream(cfg["composition"]["streams"][sname], forcing, rng)
                        all_x.append(x)
                        all_log_e.append(log_e)
                    composed = np.sum(all_log_e, axis=0)
                    feats = compute_features(composed)
                    label, _ = classify_v4(feats, thresholds)
                    rows_list.append({"condition": "mixed", "rep": rep, "signal": "composed",
                                      "label": label, "gt": "mixed"})
                seed += 100
                summary = summarize_cell(rows_list, "mixed", pair_name, f"{r:.1f}:{1-r:.1f}", is_mixed=True)
            else:
                cond1, cond2 = parts_spec[0][0], parts_spec[1][0]
                mixed_parts = [(cond1, r), (cond2, 1-r)]
                rows = run_cell(cfg, thresholds, "mixed", pair_name,
                                f"{r:.1f}:{1-r:.1f}", ["null"], V3_AMPS, 100, seed,
                                mixed_parts=mixed_parts)
                seed += 10000
                summary = summarize_cell(rows, "mixed", pair_name, f"{r:.1f}:{1-r:.1f}", is_mixed=True)
            all_results.extend(summary)
            counts = [s["label_counts"] for s in summary if s["signal"] == "composed"]
            print(f"  {pair_name} {r:.1f}:{1-r:.1f}: {counts[0] if counts else '?'}")

    # E-VALUE VALIDITY
    print("\nE-value validity...")
    exceedance_count = 0
    max_mean_et = 0.0
    validity_reps = 1000
    for rep in range(validity_reps):
        rng = np.random.default_rng(seed + rep)
        n = cfg["composition"]["n_observations"]
        forcing = np.zeros(n)
        all_log_e = []
        for sname in cfg["composition"]["streams"]:
            s_cfg = cfg["composition"]["streams"][sname]
            if s_cfg["distribution"] == "normal":
                df = 3
                sigma = s_cfg["null_params"]["sigma"]
                scale = sigma * math.sqrt(df / (df - 2))
                x = rng.standard_t(df, n) * scale
                lam = s_cfg["alt_params"]["lambda"]
                log_e = lam * x - lam**2 / 2
            else:
                x, log_e = generate_stream(s_cfg, forcing, rng)
            all_log_e.append(log_e)
        composed = np.sum(all_log_e, axis=0)
        cum_log_e = np.cumsum(composed)
        max_E = np.exp(np.max(cum_log_e))
        if max_E > 20:
            exceedance_count += 1
        # Track mean e_t at each time
        if rep == 0:
            mean_et_accum = np.exp(composed)
        else:
            mean_et_accum += np.exp(composed)
    seed += validity_reps

    type_i_rate = exceedance_count / validity_reps
    mean_et = mean_et_accum / validity_reps
    max_mean = float(np.max(mean_et))
    all_results.append({
        "category": "validity", "perturbation": "t_null_type_i",
        "severity": "df=3", "signal": "composed",
        "macro_f1": "", "ci_lo": "", "ci_hi": "",
        "null_fpr": f"type_i={type_i_rate:.3f}",
        "label_entropy": "",
        "label_counts": f"max_mean_et={max_mean:.4f},exceedance={exceedance_count}/{validity_reps}",
        "detection_rate": "", "paired_delta": "",
    })
    print(f"  Type-I rate: {type_i_rate:.1%}, max mean(e_t): {max_mean:.4f}")

    # SAVE
    fieldnames = ["category", "perturbation", "severity", "signal",
                  "macro_f1", "ci_lo", "ci_hi", "null_fpr",
                  "label_entropy", "label_counts", "detection_rate", "paired_delta"]
    out_path = OUT_DIR / "v4_grid.csv"
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(all_results)
    print(f"\nSaved {len(all_results)} rows to {out_path}")
    print("Done.")


if __name__ == "__main__":
    main()
