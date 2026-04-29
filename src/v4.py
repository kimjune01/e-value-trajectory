"""V4 robustness experiment for four-bin e-value trajectory classification.

Pre-registered at PREREGISTRATION_V4.md. Implements the full perturbation grid:
amplitude sweeps, ablations, degradations, mixed dynamics, E-value validity,
and standardized-sum baseline comparisons.

Run:  cd ~/e-value-trajectory && .venv/bin/python src/v4.py
"""

import hashlib
import json
import csv
import sys
from pathlib import Path
from itertools import product as iterproduct
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

# ── Imports from fourbin.py ─────────────────────────────────────────
ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT / "src"))

from fourbin import (
    load_config,
    run_null_calibration,
    generate_stream,
    make_forcing,
    compute_features,
    LABELS,
    CONDITIONS,
    GROUND_TRUTH,
)

RESULTS = ROOT / "results"
PLOTS = RESULTS / "plots"
TABLES = RESULTS / "tables"

# ── V3 base amplitudes (from fourbin.py main) ──────────────────────
BASE_AMPLITUDES = {
    "null": 0.0,
    "divergence": 1.5,
    "convergence": 1.5,
    "oscillation": 0.1,
    "aperiodic": 0.5,
}


# ====================================================================
# V4 CLASSIFIER (reimplemented with tunable knobs)
# ====================================================================

def classify_v4(features, thresholds, period_min=10, curve_mult=1.0):
    """Decision tree identical to V3 but with period_min and curve_mult knobs."""
    f = features
    th = thresholds
    N = 10000  # normative

    # TEST 1: Monotone trend
    if (abs(f["S"]) > th["S_q99"]
            and abs(f["Z_MK"]) > 2.58
            and abs(f["med_last"] - f["med_first"]) > 3 * f["MAD_dx"]):
        # TEST 2: Curvature
        if f["R_curve"] < th["R_curve_q05"] * curve_mult:
            return "convergent"
        else:
            return "divergent"

    # TEST 3: Periodicity
    if f["G_spec"] > th["G_spec_q99"] and f["Q_spec"] > 5:
        period = f["peak_period"]
        if period_min < period < (float("inf") if period_min == 0 else N / 8):
            return "oscillatory"

    # TEST 4: Aperiodic structure
    if f["K_01"] > 0.8 and 0.55 <= f["PE"] <= 0.95:
        return "aperiodic"

    return "null"


# ====================================================================
# MACRO-F1 AND BOOTSTRAP
# ====================================================================

def compute_macro_f1(y_true, y_pred):
    """Standard macro-averaged F1 over LABELS."""
    f1s = []
    for label in LABELS:
        tp = sum(1 for t, p in zip(y_true, y_pred) if t == label and p == label)
        fp = sum(1 for t, p in zip(y_true, y_pred) if t != label and p == label)
        fn = sum(1 for t, p in zip(y_true, y_pred) if t == label and p != label)
        prec = tp / (tp + fp) if (tp + fp) > 0 else 0
        rec = tp / (tp + fn) if (tp + fn) > 0 else 0
        f1 = 2 * prec * rec / (prec + rec) if (prec + rec) > 0 else 0
        f1s.append(f1)
    return float(np.mean(f1s))


def bootstrap_f1_ci(y_true, y_pred, category, perturbation, severity, signal,
                     n_boot=1000):
    """Stratified percentile bootstrap for macro-F1. Returns (point, ci_lo, ci_hi)."""
    point = compute_macro_f1(y_true, y_pred)

    # Group indices by true class (stratified)
    class_indices = {}
    for i, t in enumerate(y_true):
        class_indices.setdefault(t, []).append(i)

    boot_f1s = []
    for b in range(n_boot):
        seed_str = f"{category},{perturbation},{severity},{signal},{b}"
        seed = int(hashlib.sha256(seed_str.encode()).hexdigest()[:8], 16)
        rng = np.random.default_rng(seed)

        boot_idx = []
        for cls in LABELS:
            idx = class_indices.get(cls, [])
            if len(idx) > 0:
                boot_idx.extend(rng.choice(idx, size=len(idx), replace=True))

        bt = [y_true[i] for i in boot_idx]
        bp = [y_pred[i] for i in boot_idx]
        boot_f1s.append(compute_macro_f1(bt, bp))

    ci_lo = float(np.percentile(boot_f1s, 2.5))
    ci_hi = float(np.percentile(boot_f1s, 97.5))
    return point, ci_lo, ci_hi


def null_fpr_from_preds(y_true, y_pred):
    """False positive rate for null class: fraction of null trajectories classified
    as non-null."""
    null_mask = [t == "null" for t in y_true]
    if not any(null_mask):
        return float("nan")
    n_null = sum(null_mask)
    fp = sum(1 for m, p in zip(null_mask, y_pred) if m and p != "null")
    return fp / n_null


def label_distribution(y_pred):
    """Returns dict of label -> count and Shannon entropy (log base 2)."""
    counts = {l: 0 for l in LABELS}
    for p in y_pred:
        counts[p] = counts.get(p, 0) + 1
    total = len(y_pred)
    if total == 0:
        return counts, 0.0
    entropy = 0.0
    for c in counts.values():
        if c > 0:
            p = c / total
            entropy -= p * np.log2(p)
    return counts, entropy


# ====================================================================
# SEED SCHEDULE
# ====================================================================

def eval_seed(rep, condition):
    """Evaluation seed per spec: 99999 + rep + CONDITIONS.index(condition) * 10000."""
    return 99999 + rep + CONDITIONS.index(condition) * 10000


# ====================================================================
# CORE EVALUATION LOOP
# ====================================================================

def run_single_rep(cfg, condition, rep, thresholds, thresholds_zsum,
                   amplitude_overrides=None, degradation=None, deg_params=None,
                   period_min=10, curve_mult=1.0, feature_overrides=None,
                   mixed_forcing_fn=None):
    """Run one rep for one condition. Returns (composed_pred, zsum_pred,
    composed_features, zsum_features, null_fpr_data, log_e_composed, x_streams).

    amplitude_overrides: dict condition -> amplitude (overrides BASE_AMPLITUDES)
    degradation: one of None, "ar1", "correlated", "missing", "t_misspec", "drift"
    deg_params: dict of degradation parameters
    feature_overrides: dict of feature name -> value (for ablations)
    mixed_forcing_fn: callable(n, rep, rng) -> forcing array (for mixed dynamics)
    """
    comp = cfg["composition"]
    streams = comp["streams"]
    stream_names = list(streams.keys())
    N = comp["n_observations"]

    seed = eval_seed(rep, condition)
    rng = np.random.default_rng(seed)

    # -- Determine amplitude --
    if amplitude_overrides and condition in amplitude_overrides:
        amp = amplitude_overrides[condition]
    else:
        amp = BASE_AMPLITUDES.get(condition, 0.0)

    # -- Build forcing --
    if mixed_forcing_fn is not None:
        forcing = mixed_forcing_fn(N, rep, rng)
    else:
        forcing = make_forcing(condition, N, rep) * amp

    # -- Apply degradation to forcing (pre-stream) --
    if degradation == "ar1" and deg_params:
        phi = deg_params["phi"]
        eps = rng.normal(0, 1, N)
        ar = np.zeros(N)
        scale = np.sqrt(1 - phi**2)
        for t in range(1, N):
            ar[t] = phi * ar[t - 1] + eps[t] * scale
        forcing = forcing + 0.1 * ar

    if degradation == "drift" and deg_params:
        drift_coef = deg_params["drift_coef"]
        t_arr = np.arange(N)
        forcing = forcing + drift_coef * t_arr / N

    # -- Correlated streams: shared latent --
    z_shared = None
    if degradation == "correlated" and deg_params:
        z_shared = rng.normal(0, 1, N)

    # -- Step 3: Generate ALL streams in order (normative RNG sequence) --
    all_log_e = []
    all_x = []
    for sname in stream_names:
        if degradation == "correlated" and deg_params:
            rho = deg_params["rho"]
            forcing_corr = forcing + rho * 0.3 * z_shared
        else:
            forcing_corr = forcing
        x, log_e = generate_stream(streams[sname], forcing_corr, rng)
        all_x.append(x)
        all_log_e.append(log_e)

    log_e_arr = np.array(all_log_e)  # (K, N)
    x_arr = np.array(all_x)          # (K, N)

    # -- Step 4: Missing data (MCAR masks drawn AFTER all streams) --
    if degradation == "missing" and deg_params:
        frac = deg_params["frac"]
        for k in range(len(stream_names)):
            mask = rng.random(N) < frac
            log_e_arr[k][mask] = 0.0      # neutral evidence
            x_arr[k][mask] = np.nan        # NaN for standardized sum

    # -- Step 5: t-noise replacement (LAST, after streams + masks) --
    if degradation == "t_misspec" and deg_params:
        temp_idx = stream_names.index("temperature")
        df = deg_params["df"]
        sigma = streams["temperature"]["null_params"]["sigma"]
        lam = streams["temperature"]["alt_params"]["lambda"]
        t_noise = rng.standard_t(df, N)
        if df > 2:
            t_noise *= np.sqrt((df - 2) / df)
        if degradation == "correlated" and deg_params:
            rho = deg_params["rho"]
            f_corr = forcing + rho * 0.3 * z_shared
        else:
            f_corr = forcing
        x_new = f_corr + t_noise * sigma
        x_arr[temp_idx] = x_new
        log_e_arr[temp_idx] = lam * x_new - lam**2 / 2

    # -- Composed e-value --
    composed = np.sum(log_e_arr, axis=0)
    feats_composed = compute_features(composed)

    # -- Apply feature overrides (ablations) --
    if feature_overrides:
        for k, v in feature_overrides.items():
            feats_composed[k] = v

    pred_composed = classify_v4(feats_composed, thresholds,
                                period_min=period_min, curve_mult=curve_mult)

    # -- Standardized sum --
    z_scores = []
    for k in range(len(stream_names)):
        xk = x_arr[k]
        mu = np.nanmean(xk)
        sd = np.nanstd(xk)
        if sd < 1e-10:
            sd = 1e-10
        zk = (xk - mu) / sd
        # Missing: fill NaN z-scores with 0
        zk = np.where(np.isnan(zk), 0.0, zk)
        z_scores.append(zk)
    z_sum = np.sum(z_scores, axis=0)
    feats_zsum = compute_features(z_sum)

    # Apply same ablation overrides to zsum features
    if feature_overrides:
        for k, v in feature_overrides.items():
            feats_zsum[k] = v

    pred_zsum = classify_v4(feats_zsum, thresholds_zsum,
                            period_min=period_min, curve_mult=curve_mult)

    return pred_composed, pred_zsum, feats_composed, feats_zsum, composed, x_arr


# ====================================================================
# NULL CALIBRATION FOR STANDARDIZED SUM
# ====================================================================

def run_null_calibration_zsum(cfg, n_reps=1000):
    """Calibrate thresholds on standardized sum under null."""
    comp = cfg["composition"]
    streams = comp["streams"]
    stream_names = list(streams.keys())
    N = comp["n_observations"]

    print(f"Null calibration (standardized sum): {n_reps} reps...")
    null_features = []

    for rep in range(n_reps):
        seed = 80000 + rep
        rng = np.random.default_rng(seed)
        forcing = np.zeros(N)

        all_x = []
        for sname in stream_names:
            x, _ = generate_stream(streams[sname], forcing, rng)
            all_x.append(x)

        z_scores = []
        for xk in all_x:
            mu = np.mean(xk)
            sd = max(np.std(xk), 1e-10)
            z_scores.append((xk - mu) / sd)
        z_sum = np.sum(z_scores, axis=0)

        feats = compute_features(z_sum)
        null_features.append(feats)

        if (rep + 1) % 200 == 0:
            print(f"  zsum null calibration: {rep + 1}/{n_reps}")

    abs_S = [abs(f["S"]) for f in null_features]
    r_curves = [f["R_curve"] for f in null_features]
    g_specs = [f["G_spec"] for f in null_features]

    thresholds = {
        "S_q99": float(np.percentile(abs_S, 99)),
        "R_curve_q05": float(np.percentile(r_curves, 5)),
        "R_curve_q95": float(np.percentile(r_curves, 95)),
        "G_spec_q99": float(np.percentile(g_specs, 99)),
    }

    print(f"  Z-sum thresholds: {thresholds}")
    return thresholds


# ====================================================================
# PERTURBATION RUNNERS
# ====================================================================

def run_grid_cell(cfg, thresholds, thresholds_zsum,
                  category, perturbation, severity,
                  n_reps, amplitude_overrides=None,
                  degradation=None, deg_params=None,
                  period_min=10, curve_mult=1.0,
                  feature_overrides=None,
                  mixed_forcing_fn=None,
                  is_mixed=False,
                  sweep_condition=None,
                  return_raw=False):
    """Run one cell of the perturbation grid across all conditions.
    Returns a list of result dicts (one per signal type).
    If return_raw=True, also returns (y_true_comp, y_pred_comp) for paired comparisons."""

    y_true_comp = []
    y_pred_comp = []
    y_true_zsum = []
    y_pred_zsum = []
    sweep_correct = 0
    sweep_total = 0

    conditions_to_run = CONDITIONS
    if is_mixed:
        # Mixed dynamics: just run "null" condition slot but with custom forcing
        conditions_to_run = ["null"]  # placeholder, no ground truth

    for cond in conditions_to_run:
        gt = GROUND_TRUTH[cond]
        for rep in range(n_reps):
            pc, pz, _, _, _, _ = run_single_rep(
                cfg, cond, rep, thresholds, thresholds_zsum,
                amplitude_overrides=amplitude_overrides,
                degradation=degradation, deg_params=deg_params,
                period_min=period_min, curve_mult=curve_mult,
                feature_overrides=feature_overrides,
                mixed_forcing_fn=mixed_forcing_fn,
            )
            if not is_mixed:
                y_true_comp.append(gt)
                y_pred_comp.append(pc)
                y_true_zsum.append(gt)
                y_pred_zsum.append(pz)

                # Sweep detection rate
                if sweep_condition and cond == sweep_condition:
                    if pc == gt:
                        sweep_correct += 1
                    sweep_total += 1
            else:
                y_pred_comp.append(pc)

    results = []

    if is_mixed:
        counts, entropy = label_distribution(y_pred_comp)
        results.append({
            "category": category,
            "perturbation": perturbation,
            "severity": str(severity),
            "signal": "composed",
            "macro_f1": "",
            "ci_lo": "",
            "ci_hi": "",
            "null_fpr": "",
            "label_entropy": f"{entropy:.4f}",
            "label_counts": json.dumps(counts),
            "detection_rate": "",
            "paired_delta": "",
        })
        return results

    # Composed
    f1_c, ci_lo_c, ci_hi_c = bootstrap_f1_ci(
        y_true_comp, y_pred_comp,
        category, perturbation, str(severity), "composed")
    nfpr_c = null_fpr_from_preds(y_true_comp, y_pred_comp)
    counts_c, entropy_c = label_distribution(y_pred_comp)

    # Standardized sum
    f1_z, ci_lo_z, ci_hi_z = bootstrap_f1_ci(
        y_true_zsum, y_pred_zsum,
        category, perturbation, str(severity), "standardized_sum")
    nfpr_z = null_fpr_from_preds(y_true_zsum, y_pred_zsum)
    counts_z, entropy_z = label_distribution(y_pred_zsum)

    det_rate = f"{sweep_correct / sweep_total:.4f}" if sweep_total > 0 else ""
    paired_delta = f"{f1_c - f1_z:.4f}"

    results.append({
        "category": category,
        "perturbation": perturbation,
        "severity": str(severity),
        "signal": "composed",
        "macro_f1": f"{f1_c:.4f}",
        "ci_lo": f"{ci_lo_c:.4f}",
        "ci_hi": f"{ci_hi_c:.4f}",
        "null_fpr": f"{nfpr_c:.4f}",
        "label_entropy": f"{entropy_c:.4f}",
        "label_counts": json.dumps(counts_c),
        "detection_rate": det_rate,
        "paired_delta": paired_delta,
    })
    results.append({
        "category": category,
        "perturbation": perturbation,
        "severity": str(severity),
        "signal": "standardized_sum",
        "macro_f1": f"{f1_z:.4f}",
        "ci_lo": f"{ci_lo_z:.4f}",
        "ci_hi": f"{ci_hi_z:.4f}",
        "null_fpr": f"{nfpr_z:.4f}",
        "label_entropy": f"{entropy_z:.4f}",
        "label_counts": json.dumps(counts_z),
        "detection_rate": "",
        "paired_delta": "",
    })
    if return_raw:
        return results, y_true_comp, y_pred_comp
    return results


# ====================================================================
# AMPLITUDE SWEEP
# ====================================================================

def run_amplitude_sweep(cfg, thresholds, thresholds_zsum):
    """Sweep each non-null bin from 0 to 2x V3 amplitude, 21 steps."""
    print("\n--- Amplitude sweep ---")
    all_rows = []
    sweep_data = {}  # cond -> list of (amp, det_rate)

    non_null = [c for c in CONDITIONS if c != "null"]
    for cond in non_null:
        base_amp = BASE_AMPLITUDES[cond]
        amps = np.linspace(0, 2 * base_amp, 21)
        sweep_data[cond] = []

        for amp in amps:
            overrides = dict(BASE_AMPLITUDES)
            overrides[cond] = amp
            # Other bins stay at V3 amplitude

            rows = run_grid_cell(
                cfg, thresholds, thresholds_zsum,
                category="amplitude", perturbation=f"sweep_{cond}",
                severity=f"{amp:.4f}", n_reps=100,
                amplitude_overrides=overrides,
                sweep_condition=cond,
            )
            all_rows.extend(rows)

            # Extract detection rate
            comp_row = [r for r in rows if r["signal"] == "composed"][0]
            det_str = comp_row["detection_rate"]
            det_rate = float(det_str) if det_str else 0.0
            sweep_data[cond].append((amp, det_rate))

            print(f"  {cond} amp={amp:.4f}: det={det_rate:.3f}  F1={comp_row['macro_f1']}")

    return all_rows, sweep_data


# ====================================================================
# EQUALIZED DIFFICULTY
# ====================================================================

def run_equalized_difficulty(cfg, thresholds, thresholds_zsum, sweep_data):
    """Find amplitude closest to 80% detection per bin, run at those amplitudes."""
    print("\n--- Equalized difficulty ---")
    eq_amps = {}
    for cond, data in sweep_data.items():
        # Find amplitude closest to 80% detection
        best_amp = data[0][0]
        best_diff = abs(data[0][1] - 0.8)
        for amp, det in data:
            diff = abs(det - 0.8)
            if diff < best_diff:
                best_diff = diff
                best_amp = amp
        eq_amps[cond] = best_amp
        print(f"  {cond}: equalized amp = {best_amp:.4f}")

    overrides = {"null": 0.0}
    overrides.update(eq_amps)

    rows = run_grid_cell(
        cfg, thresholds, thresholds_zsum,
        category="amplitude", perturbation="equalized_difficulty",
        severity="1", n_reps=100,
        amplitude_overrides=overrides,
    )
    return rows


# ====================================================================
# ABLATIONS
# ====================================================================

def paired_bootstrap_delta_ci(y_true_base, y_pred_base, y_true_abl, y_pred_abl,
                               category, perturbation, n_boot=1000):
    """Paired bootstrap CI on F1 delta. Resamples by shared (condition, rep)."""
    class_indices = {}
    for i, t in enumerate(y_true_base):
        class_indices.setdefault(t, []).append(i)

    deltas = []
    for b in range(n_boot):
        seed_str = f"{category},{perturbation},delta,composed,{b}"
        seed = int(hashlib.sha256(seed_str.encode()).hexdigest()[:8], 16)
        rng = np.random.default_rng(seed)
        boot_idx = []
        for cls in LABELS:
            idx = class_indices.get(cls, [])
            if len(idx) > 0:
                boot_idx.extend(rng.choice(idx, size=len(idx), replace=True))
        f1_base = compute_macro_f1([y_true_base[i] for i in boot_idx],
                                    [y_pred_base[i] for i in boot_idx])
        f1_abl = compute_macro_f1([y_true_abl[i] for i in boot_idx],
                                   [y_pred_abl[i] for i in boot_idx])
        deltas.append(f1_abl - f1_base)
    return float(np.percentile(deltas, 2.5)), float(np.percentile(deltas, 97.5))


def run_ablations(cfg, thresholds, thresholds_zsum):
    """Run all ablation experiments with paired raw predictions."""
    print("\n--- Ablations ---")
    all_rows = []
    raw_preds = {}  # name -> (y_true, y_pred)

    # Baseline (no ablation) for paired comparison
    baseline_result = run_grid_cell(
        cfg, thresholds, thresholds_zsum,
        category="ablation", perturbation="baseline",
        severity="0", n_reps=100, return_raw=True,
    )
    baseline_rows, base_true, base_pred = baseline_result
    all_rows.extend(baseline_rows)
    raw_preds["baseline"] = (base_true, base_pred)

    ablations = [
        ("remove_monotone", {"S": float("inf"), "Z_MK": float("inf"),
                            "med_last": 100.0, "med_first": 0.0, "MAD_dx": 0.001}, 10, 1.0),
        ("remove_mann_kendall", {"Z_MK": float("inf")}, 10, 1.0),
        ("remove_curvature", {"R_curve": 0.0}, 10, 1.0),
        ("remove_period_filter", {}, 0, 1.0),
        ("remove_01_test", {"K_01": 1.0}, 10, 1.0),
        ("remove_pe", {"PE": 0.75}, 10, 1.0),
        ("curvature_minus20", {}, 10, 0.8),
        ("curvature_plus20", {}, 10, 1.2),
    ]

    for name, feat_overrides, pmin, cmult in ablations:
        print(f"  Ablation: {name}")
        result = run_grid_cell(
            cfg, thresholds, thresholds_zsum,
            category="ablation", perturbation=name,
            severity="1", n_reps=100,
            period_min=pmin, curve_mult=cmult,
            feature_overrides=feat_overrides if feat_overrides else None,
            return_raw=True,
        )
        rows, abl_true, abl_pred = result
        raw_preds[name] = (abl_true, abl_pred)
        all_rows.extend(rows)
        comp_row = [r for r in rows if r["signal"] == "composed"][0]
        print(f"    F1={comp_row['macro_f1']}  delta={comp_row['paired_delta']}")

    # Store raw_preds on the function for plot_ablation to use
    run_ablations._raw_preds = raw_preds
    return all_rows


# ====================================================================
# DEGRADATIONS
# ====================================================================

def _run_one_degradation(args):
    """Worker for parallel degradation runs."""
    cfg, thresholds, thresholds_zsum, perturbation, severity_str, n_reps, degradation, deg_params = args
    rows = run_grid_cell(
        cfg, thresholds, thresholds_zsum,
        category="degradation", perturbation=perturbation,
        severity=severity_str, n_reps=n_reps,
        degradation=degradation, deg_params=deg_params,
    )
    return perturbation, severity_str, rows


def run_degradations(cfg, thresholds, thresholds_zsum):
    """Run all degradation experiments."""
    print("\n--- Degradations ---")
    all_rows = []

    # Autocorrelated null (AR(1))
    for phi in [0.2, 0.5, 0.8]:
        print(f"  AR(1) phi={phi}")
        rows = run_grid_cell(
            cfg, thresholds, thresholds_zsum,
            category="degradation", perturbation="ar1",
            severity=f"{phi}", n_reps=500,
            degradation="ar1", deg_params={"phi": phi},
        )
        all_rows.extend(rows)
        comp_row = [r for r in rows if r["signal"] == "composed"][0]
        print(f"    F1={comp_row['macro_f1']}  null_fpr={comp_row['null_fpr']}")

    # Correlated streams
    for rho in [0.1, 0.3, 0.6]:
        print(f"  Correlated rho={rho}")
        rows = run_grid_cell(
            cfg, thresholds, thresholds_zsum,
            category="degradation", perturbation="correlated",
            severity=f"{rho}", n_reps=500,
            degradation="correlated", deg_params={"rho": rho},
        )
        all_rows.extend(rows)
        comp_row = [r for r in rows if r["signal"] == "composed"][0]
        print(f"    F1={comp_row['macro_f1']}  delta={comp_row['paired_delta']}")

    # Missing data
    for frac in [0.05, 0.10, 0.25]:
        print(f"  Missing frac={frac}")
        rows = run_grid_cell(
            cfg, thresholds, thresholds_zsum,
            category="degradation", perturbation="missing",
            severity=f"{frac}", n_reps=500,
            degradation="missing", deg_params={"frac": frac},
        )
        all_rows.extend(rows)
        comp_row = [r for r in rows if r["signal"] == "composed"][0]
        print(f"    F1={comp_row['macro_f1']}")

    # Misspecified distribution (Normal -> t)
    for df in [3, 5, 10]:
        print(f"  t-misspec df={df}")
        rows = run_grid_cell(
            cfg, thresholds, thresholds_zsum,
            category="degradation", perturbation="t_misspec",
            severity=f"{df}", n_reps=500,
            degradation="t_misspec", deg_params={"df": df},
        )
        all_rows.extend(rows)
        comp_row = [r for r in rows if r["signal"] == "composed"][0]
        print(f"    F1={comp_row['macro_f1']}")

    # Nonstationary baseline
    for dc in [0.005, 0.01, 0.02]:
        print(f"  Drift coef={dc}")
        rows = run_grid_cell(
            cfg, thresholds, thresholds_zsum,
            category="degradation", perturbation="drift",
            severity=f"{dc}", n_reps=500,
            degradation="drift", deg_params={"drift_coef": dc},
        )
        all_rows.extend(rows)
        comp_row = [r for r in rows if r["signal"] == "composed"][0]
        print(f"    F1={comp_row['macro_f1']}  null_fpr={comp_row['null_fpr']}")

    return all_rows


# ====================================================================
# PERIOD FILTER SWEEP
# ====================================================================

def run_filter_sweep(cfg, thresholds, thresholds_zsum):
    """Sweep period_min from 2 to 50 in steps of 2."""
    print("\n--- Period filter sweep ---")
    all_rows = []

    for pmin in range(2, 52, 2):
        rows = run_grid_cell(
            cfg, thresholds, thresholds_zsum,
            category="filter", perturbation="period_threshold",
            severity=f"{pmin}", n_reps=100,
            period_min=pmin,
        )
        all_rows.extend(rows)
        comp_row = [r for r in rows if r["signal"] == "composed"][0]
        print(f"  period_min={pmin}: F1={comp_row['macro_f1']}")

    return all_rows


# ====================================================================
# MIXED DYNAMICS
# ====================================================================

def run_mixed_dynamics(cfg, thresholds, thresholds_zsum):
    """Run mixed dynamics experiments (descriptive only)."""
    print("\n--- Mixed dynamics ---")
    all_rows = []
    ratios = [(0.2, 0.8), (0.4, 0.6), (0.5, 0.5), (0.6, 0.4), (0.8, 0.2)]
    N = cfg["composition"]["n_observations"]

    mixed_configs = [
        ("trend_oscillation", "divergence", "oscillation"),
        ("oscillation_aperiodic", "oscillation", "aperiodic"),
        ("convergence_oscillation", "convergence", "oscillation"),
    ]

    for name, cond1, cond2 in mixed_configs:
        a1 = BASE_AMPLITUDES[cond1]
        a2 = BASE_AMPLITUDES[cond2]
        for w1, w2 in ratios:
            def _make_mixed_forcing(n, rep, rng, _c1=cond1, _c2=cond2,
                                     _a1=a1, _a2=a2, _w1=w1, _w2=w2):
                f1 = make_forcing(_c1, n, rep)
                f2 = make_forcing(_c2, n, rep)
                return _w1 * _a1 * f1 + _w2 * _a2 * f2

            print(f"  {name} w1={w1} w2={w2}")
            rows = run_grid_cell(
                cfg, thresholds, thresholds_zsum,
                category="mixed", perturbation=name,
                severity=f"{w1}:{w2}", n_reps=100,
                mixed_forcing_fn=_make_mixed_forcing,
                is_mixed=True,
            )
            all_rows.extend(rows)
            comp_row = rows[0]
            print(f"    entropy={comp_row['label_entropy']}  counts={comp_row['label_counts']}")

    # Weak divergence + heavy-tailed noise (df=3)
    a_div = BASE_AMPLITUDES["divergence"]
    for w1, w2 in ratios:
        def _make_div_t(n, rep, rng, _a_div=a_div, _w1=w1):
            f_div = make_forcing("divergence", n, rep)
            t_noise = rng.standard_t(3, n)
            return _w1 * _a_div * f_div + (1 - _w1) * 0.3 * t_noise

        print(f"  div_heavy_tail w1={w1}")
        rows = run_grid_cell(
            cfg, thresholds, thresholds_zsum,
            category="mixed", perturbation="div_heavy_tail",
            severity=f"{w1}:{1-w1}", n_reps=100,
            mixed_forcing_fn=_make_div_t,
            is_mixed=True,
        )
        all_rows.extend(rows)
        comp_row = rows[0]
        print(f"    entropy={comp_row['label_entropy']}  counts={comp_row['label_counts']}")

    return all_rows


# ====================================================================
# E-VALUE VALIDITY
# ====================================================================

def run_evalue_validity(cfg, thresholds):
    """E-value validity check under t(df=3) misspecification for null."""
    print("\n--- E-value validity ---")
    comp = cfg["composition"]
    streams = comp["streams"]
    stream_names = list(streams.keys())
    N = comp["n_observations"]

    n_reps = 1000
    max_Es = []
    # For supermartingale: accumulate mean E_t at each time point
    cum_log_e_all = np.zeros((n_reps, N))

    for rep in range(n_reps):
        # Separate seeds from degradation reps: use 70000 + rep
        seed = 70000 + rep
        rng = np.random.default_rng(seed)
        forcing = np.zeros(N)

        # Step 3: generate ALL streams first (normative RNG order)
        all_log_e = []
        for sname in stream_names:
            _, log_e = generate_stream(streams[sname], forcing, rng)
            all_log_e.append(log_e)
        # Step 5: t-noise replacement LAST
        temp_idx = stream_names.index("temperature")
        df = 3
        sigma = streams["temperature"]["null_params"]["sigma"]
        lam = streams["temperature"]["alt_params"]["lambda"]
        t_noise = rng.standard_t(df, N)
        t_noise *= np.sqrt((df - 2) / df)
        x_new = forcing + t_noise * sigma
        all_log_e[temp_idx] = lam * x_new - lam**2 / 2

        composed = np.sum(all_log_e, axis=0)
        cum_log_e = np.cumsum(composed)
        cum_log_e_all[rep] = cum_log_e
        max_E = np.exp(np.max(cum_log_e))
        max_Es.append(max_E)

        if (rep + 1) % 200 == 0:
            print(f"  validity: {rep + 1}/{n_reps}")

    # Anytime exceedance
    exceedance_rate = np.mean([e > 20 for e in max_Es])

    # Supermartingale: mean of E_t = exp(cumsum) at each time
    # Use log-sum-exp for numerical stability instead of clipping
    # mean(exp(x)) = exp(log(mean(exp(x)))) = exp(logsumexp(x) - log(n))
    from scipy.special import logsumexp
    log_mean_Et = np.array([
        logsumexp(cum_log_e_all[:, t]) - np.log(n_reps)
        for t in range(N)
    ])
    mean_Et = np.exp(log_mean_Et)
    max_mean_Et = float(np.max(mean_Et))
    if not np.isfinite(max_mean_Et):
        max_mean_Et = float("inf")

    print(f"  Anytime exceedance rate (>20): {exceedance_rate:.4f}")
    print(f"  Max mean E_t (supermartingale): {max_mean_Et:.4f}")
    print(f"  Inflated? {'YES' if max_mean_Et > 1.05 else 'no'}")

    return {
        "exceedance_rate": exceedance_rate,
        "max_mean_Et": max_mean_Et,
    }


# ====================================================================
# PLOTTING
# ====================================================================

def plot_sensitivity(sweep_data):
    """4 subplots (one per bin), x=amplitude, y=detection rate."""
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    axes = axes.flatten()

    non_null = [c for c in CONDITIONS if c != "null"]
    for i, cond in enumerate(non_null):
        ax = axes[i]
        data = sweep_data[cond]
        amps = [d[0] for d in data]
        rates = [d[1] for d in data]

        # 95% CI via Wilson interval for proportions
        n_reps = 100
        ci_lo = [max(0, r - 1.96 * np.sqrt(r * (1 - r) / n_reps)) for r in rates]
        ci_hi = [min(1, r + 1.96 * np.sqrt(r * (1 - r) / n_reps)) for r in rates]
        ax.fill_between(amps, ci_lo, ci_hi, alpha=0.2)
        ax.plot(amps, rates, "o-", markersize=4, linewidth=1.5)
        ax.axhline(0.8, color="orange", linestyle="--", alpha=0.7, label="80%")
        ax.axhline(0.9, color="green", linestyle="--", alpha=0.7, label="90%")
        ax.set_xlabel("Amplitude")
        ax.set_ylabel("Detection rate")
        ax.set_title(f"{cond} (base={BASE_AMPLITUDES[cond]:.3f})")
        ax.set_ylim(-0.05, 1.05)
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)

    fig.suptitle("V4 Amplitude Sensitivity", fontsize=14)
    fig.tight_layout()
    fig.savefig(PLOTS / "v4_sensitivity.png", dpi=150)
    plt.close(fig)
    print("  v4_sensitivity.png saved")


def plot_degraded(all_rows):
    """Heatmap of degradation results: rows=type, cols=severity, cell=F1."""
    deg_types = ["ar1", "correlated", "missing", "t_misspec", "drift"]
    severities_map = {
        "ar1": ["0.2", "0.5", "0.8"],
        "correlated": ["0.1", "0.3", "0.6"],
        "missing": ["0.05", "0.1", "0.25"],
        "t_misspec": ["3", "5", "10"],
        "drift": ["0.005", "0.01", "0.02"],
    }

    # Build matrix
    f1_matrix = np.zeros((len(deg_types), 3))
    for i, dt in enumerate(deg_types):
        sevs = severities_map[dt]
        for j, sev in enumerate(sevs):
            match = [r for r in all_rows
                     if r["category"] == "degradation"
                     and r["perturbation"] == dt
                     and r["severity"] == sev
                     and r["signal"] == "composed"]
            if match:
                f1_val = float(match[0]["macro_f1"]) if match[0]["macro_f1"] else 0
                f1_matrix[i, j] = f1_val

    # Custom colormap: red(<0.6) -> yellow(0.8) -> green(>0.9)
    cmap = LinearSegmentedColormap.from_list("ryg", [
        (0.0, "red"), (0.6, "red"), (0.8, "yellow"), (0.9, "green"), (1.0, "green")
    ])

    fig, ax = plt.subplots(figsize=(8, 6))
    im = ax.imshow(f1_matrix, cmap=cmap, vmin=0, vmax=1, aspect="auto")

    # Labels
    sev_labels = []
    for dt in deg_types:
        sev_labels.append(severities_map[dt])

    ax.set_yticks(range(len(deg_types)))
    ax.set_yticklabels([f"{dt}" for dt in deg_types])
    ax.set_xticks(range(3))
    ax.set_xticklabels(["Mild", "Medium", "Harsh"])
    ax.set_xlabel("Severity")
    ax.set_ylabel("Degradation")

    # Annotate
    for i in range(len(deg_types)):
        for j in range(3):
            val = f1_matrix[i, j]
            color = "white" if val < 0.5 else "black"
            sev_str = severities_map[deg_types[i]][j]
            ax.text(j, i, f"{val:.3f}\n({sev_str})",
                    ha="center", va="center", fontsize=9, color=color)

    fig.colorbar(im, ax=ax, label="Composed macro-F1")
    ax.set_title("V4 Degradation Robustness", fontsize=14)
    fig.tight_layout()
    fig.savefig(PLOTS / "v4_degraded.png", dpi=150)
    plt.close(fig)
    print("  v4_degraded.png saved")


def plot_ablation(ablation_rows):
    """Horizontal bar chart: one bar per ablation, x=paired delta F1 vs baseline."""
    # Get baseline F1
    baseline = [r for r in ablation_rows
                if r["perturbation"] == "baseline" and r["signal"] == "composed"]
    if not baseline:
        print("  WARNING: no baseline row for ablation plot")
        return

    baseline_f1 = float(baseline[0]["macro_f1"])

    ablation_names = []
    deltas = []
    ci_los = []
    ci_his = []

    raw_preds = getattr(run_ablations, '_raw_preds', {})
    base_true, base_pred = raw_preds.get("baseline", ([], []))

    for r in ablation_rows:
        if r["signal"] != "composed" or r["perturbation"] == "baseline":
            continue
        name = r["perturbation"]
        f1 = float(r["macro_f1"])
        delta = f1 - baseline_f1
        ablation_names.append(name)
        deltas.append(delta)

        # Paired bootstrap delta CI
        abl_true, abl_pred = raw_preds.get(name, ([], []))
        if base_true and abl_true:
            ci_lo_d, ci_hi_d = paired_bootstrap_delta_ci(
                base_true, base_pred, abl_true, abl_pred, "ablation", name)
        else:
            ci_lo_d, ci_hi_d = delta, delta
        ci_los.append(delta - ci_lo_d)
        ci_his.append(ci_hi_d - delta)

    # Sort by magnitude
    order = np.argsort(np.abs(deltas))[::-1]
    ablation_names = [ablation_names[i] for i in order]
    deltas = [deltas[i] for i in order]
    ci_los = [ci_los[i] for i in order]
    ci_his = [ci_his[i] for i in order]

    fig, ax = plt.subplots(figsize=(10, 6))
    y_pos = range(len(ablation_names))
    colors = ["red" if d < 0 else "green" for d in deltas]
    ax.barh(y_pos, deltas, xerr=[ci_los, ci_his], color=colors, alpha=0.7,
            capsize=3)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(ablation_names, fontsize=9)
    ax.set_xlabel("Delta macro-F1 vs baseline")
    ax.axvline(0, color="black", linewidth=0.8)
    ax.set_title("V4 Ablation Impact", fontsize=14)
    ax.grid(True, axis="x", alpha=0.3)
    fig.tight_layout()
    fig.savefig(PLOTS / "v4_ablation.png", dpi=150)
    plt.close(fig)
    print("  v4_ablation.png saved")


def plot_mixed(mixed_rows):
    """Stacked bar chart per mixed condition, faceted by ratio."""
    # Group by perturbation
    perturbations = []
    seen = set()
    for r in mixed_rows:
        p = r["perturbation"]
        if p not in seen:
            perturbations.append(p)
            seen.add(p)

    n_perturb = len(perturbations)
    fig, axes = plt.subplots(1, n_perturb, figsize=(5 * n_perturb, 6), sharey=True)
    if n_perturb == 1:
        axes = [axes]

    colors_map = {
        "null": "#999999",
        "divergent": "#e74c3c",
        "convergent": "#2ecc71",
        "oscillatory": "#3498db",
        "aperiodic": "#9b59b6",
    }

    for ax_idx, perturb in enumerate(perturbations):
        ax = axes[ax_idx]
        rows = [r for r in mixed_rows if r["perturbation"] == perturb]

        severities = [r["severity"] for r in rows]
        bottoms = np.zeros(len(severities))

        for label in LABELS:
            counts = []
            for r in rows:
                lc = json.loads(r["label_counts"])
                counts.append(lc.get(label, 0))
            ax.bar(range(len(severities)), counts, bottom=bottoms,
                   label=label, color=colors_map.get(label, "#cccccc"))
            bottoms += np.array(counts)

        ax.set_xticks(range(len(severities)))
        ax.set_xticklabels(severities, rotation=45, fontsize=8, ha="right")
        ax.set_title(perturb, fontsize=10)
        ax.set_xlabel("Ratio (w1:w2)")
        if ax_idx == 0:
            ax.set_ylabel("Count")

    # Shared legend
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper right", fontsize=8)
    fig.suptitle("V4 Mixed Dynamics: Label Distribution", fontsize=14)
    fig.tight_layout()
    fig.savefig(PLOTS / "v4_mixed.png", dpi=150)
    plt.close(fig)
    print("  v4_mixed.png saved")


# ====================================================================
# MAIN
# ====================================================================

def main():
    cfg = load_config()
    PLOTS.mkdir(parents=True, exist_ok=True)
    TABLES.mkdir(parents=True, exist_ok=True)

    all_csv_rows = []

    # ── Phase 0: Null calibration ──────────────────────────────────
    print("=" * 60)
    print("PHASE 0: Null calibration")
    print("=" * 60)
    thresholds, _ = run_null_calibration(cfg, n_reps=1000)
    thresholds_zsum = run_null_calibration_zsum(cfg, n_reps=1000)

    # ── Phase 1: Amplitude sweep ───────────────────────────────────
    print("\n" + "=" * 60)
    print("PHASE 1: Amplitude sweep (21 steps x 4 bins x 100 reps)")
    print("=" * 60)
    amp_rows, sweep_data = run_amplitude_sweep(cfg, thresholds, thresholds_zsum)
    all_csv_rows.extend(amp_rows)

    # ── Phase 1b: Equalized difficulty ─────────────────────────────
    print("\n" + "=" * 60)
    print("PHASE 1b: Equalized difficulty")
    print("=" * 60)
    eq_rows = run_equalized_difficulty(cfg, thresholds, thresholds_zsum, sweep_data)
    all_csv_rows.extend(eq_rows)

    # ── Phase 2: Ablations ─────────────────────────────────────────
    print("\n" + "=" * 60)
    print("PHASE 2: Ablations (8 conditions x 100 reps)")
    print("=" * 60)
    abl_rows = run_ablations(cfg, thresholds, thresholds_zsum)
    all_csv_rows.extend(abl_rows)

    # ── Phase 3: Degradations ──────────────────────────────────────
    print("\n" + "=" * 60)
    print("PHASE 3: Degradations (5 types x 3 severities x 500 reps)")
    print("=" * 60)
    deg_rows = run_degradations(cfg, thresholds, thresholds_zsum)
    all_csv_rows.extend(deg_rows)

    # ── Phase 4: Period filter sweep ───────────────────────────────
    print("\n" + "=" * 60)
    print("PHASE 4: Period filter sweep (25 thresholds x 100 reps)")
    print("=" * 60)
    filt_rows = run_filter_sweep(cfg, thresholds, thresholds_zsum)
    all_csv_rows.extend(filt_rows)

    # ── Phase 5: Mixed dynamics ────────────────────────────────────
    print("\n" + "=" * 60)
    print("PHASE 5: Mixed dynamics (4 types x 5 ratios x 100 reps)")
    print("=" * 60)
    mix_rows = run_mixed_dynamics(cfg, thresholds, thresholds_zsum)
    all_csv_rows.extend(mix_rows)

    # ── Phase 6: E-value validity ──────────────────────────────────
    print("\n" + "=" * 60)
    print("PHASE 6: E-value validity (1000 reps)")
    print("=" * 60)
    validity = run_evalue_validity(cfg, thresholds)

    # Validity uses its own semantics — label_counts carries the diagnostics,
    # other fields left empty to avoid schema abuse
    all_csv_rows.append({
        "category": "validity",
        "perturbation": "anytime_exceedance",
        "severity": "df=3",
        "signal": "composed",
        "macro_f1": "",
        "ci_lo": "",
        "ci_hi": "",
        "null_fpr": "",
        "label_entropy": "",
        "label_counts": "",
        "detection_rate": f"exceedance={validity['exceedance_rate']:.4f},max_mean_Et={validity['max_mean_Et']:.4f}",
        "paired_delta": "",
    })
    all_csv_rows.append({
        "category": "validity",
        "perturbation": "supermartingale",
        "severity": "df=3",
        "signal": "composed",
        "macro_f1": "",
        "ci_lo": "",
        "ci_hi": "",
        "null_fpr": "",
        "label_entropy": "",
        "label_counts": "",
        "detection_rate": f"max_mean_Et={validity['max_mean_Et']:.4f},inflated={'yes' if validity['max_mean_Et'] > 1.05 else 'no'}",
        "paired_delta": "",
    })

    # ── Write CSV ──────────────────────────────────────────────────
    print("\n" + "=" * 60)
    print("Writing results")
    print("=" * 60)

    fieldnames = [
        "category", "perturbation", "severity", "signal",
        "macro_f1", "ci_lo", "ci_hi", "null_fpr",
        "label_entropy", "label_counts", "detection_rate", "paired_delta",
    ]
    csv_path = TABLES / "v4_grid.csv"
    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(all_csv_rows)
    print(f"  {csv_path} ({len(all_csv_rows)} rows)")

    # ── Plots ──────────────────────────────────────────────────────
    print("\nGenerating plots...")
    plot_sensitivity(sweep_data)
    plot_degraded(deg_rows)
    plot_ablation(abl_rows)
    plot_mixed(mix_rows)

    # ── Summary ────────────────────────────────────────────────────
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)

    # Report pass/fail for mildest degradations
    mild_degs = [
        ("ar1", "0.2"),
        ("correlated", "0.1"),
        ("missing", "0.05"),
        ("t_misspec", "10"),
        ("drift", "0.005"),
    ]
    for pert, sev in mild_degs:
        match = [r for r in all_csv_rows
                 if r["category"] == "degradation"
                 and r["perturbation"] == pert
                 and r["severity"] == sev
                 and r["signal"] == "composed"]
        if match:
            f1 = float(match[0]["macro_f1"])
            ci_lo = float(match[0]["ci_lo"])
            if ci_lo > 0.8:
                verdict = "PASS"
            elif f1 >= 0.6:
                verdict = "INCONCLUSIVE"
            else:
                verdict = "FAIL"
            print(f"  {pert} ({sev}): F1={f1:.3f} CI=[{ci_lo:.3f}, {float(match[0]['ci_hi']):.3f}] -> {verdict}")

    print(f"\n  E-value exceedance rate: {validity['exceedance_rate']:.4f} (target: <0.05)")
    print(f"  Supermartingale max mean E_t: {validity['max_mean_Et']:.4f} (target: <1.05)")

    print("\nDone.")


if __name__ == "__main__":
    main()
