"""Compute e-values for all generated data. Reads configs/conditions.yaml."""

import numpy as np
import yaml
from pathlib import Path

ROOT = Path(__file__).parent.parent
CONFIG = ROOT / "configs" / "conditions.yaml"
DATA = ROOT / "data"


def compute_evalues(x, lam):
    log_e = lam * x - lam ** 2 / 2
    cum_log_e = np.cumsum(log_e)
    return log_e, cum_log_e


def main():
    with open(CONFIG) as f:
        cfg = yaml.safe_load(f)

    lam = cfg["evalue"]["lambda"]

    all_conditions = list(cfg["conditions"].keys()) + list(cfg["regime_switch_conditions"].keys())

    for cond_key in all_conditions:
        cond_dir = DATA / cond_key
        if not cond_dir.exists():
            print(f"  skip {cond_key} (no data)")
            continue

        for csv_path in sorted(cond_dir.glob("*.csv")):
            if csv_path.stem.endswith("_evalues"):
                continue
            x = np.loadtxt(csv_path)
            log_e, cum_log_e = compute_evalues(x, lam)

            out = np.column_stack([log_e, cum_log_e])
            np.savetxt(
                csv_path.parent / f"{csv_path.stem}_evalues.csv",
                out, header="log_e_t,cum_log_E_t", comments=""
            )

        print(f"  {cond_key}: e-values computed")

    print("Done.")


if __name__ == "__main__":
    main()
