"""Generate synthetic data for all conditions. Reads configs/conditions.yaml."""

import numpy as np
import yaml
import os
from pathlib import Path

ROOT = Path(__file__).parent.parent
CONFIG = ROOT / "configs" / "conditions.yaml"
DATA = ROOT / "data"


def load_config():
    with open(CONFIG) as f:
        return yaml.safe_load(f)


def generate_iid_normal(n, mu, sigma, rng):
    return rng.normal(mu, sigma, n)


def generate_sinusoidal(n, offset, amplitude, period, sigma, rng):
    t = np.arange(n)
    mu_t = offset + amplitude * np.sin(2 * np.pi * t / period)
    return mu_t + rng.normal(0, sigma, n)


def generate_ar1(n, phi, c, sigma, rng):
    x = np.zeros(n)
    x[0] = c / (1 - phi)
    for t in range(1, n):
        x[t] = phi * x[t - 1] + c + rng.normal(0, sigma)
    return x


def _run_lotka_volterra(n_steps, alpha, beta, delta, gamma, x0, y0, dt,
                        process_sigma, rng):
    x = np.zeros(n_steps)
    y = np.zeros(n_steps)
    x[0], y[0] = x0, y0
    for t in range(1, n_steps):
        dx = (alpha * x[t - 1] - beta * x[t - 1] * y[t - 1]) * dt
        dy = (delta * x[t - 1] * y[t - 1] - gamma * y[t - 1]) * dt
        if process_sigma > 0:
            dx += process_sigma * x[t - 1] * rng.normal() * np.sqrt(dt)
            dy += process_sigma * y[t - 1] * rng.normal() * np.sqrt(dt)
        x[t] = max(x[t - 1] + dx, 0.001)
        y[t] = max(y[t - 1] + dy, 0.001)
    return x, y


def generate_lotka_volterra(n, alpha, beta, delta, gamma, x0, y0, dt,
                            obs_sigma, process_sigma, rng,
                            burn_in_mean=None, burn_in_std=None):
    steps_per_obs = max(1, int(1.0 / dt))
    total_steps = n * steps_per_obs
    x_raw, _ = _run_lotka_volterra(total_steps, alpha, beta, delta, gamma,
                                    x0, y0, dt, process_sigma, rng)
    x_sampled = x_raw[::steps_per_obs][:n]

    if burn_in_mean is not None and burn_in_std is not None:
        x_centered = (x_sampled - burn_in_mean) / burn_in_std + 0.3
    else:
        x_centered = (x_sampled - x_sampled.mean()) / x_sampled.std() + 0.3

    return x_centered + rng.normal(0, obs_sigma, n)


def generate_hmm(n, mu_states, transition_prob, sigma, rng):
    states = np.zeros(n, dtype=int)
    obs = np.zeros(n)
    states[0] = 0
    obs[0] = rng.normal(mu_states[states[0]], sigma)
    for t in range(1, n):
        if rng.random() < transition_prob:
            states[t] = 1 - states[t - 1]
        else:
            states[t] = states[t - 1]
        obs[t] = rng.normal(mu_states[states[t]], sigma)
    return obs


def generate_regime_switch(n, mu_before, mu_after, switch_at, sigma, rng):
    obs = np.zeros(n)
    for t in range(n):
        mu = mu_before if t < switch_at else mu_after
        obs[t] = rng.normal(mu, sigma)
    return obs


def compute_burn_in_stats(config_params):
    """Run a deterministic LV to get centering stats."""
    rng = np.random.default_rng(0)
    p = config_params
    steps_per_obs = max(1, int(1.0 / p["dt"]))
    total_steps = 20000 * steps_per_obs
    x_raw, _ = _run_lotka_volterra(
        total_steps, p["alpha"], p["beta"], p["delta"], p["gamma"],
        p["x0"], p["y0"], p["dt"], 0.0, rng
    )
    x_sampled = x_raw[::steps_per_obs]
    return float(x_sampled.mean()), float(x_sampled.std())


GENERATORS = {
    "iid_normal": lambda n, p, rng: generate_iid_normal(n, p["mu"], p["sigma"], rng),
    "sinusoidal": lambda n, p, rng: generate_sinusoidal(
        n, p["offset"], p["amplitude"], p["period"], p["sigma"], rng),
    "ar1": lambda n, p, rng: generate_ar1(n, p["phi"], p["c"], p["sigma"], rng),
    "hmm": lambda n, p, rng: generate_hmm(n, p["mu_states"], p["transition_prob"], p["sigma"], rng),
}


def main():
    cfg = load_config()
    n = cfg["experiment"]["n_observations"]
    n_reps = cfg["experiment"]["n_replications"]
    base_seed = cfg["experiment"]["base_seed"]

    # Compute LV burn-in stats once
    lv_params = cfg["conditions"]["b_lotka_volterra"]["params"]
    burn_mean, burn_std = compute_burn_in_stats(lv_params)
    print(f"LV burn-in: mean={burn_mean:.3f}, std={burn_std:.3f}")

    # Generate Claim 1 & 2 conditions
    for cond_key, cond in cfg["conditions"].items():
        cond_dir = DATA / cond_key
        cond_dir.mkdir(parents=True, exist_ok=True)
        dgp = cond["dgp"]
        params = cond["params"]

        for rep in range(n_reps):
            seed = base_seed + hash(cond_key) % 10000 * 1000 + rep
            rng = np.random.default_rng(seed)

            if dgp == "lotka_volterra":
                x = generate_lotka_volterra(
                    n, **params, rng=rng,
                    burn_in_mean=burn_mean, burn_in_std=burn_std
                )
            else:
                x = GENERATORS[dgp](n, params, rng)

            np.savetxt(cond_dir / f"{rep:04d}.csv", x)

        print(f"  {cond['label']}: {n_reps} reps")

    # Generate Claim 3 regime-switch conditions
    for cond_key, cond in cfg["regime_switch_conditions"].items():
        cond_dir = DATA / cond_key
        cond_dir.mkdir(parents=True, exist_ok=True)

        for rep in range(n_reps):
            seed = base_seed + hash(cond_key) % 10000 * 1000 + rep
            rng = np.random.default_rng(seed)
            x = generate_regime_switch(
                n, cond["mu_before"], cond["mu_after"],
                cond["switch_at"], 1.0, rng
            )
            np.savetxt(cond_dir / f"{rep:04d}.csv", x)

        print(f"  {cond['label']}: {n_reps} reps")

    print("Done.")


if __name__ == "__main__":
    main()
