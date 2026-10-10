"""Focused follow-up checks and plots from saved results; no model fitting."""
import json
from pathlib import Path
import numpy as np
from .config import Config
from . import analysis, plotting, statistics


def main():
    out = Path("results")
    cfg = Config()
    fig, ax = plotting.figure()
    noisy = {}
    for count in (1, 2, 4, 6):
        d = np.load(out/f"meta_nc{count}.npz")
        t = d["inversion_t"]
        if count == 1:
            ax.plot(t, statistics.meta_ccdf(d["analytical_p"], t), label="Same-sample empirical CCDF")
            ax.plot(t, d["inversion_200"], "--", label="Gil–Pelaez u_max=200")
            ax.plot(t, d["inversion_400"], ":", label="Gil–Pelaez u_max=400")
        grid = np.linspace(0, 1, 101)
        indicators = d["estimated_p"][:, None] >= grid
        noisy[str(count)] = {"t": grid.tolist(),
            "noisy_estimate_ccdf": indicators.mean(0).tolist(),
            "noisy_estimate_ccdf_wilson_ci": statistics.wilson(indicators.sum(0), len(indicators)).tolist(),
            "note": "Interval describes the distribution of finite-inner-MC estimates, not exact PS(R)."}
    ax.set(xlabel="t", ylabel="Raw meta CCDF", title="Monostatic inversion diagnostic — no clipping")
    ax.legend()
    plotting.save(fig, out/"inversion_check.png")
    (out/"meta_uncertainty.json").write_text(json.dumps(noisy, indent=2)+"\n")
    # Independent higher-count reflected-only check, motivated by a 95% interval miss.
    s = np.array([1e7, 1e9, 1e11])
    rng = np.random.default_rng(cfg.seed+1801)
    samples = np.empty((200000, len(s)))
    for j in range(len(samples)):
        r = np.sqrt(rng.uniform(80**2, 1500**2,
                    rng.poisson(np.pi*cfg.density_m2*(1500**2-80**2))))
        gate = rng.random(len(r)) < np.exp(-cfg.blockage_m_inv*r)*cfg.beams*cfg.beamwidth_rad/(2*np.pi)
        interference = cfg.radar_constant/80**cfg.alpha_l*np.sum(
            rng.exponential(cfg.rcs_bi_m2, len(r))*r**(-cfg.alpha_l)*gate)
        samples[j] = np.exp(-s*interference)
    mean, ci = statistics.mean_interval(samples)
    a = s*cfg.radar_constant*cfg.rcs_bi_m2/80**cfg.alpha_l
    exact = np.exp(-cfg.density_m2*cfg.beams*cfg.beamwidth_rad*analysis.reflected_integral(
        a, 80, cfg.blockage_m_inv))
    report = {"status": "independent follow-up", "seed": cfg.seed+1801, "samples": len(samples),
              "R1_m": 80, "radius_m": 1500, "s": s.tolist(), "mean": mean.tolist(),
              "ci95": ci.tolist(), "exact_original_PPP": exact.tolist(),
              "standardized_error": ((mean-exact)/samples.std(0, ddof=1)*np.sqrt(len(samples))).tolist()}
    (out/"reflected_followup.json").write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
