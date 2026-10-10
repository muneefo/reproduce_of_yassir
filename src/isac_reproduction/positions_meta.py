"""Distinguish conditioning on BS positions alone from persistent spatial marks."""
import json
from pathlib import Path
import numpy as np
from .config import Config
from . import plotting, simulation, statistics


def main():
    cfg = Config()
    out = Path("results")
    report = []
    fig, ax = plotting.figure()
    grid = np.linspace(0, 1, 101)
    for count in (1, 2, 4, 6):
        prior = np.load(out/f"meta_nc{count}.npz")
        samples = simulation.conditional_simulation(cfg, prior["conditional_radii"], 128,
                                                    cfg.seed+count+1901, "positions")
        p = (samples > 10**.5).mean(1)
        mean, ci = statistics.mean_interval(p)
        report.append({"cluster": count, "positions_only_mean": float(mean),
                       "outer_group_mean_ci95": ci.tolist(), "fraction_estimated_p_ge_09": float((p >= .9).mean()),
                       "outer": len(p), "inner": samples.shape[1],
                       "note": "Finite-inner-MC estimates; positions fixed, beam/blockage/fading/RCS marks resampled."})
        np.savez_compressed(out/f"positions_meta_nc{count}.npz", sinr=samples, estimated_p=p,
                            conditioned_radii=prior["conditional_radii"])
        if count == 1:
            for name, values, style in (("Ordered cluster distances", prior["estimated_p"], "-"),
                                        ("All BS positions", p, "--"),
                                        ("BS positions + beams + blockages", prior["estimated_full_p"], ":")):
                ax.plot(grid, statistics.meta_ccdf(values, grid), style, label=name)
        print(report[-1], flush=True)
    ax.set(xlabel="Reliability threshold t", ylabel="Estimated meta CCDF",
           title="Monostatic conditioning audit at 5 dB — preliminary nested MC")
    ax.legend(fontsize=8)
    plotting.save(fig, out/"positions_conditioning.png")
    (out/"positions_conditioning.json").write_text(json.dumps(report, indent=2)+"\n")


if __name__ == "__main__":
    main()
