"""Reproducible commands. Run python -m src.isac_reproduction.cli --help."""
import argparse
from dataclasses import fields
import hashlib
import json
from pathlib import Path
import time
import numpy as np
from .config import Config
from . import analysis, communication, monostatic_audit, plotting, simulation, statistics


def write_json(path, data):
    def serialize(v):
        if isinstance(v, Path):
            return str(v)
        if isinstance(v, (np.ndarray, np.generic)):
            return v.tolist()
        raise TypeError(f"Unsupported JSON type: {type(v).__name__}")
    path.write_text(json.dumps(data, indent=2, default=serialize, allow_nan=False)+"\n")


def provenance(cfg, args, seconds):
    return {"status": "preliminary", "config": cfg.to_dict(), "command": vars(args),
            "elapsed_seconds": seconds,
            "source_sha256": {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                              for p in sorted(Path(__file__).parent.glob("*.py"))}}


def run_coverage(cfg, args, out):
    start = time.monotonic()
    clusters = (1,) if args.command == "baseline" else (1, 2, 4, 6)
    db = np.array([-4, -2, 0, 2, 4, 5, 6, 8, 10.])
    result = simulation.simulate(cfg, args.samples, clusters)
    cov, ci = statistics.coverage(result["sinr"], db)
    print("Simulation coverage at 5 dB:", cov[:, db == 5].ravel(), flush=True)
    radii = analysis.distance_qmc(cfg, args.qmc_power, max(clusters))
    lt = analysis.DirectLT(cfg)
    probabilities, diagnostics = [], []
    for count in clusters:
        p, d = analysis.conditional_coverage(radii[:, :count], 10**(db/10), cfg, lt)
        probabilities.append(p)
        diagnostics.append(d)
        print(f"Printed analysis Nc={count}: {p.mean(0)}; {d}", flush=True)
    ana = np.stack([p.mean(0) for p in probabilities])
    joint_mono = monostatic_audit.conditional_coverage(radii[:, 0], 10**(db/10), cfg).mean(0)
    stem = "baseline" if len(clusters) == 1 else "figure3a"
    np.savez_compressed(out/f"{stem}.npz", **result, threshold_db=db, coverage=cov,
                        coverage_ci=ci, analytical=ana, ordered_radii=radii,
                        conditional_probabilities=np.stack(probabilities), joint_monostatic=joint_mono)
    rate, rate_ci = statistics.mean_interval(cfg.beams*np.log1p(result["sinr"]))
    report = provenance(cfg, args, time.monotonic()-start)
    report.update({"clusters": clusters, "threshold_db": db, "coverage": cov,
                   "coverage_ci": ci, "printed_analysis": ana,
                   "analysis_minus_simulation": ana-cov, "analysis_diagnostics": diagnostics,
                   "joint_PGFL_monostatic_audit": joint_mono,
                   "sensing_rate": rate, "sensing_rate_ci": rate_ci,
                   "interference_pearson_correlation": [float(np.corrcoef(result["direct_w"],
                       result["reflected_w"][:, j])[0, 1]) if result["direct_w"].std() > 0
                       and result["reflected_w"][:, j].std() > 0 else None
                       for j in range(len(clusters))]})
    write_json(out/f"{stem}.json", report)
    fig, ax = plotting.figure()
    for j, count in enumerate(clusters):
        label = "Monostatic" if count == 1 else f"Nc={count}"
        color = plotting.COLORS[j]
        ax.plot(db, ana[j], color=color, label=label+" — printed analysis")
        ax.errorbar(db, cov[j], yerr=np.stack((cov[j]-ci[j, :, 0], ci[j, :, 1]-cov[j])),
                    color=color, fmt="o", markersize=3, capsize=2)
    ax.set(xlabel="Sensing SINR threshold (dB)", ylabel="Coverage probability",
           title=f"Figure 3(a) study — preliminary ({args.samples:,} geometries)")
    ax.legend(fontsize=8)
    ax.plot(db, joint_mono, color="#777777", linestyle=":", label="Monostatic shared-geometry audit")
    ax.legend(fontsize=8)
    plotting.save(fig, out/f"{stem}.png")


def run_meta(cfg, args, out):
    start = time.monotonic()
    grid = np.linspace(0, 1, 101)
    inversion_t = np.linspace(.02, .98, 49)
    all_p, all_estimates, all_full, report = [], [], [], []
    lt = analysis.DirectLT(cfg)
    for count in (1, 2, 4, 6):
        radii = analysis.distance_qmc(cfg, args.qmc_power, count, cfg.seed+count)
        p, d = analysis.conditional_coverage(radii, [10**.5], cfg, lt)
        p = p[:, 0]
        conditional_r = simulation.ordered_radii(np.random.default_rng(cfg.seed+count+100),
                                                 args.outer, count, cfg.density_m2)
        samples = simulation.conditional_simulation(cfg, conditional_r, args.slots,
                                                     cfg.seed+count+200, "distances")
        full = simulation.conditional_simulation(cfg, conditional_r, args.slots,
                                                  cfg.seed+count+300, "full")
        estimated = (samples > 10**.5).mean(axis=1)
        estimated_full = (full > 10**.5).mean(axis=1)
        inversion = statistics.gil_pelaez(p, inversion_t, 200, 4001)
        inversion_fine = statistics.gil_pelaez(p, inversion_t, 400, 8001)
        empirical = statistics.meta_ccdf(p, inversion_t)
        average, average_ci = statistics.mean_interval(estimated)
        # Equality is an algebraic identity for one nested sample, not independent validation.
        identity_error = float(abs(estimated.mean()-(samples > 10**.5).mean()))
        summary = {"cluster": count, "printed_analysis_first_moment": p.mean(),
                   "printed_analysis_second_moment": (p*p).mean(),
                   "distance_mc_first_moment": average, "distance_mc_first_moment_ci": average_ci,
                   "full_geometry_mc_first_moment": statistics.mean_interval(estimated_full),
                   "nested_pool_identity_error": identity_error,
                   "nested_inner_max_95pct_halfwidth_bound": 1.96/(2*np.sqrt(args.slots)),
                   "inversion_200_vs_400_max_error": np.max(np.abs(inversion-inversion_fine)),
                   "inversion_400_vs_empirical_max_error": np.max(np.abs(inversion_fine-empirical)),
                   "inversion_raw_min_max": [inversion_fine.min(), inversion_fine.max()],
                   "analysis_diagnostics": d}
        report.append(summary)
        all_p.append(p); all_estimates.append(estimated); all_full.append(estimated_full)
        np.savez_compressed(out/f"meta_nc{count}.npz", analytical_p=p, radii=radii,
                            conditional_radii=conditional_r, distance_sinr=samples, full_sinr=full,
                            estimated_p=estimated, estimated_full_p=estimated_full,
                            inversion_t=inversion_t, inversion_200=inversion, inversion_400=inversion_fine)
        print("Meta-distribution:", summary, flush=True)
    fig, ax = plotting.figure()
    for j, count in enumerate((1, 2, 4, 6)):
        ax.plot(grid, statistics.meta_ccdf(all_p[j], grid), color=plotting.COLORS[j],
                label=f"Nc={count}, printed distance analysis")
        ax.plot(grid, statistics.meta_ccdf(all_estimates[j], grid), color=plotting.COLORS[j],
                linestyle="--", alpha=.6)
    ax.set(xlabel="Reliability threshold t", ylabel="Fraction with conditional success ≥ t",
           title="Figure 3(b) study — 5 dB, distance conditioning (preliminary)")
    ax.legend(fontsize=8)
    plotting.save(fig, out/"figure3b.png")
    fig, ax = plotting.figure()
    for j, count in enumerate((1, 2, 4, 6)):
        ax.plot(grid, statistics.meta_ccdf(all_estimates[j], grid), color=plotting.COLORS[j],
                label=f"Nc={count}, distances")
        ax.plot(grid, statistics.meta_ccdf(all_full[j], grid), color=plotting.COLORS[j], linestyle=":")
    ax.set(xlabel="Reliability threshold t", ylabel="Estimated meta CCDF",
           title="Conditioning comparison — solid: distances; dotted: full marked geometry")
    ax.legend(fontsize=8)
    plotting.save(fig, out/"conditioning_comparison.png")
    manifest = provenance(cfg, args, time.monotonic()-start)
    manifest["checks"] = report
    manifest["warning"] = "Nested MC tail estimates include binomial noise; not original Figure 3 data."
    write_json(out/"figure3b.json", manifest)


def run_rates(cfg, args, out):
    start = time.monotonic()
    clusters = (1, 2, 4, 6, 8)
    density = np.array([10, 30, 50, 70, 90, 110, 120.])
    gamma = np.array([0, .1, .3, .5, .7, .9, .99, 1.])
    summary = {}
    for name, grid in (("density", density), ("gamma", gamma)):
        sensing, sensing_ci, comm, comm_ci, ref3, ref3_ci, only = [], [], [], [], [], [], []
        for value in grid:
            c = cfg.changed(density_m2=value/1e6) if name == "density" else cfg.changed(gamma=value)
            r = simulation.simulate(c, args.samples, clusters)
            mean, ci = statistics.mean_interval(c.beams*np.log1p(r["sinr"]))
            sensing.append(mean); sensing_ci.append(ci)
            cs = simulation.communication_simulation(c, args.samples)
            factor = c.beams*(c.slot_s-c.pulse_s)/c.slot_s
            cm, cci = statistics.mean_interval(factor*np.log1p(cs))
            comm.append(cm); comm_ci.append(cci)
            alternate = simulation.communication_simulation(c, args.samples, mode="ref3_single")
            am, aci = statistics.mean_interval(factor*np.log1p(alternate))
            ref3.append(am); ref3_ci.append(aci)
            # Pure communication: Pc=1 W, no pulse duration penalty.
            pure = c.changed(gamma=0, average_power_w=c.average_power_w*(c.slot_s-c.pulse_s)/c.slot_s)
            pure_samples = simulation.communication_simulation(pure, args.samples)
            only.append(c.beams*np.log1p(pure_samples).mean())
            print(name, value, "sensing", mean, "communication SDMA", cm,
                  "communication ref3 single-frequency", am, flush=True)
        result = {"grid": grid, "sensing": np.asarray(sensing), "sensing_ci": np.asarray(sensing_ci),
                  "communication": np.asarray(comm), "communication_ci": np.asarray(comm_ci),
                  "ref3_single": np.asarray(ref3), "ref3_single_ci": np.asarray(ref3_ci),
                  "communication_only": np.asarray(only), "clusters": np.asarray(clusters)}
        summary[name] = result
        np.savez_compressed(out/f"rates_{name}.npz", **result)
        fig, ax = plotting.figure()
        if name == "density":
            for j, count in enumerate(clusters):
                ax.plot(grid, result["sensing"][:, j], color=plotting.COLORS[j], label=f"Nc={count} sensing")
                ax.fill_between(grid, result["sensing_ci"][:, j, 0], result["sensing_ci"][:, j, 1],
                                color=plotting.COLORS[j], alpha=.12)
            ax.plot(grid, comm, "k--", label="Communication: shared bandwidth")
            ax.plot(grid, ref3, "k:", label="Communication: [3] frequency-separated")
            xlabel, stem = "BS density (BS/km²)", "figure3c"
        else:
            ax.plot(grid, result["sensing"][:, 0], label="Monostatic sensing")
            ax.plot(grid, result["sensing"][:, 3], label="Nc=6 sensing")
            ax.plot(grid, comm, label="Communication: shared bandwidth")
            ax.plot(grid, result["sensing"][:, 0]+comm, label="Monostatic + communication")
            ax.plot(grid, result["sensing"][:, 3]+comm, label="Nc=6 + communication")
            ax.plot(grid, only, "k--", label="Communication-only")
            ax.plot(grid, ref3, "k:", label="Communication: [3] frequency-separated")
            xlabel, stem = "Sensing energy fraction γ", "figure3d"
        ax.set(xlabel=xlabel, ylabel="Rate (nats/s/Hz per BS)",
               title=f"Figure 3({'c' if name == 'density' else 'd'}) study — preliminary; communication interpretation unresolved")
        ax.legend(fontsize=7)
        plotting.save(fig, out/f"{stem}.png")
    report = provenance(cfg, args, time.monotonic()-start)
    report["results"] = summary
    report["note"] = "Rates estimated directly from log(1+SINR), Eqs. (17), (22). No digitization or numerical agreement claim."
    write_json(out/"rates.json", report)


def run_convergence(cfg, args, out):
    start = time.monotonic()
    reports = []
    for area in (6.25e6, 25e6, 100e6):
        c = cfg.changed(area_m2=area, seed=cfg.seed+round(area/1e6))
        result = simulation.simulate(c, args.samples, (1, 2, 4, 6))
        cov, ci = statistics.coverage(result["sinr"], [5])
        rates, rci = statistics.mean_interval(c.beams*np.log1p(result["sinr"]))
        prefix_cov, prefix_ci = statistics.coverage(result["sinr"][:args.samples//4], [5])
        reports.append({"area_km2": area/1e6, "samples": args.samples,
                        "coverage": cov[:, 0], "coverage_ci": ci[:, 0], "rates": rates,
                        "rate_ci": rci, "quarter_sample_coverage": prefix_cov[:, 0],
                        "quarter_sample_ci": prefix_ci[:, 0]})
        print("Area sensitivity:", reports[-1], flush=True)
        np.savez_compressed(out/f"area_{area/1e6:g}.npz", **result)
    lt = analysis.DirectLT(cfg)
    spatial = []
    for exponent in (6, 8, 10):
        for count in (1, 6):
            r = analysis.distance_qmc(cfg, exponent, count)
            p, d = analysis.conditional_coverage(r, [10**.5], cfg, lt)
            spatial.append({"qmc_power": exponent, "cluster": count,
                            "coverage": p.mean(), "diagnostics": d})
            print("Spatial quadrature:", spatial[-1], flush=True)
    manifest = provenance(cfg, args, time.monotonic()-start)
    manifest.update({"area_and_count": reports, "spatial_quadrature": spatial})
    write_json(out/"convergence.json", manifest)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=["baseline", "coverage", "meta", "rates", "convergence", "validate"])
    parser.add_argument("--samples", type=int, default=10000)
    parser.add_argument("--qmc-power", type=int, default=8, help="2**power analytical distance samples")
    parser.add_argument("--outer", type=int, default=128, help="Independent nested-MC distance vectors")
    parser.add_argument("--slots", type=int, default=128, help="Inner resamples per distance vector")
    parser.add_argument("--config", type=Path, help="JSON overrides of Config SI-valued fields")
    parser.add_argument("--output", type=Path, default=Path("results"))
    args = parser.parse_args()
    if args.samples < 4 or args.outer < 2 or args.slots < 2:
        parser.error("At least 4 samples and 2 nested samples per level required")
    if args.config:
        data = json.loads(args.config.read_text())
        unknown = set(data)-{f.name for f in fields(Config)}
        if unknown:
            raise ValueError(f"Unknown config fields: {unknown}")
        cfg = Config(**data)
    else:
        cfg = Config()
    args.output.mkdir(parents=True, exist_ok=True)
    if args.command in ("baseline", "coverage"):
        run_coverage(cfg, args, args.output)
    elif args.command == "meta":
        run_meta(cfg, args, args.output)
    elif args.command == "rates":
        run_rates(cfg, args, args.output)
    elif args.command == "convergence":
        run_convergence(cfg, args, args.output)
    else:
        from .validation import run
        run(cfg, args.output)


if __name__ == "__main__":
    main()
