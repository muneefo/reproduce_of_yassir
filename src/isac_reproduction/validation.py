"""Numerical/Monte Carlo evidence; failures raise rather than being suppressed."""
import json
from pathlib import Path
import time
import numpy as np
from scipy.integrate import quad
from scipy.special import gammaincc, gammaln
from . import analysis, communication, monostatic_audit, simulation, statistics


def alzer_series(k, x, terms):
    eta = np.exp(-gammaln(1+k)/k)
    value = np.zeros_like(np.asarray(x, dtype=float))
    coefficient = k
    for z in range(1, terms+1):
        value += coefficient*np.exp(-z*eta*x)
        coefficient *= -(k-z)/(z+1)
    return value


def independent_guard_mc(cfg, s, samples=8000):
    """Surrogate behind Eq. (11), omitting nearest BS, NOT the exact network."""
    rng = np.random.default_rng(cfg.seed+500)
    radius = 1500.
    result = np.empty((samples, len(s)))
    alignment = cfg.beams*(cfg.beamwidth_rad/(2*np.pi))**2
    for j in range(samples):
        guard_sq = rng.exponential()/(np.pi*cfg.density_m2)
        if guard_sq >= radius**2:
            raise ValueError("Increase guard-MC disk size")
        n = rng.poisson(np.pi*cfg.density_m2*alignment*(radius**2-guard_sq))
        r = np.sqrt(rng.uniform(guard_sq, radius**2, n))
        los = rng.random(n) < np.exp(-cfg.blockage_m_inv*r)
        m = np.where(los, cfg.m_l, cfg.m_n)
        path = np.where(los, cfg.intercept_l*r**(-cfg.alpha_l),
                        cfg.intercept_n*r**(-cfg.alpha_n))
        interference = cfg.communication_w*cfg.gain**2*np.sum(rng.gamma(m, 1/m)*path)
        result[j] = np.exp(-s*interference)
    return statistics.mean_interval(result)


def run(cfg, out):
    start = time.monotonic()
    report = {"status": "component validation; not a reproduction claim", "checks": {}}
    r = np.array([[80., 110., 160., 220., 300., 350.]])
    mean, variance, k, scale = analysis.echo_moments(r, cfg)
    rng = np.random.default_rng(cfg.seed+501)
    n = 200000
    mono = cfg.radar_constant*cfg.rcs_mono_m2/r[0, 0]**(2*cfg.alpha_l)
    bi = cfg.radar_constant*cfg.rcs_bi_m2/(r[0, 0]**cfg.alpha_l*r[0, 1:]**cfg.alpha_l)
    probability = np.exp(-cfg.blockage_m_inv*r[0, 1:])*cfg.beams*cfg.beamwidth_rad/(2*np.pi)
    echoes = rng.exponential(mono, n)+np.sum((rng.random((n, 5)) < probability)*rng.exponential(bi, (n, 5)), axis=1)
    mean_error = abs(echoes.mean()/mean[0]-1)
    variance_error = abs(echoes.var()/variance[0]-1)
    if mean_error > .02 or variance_error > .04:
        raise AssertionError("Eq. (5) empirical moment check failed")
    x = np.linspace(.01, 8, 150)
    empirical_sf = (echoes[:, None] > x[None, :]*scale[0]).mean(axis=0)
    gamma_sf = gammaincc(k[0], x)
    eta = np.exp(-gammaln(1+k[0])/k[0])
    alzer = -np.expm1(k[0]*np.log(-np.expm1(-eta*x)))
    report["checks"]["echo_moments_and_approximations"] = {
        "samples": n, "radii_m": r[0].tolist(), "mean_relative_error": mean_error,
        "variance_relative_error": variance_error, "shape": k[0],
        "gamma_vs_original_max_sf_error": float(np.max(np.abs(gamma_sf-empirical_sf))),
        "alzer_vs_gamma_max_sf_error": float(np.max(np.abs(alzer-gamma_sf)))}
    e1 = []
    for a in (1e-4, 1, 1e4, 1e8):
        exact, error = quad(lambda rr: np.exp(-cfg.blockage_m_inv*rr)*a*rr/(rr*rr+a),
                            80, np.inf, epsabs=1e-10, epsrel=1e-10)
        value = float(analysis.reflected_integral(a, 80, cfg.blockage_m_inv))
        if abs(value-exact) > max(1e-9, abs(exact)*1e-8):
            raise AssertionError("Reflected E1 integral failed quadrature comparison")
        e1.append({"a": a, "E1": value, "quadrature": exact, "quadrature_error": error})
    report["checks"]["reflected_integral"] = e1
    s = np.array([1e4, 1e8, 1e12])
    direct = analysis.DirectLT(cfg)
    fine = analysis.DirectLT(cfg, radial_nodes=4096, outer_nodes=96, s_nodes=4096)
    quad_values = np.array([analysis.direct_lt_quad(ss, cfg, outer_nodes=48) for ss in s])
    if np.max(np.abs(direct(s)-quad_values)) > 2e-4:
        raise AssertionError("Direct LT failed independent quadrature comparison")
    surrogate, surrogate_ci = independent_guard_mc(cfg, s)
    if np.any(np.abs(surrogate-direct(s)) > (surrogate_ci[:, 1]-surrogate_ci[:, 0])+5e-4):
        raise AssertionError("Direct LT disagrees with its surrogate Monte Carlo")
    report["checks"]["direct_lt"] = {"s": s.tolist(), "table": direct(s).tolist(),
        "quad": quad_values.tolist(), "doubled_resolution": fine(s).tolist(),
        "surrogate_mc": surrogate.tolist(), "surrogate_mc_ci": surrogate_ci.tolist(),
        "surrogate_disk_radius_m": 1500}
    # Literal Eq. (13) vs the correct conditional PPP outside the cluster radius.
    s_ref = np.array([1e7, 1e9, 1e11])
    a = s_ref*cfg.radar_constant*cfg.rcs_bi_m2/80**cfg.alpha_l
    true_mono = np.exp(-cfg.density_m2*cfg.beams*cfg.beamwidth_rad*
                       analysis.reflected_integral(a, 80, cfg.blockage_m_inv))
    literal = analysis.reflected_lt(s_ref, np.full(3, 80.), 1, cfg)
    ref_fine = analysis.reflected_lt(s_ref, np.full(3, 80.), 1, cfg, 64)
    report["checks"]["reflected_exclusion_audit"] = {"s": s_ref.tolist(),
        "R1_m": 80, "literal_eq13": literal.tolist(), "eq13_64_nodes": ref_fine.tolist(),
        "exact_conditional_PPP_outside_R1": true_mono.tolist(),
        "note": "Different models: Eq. (13) excludes the first outside point contribution."}
    # Original marked PPP reflected interference at fixed R1.
    ref_samples = np.empty((50000, 3))
    for j in range(len(ref_samples)):
        radii = np.sqrt(rng.uniform(80**2, 1500**2,
                     rng.poisson(np.pi*cfg.density_m2*(1500**2-80**2))))
        gate = rng.random(len(radii)) < np.exp(-cfg.blockage_m_inv*radii)*cfg.beams*cfg.beamwidth_rad/(2*np.pi)
        interference = cfg.radar_constant/80**cfg.alpha_l*np.sum(
            rng.exponential(cfg.rcs_bi_m2, len(radii))*radii**(-cfg.alpha_l)*gate)
        ref_samples[j] = np.exp(-s_ref*interference)
    ref_mean, ref_ci = statistics.mean_interval(ref_samples)
    if np.any(np.abs(ref_mean-true_mono) > (ref_ci[:, 1]-ref_ci[:, 0])+5e-4):
        raise AssertionError("Original reflected PPP LT check failed")
    report["checks"]["reflected_exclusion_audit"].update(
        {"original_ppp_mc": ref_mean.tolist(), "original_ppp_mc_ci": ref_ci.tolist()})
    series = []
    xx = np.array([.01, .1, 1., 10.])
    for shape in (1., 1.15, 1.4, 2., 2.5, 4.):
        eta = np.exp(-gammaln(1+shape)/shape)
        stable = -np.expm1(shape*np.log(-np.expm1(-eta*xx)))
        errors = [float(np.max(np.abs(alzer_series(shape, xx, z)-stable))) for z in (64, 256, 1024)]
        if errors[-1] > 1e-6:
            raise AssertionError("Generalized-binomial series failed convergence")
        series.append({"shape": shape, "max_errors_64_256_1024": errors,
                       "alzer_vs_gamma_max_error": float(np.max(np.abs(stable-gammaincc(shape, xx))))})
    report["checks"]["alzer_noninteger_series"] = series
    rr = analysis.distance_qmc(cfg, 5, 6)
    p1, d1 = analysis.conditional_coverage(rr, [10**.5], cfg, direct, tolerance=2e-4)
    p2, d2 = analysis.conditional_coverage(rr, [10**.5], cfg, direct, tolerance=2e-5, max_terms=4096)
    p3, d3 = analysis.conditional_coverage(rr, [10**.5], cfg, fine, tolerance=2e-4, outer_nodes=64)
    report["checks"]["coverage_numerics"] = {
        "tightened_series_max_difference": float(np.max(np.abs(p1-p2))),
        "doubled_quadrature_max_difference": float(np.max(np.abs(p1-p3))),
        "default": d1, "tight_series": d2, "fine_quadrature": d3}
    # Limits and PPP order statistics.
    zero_direct = cfg.changed(gamma=1)
    if not np.allclose(analysis.DirectLT(zero_direct)(s), 1):
        raise AssertionError("Zero communication power LT limit failed")
    zero = cfg.changed(gamma=0)
    sinr = simulation.simulate(zero, 50, (1, 6))["sinr"]
    if np.any(sinr != 0):
        raise AssertionError("Zero sensing-power limit failed")
    ro = simulation.ordered_radii(rng, 100000, 6, cfg.density_m2)
    area = np.pi*cfg.density_m2*ro**2
    increments = np.diff(np.column_stack((np.zeros(len(area)), area)), axis=1)
    if np.max(np.abs(increments.mean(0)-1)) > .015:
        raise AssertionError("PPP area increment check failed")
    report["checks"]["limits_and_PPP"] = {"area_increment_means": increments.mean(0).tolist(),
        "zero_sensing": "passed", "zero_communication_direct_LT": "passed"}
    uniform_p = (np.arange(4096)+.5)/4096
    ts = np.array([.2, .5, .8])
    inversion = statistics.gil_pelaez(uniform_p, ts, 400, 8001)
    if np.max(np.abs(inversion-(1-ts))) > .003:
        raise AssertionError("Gil-Pelaez inversion failed uniform limiting case")
    report["checks"]["moments_and_inversion"] = {"uniform_M1": statistics.moments(uniform_p, 1).item(),
        "uniform_M2": statistics.moments(uniform_p, 2).item(), "t": ts.tolist(),
        "inversion": inversion.tolist(), "exact": (1-ts).tolist(),
        "distance_M1_equals_spatial_average": float(statistics.moments(p1[:, 0], 1)-p1.mean())}
    # Eq. (17) quadrature versus MC rate at the default setting (same analysis discrepancy persists).
    rate_checks = []
    for count in (1, 6):
        rrate = analysis.distance_qmc(cfg, 5, count)
        rate, tail, d = analysis.sensing_rate(rrate, cfg, t_max=30, points=61, direct=direct)
        fine_rate, fine_tail, _ = analysis.sensing_rate(rrate, cfg, t_max=40, points=161, direct=direct)
        rate_checks.append({"cluster": count, "t30_step_half": rate, "t40_step_quarter": fine_rate,
                            "difference": fine_rate-rate, "coverage_at_t30": tail,
                            "coverage_at_t40": fine_tail, "diagnostic": d})
    report["checks"]["sensing_rate_quadrature"] = rate_checks
    # [3] communication transfer: compare both explicit beam multiplicities.
    comm_checks = []
    r0 = analysis.distance_qmc(cfg, 9, 1)[:, 0]
    for mode in ("same_band", "ref3_single"):
        raw = simulation.communication_simulation(cfg, 20000, cfg.seed+601, mode)
        thresholds = np.array([1., 10., 100.])
        mc = (raw[:, None] > thresholds).mean(0)
        ci = statistics.wilson((raw[:, None] > thresholds).sum(0), len(raw))
        predicted = communication.conditional_coverage(r0, thresholds, cfg, mode).mean(0)
        factor = cfg.beams*(cfg.slot_s-cfg.pulse_s)/cfg.slot_s
        rate_mc, rate_ci = statistics.mean_interval(factor*np.log1p(raw))
        rate_an, tail = communication.rate(r0[:64], cfg, mode, points=81, t_max=40)
        comm_checks.append({"mode": mode, "thresholds": thresholds.tolist(), "mc": mc.tolist(),
                            "mc_ci": ci.tolist(), "analysis": predicted.tolist(),
                            "analysis_minus_mc": (predicted-mc).tolist(), "rate_mc": rate_mc,
                            "rate_mc_ci": rate_ci.tolist(), "rate_analysis_64_QMC": rate_an,
                            "rate_tail_at40": tail,
                            "note": "Alzer approximation; transfer of [3] to shared bandwidth unresolved."})
    report["checks"]["communication"] = comm_checks
    rjoint = analysis.distance_qmc(cfg, 10, 1)[:, 0]
    taus = 10**(np.array([-4., 0., 5., 10.])/10)
    joint = monostatic_audit.conditional_coverage(rjoint, taus, cfg)
    joint_fine = monostatic_audit.conditional_coverage(rjoint, taus, cfg, 48, 192)
    rate_t = np.linspace(0, 40, 161)
    from scipy.integrate import trapezoid
    rate_joint = cfg.beams*trapezoid(monostatic_audit.conditional_coverage(
        rjoint[:256], np.expm1(rate_t), cfg).mean(0), rate_t)
    mono_sim = simulation.simulate(cfg, 20000, (1,), cfg.seed+711)["sinr"]
    mc, mc_ci = statistics.coverage(mono_sim, [-4, 0, 5, 10])
    rate_mc, rate_ci = statistics.mean_interval(cfg.beams*np.log1p(mono_sim))
    if np.max(np.abs(joint-joint_fine)) > 1e-4:
        raise AssertionError("Shared-geometry monostatic quadrature failed refinement")
    if np.any(np.abs(mc[0]-joint.mean(0)) > (mc_ci[0, :, 1]-mc_ci[0, :, 0])):
        raise AssertionError("Joint monostatic PGFL disagrees with simulation")
    report["checks"]["shared_geometry_monostatic_audit"] = {
        "threshold_db": [-4, 0, 5, 10], "joint_PGFL": joint.mean(0).tolist(),
        "mc": mc[0].tolist(), "mc_ci": mc_ci[0].tolist(),
        "quadrature_refinement_max_difference": float(abs(joint-joint_fine).max()),
        "rate_joint": rate_joint, "rate_mc": rate_mc[0], "rate_mc_ci": rate_ci[0].tolist(),
        "note": "Independent derivation, not the paper's Eq. (14)."}
    report["elapsed_seconds"] = time.monotonic()-start
    Path(out/"validation.json").write_text(json.dumps(report, indent=2, default=float)+"\n")
    print(json.dumps(report, indent=2, default=float), flush=True)
