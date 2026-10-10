"""Shared-geometry Monte Carlo, using exponential echoes, never Gamma matching."""
import numpy as np
from .config import Config


def beam_offset(angle, rotation, beams):
    """Nearest beam in a randomly rotated, equally spaced beam set (a choice)."""
    period = 2 * np.pi / beams
    return (angle - rotation + period/2) % period - period/2


def ordered_radii(rng, n, count, density):
    """Eq. (15): PPP squared-area increments are independent Exp(1)."""
    return np.sqrt(np.cumsum(rng.exponential(size=(n, count)), axis=1)/(np.pi*density))


def geometry(rng, cfg, radii=None):
    """Target at square center. No wrapping; truncation is tested separately."""
    side = np.sqrt(cfg.area_m2)
    points = rng.uniform(-side/2, side/2,
                         (rng.poisson(cfg.density_m2*cfg.area_m2), 2))
    if radii is not None:
        if radii[-1] >= side/2:
            raise ValueError("Conditioned cluster reaches square boundary; increase area")
        points = points[np.linalg.norm(points, axis=1) > radii[-1]]
        phi = rng.uniform(0, 2*np.pi, len(radii))
        inside = radii[:, None]*np.column_stack((np.cos(phi), np.sin(phi)))
        points = np.vstack((inside, points))
    if not len(points):
        raise ValueError("Empty finite PPP; increase area rather than censoring this event")
    return points[np.argsort(np.linalg.norm(points, axis=1))]


def spatial_marks(rng, points, cfg):
    """Shared rotations couple direct illumination and target illumination."""
    bs = points[0]
    others = points[1:]
    r = np.linalg.norm(points, axis=1)
    delta = others-bs
    d = np.linalg.norm(delta, axis=1)
    to_target = np.arctan2(-others[:, 1], -others[:, 0])
    to_bs = np.arctan2(-delta[:, 1], -delta[:, 0])
    rotation = rng.uniform(0, 2*np.pi/cfg.beams, len(others))
    target_hit = np.abs(beam_offset(to_target, rotation, cfg.beams)) <= cfg.beamwidth_rad/2
    tx_hit = np.abs(beam_offset(to_bs, rotation, cfg.beams)) <= cfg.beamwidth_rad/2
    receiver_direction = np.arctan2(-bs[1], -bs[0])
    arrival_direction = np.arctan2(delta[:, 1], delta[:, 0])
    difference = (arrival_direction-receiver_direction+np.pi) % (2*np.pi)-np.pi
    rx_hit = np.abs(difference) <= cfg.beamwidth_rad/2
    target_los = rng.random(len(others)) < np.exp(-cfg.blockage_m_inv*r[1:])
    direct_los = rng.random(len(others)) < np.exp(-cfg.blockage_m_inv*d)
    return r, d, target_hit & target_los, tx_hit & rx_hit, direct_los


def sensing_slot(rng, points, cfg, clusters=(1,), marks=None, components=False):
    """Eqs. (3), (7)–(10). All cluster sizes use the same echo/interference draws."""
    if len(points) < max(clusters):
        raise ValueError("Too few BSs; increase finite area")
    if marks is None:
        marks = spatial_marks(rng, points, cfg)
    r, d, target_gate, direct_gate, los = marks
    mono = cfg.radar_constant*rng.exponential(cfg.rcs_mono_m2)/r[0]**(2*cfg.alpha_l)
    reflected = (cfg.radar_constant/r[0]**cfg.alpha_l *
                 rng.exponential(cfg.rcs_bi_m2, len(d))*r[1:]**(-cfg.alpha_l)*target_gate)
    shape = np.where(los, cfg.m_l, cfg.m_n)
    fading = rng.gamma(shape, 1/shape)
    path = np.where(los, cfg.intercept_l*d**(-cfg.alpha_l),
                    cfg.intercept_n*d**(-cfg.alpha_n))
    direct = cfg.communication_w*cfg.gain**2*np.sum(direct_gate*fading*path)
    combined = np.concatenate(([0.], np.cumsum(reflected)))
    desired = mono+combined[np.asarray(clusters)-1]
    # Numerical cancellation is avoidable: sum the excluded tail directly.
    reflected_interference = np.array([reflected[n-1:].sum() for n in clusters])
    denominator = reflected_interference+direct+cfg.noise_w+cfg.communication_w*cfg.residual_si
    sinr = desired/denominator
    if not np.all(np.isfinite(sinr)) or np.any(sinr < 0):
        raise FloatingPointError("Invalid sensing SINR")
    if components:
        return sinr, direct, reflected_interference, desired
    return sinr


def simulate(cfg: Config, n=1000, clusters=(1, 2, 4, 6), seed=None):
    rng = np.random.default_rng(cfg.seed if seed is None else seed)
    samples = np.empty((n, len(clusters)))
    direct, reflected = np.empty(n), np.empty((n, len(clusters)))
    for j in range(n):
        samples[j], direct[j], reflected[j], _ = sensing_slot(
            rng, geometry(rng, cfg), cfg, clusters, components=True)
    return {"sinr": samples, "direct_w": direct, "reflected_w": reflected,
            "clusters": np.asarray(clusters)}


def conditional_simulation(cfg, radii, slots=128, seed=None, mode="distances"):
    """Distances mode resamples geometry; full mode fixes positions/beams/blockage.

    The latter is a fading/RCS-only success probability, a different conditioning.
    """
    if mode not in ("distances", "positions", "full"):
        raise ValueError(mode)
    rng = np.random.default_rng(cfg.seed if seed is None else seed)
    result = np.empty((len(radii), slots))
    count = radii.shape[1]
    for j, distances in enumerate(radii):
        if mode in ("positions", "full"):
            points = geometry(rng, cfg, distances)
            marks = spatial_marks(rng, points, cfg) if mode == "full" else None
        for k in range(slots):
            if mode == "distances":
                points = geometry(rng, cfg, distances)
                marks = None
            result[j, k] = sensing_slot(rng, points, cfg, (count,), marks)[0]
    return result


def communication_simulation(cfg, n=1000, seed=None, mode="same_band"):
    """Eq. (21), with explicit same-band vs [3] frequency-separated ambiguity."""
    if mode not in ("same_band", "ref3_single"):
        raise ValueError(mode)
    rng = np.random.default_rng(cfg.seed+17 if seed is None else seed)
    result = np.empty(n)
    q = np.pi/cfg.beamwidth_rad  # Derived from half-power points, not stated q.
    for j in range(n):
        points = geometry(rng, cfg)
        r = np.linalg.norm(points, axis=1)
        extent = np.pi/cfg.beams if mode == "same_band" else np.pi
        theta = rng.uniform(-extent, extent, len(r)-1)
        gain = np.where(np.abs(theta) <= np.pi/q,
                        cfg.gain*np.cos(q*theta/2)**2, 0)
        los = rng.random(len(gain)) < np.exp(-cfg.blockage_m_inv*r[1:])
        m = np.where(los, cfg.m_l, cfg.m_n)
        fading = rng.gamma(m, 1/m)
        path = np.where(los, cfg.intercept_l*r[1:]**(-cfg.alpha_l),
                        cfg.intercept_n*r[1:]**(-cfg.alpha_n))
        desired = (cfg.communication_w*cfg.gain*cfg.intercept_l*r[0]**(-cfg.alpha_l)*
                   rng.gamma(cfg.m_l, 1/cfg.m_l))
        result[j] = desired/(cfg.communication_w*np.sum(gain*fading*path)+cfg.noise_w)
    return result
