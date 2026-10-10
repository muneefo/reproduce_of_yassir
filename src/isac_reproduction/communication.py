"""[3] §IV, Lemma 7 / Theorem 4, with explicit shared-band adaptation.

No misalignment is included: the primary paper Eq. (21) uses G_m.
"""
import math
import numpy as np
from scipy.integrate import cumulative_trapezoid, trapezoid
from scipy.special import gammaln, roots_legendre


def interference_lt(s, serving_distance, cfg, mode="same_band", radial_nodes=1536, angle_nodes=24):
    if mode not in ("same_band", "ref3_single"):
        raise ValueError(mode)
    s, serving_distance = np.broadcast_arrays(s, serving_distance)
    shape = s.shape
    s, serving_distance = s.ravel(), serving_distance.ravel()
    r = np.geomspace(1e-4, 1e12, radial_nodes)
    nodes, weights = roots_legendre(angle_nodes)
    q = np.pi/cfg.beamwidth_rad
    theta = nodes*np.pi/q
    weights = weights*np.pi/q
    gains = cfg.gain*np.cos(q*theta/2)**2
    p = np.exp(-cfg.blockage_m_inv*r)
    multiplicity = cfg.beams if mode == "same_band" else 1
    output = np.empty(len(s))
    for start in range(0, len(s), 32):
        ss = s[start:start+32, None, None]
        l = ss*cfg.communication_w*gains[None, :, None]*cfg.intercept_l*r**(-cfg.alpha_l)/cfg.m_l
        n = ss*cfg.communication_w*gains[None, :, None]*cfg.intercept_n*r**(-cfg.alpha_n)/cfg.m_n
        f = (p*(-np.expm1(-cfg.m_l*np.log1p(l)))+
             (1-p)*(-np.expm1(-cfg.m_n*np.log1p(n))))*r
        f = np.sum(f*weights[None, :, None], axis=1)
        reverse = -cumulative_trapezoid(f[:, ::-1], r[::-1], axis=1, initial=0)[:, ::-1]
        tail = (s[start:start+32]*cfg.communication_w*np.dot(weights, gains)*cfg.intercept_n*
                r[-1]**(2-cfg.alpha_n)/(cfg.alpha_n-2))
        integral = np.array([np.interp(np.log(ro), np.log(r), row)
                             for ro, row in zip(serving_distance[start:start+32], reverse)])+tail
        output[start:start+32] = np.exp(-cfg.density_m2*multiplicity*integral)
    return output.reshape(shape)


def conditional_coverage(r0, thresholds, cfg, mode="same_band", **kwargs):
    r0 = np.atleast_1d(r0)
    thresholds = np.atleast_1d(thresholds)
    if cfg.communication_w == 0:
        return np.zeros((len(r0), len(thresholds)))
    m = int(cfg.m_l)
    if m != cfg.m_l:
        raise NotImplementedError("[3] finite Alzer sum requires integer m_L")
    k = m*np.exp(-gammaln(1+m)/m)
    base = k*thresholds[None, :]*r0[:, None]**cfg.alpha_l/(cfg.communication_w*cfg.gain*cfg.intercept_l)
    result = np.zeros(base.shape)
    for z in range(1, m+1):
        s = z*base
        result += (-1)**(z+1)*math.comb(m, z)*interference_lt(
            s, np.broadcast_to(r0[:, None], s.shape), cfg, mode, **kwargs)*np.exp(-s*cfg.noise_w)
    if np.any(result < 0) or np.any(result > 1+1e-12) or np.any(~np.isfinite(result)):
        raise FloatingPointError("Invalid communication coverage; no clipping")
    return result


def rate(r0, cfg, mode="same_band", t_max=30, points=121, **kwargs):
    t = np.linspace(0, t_max, points)
    probability = conditional_coverage(r0, np.expm1(t), cfg, mode, **kwargs).mean(axis=0)
    factor = cfg.beams*(cfg.slot_s-cfg.pulse_s)/cfg.slot_s
    return factor*trapezoid(probability, t), probability[-1]
