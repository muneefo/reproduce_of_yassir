"""Independent exact-PPP monostatic audit under our stated beam/blockage choices.

This is NOT Eq. (14). It retains target-centered exclusion and shared-beam
direct/reflected correlation, using an exponential desired echo exactly.
"""
import numpy as np
from scipy.special import roots_legendre
from .analysis import reflected_integral


def conditional_coverage(r1, thresholds, cfg, angle_nodes=24, radial_nodes=96):
    if cfg.alpha_l != 2:
        raise NotImplementedError("Reflected E1 reduction assumes alpha_L=2")
    r1 = np.atleast_1d(r1)
    tau = np.atleast_1d(thresholds)
    if cfg.sensing_w == 0:
        return np.zeros((len(r1), len(tau)))
    mono_mean = cfg.radar_constant*cfg.rcs_mono_m2/r1**(2*cfg.alpha_l)
    s = tau[None, :]/mono_mean[:, None]
    reflected_a = s*cfg.radar_constant*cfg.rcs_bi_m2/r1[:, None]**cfg.alpha_l
    reflected_exponent = cfg.density_m2*cfg.beams*cfg.beamwidth_rad*reflected_integral(
        reflected_a, r1[:, None], cfg.blockage_m_inv)
    an, aw = roots_legendre(angle_nodes)
    rn, rw = roots_legendre(radial_nodes)
    phi = an*cfg.beamwidth_rad/2
    aw = aw*cfg.beamwidth_rad/2
    u = (rn+1)/2
    rw = rw/2
    period = 2*np.pi/cfg.beams
    pbeam = cfg.beamwidth_rad/period
    result = np.empty(s.shape)
    for j, distance in enumerate(r1):
        # Ray exits the target's empty disk at 2 R1 cos(phi).
        lower = 2*distance*np.cos(phi)
        scale = max(distance, 1/cfg.blockage_m_inv)
        d = lower[:, None]+scale*u[None, :]/(1-u[None, :])
        jacobian = scale/(1-u)**2
        x = distance-d*np.cos(phi)[:, None]
        y = -d*np.sin(phi)[:, None]
        target_distance = np.hypot(x, y)
        to_target = np.arctan2(-y, -x)
        to_receiver = np.arctan2(-y, distance-x)
        difference = (to_target-to_receiver+period/2) % period-period/2
        overlap = np.maximum(cfg.beamwidth_rad-np.abs(difference), 0)/period
        # max here is the geometric intersection length, not probability clipping.
        los = np.exp(-cfg.blockage_m_inv*d)
        ss = s[j, :, None, None]
        xl = ss*cfg.communication_w*cfg.gain**2*cfg.intercept_l*d**(-cfg.alpha_l)/cfg.m_l
        xn = ss*cfg.communication_w*cfg.gain**2*cfg.intercept_n*d**(-cfg.alpha_n)/cfg.m_n
        fd = los*(-np.expm1(-cfg.m_l*np.log1p(xl)))+(1-los)*(-np.expm1(-cfg.m_n*np.log1p(xn)))
        a = reflected_a[j, :, None, None]
        ft = np.exp(-cfg.blockage_m_inv*target_distance)*a/(target_distance**2+a)
        integrand = (pbeam*fd-overlap*ft*fd)*d*jacobian*aw[:, None]*rw[None, :]
        extra = cfg.density_m2*integrand.sum(axis=(1, 2))
        result[j] = np.exp(-reflected_exponent[j]-extra-s[j]*(cfg.noise_w+cfg.communication_w*cfg.residual_si))
    if np.any(result < 0) or np.any(result > 1+1e-12) or np.any(~np.isfinite(result)):
        raise FloatingPointError("Monostatic joint-PGFL numerical error")
    return result
