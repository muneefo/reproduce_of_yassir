"""Printed Eqs. (5), (6), (11)–(20), deliberately separate from simulation.

The printed exclusion/independence approximations are retained, not silently repaired.
"""
import numpy as np
from scipy.integrate import cumulative_trapezoid, quad, trapezoid
from scipy.special import exp1, gammaln, roots_genlaguerre, roots_laguerre, logsumexp
from scipy.stats import qmc


def distance_qmc(cfg, exponent=8, count=6, seed=None):
    u = qmc.Sobol(count, scramble=True, seed=cfg.seed+41 if seed is None else seed).random_base2(exponent)
    if np.any(u <= 0) or np.any(u >= 1):
        raise FloatingPointError("Unexpected QMC endpoint")
    return np.sqrt(np.cumsum(-np.log(u), axis=1)/(np.pi*cfg.density_m2))


def echo_moments(radii, cfg):
    """Bernoulli-exponential moments (5), Gamma shape/scale (6)."""
    r = np.atleast_2d(radii)
    mono = cfg.radar_constant*cfg.rcs_mono_m2/r[:, 0]**(2*cfg.alpha_l)
    bi = cfg.radar_constant*cfg.rcs_bi_m2/(r[:, 0, None]**cfg.alpha_l*r[:, 1:]**cfg.alpha_l)
    p = np.exp(-cfg.blockage_m_inv*r[:, 1:])*cfg.beams*cfg.beamwidth_rad/(2*np.pi)
    mean = mono+np.sum(p*bi, axis=1)
    variance = mono**2+np.sum(p*(2-p)*bi**2, axis=1)
    if cfg.sensing_w == 0:
        return mean, variance, np.full(len(r), np.nan), np.full(len(r), np.nan)
    return mean, variance, mean**2/variance, variance/mean


def reflected_integral(a, lower, blockage):
    """Integral in (13) for alpha_L=2 via complex E1, checked by quadrature.

    ∫_lower^∞ exp(-eta*r) a*r/(r²+a) dr.
    """
    a, lower = np.broadcast_arrays(np.asarray(a), np.asarray(lower))
    if np.any(a < 0) or np.any(lower < 0):
        raise ValueError("Nonnegative Laplace argument and radius required")
    result = np.zeros(a.shape)
    nonzero = a > 0
    root = np.sqrt(a[nonzero])
    result[nonzero] = a[nonzero]*np.real(
        np.exp(-1j*blockage*root)*exp1(blockage*(lower[nonzero]-1j*root)))
    if np.any(~np.isfinite(result)) or np.any(result < 0):
        raise FloatingPointError("E1 reflected integral failed; no clipping allowed")
    return result


def reflected_lt(s, r1, count, cfg, outer_nodes=32):
    """Literal (12)–(13): average R_(Nc+1)|R1; omit its point contribution."""
    if cfg.alpha_l != 2:
        raise NotImplementedError("E1 reduction uses the paper's alpha_L=2")
    s, r1 = np.broadcast_arrays(np.asarray(s), np.asarray(r1))
    nodes, weights = roots_genlaguerre(outer_nodes, count-1)
    weights = weights/np.exp(gammaln(count))
    lower = np.sqrt(r1[..., None]**2+nodes/(np.pi*cfg.density_m2))
    a = s*cfg.radar_constant*cfg.rcs_bi_m2/r1**cfg.alpha_l
    integral = reflected_integral(a[..., None], lower, cfg.blockage_m_inv)
    return np.sum(weights*np.exp(-cfg.density_m2*cfg.beams*cfg.beamwidth_rad*integral), axis=-1)


class DirectLT:
    """Eq. (11), logarithmic radial quadrature plus bounded log-s interpolation.

    Zero argument is exact; out-of-table arguments raise, not extrapolate.
    The NLoS tail beyond r_max uses its small-argument first-order integral.
    """
    def __init__(self, cfg, radial_nodes=2048, outer_nodes=48, s_nodes=2048):
        self.cfg = cfg
        self.log_s = np.linspace(np.log(1e-10), np.log(1e32), s_nodes)
        r = np.geomspace(1e-5, 1e12, radial_nodes)
        nodes, weights = roots_laguerre(outer_nodes)
        lower = np.sqrt(nodes/(np.pi*cfg.density_m2))
        los = np.exp(-cfg.blockage_m_inv*r)
        xi = 2*np.pi*cfg.density_m2*cfg.beams*(cfg.beamwidth_rad/(2*np.pi))**2
        log_values = []
        for start in range(0, s_nodes, 128):
            s = np.exp(self.log_s[start:start+128])[:, None]
            l = s*cfg.communication_w*cfg.gain**2*cfg.intercept_l*r**(-cfg.alpha_l)/cfg.m_l
            n = s*cfg.communication_w*cfg.gain**2*cfg.intercept_n*r**(-cfg.alpha_n)/cfg.m_n
            integrand = (los*(-np.expm1(-cfg.m_l*np.log1p(l)))+
                         (1-los)*(-np.expm1(-cfg.m_n*np.log1p(n))))*r
            tail = (s[:, 0]*cfg.communication_w*cfg.gain**2*cfg.intercept_n*
                    r[-1]**(2-cfg.alpha_n)/(cfg.alpha_n-2))
            # Integrate backwards without total-minus-prefix cancellation.
            cumulative = -cumulative_trapezoid(integrand[:, ::-1], r[::-1], axis=1, initial=0)[:, ::-1]
            at_lower = np.stack([np.interp(np.log(lower), np.log(r), row)
                                 for row in cumulative])+tail[:, None]
            log_value = logsumexp(np.log(weights)[None, :]-xi*at_lower, axis=1)
            if np.any(~np.isfinite(log_value)) or np.any(log_value > 1e-12):
                raise FloatingPointError("Invalid Eq. (11) quadrature")
            log_values.extend(log_value)
        self.log_values = np.asarray(log_values)

    def __call__(self, s):
        s = np.asarray(s)
        if np.any(s < 0):
            raise ValueError("Negative Laplace argument")
        out = np.ones(s.shape)
        nz = s > 0
        log_s = np.log(s[nz])
        if np.any(log_s < self.log_s[0]) or np.any(log_s > self.log_s[-1]):
            raise ValueError("Laplace argument outside validated table; extend it explicitly")
        out[nz] = np.exp(np.interp(log_s, self.log_s, self.log_values))
        return out


def conditional_coverage(radii, thresholds, cfg, direct=None, max_terms=2048,
                         tolerance=2e-4, outer_nodes=32, require_convergence=True):
    """Eq. (14) with real-shape binomial recurrence and an absolute tail bound.

    Returned diagnostics expose truncation errors. Out-of-range probabilities raise.
    """
    r = np.atleast_2d(radii)
    tau = np.atleast_1d(thresholds)
    if np.any(tau < 0):
        raise ValueError("Negative SINR threshold")
    if cfg.sensing_w == 0:
        return np.zeros((len(r), len(tau))), {"terms": 0, "tail_bound": 0.}
    _, _, shape, scale = echo_moments(r, cfg)
    eta = np.exp(-gammaln(1+shape)/shape)
    base = eta[:, None]*tau[None, :]/scale[:, None]
    direct = DirectLT(cfg) if direct is None else direct
    result = np.zeros(base.shape)
    coefficient = shape.copy()
    bound = np.full(base.shape, np.inf)
    positive = tau > 0
    for z in range(1, max_terms+1):
        s = z*base
        lt = (direct(s)*reflected_lt(s, r[:, 0, None], r.shape[1], cfg, outer_nodes)*
              np.exp(-s*(cfg.noise_w+cfg.communication_w*cfg.residual_si)))
        result += coefficient[:, None]*lt
        if z >= np.ceil(shape.max()):
            # All remaining coefficients have the same sign; LT is decreasing.
            mass = np.abs(coefficient)*(z-shape)/shape
            bound = mass[:, None]*lt
            if not positive.any() or np.max(bound[:, positive]) <= tolerance:
                break
        coefficient *= -(shape-z)/(z+1)
    result[:, ~positive] = 1
    bound[:, ~positive] = 0
    worst = float(bound.max())
    if require_convergence and worst > tolerance:
        raise RuntimeError(f"Alzer series unconverged: {z} terms, tail bound {worst:g}")
    if np.any(~np.isfinite(result)) or np.any(result < 0) or np.any(result > 1+1e-12):
        raise FloatingPointError("Invalid conditional probability; investigate, do not clip")
    return result, {"terms": z, "tail_bound": worst, "shape_min": float(shape.min()),
                    "shape_max": float(shape.max())}


def sensing_rate(radii, cfg, t_max=30, points=121, **kwargs):
    t = np.linspace(0, t_max, points)
    p, diagnostic = conditional_coverage(radii, np.expm1(t), cfg, **kwargs)
    average = p.mean(axis=0)
    return cfg.beams*trapezoid(average, t), average[-1], diagnostic


def direct_lt_quad(s, cfg, outer_nodes=32):
    """Slow independent integral check of Eq. (11), not used for plotting."""
    nodes, weights = roots_laguerre(outer_nodes)
    xi = 2*np.pi*cfg.density_m2*cfg.beams*(cfg.beamwidth_rad/(2*np.pi))**2
    total = 0.
    for node, weight in zip(nodes, weights):
        lower = np.sqrt(node/(np.pi*cfg.density_m2))
        def integrand(r):
            p = np.exp(-cfg.blockage_m_inv*r)
            xl = s*cfg.communication_w*cfg.gain**2*cfg.intercept_l*r**(-cfg.alpha_l)/cfg.m_l
            xn = s*cfg.communication_w*cfg.gain**2*cfg.intercept_n*r**(-cfg.alpha_n)/cfg.m_n
            return (p*(-np.expm1(-cfg.m_l*np.log1p(xl)))+
                    (1-p)*(-np.expm1(-cfg.m_n*np.log1p(xn))))*r
        value, error = quad(integrand, lower, np.inf, epsabs=1e-5, epsrel=1e-8, limit=200)
        if error > max(1e-4, abs(value)*1e-6):
            raise RuntimeError("Independent direct integral did not converge")
        total += weight*np.exp(-xi*value)
    return total
