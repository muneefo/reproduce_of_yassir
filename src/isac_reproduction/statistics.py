"""Intervals and meta-distribution calculations, with no probability clipping."""
import numpy as np
from scipy.integrate import trapezoid
from scipy.stats import norm, t as student_t


def wilson(successes, trials, confidence=.95):
    p = np.asarray(successes)/trials
    z = norm.ppf((1+confidence)/2)
    # Rationalized Wilson root avoids cancellation at zero/all successes.
    # This is algebraically the usual interval, not clipping an invalid result.
    def lower(q):
        return 2*trials*q*q/(2*trials*q+z*z+z*np.sqrt(z*z+4*trials*q*(1-q)))
    return np.stack((lower(p), 1-lower(1-p)), axis=-1)


def mean_interval(samples, axis=0):
    samples = np.asarray(samples)
    n = samples.shape[axis]
    mean = samples.mean(axis=axis)
    half = student_t.ppf(.975, n-1)*samples.std(axis=axis, ddof=1)/np.sqrt(n)
    return mean, np.stack((mean-half, mean+half), axis=-1)


def coverage(samples, threshold_db):
    thresholds = 10**(np.asarray(threshold_db)/10)
    hits = samples[:, :, None] > thresholds
    return hits.mean(axis=0), wilson(hits.sum(axis=0), len(samples))


def moments(probabilities, b):
    p = np.asarray(probabilities)
    if np.any(p <= 0) or np.any(p > 1) or np.any(~np.isfinite(p)):
        raise ValueError("Complex moments require finite probabilities in (0,1]")
    return np.exp(np.asarray(b)[..., None]*np.log(p)).mean(axis=-1)


def meta_ccdf(probabilities, t):
    return (np.asarray(probabilities)[:, None] >= np.asarray(t)).mean(axis=0)


def gil_pelaez(probabilities, t, u_max=200, steps=4001):
    """Eq. (20) truncated inversion. Returns raw values, including ringing."""
    p = np.asarray(probabilities)
    t = np.asarray(t)
    if np.any(t <= 0) or np.any(t >= 1):
        raise ValueError("Inversion is evaluated on 0<t<1; endpoints handled separately")
    u = np.linspace(0, u_max, steps)
    characteristic = moments(p, 1j*u)
    integral = np.empty((len(t), len(u)))
    integral[:, 0] = np.log(p).mean()-np.log(t)
    integral[:, 1:] = np.imag(np.exp(-1j*np.log(t)[:, None]*u[None, 1:])*
                                     characteristic[None, 1:])/u[None, 1:]
    return .5+trapezoid(integral, u, axis=1)/np.pi
