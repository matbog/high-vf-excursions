"""Area-Palm, small-margin normal film-width law in three dimensions.

For iid unit Gaussian components with C(r)=exp(-(r/lc)**2), the limiting
linearized width is delta*lc/R, R~chi_4. This is not an exact finite-margin
width distribution or a lower bound on physical gaps.
"""
from __future__ import annotations

import numpy as np
from scipy.special import gamma
from scipy.stats import chi

from theory.check_margin_tubes import expected_maximum


def width_scale(delta, lc):
    if not np.isfinite(delta) or delta <= 0 or not np.isfinite(lc) or lc <= 0:
        raise ValueError("delta and lc must be finite and positive")
    return delta*lc


def width_cdf(width, delta, lc):
    """Area fraction whose linearized normal width is <= width."""
    scale = width_scale(delta, lc)
    width = np.asarray(width, dtype=float)
    if np.any(np.isnan(width)):
        raise ValueError("Width cannot be NaN")
    ratio = np.divide(scale, width, out=np.full_like(width, np.inf), where=width > 0)
    return chi.sf(ratio, df=4)


def width_quantile(probability, delta, lc):
    scale = width_scale(delta, lc)
    p = np.asarray(probability, dtype=float)
    if np.any(~np.isfinite(p)) or np.any((p <= 0) | (p >= 1)):
        raise ValueError("Quantile probability must be strictly between zero and one")
    return scale/chi.isf(p, df=4)


def width_moment(power, delta, lc):
    scale = width_scale(delta, lc)
    if not np.isfinite(power) or power < 0:
        raise ValueError("Require a finite nonnegative moment order")
    if power >= 4:
        return float("inf")
    return float(scale**power * 2**(-power/2) * gamma((4-power)/2))


def interface_area_density(q, lc):
    """Mean area per volume of zero-margin winner ties, each face counted once."""
    if not isinstance(q, (int, np.integer)) or q < 2:
        raise ValueError("Require integer q>=2")
    width_scale(1., lc)
    return float(2*np.sqrt(2/np.pi)*expected_maximum(q)/lc)


def weighted_quantile(values, weights, probabilities):
    values, weights = np.asarray(values), np.asarray(weights)
    if len(values) == 0 or np.any(weights <= 0):
        raise ValueError("Require nonempty values and positive weights")
    order = np.argsort(values)
    cdf = np.cumsum(weights[order])/np.sum(weights)
    return np.interp(probabilities, cdf, values[order])


def weighted_cdf_distance(values, weights, cdf):
    """Descriptive weighted KS distance; not an iid hypothesis-test p-value."""
    order = np.argsort(values)
    weights = np.asarray(weights)[order]
    theory = cdf(np.asarray(values)[order])
    upper = np.cumsum(weights)/np.sum(weights)
    lower = upper-weights/np.sum(weights)
    return float(max(np.max(abs(upper-theory)), np.max(abs(lower-theory))))
