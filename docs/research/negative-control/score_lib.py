"""Scoring computations of the negative-control pre-registration (docs/research/negative-control-preregistration.md,
commit 76c3f99). Shared by the exploratory caller-only sweep (#68) and meant for reuse by the scoring ticket (#67).

The call rule is never re-implemented: `call` is imported from bin/long_reads/genotyping.py at genotype-v3@61d8809,
whose directory is passed as GENOTYPING_DIR (environment variable) or added to sys.path by the caller.
"""
import os
import sys

import numpy as np
from scipy.stats import beta, binom

if os.environ.get("GENOTYPING_DIR"):
    sys.path.insert(0, os.environ["GENOTYPING_DIR"])
from genotyping import call  # noqa: E402  (genotype-v3@61d8809)

SEEDS = [2119795751, 703054139, 967594051, 2137092682, 479082201,
         1922421208, 640385892, 828295203, 1190487198, 31658864,
         1178418751, 1637046186, 1438843744, 2114041616, 913149788,
         174283017, 1834471605, 912037160, 1502725731, 2038317252]
MIN_B_CELLS, MIN_A_UMIS, MIN_EVAL_SEEDS = 100, 200, 10


def split(neg_sorted, seed):
    """Held-out split of the lexicographically sorted negatives: (A, B) as sets (pre-registration §3)."""
    perm = np.random.default_rng(seed).permutation(len(neg_sorted))
    half = len(neg_sorted) // 2
    return {neg_sorted[i] for i in perm[:half]}, {neg_sorted[i] for i in perm[half:]}


def pooled_rate(mut, wt):
    """c = Σ MUT / Σ (WT + MUT), 0 if the denominator is 0 (as genotype-v3's target_contamination)."""
    d = int(np.sum(mut) + np.sum(wt))
    return (float(np.sum(mut)) / d if d else 0.0), int(np.sum(mut)), d


def cp_upper(x, n, conf=0.95):
    """One-sided Clopper-Pearson upper bound of a binomial proportion; 0 if n = 0 (no information, as c = 0)."""
    if n == 0:
        return 0.0
    return 1.0 if x >= n else float(beta.ppf(conf, x + 1, n - x))


def wilson(x, n, z=1.96):
    """Wilson score interval (lo, hi); (nan, nan) if n = 0."""
    if n == 0:
        return float("nan"), float("nan")
    p = x / n
    den = 1 + z * z / n
    mid = (p + z * z / (2 * n)) / den
    half = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return float(mid - half), float(mid + half)


_REGION = {}


def mut_region(n, theta, m_min):
    """Boolean mask over m = 0..n of the counts the imported call labels MUT, for n labelled UMIs."""
    key = (int(n), float(theta), int(m_min))
    if key not in _REGION:
        m = np.arange(n + 1)
        _, calls = call(m, n - m, np.full(n + 1, theta), m_min)
        _REGION[key] = calls == "MUT"
    return _REGION[key]


def p_mut(n, theta, m_min, rho):
    """P(MUT | n_i; θ_i, m_min, ρ_i) per cell: P(X ∈ MUT region), X ~ Bin(n_i, ρ_i) (pre-registration §2)."""
    n = np.asarray(n, dtype=int)
    theta = np.broadcast_to(np.asarray(theta, dtype=float), n.shape)
    rho = np.broadcast_to(np.minimum(np.asarray(rho, dtype=float), 1.0), n.shape)
    if n.size == 0:
        return np.zeros(0)
    keys, inv = np.unique(np.stack([n, theta, rho], axis=1), axis=0, return_inverse=True)
    vals = np.empty(len(keys))
    for j, (ni, ti, ri) in enumerate(keys):
        ni = int(ni)
        mask = mut_region(ni, ti, m_min)
        vals[j] = binom.pmf(np.nonzero(mask)[0], ni, ri).sum() if mask.any() else 0.0
    return vals[np.ravel(inv)]


def poisson_binomial_sf(probs, k):
    """P(S >= k) for S a sum of independent Bernoulli(probs), by dynamic programming on P(S = j), j < k."""
    if k <= 0:
        return 1.0
    dist = np.zeros(k)
    dist[0] = 1.0
    for p in probs:
        dist[1:] = dist[1:] * (1 - p) + dist[:-1] * p
        dist[0] *= 1 - p
    return float(max(0.0, 1.0 - dist.sum()))


def ratio(o, e):
    """O / E with the pre-registered conventions: ∞ if E = 0 < O, 1 if O = E = 0."""
    if e == 0:
        return float("inf") if o > 0 else 1.0
    return o / e


def em_mixture(m, n, rho0, tol=1e-8, max_iter=10_000):
    """ML fit of m_i ~ (1-π) Bin(n_i, ρ0) + π Bin(n_i, f), ρ0 fixed, f ∈ (ρ0, 1], by EM from π = f = 0.5."""
    m, n = np.asarray(m, float), np.asarray(n, float)
    pi, f = 0.5, max(0.5, rho0 + 1e-6)
    l0 = binom.logpmf(m, n, min(rho0, 1.0))
    last = -np.inf
    for _ in range(max_iter):
        l1 = binom.logpmf(m, n, f)
        with np.errstate(divide="ignore"):
            a = np.log(1 - pi) + l0
            b = np.log(pi) + l1
        mx = np.maximum(a, b)
        s = np.exp(a - mx) + np.exp(b - mx)
        ll = float(np.sum(mx + np.log(s)))
        w = np.exp(b - mx) / s
        pi = float(w.mean())
        f = float(np.clip((w * m).sum() / (w * n).sum(), rho0 + 1e-9, 1.0)) if (w * n).sum() > 0 else f
        if abs(ll - last) < tol:
            break
        last = ll
    return pi, f


def efficiency(o_t, e_t, pi, f, n_cov, rho0):
    """Detection efficiency (O_T − (1 − π̂) E_T) / (π̂ n_cov,T); NaN when π̂ < 0.05 or f̂ < ρ0 + 0.05."""
    if pi < 0.05 or f < rho0 + 0.05 or n_cov == 0:
        return float("nan")
    return (o_t - (1 - pi) * e_t) / (pi * n_cov)


def seed_summary(values):
    """Median and 5–95 % range over seeds (numpy percentile, linear), NaN-free values only."""
    v = np.asarray([x for x in values if not np.isnan(x)], float)
    if v.size == 0:
        return float("nan"), float("nan"), float("nan")
    return tuple(float(x) for x in np.percentile(v, [50, 5, 95]))
