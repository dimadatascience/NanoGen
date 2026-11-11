"""Compare guards for a pooled calibration test against a single-site miss.

O_i ~ Poisson(r_i * E_i) per target. Aggregate output only.
"""
import sys, json
import numpy as np
from scipy import stats

rng = np.random.default_rng(20261111)
NSIM = int(sys.argv[2]) if len(sys.argv) > 2 else 100000
NCAL = 4000  # multinomial draws per N for conditional MC critical values

rows = [l.rstrip("\n").split("\t") for l in open(sys.argv[1])][1:]
pat = np.array([r[0] for r in rows])
E_all = {"c1": np.array([float(r[2]) for r in rows]), "c2": np.array([float(r[3]) for r in rows])}
keep = E_all["c1"] > 0  # E = 0 targets carry no information
pat = pat[keep]
for k in E_all:
    E_all[k] = E_all[k][keep]
upat, pidx = np.unique(pat, return_inverse=True)
K = keep.sum()


def pois_sf(O, E):
    # P(X >= O)
    return stats.poisson.sf(O - 1, E)


def G_stat(O, Eexp):
    # LR statistic for common O/E; Eexp already scaled to sum O
    with np.errstate(divide="ignore", invalid="ignore"):
        t = np.where(O > 0, O * np.log(O / Eexp), 0.0)
    return 2 * t.sum(-1)


def X2_stat(O, Eexp):
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(Eexp > 0, (O - Eexp) ** 2 / Eexp, 0.0).sum(-1)


def het(O, E, crit_G, crit_X):
    N = O.sum(-1)
    w = E / E.sum()
    Eexp = N[:, None] * w[None, :]
    G = G_stat(O, Eexp)
    X = X2_stat(O, Eexp)
    df = O.shape[1] - 1
    asym = stats.chi2.sf(G, df) < 0.05
    Nc = np.minimum(N, len(crit_G) - 1).astype(int)
    return G > crit_G[Nc], X > crit_X[Nc], asym


def calib(E, Nmax):
    w = E / E.sum()
    cG = np.full(Nmax + 1, np.inf)
    cX = np.full(Nmax + 1, np.inf)
    for N in range(1, Nmax + 1):
        O = rng.multinomial(N, w, size=NCAL).astype(float)
        Eexp = N * w[None, :]
        cG[N] = np.quantile(G_stat(O, Eexp), 0.95)
        cX[N] = np.quantile(X2_stat(O, Eexp), 0.95)
    return cG, cX


def holm(p, a):
    m = p.shape[1]
    o = np.argsort(p, 1)
    ps = np.take_along_axis(p, o, 1)
    thr = a / (m - np.arange(m))
    ok = ps <= thr
    # reject first j where all previous ok
    run = np.cumprod(ok, 1).astype(bool)
    rej = np.zeros_like(run)
    np.put_along_axis(rej, o, run, 1)
    return rej


def bh(p, a):
    m = p.shape[1]
    o = np.argsort(p, 1)
    ps = np.take_along_axis(p, o, 1)
    thr = a * (np.arange(1, m + 1)) / m
    ok = ps <= thr
    kmax = np.where(ok.any(1), m - np.argmax(ok[:, ::-1], 1), 0)
    run = np.arange(m)[None, :] < kmax[:, None]
    rej = np.zeros_like(run)
    np.put_along_axis(rej, o, run, 1)
    return rej


out = {"K_targets_nonzero": int(K)}
for ck in ["c1", "c2"]:
    for scale in [1.0, 1.27]:
        E = E_all[ck] * scale
        Ep = np.bincount(pidx, E)
        Etot = E.sum()
        imax = int(np.argmax(E))
        order = np.argsort(E)
        imed = int(order[len(E) // 2])
        Nmax = int(Etot * 12 + 60)
        cGt, cXt = calib(E, Nmax)
        cGp, cXp = calib(Ep, Nmax)
        cfg = f"{ck}_x{scale}"
        res = {"E_pool": Etot, "E_max": E[imax], "E_med": E[imed],
               "share_max": E[imax] / Etot, "share_med": E[imed] / Etot,
               "E_patients_sorted": sorted(Ep.round(2).tolist()),
               "n_E_lt1": int((E < 1).sum()), "scen": {}}
        scen = [("null", None, 1.0), ("unif1.5", None, 1.5)]
        for rs in [1.7, 2, 3, 5, 10]:
            scen += [(f"max_r{rs}", imax, rs), (f"med_r{rs}", imed, rs)]
        for name, site, rs in scen:
            r = np.ones(K) * (rs if site is None else 1.0)
            if site is not None:
                r[site] = rs
            O = rng.poisson(r * E, size=(NSIM, K)).astype(float)
            Ot = O.sum(1)
            p_pool = pois_sf(Ot, Etot)
            pool_fail = (p_pool < 0.05) | (Ot / Etot > 1.5)
            p_t = pois_sf(O, E[None, :])
            d = {"pool_fail": pool_fail.mean()}
            for meth, f in [("holm", holm), ("bh", bh)]:
                for a in [0.1, 0.05]:
                    rej = f(p_t, a)
                    anyf = rej.any(1)
                    key = f"{meth}{a}"
                    d[key + "_anyflag"] = anyf.mean()
                    d[key + "_pool_or_flag"] = (pool_fail | anyf).mean()
                    if site is not None:
                        d[key + "_siteflag"] = rej[:, site].mean()
                        d[key + "_pool_or_siteflag"] = (pool_fail | rej[:, site]).mean()
            # unadjusted per-target p<0.05 for reference
            d["unadj0.05_anyflag"] = (p_t < 0.05).any(1).mean()
            gT, xT, aT = het(O, E, cGt, cXt)
            Op = np.zeros((NSIM, len(Ep)))
            for j in range(len(Ep)):
                Op[:, j] = O[:, pidx == j].sum(1)
            gP, xP, aP = het(Op, Ep, cGp, cXp)
            d.update({"hetT_G_mc": gT.mean(), "hetT_X2_mc": xT.mean(), "hetT_G_asym": aT.mean(),
                      "hetP_G_mc": gP.mean(), "hetP_X2_mc": xP.mean(), "hetP_G_asym": aP.mean(),
                      "pool_or_hetT_G": (pool_fail | gT).mean(), "pool_or_hetT_X2": (pool_fail | xT).mean(),
                      "pool_or_hetP_G": (pool_fail | gP).mean()})
            res["scen"][name] = {k: round(float(v), 4) for k, v in d.items()}
        out[cfg] = res
        print(cfg, "done", flush=True)

# Minimum E for a per-target one-sided Poisson test to detect r with 80% power
def min_E(r, a, target=0.8):
    grid = np.round(np.arange(0.05, 200, 0.05), 3)
    for e in grid:
        # smallest O with sf <= a
        o = stats.poisson.isf(a, e) + 1
        while stats.poisson.sf(o - 2, e) <= a:
            o -= 1
        if stats.poisson.sf(o - 1, r * e) >= target:
            return float(e)
    return None

me = {}
for r in [1.7, 2, 3, 5, 10]:
    for lab, a in [("a0.05_unadj", 0.05), ("a0.1_over_K", 0.1 / K), ("a0.05_over_K", 0.05 / K)]:
        me[f"r{r}_{lab}"] = min_E(r, a)
out["minE_80"] = me
# power of the per-target test at the largest / median E (Bonferroni-first-step alpha) for r=3
json.dump(out, open("pool_guard_out.json", "w"), indent=1)
print(json.dumps(out["minE_80"], indent=1))
