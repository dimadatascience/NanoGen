"""Measure the hybrid UMI merge rule on PALM14E (wayfinder ticket #13).

Hybrid rule (UMI grouping redesign, #12): within a cell x target bundle of primary records,
UMI a absorbs b if LD(a, b) <= 2 and either n_a >= 2 n_b - 1 (directional ratio) or
n_b <= T (small clause). Clusters grow downhill from the most abundant UMI, so the small
clause also needs n_a >= n_b; the "literal" scheme drops that guard to show whether it matters.

Reads the per-read table cached by followup.py (records.pkl; stays on the HPC) and writes
AGGREGATE tables only (targets as T1..Tn):
  H1 merges: clusters vs plain directional, UBs/reads absorbed by the small clause,
     molecules (clusters with >= n_min reads) per cell x target
  H2 purity: minority-allele rate of reads absorbed through each clause vs intrinsic
     (exact UB, n >= 2; same definition as F1 / check 4)
  H3 share of clusters and reads below n_min
  H4 main's caller on the n_min consensus: calls per target and transitions vs plain
"""
import argparse
import os
from collections import defaultdict

import numpy as np
import pandas as pd

from measure_umi import FMIN, call_main, close_pairs, levenshtein, load_targets

MAXD = 2
NMIN = 3
# name -> (small-clause threshold T, require n_a >= n_b for the small clause)
SCHEMES = {"plain": (0, True), "T2": (2, True), "T3": (3, True), "T5": (5, True),
           "T3-literal": (3, False)}


def hybrid(umis, nbrs, small, downhill):
    """umis: list of (umi, count) sorted by (-count, umi); nbrs: index -> neighbour indices.
    Returns (rep, via): representative index and absorbing clause per index."""
    rep, via = {}, {}
    for i in range(len(umis)):
        if i in rep:
            continue
        rep[i], via[i] = i, "seed"
        stack = [i]
        while stack:
            x = stack.pop()
            na = umis[x][1]
            for y in nbrs[x]:
                if y in rep:
                    continue
                nb = umis[y][1]
                if na >= 2 * nb - 1:
                    v = "ratio"
                elif nb <= small and (na >= nb or not downhill):
                    v = "small"
                else:
                    continue
                rep[y], via[y] = i, v
                stack.append(y)
    return rep, via


def assign(R1):
    """Add one cluster column and one via column per scheme."""
    out = {f"cl_{k}": np.empty(len(R1), dtype=object) for k in SCHEMES}
    out.update({f"via_{k}": np.empty(len(R1), dtype=object) for k in SCHEMES})
    for (cell, tid), idx in R1.groupby(["cell", "tid"]).indices.items():
        ubs = R1.ub.values[idx]
        vc = pd.Series(ubs).value_counts()
        umis = sorted(vc.items(), key=lambda x: (-x[1], x[0]))
        us = [u for u, _ in umis]
        nbrs = defaultdict(list)
        for i, j, _ in close_pairs(us, MAXD, levenshtein):
            nbrs[i].append(j)
            nbrs[j].append(i)
        pos = {u: i for i, u in enumerate(us)}
        ui = np.array([pos[u] for u in ubs])
        for k, (small, downhill) in SCHEMES.items():
            rep, via = hybrid(umis, nbrs, small, downhill)
            out[f"cl_{k}"][idx] = [us[rep[i]] for i in ui]
            out[f"via_{k}"][idx] = [via[i] for i in ui]
    return R1.assign(**out)


def consensus(R1, col):
    """Per (cell, tid, cluster): n reads, labelled counts, uniform-rule label at NMIN / FMIN."""
    g = R1.groupby(["cell", "tid", col])
    n = g.size().rename("n")
    cnt = R1[R1.allele.notna()].groupby(["cell", "tid", col, "allele"]).size().unstack(fill_value=0)
    cnt = cnt.reindex(columns=["REF", "ALT", "OTH"], fill_value=0).reindex(n.index, fill_value=0)
    tot = cnt.sum(axis=1).values
    frac = cnt.values.max(axis=1) / np.maximum(tot, 1)
    lab = np.array(["WT", "MUT", "OTH"])[cnt.values.argmax(axis=1)]
    ok = (n.values >= NMIN) & (tot > 0) & (frac >= FMIN)
    return pd.DataFrame({"n": n.values, "tot": tot, "label": np.where(ok, lab, "drop")}, index=n.index)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bed", required=True)
    ap.add_argument("--cache", default="records.pkl")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    T = load_targets(a.bed)
    R = pd.read_pickle(a.cache)
    R1 = R[R.kind == "prim"].drop_duplicates(["cell", "qname", "pos", "tid"]).reset_index(drop=True)
    R1 = assign(R1)
    lab = R1[R1.allele.isin(["REF", "ALT"])]

    # ---------------- H1 merges and molecules per cell x target
    h1 = []
    mol = {}
    for k in SCHEMES:
        cl, via = f"cl_{k}", f"via_{k}"
        sz = R1.groupby(["cell", "tid", cl]).size()
        per_ct = (sz >= NMIN).groupby(level=["cell", "tid"]).sum()
        mol[k] = per_ct
        for tid in T.tid:
            s = R1[R1.tid == tid]
            z = sz[sz.index.get_level_values("tid") == tid]
            p = per_ct[per_ct.index.get_level_values("tid") == tid]
            ubv = s.drop_duplicates(["cell", "ub"])[via]
            h1.append({"scheme": k, "tid": tid, "ubs": len(ubv), "clusters": len(z),
                       "ubs_via_ratio": int((ubv == "ratio").sum()), "ubs_via_small": int((ubv == "small").sum()),
                       "reads": len(s), "reads_via_small": int((s[via] == "small").sum()),
                       "clusters_ge_nmin": int((z >= NMIN).sum()),
                       "mol_per_ct_mean": p[p > 0].mean(), "mol_per_ct_median": p[p > 0].median(),
                       "mol_per_ct_p90": p[p > 0].quantile(.9), "mol_per_ct_max": p.max()})
    h1 = pd.DataFrame(h1)
    base = h1[h1.scheme == "plain"].set_index("tid")
    h1["extra_merges_vs_plain"] = [base.clusters[t] - c for t, c in zip(h1.tid, h1.clusters)]
    h1["delta_clusters_ge_nmin_vs_plain"] = [c - base.clusters_ge_nmin[t] for t, c in zip(h1.tid, h1.clusters_ge_nmin)]
    for k in SCHEMES:
        j = pd.concat([mol["plain"], mol[k]], axis=1, keys=["p", "h"]).fillna(0)
        j = j[(j.p > 0) | (j.h > 0)]
        for tid in T.tid:
            jj = j[j.index.get_level_values("tid") == tid]
            m = (h1.scheme == k) & (h1.tid == tid)
            h1.loc[m, "ct_covered"] = len(jj)
            h1.loc[m, "ct_mol_changed"] = int((jj.p != jj.h).sum())
            h1.loc[m, "ct_mol_up"] = int((jj.h > jj.p).sum())
            h1.loc[m, "ct_mol_down"] = int((jj.h < jj.p).sum())
    h1.to_csv(os.path.join(a.out, "H1_merges.csv"), index=False)

    # ---------------- H2 purity
    g = lab.groupby(["cell", "tid", "ub"]).allele
    L = lab.join(g.agg(lambda x: x.value_counts().idxmax()).rename("ub_maj"), on=["cell", "tid", "ub"]) \
           .join(g.size().rename("ub_n"), on=["cell", "tid", "ub"])
    h2 = []
    for tid in T.tid:
        s = L[(L.tid == tid) & (L.ub_n >= 2)]
        h2.append({"scheme": "-", "tid": tid, "class": "intrinsic (exact UB, n>=2)", "reads": len(s),
                   "minority_frac": (s.allele != s.ub_maj).mean() if len(s) else np.nan})
    for k in SCHEMES:
        cl, via = f"cl_{k}", f"via_{k}"
        cc = lab.groupby(["cell", "tid", cl, "allele"]).size().unstack(fill_value=0) \
                .reindex(columns=["REF", "ALT"], fill_value=0)
        uc = lab.groupby(["cell", "tid", "ub", "allele"]).size().unstack(fill_value=0) \
                .reindex(columns=["REF", "ALT"], fill_value=0)
        L2 = lab.join(cc, on=["cell", "tid", cl]).join(uc, on=["cell", "tid", "ub"], rsuffix="_ub")
        # cluster majority with the read's own UB, and leave-own-UB-out (ties excluded)
        L2["maj"] = np.where(L2.REF >= L2.ALT, "REF", "ALT")
        rr, ra = L2.REF - L2.REF_ub, L2.ALT - L2.ALT_ub
        L2["maj_loo"] = np.where(rr > ra, "REF", np.where(ra > rr, "ALT", None))
        for tid in T.tid:
            for v in ["ratio", "small"]:
                s = L2[(L2.tid == tid) & (L2[via] == v)]
                sl = s[s.maj_loo.notna()]
                h2.append({"scheme": k, "tid": tid, "class": f"absorbed via {v}", "reads": len(s),
                           "minority_frac": (s.allele != s.maj).mean() if len(s) else np.nan,
                           "reads_loo": len(sl),
                           "minority_frac_loo": (sl.allele != sl.maj_loo).mean() if len(sl) else np.nan})
    pd.DataFrame(h2).to_csv(os.path.join(a.out, "H2_purity.csv"), index=False)

    # ---------------- H3 below n_min, H4 calls
    h3, calls, ct = [], [], {}
    cells = R1.cell.unique()
    grid = pd.MultiIndex.from_product([cells, T.tid], names=["cell", "tid"]).to_frame(index=False)
    for k in SCHEMES:
        C = consensus(R1, f"cl_{k}")
        for tid in list(T.tid) + ["all"]:
            c = C if tid == "all" else C[C.index.get_level_values("tid") == tid]
            labs = c.label.value_counts()
            h3.append({"scheme": k, "tid": tid, "clusters": len(c),
                       "frac_clusters_below_nmin": (c.n < NMIN).mean(),
                       "frac_reads_below_nmin": c.n[c.n < NMIN].sum() / c.n.sum(),
                       "WT_umis": int(labs.get("WT", 0)), "MUT_umis": int(labs.get("MUT", 0)),
                       "OTH_umis": int(labs.get("OTH", 0)), "dropped_fmin": int(((c.n >= NMIN) & (c.label == "drop")).sum())})
        tab = C[C.label != "drop"].groupby(["cell", "tid", "label"]).size().unstack(fill_value=0) \
                                  .reindex(columns=["WT", "MUT", "OTH"], fill_value=0).reset_index()
        tab = grid.merge(tab, on=["cell", "tid"], how="left").fillna(0)
        ct[k] = call_main(tab).set_index(["cell", "tid"])
        for tid, r in ct[k].groupby(level="tid").genotype.value_counts().unstack(fill_value=0).iterrows():
            calls.append({"scheme": k, "tid": tid, **r.to_dict()})
    pd.DataFrame(h3).to_csv(os.path.join(a.out, "H3_nmin.csv"), index=False)
    pd.DataFrame(calls).fillna(0).to_csv(os.path.join(a.out, "H4_calls_per_target.csv"), index=False)
    tr = []
    for k in SCHEMES:
        j = ct["plain"][["genotype"]].join(ct[k][["genotype"]], rsuffix="_new")
        for (tid, g0, g1), n in j.groupby([j.index.get_level_values("tid"), "genotype", "genotype_new"]).size().items():
            tr.append({"scheme": k, "tid": tid, "from_plain": g0, "to": g1, "n": int(n)})
    pd.DataFrame(tr).to_csv(os.path.join(a.out, "H4_call_transitions.csv"), index=False)
    print("done")


if __name__ == "__main__":
    main()
