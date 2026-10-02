"""Follow-up to measure_umi.py: are directional-clustering merges genuine siblings?

F1 read-level purity: share of WT/MUT reads disagreeing with their group's majority, for
   exact-UB groups (intrinsic read error) vs reads absorbed from merged-in child UBs.
F2 parent-child label concordance for children with >= 2 labelled reads, by distance.
F3 reads per molecule and molecules per (cell, target) under B1 and B1+C.
F4 current-scheme inflation: passing current groups per passing B1+C(LD<=2) cluster.
Caches the per-read table (stays on the HPC; patient data) as records.pkl.
"""
import argparse
import glob
import os
from collections import Counter

import numpy as np
import pandas as pd
import pysam

from measure_umi import MAPQ, WIN, directional, hamming, levenshtein, load_targets, site_alleles


def build_records(cells, T):
    rows = []
    for cd in sorted(glob.glob(os.path.join(cells, "*"))):
        cell = os.path.basename(cd)
        bam = os.path.join(cd, "filtered_input.bam")
        if not os.path.exists(bam):
            continue
        if not os.path.exists(bam + ".bai"):
            bam = os.path.join(cd, "sorted_fi.bam")
        with pysam.AlignmentFile(bam, "rb") as f:
            for t in T.itertuples(index=False):
                if t.chr not in f.references:
                    continue
                for r, al in site_alleles(f, t):
                    if r.mapping_quality < MAPQ:
                        continue
                    kind = "supp" if r.is_supplementary else ("sec" if r.is_secondary else "prim")
                    ub = r.get_tag("UB") if r.has_tag("UB") else None
                    bx = r.get_tag("BX") if r.has_tag("BX") else ub
                    rows.append((cell, t.tid, r.query_name, ub, bx, r.reference_start + 1, kind, al))
    return pd.DataFrame(rows, columns=["cell", "tid", "qname", "ub", "bx", "pos", "kind", "allele"])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cells", required=True)
    ap.add_argument("--bed", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--cache", default="records.pkl")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    T = load_targets(a.bed)
    if os.path.exists(a.cache):
        R = pd.read_pickle(a.cache)
    else:
        R = build_records(a.cells, T)
        R.to_pickle(a.cache)
    R1 = R[R.kind == "prim"].copy()

    # cluster representatives per scheme, plus parent-of-child edges
    schemes = {"H1": (1, hamming), "LD2": (2, levenshtein), "LD3": (3, levenshtein)}
    reps = {k: {} for k in schemes}
    for (cell, tid), sub in R1.groupby(["cell", "tid"]):
        u = sub.ub.value_counts().to_dict()
        for k, (d, fn) in schemes.items():
            for x, rep in directional(u, d, fn).items():
                reps[k][(cell, tid, x)] = rep
    for k in schemes:
        R1[k] = [reps[k][(c, t, x)] for c, t, x in zip(R1.cell, R1.tid, R1.ub)]

    lab = R1[R1.allele.isin(["REF", "ALT"])]

    # F1 read-level purity
    f1 = []
    g = lab.groupby(["cell", "tid", "ub"]).allele
    maj = g.agg(lambda x: x.value_counts().idxmax()).rename("ub_maj")
    n = g.size().rename("ub_n")
    L = lab.join(maj, on=["cell", "tid", "ub"]).join(n, on=["cell", "tid", "ub"])
    for tid, sub in L.groupby("tid"):
        s2 = sub[sub.ub_n >= 2]
        f1.append({"tid": tid, "scheme": "exact UB (n>=2): intrinsic",
                   "reads": len(s2), "minority_frac": (s2.allele != s2.ub_maj).mean()})
    for k in schemes:
        cm = lab.groupby(["cell", "tid", k]).allele.agg(lambda x: x.value_counts().idxmax()).rename("cl_maj")
        L2 = lab.join(cm, on=["cell", "tid", k])
        child = L2[L2.ub != L2[k]]
        par = L2[L2.ub == L2[k]]
        for tid in T.tid:
            c = child[child.tid == tid]
            p = par[par.tid == tid]
            f1.append({"tid": tid, "scheme": f"{k}: reads from merged-in children",
                       "reads": len(c), "minority_frac": (c.allele != c.cl_maj).mean() if len(c) else np.nan})
            f1.append({"tid": tid, "scheme": f"{k}: reads of cluster parent UB",
                       "reads": len(p), "minority_frac": (p.allele != p.cl_maj).mean() if len(p) else np.nan})
    pd.DataFrame(f1).to_csv(os.path.join(a.out, "F1_read_purity.csv"), index=False)

    # F2 parent-child label concordance (both with >= 2 labelled reads)
    ublab = pd.concat([maj, n], axis=1).reset_index()
    ublab = ublab[ublab.ub_n >= 2].set_index(["cell", "tid", "ub"]).ub_maj.to_dict()
    f2 = []
    for k, (d, fn) in schemes.items():
        acc = Counter()
        for (cell, tid, ub), rep in reps[k].items():
            if ub == rep:
                continue
            if (cell, tid, ub) not in ublab or (cell, tid, rep) not in ublab:
                continue
            dist = levenshtein(ub, rep)
            disc = ublab[(cell, tid, ub)] != ublab[(cell, tid, rep)]
            acc[(tid, dist, "n")] += 1
            acc[(tid, dist, "disc")] += disc
        for (tid, dist, w), v in list(acc.items()):
            if w == "n":
                f2.append({"scheme": k, "tid": tid, "child_parent_dist": dist, "pairs": v,
                           "discordant": acc[(tid, dist, "disc")], "disc_rate": acc[(tid, dist, "disc")] / v})
    pd.DataFrame(f2, columns=["scheme", "tid", "child_parent_dist", "pairs", "discordant", "disc_rate"]).sort_values(["scheme", "tid", "child_parent_dist"]).to_csv(
        os.path.join(a.out, "F2_parent_child.csv"), index=False)

    # F3 reads per molecule / molecules per cell x target
    f3 = []
    for k in ["ub", "H1", "LD2", "LD3"]:
        sz = R1.groupby(["cell", "tid", k]).size()
        for tid in T.tid:
            s = sz[sz.index.get_level_values("tid") == tid]
            per_ct = s.groupby(level=["cell", "tid"]).size()
            per_ct2 = (s >= 2).groupby(level=["cell", "tid"]).sum()
            f3.append({"grouping": k, "tid": tid, "molecules": len(s), "frac_singletons": (s == 1).mean(),
                       "reads_in_singletons": s[s == 1].sum() / s.sum(),
                       "median_reads_mol_ge2": s[s >= 2].median(), "p90_reads_mol_ge2": s[s >= 2].quantile(.9),
                       "median_mol_ge2_per_cell": per_ct2[per_ct2 > 0].median(),
                       "p90_mol_ge2_per_cell": per_ct2[per_ct2 > 0].quantile(.9),
                       "max_mol_ge2_per_cell": per_ct2.max()})
    pd.DataFrame(f3).to_csv(os.path.join(a.out, "F3_molecule_sizes.csv"), index=False)

    # F4 inflation of current scheme: passing current groups (BX x win100, n>=2) per LD2 cluster
    Rc = R.copy()
    Rc["win"] = (Rc.pos // WIN) * WIN
    cur = Rc.groupby(["cell", "tid", "bx", "win"]).size()
    cur = cur[cur >= 2].reset_index()
    # map each current group to the LD2 cluster of its BX (= UB for 98% of reads)
    cur["cl"] = [reps["LD2"].get((c, t, b)) for c, t, b in zip(cur.cell, cur.tid, cur.bx)]
    f4 = cur.dropna(subset=["cl"]).groupby(["tid", "cell", "cl"]).size().groupby("tid").agg(
        ["mean", "median", lambda x: (x > 1).mean()])
    f4.columns = ["current_groups_per_LD2_cluster_mean", "median", "frac_clusters_counted_more_than_once"]
    f4.to_csv(os.path.join(a.out, "F4_current_inflation.csv"))
    print("done")


if __name__ == "__main__":
    main()
