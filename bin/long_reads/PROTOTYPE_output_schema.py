#!/usr/bin/env python
"""
PROTOTYPE - throwaway. Not production code; do not import or wire into Nextflow.

Answers wayfinder ticket "Output schema: OTH counts, per-target error and
limitation reporting" (dimadatascience/NanoGen#7): what <sample>.csv and
genotype.log look like after genotype-v2.

Runs on a tiny synthetic counts table (no patient data) chosen to hit the edge
cases: m=1 with few/many WT UMIs, m>=2, m=0 with k<2, N=0, an indel target with
WT->ALT leakage, a target with w=0 (floor applied).

    python bin/long_reads/PROTOTYPE_output_schema.py            # proposed schema
    python bin/long_reads/PROTOTYPE_output_schema.py --mis-alias # + MIS compat column
    python bin/long_reads/PROTOTYPE_output_schema.py --outdir x # also write the files
"""

import argparse
import os

import numpy as np
import pandas as pd
from scipy.stats import binom

EPS = 1e-4
P_MUT, P_WT = 0.1, 0.9

# --- Proposed per-cell table from count_consensus.py (one row per cell x target) ---
# New vs main: MIS -> OTH; three read-level columns for the indel ALT-side diagnostic,
# filled for every target (cheap) but only *reported* for indels.
#   wt_reads          reads (bq>=20) inside WT-labelled UMIs
#   wt_altlike_reads  of those, reads matching ALT under the tolerant indel match
#   mixed_discarded   UMIs with >= nmin reads, discarded because REF and ALT both present
#                     and neither reached fmin
TARGETS = [
    # chr,   start,     end,       ref, alt,    gene
    ("chr4", 54733155, 54733155, "A", "T", "KIT"),     # SNV transversion
    ("chr4", 54733160, 54733162, "GAT", "-GAT", "KIT"),  # 3-bp deletion
    ("chr12", 25245350, 25245350, "C", "T", "KRAS"),   # SNV, w = 0 -> floor
]
# cell -> per target (WT, MUT, OTH, wt_reads, wt_altlike_reads, mixed_discarded)
CELLS = {
    "AAAC-1": [(1, 1, 0, 6, 0, 0), (2, 0, 0, 14, 1, 0), (0, 0, 0, 0, 0, 0)],
    "AAAG-1": [(0, 3, 0, 0, 0, 0), (0, 2, 0, 0, 0, 1), (4, 0, 0, 30, 0, 0)],
    "AAAT-1": [(2, 0, 1, 11, 0, 0), (1, 0, 0, 5, 1, 0), (1, 0, 0, 4, 0, 0)],
    "AACA-1": [(5000, 1, 3, 40000, 2, 0), (900, 0, 1, 7000, 210, 12), (300, 0, 0, 2400, 0, 0)],
    "AACC-1": [(0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0)],
    "AACG-1": [(10, 2, 0, 70, 0, 0), (3, 1, 0, 20, 2, 0), (0, 1, 0, 0, 0, 0)],
}


def counts_table():
    rows = []
    for cell, per_target in CELLS.items():
        for (chrom, start, end, ref, alt, gene), vals in zip(TARGETS, per_target):
            rows.append(dict(chr=chrom, start=start, end=end, ref=ref, alt=alt, gene=gene,
                             WT=vals[0], MUT=vals[1], OTH=vals[2], wt_reads=vals[3],
                             wt_altlike_reads=vals[4], mixed_discarded=vals[5], barcode=cell))
    return pd.DataFrame(rows)


def target_id(df):
    # 1-based bam-readcount position as in the target file; indel alt kept verbatim (-GAT / +GAT)
    return df["chr"] + ":" + df["start"].astype(str) + ":" + df["ref"] + ">" + df["alt"]


def call(m, k, p, min_m):
    if m >= min_m and p < P_MUT:
        return "MUT"
    if (m == 0 and k >= 2) or (m > 0 and p >= P_WT):
        return "WT"
    return np.nan  # undetermined (incl. N=0, m=0 & k<2, 0.1<=p<0.9, and m<min_m)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mis-alias", action="store_true")
    ap.add_argument("--outdir")
    args = ap.parse_args()

    c = counts_table()
    c["target_id"] = target_id(c)
    c["type"] = np.where(c["alt"].str[0].isin(["-", "+"]), "INDEL", "SNV")

    # ---------------- per-target error (pooled over cells) ----------------
    t = c.groupby(["target_id", "gene", "type"], sort=False).agg(
        N=("WT", lambda s: 0), w=("OTH", "sum"), wt_reads=("wt_reads", "sum"),
        wt_altlike_reads=("wt_altlike_reads", "sum"), mixed_discarded=("mixed_discarded", "sum"),
    ).reset_index()
    t["N"] = c.groupby("target_id", sort=False)[["WT", "MUT", "OTH"]].sum().sum(axis=1).values
    t["d"] = np.where(t["type"] == "SNV", 2, 1)
    t["q_hat"] = t["w"] / (t["d"] * t["N"]).replace(0, np.nan)
    t["eps"] = EPS
    t["P_s"] = np.maximum(t["q_hat"].fillna(0), EPS)
    t["floor_applied"] = t["q_hat"].fillna(0) < EPS
    t["wt_altlike_frac"] = t["wt_altlike_reads"] / t["wt_reads"].replace(0, np.nan)
    snv = t["type"] == "SNV"
    t.loc[snv, ["wt_reads", "wt_altlike_reads", "wt_altlike_frac", "mixed_discarded"]] = np.nan

    # ---------------- per cell x target calls ----------------
    s = c.merge(t[["target_id", "P_s"]], on="target_id")
    m, k = s["MUT"].values, s["WT"].values
    # upper tail P(X >= m), X ~ Bin(m+k, P_s); P(X>=0)=1
    s["p"] = np.where(m + k > 0, binom.sf(m - 1, m + k, s["P_s"]), np.nan)
    s["genotype"] = [call(mi, ki, pi, 1) for mi, ki, pi in zip(m, k, s["p"])]
    s["genotype_m2"] = [call(mi, ki, pi, 2) for mi, ki, pi in zip(m, k, s["p"])]
    s["het_detect_lb"] = np.where(s["genotype"] == "WT", 1 - 0.5 ** k, np.nan)
    s = s.rename(columns={"barcode": "cell"})
    cols = ["cell", "target_id", "gene", "type", "WT", "MUT", "OTH", "P_s", "p",
            "genotype", "genotype_m2", "het_detect_lb"]
    if args.mis_alias:
        s["MIS"] = s["OTH"]
        cols.insert(7, "MIS")
    sample = s[cols]

    # per-target call summary goes into the targets table too
    g = sample.assign(genotype=sample["genotype"].fillna("undet"))
    summ = g.pivot_table(index="target_id", columns="genotype", values="cell",
                         aggfunc="count", fill_value=0)
    t = t.merge(summ.rename(columns=lambda x: f"n_{x}"), on="target_id", how="left")
    single = g[(g["genotype"] == "MUT") & (g["MUT"] == 1)].groupby("target_id").size()
    t["n_MUT_single_umi"] = t["target_id"].map(single).fillna(0).astype(int)
    t = t[["target_id", "gene", "type", "N", "w", "d", "q_hat", "eps", "P_s", "floor_applied",
           "wt_reads", "wt_altlike_reads", "wt_altlike_frac", "mixed_discarded",
           "n_MUT", "n_WT", "n_undet", "n_MUT_single_umi"]]

    # ---------------- genotype.log ----------------
    tot = s.assign(u=s["WT"] + s["MUT"] + s["OTH"])
    numis = tot.groupby("cell")["u"].sum()
    ngenes = tot.groupby("cell")["u"].apply(lambda x: int(np.sum(x > 0)))
    mut = sample[sample["genotype"] == "MUT"]
    wt = sample[sample["genotype"] == "WT"]
    lb = wt["het_detect_lb"]
    log = [
        "# cells / UMIs (unchanged from main)",
        f"Total number of cells: {numis.size}",
        f"Total number of UMIs: {numis.sum()}",
        f"Median UMI per cell: {np.median(numis)}",
        f"Median genes covered per cell: {np.median(ngenes)}",
        f"Fraction of cells with UMIs: {np.mean(numis > 0):.3f}",
        "",
        "# calls (cell x target)",
        f"MUT: {len(mut)}  WT: {len(wt)}  undetermined: {sample['genotype'].isna().sum()}",
        f"MUT supported by a single UMI (m=1): {int((mut['MUT'] == 1).sum())} "
        f"({(mut['MUT'] == 1).mean():.1%} of MUT)",
        f"MUT with m>=2 required (genotype_m2): {int((sample['genotype_m2'] == 'MUT').sum())}",
        f"WT calls with k=2 / 3 / >=4 WT UMIs: {int((wt['WT'] == 2).sum())} / "
        f"{int((wt['WT'] == 3).sum())} / {int((wt['WT'] >= 4).sum())}",
        f"Het-detection lower bound (1-0.5^k) over WT calls: median {lb.median():.3f}, "
        f"min {lb.min():.3f}",
        "",
        "# per-target error: see <sample>_targets.csv",
        f"Targets with floor applied (q_hat < {EPS:g}): "
        + (", ".join(t.loc[t["floor_applied"], "target_id"]) or "none"),
    ]

    pd.set_option("display.width", 250, "display.max_columns", 30)
    print("=" * 30, "<sample>.csv  (cell x target, long)", "=" * 30)
    print(sample.to_string(index=False))
    print("\n" + "=" * 30, "<sample>_targets.csv  (per target)", "=" * 30)
    print(t.to_string(index=False))
    print("\n" + "=" * 30, "genotype.log", "=" * 30)
    print("\n".join(log))

    if args.outdir:
        os.makedirs(args.outdir, exist_ok=True)
        sample.to_csv(os.path.join(args.outdir, "PROTOTYPE_sample.csv"), index=False)
        t.to_csv(os.path.join(args.outdir, "PROTOTYPE_sample_targets.csv"), index=False)
        with open(os.path.join(args.outdir, "PROTOTYPE_genotype.log"), "w") as f:
            f.write("\n".join(log) + "\n")


if __name__ == "__main__":
    main()
