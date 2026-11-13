#!/usr/bin/env python3
"""Ticket #91, §10 indel-arm diagnostics: tol_reads, tol_unseen_reads, brc_mismatch_umis summed per indel target and
cell set, from the prototype run's per-cell CONSENSUS tables (successful tasks only: .exitcode == 0; one table per
barcode). Aggregate output only. Usage: diag91.py OUTDIR
"""
import glob
import os
import sys

import pandas as pd

import common89 as C
from score91 import load_arm

COLS = ["tol_reads", "tol_unseen_reads", "brc_mismatch_umis", "WT", "MUT"]


def main():
    out = sys.argv[1]
    rows, meta = [], []
    for S in C.PATIENTS:
        d, tg, neg = load_arm(S, "proto")
        indel = {tid: T for T, tid in tg["target_id"].items() if tg.loc[T, "type"] == "INDEL"}
        if not indel:
            continue
        sets = d.drop_duplicates("bc").set_index("bc")["set"]
        seen, dup, bad = {}, 0, 0
        for f in glob.glob(f"{C.RUNS}/proto_{S}/work/*/*/*_table.csv"):
            wd = os.path.dirname(f)
            try:
                ok = open(f"{wd}/.exitcode").read().strip() == "0"
            except OSError:
                ok = False
            if not ok:
                bad += 1
                continue
            bc = os.path.basename(f)[:-len("_table.csv")]
            if bc in seen:
                dup += 1
                continue
            seen[bc] = f
        parts = []
        for bc, f in seen.items():
            if bc not in sets.index:
                continue
            t = pd.read_csv(f)
            t["tid"] = t["chr"].astype(str) + ":" + t["start"].astype(str) + ":" + t["ref"] + ">" + t["alt"]
            t = t[t["tid"].isin(indel)]
            t["set"] = sets[bc]
            parts.append(t[["tid", "set"] + COLS])
        a = pd.concat(parts)
        a["T"] = a["tid"].map(indel)
        g = a.groupby(["T", "set"])[COLS].sum().reset_index()
        g["cells_tol"] = a[a["tol_reads"] > 0].groupby(["T", "set"]).size().reindex(
            pd.MultiIndex.from_frame(g[["T", "set"]])).fillna(0).astype(int).values
        g.insert(0, "sample", S)
        rows.append(g)
        meta.append(dict(sample=S, qc_cells_with_table=int(sum(1 for b in seen if b in sets.index)),
                         qc_cells_in_output=int(len(sets)), duplicate_ok_tables=dup, failed_task_tables=bad))
    pd.concat(rows).to_csv(f"{out}/diag_by_set.tsv", sep="\t", index=False)
    pd.DataFrame(meta).to_csv(f"{out}/diag_meta.tsv", sep="\t", index=False)
    tot = pd.concat(rows).groupby(["sample", "T"])[COLS + ["cells_tol"]].sum()
    tot.to_csv(f"{out}/diag.tsv", sep="\t")
    pd.set_option("display.width", 250)
    print(pd.DataFrame(meta).to_string())
    print(tot.to_string())


if __name__ == "__main__":
    main()
