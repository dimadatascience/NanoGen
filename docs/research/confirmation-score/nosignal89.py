#!/usr/bin/env python3
"""Stage 1 of #89 (pre-registration @ e6f629f §3 checks, §5.4 flag, §6.2 no-signal list).

Computes, from the v3 run's counts only and before any contender's tumour outcome:
  - the implementation checks (abort on failure);
  - the outlier flag and in-sample contamination for every target;
  - the 3b mixture fit (π̂, f̂) on every tumour-evaluable, unflagged target, and the no-signal list.
No calls are made on tumour cells here. Output: out/stage1_flags.tsv, out/stage1_nosignal.tsv, out/stage1_checks.json.
Usage: nosignal89.py OUTDIR
"""
import json
import sys

import numpy as np
import pandas as pd

import common89 as C
import score_lib as L


def main():
    out = sys.argv[1]
    C.check_fixed_sets()
    checks, flags, ns = {}, [], []
    for S in C.PATIENTS:
        d, tg, neg, meta = C.load(S)
        full = pd.read_csv(f"{C.RUNS}/v3_{S}/results/{S}/{S}.csv")
        _, rc = L.call(full["MUT"].values, full["WT"].values, full["P_null"].values, 1)
        same = bool((pd.Series(rc).fillna("NA").values == full["genotype"].fillna("NA").values).all())
        max_dc = 0.0
        for T in tg.index:
            st = C.target_stats(d, tg, T)
            max_dc = max(max_dc, abs(st["c_s"] - st["c_s_run"]))
            row = dict(sample=S, target=T, in_pool=T in C.POOL[S], tumour_evaluable=T in C.TUMOUR.get(S, []),
                       **{k: v for k, v in st.items() if k != "c_s_run"})
            flags.append(row)
            if T in C.TUMOUR.get(S, []):
                lk = d[(d["T"] == T) & (d["set"] == "leuk") & (d["n"] >= 1)]
                rho0 = min(1.0, st["P_s"] + st["c_s"])
                if st["flag"]:
                    ns.append(dict(sample=S, target=T, flagged=True, rho0=rho0, pi_hat=np.nan, f_hat=np.nan,
                                   leuk_covered=len(lk), no_signal=None))
                    continue
                pi, f = L.em_mixture(lk["MUT"].values, lk["n"].values, rho0)
                ns.append(dict(sample=S, target=T, flagged=False, rho0=rho0, pi_hat=pi, f_hat=f,
                               leuk_covered=len(lk), no_signal=bool(pi < 0.05 or f < rho0 + 0.05)))
        checks[S] = dict(recall_identical=same, max_abs_c_s_diff=max_dc, **meta)
        if not same or max_dc > 1e-12:
            json.dump(checks, open(f"{out}/stage1_checks.json", "w"), indent=1)
            raise SystemExit(f"implementation check failed for {S}: {checks[S]}")
    checks["numpy"], checks["pandas"] = np.__version__, pd.__version__
    import scipy
    checks["scipy"] = scipy.__version__
    json.dump(checks, open(f"{out}/stage1_checks.json", "w"), indent=1)
    pd.DataFrame(flags).to_csv(f"{out}/stage1_flags.tsv", sep="\t", index=False)
    pd.DataFrame(ns).to_csv(f"{out}/stage1_nosignal.tsv", sep="\t", index=False)
    print(json.dumps(checks, indent=1))
    print(pd.DataFrame(ns).to_string())
    fl = pd.DataFrame(flags)
    print(fl[fl["flag"]][["sample", "target", "in_pool", "tumour_evaluable", "c_s", "q_T", "flag_q_clause_only",
                          "n_s", "leuk_covered"]].to_string())


if __name__ == "__main__":
    main()
