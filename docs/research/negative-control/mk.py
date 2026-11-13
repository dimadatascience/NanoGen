"""Synthetic v3-run output for the smoke tests of sweep.py (#68) and score.py (#67)."""
import numpy as np, pandas as pd, sys, random
sys.path.insert(0, sys.argv[2]); from genotyping import call
rng = np.random.default_rng(1); out = sys.argv[1]
_r = random.Random(1)  # one generator: 3000 distinct barcodes (#68 re-seeded per base, giving 4)
bcs = sorted({"".join(_r.choice("ACGT") for _ in range(16)) for _ in range(3100)})[:3000]
ct = np.array(["Leukemic"] * 1500 + ["T/NK_cells"] * 900 + ["B_cells"] * 300 + ["Ambiguous"] * 300)
agg = np.where(ct == "T/NK_cells", np.where(np.arange(3000) % 5 == 0, "NK_cells", "T_CD4"),
               np.where(ct == "B_cells", "B_mature", "Early_myeloid"))
pd.DataFrame(dict(cell16=bcs, ct=ct, aggregated_ct=agg)).to_csv(f"{out}/meta.csv", index=False)
open(f"{out}/qc.txt", "w").write("\n".join(bcs))
neg = [b for b, c in zip(bcs, ct) if c in ("T/NK_cells", "B_cells")]
negs = set(neg)
pd.DataFrame(dict(sample="X", cell=neg)).to_csv(f"{out}/neg.csv", index=False)
rows = []; trows = []
for t, (g, cont, vaf) in enumerate([("GA", 0.01, 0.5), ("GB", 0.03, 0.4)]):
    tid = f"chr1:{t}:A>T"; Ps = 1e-4
    sub = []
    for b, c in zip(bcs, ct):
        n = rng.poisson(3)
        rate = vaf if c == "Leukemic" else cont
        m = rng.binomial(n, rate); sub.append((b, c, m, n - m))
    nm = sum(m for b, c, m, k in sub if b in negs); nn = sum(m + k for b, c, m, k in sub if b in negs)
    cs = nm / nn; Pn = Ps + cs
    M = np.array([s[2] for s in sub]); K = np.array([s[3] for s in sub])
    p, cl = call(M, K, np.full(len(M), Pn), 1)
    for (b, c, m, k), pp, cc in zip(sub, p, cl):
        rows.append(dict(cell=f"X_{b}-1", target_id=tid, gene=g, type="SNV", WT=k, MUT=m, OTH=0, P_s=Ps, c_s=cs,
                         P_null=Pn, p=pp, genotype=cc, het_detect_lb=np.nan, negative_control=b in negs))
    trows.append(dict(target_id=tid, gene=g, P_s=Ps, c_s=cs))
pd.DataFrame(rows).to_csv(f"{out}/res/X.csv", index=False)
pd.DataFrame(trows).to_csv(f"{out}/res/X_targets.csv", index=False)
open(f"{out}/key.txt", "w").write("GA\nGB\n")
