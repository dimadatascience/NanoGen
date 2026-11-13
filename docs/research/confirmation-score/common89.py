"""Shared inputs and fixed sets for scoring the confirmation pre-registration
(docs/research/confirmation-preregistration.md @ e6f629f), ticket #89, map #72.

score_lib.py is the frozen P64 scoring library (negative-control-palm22e-palm09e, unmodified); `call` is imported
from bin/long_reads/genotyping.py at genotype-v3@61d8809 (GENOTYPING_DIR). Aggregate output only: targets T1..Tn.
"""
import csv
import os
import re

import numpy as np
import pandas as pd
from scipy.stats import beta

import score_lib as L

W = "/beegfs/scratch/ieo6129/scmseq"
RUNS = f"{W}/q79_full/runs"
PATIENTS = ["PALM07_E", "PALM08_E", "PALM11_E", "PALM14_E", "PALM21_E", "AML4", "AML5", "sAML1"]

# §4: 20 new seeds = random.Random(20261112).sample(range(1, 2**31), 20); re-derived and asserted below.
SEEDS = [452483223, 1695653105, 1085084613, 441224510, 821422349,
         1758734097, 1576893678, 1902659721, 536475717, 2055609293,
         1297200436, 2099741074, 248750689, 1900119910, 1632811357,
         522860322, 934946355, 760461626, 377745768, 1704488278]


def _rng(a, b):
    return [f"T{i}" for i in range(a, b + 1)]


# §5: the 57 pool targets (pool86.tsv @ 546386e, in_pool = True); asserted against the table.
POOL = {
    "PALM07_E": _rng(1, 5),
    "PALM08_E": ["T1"] + _rng(3, 7),
    "PALM11_E": ["T1"] + _rng(3, 5),
    "PALM14_E": _rng(1, 7),
    "PALM21_E": _rng(1, 5) + _rng(7, 12),
    "AML4": _rng(1, 4) + ["T6", "T9", "T10"],
    "AML5": _rng(1, 4) + ["T7", "T10", "T11", "T16"],
    "sAML1": ["T1", "T2", "T4"] + _rng(7, 12),
}
# §6.1: tumour-evaluable targets (17 targets = 15 loci in 5 patients); phased locus PALM14E {T3-T5}.
TUMOUR = {
    "PALM08_E": ["T3", "T4"],
    "PALM11_E": ["T1"],
    "PALM14_E": ["T1", "T2", "T3", "T4", "T5", "T6", "T7"],
    "AML5": ["T1", "T3"],
    "sAML1": ["T7", "T9", "T10", "T11", "T12"],
}
PHASED = ("PALM14_E", ("T3", "T4", "T5"))
# WES-clonal sets (#77, #87), for listing under-covered clonal targets (§6.1).
CLONAL = {
    "PALM08_E": _rng(1, 4), "PALM11_E": ["T1"], "PALM14_E": _rng(1, 7),
    "AML4": _rng(1, 4) + _rng(6, 10), "AML5": _rng(1, 3) + _rng(5, 7) + ["T9", "T10", "T12"],
    "sAML1": _rng(2, 5) + ["T7"] + _rng(9, 14),
}
ALPHA_HOLM = 0.10
FRAG_MIN_NA = 11

BARE = re.compile(r"[ACGTN]{16}")


def bare(c):
    m = BARE.search(c)
    return m.group(0) if m else None


def read_list(path):  # verbatim from q76 power_check.py @ 9a2d94f
    out = set()
    with open(path) as f:
        rows = list(csv.reader(f))
    if not rows:
        return out
    head = [h.strip().strip('"') for h in rows[0]]
    col = next((head.index(k) for k in ("cell", "cell16") if k in head), 0)
    for r in rows[1:] if any(k in head for k in ("cell", "cell16")) else rows:
        if len(r) > col and bare(r[col]):
            out.add(bare(r[col]))
    return out


def read_labels(path):  # verbatim from q76 power_check.py @ 9a2d94f
    lab = {}
    with open(path) as f:
        rows = list(csv.reader(f))
    head = [h.strip().strip('"') for h in rows[0]]
    if "ct" in head:
        ic, it = head.index("cell16"), head.index("ct")
        body = rows[1:]
    else:
        ic, it, body = 0, 1, rows
    for r in body:
        b = bare(r[ic])
        if b:
            lab[b] = r[it]
    return lab


def label_path(S):
    if S == "PALM07_E":
        return f"{W}/q88_lib2/derived/PALM07_E_lib2_ct.csv"
    if S in ("AML4", "AML5", "sAML1"):
        return f"{W}/q84_restart/derived/{S}_ct.csv"
    return f"{W}/q74_intake/out/annot_{S}.csv"


def cpu(x, n):
    """§3: CPU(x, n) = beta.ppf(0.95, x + 1, n − x); 1 if x = n or n = 0 (pre-registration definition)."""
    if n == 0 or x >= n:
        return 1.0
    return float(beta.ppf(0.95, x + 1, n - x))


def bed_keys(bed):
    """BED order -> T1..Tn; returns list of (chrom, start, end, ref, alt)."""
    out = []
    for line in open(bed):
        f = line.rstrip("\r\n").split("\t")
        if len(f) < 5 or not f[0] or f[0].startswith(("#", "track")):
            continue
        out.append((f[0], int(f[1]), int(f[2]), f[3], f[4]))
    return out


def load(S):
    """Patient S: v3-run rows restricted to QC cells, with T labels, n, and cell-set membership."""
    inputs = pd.read_csv(f"{W}/q79_full/sheets/inputs.tsv", sep="\t").set_index("sample")
    d = pd.read_csv(f"{RUNS}/v3_{S}/results/{S}/{S}.csv")
    tg = pd.read_csv(f"{RUNS}/v3_{S}/results/{S}/{S}_targets.csv")
    keys = bed_keys(inputs.loc[S, "bed"])
    t_of = {}
    for i, (c, s, e, r, a) in enumerate(keys):
        cand = {f"{c}:{p}:{r}>{a}" for p in (s, s + 1, e)}
        hit = [t for t in tg["target_id"] if t in cand]
        assert len(hit) == 1, f"{S}: BED row {i + 1} matches {len(hit)} target ids"
        t_of[hit[0]] = f"T{i + 1}"
    assert len(t_of) == len(tg) == len(set(t_of.values())), f"{S}: target map not a bijection"
    qc = read_list(inputs.loc[S, "qc"])
    neg = read_list(inputs.loc[S, "negatives"]) & qc
    lab = read_labels(label_path(S))
    leuk = {b for b, c in lab.items() if c == "Leukemic"} & qc
    amb = {b for b, c in lab.items() if c == "Ambiguous"} & qc
    d["bc"] = d["cell"].map(bare)
    d = d[d["bc"].isin(qc)].copy()
    d["n"] = d["WT"].astype(int) + d["MUT"].astype(int)
    d["T"] = d["target_id"].map(t_of)
    d["set"] = np.where(d["bc"].isin(neg), "neg", np.where(d["bc"].isin(leuk), "leuk",
                                                           np.where(d["bc"].isin(amb), "amb", "other")))
    tg["T"] = tg["target_id"].map(t_of)
    meta = dict(qc=len(qc), neg=len(neg), leuk=len(leuk), amb=len(amb), neg_leuk_overlap=len(neg & leuk),
                cells_in_run=int(d["bc"].nunique()))
    return d, tg.set_index("T"), sorted(neg), meta


def target_stats(d, tg, T):
    """In-sample per-target quantities from the v3 run (§5.4, §10)."""
    sub = d[d["T"] == T]
    ng, lk = sub[sub["set"] == "neg"], sub[(sub["set"] == "leuk") & (sub["n"] >= 1)]
    c_s, x_s, n_s = L.pooled_rate(ng["MUT"].values, ng["WT"].values)
    q_T, _, _ = L.pooled_rate(lk["MUT"].values, lk["WT"].values)
    flag = (c_s >= 0.25 or c_s >= 0.5 * q_T) and not (c_s == 0 and q_T == 0)
    P_s = float(tg.loc[T, "P_s"])
    return dict(P_s=P_s, c_s=c_s, x_s=x_s, n_s=n_s, cpu_s=cpu(x_s, n_s), q_T=q_T, flag=bool(flag),
                flag_q_clause_only=bool(flag and c_s < 0.25), c_s_run=float(tg.loc[T, "c_s"]),
                leuk_covered=int(len(lk)), neg_covered=int((ng["n"] >= 1).sum()))


def check_fixed_sets():
    import random
    assert random.Random(20261112).sample(range(1, 2 ** 31), 20) == SEEDS, "seed list mismatch"
    here = os.path.dirname(os.path.abspath(__file__))
    p = pd.read_csv(f"{here}/pool86.tsv", sep="\t")
    p = p[p["in_pool"].astype(str) == "True"]
    got = {S: sorted(p.loc[p["sample"] == S, "target"], key=lambda t: int(t[1:])) for S in PATIENTS}
    assert got == POOL and sum(map(len, POOL.values())) == 57, f"pool mismatch {got}"
    pl = pd.read_csv(f"{here}/pool86.tsv", sep="\t")
    ev = {}
    for S, ts in CLONAL.items():
        for T in ts:
            r = pl[(pl["sample"] == S) & (pl["target"] == T)]
            if len(r) and int(r["leuk_covered"].iloc[0]) >= 100:
                ev.setdefault(S, []).append(T)
    assert ev == TUMOUR, f"tumour-evaluable mismatch {ev}"
