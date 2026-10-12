# Writes the tiny reference, targets and per-cell SAMs for tests/consensus_tail.nf.test.
# Usage: python make_data.py tests/data/consensus_tail
import random
import sys

out = sys.argv[1]
random.seed(28)
L = 80
ref = "".join(random.choice("ACGT") for _ in range(L))
SNV, DEL1, DEL3 = 20, 40, 55  # 1-based first altered base


def other(b, k=0):
    return "ACGT".replace(b, "")[k % 3]


open(f"{out}/ref.fa", "w").write(f">chrT\n{ref}\n")
# bam-readcount reports a deletion on its first deleted base
bed = [
    ("chrT", SNV, SNV, ref[SNV - 1], other(ref[SNV - 1]), "SNV"),
    ("chrT", DEL1, DEL1, ref[DEL1 - 1], "-" + ref[DEL1 - 1], "DEL1"),
    ("chrT", DEL3, DEL3, ref[DEL3 - 1], "-" + ref[DEL3 - 1 : DEL3 + 2], "DEL3"),
]
open(f"{out}/targets.bed", "w").write(
    "".join("\t".join(map(str, r)) + "\n" for r in bed)
)


def read(name, cb, ub, snv=None, del1=False, del3=False, lowq_snv=False):
    seq, cig, qual = [], [], []
    i = 1
    while i <= L:
        if del1 and i == DEL1:
            cig.append("D")
            i += 1
            continue
        if del3 and i == DEL3:
            cig += ["D"] * 3
            i += 3
            continue
        b = ref[i - 1]
        if i == SNV and snv is not None:
            b = snv
        seq.append(b)
        cig.append("M")
        qual.append("+" if (i == SNV and lowq_snv) else "I")
        i += 1
    c, prev, n = "", cig[0], 0
    for x in cig:
        if x == prev:
            n += 1
        else:
            c += f"{n}{prev}"
            prev, n = x, 1
    c += f"{n}{prev}"
    return (
        f"{name}\t0\tchrT\t1\t60\t{c}\t*\t0\t0\t{''.join(seq)}\t{''.join(qual)}"
        f"\tCB:Z:{cb}\tUB:Z:{ub}\n"
    )


hdr = f"@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chrT\tLN:{L}\n"
alt = bed[0][4]
oth = other(ref[SNV - 1], 1)
cells = {
    "AAACCCAAGAAACACT": [
        *[("mut", "ACGTACGTACGT", dict(snv=alt, del1=True))] * 5,
        ("mut_edit1", "ACGTACGTACGA", dict(snv=alt, del1=True)),
        *[("wt", "TTTTGGGGCCCC", {})] * 4,
        *[("longdel", "GGGGAAAATTTT", dict(del3=True))] * 4,
    ],
    "AAACCCAAGAAACCAT": [
        *[("oth", "CCCCAAAAGGGG", dict(snv=oth))] * 3,
        *[("wt", "GAGAGAGATCTC", {})] * 3,
        ("wt_lowq", "GAGAGAGATCTC", dict(snv=alt, lowq_snv=True)),
        *[("small", "ATATATATCGCG", dict(snv=alt))] * 1,
        # Shortcut MUT where the general rule would say WT:
        # 3 of 9 reads carry the 1 bp deletion, 4 of 14 the 3 bp deletion
        *[("shortcut_del1", "TGCATGCATGCA", dict(del1=True))] * 3,
        *[("shortcut_del1_ref", "TGCATGCATGCA", {})] * 6,
        *[("shortcut_del3", "CATGCATGCATG", dict(del3=True))] * 4,
        *[("shortcut_del3_ref", "CATGCATGCATG", {})] * 10,
    ],
}
for cb, reads in cells.items():
    with open(f"{out}/{cb}.sam", "w") as f:
        f.write(hdr)
        for k, (name, ub, kw) in enumerate(reads):
            f.write(read(f"{name}_{k}", cb, ub, **kw))
