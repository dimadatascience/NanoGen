#!/usr/bin/env python
"""
Seeded synthetic fixtures for NanoGen tests. No real or patient data.

Writes, under tests/fixtures/synthetic/:
    contig.fa, contig.fa.fai  one random contig with planted indel sites
    targets.bed               SNVs, a deletion, an insertion, a homopolymer indel
    cells/<CB>.bam            per-cell BAMs tagged CB/UB, as SPLIT_BAM emits them

Every read is named <scenario>:<molecule>:<read>, so a test can pick out the
records of one scenario. Reads match the contig exactly except for the planted
allele; a molecule is one UB.

    python tests/fixtures/make_fixtures.py           # regenerate
    python tests/fixtures/make_fixtures.py --check   # fail if committed outputs drift

The check compares `samtools view -h` text, not BAM bytes.
"""

import argparse
import difflib
from array import array
import random
import sys
import tempfile
from pathlib import Path

import pysam

SEED = 27
CONTIG = "chrF"
CONTIG_LEN = 2000
READ_LEN = 150
UMI_LEN = 12
MIN_UMI_DIST = 5  # Levenshtein distance between unrelated molecules in a cell
QUAL = 30
LOW_QUAL = 5
MAPQ = 60
SECONDARY_OR_SUPPLEMENTARY = 0x900
DEFAULT_OUTDIR = Path(__file__).resolve().parent / "synthetic"

# Planted sequence (0-based start -> bases), so indel sites have a single,
# left-aligned representation: DEL anchor G + TTACA, INS anchor A + T,
# HP anchor G + a 7 bp T homopolymer.
PLANTS = {1000: "GTTACA", 1300: "AT", 1600: "GTTTTTTTC"}

# gene -> (1-based position, alt); ref is read off the contig. Positions follow
# bam-readcount: a deletion sits on its first deleted base, an insertion on its
# anchor. SNV alts are drawn from the seed.
INDEL_TARGETS = {"DEL": (1002, "-TTA"), "INS": (1301, "+GTC"), "HP": (1602, "-T")}
SNV_POSITIONS = {"SNV1": 301, "SNV2": 701}

CELL_SNV = "AAACCCAAGAAACACT"
CELL_MIX = "AAACCCAAGAAACCAT"
CELL_DEL = "AAACCCAAGAAACCCA"

SCENARIOS = (
    "wt",  # a REF molecule at every target
    "ld2_pair",  # UMIs at Levenshtein 2 (Hamming > 2), 10:3 reads, both ALT at SNV1
    "straddle",  # one ALT molecule at SNV2 whose reads start either side of a 100 bp window
    "subnmin",  # a REF molecule at SNV1 with 2 reads (< n_min = 3)
    "oth",  # a molecule at SNV1 carrying neither REF nor ALT
    "sec_supp",  # 2 primary ALT reads at SNV1 plus a secondary and a supplementary record
    "del_k",  # one DEL molecule: 2 exact, 2 one base shorter, 1 one base longer, 1 starting 1 bp right
    "del_lowq",  # an exact DEL molecule whose anchor base has quality 5
    "ins",  # an INS molecule
    "hp_del",  # a homopolymer 1 bp deletion molecule
)


##


def levenshtein(a, b):
    prev = list(range(len(b) + 1))
    for i, x in enumerate(a, 1):
        cur = [i]
        for j, y in enumerate(b, 1):
            cur.append(min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (x != y)))
        prev = cur
    return prev[-1]


def contig_sequence(rng):
    seq = [rng.choice("ACGT") for _ in range(CONTIG_LEN)]
    for start, bases in PLANTS.items():
        seq[start:start + len(bases)] = bases
    return "".join(seq)


def targets(contig, rng):
    """BED rows (chr, start, end, ref, alt, gene), sorted by position, plus the OTH base per SNV."""
    rows, oth = [], {}
    for gene, pos in SNV_POSITIONS.items():
        ref = contig[pos - 1]
        alt, oth[gene] = rng.sample([b for b in "ACGT" if b != ref], 2)
        rows.append((CONTIG, pos, pos, ref, alt, gene))
    for gene, (pos, alt) in INDEL_TARGETS.items():
        rows.append((CONTIG, pos, pos, contig[pos - 1], alt, gene))
    return sorted(rows, key=lambda r: r[1]), oth


class Umis:
    """Seeded UMIs, at least MIN_UMI_DIST from every other UMI of the same cell."""

    def __init__(self, rng):
        self.rng = rng
        self.taken = {}

    def new(self, cell):
        taken = self.taken.setdefault(cell, [])
        while True:
            umi = "".join(self.rng.choice("ACGT") for _ in range(UMI_LEN))
            if all(levenshtein(umi, u) >= MIN_UMI_DIST for u in taken):
                taken.append(umi)
                return umi

    def ld2_partner(self, cell, umi):
        """umi with one base deleted and one appended: Levenshtein 2, Hamming > 2."""
        for i in range(2, UMI_LEN - 2):
            for b in "ACGT":
                partner = umi[:i] + umi[i + 1:] + b
                hamming = sum(x != y for x, y in zip(umi, partner))
                others = [u for u in self.taken[cell] if u != umi]
                if (levenshtein(umi, partner) == 2 and hamming > 2
                        and all(levenshtein(partner, u) >= MIN_UMI_DIST for u in others)):
                    self.taken[cell].append(partner)
                    return partner
        raise ValueError(f"no LD-2 partner for {umi}")


def make_read(header, contig, name, cell, umi, start, snv=None, indel=None,
              anchor_qual=None, flag=0):
    """A read over contig[start:start + READ_LEN], with an SNV (pos0, base) or an
    indel (anchor0, '-SEQ' / '+SEQ') planted after the 0-based anchor."""
    seq = contig[start:start + READ_LEN]
    cigar = [(0, READ_LEN)]
    if snv:
        pos, base = snv
        seq = seq[:pos - start] + base + seq[pos - start + 1:]
    if indel:
        anchor, alt = indel
        n, size = anchor - start + 1, len(alt) - 1
        if alt[0] == "-":
            assert contig[anchor + 1:anchor + 1 + size] == alt[1:]
            seq = seq[:n] + seq[n + size:]
            cigar = [(0, n), (2, size), (0, READ_LEN - n - size)]
        else:
            seq = seq[:n] + alt[1:] + seq[n:]
            cigar = [(0, n), (1, size), (0, READ_LEN - n)]
    quals = [QUAL] * len(seq)
    if anchor_qual is not None:
        quals[indel[0] - start] = anchor_qual

    read = pysam.AlignedSegment(header)
    read.query_name = name
    read.flag = flag
    read.reference_id = 0
    read.reference_start = start
    read.mapping_quality = MAPQ
    read.cigartuples = cigar
    read.query_sequence = seq
    read.query_qualities = array("B", quals)
    read.set_tags([("CB", cell), ("UB", umi)])
    return read


def molecules(contig, bed, oth, umis):
    """Yield (cell, scenario, molecule, umi, reads); each read is a make_read kwargs dict."""
    site = {gene: (pos - 1, ref, alt) for _, pos, _, ref, alt, gene in bed}

    def starts(gene, n, step=3):
        return [site[gene][0] - 60 + step * j for j in range(n)]

    def snv_reads(gene, base, n):
        return [dict(start=s, snv=(site[gene][0], base)) for s in starts(gene, n)]

    def indel_reads(gene, alts, **kw):
        anchor = anchor_of(gene)
        return [dict(start=s, indel=(anchor, a), **kw)
                for s, a in zip(starts(gene, len(alts)), alts)]

    def anchor_of(gene):
        pos, _, alt = site[gene]
        return pos - 1 if alt[0] == "-" else pos

    def ref_reads(gene, n):
        return [dict(start=s) for s in starts(gene, n)]

    # Cell with the SNV scenarios at SNV1
    _, ref1, alt1 = site["SNV1"]
    yield CELL_SNV, "wt", "snv1", umis.new(CELL_SNV), ref_reads("SNV1", 4)
    major = umis.new(CELL_SNV)
    yield CELL_SNV, "ld2_pair", "major", major, snv_reads("SNV1", alt1, 10)
    minor = umis.ld2_partner(CELL_SNV, major)
    yield CELL_SNV, "ld2_pair", "minor", minor, snv_reads("SNV1", alt1, 3)
    yield CELL_SNV, "subnmin", "snv1", umis.new(CELL_SNV), ref_reads("SNV1", 2)
    yield CELL_SNV, "oth", "snv1", umis.new(CELL_SNV), snv_reads("SNV1", oth["SNV1"], 4)
    primary = snv_reads("SNV1", alt1, 2)
    # Full copies of the primaries at MAPQ 60, so only the flag filter removes them
    extra = [dict(r, flag=f) for r, f in zip(snv_reads("SNV1", alt1, 2), (0x100, 0x800))]
    yield CELL_SNV, "sec_supp", "snv1", umis.new(CELL_SNV), primary + extra

    # Cell with the window-straddling molecule, the insertion and the homopolymer
    pos2, _, alt2 = site["SNV2"]
    straddle = [dict(start=s, snv=(pos2, alt2)) for s in (590, 594, 598, 600, 604, 608)]
    yield CELL_MIX, "straddle", "snv2", umis.new(CELL_MIX), straddle
    yield CELL_MIX, "wt", "snv2", umis.new(CELL_MIX), ref_reads("SNV2", 4)
    yield CELL_MIX, "ins", "ins", umis.new(CELL_MIX), indel_reads("INS", [site["INS"][2]] * 4)
    yield CELL_MIX, "wt", "ins", umis.new(CELL_MIX), ref_reads("INS", 4)
    yield CELL_MIX, "hp_del", "hp", umis.new(CELL_MIX), indel_reads("HP", [site["HP"][2]] * 4)
    yield CELL_MIX, "wt", "hp", umis.new(CELL_MIX), ref_reads("HP", 4)

    # Cell with the deletion scenarios
    anchor, alt = anchor_of("DEL"), site["DEL"][2]
    longer = alt + contig[anchor + len(alt)]
    shifted = "-" + contig[anchor + 2:anchor + 1 + len(alt)]  # same length, 1 bp right
    yield CELL_DEL, "wt", "del", umis.new(CELL_DEL), ref_reads("DEL", 4)
    del_k = indel_reads("DEL", [alt, alt, alt[:-1], alt[:-1], longer])
    del_k.append(dict(start=del_k[-1]["start"] + 3, indel=(anchor + 1, shifted)))
    yield CELL_DEL, "del_k", "del", umis.new(CELL_DEL), del_k
    yield CELL_DEL, "del_lowq", "del", umis.new(CELL_DEL), indel_reads(
        "DEL", [alt] * 4, anchor_qual=LOW_QUAL)


def build(outdir):
    """Write every fixture file under outdir."""
    rng = random.Random(SEED)
    contig = contig_sequence(rng)
    bed, oth = targets(contig, rng)
    outdir = Path(outdir)
    cells_dir = outdir / "cells"
    cells_dir.mkdir(parents=True, exist_ok=True)
    for stale in cells_dir.glob("*.bam"):
        stale.unlink()

    fasta = outdir / "contig.fa"
    lines = [contig[i:i + 60] for i in range(0, CONTIG_LEN, 60)]
    fasta.write_text(f">{CONTIG}\n" + "\n".join(lines) + "\n")
    pysam.faidx(str(fasta))
    (outdir / "targets.bed").write_text("".join("\t".join(map(str, r)) + "\n" for r in bed))

    header = pysam.AlignmentHeader.from_dict(
        {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": CONTIG, "LN": CONTIG_LEN}]})
    by_cell = {}
    for cell, scenario, molecule, umi, reads in molecules(contig, bed, oth, Umis(rng)):
        assert scenario in SCENARIOS
        # Secondary and supplementary records reuse the names of the primary reads
        index = {0: 0, 1: 0}
        for spec in reads:
            extra = bool(spec.get("flag", 0) & SECONDARY_OR_SUPPLEMENTARY)
            index[extra] += 1
            name = f"{scenario}:{molecule}:{index[extra]:02d}"
            by_cell.setdefault(cell, []).append(make_read(header, contig, name, cell, umi, **spec))

    for cell, reads in by_cell.items():
        reads.sort(key=lambda r: (r.reference_start, r.query_name, r.flag))
        with pysam.AlignmentFile(str(cells_dir / f"{cell}.bam"), "wb", header=header) as f:
            for read in reads:
                f.write(read)


def _text(path):
    if path.suffix == ".bam":
        return pysam.view("-h", "--no-PG", str(path))
    return path.read_text()


def check(committed):
    """Regenerate into a temporary directory; return one message per drifted file."""
    committed = Path(committed)
    with tempfile.TemporaryDirectory() as tmp:
        fresh = Path(tmp)
        build(fresh)
        files = lambda d: {p.relative_to(d) for p in d.rglob("*")
                           if p.is_file() and not p.name.startswith(".")}
        ours, theirs = files(fresh), files(committed)
        diffs = [f"{p}: missing from {committed}" for p in sorted(ours - theirs)]
        diffs += [f"{p}: not produced by the generator" for p in sorted(theirs - ours)]
        for p in sorted(ours & theirs):
            a, b = _text(committed / p), _text(fresh / p)
            if a != b:
                diff = difflib.unified_diff(a.splitlines(), b.splitlines(), f"committed/{p}",
                                            f"regenerated/{p}", lineterm="", n=1)
                diffs.append(f"{p}: differs\n" + "\n".join(list(diff)[:40]))
    return diffs


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--check", action="store_true",
                        help="Regenerate into a temporary directory and fail if outdir differs.")
    args = parser.parse_args()

    if not args.check:
        build(args.outdir)
        return 0
    diffs = check(args.outdir)
    for d in diffs:
        print(d, file=sys.stderr)
    if diffs:
        print("Fixtures drifted: rerun tests/fixtures/make_fixtures.py and commit.", file=sys.stderr)
        return 1
    print(f"Fixtures in {args.outdir} match the generator.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
