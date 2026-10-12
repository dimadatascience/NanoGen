# The committed synthetic fixtures (tests/fixtures/synthetic/) carry every
# named scenario they claim to, and the drift check catches any change to them.
# CI runs the drift check itself (make_fixtures.py --check) before pytest.

import functools
import shutil
import subprocess
from collections import defaultdict

import pysam
import pytest

import make_fixtures

DATA = make_fixtures.DEFAULT_OUTDIR
WINDOW = 100  # main's window width, the boundary the straddling molecule crosses


def _levenshtein(a, b):
    # Kept independent of make_fixtures.levenshtein, so the tests do not trust the generator
    prev = list(range(len(b) + 1))
    for i, x in enumerate(a, 1):
        cur = [i]
        for j, y in enumerate(b, 1):
            cur.append(min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (x != y)))
        prev = cur
    return prev[-1]


@functools.lru_cache
def _contig():
    with pysam.FastaFile(str(DATA / "contig.fa")) as fa:
        return fa.fetch(make_fixtures.CONTIG)


def _targets():
    rows = [line.split("\t") for line in (DATA / "targets.bed").read_text().splitlines()]
    return {gene: (chrom, int(start), ref, alt) for chrom, start, _, ref, alt, gene in rows}


def _records():
    """Every record in every cell BAM, as (cell file stem, record)."""
    for bam in sorted((DATA / "cells").glob("*.bam")):
        with pysam.AlignmentFile(str(bam)) as f:
            for read in f:
                yield bam.stem, read


def _scenario(name):
    return [r for _, r in _records() if r.query_name.split(":")[0] == name]


def _primary(reads):
    return [r for r in reads if not r.flag & 0x904]


def _site_base(read, pos1):
    """Base and quality the read carries at 1-based position pos1, or None."""
    for q, ref in read.get_aligned_pairs(matches_only=True):
        if ref == pos1 - 1:
            return read.query_sequence[q], read.query_qualities[q]
    return None


def _indel_at(read, pos1):
    """The indel the read carries at a target at 1-based pos1, as '-SEQ' / '+SEQ'.
    As in bam-readcount, a deletion sits on its first deleted base and an
    insertion on its anchor."""
    ref_pos, q_pos = read.reference_start, 0
    contig = _contig()
    for op, n in read.cigartuples:
        if op == 2 and ref_pos == pos1 - 1:
            return "-" + contig[ref_pos:ref_pos + n]
        if op == 1 and ref_pos == pos1:
            return "+" + read.query_sequence[q_pos:q_pos + n]
        if op in (0, 2):
            ref_pos += n
        if op in (0, 1):
            q_pos += n
    return None


def _umis(reads):
    by_umi = defaultdict(int)
    for r in reads:
        by_umi[r.get_tag("UB")] += 1
    return dict(by_umi)


# --- files -------------------------------------------------------------------


def test_outputs_are_small():
    for path in DATA.rglob("*"):
        if path.is_file():
            assert path.stat().st_size < 64 * 1024, path


def test_bed_has_snvs_an_insertion_a_deletion_and_a_homopolymer_indel():
    targets = _targets()
    contig = _contig()
    snvs = [t for t in targets.values() if not t[3][0] in "+-"]
    assert len(snvs) >= 2
    assert targets["INS"][3].startswith("+")
    assert targets["DEL"][3].startswith("-")
    _, pos, _, alt = targets["HP"]
    run = contig[pos - 1:pos + 4]
    assert alt.startswith("-") and set(alt[1:]) == set(run) and len(set(run)) == 1


def test_bed_ref_matches_the_contig():
    contig = _contig()
    for chrom, pos, ref, _ in _targets().values():
        assert chrom == make_fixtures.CONTIG
        assert contig[pos - 1] == ref


def test_every_read_name_carries_a_known_scenario():
    names = {r.query_name.split(":")[0] for _, r in _records()}
    assert names == {"wt", "ld2_pair", "straddle", "subnmin", "oth", "sec_supp",
                     "del_k", "del_lowq", "ins", "hp_del", "snv_lowq"}
    assert names == set(make_fixtures.SCENARIOS)


def test_cell_bams_are_tagged_and_named_by_cell_barcode():
    for cell, read in _records():
        assert read.get_tag("CB") == cell
        assert len(read.get_tag("UB")) == 12


def test_tagged_bam_holds_every_cell_bam_record():
    # The sample-level BAM as FIXTAGS emits it: SPLIT_BAM on it gives back the cells
    with pysam.AlignmentFile(str(DATA / "tagged" / "S1.bam")) as f:
        tagged = sorted((r.get_tag("CB"), r.to_string()) for r in f)
    cells = sorted((cell, r.to_string()) for cell, r in _records())
    assert tagged == cells
    assert len({cell for cell, _ in tagged}) == 3


def test_distinct_molecules_in_a_cell_have_far_apart_umis():
    by_cell = defaultdict(set)
    for cell, r in _records():
        by_cell[cell].add(r.get_tag("UB"))
    pair = set(_umis(_scenario("ld2_pair")))
    for umis in by_cell.values():
        umis = sorted(umis)
        for i, a in enumerate(umis):
            for b in umis[i + 1:]:
                if {a, b} != pair:
                    assert _levenshtein(a, b) > 4


# --- scenarios ---------------------------------------------------------------


def test_ld2_pair_is_two_umis_at_levenshtein_2_with_10_and_3_reads():
    reads = _scenario("ld2_pair")
    counts = _umis(reads)
    assert sorted(counts.values()) == [3, 10]
    a, b = counts
    assert _levenshtein(a, b) == 2
    assert sum(x != y for x, y in zip(a, b)) > 2  # Hamming alone would not merge them
    assert len({r.get_tag("CB") for r in reads}) == 1


def test_straddle_is_one_molecule_crossing_a_100bp_window():
    reads = _scenario("straddle")
    _, pos, _, _ = _targets()["SNV2"]
    assert len(_umis(reads)) == 1
    windows = {(r.reference_start + 1) // WINDOW for r in reads}  # main bins 1-based POS
    assert len(windows) == 2
    assert all(_site_base(r, pos) for r in reads)


def test_subnmin_cluster_has_fewer_than_n_min_usable_reads():
    _, pos, ref, _ = _targets()["SNV1"]
    reads = _scenario("subnmin")
    assert len(_umis(reads)) == 1
    usable = [r for r in reads if _site_base(r, pos) and _site_base(r, pos)[1] >= 20]
    assert len(usable) == 2  # n_min = 3
    assert {_site_base(r, pos)[0] for r in usable} == {ref}


def test_oth_cluster_carries_neither_ref_nor_alt():
    _, pos, ref, alt = _targets()["SNV1"]
    reads = _scenario("oth")
    bases = {_site_base(r, pos)[0] for r in reads}
    assert len(bases) == 1 and bases.isdisjoint({ref, alt})
    assert len(reads) >= 3


def test_del_k_carries_deletions_within_one_bp_of_the_target_length():
    _, pos, _, alt = _targets()["DEL"]
    reads = _scenario("del_k")
    assert len(_umis(reads)) == 1
    dels = {_indel_at(r, pos) for r in reads} - {None}
    assert alt in dels
    assert {len(d) - len(alt) for d in dels} == {-1, 0, 1}
    # One read deletes the same length 1 bp right of the target: the ±k bundle case
    shifted = "-" + _contig()[pos:pos + len(alt) - 1]
    moved = [r for r in reads if _indel_at(r, pos + 1) == shifted]
    assert len(moved) == 1 and _indel_at(moved[0], pos) is None


def test_del_lowq_deletions_have_a_low_quality_anchor():
    _, pos, _, alt = _targets()["DEL"]
    for r in _scenario("del_lowq"):
        assert _indel_at(r, pos) == alt
        base, qual = _site_base(r, pos - 1)  # the anchor, before the first deleted base
        assert base == _contig()[pos - 2] and qual < 20


def test_snv_lowq_alt_bases_have_low_quality():
    _, pos, _, alt = _targets()["SNV2"]
    reads = _scenario("snv_lowq")
    assert len(_umis(reads)) == 1 and len(reads) >= 3
    assert {_site_base(r, pos)[0] for r in reads} == {alt}
    assert all(_site_base(r, pos)[1] < 20 for r in reads)
    # Alone at SNV2 in its cell, so the SNV base-quality filter leaves it no reads
    assert {r.get_tag("CB") for r in reads}.isdisjoint(
        {r.get_tag("CB") for r in _scenario("straddle") + _scenario("wt")
         if _site_base(r, pos)})


def test_ins_and_hp_del_carry_their_target_indel():
    targets = _targets()
    for scenario, gene in [("ins", "INS"), ("hp_del", "HP")]:
        _, pos, _, alt = targets[gene]
        reads = _scenario(scenario)
        assert reads
        assert {_indel_at(r, pos) for r in reads} == {alt}


def test_sec_supp_has_secondary_and_supplementary_records_of_one_molecule():
    reads = _scenario("sec_supp")
    flags = [r.flag for r in reads]
    assert any(f & 0x100 for f in flags) and any(f & 0x800 for f in flags)
    assert len(_umis(reads)) == 1
    assert len(_primary(reads)) == 2  # sub-n_min once the extra records are dropped
    assert all(r.mapping_quality >= 20 for r in reads)  # only the flag filter removes them


def test_wt_molecules_cover_every_target():
    targets = _targets()
    covered = set()
    for r in _scenario("wt"):
        for gene, (_, pos, ref, _) in targets.items():
            hit = _site_base(r, pos)
            if hit and hit[0] == ref and not _indel_at(r, pos):
                covered.add(gene)
    assert covered == set(targets)


# --- bam-readcount sees the scenarios as intended --------------------------------


def _readcount(tmp_path, cell, min_base_qual):
    if shutil.which("bam-readcount") is None:
        pytest.skip("bam-readcount not on PATH")
    bam = tmp_path / f"{cell}.bam"
    shutil.copy(DATA / "cells" / f"{cell}.bam", bam)
    pysam.index(str(bam))
    out = subprocess.run(
        ["bam-readcount", "-f", str(DATA / "contig.fa"), "-l", str(DATA / "targets.bed"),
         "-b", str(min_base_qual), "-w", "1", str(bam)],
        capture_output=True, text=True, check=True,
    ).stdout
    counts = {}
    for line in out.splitlines():
        fields = line.split("\t")
        counts[int(fields[1])] = {e.split(":")[0]: int(e.split(":")[1]) for e in fields[4:]}
    return counts


def _cell_of(scenario):
    return next(cell for cell, r in _records() if r.query_name.startswith(scenario + ":"))


def test_readcount_reports_the_shifted_deletions_at_the_anchor(tmp_path):
    _, pos, _, alt = _targets()["DEL"]
    counts = _readcount(tmp_path, _cell_of("del_k"), 0)[pos]
    assert counts[alt] == 2 + 4  # del_k exact + del_lowq
    assert counts[alt[:-1]] == 2
    assert counts[alt + _contig()[pos + len(alt) - 2]] == 1
    assert counts[_contig()[pos - 1]] == 4 + 1  # wt + the 1 bp right-shifted read


def test_readcount_base_quality_filter_drops_the_low_quality_anchor(tmp_path):
    _, pos, _, alt = _targets()["DEL"]
    cell = _cell_of("del_lowq")
    assert _readcount(tmp_path, cell, 0)[pos][alt] - _readcount(tmp_path, cell, 20)[pos][alt] == 4


def test_readcount_base_quality_filter_drops_the_low_quality_snv_base(tmp_path):
    _, pos, _, alt = _targets()["SNV2"]
    cell = _cell_of("snv_lowq")
    assert _readcount(tmp_path, cell, 0)[pos][alt] == 4
    assert _readcount(tmp_path, cell, 20)[pos][alt] == 0


def test_readcount_reports_insertion_and_homopolymer_deletion(tmp_path):
    targets = _targets()
    for scenario, gene in [("ins", "INS"), ("hp_del", "HP")]:
        _, pos, _, alt = targets[gene]
        assert _readcount(tmp_path, _cell_of(scenario), 0)[pos][alt] == 4


# --- drift check -------------------------------------------------------------


def test_check_reports_drift(tmp_path):
    copy = tmp_path / "synthetic"
    shutil.copytree(DATA, copy)
    with open(copy / "targets.bed", "a") as f:
        f.write("chrF\t5\t5\tA\tC\tEXTRA\n")
    (copy / "cells" / sorted(p.name for p in (copy / "cells").glob("*.bam"))[0]).unlink()

    diffs = make_fixtures.check(copy)

    assert any("targets.bed" in d for d in diffs)
    assert any(".bam" in d for d in diffs)
