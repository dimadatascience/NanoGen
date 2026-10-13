# Spec tests for the UMI grouping rule (#38): banded Levenshtein and the
# hybrid merge within one bundle.

import random
from pathlib import Path

import pysam
import pytest

import umi_grouping as ug

REPO = Path(__file__).resolve().parents[2]

THRESHOLDS = dict(max_ld=2, ratio=2, small_reads=3)


def cluster(counts, **overrides):
    return ug.cluster_umis(counts, **{**THRESHOLDS, **overrides})


def same_cluster(clusters, a, b):
    return clusters[a] == clusters[b]


def full_levenshtein(a, b):
    """Unbanded reference DP."""
    prev = list(range(len(b) + 1))
    for i, ca in enumerate(a, 1):
        cur = [i]
        for j, cb in enumerate(b, 1):
            cur.append(min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (ca != cb)))
        prev = cur
    return prev[-1]


# Banded Levenshtein


def test_levenshtein_identical_is_zero():
    assert ug.levenshtein("ACGTACGTACGT", "ACGTACGTACGT", 2) == 0


def test_levenshtein_counts_a_substitution():
    assert ug.levenshtein("ACGTACGTACGT", "ACGTACCTACGT", 2) == 1


def test_levenshtein_counts_an_insertion_and_a_deletion():
    assert ug.levenshtein("ACGTACGTACGT", "ACGTAACGTACGT", 2) == 1
    assert ug.levenshtein("ACGTACGTACGT", "ACGACGTACGT", 2) == 1


def test_levenshtein_sees_a_shifted_umi_as_two_edits():
    # A deleted base early and a gained base at the end shift every base in
    # between: Hamming distance is large, Levenshtein distance is 2.
    a, b = "ACGTACGTACGT", "ACTACGTACGTA"
    assert sum(x != y for x, y in zip(a, b)) > 2
    assert ug.levenshtein(a, b, 2) == 2


def test_levenshtein_reports_max_d_plus_one_beyond_the_band():
    assert ug.levenshtein("AAAAAAAAAAAA", "AAAAAAAAACCC", 2) == 3
    assert ug.levenshtein("AAAAAAAAAAAA", "CCCCCCCCCCCC", 2) == 3


def test_levenshtein_length_difference_beyond_the_band_is_out():
    assert ug.levenshtein("ACGTACGTACGT", "ACGTACGTA", 2) == 3


def test_levenshtein_with_empty_strings():
    assert ug.levenshtein("", "", 2) == 0
    assert ug.levenshtein("", "AC", 2) == 2
    assert ug.levenshtein("ACG", "", 2) == 3


def test_levenshtein_matches_the_full_dp_within_the_band():
    rng = random.Random(38)
    for _ in range(2000):
        a = "".join(rng.choice("ACGT") for _ in range(rng.randint(8, 12)))
        b = list(a)
        for _ in range(rng.randint(0, 4)):
            op = rng.choice("sid")
            p = rng.randrange(len(b) + (op == "i"))
            if op == "s" and b:
                b[p % len(b)] = rng.choice("ACGT")
            elif op == "i":
                b.insert(p, rng.choice("ACGT"))
            elif b:
                del b[p % len(b)]
        b = "".join(b)
        for k in (0, 1, 2, 3):
            assert ug.levenshtein(a, b, k) == min(full_levenshtein(a, b), k + 1)
            assert ug.levenshtein(b, a, k) == ug.levenshtein(a, b, k)


# Candidate pairs


def test_close_pairs_finds_every_pair_within_max_d():
    rng = random.Random(7)
    umis = sorted({"".join(rng.choice("ACGT") for _ in range(6)) for _ in range(300)})
    want = {
        (i, j, full_levenshtein(umis[i], umis[j]))
        for i in range(len(umis))
        for j in range(i + 1, len(umis))
        if full_levenshtein(umis[i], umis[j]) <= 2
    }
    got = {(min(i, j), max(i, j), d) for i, j, d in ug.close_pairs(umis, 2)}
    assert got == want


# Hybrid merge


def test_ten_to_three_at_ld_two_merges_via_the_ratio():
    clusters = cluster({"ACGTACGTACGT": 10, "ACGTACGTACCA": 3})
    assert ug.levenshtein("ACGTACGTACGT", "ACGTACGTACCA", 2) == 2
    assert same_cluster(clusters, "ACGTACGTACGT", "ACGTACGTACCA")


def test_two_large_comparable_umis_at_ld_one_stay_apart():
    clusters = cluster({"ACGTACGTACGT": 20, "ACGTACGTACGA": 15})
    assert not same_cluster(clusters, "ACGTACGTACGT", "ACGTACGTACGA")


def test_two_plus_two_at_ld_one_merges_via_the_small_clause():
    counts = {"ACGTACGTACGT": 2, "ACGTACGTACGA": 2}
    assert same_cluster(cluster(counts), "ACGTACGTACGT", "ACGTACGTACGA")
    # Without the small clause the ratio alone (2 >= 2*2 - 1) does not merge
    assert not same_cluster(cluster(counts, small_reads=0), "ACGTACGTACGT", "ACGTACGTACGA")


def test_umis_beyond_max_ld_never_merge():
    clusters = cluster({"ACGTACGTACGT": 50, "ACGTACGTATTA": 1})
    assert not same_cluster(clusters, "ACGTACGTACGT", "ACGTACGTATTA")


def test_shifted_umi_merges_where_hamming_would_not():
    clusters = cluster({"ACGTACGTACGT": 10, "ACTACGTACGTA": 1})
    assert same_cluster(clusters, "ACGTACGTACGT", "ACTACGTACGTA")


def test_downhill_guard_stops_a_one_read_child_pulling_in_a_three_read_umi():
    seed, child, big = "AAAAAAAAAAAA", "AAAAAAAAAACC", "AAAAAAAACCCC"
    assert ug.levenshtein(seed, child, 2) == 2
    assert ug.levenshtein(child, big, 2) == 2
    assert ug.levenshtein(seed, big, 2) == 3
    clusters = cluster({seed: 10, child: 1, big: 3})
    assert same_cluster(clusters, seed, child)
    assert not same_cluster(clusters, child, big)


def test_clusters_grow_through_children_by_the_ratio():
    # The 1-read UMI is beyond max_ld of the seed but within it of the child
    seed, child, grandchild = "AAAAAAAAAAAA", "AAAAAAAAAACC", "AAAAAAAACCCC"
    clusters = cluster({seed: 20, child: 5, grandchild: 1})
    assert same_cluster(clusters, seed, child)
    assert same_cluster(clusters, child, grandchild)


def test_ratio_threshold_is_an_argument():
    counts = {"ACGTACGTACGT": 10, "ACGTACGTACGA": 4}
    assert same_cluster(cluster(counts, ratio=2), "ACGTACGTACGT", "ACGTACGTACGA")
    assert not same_cluster(cluster(counts, ratio=3), "ACGTACGTACGT", "ACGTACGTACGA")


def test_max_ld_threshold_is_an_argument():
    counts = {"ACGTACGTACGT": 10, "ACGTACGTACCA": 1}
    assert not same_cluster(cluster(counts, max_ld=1), "ACGTACGTACGT", "ACGTACGTACCA")


def test_tie_breaks_by_sequence_and_ignores_input_order():
    # Equal counts: the lexicographically smaller UMI seeds cluster 0, so it
    # claims the shared neighbour; the other one is left on its own.
    a, b, shared = "AAAAAAAAAAAA", "AAAAAAAACCCC", "AAAAAAAAAACC"
    counts = {b: 4, shared: 1, a: 4}
    first = cluster(counts)
    assert first[a] == 0
    assert same_cluster(first, a, shared)
    assert not same_cluster(first, b, shared)
    for _ in range(20):
        items = list(counts.items())
        random.shuffle(items)
        assert cluster(dict(items)) == first


def test_cluster_ids_follow_seed_order_by_count_then_sequence():
    clusters = cluster({"TTTTTTTTTTTT": 5, "CCCCCCCCCCCC": 9, "GGGGGGGGGGGG": 5})
    assert clusters == {"CCCCCCCCCCCC": 0, "GGGGGGGGGGGG": 1, "TTTTTTTTTTTT": 2}


def test_output_maps_every_raw_ub_to_a_cluster_id():
    counts = {"ACGTACGTACGT": 10, "ACGTACGTACGA": 1, "TTTTTTTTTTTT": 2}
    clusters = cluster(counts)
    assert set(clusters) == set(counts)
    assert sorted(set(clusters.values())) == [0, 1]


def test_empty_bundle_gives_an_empty_map():
    assert cluster({}) == {}


@pytest.mark.parametrize("bad", [{"ACGT": 0}, {"ACGT": -1}])
def test_counts_must_be_positive(bad):
    with pytest.raises(ValueError):
        cluster(bad)


# Reads and bundles (#39): primary mapped records with MAPQ >= min_mapq,
# bundled per cell x target with no position term.

HEADER = pysam.AlignmentHeader.from_dict({"SQ": [{"SN": "chr1", "LN": 10000}]})
GROUPING = dict(min_mapq=20, **THRESHOLDS)
A, B = "ACGTACGTACGT", "ACGTACGTACGA"  # Levenshtein 1


def read(start, cigar="100M", ub=A, flag=0, mapq=60, name="r"):
    """A record on chr1 from 0-based start."""
    r = pysam.AlignedSegment(HEADER)
    r.query_name = name
    r.flag = flag
    r.reference_id = 0
    r.reference_start = start
    r.mapping_quality = mapq
    r.cigarstring = cigar
    r.query_sequence = "A" * r.infer_query_length()
    if ub is not None:
        r.set_tag("UB", ub)
    return r


def reads(n, start, ub, prefix, **kw):
    return [read(start, ub=ub, name=f"{prefix}{i}", **kw) for i in range(n)]


def target(start, ref, alt, chrom="chr1"):
    """A BED target; start is the 1-based bam-readcount position."""
    return ug.Target.from_bed(chrom, start, ref, alt)


def assignments(records, targets, **overrides):
    """(target_id, read name) -> cluster ID of every grouped record."""
    return {(tid, r.query_name): c
            for tid, c, r in ug.group_reads(records, targets, **{**GROUPING, **overrides})}


def test_target_id_is_chr_start_ref_alt():
    assert target(1002, "T", "-TTA").target_id == "chr1:1002:T>-TTA"


@pytest.mark.parametrize("flag", [0x4, 0x100, 0x800, 0x904])
def test_unmapped_secondary_and_supplementary_records_are_dropped(flag):
    assert not ug.is_primary_mapped(read(0, flag=flag), min_mapq=20)


def test_reverse_strand_primary_records_are_kept():
    assert ug.is_primary_mapped(read(0, flag=0x10), min_mapq=20)


def test_mapq_below_the_minimum_is_dropped_and_equal_is_kept():
    assert not ug.is_primary_mapped(read(0, mapq=19), min_mapq=20)
    assert ug.is_primary_mapped(read(0, mapq=20), min_mapq=20)


def test_snv_bundle_holds_reads_whose_alignment_covers_the_site():
    snv = target(501, "A", "G")  # 0-based 500
    assert ug.in_bundle(read(401), snv)  # 401..500
    assert ug.in_bundle(read(500), snv)
    assert not ug.in_bundle(read(400), snv)  # ends on 499
    assert not ug.in_bundle(read(501), snv)


def test_reads_on_another_chromosome_are_not_in_the_bundle():
    assert not ug.in_bundle(read(450), target(501, "A", "G", chrom="chr2"))


def test_soft_clips_are_not_in_the_alignment():
    assert not ug.in_bundle(read(400, "100M50S"), target(501, "A", "G"))


def test_deletion_bundle_is_the_deleted_span_plus_minus_k():
    # -TTA deletes 0-based 1001..1003; with k = 1 the bundle spans 1000..1004
    deletion = target(1002, "T", "-TTA")
    assert ug.INDEL_K == 1
    assert ug.in_bundle(read(901), deletion)  # ends on 1000
    assert not ug.in_bundle(read(900), deletion)  # ends on 999
    assert ug.in_bundle(read(1004), deletion)
    assert not ug.in_bundle(read(1005), deletion)


def test_insertion_bundle_is_the_anchor_plus_minus_k():
    # +GTC is anchored on 0-based 1300; with k = 1 the bundle spans 1299..1301
    insertion = target(1301, "A", "+GTC")
    assert ug.in_bundle(read(1200), insertion)  # ends on 1299
    assert not ug.in_bundle(read(1199), insertion)
    assert ug.in_bundle(read(1301), insertion)
    assert not ug.in_bundle(read(1302), insertion)


def test_indel_bundle_counts_deleted_bases_in_the_reference_span():
    # 50 aligned bases end on 949; a deletion carries the span on, so a read
    # whose last aligned bases stop short of the target can still reach it
    deletion = target(1002, "T", "-TTA")
    assert not ug.in_bundle(read(900, "50M49D1M"), deletion)  # span ends on 999
    assert ug.in_bundle(read(900, "50M50D1M"), deletion)  # span ends on 1000
    assert ug.in_bundle(read(900, "50M200D50M"), deletion)  # deletion over the target


def test_a_molecule_is_one_cluster_however_far_apart_its_reads_start():
    # The reads start either side of a former 100 bp window boundary
    got = assignments([read(s, "150M", name=f"m{s}") for s in (590, 594, 598, 600, 604, 608)],
                      [target(701, "A", "G")])
    assert len(got) == 6
    assert set(got.values()) == {0}


def test_ten_to_three_ld_two_pair_is_one_cluster():
    major, minor = "ACGTACGTACGT", "ACTACGTACGTA"
    got = assignments(reads(10, 250, major, "major") + reads(3, 250, minor, "minor"),
                      [target(301, "A", "G")])
    assert len(got) == 13
    assert set(got.values()) == {0}


def test_secondary_supplementary_and_low_mapq_records_are_not_grouped():
    records = [read(250, name="primary"), read(250, flag=0x100, name="secondary"),
               read(250, flag=0x800, name="supplementary"), read(250, mapq=5, name="lowmapq")]
    assert assignments(records, [target(301, "A", "G")]) == {("chr1:301:A>G", "primary"): 0}


def test_dropped_records_do_not_count_as_reads_per_ub():
    # 4 vs 4 primaries stay apart; A's 4 secondaries, if counted, would make
    # it 8 >= 2*4 - 1 and merge B
    records = reads(4, 250, A, "a") + reads(4, 250, A, "a_sec", flag=0x100) + reads(4, 250, B, "b")
    got = assignments(records, [target(301, "A", "G")])
    assert got[("chr1:301:A>G", "a0")] != got[("chr1:301:A>G", "b0")]


def test_reads_without_a_ub_are_not_grouped():
    assert assignments([read(250, ub=None)], [target(301, "A", "G")]) == {}


def test_reads_outside_every_bundle_are_not_grouped():
    assert assignments([read(5000)], [target(301, "A", "G")]) == {}


def test_a_read_joins_the_bundle_of_every_target_it_reaches():
    near, far = target(301, "A", "G"), target(341, "C", "T")  # 0-based 300, 340
    records = reads(4, 250, A, "both")  # 250..349
    records += reads(4, 210, B, "near")  # 210..309
    got = assignments(records, [near, far])
    assert set(got) == ({(near.target_id, f"both{i}") for i in range(4)}
                        | {(near.target_id, f"near{i}") for i in range(4)}
                        | {(far.target_id, f"both{i}") for i in range(4)})


def test_reads_per_ub_are_counted_within_the_bundle():
    # In the first bundle B has 1 read and is absorbed by A's 4; in the second
    # B has 4 reads, as many as A, and the two stay apart
    first, second = target(301, "A", "G"), target(601, "C", "T")
    records = reads(4, 250, A, "a") + reads(1, 250, B, "b")
    records += reads(4, 550, A, "x") + reads(4, 550, B, "y")
    got = assignments(records, [first, second])
    assert got[(first.target_id, "a0")] == got[(first.target_id, "b0")]
    assert got[(second.target_id, "x0")] != got[(second.target_id, "y0")]


def test_cluster_ids_restart_per_target_in_seed_order():
    first, second = target(301, "A", "G"), target(601, "C", "T")
    records = reads(1, 250, "TTTTTTTTTTTT", "t") + reads(3, 250, "CCCCCCCCCCCC", "c")
    records += reads(2, 550, "GGGGGGGGGGGG", "g")
    got = assignments(records, [first, second])
    assert got == {
        (first.target_id, "c0"): 0, (first.target_id, "c1"): 0, (first.target_id, "c2"): 0,
        (first.target_id, "t0"): 1, (second.target_id, "g0"): 0, (second.target_id, "g1"): 0,
    }


def test_grouping_thresholds_are_arguments():
    t = target(301, "A", "G")
    records = reads(10, 250, A, "a") + reads(1, 250, B, "b")
    assert len(set(assignments(records, [t]).values())) == 1
    assert len(set(assignments(records, [t], max_ld=0).values())) == 2


# Grouped BAM on the synthetic fixtures: a record per bundle, one read group
# per cluster named <target_id>_<cluster ID>, carried as RG and MI.

FIXTURES = REPO / "tests" / "fixtures" / "synthetic"


@pytest.fixture(scope="module")
def grouped(tmp_path_factory):
    """Read name -> {target_id: cluster ID} per fixture cell, plus the BAMs."""
    out = {}
    for bam in sorted((FIXTURES / "cells").glob("*.bam")):
        path = tmp_path_factory.mktemp(bam.stem) / f"{bam.stem}.bam"
        ug.write_grouped_bam(bam, FIXTURES / "targets.bed", path, **GROUPING)
        out[bam.stem] = path
    return out


def records_of(path):
    with pysam.AlignmentFile(str(path)) as f:
        header = f.header.to_dict()
        return header, list(f)


def test_grouped_bam_has_a_read_group_per_cluster(grouped):
    for path in grouped.values():
        header, records = records_of(path)
        groups = {rg["ID"]: rg for rg in header["RG"]}
        for r in records:
            name = r.get_tag("RG")
            assert r.get_tag("MI") == name
            assert groups[name]["LB"] == name
        assert set(groups) == {r.get_tag("RG") for r in records}


def test_grouped_bam_is_sorted_and_indexed(grouped):
    for path in grouped.values():
        header, records = records_of(path)
        assert header["HD"]["SO"] == "coordinate"
        assert [r.reference_start for r in records] == sorted(r.reference_start for r in records)
        assert Path(f"{path}.bai").exists()


def clusters_of(path, scenario):
    _, records = records_of(path)
    return {r.get_tag("RG") for r in records if r.query_name.startswith(f"{scenario}:")}


def test_fixture_window_straddling_molecule_is_one_cluster(grouped):
    assert clusters_of(grouped["AAACCCAAGAAACCAT"], "straddle") == {"chrF:701:A>G_0"}


def test_fixture_ld2_pair_is_one_cluster(grouped):
    assert clusters_of(grouped["AAACCCAAGAAACACT"], "ld2_pair") == {"chrF:301:A>G_0"}


def test_fixture_secondary_and_supplementary_records_are_dropped(grouped):
    _, records = records_of(grouped["AAACCCAAGAAACACT"])
    sec_supp = [r for r in records if r.query_name.startswith("sec_supp:")]
    assert len(sec_supp) == 2
    assert not any(r.flag & 0x904 for r in sec_supp)
