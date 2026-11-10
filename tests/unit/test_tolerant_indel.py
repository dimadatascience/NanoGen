# PROTOTYPE tests (position-tolerant indel matching, never merged): the read
# classifier in umi_grouping.py and the move to ALT in count_consensus.py.

import pysam
import pytest

import count_consensus as cc
import umi_grouping as ug

HEADER = pysam.AlignmentHeader.from_dict({"SQ": [{"SN": "chr1", "LN": 10000}]})
# +TCTG anchored on 1-based 1001 (0-based 1000); -TTA from 1-based 2002 (0-based 2001)
INS = ug.IndelTarget("chr1:1001:C>+TCTG", "chr1", 1000, "C", "+", "TCTG")
DEL = ug.IndelTarget("chr1:2002:T>-TTA", "chr1", 2001, "T", "-", "TTA")


def read(start, cigar, seq=None, ub="ACGTACGTACGT", name="r"):
    """A record on chr1 from 0-based start; seq defaults to all C."""
    r = pysam.AlignedSegment(HEADER)
    r.query_name = name
    r.reference_id = 0
    r.reference_start = start
    r.mapping_quality = 60
    r.cigarstring = cigar
    r.query_sequence = seq or "C" * r.infer_query_length()
    r.set_tag("UB", ub)
    return r


def insertion_read(offset, ins="TCTG"):
    """50 M to 0-based 999, then the insertion after reference base 1000 + offset."""
    left = 1000 + offset - 950 + 1
    seq = "C" * left + ins + "C" * 50
    return read(950, f"{left}M{len(ins)}I50M", seq)


def deletion_read(offset, length=3):
    """Deletion of `length` bases from 0-based 2001 + offset."""
    left = 2001 + offset - 1950
    return read(1950, f"{left}M{length}D50M")


def classes(t, r, window=ug.INDEL_WINDOW):
    indels, base = ug.read_indels(r, t.p0)
    return ug.v3_anchor_class(t, indels, base), ug.is_tolerant(t, indels, window)


# Read classifier


def test_exact_insertion_at_the_anchor_is_v3_alt_and_tolerant():
    assert classes(INS, insertion_read(0)) == (ug.ANCHOR_ALT, True)


def test_reference_read_is_ref_and_not_tolerant():
    assert classes(INS, read(950, "100M")) == (ug.ANCHOR_REF, False)


@pytest.mark.parametrize("offset", [-12, -2, 3, 12])
def test_shifted_insertion_within_the_window_is_tolerant_and_ref_at_the_anchor(offset):
    assert classes(INS, insertion_read(offset)) == (ug.ANCHOR_REF, True)


@pytest.mark.parametrize("offset", [-13, 13])
def test_shifted_insertion_beyond_the_window_is_not_tolerant(offset):
    assert classes(INS, insertion_read(offset)) == (ug.ANCHOR_REF, False)


def test_any_sequence_of_a_compatible_length_is_tolerant():
    # v3 needs the prefix to agree: at the anchor this read is OTH, and tolerant
    assert classes(INS, insertion_read(0, "GGGG")) == (ug.ANCHOR_OTH, True)
    assert classes(INS, insertion_read(5, "GGG")) == (ug.ANCHOR_REF, True)


def test_insertion_length_beyond_k_is_not_tolerant():
    assert classes(INS, insertion_read(0, "TCTGTC")) == (ug.ANCHOR_OTH, False)
    assert classes(INS, insertion_read(4, "TC")) == (ug.ANCHOR_REF, False)


def test_a_deletion_does_not_match_an_insertion_target():
    _, tol = classes(INS, read(950, "51M4D50M"))
    assert not tol


def test_window_minus_one_switches_the_rule_off():
    assert classes(INS, insertion_read(0, "GGGG"), window=-1) == (ug.ANCHOR_OTH, False)
    assert classes(INS, insertion_read(3), window=-1) == (ug.ANCHOR_REF, False)
    # v3's own ALT is unaffected
    assert classes(INS, insertion_read(0), window=-1)[0] == ug.ANCHOR_ALT


def test_window_zero_accepts_any_sequence_at_the_anchor_only():
    assert classes(INS, insertion_read(0, "GGGG"), window=0) == (ug.ANCHOR_OTH, True)
    assert classes(INS, insertion_read(1), window=0) == (ug.ANCHOR_REF, False)


def test_deletion_at_its_first_deleted_base_is_v3_alt():
    assert classes(DEL, deletion_read(0)) == (ug.ANCHOR_ALT, True)
    assert classes(DEL, deletion_read(0, length=4)) == (ug.ANCHOR_ALT, True)


def test_deletion_length_beyond_k_at_the_anchor_is_oth():
    assert classes(DEL, deletion_read(0, length=5)) == (ug.ANCHOR_OTH, False)


def test_shifted_deletion_over_the_anchor_is_not_usable_there_and_tolerant():
    # From 0-based 2000, the deletion covers the anchor: bam-readcount has no base there
    assert classes(DEL, deletion_read(-1)) == (ug.ANCHOR_NONE, True)


def test_shifted_deletion_after_the_anchor_is_ref_there_and_tolerant():
    r = read(1950, "60M3D50M", "C" * 51 + "T" + "C" * 58)  # anchor base T, deletion from 2010
    assert classes(DEL, r) == (ug.ANCHOR_REF, True)


def test_non_reference_anchor_base_is_oth():
    r = read(1000, "100M", "A" + "C" * 99)  # A on the anchor, whose reference base is C
    assert ug.v3_anchor_class(INS, *ug.read_indels(r, INS.p0)) == ug.ANCHOR_OTH


def test_insertion_at_the_anchor_after_an_n_base_is_still_classed_by_the_entry():
    seq = "C" * 50 + "N" + "TCTG" + "C" * 50
    assert classes(INS, read(950, "51M4I50M", seq)) == (ug.ANCHOR_ALT, True)


# Per-cluster counts


def test_tolerant_counts_split_tolerant_only_reads_by_anchor_class():
    tid = INS.target_id
    grouped = [(tid, 0, insertion_read(0)), (tid, 0, insertion_read(3)),
               (tid, 0, insertion_read(0, "GGGG")), (tid, 0, read(950, "100M")),
               (tid, 1, insertion_read(-5)),
               ("chr1:500:A>G", 0, insertion_read(3))]  # an SNV bundle: ignored
    got = ug.tolerant_counts(grouped, [INS])
    assert set(got) == {f"{tid}_0", f"{tid}_1"}
    assert got[f"{tid}_0"] == dict(shift_ref=1, shift_oth=1, shift_none=0,
                                   v3_ref=2, v3_alt=1, v3_oth=1)
    assert got[f"{tid}_1"] == dict(shift_ref=1, shift_oth=0, shift_none=0,
                                   v3_ref=1, v3_alt=0, v3_oth=0)


def test_read_indel_targets_keeps_indels_only(tmp_path):
    bed = tmp_path / "t.bed"
    bed.write_text("chr1\t500\t500\tA\tG\tg1\nchr1\t1001\t1001\tC\t+TCTG\tg2\n"
                   "chr1\t2002\t2002\tT\t-TTA\tg3\n")
    assert ug.read_indel_targets(bed) == [INS, DEL]


def test_tolerant_tsv_round_trips(tmp_path):
    counts = {"chr1:1001:C>+TCTG_0": dict(shift_ref=2, shift_oth=1, shift_none=0,
                                          v3_ref=5, v3_alt=1, v3_oth=1)}
    path = tmp_path / "tol.tsv"
    ug.write_tolerant_counts(counts, path)
    assert cc.read_tolerant(path) == counts


# The move to ALT in count_consensus


def tol(shift_ref=0, shift_oth=0, shift_none=0, v3=(0, 0, 0)):
    return dict(shift_ref=shift_ref, shift_oth=shift_oth, shift_none=shift_none,
                v3_ref=v3[0], v3_alt=v3[1], v3_oth=v3[2])


def test_apply_tolerant_moves_ref_and_oth_reads_and_adds_unusable_ones():
    assert cc.apply_tolerant(5, 1, 2, tol(2, 1, 1, v3=(5, 1, 2))) == (3, 5, 1, 4, False)


def test_apply_tolerant_flags_a_classifier_mismatch_and_never_goes_negative():
    r, a, o, moved, mismatch = cc.apply_tolerant(1, 0, 0, tol(3, 0, 0, v3=(3, 0, 0)))
    assert (r, a, o, moved, mismatch) == (0, 3, 0, 3, True)


def test_tally_turns_a_ref_cluster_with_shifted_insertions_into_mut():
    tid = INS.target_id
    entries = [("C", 4), ("+TCTG", 1)]  # 1 exact, 3 REF reads after the insertion's anchor
    clusters, names = [entries], [f"{tid}_0"]
    v3 = cc.tally_clusters(clusters, "C", "+TCTG", names=names)
    assert v3["MUT"] == 0
    got = cc.tally_clusters(clusters, "C", "+TCTG", names=names,
                            tolerant={f"{tid}_0": tol(2, 0, 0, v3=(3, 1, 0))})
    assert got["MUT"] == 1 and got["tol_reads"] == 2 and got["brc_mismatch_umis"] == 0


def test_tally_ignores_tolerant_counts_at_snv_targets():
    got = cc.tally_clusters([[("A", 3)]], "A", "G", names=["x_0"],
                            tolerant={"x_0": tol(3, 0, 0, v3=(3, 0, 0))})
    assert got["WT"] == 1 and got["tol_reads"] == 0


def test_count_targets_reports_tolerant_reads_of_clusters_not_at_the_anchor():
    import pandas as pd
    tid = INS.target_id
    target = pd.DataFrame([["chr1", 1001, 1001, "C", "+TCTG", "g"]],
                          columns=["chr", "start", "end", "ref", "alt", "gene"])
    stats = "60.00:40.00:60.00:1:0:0.50:0.00:40.00:1:0.70:79.00:0.70"
    blk = "{\t" + f"C:4:{stats}\t+TCTG:1:{stats}" + "\t}"
    lines = ["\t".join(["chr1", "1001", "C", "0", f"{tid}_0", blk]) + "\n"]
    tolerant = {f"{tid}_0": tol(2, 0, 0, v3=(3, 1, 0)), f"{tid}_7": tol(0, 0, 3)}
    got = cc.count_targets(target, lines, tolerant=tolerant).iloc[0]
    assert (got["MUT"], got["tol_reads"], got["tol_unseen_reads"], got["brc_mismatch_umis"]) == (1, 2, 3, 0)
    plain = cc.count_targets(target, lines).iloc[0]
    assert (plain["MUT"], plain["tol_reads"]) == (0, 0)
