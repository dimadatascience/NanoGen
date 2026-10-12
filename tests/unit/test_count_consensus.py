# Spec tests for consensus labelling (rule A, equation (1), on usable reads) and
# the sub-n_min / labelled / discarded partition. The indel tests still pin
# the legacy allele matching and shortcuts until the tolerant indel ALT ticket
# replaces them.

import itertools

import pandas as pd
import pytest

import count_consensus as cc

STATS = "60.00:40.00:60.00:1:0:0.50:0.00:40.00:1:0.70:79.00:0.70"


def block(*entries):
    """One bam-readcount -p cluster block: entries are (base, read count)."""
    fields = ["=:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00"]
    fields += [f"{base}:{n}:{STATS}" for base, n in entries]
    return "{\t" + "\t".join(fields) + "\t}"


def line(chrom, pos, ref, *blocks):
    """One bam-readcount -p line with a named cluster per block."""
    parts = [chrom, str(pos), ref, "0"]
    for i, b in enumerate(blocks):
        parts += [f"'UMI{i}_{chrom}_0'", b]
    return "\t".join(parts) + "\n"


def classify(ref, alt, *entries, min_fraction=0.6, min_reads=3):
    """Class of one cluster: a consensus label, sub-n_min or discarded."""
    r, a, o = cc.allele_counts(entries, ref, alt)
    return cc.classify_counts(r, a, o, alt, min_fraction, min_reads)


# bam-readcount line parser


def test_parse_readcount_line_splits_site_and_clusters():
    text = line(
        "chr1", 120, "A",
        block(("A", 5), ("C", 1), ("N", 0)),
        block(("A", 0), ("-T", 3)),
    )

    chrom, pos, ref, clusters = cc.parse_readcount_line(text)

    assert (chrom, pos, ref) == ("chr1", 120, "A")
    assert clusters == [
        [("=", 0), ("A", 5), ("C", 1), ("N", 0)],
        [("=", 0), ("A", 0), ("-T", 3)],
    ]


def test_parse_clusters_without_braces_reads_the_whole_line_as_one_cluster():
    # Fields without a colon (chrom, pos, ref, depth) are skipped
    assert cc.parse_clusters("chr1\t120\tA\t6\tA:5:60.00\tC:1:60.00\n") == [
        [("A", 5), ("C", 1)]
    ]


def test_parse_readcount_line_on_real_output():
    text = (
        "chrT\t40\tA\t14\t'ACGTACGTACGT_chrT_0'\t{\t"
        "=:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00\t"
        "A:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00\t"
        "N:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00\t"
        "-A:6:60.00:40.00:60.00:6:0:0.48:0.00:40.00:6:0.73:79.00:0.73\t}\n"
    )

    assert cc.parse_readcount_line(text) == (
        "chrT", 40, "A", [[("=", 0), ("A", 0), ("N", 0), ("-A", 6)]]
    )


# Usable reads at an SNV. Base quality is filtered upstream, by bam-readcount -b.


def test_snv_usable_reads_split_into_ref_alt_and_other_bases():
    entries = [("=", 0), ("A", 5), ("C", 2), ("G", 1), ("T", 1)]

    assert cc.allele_counts(entries, "A", "C") == (5, 2, 2)


def test_snv_n_reads_are_not_usable():
    assert cc.allele_counts([("A", 3), ("N", 10)], "A", "C") == (3, 0, 0)


def test_snv_reads_with_a_deletion_over_the_site_are_not_usable():
    # bam-readcount reports a deletion on its deleted bases, outside the base counts
    assert cc.allele_counts([("A", 3), ("-A", 4), ("-AG", 1)], "A", "C") == (3, 0, 0)


def test_snv_insertion_entries_are_not_counted_twice():
    # A read with an insertion after the site is already in its base count
    assert cc.allele_counts([("A", 3), ("C", 1), ("+GT", 2)], "A", "C") == (3, 1, 0)


def test_snv_repeated_entries_add_up():
    assert cc.allele_counts([("A", 2), ("A", 1), ("G", 1), ("G", 1)], "A", "C") == (3, 0, 2)


# Rule A, equation (1): label X iff count_X >= n_min and count_X / (r + a + o) >= f_min


@pytest.mark.parametrize(
    "r, a, o, label",
    [
        (3, 0, 2, "WT"),  # 3 reads at 3 / 5 = 0.6: both thresholds are inclusive
        (0, 3, 2, "MUT"),
        (2, 0, 3, "OTH"),
        (1, 0, 9, "OTH"),  # OTH is tested on o itself
        (3, 2, 1, None),  # 3 / 6 < f_min
        (1, 2, 0, None),  # 2 / 3 >= f_min but 2 < n_min
        (0, 0, 0, None),
    ],
)
def test_consensus_label_follows_rule_a(r, a, o, label):
    assert cc.consensus_label(r, a, o, min_fraction=0.6, min_reads=3) == label


def test_consensus_label_uses_f_min():
    assert cc.consensus_label(3, 1, 0, min_fraction=0.75, min_reads=3) == "WT"
    assert cc.consensus_label(3, 1, 1, min_fraction=0.75, min_reads=3) is None


# Partition of the clusters covering a target


def test_cluster_below_n_min_usable_reads_is_sub_n_min_even_if_unanimous():
    assert classify("A", "C", ("A", 2)) == cc.SUBNMIN


def test_reads_that_are_not_usable_do_not_lift_a_cluster_out_of_sub_n_min():
    assert classify("A", "C", ("A", 2), ("N", 5), ("-A", 3)) == cc.SUBNMIN


def test_covering_cluster_with_no_usable_read_is_sub_n_min():
    assert classify("A", "C", ("=", 0), ("A", 0), ("C", 0), ("N", 0)) == cc.SUBNMIN


def test_two_alt_one_ref_at_n_min_3_is_discarded():
    assert classify("A", "C", ("A", 1), ("C", 2), min_reads=3) == cc.DISCARDED


def test_cluster_with_enough_usable_reads_but_no_majority_is_discarded():
    assert classify("A", "C", ("A", 3), ("C", 2), ("G", 1)) == cc.DISCARDED


@pytest.mark.parametrize("min_reads", [1, 2, 3, 4])
@pytest.mark.parametrize("min_fraction", [0.6, 0.75, 1.0])
def test_every_cluster_is_exactly_one_of_sub_n_min_labelled_or_discarded(
    min_reads, min_fraction
):
    for r, a, o in itertools.product(range(7), repeat=3):
        cls = cc.classify_counts(r, a, o, "C", min_fraction, min_reads)
        usable = r + a + o
        qualifies = {
            label
            for label, n in (("WT", r), ("MUT", a), ("OTH", o))
            if usable and n >= min_reads and n / usable >= min_fraction
        }
        if usable < min_reads:
            assert cls == cc.SUBNMIN
        elif qualifies:
            assert {cls} == qualifies
        else:
            assert cls == cc.DISCARDED


# Indel targets: the legacy allele matching and shortcuts are kept until the
# tolerant indel ALT ticket replaces them; rule A and the partition apply.


def test_short_indel_mut_needs_more_than_two_reads_and_a_quarter():
    # 3 / (8 + 3) = 0.27 > 0.25, with no f_min test
    assert classify("A", "-A", ("A", 8), ("-A", 3), min_reads=2) == "MUT"


def test_short_indel_at_exactly_a_quarter_falls_back_to_rule_a():
    # 3 / 12 = 0.25, then A at 9 / 12 = 0.75 passes f_min
    assert classify("A", "-A", ("A", 9), ("-A", 3), min_reads=2) == "WT"


def test_short_indel_with_two_reads_falls_back_to_rule_a():
    assert classify("A", "-A", ("-A", 2), min_reads=2) == "MUT"
    assert classify("A", "-A", ("-A", 2), min_reads=3) == cc.SUBNMIN


def test_short_indel_other_indels_count_as_oth():
    # -AT is not the alt, so it counts as OTH, which rule A tests on o itself
    assert classify("A", "-A", ("A", 1), ("-AT", 9), min_reads=2) == "OTH"


def test_three_char_indel_matches_by_substring_and_takes_the_max():
    # -AC is in -ACG: alt count = max(2, 4) = 4, rule A only
    assert classify("A", "-AC", ("A", 1), ("-AC", 2), ("-ACG", 4), min_reads=2) == "MUT"


def test_three_char_indel_has_no_shortcut():
    # 4 reads at 4 / 9 would pass the len(alt) == 2 shortcut; here neither passes f_min
    assert classify("A", "-AC", ("A", 5), ("-AC", 4), min_reads=2) == cc.DISCARDED


def test_long_indel_mut_needs_more_than_three_reads_and_no_fraction():
    assert classify("A", "-ACC", ("A", 100), ("-ACC", 4), min_reads=2) == "MUT"


def test_long_indel_with_three_reads_falls_back_to_rule_a():
    assert classify("A", "-ACC", ("A", 100), ("-ACC", 3), min_reads=2) == "WT"


def test_long_indel_substring_match_takes_the_max_count():
    # -AC and -ACCT both overlap -ACC; the larger, 4, enters the shortcut
    assert classify("A", "-ACC", ("A", 100), ("-AC", 4), ("-ACCT", 2), min_reads=2) == "MUT"


def test_long_indel_short_entries_count_as_oth():
    # +T and C are too short to match: r = 3, o = 3
    assert classify("A", "-ACC", ("A", 3), ("+T", 2), ("C", 1), min_reads=2) == cc.DISCARDED


def test_indel_shortcut_never_labels_a_sub_n_min_cluster():
    assert classify("A", "-ACC", ("-ACC", 4), min_reads=5) == cc.SUBNMIN


# Cluster tally


def test_tally_puts_every_covering_cluster_in_exactly_one_class():
    clusters = [
        [("A", 5)],  # WT
        [("C", 4)],  # MUT
        [("G", 3), ("T", 1)],  # OTH
        [("A", 2)],  # sub-n_min, 2 usable reads
        [("A", 0), ("N", 4)],  # sub-n_min, no usable read
        [("A", 1), ("C", 2)],  # discarded
    ]

    counts = cc.tally_clusters(clusters, "A", "C", min_fraction=0.6, min_reads=3)

    assert counts == {
        "WT": 1, "MUT": 1, "OTH": 1, "discarded": 1, "subnmin_umis": 2, "subnmin_reads": 2,
    }


def test_consensus_on_umi_tallies_a_bam_readcount_line():
    text = line(
        "chr1", 120, "A",
        block(("A", 5)),
        block(("C", 4)),
        block(("C", 3)),
        block(("A", 1)),
    )

    assert cc.consensus_on_umi(text, "A", "C", 0.6, 2) == {
        "WT": 1, "MUT": 2, "OTH": 0, "discarded": 0, "subnmin_umis": 1, "subnmin_reads": 1,
    }


# Per-cell target counting


def targets():
    return pd.DataFrame(
        [
            ["chr1", 120, 120, "A", "C", "G1"],
            ["chr1", 120, 120, "A", "T", "G1"],  # multiallelic site
            ["chr1", 300, 300, "G", "-G", "G2"],
            ["chr2", 50, 50, "T", "A", "G3"],  # never covered
        ],
        columns=["chr", "start", "end", "ref", "alt", "gene"],
    )


def test_count_targets_scores_every_alt_at_a_covered_site():
    lines = [
        line(
            "chr1", 120, "A",
            block(("A", 5)), block(("C", 4)), block(("T", 4)), block(("A", 1)),
        ),
        line("chr1", 300, "G", block(("G", 8), ("-G", 3))),
    ]

    table = cc.count_targets(targets(), lines, min_fraction=0.6, min_reads=2)

    assert table.columns.tolist() == [
        "chr", "start", "end", "ref", "alt", "gene",
        "WT", "MUT", "OTH", "subnmin_umis", "subnmin_reads",
    ]
    assert table[["WT", "MUT", "OTH", "subnmin_umis", "subnmin_reads"]].values.tolist() == [
        # alt C: the T cluster is OTH; alt T: the C cluster is OTH
        [1, 1, 1, 1, 1],
        [1, 1, 1, 1, 1],
        [0, 1, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ]


def test_count_targets_skips_lines_whose_ref_does_not_match():
    lines = [line("chr1", 120, "G", block(("A", 5)))]

    table = cc.count_targets(targets(), lines, min_fraction=0.6, min_reads=2)

    assert table[["WT", "MUT", "OTH", "subnmin_umis", "subnmin_reads"]].values.sum() == 0


def test_count_targets_does_not_mutate_its_input():
    t = targets()

    cc.count_targets(t, [], min_fraction=0.6, min_reads=2)

    assert "WT" not in t.columns
