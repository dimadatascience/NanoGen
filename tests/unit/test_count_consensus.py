# Spec tests for consensus labelling (rule A, equation (1), on usable reads) and
# the sub-n_min / labelled / discarded partition, with tolerant indel ALT.

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
    return cc.classify_counts(r, a, o, min_fraction, min_reads)


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
        cls = cc.classify_counts(r, a, o, min_fraction, min_reads)
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


# Indel targets: REF is the reference base minus every '+' entry, ALT the sum of
# same-type entries within +-k bp that agree with the target over their shared
# part (tolerant match), OTH every other non-N read. Rule A applies unchanged.


@pytest.mark.parametrize(
    "entry, matches",
    [
        ("-TT", True),  # 1 bp shorter
        ("-TTA", True),  # exact
        ("-TTAC", True),  # 1 bp longer
        ("-GG", False),  # disagrees over the shared part
        ("-TTACA", False),  # 2 bp longer
        ("-T", False),  # 2 bp shorter
        ("-TAA", False),  # same length, disagrees
        ("+TTA", False),  # other type
    ],
)
def test_tolerant_indel_match(entry, matches):
    assert cc.indel_matches(entry, "-TTA") is matches


def test_deletion_alt_sums_every_tolerant_entry():
    entries = [("T", 1), ("-TT", 2), ("-TTA", 3), ("-TTAC", 1), ("-GG", 4)]

    # -GG is OTH; T reads keep their REF count at a deletion target
    assert cc.allele_counts(entries, "T", "-TTA") == (1, 6, 4)


def test_insertion_reads_are_removed_from_the_anchor_ref_count():
    # bam-readcount counts an insertion read on its anchor base and in its '+' entry
    entries = [("A", 9), ("+GTC", 3), ("+GT", 1), ("+CC", 2)]

    assert cc.allele_counts(entries, "A", "+GTC") == (3, 4, 2)


def test_every_plus_entry_is_removed_from_ref_at_a_deletion_target_too():
    assert cc.allele_counts([("T", 5), ("+G", 2), ("-TTA", 3)], "T", "-TTA") == (3, 3, 2)


def test_insertion_reads_on_a_non_ref_anchor_are_counted_once():
    # The 2 insertion reads carry C at the anchor: ALT, not also OTH
    assert cc.allele_counts([("C", 2), ("+GTC", 2)], "A", "+GTC") == (0, 2, 0)
    assert cc.allele_counts([("A", 1), ("C", 2), ("+GTC", 2)], "A", "+GTC") == (0, 2, 1)


def test_indel_other_bases_count_as_oth_and_n_is_not_usable():
    entries = [("=", 0), ("A", 4), ("C", 1), ("G", 1), ("N", 7), ("+T", 1)]

    assert cc.allele_counts(entries, "A", "-AC") == (3, 0, 3)


def test_indel_repeated_entries_add_up():
    assert cc.allele_counts([("A", 2), ("A", 1), ("-AC", 1), ("-AC", 2)], "A", "-AC") == (3, 3, 0)


def test_homopolymer_deletion_uses_the_same_rules():
    # -T in a T run: -TT is 1 bp longer and agrees, -TTT is 2 bp longer
    entries = [("T", 2), ("-T", 3), ("-TT", 2), ("-TTT", 1), ("-G", 1)]

    assert cc.allele_counts(entries, "T", "-T") == (2, 5, 2)
    assert classify("T", "-T", *entries) == cc.DISCARDED  # 5 / 9 < f_min
    assert classify("T", "-T", ("T", 1), ("-T", 3), ("-TT", 2)) == "MUT"


@pytest.mark.parametrize(
    "alt, entries, label",
    [
        # 3 ALT reads at 3 / 11 do not pass f_min, however short the alt
        ("-A", [("A", 8), ("-A", 3)], "WT"),
        # 4 ALT reads at 4 / 104 do not pass f_min, however long the alt
        ("-ACC", [("A", 100), ("-ACC", 4)], "WT"),
        ("-ACC", [("A", 5), ("-ACC", 4)], cc.DISCARDED),
    ],
)
def test_indel_labels_follow_rule_a_with_no_shortcut(alt, entries, label):
    assert classify("A", alt, *entries, min_reads=2) == label


def test_indel_entries_beyond_k_count_as_oth_even_if_they_overlap_the_alt():
    # -AC is in -ACCTG but 3 bp shorter: OTH, not ALT
    assert cc.allele_counts([("A", 1), ("-AC", 4), ("-ACCTG", 1)], "A", "-ACCTG") == (1, 1, 4)


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
        "WT": 1, "MUT": 1, "OTH": 1, "discarded": 1,
        "wt_reads": 5, "wt_altlike_reads": 0, "mixed_discarded": 1,
        "subnmin_umis": 2, "subnmin_reads": 2,
    }


def test_tally_reports_alt_like_reads_in_wt_clusters_and_mixed_discarded_clusters():
    clusters = [
        [("T", 4), ("-TT", 1)],  # WT, 1 ALT-like read of 5
        [("T", 3), ("-TTA", 1), ("G", 1)],  # WT at 3 / 5, 1 ALT-like read of 5
        [("T", 2), ("-TTAC", 2)],  # discarded, REF and ALT: mixed
        [("T", 2), ("C", 2)],  # discarded, no ALT: not mixed
        [("-TTA", 3)],  # MUT
    ]

    counts = cc.tally_clusters(clusters, "T", "-TTA", min_fraction=0.6, min_reads=3)

    assert counts == {
        "WT": 2, "MUT": 1, "OTH": 0, "discarded": 2,
        "wt_reads": 10, "wt_altlike_reads": 2, "mixed_discarded": 1,
        "subnmin_umis": 0, "subnmin_reads": 0,
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
        "WT": 1, "MUT": 2, "OTH": 0, "discarded": 0,
        "wt_reads": 5, "wt_altlike_reads": 0, "mixed_discarded": 0,
        "subnmin_umis": 1, "subnmin_reads": 1,
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
        "WT", "MUT", "OTH", "wt_reads", "wt_altlike_reads", "mixed_discarded",
        "subnmin_umis", "subnmin_reads",
    ]
    assert table[cc.COUNT_COLUMNS].values.tolist() == [
        # alt C: the T cluster is OTH; alt T: the C cluster is OTH
        [1, 1, 1, 5, 0, 0, 1, 1],
        [1, 1, 1, 5, 0, 0, 1, 1],
        # 8 REF : 3 ALT is WT by rule A
        [1, 0, 0, 11, 3, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
    ]


def test_count_targets_skips_lines_whose_ref_does_not_match():
    lines = [line("chr1", 120, "G", block(("A", 5)))]

    table = cc.count_targets(targets(), lines, min_fraction=0.6, min_reads=2)

    assert table[cc.COUNT_COLUMNS].values.sum() == 0


def test_count_targets_does_not_mutate_its_input():
    t = targets()

    cc.count_targets(t, [], min_fraction=0.6, min_reads=2)

    assert "WT" not in t.columns


def test_count_targets_scores_snv_alts_on_the_first_line_and_indel_alts_on_the_last():
    # A site with both an SNV and an indel target has one line per bam-readcount
    # pass, the SNV pass (base-quality filtered) before the indel pass (unfiltered)
    t = pd.DataFrame(
        [["chr1", 300, 300, "G", "C", "S"], ["chr1", 300, 300, "G", "-G", "D"]],
        columns=["chr", "start", "end", "ref", "alt", "gene"],
    )
    lines = [
        # Low-quality C bases and deletion anchors are dropped
        line("chr1", 300, "G", block(("G", 4)), block(), block()),
        line("chr1", 300, "G", block(("G", 4)), block(("C", 4)), block(("-G", 4))),
    ]

    table = cc.count_targets(t, lines, min_fraction=0.6, min_reads=2)

    # At the SNV, the two clusters emptied by the filter are sub-n_min; at the
    # indel, the C cluster is OTH under rule A
    assert table[["WT", "MUT", "OTH", "subnmin_umis", "subnmin_reads"]].values.tolist() == [
        [1, 0, 0, 2, 0],
        [1, 1, 1, 0, 0],
    ]
