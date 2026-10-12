# Characterisation tests: they pin current main behaviour of count_consensus.py,
# defects included, so the pure-function refactor provably changes nothing.
# Later tickets replace them with spec tests.

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


def label(ref, alt, *entries, min_fraction=0.6, min_consensus=2):
    return cc.label_cluster(entries, ref, alt, min_fraction, min_consensus)


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


# Per-cluster consensus labeller: SNV


def test_snv_ref_majority_is_wt():
    assert label("A", "C", ("A", 5), ("C", 1)) == "WT"


def test_snv_alt_majority_is_mut():
    assert label("A", "C", ("A", 1), ("C", 4)) == "MUT"


def test_snv_below_min_consensus_is_unlabelled():
    assert label("A", "C", ("A", 2), min_consensus=3) is None


def test_snv_min_consensus_and_min_fraction_are_inclusive():
    # 3 of 5 reads = 0.6
    assert label("A", "C", ("A", 3), ("G", 2), min_consensus=3) == "WT"


def test_snv_below_min_fraction_is_unlabelled():
    assert label("A", "C", ("A", 3), ("C", 2), ("G", 1)) is None


def test_snv_n_reads_are_ignored():
    assert label("A", "C", ("A", 3), ("N", 10)) == "WT"


def test_snv_other_bases_are_pooled_in_the_fraction():
    # G and T together make 4 of 7 reads, so A at 3/7 fails f_min
    assert label("A", "C", ("A", 3), ("G", 2), ("T", 2)) is None


def test_snv_oth_majority_is_never_labelled_mis():
    # The count tested is max(ref, alt), so an OTH-dominated cluster is dropped
    assert label("A", "C", ("A", 0), ("G", 9)) is None


def test_snv_ref_alt_tie_passing_min_fraction_raises():
    # Defect pinned as-is: a tie yields two argmax keys and an ambiguous truth value
    with pytest.raises(ValueError):
        label("A", "C", ("A", 3), ("C", 3), min_fraction=0.5)


def test_snv_last_entry_for_ref_wins_and_other_bases_add_up():
    assert label("A", "C", ("A", 9), ("A", 1), ("G", 1), ("G", 1)) is None


# Per-cluster consensus labeller: len(alt) == 2 indel shortcut


def test_short_indel_mut_needs_more_than_two_reads_and_a_quarter():
    # 3 / (8 + 3) = 0.27 > 0.25, with no f_min test
    assert label("A", "-A", ("A", 8), ("-A", 3)) == "MUT"


def test_short_indel_at_exactly_a_quarter_falls_back_to_the_general_rule():
    # 3 / 12 = 0.25, then A at 9 / 12 = 0.75 passes f_min
    assert label("A", "-A", ("A", 9), ("-A", 3)) == "WT"


def test_short_indel_with_two_reads_falls_back_to_the_general_rule():
    assert label("A", "-A", ("-A", 2)) == "MUT"
    assert label("A", "-A", ("-A", 2), min_consensus=3) is None


def test_short_indel_ignores_other_indels_in_the_shortcut():
    # -AT is not the alt, so it only counts as OTH
    assert label("A", "-A", ("A", 1), ("-AT", 9)) is None


# Per-cluster consensus labeller: longer indels


def test_three_char_indel_matches_by_substring_and_takes_the_max():
    # -AC is in -ACG: alt count = max(4, 2) = 4, general rule only
    assert label("A", "-AC", ("A", 1), ("-AC", 2), ("-ACG", 4)) == "MUT"


def test_three_char_indel_has_no_shortcut():
    # 4 reads at 4 / 9 would pass the len(alt) == 2 shortcut; here neither passes f_min
    assert label("A", "-AC", ("A", 5), ("-AC", 4)) is None


def test_long_indel_mut_needs_more_than_three_reads_and_no_fraction():
    assert label("A", "-ACC", ("A", 100), ("-ACC", 4)) == "MUT"


def test_long_indel_with_three_reads_falls_back_to_the_general_rule():
    assert label("A", "-ACC", ("A", 100), ("-ACC", 3)) == "WT"


def test_long_indel_substring_match_takes_the_max_count():
    # -AC and -ACCT both overlap -ACC; the larger, 4, enters the shortcut
    assert label("A", "-ACC", ("A", 100), ("-AC", 4), ("-ACCT", 2)) == "MUT"


def test_long_indel_short_entries_count_as_oth():
    # +T and C are too short to match, so they only add to the denominator
    assert label("A", "-ACC", ("A", 3), ("+T", 2), ("C", 1)) is None


# Cluster tally


def test_consensus_on_umi_tallies_labels_over_clusters():
    text = line(
        "chr1", 120, "A",
        block(("A", 5)),
        block(("C", 4)),
        block(("C", 3)),
        block(("A", 1)),  # below min_consensus: not counted
    )

    assert cc.consensus_on_umi(text, "A", "C", 0.6, 2) == {"WT": 1, "MUT": 2, "MIS": 0}


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
        line("chr1", 120, "A", block(("A", 5)), block(("C", 4)), block(("T", 4))),
        line("chr1", 300, "G", block(("G", 8), ("-G", 3))),
    ]

    table = cc.count_targets(targets(), lines, min_fraction=0.6, min_consensus=2)

    assert table.columns.tolist() == [
        "chr", "start", "end", "ref", "alt", "gene", "WT", "MUT", "MIS",
    ]
    assert table[["WT", "MUT", "MIS"]].values.tolist() == [
        # alt C: the T cluster is OTH at 4/4 but has max(ref, alt) = 0
        [1, 1, 0],
        [1, 1, 0],
        [0, 1, 0],
        [0, 0, 0],
    ]


def test_count_targets_skips_lines_whose_ref_does_not_match():
    lines = [line("chr1", 120, "G", block(("A", 5)))]

    table = cc.count_targets(targets(), lines, min_fraction=0.6, min_consensus=2)

    assert table[["WT", "MUT", "MIS"]].values.sum() == 0


def test_count_targets_does_not_mutate_its_input():
    t = targets()

    cc.count_targets(t, [], min_fraction=0.6, min_consensus=2)

    assert "WT" not in t.columns
