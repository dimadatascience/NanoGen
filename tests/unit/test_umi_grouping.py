# Spec tests for the UMI grouping rule (#38): banded Levenshtein and the
# hybrid merge within one bundle.

import random

import pytest

import umi_grouping as ug

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
