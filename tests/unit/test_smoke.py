# Smoke test: the conftest shim makes the bin/long_reads/ scripts importable.

import count_consensus


def test_consensus_on_umi_labels_a_ref_cluster_wt():
    # One UMI cluster: 5 reads carry the reference base A, none the alt C.
    line = "{A:5:60.00 C:0:0.00}"

    counts = count_consensus.consensus_on_umi(
        line, ref="A", alt="C", min_fraction=0.6, min_reads=2
    )

    assert (counts["WT"], counts["MUT"], counts["OTH"]) == (1, 0, 0)
