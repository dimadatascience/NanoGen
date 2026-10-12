# Spec tests for genotyping.py: the per-target error pooled over cells, the
# binomial caller and its decision rule (genotyping.tex, sections 3-5), and the
# <sample>.csv and <sample>_targets.csv contracts.

import math
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import genotyping as gt

SPEC = Path(__file__).resolve().parents[2] / "docs" / "genotyping" / "genotyping.tex"


def collapsed(*rows):
    """Collapsed consensus counts as COLLAPSE writes them: rows are
    (chr, start, ref, alt, gene, WT, MUT, OTH, barcode)."""
    table = pd.DataFrame(
        rows,
        columns=["chr", "start", "ref", "alt", "gene", "WT", "MUT", "OTH", "barcode"],
    )
    table.insert(2, "end", table["start"])
    return table


def snv(barcode, wt, mut, oth):
    return ("chrF", 301, "A", "G", "SNV1", wt, mut, oth, barcode)


def deletion(barcode, wt, mut, oth):
    return ("chrF", 1002, "T", "-TTA", "DEL", wt, mut, oth, barcode)


def is_na(x):
    return isinstance(x, float) and math.isnan(x)


# Spelled out rather than read from gt, so the contracts are pinned
OUTPUT_COLUMNS = [
    "cell", "target_id", "gene", "type", "WT", "MUT", "OTH",
    "P_s", "p", "genotype", "het_detect_lb",
]
TARGET_COLUMNS = [
    "target_id", "gene", "type", "N", "w", "d", "q_hat", "eps", "P_s", "floor_applied",
]


# Per-target error


def errors(*rows, error_floor=1e-4):
    return gt.target_error(gt.annotate(collapsed(*rows)), error_floor)


def test_n_and_w_pool_labelled_umis_over_cells():
    table = errors(snv("c1", 3, 1, 1), snv("c2", 5, 0, 2), snv("c3", 0, 0, 0))

    assert table[["N", "w"]].values.tolist() == [[12, 3]]


def test_snv_error_halves_the_oth_fraction():
    table = errors(snv("c1", 3, 1, 1), snv("c2", 5, 0, 2))

    assert table["d"].tolist() == [2]
    assert table["q_hat"].tolist() == pytest.approx([3 / 24])


def test_indel_error_is_the_oth_fraction():
    table = errors(deletion("c1", 3, 1, 1), deletion("c2", 5, 0, 2))

    assert table["d"].tolist() == [1]
    assert table["q_hat"].tolist() == pytest.approx([3 / 12])


def test_error_above_the_floor_is_used():
    table = errors(snv("c1", 2, 1, 1), error_floor=1e-4)

    assert table["P_s"].tolist() == pytest.approx([1 / 8])
    assert table["floor_applied"].tolist() == [False]


def test_error_below_the_floor_is_raised_to_it():
    table = errors(snv("c1", 4, 0, 0), error_floor=0.01)

    assert table["q_hat"].tolist() == [0.0]
    assert table["eps"].tolist() == [0.01]
    assert table["P_s"].tolist() == [0.01]
    assert table["floor_applied"].tolist() == [True]


def test_no_minimum_n_a_single_umi_estimates_the_error():
    table = errors(deletion("c1", 0, 0, 1), error_floor=1e-4)

    assert table[["N", "w", "q_hat", "P_s"]].values.tolist() == [[1, 1, 1.0, 1.0]]


def test_target_without_labelled_umis_has_no_estimate_and_sits_at_the_floor():
    table = errors(snv("c1", 0, 0, 0), snv("c2", 0, 0, 0), error_floor=1e-4)

    assert table["N"].tolist() == [0]
    assert is_na(table["q_hat"].tolist()[0])
    assert table["P_s"].tolist() == [1e-4]
    assert table["floor_applied"].tolist() == [True]


def test_targets_table_has_one_row_per_target_in_input_order():
    table = errors(
        deletion("c1", 1, 0, 0), snv("c1", 1, 0, 0),
        deletion("c2", 1, 0, 0), snv("c2", 1, 0, 0),
    )

    assert table.columns.tolist() == TARGET_COLUMNS
    assert table[["target_id", "gene", "type"]].values.tolist() == [
        ["chrF:1002:T>-TTA", "DEL", "INDEL"],
        ["chrF:301:A>G", "SNV1", "SNV"],
    ]


def test_errors_are_estimated_per_target_not_across_targets():
    table = errors(snv("c1", 3, 0, 1), deletion("c1", 4, 0, 0))

    assert table["q_hat"].tolist() == pytest.approx([1 / 8, 0.0])


# Caller


def call(mut, wt, p_s, min_mut_umis=1):
    p, called = gt.call(np.array([mut]), np.array([wt]), np.array([p_s]), min_mut_umis)
    return p[0], called[0]


def test_thresholds_are_named_constants():
    assert gt.MUT_P_MAX == 0.1
    assert gt.WT_P_MIN == 0.9


def test_p_is_the_binomial_upper_tail_of_mut_given_m_plus_k():
    # P(X >= 2), X ~ Bin(3, 0.5) = 3/8 + 1/8
    assert call(2, 1, 0.5)[0] == pytest.approx(0.5)


def test_tail_starts_at_j_0_so_only_mut_umis_is_not_p_0():
    p, called = call(1, 0, 1e-4)

    assert p == pytest.approx(1e-4)
    assert called == "MUT"


def test_no_wt_or_mut_umis_gives_nan_p_and_na_call():
    p, called = call(0, 0, 1e-4)
    assert math.isnan(p) and is_na(called)


def test_only_oth_umis_gives_nan_p_and_na_call():
    p, called = gt.call(np.array([0]), np.array([0]), np.array([0.2]), 1)
    assert math.isnan(p[0]) and is_na(called[0])


def test_no_mut_umis_and_two_wt_is_wt():
    p, called = call(0, 2, 1e-4)
    assert p == 1.0 and called == "WT"


def test_no_mut_umis_and_one_wt_is_undetermined():
    p, called = call(0, 1, 1e-4)
    assert p == 1.0 and is_na(called)


def test_p_below_0_1_is_mut():
    p, called = call(1, 2, 1e-4)
    assert p < 0.1 and called == "MUT"


def test_p_at_least_0_9_with_mut_umis_is_wt():
    p, called = call(1, 30000, 1e-4)
    assert p >= 0.9 and called == "WT"


def test_p_between_0_1_and_0_9_is_undetermined():
    p, called = call(1, 5000, 1e-4)
    assert 0.1 <= p < 0.9 and is_na(called)


def test_min_mut_umis_2_makes_a_single_mut_umi_undetermined():
    p, called = call(1, 2, 1e-4, min_mut_umis=2)
    assert p < 0.1 and is_na(called)


def test_min_mut_umis_2_calls_two_mut_umis_mut():
    assert call(2, 2, 1e-4, min_mut_umis=2)[1] == "MUT"


def test_min_mut_umis_2_keeps_wt_calls():
    assert call(0, 2, 1e-4, min_mut_umis=2)[1] == "WT"


def test_het_detect_lb_is_filled_for_wt_calls_only():
    lb = gt.het_detect_lb(
        np.array([2, 4, 3, 1]), np.array(["WT", "WT", "MUT", np.nan], dtype=object)
    )

    assert lb[:2].tolist() == [0.75, 0.9375]
    assert np.isnan(lb[2:]).all()


# The document's Table 1: p for a single MUT UMI at P_s = 1e-4, and its call


def _table_1():
    tex = SPEC.read_text()
    body = tex[tex.index("\\midrule"):tex.index("\\bottomrule")]
    rows = []
    for line in body.splitlines():
        cells = [c.strip() for c in line.rstrip("\\ ").split("&")]
        if len(cells) == 4:
            rows.append((int(cells[0]), int(cells[1]), cells[2], cells[3]))
    return rows


TABLE_1 = _table_1()


def test_table_1_is_read_from_the_document():
    assert [(m, k) for m, k, _, _ in TABLE_1] == [
        (1, 1), (1, 10), (1, 100), (1, 1000), (1, 5000),
    ]


@pytest.mark.parametrize("m, k, p_shown, call_shown", TABLE_1)
def test_table_1(m, k, p_shown, call_shown):
    p, called = call(m, k, 1e-4)

    decimals = len(re.sub(r"^\d*\.", "", p_shown))
    assert p == pytest.approx(float(p_shown), abs=0.5 * 10**-decimals)
    assert (called if isinstance(called, str) else "undetermined") == call_shown


# <sample>.csv


def test_target_id_is_chr_position_ref_and_alt_as_in_the_bed():
    assert gt.target_id("chr4", 54733160, "GAT", "-GAT") == "chr4:54733160:GAT>-GAT"
    assert gt.target_id("chrF", 301, "A", "G") == "chrF:301:A>G"


@pytest.mark.parametrize(
    "alt, expected",
    [("G", "SNV"), ("-T", "INDEL"), ("-TTA", "INDEL"), ("+GTC", "INDEL")],
)
def test_type_is_indel_for_bam_readcount_indel_alts_else_snv(alt, expected):
    assert gt.variant_type(alt) == expected


def test_output_table_has_the_v2_columns_one_row_per_cell_and_target():
    table, _ = gt.call_targets(
        collapsed(
            snv("c1", 2, 3, 0),
            deletion("c1", 0, 0, 0),
            snv("c2", 0, 0, 0),
        ),
        error_floor=1e-4,
    )

    assert table.columns.tolist() == OUTPUT_COLUMNS
    assert table[["cell", "target_id", "gene", "type"]].values.tolist() == [
        ["c1", "chrF:301:A>G", "SNV1", "SNV"],
        ["c1", "chrF:1002:T>-TTA", "DEL", "INDEL"],
        ["c2", "chrF:301:A>G", "SNV1", "SNV"],
    ]
    assert table[["WT", "MUT", "OTH"]].values.tolist() == [[2, 3, 0], [0, 0, 0], [0, 0, 0]]
    assert table["genotype"].tolist()[0] == "MUT"


def test_each_row_uses_its_targets_pooled_p_s():
    # SNV: w / 2N = 2 / 16; deletion: w / N = 2 / 8
    table, targets = gt.call_targets(
        collapsed(
            snv("c1", 4, 0, 1), snv("c2", 2, 0, 1),
            deletion("c1", 4, 0, 1), deletion("c2", 2, 0, 1),
        ),
        error_floor=1e-4,
    )

    assert table["P_s"].tolist() == pytest.approx([1 / 8, 1 / 8, 1 / 4, 1 / 4])
    assert targets["P_s"].tolist() == pytest.approx([1 / 8, 1 / 4])


def test_output_table_fills_het_detect_lb_for_wt_calls():
    table, _ = gt.call_targets(
        collapsed(snv("c1", 3, 0, 0), snv("c2", 2, 1, 0), snv("c3", 1, 0, 0)),
        error_floor=1e-4,
    )

    assert table["genotype"].tolist()[:2] == ["WT", "MUT"]
    assert table["het_detect_lb"].tolist()[0] == 0.875
    assert np.isnan(table["het_detect_lb"].tolist()[1:]).all()


def test_call_targets_passes_min_mut_umis_to_the_caller():
    table, _ = gt.call_targets(
        collapsed(snv("c1", 2, 1, 0), snv("c2", 2, 2, 0)),
        error_floor=1e-4,
        min_mut_umis=2,
    )

    assert is_na(table["genotype"].tolist()[0])
    assert table["genotype"].tolist()[1] == "MUT"


def test_call_targets_does_not_mutate_its_input():
    c = collapsed(snv("c1", 1, 2, 0))
    before = c.copy()

    gt.call_targets(c, error_floor=1e-4)

    pd.testing.assert_frame_equal(c, before)


# genotype.log


def test_summary_reports_cells_umis_and_genes():
    table, _ = gt.call_targets(
        collapsed(
            snv("c1", 2, 1, 0),
            deletion("c1", 3, 0, 1),
            snv("c2", 0, 0, 0),
            deletion("c2", 1, 0, 0),
            snv("c3", 0, 0, 0),
        ),
        error_floor=1e-4,
    )

    # UMIs per cell: c1 = 7, c2 = 1, c3 = 0; genes covered: 2, 1, 0
    assert gt.summarise(table) == [
        "Total number of cells: 3",
        "Total number of UMIs: 8",
        "Median UMI per cell: 1.0",
        "Median genes covered per cell: 1.0",
        "Fraction of cells with UMIs: 0.6666666666666666",
    ]


# CLI


def run_cli(tmp_path, *extra):
    collapsed(
        snv("c1", 2, 3, 0),
        ("chrF", 1602, "T", "-T", "HP", 1, 1, 0, "c1"),
        snv("c2", 2, 1, 0),
    ).to_csv(tmp_path / "counts.tsv", index=False)

    subprocess.run(
        [
            sys.executable, str(Path(gt.__file__)),
            "--input", "counts.tsv",
            "--output", "S1.csv",
            "--targets_output", "S1_targets.csv",
            *extra,
        ],
        cwd=tmp_path,
        check=True,
    )
    return (
        (tmp_path / "S1.csv").read_text().splitlines(),
        (tmp_path / "S1_targets.csv").read_text().splitlines(),
    )


def test_cli_writes_the_sample_csv_without_an_index(tmp_path):
    lines, _ = run_cli(tmp_path)

    assert lines[0] == ",".join(OUTPUT_COLUMNS)
    assert lines[1].startswith("c1,chrF:301:A>G,SNV1,SNV,2,3,0,")
    assert lines[2].startswith("c1,chrF:1602:T>-T,HP,INDEL,1,1,0,")


def test_cli_writes_the_targets_csv_without_an_index(tmp_path):
    _, lines = run_cli(tmp_path)

    assert lines == [
        ",".join(TARGET_COLUMNS),
        "chrF:301:A>G,SNV1,SNV,8,0,2,0.0,0.0001,0.0001,True",
        "chrF:1602:T>-T,HP,INDEL,2,0,1,0.0,0.0001,0.0001,True",
    ]


def test_cli_min_mut_umis(tmp_path):
    lines, _ = run_cli(tmp_path, "--min_mut_umis", "2")

    genotypes = [line.split(",")[9] for line in lines[1:]]
    assert genotypes == ["MUT", "", ""]
