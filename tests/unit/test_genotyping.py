# Characterisation tests: they pin current main behaviour of genotyping.py,
# defects included, so the pure-function refactor provably changes nothing.
# Later tickets replace them with spec tests.

import math
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import genotyping as gt


def counts(*rows):
    """Collapsed consensus counts: rows are (barcode, gene, MUT, WT, OTH)."""
    return pd.DataFrame(rows, columns=["barcode", "gene", "MUT", "WT", "OTH"])


# Per-cell error estimate


def test_cell_error_is_half_the_oth_fraction_over_all_targets_of_the_cell():
    table = gt.estimate_cell_error(
        counts(
            ("c1", "G1", 1, 2, 0),
            ("c1", "G2", 0, 3, 2),
            ("c2", "G1", 0, 4, 0),
        ).rename(columns={"barcode": "cell"})
    )

    # c1: 2 OTH of 8 UMIs, halved
    assert table["err"].tolist() == pytest.approx([0.125, 0.125, 0.0])


def test_cell_error_stays_zero_for_a_cell_without_oth():
    # The zero-replacement uses the overall min, which is then 0 itself
    table = gt.estimate_cell_error(
        counts(("c1", "G1", 0, 3, 1), ("c2", "G1", 0, 4, 0)).rename(
            columns={"barcode": "cell"}
        )
    )

    assert table["err"].tolist()[1] == 0.0


def test_cell_error_of_a_cell_without_umis_is_zero():
    table = gt.estimate_cell_error(
        counts(("c1", "G1", 0, 0, 0)).rename(columns={"barcode": "cell"})
    )

    assert table["err"].tolist() == [0.0]


# Caller


def call(mut, wt, err, error_floor=1e-4):
    p, called = gt.call(
        np.array([mut]), np.array([wt]), np.array([err]), error_floor
    )
    return p[0], called[0]


def is_na(x):
    return isinstance(x, float) and math.isnan(x)


def test_no_umis_is_na():
    p, called = call(0, 0, 0.0)
    assert math.isnan(p) and is_na(called)


def test_no_mut_umis_and_at_least_two_wt_is_wt():
    p, called = call(0, 2, 0.0)
    assert math.isnan(p) and called == "WT"


def test_no_mut_umis_and_one_wt_is_na():
    p, called = call(0, 1, 0.0)
    assert math.isnan(p) and is_na(called)


def test_p_is_the_negative_binomial_cdf_of_wt_given_mut():
    # n = 1: P(at most wt failures before 1 success) = 1 - (1 - err)^(wt + 1)
    p, _ = call(1, 2, 0.5)
    assert p == pytest.approx(1 - 0.5**3)


def test_p_below_0_1_is_mut():
    p, called = call(1, 2, 1e-4)
    assert p == pytest.approx(1 - (1 - 1e-4) ** 3)
    assert called == "MUT"


def test_p_at_least_0_9_is_wt():
    p, called = call(1, 30000, 1e-4)
    assert p == pytest.approx(1 - (1 - 1e-4) ** 30001)
    assert called == "WT"


def test_p_between_0_1_and_0_9_is_na():
    p, called = call(1, 5000, 1e-4)
    assert 0.1 < p < 0.9
    assert is_na(called)


def test_err_below_the_error_floor_is_raised_to_it():
    assert call(1, 2, 0.0, error_floor=0.5)[0] == pytest.approx(1 - 0.5**3)


def test_err_above_the_error_floor_is_used():
    assert call(1, 2, 0.5, error_floor=1e-4)[0] == pytest.approx(1 - 0.5**3)


# Output table and log


def test_call_table_adds_err_p_and_genotype_per_row():
    table = gt.call_targets(
        counts(("c1", "G1", 1, 2, 0), ("c1", "G2", 0, 3, 0), ("c2", "G1", 0, 1, 0)),
        error_floor=1e-4,
    )

    assert table.columns.tolist() == [
        "cell", "gene", "MUT", "WT", "OTH", "err", "p", "genotype",
    ]
    assert table["err"].tolist() == [0.0, 0.0, 0.0]
    assert table["genotype"].tolist()[:2] == ["MUT", "WT"]
    assert is_na(table["genotype"].tolist()[2])


def test_call_targets_does_not_mutate_its_input():
    c = counts(("c1", "G1", 1, 2, 0))

    gt.call_targets(c, error_floor=1e-4)

    assert c.columns.tolist() == ["barcode", "gene", "MUT", "WT", "OTH"]


def test_summary_reports_cells_umis_and_genes():
    table = gt.call_targets(
        counts(
            ("c1", "G1", 1, 2, 0),
            ("c1", "G2", 0, 3, 1),
            ("c2", "G1", 0, 0, 0),
            ("c2", "G2", 0, 1, 0),
            ("c3", "G1", 0, 0, 0),
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


# Output contract: <sample>.csv


def collapsed(*rows):
    """Collapsed consensus counts as COLLAPSE writes them: rows are
    (chr, start, ref, alt, gene, WT, MUT, OTH, barcode)."""
    table = pd.DataFrame(
        rows,
        columns=["chr", "start", "ref", "alt", "gene", "WT", "MUT", "OTH", "barcode"],
    )
    table.insert(2, "end", table["start"])
    return table


OUTPUT_COLUMNS = [
    "cell", "target_id", "gene", "type", "WT", "MUT", "OTH", "err", "p", "genotype",
]


def test_target_id_is_chr_position_ref_and_alt_as_in_the_bed():
    assert gt.target_id("chr4", 54733160, "GAT", "-GAT") == "chr4:54733160:GAT>-GAT"
    assert gt.target_id("chrF", 301, "A", "G") == "chrF:301:A>G"


@pytest.mark.parametrize(
    "ref, alt, expected",
    [
        ("A", "G", "SNV"),
        ("T", "-T", "INDEL"),
        ("T", "-TTA", "INDEL"),
        ("A", "+GTC", "INDEL"),
        ("GAT", "-GAT", "INDEL"),
    ],
)
def test_type_is_indel_for_bam_readcount_indel_alts_else_snv(ref, alt, expected):
    assert gt.variant_type(ref, alt) == expected


def test_output_table_has_the_v2_columns_one_row_per_cell_and_target():
    table = gt.output_table(
        gt.call_targets(
            collapsed(
                ("chrF", 301, "A", "G", "SNV1", 2, 3, 0, "c1"),
                ("chrF", 1002, "T", "-TTA", "DEL", 0, 0, 0, "c1"),
                ("chrF", 301, "A", "G", "SNV1", 0, 0, 0, "c2"),
            ),
            error_floor=1e-4,
        )
    )

    assert table.columns.tolist() == OUTPUT_COLUMNS
    assert table[["cell", "target_id", "gene", "type"]].values.tolist() == [
        ["c1", "chrF:301:A>G", "SNV1", "SNV"],
        ["c1", "chrF:1002:T>-TTA", "DEL", "INDEL"],
        ["c2", "chrF:301:A>G", "SNV1", "SNV"],
    ]
    assert table[["WT", "MUT", "OTH"]].values.tolist() == [[2, 3, 0], [0, 0, 0], [0, 0, 0]]
    assert table["genotype"].tolist()[0] == "MUT"


def test_cli_writes_the_sample_csv_without_an_index(tmp_path):
    collapsed(
        ("chrF", 301, "A", "G", "SNV1", 2, 3, 0, "c1"),
        ("chrF", 1602, "T", "-T", "HP", 1, 1, 0, "c1"),
    ).to_csv(tmp_path / "counts.tsv", index=False)

    subprocess.run(
        [
            sys.executable, str(Path(gt.__file__)),
            "--input", "counts.tsv", "--output", "S1.csv",
        ],
        cwd=tmp_path,
        check=True,
    )

    lines = (tmp_path / "S1.csv").read_text().splitlines()
    assert lines[0] == ",".join(OUTPUT_COLUMNS)
    assert lines[1].startswith("c1,chrF:301:A>G,SNV1,SNV,2,3,0,")
    assert lines[2].startswith("c1,chrF:1602:T>-T,HP,INDEL,1,1,0,")
