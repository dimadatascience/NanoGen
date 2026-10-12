# Characterisation tests: they pin current main behaviour of genotyping.py,
# defects included, so the pure-function refactor provably changes nothing.
# Later tickets replace them with spec tests.

import math

import numpy as np
import pandas as pd
import pytest

import genotyping as gt


def counts(*rows):
    """Collapsed consensus counts: rows are (barcode, gene, MUT, WT, MIS)."""
    return pd.DataFrame(rows, columns=["barcode", "gene", "MUT", "WT", "MIS"])


# Per-cell error estimate


def test_cell_error_is_half_the_mis_fraction_over_all_targets_of_the_cell():
    table = gt.estimate_cell_error(
        counts(
            ("c1", "G1", 1, 2, 0),
            ("c1", "G2", 0, 3, 2),
            ("c2", "G1", 0, 4, 0),
        ).rename(columns={"barcode": "cell"})
    )

    # c1: 2 MIS of 8 UMIs, halved; the 1e-8 guard shifts it negligibly
    assert table["err"].tolist() == pytest.approx([0.125, 0.125, 0.0])
    assert table["err"].iloc[0] < 0.125


def test_cell_error_stays_zero_for_a_cell_without_mis():
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


def call(mut, wt, err, baseline_error=1e-4):
    p, genotype = gt.call_genotypes(
        np.array([mut]), np.array([wt]), np.array([err]), baseline_error
    )
    return p[0], genotype[0]


def is_na(x):
    return isinstance(x, float) and math.isnan(x)


def test_no_umis_is_na():
    p, genotype = call(0, 0, 0.0)
    assert math.isnan(p) and is_na(genotype)


def test_no_mut_umis_and_at_least_two_wt_is_wt():
    p, genotype = call(0, 2, 0.0)
    assert math.isnan(p) and genotype == "WT"


def test_no_mut_umis_and_one_wt_is_na():
    p, genotype = call(0, 1, 0.0)
    assert math.isnan(p) and is_na(genotype)


def test_p_is_the_negative_binomial_cdf_of_wt_given_mut():
    # n = 1: P(at most wt failures before 1 success) = 1 - (1 - err)^(wt + 1)
    p, _ = call(1, 2, 0.5)
    assert p == pytest.approx(1 - 0.5**3)


def test_p_below_0_1_is_mut():
    p, genotype = call(1, 2, 1e-4)
    assert p == pytest.approx(1 - (1 - 1e-4) ** 3)
    assert genotype == "MUT"


def test_p_at_least_0_9_is_wt():
    p, genotype = call(1, 30000, 1e-4)
    assert p == pytest.approx(1 - (1 - 1e-4) ** 30001)
    assert genotype == "WT"


def test_p_between_0_1_and_0_9_is_na():
    p, genotype = call(1, 5000, 1e-4)
    assert 0.1 < p < 0.9
    assert is_na(genotype)


def test_err_below_the_baseline_is_raised_to_it():
    assert call(1, 2, 0.0, baseline_error=0.5)[0] == pytest.approx(1 - 0.5**3)


def test_err_above_the_baseline_is_used():
    assert call(1, 2, 0.5, baseline_error=1e-4)[0] == pytest.approx(1 - 0.5**3)


# Output table and log


def test_genotype_table_adds_err_p_and_genotype_per_row():
    table = gt.genotype(
        counts(("c1", "G1", 1, 2, 0), ("c1", "G2", 0, 3, 0), ("c2", "G1", 0, 1, 0)),
        baseline_error=1e-4,
    )

    assert table.columns.tolist() == [
        "cell", "gene", "MUT", "WT", "MIS", "err", "p", "genotype",
    ]
    assert table["err"].tolist() == [0.0, 0.0, 0.0]
    assert table["genotype"].tolist()[:2] == ["MUT", "WT"]
    assert is_na(table["genotype"].tolist()[2])


def test_genotype_does_not_mutate_its_input():
    c = counts(("c1", "G1", 1, 2, 0))

    gt.genotype(c, baseline_error=1e-4)

    assert c.columns.tolist() == ["barcode", "gene", "MUT", "WT", "MIS"]


def test_summary_reports_cells_umis_and_genes():
    table = gt.genotype(
        counts(
            ("c1", "G1", 1, 2, 0),
            ("c1", "G2", 0, 3, 1),
            ("c2", "G1", 0, 0, 0),
            ("c2", "G2", 0, 1, 0),
            ("c3", "G1", 0, 0, 0),
        ),
        baseline_error=1e-4,
    )

    # UMIs per cell: c1 = 7, c2 = 1, c3 = 0; genes covered: 2, 1, 0
    assert gt.summarise(table) == [
        "Total number of cells: 3",
        "Total number of UMIs: 8",
        "Median UMI per cell: 1.0",
        "Median genes covered per cell: 1.0",
        "Fraction of cells with UMIs: 0.6666666666666666",
    ]
