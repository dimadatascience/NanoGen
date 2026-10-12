#!/usr/bin/python

import argparse

import numpy as np
import pandas as pd
from scipy.stats import nbinom


def _parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        default=None,
        required=True,
        help="File containing consensus counts for each gene in the panel. Required columns: [MUT, WT, MIS, barcode]",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        required=True,
        help="Output name.",
    )
    parser.add_argument(
        "-b",
        "--baseline_error",
        type=float,
        default=1e-4,
        help="Baseline error accounting for any form of errors. For instance contamination.",
    )

    args = parser.parse_args()
    return args


def estimate_cell_error(barcodes):
    """
    Estimate the error rate of each cell from its MIS consensus counts.

    Parameters:
        barcodes (pd.DataFrame): Consensus counts with columns [MUT, WT, MIS, cell].

    Returns:
        pd.DataFrame: barcodes with an 'err' column: half the MIS fraction of
                      all the cell's UMIs.
    """
    grouped = barcodes.groupby("cell")
    barcodes = barcodes.merge(
        (
            (
                grouped["MIS"].sum()
                / (grouped[["MUT", "WT", "MIS"]].sum().sum(axis=1) + 1e-8)
                / 2
            )
            .to_frame("err")
            .reset_index()
        ),
        on="cell",
    )
    barcodes.loc[barcodes["err"] == 0, "err"] = barcodes["err"].min()
    return barcodes


def call(mut, wt, experimental_error, baseline_error=1e-4):
    """
    Make calls from MUT and WT consensus counts.

    Parameters:
        mut (np.ndarray): MUT consensus counts.
        wt (np.ndarray): WT consensus counts.
        experimental_error (np.ndarray): Error rate of each count pair.
        baseline_error (float, optional): Floor for the error rate. Defaults to 1e-4.

    Returns:
        Tuple[np.ndarray, np.ndarray]: The negative binomial p of each pair, and
                                       its call: 'WT', 'MUT' or NaN.
    """
    experimental_error = np.maximum(experimental_error, baseline_error)

    p = nbinom.cdf(k=wt, n=mut, p=experimental_error)

    codes = np.zeros(p.size)
    codes[(wt == 0) & (mut == 0)] = (
        np.nan
    )  # Real NAs: counts mut==0 and counts wt==0
    codes[p < 0.1] = 1  # MUT: p WT very low (<.1)
    codes[(np.isnan(p)) & (wt >= 2)] = (
        0  # WT: p NA for counts mut == 0 and count wt>=3
    )
    codes[p >= 0.9] = 0  # WT: p WT very high (>=.9)
    codes[(np.isnan(p)) & (wt < 2)] = (
        np.nan
    )  # NAs (impossible to assign): p NA for counts mut == 0 but count wt<3
    codes[(p >= 0.1) & (p < 0.9)] = (
        np.nan
    )  # NAs (impossible to assign): p WT between .1 and .9

    calls = np.select(
        [codes == 0, codes == 1], ["WT", "MUT"], default="-"
    ).astype(object)
    calls[calls == "-"] = np.nan
    return p, calls


def call_targets(counts, baseline_error=1e-4):
    """
    Make the call of every cell x target.

    Parameters:
        counts (pd.DataFrame): Consensus counts with columns [MUT, WT, MIS, barcode].
        baseline_error (float, optional): Floor for the error rate. Defaults to 1e-4.

    Returns:
        pd.DataFrame: counts with 'barcode' renamed 'cell', and 'err', 'p' and
                      'genotype' columns.
    """
    barcodes = estimate_cell_error(counts.rename(columns={"barcode": "cell"}))
    p, calls = call(
        barcodes["MUT"].values,
        barcodes["WT"].values,
        barcodes["err"].values,
        baseline_error,
    )
    barcodes["p"] = p
    barcodes["genotype"] = calls
    return barcodes


def summarise(barcodes):
    """
    Summarise the called table for genotype.log.

    Returns:
        List[str]: One line per statistic.
    """
    total_umis = barcodes["WT"] + barcodes["MUT"] + barcodes["MIS"]
    numis = total_umis.groupby(barcodes["cell"]).sum()
    ngenes = total_umis.groupby(barcodes["cell"]).apply(lambda x: np.sum(x > 0))
    ncells = len(barcodes["cell"].unique())

    return [
        f"Total number of cells: {ncells}",
        f"Total number of UMIs: {np.sum(numis)}",
        f"Median UMI per cell: {np.median(numis)}",
        f"Median genes covered per cell: {np.median(ngenes)}",
        f"Fraction of cells with UMIs: {np.mean(numis>0)}",
    ]


def main():

    args = _parse_args()

    barcodes = call_targets(pd.read_csv(args.input), args.baseline_error)
    barcodes.to_csv(args.output)

    # collect statistics
    with open("genotype.log", "w") as f:
        for line in summarise(barcodes):
            print(line, file=f)


# Run
if __name__ == "__main__":
    main()
