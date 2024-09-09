#!/usr/bin/python

import argparse
import sys

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


def main():

    args = _parse_args()

    barcodes = pd.read_csv(args.input).rename(columns={"barcode": "cell"})
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

    mut = barcodes["MUT"].values
    wt = barcodes["WT"].values
    experimental_error = barcodes["err"].values
    experimental_error = np.maximum(experimental_error, args.baseline_error)

    p = nbinom.cdf(k=wt, n=mut, p=experimental_error)
    barcodes["p"] = p

    genotypes = np.zeros(p.size)
    genotypes[(wt == 0) & (mut == 0)] = (
        np.nan
    )  # Real NAs: counts mut==0 and counts wt==0
    genotypes[p < 0.1] = 1  # MUT: p WT very low (<.1)
    genotypes[(np.isnan(p)) & (wt >= 3)] = (
        0  # WT: p NA for counts mut == 0 and count wt>=3
    )
    genotypes[p >= 0.9] = 0  # WT: p WT very high (>=.9)
    genotypes[(np.isnan(p)) & (wt < 3)] = (
        np.nan
    )  # NAs (impossible to assign): p NA for counts mut == 0 but count wt<3
    genotypes[(p >= 0.1) & (p < 0.9)] = (
        np.nan
    )  # NAs (impossible to assign): p WT between .1 and .9

    barcodes["genotype"] = np.select(
        [genotypes == 0, genotypes == 1], ["WT", "MUT"], default="-"
    )
    barcodes.loc[barcodes["genotype"] == "-", "genotype"] = np.nan
    barcodes.to_csv(args.output)

    # collect statistics
    barcodes["total_umis"] = barcodes["WT"] + barcodes["MUT"] + barcodes["MIS"]
    numis = barcodes.groupby("cell")["total_umis"].sum()
    ngenes = barcodes.groupby("cell")["total_umis"].apply(lambda x: np.sum(x > 0))
    ncells = len(barcodes["cell"].unique())

    # Open a file to write the output
    with open("genotype.log", "w") as f:
        # Redirect stdout to the file
        sys.stdout = f

        # Your print statements
        print(f"Total number of cells: {ncells}")
        print(f"Total number of UMIs: {np.sum(numis)}")
        print(f"Median UMI per cell: {np.median(numis)}")
        print(f"Median genes covered per cell: {np.median(ngenes)}")
        print(f"Fraction of cells with UMIs: {np.mean(numis>0)}")

        # Restore stdout to its original state
        sys.stdout = sys.__stdout__


# Run
if __name__ == "__main__":
    main()
