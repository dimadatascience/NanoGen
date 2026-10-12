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
        help="File containing consensus counts for each gene in the panel. Required columns: [chr, start, ref, alt, gene, WT, MUT, OTH, barcode]",
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
        "--error_floor",
        type=float,
        default=1e-4,
        help="Error floor (epsilon): lower bound on the per-target error rate.",
    )

    args = parser.parse_args()
    return args


def estimate_cell_error(barcodes):
    """
    Estimate the error rate of each cell from its OTH consensus counts.

    Parameters:
        barcodes (pd.DataFrame): Consensus counts with columns [MUT, WT, OTH, cell].

    Returns:
        pd.DataFrame: barcodes with an 'err' column: half the OTH fraction of
                      all the cell's UMIs.
    """
    grouped = barcodes.groupby("cell")
    barcodes = barcodes.merge(
        (
            (
                grouped["OTH"].sum()
                / (grouped[["MUT", "WT", "OTH"]].sum().sum(axis=1) + 1e-8)
                / 2
            )
            .to_frame("err")
            .reset_index()
        ),
        on="cell",
    )
    barcodes.loc[barcodes["err"] == 0, "err"] = barcodes["err"].min()
    return barcodes


def call(mut, wt, experimental_error, error_floor=1e-4):
    """
    Make calls from MUT and WT consensus counts.

    Parameters:
        mut (np.ndarray): MUT consensus counts.
        wt (np.ndarray): WT consensus counts.
        experimental_error (np.ndarray): Error rate of each count pair.
        error_floor (float, optional): Floor for the error rate. Defaults to 1e-4.

    Returns:
        Tuple[np.ndarray, np.ndarray]: The negative binomial p of each pair, and
                                       its call: 'WT', 'MUT' or NaN.
    """
    experimental_error = np.maximum(experimental_error, error_floor)

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


def call_targets(counts, error_floor=1e-4):
    """
    Make the call of every cell x target.

    Parameters:
        counts (pd.DataFrame): Consensus counts with columns [MUT, WT, OTH, barcode].
        error_floor (float, optional): Floor for the error rate. Defaults to 1e-4.

    Returns:
        pd.DataFrame: counts with 'barcode' renamed 'cell', and 'err', 'p' and
                      'genotype' columns.
    """
    barcodes = estimate_cell_error(counts.rename(columns={"barcode": "cell"}))
    p, calls = call(
        barcodes["MUT"].values,
        barcodes["WT"].values,
        barcodes["err"].values,
        error_floor,
    )
    barcodes["p"] = p
    barcodes["genotype"] = calls
    return barcodes


def target_id(chrom, start, ref, alt):
    """
    Identify a target as chr:start:ref>alt.

    start is the 1-based bam-readcount position, and ref and alt are as written
    in the BED, e.g. chr4:54733160:GAT>-GAT.
    """
    return f"{chrom}:{start}:{ref}>{alt}"


def variant_type(ref, alt):
    """
    Type of a target: 'INDEL' for a bam-readcount indel alt (+ins or -del),
    'SNV' otherwise.
    """
    return "INDEL" if alt.startswith(("+", "-")) else "SNV"


OUTPUT_COLUMNS = [
    "cell", "target_id", "gene", "type", "WT", "MUT", "OTH", "err", "p", "genotype",
]


def output_table(barcodes):
    """
    Shape the called table (see call_targets) as <sample>.csv.

    Returns:
        pd.DataFrame: One row per cell x target with OUTPUT_COLUMNS.
    """
    out = barcodes.copy()
    out["target_id"] = [
        target_id(*t) for t in out[["chr", "start", "ref", "alt"]].itertuples(index=False)
    ]
    out["type"] = [variant_type(r, a) for r, a in zip(out["ref"], out["alt"])]
    return out[OUTPUT_COLUMNS].reset_index(drop=True)


def summarise(barcodes):
    """
    Summarise the called table for genotype.log.

    Returns:
        List[str]: One line per statistic.
    """
    total_umis = barcodes["WT"] + barcodes["MUT"] + barcodes["OTH"]
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

    barcodes = output_table(call_targets(pd.read_csv(args.input), args.error_floor))
    barcodes.to_csv(args.output, index=False)

    # collect statistics
    with open("genotype.log", "w") as f:
        for line in summarise(barcodes):
            print(line, file=f)


# Run
if __name__ == "__main__":
    main()
