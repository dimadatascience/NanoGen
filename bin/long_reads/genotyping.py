#!/usr/bin/python

import argparse

import numpy as np
import pandas as pd
from scipy.stats import binom

# Decision rule thresholds on p (genotyping.tex, section 4); not parameters
MUT_P_MAX = 0.1
WT_P_MIN = 0.9

# Error channels per target type: q_hat = w / (d * N)
ERROR_CHANNELS = {"SNV": 2, "INDEL": 1}


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
        help="Output name of the per cell x target calls (<sample>.csv).",
    )
    parser.add_argument(
        "-t",
        "--targets_output",
        type=str,
        required=True,
        help="Output name of the per-target error table (<sample>_targets.csv).",
    )
    parser.add_argument(
        "-b",
        "--error_floor",
        type=float,
        default=1e-4,
        help="Error floor (epsilon): lower bound on the per-target error rate.",
    )
    parser.add_argument(
        "-m",
        "--min_mut_umis",
        type=int,
        default=1,
        help="Minimum number of MUT UMIs (m_min) for a MUT call.",
    )

    args = parser.parse_args()
    return args


def target_id(chrom, start, ref, alt):
    """
    Identify a target as chr:start:ref>alt.

    Parameters:
        chrom (str): Chromosome.
        start (int): 1-based bam-readcount position, i.e. the BED start column.
        ref (str): Reference allele as written in the BED.
        alt (str): Alternative allele as written in the BED.

    Returns:
        str: The target ID, e.g. chr4:54733160:GAT>-GAT.
    """
    return f"{chrom}:{start}:{ref}>{alt}"


def variant_type(alt):
    """
    Type of a target from its alt.

    Parameters:
        alt (str): Alternative allele in bam-readcount notation.

    Returns:
        str: 'INDEL' for an insertion (+) or deletion (-) alt, 'SNV' otherwise.
    """
    return "INDEL" if alt.startswith(("+", "-")) else "SNV"


def annotate(counts):
    """
    Key collapsed consensus counts by cell and target.

    Parameters:
        counts (pd.DataFrame): Collapsed consensus counts with columns
                              [chr, start, ref, alt, gene, WT, MUT, OTH, barcode].

    Returns:
        pd.DataFrame: A copy of counts with 'barcode' renamed 'cell', and
                      'target_id' and 'type' columns.
    """
    table = counts.rename(columns={"barcode": "cell"})
    table["target_id"] = [
        target_id(chrom, start, ref, alt)
        for chrom, start, ref, alt in zip(
            table["chr"], table["start"], table["ref"], table["alt"]
        )
    ]
    table["type"] = [variant_type(alt) for alt in table["alt"]]
    return table


TARGET_COLUMNS = [
    "target_id", "gene", "type", "N", "w", "d", "q_hat", "eps", "P_s", "floor_applied",
]


def target_error(table, error_floor=1e-4):
    """
    Estimate the error of each target, pooled over cells.

    N is the number of labelled (WT, MUT or OTH) UMIs and w the number of OTH
    UMIs; q_hat = w / (d * N), with d = 2 at SNVs and 1 at indels, and
    P_s = max(q_hat, error_floor). No minimum N is required; q_hat is NaN when
    N = 0, and P_s is then the floor.

    Parameters:
        table (pd.DataFrame): Output of annotate.
        error_floor (float, optional): Floor (epsilon) for the error rate.
                                       Defaults to 1e-4.

    Returns:
        pd.DataFrame: One row per target, in input order, with TARGET_COLUMNS.
    """
    targets = (
        table.assign(N=table[["WT", "MUT", "OTH"]].sum(axis=1), w=table["OTH"])
        .groupby("target_id", sort=False)
        .agg(gene=("gene", "first"), type=("type", "first"), N=("N", "sum"), w=("w", "sum"))
        .reset_index()
    )
    targets["d"] = targets["type"].map(ERROR_CHANNELS)
    channel_umis = targets["d"] * targets["N"]
    targets["q_hat"] = (targets["w"] / channel_umis.where(channel_umis > 0)).astype(float)
    targets["eps"] = error_floor
    targets["floor_applied"] = targets["q_hat"].isna() | (targets["q_hat"] < error_floor)
    targets["P_s"] = targets["q_hat"].where(~targets["floor_applied"], error_floor)
    return targets[TARGET_COLUMNS]


def call(mut, wt, p_s, min_mut_umis=1):
    """
    Make calls from MUT and WT consensus counts.

    p = P(X >= m) with X ~ Binomial(m + k, P_s), the tail summed from j = 0.
    MUT if m >= min_mut_umis and p < MUT_P_MAX; WT if m = 0 and k >= 2, or
    m >= min_mut_umis and p >= WT_P_MIN; undetermined (NaN) otherwise,
    including 0 < m < min_mut_umis. p is NaN when m + k = 0.

    Parameters:
        mut (np.ndarray): MUT consensus counts (m).
        wt (np.ndarray): WT consensus counts (k).
        p_s (np.ndarray): Error rate P_s of each count pair.
        min_mut_umis (int, optional): Minimum MUT UMIs for a MUT call. Defaults to 1.

    Returns:
        Tuple[np.ndarray, np.ndarray]: p of each pair, and its call: 'WT', 'MUT'
                                       or NaN.
    """
    p = binom.sf(mut - 1, mut + wt, p_s)
    p[mut + wt == 0] = np.nan

    calls = np.full(p.size, np.nan, dtype=object)
    enough_mut = mut >= min_mut_umis
    calls[enough_mut & (p < MUT_P_MAX)] = "MUT"
    calls[((mut == 0) & (wt >= 2)) | (enough_mut & (p >= WT_P_MIN))] = "WT"
    return p, calls


def het_detect_lb(wt, calls):
    """
    Lower bound on the probability of detecting a heterozygous cell, 1 - 0.5^k.

    Parameters:
        wt (np.ndarray): WT consensus counts (k).
        calls (np.ndarray): Calls from call.

    Returns:
        np.ndarray: 1 - 0.5^k for WT calls, NaN otherwise.
    """
    return np.where(calls == "WT", 1 - 0.5 ** wt.astype(float), np.nan)


OUTPUT_COLUMNS = [
    "cell", "target_id", "gene", "type", "WT", "MUT", "OTH",
    "P_s", "p", "genotype", "het_detect_lb",
]


def call_targets(counts, error_floor=1e-4, min_mut_umis=1):
    """
    Call every cell x target with its target's pooled error.

    Parameters:
        counts (pd.DataFrame): Collapsed consensus counts with columns
                              [chr, start, ref, alt, gene, WT, MUT, OTH, barcode].
        error_floor (float, optional): Floor (epsilon) for the error rate.
                                       Defaults to 1e-4.
        min_mut_umis (int, optional): Minimum MUT UMIs for a MUT call. Defaults to 1.

    Returns:
        Tuple[pd.DataFrame, pd.DataFrame]: <sample>.csv, one row per
            cell x target with OUTPUT_COLUMNS, and <sample>_targets.csv, one
            row per target with TARGET_COLUMNS.
    """
    table = annotate(counts)
    targets = target_error(table, error_floor)
    table["P_s"] = table["target_id"].map(targets.set_index("target_id")["P_s"])

    mut, wt = table["MUT"].values, table["WT"].values
    table["p"], table["genotype"] = call(mut, wt, table["P_s"].values, min_mut_umis)
    table["het_detect_lb"] = het_detect_lb(wt, table["genotype"].values)
    return table[OUTPUT_COLUMNS].reset_index(drop=True), targets


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

    table, targets = call_targets(
        pd.read_csv(args.input), args.error_floor, args.min_mut_umis
    )
    table.to_csv(args.output, index=False)
    targets.to_csv(args.targets_output, index=False)

    # collect statistics
    with open("genotype.log", "w") as f:
        for line in summarise(table):
            print(line, file=f)


# Run
if __name__ == "__main__":
    main()
