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
        help="File containing consensus counts for each gene in the panel. Required columns: [chr, start, ref, alt, gene, WT, MUT, OTH, wt_reads, wt_altlike_reads, mixed_discarded, subnmin_umis, subnmin_reads, barcode]",
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
    parser.add_argument(
        "-q",
        "--snv_min_base_quality",
        type=int,
        default=20,
        help="Base-quality filter of the SNV bam-readcount pass, reported as min_base_qual.",
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


ERROR_COLUMNS = [
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
        pd.DataFrame: One row per target, in input order, with ERROR_COLUMNS.
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
    return targets[ERROR_COLUMNS]


# Per-cell diagnostic counts written by CONSENSUS, pooled over cells
DIAGNOSTIC_COLUMNS = [
    "wt_reads", "wt_altlike_reads", "mixed_discarded", "subnmin_umis", "subnmin_reads",
]
# The ALT-side noise diagnostic, defined at indel targets only
INDEL_DIAGNOSTIC_COLUMNS = ["wt_reads", "wt_altlike_reads", "wt_altlike_frac"]
CALL_COUNT_COLUMNS = ["n_MUT", "n_WT", "n_undet", "n_MUT_single_umi"]
TARGET_COLUMNS = ERROR_COLUMNS + [
    "min_base_qual", "wt_reads", "wt_altlike_reads", "wt_altlike_frac",
    "mixed_discarded", "subnmin_umis", "subnmin_reads",
] + CALL_COUNT_COLUMNS


def target_diagnostics(called, snv_min_base_quality=20):
    """
    Audit the limitations of each target, pooled over cells.

    min_base_qual is the base-quality filter of the target's bam-readcount pass:
    snv_min_base_quality at SNVs, 0 at indels. The DIAGNOSTIC_COLUMNS are summed
    over cells, and wt_altlike_frac = wt_altlike_reads / wt_reads (NaN when
    wt_reads = 0); the INDEL_DIAGNOSTIC_COLUMNS are NaN at SNVs. n_MUT, n_WT and
    n_undet count the calls over cells, undetermined including cells with no
    UMIs, and n_MUT_single_umi the MUT calls with a single MUT UMI.

    Parameters:
        called (pd.DataFrame): Output of annotate with a 'genotype' column
                               (see call).
        snv_min_base_quality (int, optional): Base-quality filter at SNV targets.
                                              Defaults to 20.

    Returns:
        pd.DataFrame: One row per target, in input order, with 'target_id', then
                      the columns of TARGET_COLUMNS after ERROR_COLUMNS.
    """
    genotype = called["genotype"]
    mut = genotype == "MUT"
    targets = (
        called.assign(
            n_MUT=mut,
            n_WT=genotype == "WT",
            n_undet=genotype.isna(),
            n_MUT_single_umi=mut & (called["MUT"] == 1),
        )
        .groupby("target_id", sort=False)
        .agg(type=("type", "first"),
             **{col: (col, "sum") for col in DIAGNOSTIC_COLUMNS + CALL_COUNT_COLUMNS})
        .reset_index()
    )
    targets["min_base_qual"] = np.where(targets["type"] == "SNV", snv_min_base_quality, 0)
    targets["wt_altlike_frac"] = (
        targets["wt_altlike_reads"] / targets["wt_reads"].where(targets["wt_reads"] > 0)
    ).astype(float)

    # NaN at SNVs; the read counts stay integers as nullable Int64
    snv = targets["type"] == "SNV"
    for col in INDEL_DIAGNOSTIC_COLUMNS:
        targets[col] = targets[col].mask(snv)
    targets[["wt_reads", "wt_altlike_reads"]] = targets[["wt_reads", "wt_altlike_reads"]].astype("Int64")
    return targets[["target_id"] + TARGET_COLUMNS[len(ERROR_COLUMNS):]]


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


def call_targets(counts, error_floor=1e-4, min_mut_umis=1, snv_min_base_quality=20):
    """
    Call every cell x target with its target's pooled error.

    Parameters:
        counts (pd.DataFrame): Collapsed consensus counts with columns
                              [chr, start, ref, alt, gene, WT, MUT, OTH,
                              barcode] and the DIAGNOSTIC_COLUMNS.
        error_floor (float, optional): Floor (epsilon) for the error rate.
                                       Defaults to 1e-4.
        min_mut_umis (int, optional): Minimum MUT UMIs for a MUT call. Defaults to 1.
        snv_min_base_quality (int, optional): Base-quality filter at SNV targets,
                                              reported as min_base_qual. Defaults to 20.

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

    targets = targets.merge(target_diagnostics(table, snv_min_base_quality), on="target_id")
    return table[OUTPUT_COLUMNS].reset_index(drop=True), targets[TARGET_COLUMNS]


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


def summarise_calls(table, targets, min_mut_umis=1):
    """
    Summarise the calls and their limitations for genotype.log: call totals,
    single-UMI MUT calls, WT calls by k, het_detect_lb over WT calls, the
    targets at the error floor and min_mut_umis. Statistics over no calls are nan.

    Parameters:
        table (pd.DataFrame): <sample>.csv, from call_targets.
        targets (pd.DataFrame): <sample>_targets.csv, from call_targets.
        min_mut_umis (int, optional): m_min used for the calls. Defaults to 1.

    Returns:
        List[str]: One line per statistic.
    """
    mut = table[table["genotype"] == "MUT"]
    wt = table[table["genotype"] == "WT"]
    single_umi = int((mut["MUT"] == 1).sum())
    floored = targets.loc[targets["floor_applied"], "target_id"].tolist()
    floor_line = (
        f"Targets at the error floor (eps = {targets['eps'].iloc[0]}): "
        f"{len(floored)} of {len(targets)}"
    )
    if floored:
        floor_line += ": " + ", ".join(floored)

    return [
        f"MUT calls: {len(mut)}",
        f"WT calls: {len(wt)}",
        f"Undetermined calls: {int(table['genotype'].isna().sum())}",
        f"MUT calls with m = 1: {single_umi}",
        f"Fraction of MUT calls with m = 1: {single_umi / len(mut) if len(mut) else np.nan}",
        f"WT calls with k = 2: {int((wt['WT'] == 2).sum())}",
        f"WT calls with k = 3: {int((wt['WT'] == 3).sum())}",
        f"WT calls with k >= 4: {int((wt['WT'] >= 4).sum())}",
        f"Median het_detect_lb over WT calls: {wt['het_detect_lb'].median()}",
        f"Minimum het_detect_lb over WT calls: {wt['het_detect_lb'].min()}",
        floor_line,
        f"min_mut_umis: {min_mut_umis}",
    ]


def main():

    args = _parse_args()

    table, targets = call_targets(
        pd.read_csv(args.input), args.error_floor, args.min_mut_umis,
        args.snv_min_base_quality,
    )
    table.to_csv(args.output, index=False)
    targets.to_csv(args.targets_output, index=False)

    # collect statistics
    with open("genotype.log", "w") as f:
        for line in summarise(table) + summarise_calls(table, targets, args.min_mut_umis):
            print(line, file=f)


# Run
if __name__ == "__main__":
    main()
