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
        help="Output name of the per-target error and diagnostics table (<sample>_targets.csv).",
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
    parser.add_argument(
        "-n",
        "--negative_cells",
        type=str,
        default=None,
        help="CSV of negative-control cells (columns sample, cell): cells that cannot carry the targeted mutations, e.g. T and B cells. Their MUT UMIs estimate the contamination rate c_s of each target. Without it c_s = 0.",
    )
    parser.add_argument(
        "-s",
        "--sample",
        type=str,
        default=None,
        help="Sample name, used to select this sample's rows of --negative_cells.",
    )

    args = parser.parse_args()
    return args


def read_negative_cells(path, sample=None):
    """
    Read the negative-control cell barcodes of a sample.

    Parameters:
        path (str or None): CSV with columns 'sample' and 'cell'. None means no
                            negative-control cells.
        sample (str or None): Sample whose rows are kept; None keeps every row.

    Returns:
        set: The barcodes of the negative-control cells (empty without a file).
    """
    if path is None:
        return set()
    negatives = pd.read_csv(path, dtype=str)
    missing = {"sample", "cell"} - set(negatives.columns)
    if missing:
        raise ValueError(f"--negative_cells is missing the columns {sorted(missing)}")
    if sample is not None:
        negatives = negatives[negatives["sample"] == sample]
    return set(negatives["cell"].dropna())


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
CONTAMINATION_COLUMNS = ["n_neg_cells", "n_neg_umis", "n_neg_mut", "c_s", "P_null"]


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


def target_contamination(table, targets, negative_cells):
    """
    Estimate the contamination of each target from its negative-control cells.

    A negative-control cell cannot carry the mutation, so every MUT UMI it
    shows came from elsewhere (ambient RNA, a PCR chimera, barcode swapping).
    Pooled over the negative-control cells of a target, c_s = n_neg_mut /
    n_neg_umis, the fraction of their labelled WT+MUT UMIs that read MUT, is
    the probability that a UMI of a cell without the mutation reads MUT by
    contamination; the OTH-based error q_hat cannot see it, because a
    contaminating molecule gives a clean MUT consensus. P_null = P_s + c_s is
    the per-UMI probability used by the caller. Without negative-control cells,
    or when they show no labelled UMI at the target, c_s = 0 and P_null = P_s,
    which is the genotype-v2 behaviour.

    Parameters:
        table (pd.DataFrame): Output of annotate.
        targets (pd.DataFrame): Output of target_error.
        negative_cells (set): Barcodes of the negative-control cells.

    Returns:
        pd.DataFrame: targets with CONTAMINATION_COLUMNS appended.
    """
    negative = table[table["cell"].isin(negative_cells)]
    pooled = (
        negative.assign(
            n_neg_cells=negative["WT"] + negative["MUT"] > 0,
            n_neg_umis=negative["WT"] + negative["MUT"],
            n_neg_mut=negative["MUT"],
        )
        .groupby("target_id", sort=False)[["n_neg_cells", "n_neg_umis", "n_neg_mut"]]
        .sum()
    )
    targets = targets.copy()
    for col in ["n_neg_cells", "n_neg_umis", "n_neg_mut"]:
        targets[col] = targets["target_id"].map(pooled[col]).fillna(0).astype(int)
    targets["c_s"] = np.where(
        targets["n_neg_umis"] > 0,
        targets["n_neg_mut"] / targets["n_neg_umis"].where(targets["n_neg_umis"] > 0),
        0.0,
    )
    targets["P_null"] = (targets["P_s"] + targets["c_s"]).clip(upper=1.0)
    return targets


# Per-cell diagnostic counts written by CONSENSUS, pooled over cells
DIAGNOSTIC_COLUMNS = [
    "wt_reads", "wt_altlike_reads", "mixed_discarded", "subnmin_umis", "subnmin_reads",
]
# The ALT-side noise diagnostic, defined at indel targets only
INDEL_DIAGNOSTIC_COLUMNS = ["wt_reads", "wt_altlike_reads", "wt_altlike_frac"]
CALL_COUNT_COLUMNS = ["n_MUT", "n_WT", "n_undet", "n_MUT_single_umi", "n_MUT_neg"]
ONE_UMI_COLUMNS = ["one_umi_mut_frac_neg", "one_umi_mut_frac_other"]
LIMITATION_COLUMNS = [
    "min_base_qual", "wt_reads", "wt_altlike_reads", "wt_altlike_frac",
    "mixed_discarded", "subnmin_umis", "subnmin_reads",
] + CALL_COUNT_COLUMNS + ONE_UMI_COLUMNS
TARGET_COLUMNS = ERROR_COLUMNS + CONTAMINATION_COLUMNS + LIMITATION_COLUMNS


def target_diagnostics(called, snv_min_base_quality=20, negative_cells=frozenset()):
    """
    Audit the limitations of each target, pooled over cells.

    min_base_qual is the base-quality filter of the target's bam-readcount pass:
    snv_min_base_quality at SNVs, 0 at indels. The DIAGNOSTIC_COLUMNS are summed
    over cells, and wt_altlike_frac = wt_altlike_reads / wt_reads (NaN when
    wt_reads = 0); the INDEL_DIAGNOSTIC_COLUMNS are NaN at SNVs. n_MUT, n_WT and
    n_undet count the calls over cells, undetermined including cells with no
    UMIs, n_MUT_single_umi the MUT calls with a single MUT UMI and n_MUT_neg the
    MUT calls in negative-control cells (contamination the caller still lets
    through). The one-UMI test: among cells with exactly one labelled WT or MUT
    UMI at the target, one_umi_mut_frac_neg is the fraction reading MUT in
    negative-control cells and one_umi_mut_frac_other in the other cells; a
    single MUT UMI carries information only where the second is well above the
    first. Each is NaN without such cells.

    Parameters:
        called (pd.DataFrame): Output of annotate with a 'genotype' column
                               (see call).
        snv_min_base_quality (int, optional): Base-quality filter at SNV targets.
                                              Defaults to 20.
        negative_cells (set, optional): Barcodes of the negative-control cells.
                                        Defaults to none.

    Returns:
        pd.DataFrame: One row per target, in input order, with 'target_id' and
                      LIMITATION_COLUMNS.
    """
    genotype = called["genotype"]
    mut = genotype == "MUT"
    negative = called["cell"].isin(negative_cells)
    one_umi = called["MUT"] + called["WT"] == 1
    targets = (
        called.assign(
            n_MUT=mut,
            n_WT=genotype == "WT",
            n_undet=genotype.isna(),
            n_MUT_single_umi=mut & (called["MUT"] == 1),
            n_MUT_neg=mut & negative,
            one_neg=one_umi & negative,
            one_neg_mut=one_umi & negative & (called["MUT"] == 1),
            one_other=one_umi & ~negative,
            one_other_mut=one_umi & ~negative & (called["MUT"] == 1),
        )
        .groupby("target_id", sort=False)
        .agg(type=("type", "first"),
             **{col: (col, "sum") for col in DIAGNOSTIC_COLUMNS + CALL_COUNT_COLUMNS
                + ["one_neg", "one_neg_mut", "one_other", "one_other_mut"]})
        .reset_index()
    )
    snv = targets["type"] == "SNV"
    targets["min_base_qual"] = np.where(snv, snv_min_base_quality, 0)
    targets["wt_altlike_frac"] = (
        targets["wt_altlike_reads"] / targets["wt_reads"].where(targets["wt_reads"] > 0)
    ).astype(float)
    targets["one_umi_mut_frac_neg"] = (
        targets["one_neg_mut"] / targets["one_neg"].where(targets["one_neg"] > 0)
    ).astype(float)
    targets["one_umi_mut_frac_other"] = (
        targets["one_other_mut"] / targets["one_other"].where(targets["one_other"] > 0)
    ).astype(float)

    # NaN at SNVs; the read counts stay integers as nullable Int64
    for col in INDEL_DIAGNOSTIC_COLUMNS:
        targets[col] = targets[col].mask(snv)
    targets[["wt_reads", "wt_altlike_reads"]] = targets[["wt_reads", "wt_altlike_reads"]].astype("Int64")
    return targets[["target_id"] + LIMITATION_COLUMNS]


def call(mut, wt, p_null, min_mut_umis=1):
    """
    Make calls from MUT and WT consensus counts.

    p = P(X >= m) with X ~ Binomial(m + k, P_null), the tail summed from j = 0,
    where P_null = P_s + c_s is the probability that a UMI of a cell without
    the mutation reads MUT by error or contamination (P_null = P_s without
    negative-control cells). MUT if m >= min_mut_umis and p < MUT_P_MAX; WT if
    m = 0 and k >= 2, or m >= min_mut_umis and p >= WT_P_MIN; undetermined
    (NaN) otherwise, including 0 < m < min_mut_umis. p is NaN when m + k = 0.

    Parameters:
        mut (np.ndarray): MUT consensus counts (m).
        wt (np.ndarray): WT consensus counts (k).
        p_null (np.ndarray): Null MUT probability P_null of each count pair.
        min_mut_umis (int, optional): Minimum MUT UMIs for a MUT call. Defaults to 1.

    Returns:
        Tuple[np.ndarray, np.ndarray]: p of each pair, and its call: 'WT', 'MUT'
                                       or NaN.
    """
    p = binom.sf(mut - 1, mut + wt, p_null)
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
    "P_s", "c_s", "P_null", "p", "genotype", "het_detect_lb", "negative_control",
]


def call_targets(counts, error_floor=1e-4, min_mut_umis=1, snv_min_base_quality=20,
                 negative_cells=frozenset()):
    """
    Call every cell x target with its target's pooled error and contamination.

    Parameters:
        counts (pd.DataFrame): Collapsed consensus counts with columns
                              [chr, start, ref, alt, gene, WT, MUT, OTH,
                              barcode] and the DIAGNOSTIC_COLUMNS.
        error_floor (float, optional): Floor (epsilon) for the error rate.
                                       Defaults to 1e-4.
        min_mut_umis (int, optional): Minimum MUT UMIs for a MUT call. Defaults to 1.
        snv_min_base_quality (int, optional): Base-quality filter at SNV targets,
                                              reported as min_base_qual. Defaults to 20.
        negative_cells (set, optional): Barcodes of the negative-control cells,
                                        from read_negative_cells. Defaults to
                                        none, so c_s = 0 at every target.

    Returns:
        Tuple[pd.DataFrame, pd.DataFrame]: <sample>.csv, one row per
            cell x target with OUTPUT_COLUMNS, and <sample>_targets.csv, one
            row per target with TARGET_COLUMNS.
    """
    table = annotate(counts)
    targets = target_contamination(table, target_error(table, error_floor), negative_cells)
    by_target = targets.set_index("target_id")
    for col in ["P_s", "c_s", "P_null"]:
        table[col] = table["target_id"].map(by_target[col])
    table["negative_control"] = table["cell"].isin(negative_cells)

    mut, wt = table["MUT"].values, table["WT"].values
    table["p"], table["genotype"] = call(mut, wt, table["P_null"].values, min_mut_umis)
    table["het_detect_lb"] = het_detect_lb(wt, table["genotype"].values)

    targets = targets.merge(
        target_diagnostics(table, snv_min_base_quality, negative_cells), on="target_id"
    )
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


def summarise_contamination(table, targets):
    """
    Summarise the negative-control cells and the contamination estimates for
    genotype.log: how many negative-control cells the sample has, the targets
    with c_s > 0 and their c_s, and the MUT calls left in negative-control
    cells. Without negative-control cells the lines say so and c_s is 0.

    Parameters:
        table (pd.DataFrame): <sample>.csv, from call_targets.
        targets (pd.DataFrame): <sample>_targets.csv, from call_targets.

    Returns:
        List[str]: One line per statistic.
    """
    n_negative = int(table.loc[table["negative_control"], "cell"].nunique())
    contaminated = targets[targets["c_s"] > 0]
    c_s_line = f"Targets with c_s > 0: {len(contaminated)} of {len(targets)}"
    if len(contaminated):
        c_s_line += ": " + ", ".join(
            f"{t} ({c:.4g})" for t, c in zip(contaminated["target_id"], contaminated["c_s"])
        )
    return [
        f"Negative-control cells: {n_negative}",
        c_s_line,
        f"MUT calls in negative-control cells: {int(targets['n_MUT_neg'].sum())}",
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
    n = targets[CALL_COUNT_COLUMNS].sum()
    wt = table[table["genotype"] == "WT"]
    floored = targets.loc[targets["floor_applied"], "target_id"].tolist()
    floor_line = (
        f"Targets at the error floor (eps = {targets['eps'].iloc[0]}): "
        f"{len(floored)} of {len(targets)}"
    )
    if floored:
        floor_line += ": " + ", ".join(floored)

    return [
        f"MUT calls: {n['n_MUT']}",
        f"WT calls: {n['n_WT']}",
        f"Undetermined calls: {n['n_undet']}",
        f"MUT calls with m = 1: {n['n_MUT_single_umi']}",
        f"Fraction of MUT calls with m = 1: {n['n_MUT_single_umi'] / n['n_MUT'] if n['n_MUT'] else np.nan}",
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

    negative_cells = read_negative_cells(args.negative_cells, args.sample)
    table, targets = call_targets(
        pd.read_csv(args.input), args.error_floor, args.min_mut_umis,
        args.snv_min_base_quality, negative_cells,
    )
    table.to_csv(args.output, index=False)
    targets.to_csv(args.targets_output, index=False)

    # collect statistics
    lines = (
        summarise(table)
        + summarise_calls(table, targets, args.min_mut_umis)
        + summarise_contamination(table, targets)
    )
    with open("genotype.log", "w") as f:
        for line in lines:
            print(line, file=f)


# Run
if __name__ == "__main__":
    main()
