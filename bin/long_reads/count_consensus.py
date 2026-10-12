#!/usr/bin/env python

import argparse
import os
import re

import numpy as np
import pandas as pd


##


# Utils
def _parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        default=None,
        required=True,
        help="Bam-readcount output file run on multiple libraries.",
    )
    parser.add_argument(
        "-t",
        "--target",
        type=str,
        default=None,
        required=True,
        help="Tab delimited file with target positions in the following format: ['chr', 'start', 'end', 'ref', 'alt', 'gene']",
    )
    parser.add_argument(
        "-n",
        "--cell_barcode",
        type=str,
        required=True,
        help="Barcode of the cell.",
    )
    parser.add_argument(
        "-mf",
        "--min_fraction",
        type=float,
        default=0.75,
        required=False,
        help="The minimum fraction of reads supporting a base.",
    )
    parser.add_argument(
        "-mr",
        "--min_read",
        type=int,
        default=3,
        required=False,
        help="The minimum number of reads supporting a base required to consider it for consensus.",
    )
    args = parser.parse_args()
    return args


##


# bam-readcount parsing
def parse_clusters(line):
    """
    Extract per-cluster base counts from one bam-readcount line.

    bam-readcount -p writes one data block within curly braces per read group
    (here, per UMI cluster). Without curly braces the whole line is read as a
    single block.

    Parameters:
        line (str): A line of bam-readcount output.

    Returns:
        List[List[Tuple[str, int]]]: For each block, the (base, read count) of
                                     each entry, in file order.
    """
    # Regex to extract data within curly braces
    if "{" in line:
        data_blocks = re.findall(r"\{([^}]*)\}", line)
    else:
        data_blocks = [line]

    clusters = []
    for block in data_blocks:
        entries = []
        # Split block into components, each corresponding to a base
        for entry in block.split():
            if ":" in entry:
                parts = entry.split(":")
                entries.append((parts[0], int(parts[1])))
        clusters.append(entries)
    return clusters


def parse_readcount_line(line):
    """
    Parse one bam-readcount -p line.

    Returns:
        Tuple[str, int, str, List[List[Tuple[str, int]]]]: Chromosome, position,
                                                          reference base and the
                                                          per-cluster base counts.
    """
    fields = line.split("\t")
    return fields[0], int(fields[1]), fields[2], parse_clusters(line)


##


# Consensus
def label_cluster(entries, ref, alt, min_fraction=0.75, min_consensus=3):
    """
    Assign a consensus label to one UMI cluster.

    Evaluates whether the reference or alternative base meets consensus
    criteria based on minimum fraction and minimum consensus read counts.

    Parameters:
        entries (Iterable[Tuple[str, int]]): (base, read count) of the cluster.
        ref (str): The reference base to be evaluated for consensus.
        alt (str): The alternative base to be evaluated for consensus.
        min_fraction (float, optional): The minimum fraction of reads supporting a base
                                        necessary to declare a consensus. Defaults to 0.75.
        min_consensus (int, optional): The minimum number of reads supporting a base required
                                       to consider it for consensus. Defaults to 3.

    Returns:
        str or None: 'WT', 'MUT' or 'MIS', or None if the cluster gets no label.
    """
    counts = {ref: 0, alt: 0, "mis": 0}
    for base, read_count in entries:
        if base == ref:
            counts[ref] = read_count
        elif len(alt) <= 2 and base == alt:
            counts[alt] = read_count
        elif len(alt) > 2 and len(base) > 2 and (base in alt or alt in base):
            counts[alt] = np.maximum(read_count, counts[alt])
        elif base != "N":
            counts["mis"] += read_count

    maximum = np.max([counts[ref], counts[alt]])
    count_values = np.array(list(counts.values()))
    key = np.array(list(counts.keys()))[count_values == np.max(count_values)]

    if len(alt) == 2 and counts[alt] > 2 and counts[alt]/(counts[ref] + counts[alt]) > 0.25:
        return "MUT"
    elif len(alt) > 3 and counts[alt] > 3:
        return "MUT"
    elif (
        maximum >= min_consensus and maximum / np.sum(count_values) >= min_fraction
    ):
        if key == ref:
            return "WT"
        elif key == alt:
            return "MUT"
        else:
            return "MIS"
    return None


def consensus_on_umi(line, ref, alt, min_fraction=0.75, min_consensus=3):
    """
    Create consensus based on bam-readcount output where each read group (e.g. same UMI tag)
    has been assigned to library ID in RG within bam.

    Parameters:
        line (str): A string from a sequencing data file which contains multiple data blocks
                    enclosed in curly braces. Each block includes base counts.
        ref (str): The reference base to be evaluated for consensus.
        alt (str): The alternative base to be evaluated for consensus.
        min_fraction (float, optional): The minimum fraction of reads supporting a base
                                        necessary to declare a consensus. Defaults to 0.75.
        min_consensus (int, optional): The minimum number of reads supporting a base required
                                       to consider it for consensus. Defaults to 3.

    Returns:
        Dict: A dictionary containing counts of consensus results for each data
                       block, divided by consensus types ('WT', 'MUT', 'MIS').
    """
    consensus_counts = {"WT": 0, "MUT": 0, "MIS": 0}
    for entries in parse_clusters(line):
        label = label_cluster(entries, ref, alt, min_fraction, min_consensus)
        if label is not None:
            consensus_counts[label] += 1
    return consensus_counts


def count_targets(target, lines, min_fraction=0.75, min_read=3):
    """
    Count consensus labels per target for one cell.

    Parameters:
        target (pd.DataFrame): Targets with columns ['chr', 'start', 'end', 'ref', 'alt', 'gene'].
        lines (Iterable[str]): bam-readcount -p output lines.
        min_fraction (float, optional): See label_cluster. Defaults to 0.75.
        min_read (int, optional): See label_cluster (min_consensus). Defaults to 3.

    Returns:
        pd.DataFrame: A copy of target with 'WT', 'MUT' and 'MIS' consensus counts.
    """
    target = target.copy()
    target[["WT", "MUT", "MIS"]] = 0

    for line in lines:
        chrom, start_pos, ref, _ = parse_readcount_line(line)
        alts = target["alt"][
            (target["chr"] == chrom)
            & (target["start"] == start_pos)
            & (target["ref"] == ref)
        ]
        for alt in alts:
            counts = consensus_on_umi(line, ref, alt, min_fraction, min_read)
            for val in ["WT", "MUT", "MIS"]:
                target.loc[
                    (target["chr"] == chrom)
                    & (target["start"] == start_pos)
                    & (target["alt"] == alt),
                    val,
                ] = counts[val]
    return target


##


def main():
    """
    Generate count matrix based of consensus base at positions within target position.
    """

    # Parse input
    args = _parse_args()

    if not os.path.exists(args.input):
        raise ValueError(f"bam-readcount output file {args.input} does not exist!")
    if not os.path.exists(args.target):
        raise ValueError(f"bam-readcount output file {args.target} does not exist!")

    target = pd.read_csv(args.target, sep="\t", header=None)
    target.columns = ["chr", "start", "end", "ref", "alt", "gene"]

    with open(args.input, "r") as file:
        target = count_targets(target, file, args.min_fraction, args.min_read)

    target["barcode"] = args.cell_barcode
    target.to_csv(f"{args.cell_barcode}_table.csv", index=False)




##


# Run
if __name__ == "__main__":
    main()
