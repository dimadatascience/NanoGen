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
        help="Barcode of the cell.",
    )
    args = parser.parse_args()
    return args


##


def consensus_on_umi(line, ref, alt, min_fraction=0.75, min_consensus=3):
    """
    Create consensus based on bam-readcount output where each read group (e.g. same UMI tag)
    has been assigned to library ID in RG within bam.

    The function extracts data blocks within curly braces (counts of each read group), processes each
    base count within these blocks, and evaluates whether the reference or alternative base
    meets consensus criteria based on minimum fraction and minimum consensus read counts.

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
    # Regex to extract data within curly braces
    data_blocks = re.findall(r"\{([^}]*)\}", line)

    results = []
    consensus_counts = {"WT": 0, "MUT": 0, "MIS": 0}
    for block in data_blocks:
        counts = {ref: 0, alt: 0, "mis": 0}
        # Split block into components, each corresponding to a base
        entries = block.split()

        # Process each entry
        for entry in entries:
            parts = entry.split(":")
            base = parts[0]
            read_count = int(parts[1])

            if base == ref:
                counts[ref] = read_count
            elif base == alt:
                counts[alt] = read_count
            elif base != "N":
                counts["mis"] += read_count

        maximum = np.max([counts[ref], counts[alt]])
        count_values = np.array(list(counts.values()))

        if maximum >= min_consensus and maximum / np.sum(count_values) >= min_fraction:
            key = np.array(list(counts.keys()))[count_values == np.max(count_values)]
            if key == ref:
                consensus_counts["WT"] += 1
            elif key == alt:
                consensus_counts["MUT"] += 1
            else:
                consensus_counts["MIS"] += 1
    return consensus_counts


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
    target[["WT", "MUT", "MIS"]] = 0

    # Open the file using a with statement
    with open(args.input, "r") as file:
        # Read each line in the file one by one
        for line in file:
            chrom = line.split("\t")[0]
            start_pos = int(line.split("\t")[1])
            ref = line.split("\t")[2]
            alts = target["alt"][
                (target["chr"] == chrom)
                & (target["start"] == start_pos)
                & (target["ref"] == ref)
            ]
            for alt in alts:
                for val in ["WT", "MUT", "MIS"]:
                    target.loc[
                        (target["chr"] == chrom)
                        & (target["start"] == start_pos)
                        & (target["alt"] == alt),
                        val,
                    ] = consensus_on_umi(line, ref, alt)[val]

    target["barcode"] = args.cell_barcode
    target.to_csv(f"{args.cell_barcode}_table.csv", index=False)


##


# Run
if __name__ == "__main__":
    main()
