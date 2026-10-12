#!/usr/bin/env python

import argparse
import os
import re

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
        default=0.6,
        required=False,
        help="f_min: the minimum fraction of a cluster's usable reads supporting its label.",
    )
    parser.add_argument(
        "-mr",
        "--min_read",
        type=int,
        default=3,
        required=False,
        help="n_min: the minimum number of usable reads supporting a cluster's label.",
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
#
# bam-readcount -p lists a block for every UMI cluster with reads in the pileup
# at a site, so every block is a cluster covering the target. Each covering
# cluster is exactly one of: sub-n_min, labelled (WT, MUT or OTH) or discarded.

SUBNMIN = "sub-n_min"
DISCARDED = "discarded"
COUNT_COLUMNS = ["WT", "MUT", "OTH", "subnmin_umis", "subnmin_reads"]
BASES = ("A", "C", "G", "T")


def allele_counts(entries, ref, alt):
    """
    Count the usable reads of one UMI cluster as REF, ALT and OTH.

    At an SNV, a usable read carries A, C, G or T at the site. N reads and reads
    with a deletion over the site (reported as '-' entries) are not usable, and
    '+' entries are skipped because their reads are already in their base count.
    Base quality is filtered upstream, by bam-readcount -b.

    At an indel target, the legacy matching is kept until tolerant indel ALT
    replaces it: REF is the reference base, ALT the exact alt (len(alt) <= 2)
    or the largest entry overlapping it by substring, OTH every other non-N entry.

    Parameters:
        entries (Iterable[Tuple[str, int]]): (base, read count) of the cluster.
        ref (str): The reference base.
        alt (str): The alternative allele, in bam-readcount notation.

    Returns:
        Tuple[int, int, int]: r, a and o, the REF, ALT and OTH usable reads.
    """
    if alt.startswith(("+", "-")):
        return _indel_allele_counts(entries, ref, alt)

    r = a = o = 0
    for base, read_count in entries:
        if base == ref:
            r += read_count
        elif base == alt:
            a += read_count
        elif base in BASES:
            o += read_count
    return r, a, o


def _indel_allele_counts(entries, ref, alt):
    """Legacy indel matching; see allele_counts."""
    r = a = o = 0
    for base, read_count in entries:
        if base == ref:
            r = read_count
        elif len(alt) <= 2 and base == alt:
            a = read_count
        elif len(alt) > 2 and len(base) > 2 and (base in alt or alt in base):
            a = max(read_count, a)
        elif base != "N":
            o += read_count
    return r, a, o


def consensus_label(r, a, o, min_fraction=0.6, min_reads=3):
    """
    Rule A, equation (1): the consensus label of a cluster's usable reads.

    A label X is assigned iff count_X >= n_min and count_X / (r + a + o) >= f_min.
    OTH is tested on o itself. With f_min > 0.5 at most one allele qualifies.

    Parameters:
        r, a, o (int): REF, ALT and OTH usable reads (see allele_counts).
        min_fraction (float, optional): f_min. Defaults to 0.6.
        min_reads (int, optional): n_min. Defaults to 3.

    Returns:
        str or None: 'WT', 'MUT' or 'OTH', or None if no allele qualifies.
    """
    usable = r + a + o
    for label, count in (("WT", r), ("MUT", a), ("OTH", o)):
        if count and count >= min_reads and count / usable >= min_fraction:
            return label
    return None


def classify_counts(r, a, o, alt, min_fraction=0.6, min_reads=3):
    """
    Class of a covering cluster: sub-n_min, a consensus label or discarded.

    Parameters:
        r, a, o (int): REF, ALT and OTH usable reads (see allele_counts).
        alt (str): The alternative allele, for the legacy indel shortcuts.
        min_fraction (float, optional): f_min. Defaults to 0.6.
        min_reads (int, optional): n_min. Defaults to 3.

    Returns:
        str: SUBNMIN if r + a + o < n_min, else 'WT', 'MUT' or 'OTH' by rule A,
             else DISCARDED.
    """
    if r + a + o < min_reads:
        return SUBNMIN
    # Legacy indel shortcuts, kept until tolerant indel ALT removes them
    if len(alt) == 2 and a > 2 and a / (r + a) > 0.25:
        return "MUT"
    if len(alt) > 3 and a > 3:
        return "MUT"
    return consensus_label(r, a, o, min_fraction, min_reads) or DISCARDED


def consensus_on_umi(line, ref, alt, min_fraction=0.6, min_reads=3):
    """
    Tally the clusters of one bam-readcount -p line, where each read group
    (a UMI cluster) has its own library ID in the BAM's RG.

    Parameters:
        line (str): One bam-readcount -p line.
        ref (str): The reference base.
        alt (str): The alternative allele, in bam-readcount notation.
        min_fraction (float, optional): f_min. Defaults to 0.6.
        min_reads (int, optional): n_min. Defaults to 3.

    Returns:
        Dict: See tally_clusters.
    """
    return tally_clusters(parse_clusters(line), ref, alt, min_fraction, min_reads)


def tally_clusters(clusters, ref, alt, min_fraction=0.6, min_reads=3):
    """
    Tally the classes of parsed clusters (see parse_clusters) at one target.

    Returns:
        Dict: Number of clusters labelled 'WT', 'MUT' and 'OTH' and 'discarded',
              plus 'subnmin_umis', the number of sub-n_min clusters, and
              'subnmin_reads', their usable reads.
    """
    counts = {"WT": 0, "MUT": 0, "OTH": 0, DISCARDED: 0, "subnmin_umis": 0, "subnmin_reads": 0}
    for entries in clusters:
        r, a, o = allele_counts(entries, ref, alt)
        cls = classify_counts(r, a, o, alt, min_fraction, min_reads)
        if cls == SUBNMIN:
            counts["subnmin_umis"] += 1
            counts["subnmin_reads"] += r + a + o
        else:
            counts[cls] += 1
    return counts


def count_targets(target, lines, min_fraction=0.6, min_reads=3):
    """
    Count consensus labels and sub-n_min clusters per target for one cell.

    Parameters:
        target (pd.DataFrame): Targets with columns ['chr', 'start', 'end', 'ref', 'alt', 'gene'].
        lines (Iterable[str]): bam-readcount -p output lines, SNV pass before indel pass.
        min_fraction (float, optional): f_min. Defaults to 0.6.
        min_reads (int, optional): n_min. Defaults to 3.

    Returns:
        pd.DataFrame: A copy of target with the COUNT_COLUMNS: 'WT', 'MUT' and
                      'OTH' labelled clusters, 'subnmin_umis' and 'subnmin_reads'.
    """
    target = target.copy()
    target[COUNT_COLUMNS] = 0

    # A site with both an SNV and an indel target has one line per bam-readcount
    # pass, the SNV pass first: SNV alts are scored on the first, indel alts on the last.
    # bam-readcount writes a covered site even when -b filters every base.
    scored_snvs = set()
    for line in lines:
        chrom, start_pos, ref, clusters = parse_readcount_line(line)
        alts = target["alt"][
            (target["chr"] == chrom)
            & (target["start"] == start_pos)
            & (target["ref"] == ref)
        ]
        for alt in alts:
            if not alt.startswith(("+", "-")):
                if (chrom, start_pos, alt) in scored_snvs:
                    continue
                scored_snvs.add((chrom, start_pos, alt))
            counts = tally_clusters(clusters, ref, alt, min_fraction, min_reads)
            # Rows are matched on chr, start and alt only, as on main
            for val in COUNT_COLUMNS:
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
