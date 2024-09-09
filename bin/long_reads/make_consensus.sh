#!/bin/bash

function usage {
    echo -e "usage : split_bam_by_tag_and_window_count_consensus.sh -i INPUT -b BEDFILE -r REFERENCE -q QUALITY -c COUNT -w WINDOW -m MINBASEQUAL -d UMIDIST [-h]"
    echo -e "Use option -h|--help for more information"
}

function help {
    usage
    echo 
    echo "This script groups reads based on UMIs using UMI tools."
    echo "It splits a BAM file based on BX:Z: tags and a specified window around the mapping positions."
    echo "Reads are first filtered to only include those overlapping with regions specified in a BED file."
    echo "Reads must meet a minimum mapping and base quality threshold."
    echo ""
    echo "---------------"
    echo "OPTIONS"
    echo
    echo "   -i|--input INPUT : BAM file as input"
    echo "   -b|--bed BEDFILE : BED file with genomic coordinates for filtering in the following format: ['chr', 'start', 'end', 'ref', 'alt', 'gene']"
    echo "   -r|--reference REFERENCE : Reference FASTA file"
    echo "   -q|--quality QUALITY : Minimum mapping quality of reads, default 20"
    echo "   -m|--minbasequal MINBASEQUAL : Minimum base quality, default 20"
    echo "   -c|--count COUNT : Minimum read count supporting UMI to generate a BAM file, default 3"
    echo "   -d|--umidist UMIDIST : Maximum distance to collapse UMIs using umitools group. Default 1."
    echo "   -w|--window WINDOW : Window size (in bp) for grouping reads around their mapping positions. Default 1000"
    echo "   [-h|--help] : Displays this help message"
    exit
}

# Transform long options to short ones
for arg in "$@"; do
  shift
  case "$arg" in
    "--input") set -- "$@" "-i" ;;
    "--bed") set -- "$@" "-b" ;;
    "--reference") set -- "$@" "-r" ;;
    "--quality") set -- "$@" "-q" ;;
    "--minbasequal") set -- "$@" "-m" ;;
    "--count") set -- "$@" "-c" ;;
    "--window") set -- "$@" "-w" ;;
    "--umidist") set -- "$@" "-d" ;;
    "--help") set -- "$@" "-h" ;;
    *) set -- "$@" "$arg"
  esac
done

input=""
bedfile=""
reference=""
quality=20
minbasequal=20
count=3
window=1000
umidist=1

while getopts ":i:d:b:r:q:m:c:w:h" opt; do
    case $opt in
        i) input=$OPTARG;;
        b) bedfile=$OPTARG;;
        d) umidist=$OPTARG;;
        r) reference=$OPTARG;;
        q) quality=$OPTARG;;
        m) minbasequal=$OPTARG;;
        c) count=$OPTARG;;
        w) window=$OPTARG;;
        h) help;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            usage
            exit 1
            ;;
        :)
            echo "Option -$OPTARG requires an argument." >&2
            usage
            exit 1
            ;;
    esac
done

# Check for mandatory options
if [ -z "$input" ] || [ -z "$bedfile" ] || [ -z "$reference" ]; then
    echo "Input BAM, BED file, and reference FASTA must be specified."
    usage
    exit 1
fi

if [ ! -f "$input" ] || [ ! -f "$bedfile" ] || [ ! -f "$reference" ]; then
    echo "Specified file does not exist: $input or $bedfile or $reference"
    exit 1
fi

bedfile=`realpath $bedfile`
input=`realpath $input`
reference=`realpath $reference`

cb=`basename $input | sed 's/\.bam//g'`

# Script Logic
mkdir -p temp_$cb
samtools view -H $input > bam_header.sam

# Filter reads based on BED file
if [ -f "$bedfile" ]; then
	bedtools intersect -abam $input -b $bedfile -u > filtered_input.bam
else
	cp $input filtered_input.bam
fi
samtools sort filtered_input.bam -o sorted.bam
mv sorted.bam filtered_input.bam
samtools index filtered_input.bam

# Group using umi_tools
umi_tools group -I filtered_input.bam --group-out=grouped_reads.tsv --output-bam -S output.bam --edit-distance-threshold $umidist --method adjacency --paired --umi-tag=UB  --extract-umi-method=tag
samtools view -o output.sam output.bam
cat bam_header.sam output.sam | samtools view -bS - > filtered_input.bam 
samtools index filtered_input.bam

# Directly write reads to SAM files by UB tag and window
samtools view -q $quality filtered_input.bam | awk -v win=$window -v outdir="temp_$cb/" '{ 
    for(i=12; i<=NF; ++i) {
        if($i ~ /^BX:Z:/) {
            fname = outdir substr($i, 6) "_" $3 "_" int($4/win)*win ".sam";
            print $0 >> fname;
        }
    }
}'

# Convert SAM files to BAM if they meet the read count requirement
for samfile in temp_$cb/*.sam; do
    read_count=$(wc -l < "$samfile")
    if [ "$read_count" -ge $count ]; then
        name=`basename $samfile | sed 's/\.sam//g'`
        cat bam_header.sam "$samfile" | samtools view -bS - > "${samfile%.sam}.bam"
        samtools addreplacerg -r "@RG\tID:group1\tSM:sample1\tLB:'$name'\tPL:illumina" -o "${samfile%.sam}.lb.bam" "${samfile%.sam}.bam"
    else
        rm "$samfile" # Remove SAM file if read count is too low
    fi
done

cd temp_$cb
	allbams=""
	for file in `ls *lb.bam`; do
		allbams+=" $file"
	done
	samtools merge merged.bam $allbams
	samtools index merged.bam
	bam-readcount -f ${reference} --site-list ${bedfile} -b $minbasequal merged.bam -p  > allcounts.count
cd ..


##
