#!/bin/bash
#SBATCH -J q81_kmer
#SBATCH -p short
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=03:00:00
#SBATCH -o logs/kmer_%j.log
# Ticket 81: does the PALM09R relapse FASTQ cover T3 at all? Exact flank 20-mer counts, aggregate only.
set -o pipefail
Q=${Q:?working folder}
REF=${REF:?GRCh38 reference fasta}
NG=${NG:?pinned NanoGen singularity image}
O=${O:?analysis base folder (see local sample-paths sheet)}
FQ=$(cat $Q/fq_r.txt)
# label<TAB>chrom<TAB>1-based pos ; taken from the PALM09_E BED (T1..Tn in BED order)
awk -F'\t' '{print "T"NR"\t"$1"\t"$2}' $O/nanogen_results/bed_files/PALM09_E.bed > $Q/sites.tsv
singularity exec -B /hpcnfs -B /beegfs $NG python3 - $REF $Q/sites.tsv $Q/kmers.tsv <<'PY'
import sys, pysam
fa = pysam.FastaFile(sys.argv[1]); comp = str.maketrans("ACGT", "TGCA")
with open(sys.argv[3], "w") as out:
    for l in open(sys.argv[2]):
        t, c, p = l.split(); p = int(p)
        # 0-based pos0 = p-1; flanks avoid the site and the indel window: 20-mers at -45..-26, -25..-6, +6..+25, +26..+45
        for name, s in [("L2", p-1-45), ("L1", p-1-25), ("R1", p-1+6), ("R2", p-1+26)]:
            k = fa.fetch(c, s, s+20).upper()
            out.write(f"{t}\t{name}\t+\t{k}\n{t}\t{name}\t-\t{k.translate(comp)[::-1]}\n")
PY
cut -f4 $Q/kmers.tsv > $Q/kmers.txt
echo "== start $(date)"; ls -l "$FQ" | awk '{print $5}'
LC_ALL=C grep -o -F -f $Q/kmers.txt "$FQ" | sort | uniq -c > $Q/kmer_counts.txt
echo "== reads $(( $(wc -l < "$FQ") / 4 ))"
awk 'NR==FNR{c[$2]=$1;next}{print $1"\t"$2"\t"$3"\t"($4 in c?c[$4]:0)}' $Q/kmer_counts.txt $Q/kmers.tsv
echo "== done $(date)"
