#!/bin/bash
# Why does bam-readcount report fewer deletion reads at the KIT site than the alignments hold?
# Whole-cell bam-readcount (no per-library split) at -b 20 vs -b 0, against a pysam count, for the
# 5 cells with the most site-line reads. Aggregates printed only; scratch in kit_debug/.
D=/hpcnfs/scratch/P_DIMA_SCMSEQ/analysis/nanopore_enrichment/zhan_analysis/panel2026/PALM14E_subset
P=/hpcnfs/scratch/P_DIMA_SCMSEQ/analysis/nanopore_enrichment/zhan_analysis/panel2026/PALM14E
REF=/hpcnfs/scratch/P_DIMA_SCMSEQ/resources/ref_standard/GRCh38_full_analysis_set_plus_decoy_hla.fa
IMG=/hpcnfs/data/DIMA/singularity_cache/acox1-scmseq-1.0.img
run() { singularity exec -B /hpcnfs $IMG "$@"; }
mkdir -p $D/kit_debug && cd $D/kit_debug
grep -P '^chr4\t54723603\t' $P/target.bed | cut -f1-3 > site.bed
for f in $D/cells/*/allcounts.count; do
  n=$(grep -P '^chr4\t54723603\t' $f | cut -f4)
  [ -n "$n" ] && echo "$n $(dirname $f)"
done | sort -k1,1nr | head -5 | cut -d' ' -f2 > top_cells.txt
i=0
for cd in $(cat top_cells.txt); do
  i=$((i+1))
  run samtools view -b -q 20 -o c$i.bam $cd/filtered_input.bam && run samtools index c$i.bam
  for b in 20 0; do
    echo "cell$i -b $b: $(run bam-readcount -f $REF -l site.bed -b $b -w 1 c$i.bam 2>/dev/null | tr '\t' '\n' | grep -E '^(T|-TTA|-TT|-T):' | cut -d: -f1,2 | tr '\n' ' ')"
  done
  run python - c$i.bam <<'EOF'
import sys, collections, pysam
c = collections.Counter()
for col in pysam.AlignmentFile(sys.argv[1]).pileup("chr4", 54723601, 54723603, truncate=True, stepper="nofilter",
                                                    min_base_quality=0, max_depth=10**8):
    for pr in col.pileups:
        if col.reference_pos == 54723601 and pr.indel < 0:
            c[("del_after_anchor", pr.indel)] += 1
        if col.reference_pos == 54723602 and not pr.is_del:
            q = pr.alignment.query_qualities[pr.query_position]
            c[("base_q20" if q >= 20 else "base_lowq", pr.alignment.query_sequence[pr.query_position])] += 1
print("  pysam:", {f"{a}{b}": n for (a, b), n in c.items() if n > 5})
EOF
done
rm -f c*.bam c*.bam.bai
