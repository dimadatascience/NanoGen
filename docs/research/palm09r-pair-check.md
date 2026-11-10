# PALM09R pair check: T3 negatives at relapse

Wayfinder ticket: *PALM09R pair check: T3 negatives at relapse* (#81), map *Confirm genotype-v3 vs v3-ub on new diagnosis samples* (#72). Follows [`t3-cause-palm09e.md`](t3-cause-palm09e.md) (#73).

## Answer: not possible, because the relapse library doesn't cover T3

The check needed T3 reads in PALM09R's T/B/NK negatives. The relapse library has none to speak of, so there is nothing to genotype and nothing to pre-register or score. **T3's cause stays inconclusive**, as #73 left it, leaning towards a site-specific partial-molecule artefact.

- **The panel was redesigned at relapse.** The PALM09_R panel BED has 10 targets. It shares only T1 and T2 with PALM09_E. T3 and T4 are absent, and the 8 other targets are new.
- **The library confirms this.** The probe counts exact 20-mers in the raw FASTQ, at four flank positions on each side of each site, on both strands. They are reference sequence, which avoids the variant and the indel window.

| target | on E panel | on R panel | E flank hits (51.1 M reads) | R flank hits (56.2 M reads) |
|---|---|---|---|---|
| T1 | yes | yes | 5,698,193 | 4,283,979 |
| T2 | yes | yes | 18,786,626 | 8,291,416 |
| T3 | yes | **no** | 20,602,048 | **19** |
| T4 | yes | no | 2,995,039 | 230 |

The shared targets keep the same order of magnitude at relapse. T3 drops by six orders of magnitude, to off-target noise. A FASTQ-start genotype-v3 run of PALM09R therefore can't produce T3 data, so it was not run.

## Consequences for the map

- The *Outlier-flag revision* fog no longer waits on relapse. It rests on #73's leaning (artefact) alone, or on another sample that carries T3's gene. None of the census samples (PALM07E, 08E, 11E, 14E, 21E) or PALM05E do.
- PALM09R as a *FASTQ restart for more power* sample shares only T1 and T2 with PALM09E. Its other 8 targets would be new targets for a relapse sample, which the map has ruled out of the verdict anyway.

## Method

- Jobs 9023071 (PALM09R, 19 min) and 9023075 (PALM09E, 18 min) on `hpc_new` `short`, using the pinned NanoGen image (pysam, for reference extraction). Script: [`palm09r-pair-check/kmer_probe.sbatch`](palm09r-pair-check/kmer_probe.sbatch). The E job is the same script with E input and output names.
- Inputs: the PALM09E and PALM09R pass FASTQs listed in the local sample-paths sheet, and the PALM09_E / PALM09_R BEDs. Aggregate counts are in [`palm09r-pair-check/kmer_counts.tsv`](palm09r-pair-check/kmer_counts.tsv). K-mer sequences are not published because they would locate the targets.
- Some flank positions read zero in both libraries, for example T3 R2 and the T1 left flanks. These are at amplicon or read ends, and the same pattern appears in E and R.
