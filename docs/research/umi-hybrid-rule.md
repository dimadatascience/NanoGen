# Hybrid UMI merge rule on PALM14E

Wayfinder ticket: *Measure the hybrid UMI merge rule on PALM14E* (#13), map #3. Follows *UMI grouping redesign: key, distance and n_min* (#12) and the measurement in [`umi-measurement.md`](umi-measurement.md) (#9).

**Patient data.** Everything below is aggregate. Targets are **T1–T7** in BED order. T3 (the 3-bp deletion), T4 and T5 lie within 6 bp of each other, so they share reads and move together. Per-read and per-cell tables stay on the HPC in `PALM14E_subset/` (group `deo_p_dima_scm`).

## Setup

- **Input:** the per-read table `records.pkl` cached by `followup.py`. It holds site-covering MAPQ ≥ 20 records of the on-target subset, with alleles from the bam-readcount-emulating pileup. Only primary records are kept, deduplicated per read × target: 6.76 M read × target rows in 4,904 cells.
- **Bundle:** cell × target, as decided in #12. Indel reads are those covering the site, not the ±k overlap; the difference is negligible here.
- **Rule** (`hybrid.py`, run by `stage4_hybrid.sh`):
  - UMI *a* absorbs *b* if LD(a, b) ≤ 2 and either n_a ≥ 2·n_b − 1, or n_b ≤ T **and n_a ≥ n_b**.
  - Clusters grow from the most abundant UMI; ties break by count, then sequence.
  - The n_a ≥ n_b guard keeps growth downhill. Without it, a 1-read child already in a cluster could pull in a 3-read UMI. **T3-literal** drops the guard to show its effect.
- **Scoring:** clusters are scored with the uniform rule (n ≥ n_min = 3, top allele ≥ fmin 0.6) and called with `main`'s `genotyping.py` logic (`call_main`, as in check 5). Only the grouping changes.

Schemes: **plain** (directional, LD ≤ 2, T = 0), **T2**, **T3**, **T5**, and **T3-literal**.

## 1. Extra merges and molecules per cell × target

Summed over T1–T7. Plain has 94,309 clusters, 8,979 of them with ≥ 3 reads.

| scheme | extra merges (clusters removed) | UBs absorbed via clause | reads via clause | clusters ≥ n_min | Δ vs plain |
|---|---|---|---|---|---|
| T2 | 1,098 (1.2 %) | 7,433 | 14,866 (0.22 %) | 8,835 | −144 (−1.6 %) |
| T3 | 1,536 (1.6 %) | 9,713 | 21,802 (0.32 %) | 8,397 | −582 (−6.5 %) |
| T5 | 1,961 (2.1 %) | 11,525 | 29,682 (0.44 %) | 7,972 | −1,007 (−11 %) |
| T3-literal | 3,344 (3.5 %) | 17,507 | 39,963 (0.59 %) | 7,444 | −1,535 (−17 %) |

- **Molecules mostly go down.** The clause mostly folds 3-read clusters into a sibling. Each such merge removes a molecule that passed n_min. Few merges lift 2 + 2 fragments over n_min.
- **Per cell × target under T3:** the count of molecules with ≥ 3 reads changes in 20 % of covered cell × target pairs at T6 (253 of 1,277), 11 % at T7 and 10 % at T1, almost always downward.
- **Means shift only slightly:** T6 goes from 2.96 to 2.75 and T7 from 2.32 to 2.19. The T3–T5 locus goes from 2.9 to 2.7.

## 2. Purity of reads absorbed through the clause

Minority-allele share among REF/ALT reads. Reads absorbed through the clause are compared with the intrinsic within-molecule rate (exact UB, n ≥ 2; check 4 / F1) and with reads absorbed through the ratio (plain).

| target | T1 | T2 | T3 | T4 | T5 | T6 | T7 |
|---|---|---|---|---|---|---|---|
| intrinsic (%) | 4.9 | 0.06 | 8.8 | 7.7 | 6.5 | 12.7 | 11.1 |
| plain, via ratio (%) | 5.4 | 0.08 | 9.9 | 9.0 | 7.7 | 12.4 | 11.4 |
| T2, via clause (%) | 5.5 | 0.06 | 11.8 | 11.5 | 8.3 | 11.3 | 10.6 |
| T3, via clause (%) | 5.8 | 0.04 | 10.9 | 11.9 | 8.8 | 11.2 | 10.9 |
| T5, via clause (%) | 5.5 | 0.03 | 10.1 | 10.6 | 8.0 | 11.2 | 10.6 |
| T3-literal, via clause (%) | 5.9 | 0.04 | 12.5 | 14.0 | 11.8 | 13.3 | 11.7 |

- **Where clause reads match the baseline:** at T1, T6 and T7 they sit at or below the intrinsic rate. At T2, the near-clonal target, the measure can't detect false merges, because nearly all molecules share one allele.
- **At the T3–T5 locus** clause reads run 1–4 points above intrinsic and about 1–3 points above ratio-absorbed reads. That is about 1,000 UBs from one shared set of reads, and the excess doesn't grow with T.
- **Leave-own-UB-out majorities** give the same picture (within 0.3 points).
- **T3-literal is worse** everywhere except T2 (up to 14.0 % at T4). The downhill guard matters.

## 3. Below n_min (= 3)

| scheme | clusters below n_min | reads below n_min |
|---|---|---|
| plain | 90.5 % | 1.33 % |
| T2 | 90.5 % | 1.31 % |
| T3 | 90.9 % | 1.31 % |
| T5 | 91.4 % | 1.31 % |

The clause barely touches sub-n_min mass. Stray single reads stay stray; most are not within LD 2 of a small sibling.

## 4. Calls (`main`'s caller, n_min 3) against plain

The comparison covers 34,328 cell × target rows.

| scheme | rows changed | MUT→WT/NA | WT→NA/MUT | NA→MUT/WT |
|---|---|---|---|---|
| T2 | 16 (0.05 %) | 3 | 9 | 4 |
| T3 | 40 (0.12 %) | 9 | 25 | 6 |
| T5 | 62 (0.18 %) | 16 | 40 | 6 |
| T3-literal | 121 (0.35 %) | 27 | 88 | 6 |

- **Small effect:** MUT calls total 1,631 under plain and 1,627 under T3. Most changes are WT→NA at T6 and T7, where merging a few WT molecules drops a cell below `main`'s WT evidence. T2 is unchanged under every scheme.
- **No spurious MUT calls appear:** the worry was that merging 2 + 2 MUT fragments would create MUT calls. NA→MUT is 1–3 rows.

## Conclusion

- **T = 3 is confirmed** as the `umi_small_reads` default.
  - The clause is a small correction: 1.6 % fewer clusters, 6.5 % fewer molecules passing n_min, and 0.12 % of calls changed.
  - Its merges are close to the purity of ratio merges.
  - T = 5 doubles the effect on molecules for no purity gain. T = 2 does half as much.
- **The implementation must state the downhill guard explicitly:** the small clause applies only when n_a ≥ n_b. Without it (T3-literal), merges triple, purity drops at the het targets (up to 14 % at T4), and three times as many calls change.

## Reproduce (HPC, group `deo_p_dima_scm`)

```bash
cd /hpcnfs/scratch/P_DIMA_SCMSEQ/analysis/nanopore_enrichment/zhan_analysis/panel2026/PALM14E_subset
sbatch scripts/stage4_hybrid.sh       # records.pkl -> hybrid_out/H1..H4 (about 6 min)
```
