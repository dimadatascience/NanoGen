# main → genotype-v2 attribution ladder on PALM14E

Wayfinder ticket: *Run the main → v2 attribution ladder on PALM14E* (#21), map #3. Runs part (b) of *Regression check: main vs genotype-v2 calls on PALM14E* (#19): cumulative rungs S0 (`main` emulated) → S6 (v2 expected), plus the off-ladder caller point S0+6. Every hard flip must trace back to its rung's decided mechanism.

**Patient data.** Everything below is aggregate, with targets T1–T7 in BED order. T3 is the KIT 3-bp deletion; T4 and T5 are SNVs a few bp downstream and share its reads. Independent loci are therefore {T1}, {T2}, {T3, T4, T5}, {T6} and {T7}. Per-cell and per-flip tables stay on the HPC, in `PALM14E_subset/regression/private/` (group `deo_p_dima_scm`).

## Setup

- **Inputs:** no new BLAZE or consensus run.
  - `records.pkl` holds q20 site records with their pysam pileup allele for 4,926 cells (`followup.py`, #9).
  - The T3 site line comes from `cells/<CB>/filtered_input.bam`, via `kit_indel.read_info` (#11). It gives both `-b 20` and `-b 0` entries.
- **Run:** `regression/ladder.py`, submitted with `regression/run_ladder.sh` as a Slurm job (singularity image `acox1-scmseq-1.0.img`). The full run took 8 min and 8 GB.
  - A 600-cell smoke test ran first (`--max-cells 600`).
  - `fault_review.py` runs as a dependent job and prints the acceptance-3 cells with anonymised IDs.
- **Rungs** (cumulative):

| rung | records | grouping | n_min | consensus rule | T3 site line | caller |
|---|---|---|---|---|---|---|
| S0 | all q20 | BX × chr × 100 bp window, group ≥ 2 records | 2 | `main` (`consensus_on_umi`) | `-b 20` | `main` (`call_main_exact`) |
| S1 | primary (`-F 0x904`) | same | 2 | `main` | `-b 20` | `main` |
| S2 | primary | cell × target, LD ≤ 2 hybrid T = 3 | 2 | `main` | `-b 20` | `main` |
| S3 | primary | same | 3 | `main` | `-b 20` | `main` |
| S4 | primary | same | 3 | spec: uniform fmin 0.6, OTH tested on o, T3 tolerant k = 1, no shortcuts; gate r + a + o ≥ n_min | `-b 20` | `main` |
| S5 | primary | same | 3 | spec | `-b 0` | `main` |
| **S6** | primary | same | 3 | spec | `-b 0` | spec: per-target pooled P_s = max(w / dN, 1e-4), d = 2 SNV / 1 indel, binomial upper tail, 0.1 / 0.9 |
| S0+6 | as S0 | | | | | spec |
| S4n | as S4, but the gate counts all cluster reads (§4) | | | | | `main` |

- **"Hard flip"** means MUT ↔ WT; a **"soft flip"** moves to or from NA. A row is **covered** when it has at least one consensus UMI at either rung.

## 0. S0 is tied to `main`'s real output

| check | agreement |
|---|---|
| S0 counts vs `counts_rerun.csv` (`main` on this subset), covered rows | 98.7 % count-identical (3,877 rows); UMI totals WT 24,102 vs 24,097, MUT 23,528 vs 23,497 |
| S0 calls vs `palm14e_rerun.csv` | **34,470 / 34,482** |
| S0 calls vs the original `palm14e.csv` | 34,434 / 34,482 |
| caller port on `counts_rerun.csv` vs `palm14e_rerun.csv` | 34,482 / 34,482 |

This matches the reproduction from the UMI measurement (34,943 / 34,979 there, on the 4,997-cell grid).

## 1. Per-rung result

Hard flips summed over targets. The "against" column counts hard flips in a direction or at a target that the #19 prediction did not allow. The last column is the largest per-target share of covered rows.

| step | covered rows | MUT → WT | WT → MUT | against prediction | worst target | prediction (#19) | verdict |
|---|---|---|---|---|---|---|---|
| S0 → S1 | 3,881 | 10 | 8 | 18 | T7 0.96 % | MUT → NA at T2/T7; few hard | **explained** |
| S1 → S2 | 3,519 | 149 | 20 | 20 | T3 1.7 % (3 / 175) | counts ≈ 0.3×; MUT → NA; MUT → WT | **explained** |
| S2 → S3 | 3,402 | 25 | 0 | 25 | T4 1.3 % (2 / 157) | soft to NA; ≈ no hard | **explained** |
| S3 → S4 | 2,913 | 27 (all T3) | 9 (T6 7, T7 2) | 9 | T6 0.58 % | T3 MUT → WT only; SNVs OTH only | **explained** |
| S4 → S5 | 2,929 | 1 (T3) | 1 (T3) | 0 | – | T3 only, few, either way | **explained** |
| S5 → S6 | 2,922 | 0 | 0 | 0 | – | either way | **explained (no-op, §3)** |
| S0 → S0+6 | 3,877 | 0 | 0 | – | – | off-ladder | caller alone moves nothing |

**The formal 0.1 % tripwire fires at S1–S4, and none of these rungs blocks.** Two reasons:

- **The threshold is mis-scaled.** #19 set 0.1 % as "the rerun noise", i.e. 36 differing calls out of 34,979 rows. Those 34,979 rows are all cell × target rows, and about 89 % of them are uncovered NA rows. Over covered rows, the same rerun noise is **≈ 0.9 %** (36 / ~3,880), and S0's own emulation noise is 0.3 % (12 / 3,877). Against a covered-row threshold, the rungs sit at or near noise level.
- **Every against-prediction flip moves exactly one MUT UMI, through the rung's own mechanism:**
  - In each one, m goes 0 ↔ 1 (A3).
  - None has unchanged counts, so no flip comes from the caller or from the cross-target `err` of `main`'s caller. That `err` sits at the floor in every cell (§3).
  - At the 1e-4 floor, a single MUT UMI calls MUT unless the cell has ≳ 10³ WT UMIs (revised document, Table 2). So any rung that adds or removes one UMI moves a call. The cause is the calling rule's known single-UMI exposure, which is reported and out of scope, not a defect of the rung.

### Mechanism per rung

**(k, m, o): median [IQR] before → after.**

- **S1, primary only:**
  - T2 loses half its MUT calls to NA (481 → 243, 237 MUT → NA). T2 is supplementary-heavy, as predicted; T7 loses 47.
  - The 18 hard flips are single-UMI edges:
    - **MUT → WT:** the lone MUT UMI disappeared with its supplementary records (m 1 → 0, k unchanged).
    - **WT → MUT:** removing supplementary records from a mixed BX × window group left an ALT-majority group, e.g. k 4 → 1 with m 0 → 1 at T7.
- **S2, grouping:**
  - **Counts shrink 0.19× (WT) and 0.18× (MUT)**, more than the predicted ≈ 0.3×. This is consistent with #9's 5–8× over-counting.
  - **MUT → WT (149, predicted):** k collapses (T6 20 → 4, T7 16 → 3) and the lone MUT UMI goes (m 1 → 0). It was a fragment that the window split had separated from a WT molecule.
  - **WT → MUT (20, against):** k shrinks and a single MUT UMI appears (m 0 → 1). Bundling a molecule's reads across windows lets an ALT-majority cluster reach the gate. T3's 3 cases have k 2 → 0–1. This is the grouping mechanism itself, at single-UMI scale.
- **S3, n_min 3:**
  - All 25 hard flips are MUT → WT, each the loss of a lone 2-read MUT UMI (m 1 → 0, k ≥ 2 left). There are no WT → MUT flips.
  - The direction and mechanism are exact. "≈ none" underestimated how many MUT calls rest on one 2-read UMI: 61 % of MUT calls at S2 rest on a single MUT UMI.
- **S4, consensus rule:**
  - **T3: 27 MUT → WT** (k 3 → 4, m 1 → 0). These are the deep WT molecules that `main`'s `len > 3` shortcut had labelled MUT, as found in #11.
  - **SNVs: 9 WT → MUT** (k unchanged, m 0 → 1). This is the spec's gate, which #19's prediction left out. A cluster with 3 usable reads at 2 ALT : 1 REF passes the spec rule (r + a + o ≥ 3 and ALT fraction 0.67 ≥ 0.6). `main` drops it, because it requires max(r, a) ≥ n_min.
  - **OTH reclassification:** 1 OTH UMI at T7 and none elsewhere (§3).
- **S5, `-b 0` at T3:** 2 hard flips and 13 soft flips, all at T3. The SNVs are untouched.

## 2. Per cell, S0 → S6

- **Call changes:** 1,049 of 4,926 cells change at least one call. 836 of them change one target, 145 two, and 68 three or more.
- **Hard flips:** 171 cells have one. 188 hard flips in total: 176 MUT → WT and 12 WT → MUT.
  - Rung at which they first depart from S0: S1 8, S2 141, S3 11, S4 27, S5 1.
- **Acceptance 3:** **2 cells** have hard flips at ≥ 2 independent loci. Both were reviewed on the HPC (`fault_review.py`):
  - Both are deep WT cells. Their median S0 depth is 365 UMIs; the median cell has 7.
  - In both, `main` called MUT on a single MUT UMI against 10–247 WT UMIs (T5, T6, T7), or on T3's shortcut (18 shortcut MUT UMIs against 222 WT).
  - Each flip happens at the rung that removes that UMI: S1 for T7, S2 for T5/T6, S4 for T3. Every flip goes MUT → WT.
  - **No systematic v2 cause.** v2 makes these cells consistently WT.
- **MUT calls overall:** 2,293 → 1,476.
  - **T2 drives most of it:** 481 → 140. The bulk goes to NA, from losing supplementary records and from n_min.
  - **Single-UMI MUT calls:** 1,151 → 884. As a share of MUT calls, they rise from 50 % to 60 %.

## 3. Error table, and why S6 is a no-op on PALM14E

| target | cells covered | `main` per-cell err (median, IQR, max) | cells hitting zero → min fallback | cells at the 1e-4 floor | v2 ŵ / N at S6 | v2 P_s at S6 |
|---|---|---|---|---|---|---|
| T1–T7 | 349 / 484 / 196 / 176 / 173 / 1,343 / 1,146 | 0 (0–0), max 0 | all | all | 0 at every target | 1e-4 at every target |

- **`main`'s per-cell err is identically 0.** Its consensus rule never assigns OTH: there are 0 OTH UMIs across all S0 rows, which is the limitation the revised document corrects. Every cell falls back to the minimum, 0, and then to the floor.
- **Under v2's rule, OTH is still empty.** S6 has 0 OTH UMIs at every target. With fmin 0.6, a third allele almost never wins a consensus UMI, so w = 0 and P_s sits at the floor.
  - In the first run, where the gate counted all cluster reads, T3 and T7 had one OTH UMI each (P_s 2.3e-3 and 2.2e-4). That still flipped nothing.
- **The two callers are therefore the same function on PALM14E.** The nbinom cdf and the binomial tail are equal at P_s = 1e-4, so S5 → S6 and S0 → S0+6 have zero changes. The ladder cannot exercise the caller here. The exact caller check of the merge gate, (a2) in #19, is what covers it.

## 4. The n_min gate needs a decision

- **The ambiguity:** which reads count toward n_min is left open.
  - The revised document (section 1) counts r, a and o "among reads with base quality ≥ 20" and keeps a UMI "supported by at least n_min reads". The ladder's S4 reads this as **r + a + o ≥ n_min**.
  - The measurement helpers behind #12/#13 (`hybrid.consensus`) and #11 (`kit_indel.uniform_label` on the new grouping) gated on **all records in the cluster** instead. That includes reads with no usable base at the site (q < 20, a deletion over an SNV, N).
- **S4n measures the gap.** It is S4 with the cluster-read gate:

| | S3 → S4 (usable reads) | S3 → S4n (cluster reads) | S4 → S4n |
|---|---|---|---|
| SNV WT → MUT | 9 | 22 | 15 |
| NA → MUT (all targets) | 73 | 211 | 138 |
| MUT calls (all targets) | 1,475 | 1,628 | |

- **The cluster-read gate turns UMIs with a single usable ALT read into MUT UMIs.** Example: 1 ALT read plus 2 reads without a site base gives n = 3 and fraction 1.0.
- **It doesn't block the ladder:** both readings are explained. But the implementation must pick one, and the merge gate's emulation (a3) must use the same one. This is raised as a new ticket, *n_min gate: usable reads or all cluster reads*.
- **Resolved ([n_min gate: usable reads or all cluster reads](https://github.com/dimadatascience/NanoGen/issues/22)):** neither S4 nor S4n. v2 applies the revised document's equation (1) as written.
  - The winning allele must itself have ≥ n_min usable reads, and its fraction of r + a + o must be ≥ f_min. This is the same for SNV and indel targets.
  - S4 departed from this: it gated r + a + o ≥ n_min and then labelled 2 ALT : 1 REF UMIs. Under equation (1) those UMIs are discarded. The expected S3 → S4 change should therefore be the T3 shortcut removal only, with no SNV WT → MUT expected.
  - The merge gate's emulation (a3) must use the same rule.

## Verdict

**No rung blocks, and no decision is reopened.**

- **Rungs S1–S5** are explained by their decided mechanisms.
  - The formal 0.1 % tripwire fires only because the threshold was calibrated over all rows and applied over covered rows.
  - Every flip against the prediction is a one-UMI edge produced by the rung's own change.
- **S6 is a no-op on PALM14E**, because w = 0 everywhere. The caller is covered by the merge gate's (a2) exact check, not by this ladder.
- **Two predictions were too tight:**
  - S2's count shrink is 0.19×, not ≈ 0.3×.
  - S4 also changes SNV calls, through the gate, not only through OTH.
- **One new decision:** the n_min gate definition (§4).
