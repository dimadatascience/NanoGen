# KIT indel labelling on PALM14E

Wayfinder ticket: *Measure indel labelling change on a KIT cell subset* (#11), map #3. Measures the rule decided in *Indel alleles under the uniform fmin consensus rule* (#5) against `main`@`c442a32`.

**Patient data.** Everything below is aggregate. The KIT deletion is **T3** (BED order). It is a 3-bp deletion; T4 and T5 are SNVs 3 and 5 bp downstream. Per-read and per-cell tables stay on the HPC in `PALM14E_subset/` (group `deo_p_dima_scm`), outputs in `kit_out/`.

## Setup

- **No new consensus run was needed.** The ticket asked for a re-run that keeps each cell's `allcounts.count`. `stage1_consensus.sh` (#9) had already done that for all 4,947 cells with `main`'s baseline params: window 100, `min_read` 2, `min_fraction` 0.6. Of those cells, 352 have a T3 site line, in 6,116 UMI × window libraries.
- **`main`'s grouping:** each bam-readcount library is scored from its site-line entries.
- **New grouping:** the decided grouping (#12/#13): primary records bundled per cell × target, LD ≤ 2 hybrid T = 3, n_min = 3 reads, fmin 0.6. Clusters come from `records.pkl` (4,904 cells). Their site-line entries come from a read-level emulation of bam-readcount.
- **Rules** (`kit_indel.py`, run by `stage5_kit_indel.sh`):
  - **main**: `consensus_on_umi`, i.e. the MAX over substring matches plus the `len(alt) > 3` shortcut (MUT if more than 3 ALT reads, whatever the fraction).
  - **new_k**: the #5 rule. ALT is the SUM of `-` entries within ±k bp that agree on the shared prefix. REF is the ref base minus `+` entries; OTH is the rest except N. No shortcuts; the label is the top allele with fraction ≥ fmin.
  - **window_k1** (realignment-lite): ALT if the read loses 3 ± 1 bp net within ±3 bp of the deleted bases; REF if there is no indel there and the site base is REF.
- **Callers:**
  - **main**: an exact port of `genotyping.py`@`c442a32`. Per-cell err from OTH over all targets, floor 1e-4, then the nbinom tail.
  - **spec**: the revised document + #6. P_s = max(w/N, 1e-4) pooled over cells, with d = 1 for indels. Binomial upper tail, cut at 0.1 / 0.9.
  - The two callers agree everywhere except where OTH UMIs appear (window rule and no-quality-filter rule).
- **Validation (K0), all 100 %:**
  - Re-applying `main`'s rule to `allcounts.count` reproduces `counts_rerun.csv` in every cell.
  - The read emulation reproduces bam-readcount's entries in every library.
  - The caller port reproduces `main`'s T3 calls.

## 1. UMI labels and cell calls

T3 consensus UMIs (labelled WT / MUT). The other ~2.8 k libraries on `main`'s grouping and ~7.9 k clusters on the new grouping are discarded.

| grouping | rule | WT UMIs | MUT UMIs | MUT via shortcut | MUT calls | WT calls | NA (with UMIs) | MUT calls from 1 MUT UMI |
|---|---|---|---|---|---|---|---|---|
| main | main | 1,688 | 1,576 | 728 | **126** | 11 | 63 | 40 |
| main | new_k0 | 1,824 | 1,427 | – | 101 | 35 | 63 | 48 |
| main | new_k1 | 1,824 | 1,447 | – | **101** | **35** | 63 | 48 |
| main | new_k2 | 1,824 | 1,449 | – | 101 | 35 | 63 | 48 |
| new | main | 182 | 258 | 186 | 124 | 0 | 46 | 82 |
| new | new_k1 | 232 | 204 | – | **84** | **30** | 54 | 50 |

- **`main`'s grouping, main → new_k1:**
  - UMIs: 137 MUT → WT, 7 MUT → discarded, 15 discarded → MUT (lifted by summing `-TT` with `-TTA`), 1 WT → discarded.
  - Cells: **24 MUT → WT and 1 MUT → NA** (126 → 101 MUT calls). Both callers agree.
- **New grouping, main → new_k1:**
  - UMIs: 50 MUT → WT, 4 MUT → discarded.
  - Cells: **30 MUT → WT and 10 MUT → NA**. Under the new grouping, `main`'s rule gives no WT calls at all.
- **k = 0, 1, 2 give identical calls in both groupings.** k only moves UMIs: 1,427 / 1,447 / 1,449 MUT UMIs on `main`'s grouping, and 203 / 204 / 204 on the new one.

## 2. The flipped UMIs are deep WT molecules (K6)

| grouping | UMIs | n | median reads | median ALT fraction | ALT < 10 % | 10–20 % | 20–40 % |
|---|---|---|---|---|---|---|---|
| main | main MUT → new WT | 137 | 250 | 0.10 | 70 | 60 | 7 |
| main | MUT under both | 1,432 | 3 | 1.00 | 0 | 0 | 0 |
| new | main MUT → new WT | 50 | 1,198 | 0.10 | 24 | 22 | 4 |
| new | MUT under both | 204 | 7 | 0.91 | 0 | 0 | 0 |

- **Most site reads sit in a few huge UMI groups.** 93.5 % of site reads are in 249 libraries with ≥ 50 reads (53 have ≥ 1,000), spread over 83 cells.
- **These are real UMIs.** Their UMIs are ordinary 12-mers: none is low-complexity (0.8 % of UMIs in groups under 50 reads are).
- **Their ~10 % ALT is ordinary within-molecule noise.** It matches the within-molecule minority rate already measured at T3 (8.8 %, `umi-hybrid-rule.md`).
- **So the shortcut made false MUT calls.** It labelled deeply sequenced WT molecules MUT because more than 3 of their reads carried the deletion. Dropping it corrects false MUT calls; it does not lose real ones. The genuine MUT UMIs are almost pure ALT and survive unchanged.
- **No recall is lost.** Under #5's standing rule there is nothing to widen.

## 3. Near-miss spectrum (K3)

Indel entries on the T3 site line.

| entry | bam-readcount (`-b 20`) | primary reads, no quality filter | ALT at k = 0 / 1 / 2 | `main` ALT |
|---|---|---|---|---|
| `T` (REF) | 171,641 | 243,558 | – | – |
| `-TTA` | 111,236 | 171,787 | ✓ ✓ ✓ | ✓ |
| `-TT` | 974 | 4,412 | ✗ ✓ ✓ | ✓ |
| `-T` | 82 | 2,570 | ✗ ✗ ✓ | ✗ |
| `-TTACG` | 6 | 156 | ✗ ✗ ✓ | ✓ |
| `-TTACGA`, `-TTACGACA`, longer | ≤ 5 each | ≤ 78 each | ✗ | ✓ |
| `+…` (any) | ≤ 8 each | ≤ 739 each | ✗ | ✗ |

- **k = 1 is adequate.**
  - The only near-miss of weight is `-TT`: 0.9 % of `-TTA` on the site line, 2.6 % before quality filtering.
  - k = 2 would add `-T`, a single-base deletion inside a TT run, which is the commonest ONT deletion error. It gains no calls.
  - The reference context (`…GAC|TTA|CGAC…`) leaves no equivalent left or right placement of `-TTA`.

## 4. Blind spot of the per-site line (K4)

Primary MAPQ ≥ 20 reads spanning ±3 bp around the deleted bases: 547,406 reads. 191,570 of them are deletion-like (net loss of 2–4 bp, with a deletion).

| where the aligner put it | share of deletion-like reads | site line records it as |
|---|---|---|
| deletion starting at the site (`-TTA`, `-TT`, `-TTAC`) | 91.8 % | tolerant ALT |
| shifted left (−2: 4,437; −1: 793; −4/−5: 841), or split, covering the site | 5.6 % | absent |
| shifted right (+1: 1,637; +2: 574) or split, site base intact | 2.0 % | REF |
| other | 0.5 % | OTH |

- **Labels and calls barely move.** Switching from the site rule to **window_k1** (reads judged by net deletion in the window):

  | grouping | UMI changes | cell changes | MUT ↔ WT flips |
  |---|---|---|---|
  | `main`'s | 111 (of ~3.3 k labelled), mostly discarded ↔ labelled | 7 MUT → NA, 1 WT → MUT (spec caller: 10 MUT → NA, since 13 OTH UMIs raise P_s to 0.004) | 0 |
  | new | 13 (of 436) | 5 of 168 (2 MUT → NA, 1 NA → MUT, 2 WT → NA) | 0 |

- **Answer: no escalation to haplotype realignment for KIT.**
  - Deletion-like reads off the line are 8 % at the read level. Among reads passing the anchor-quality filter the pipeline applies, it is under 4 %.
  - They never flip a cell between MUT and WT.

## 5. New finding: bam-readcount's `-b 20` also filters deletions

- **How the filter works.** bam-readcount applies `-b` (min base quality 20) to a deletion through the quality of its **anchor base**, the read base just before the deletion.
  - Checked by running bam-readcount with `-b 20` and `-b 0` on the five deepest cells (`brc_debug.sh`). At `-b 0` it matches pysam exactly.
  - At `-b 20` the anchor-quality rule reproduces the counts exactly (`brc_delq.py`). In two cells, `-TTA` gives 1,909 of 2,408 and 21,037 of 34,309.
- **It is asymmetric.** At T3 it removes 37 % of deletion-like reads but 28 % of REF reads, and near-misses far more: `-TT` 78 %, `-T` 97 %. In practice the k tolerance acts on an already filtered set.
- **Dropping it for deletions changes some calls:**
  - `main`'s grouping: about 1.3 k discarded libraries gain labels, because their top-allele count reaches 2. Cells: 20 NA → MUT, 6 NA → WT, 3 WT → MUT, 2 MUT → WT.
  - New grouping: 13 of 168 cells change, in both directions (7 MUT → NA, 2 MUT → WT, 2 NA → MUT, 1 WT → MUT, 1 WT → NA).
- **It persists in genotype-v2.** #12 keeps bam-readcount, so this filter carries over. Whether indel alleles should be quality-filtered this way is an open decision. It is not a k or realignment question.

## 6. Caveat for earlier measurements

- **The earlier caller helper differs from `genotyping.py`.** The `call_main` helper in `measure_umi.py`, used for the call numbers in #9 and #13, replaces zero per-cell error with the **minimum positive** error. `genotyping.py` uses `.min()` including zeros, so zero-error cells end up at the floor of 1e-4.
- **The difference matters when OTH UMIs are sparse.** The substitute error can exceed 0.1 and make single-UMI MUT cells NA. On the new grouping it flipped T3 MUT calls between 84 and 36 across rules that differ by 3 UMIs.
- **Where the comparisons stand.** Relative comparisons within #9/#13 used the helper consistently, but their absolute call counts may be off. `kit_indel.py` uses an exact port (`call_main_exact`, validated). The regression check should use it too.

## Files

- `kit-indel/kit_indel.py`: all tables K0–K6.
- `kit-indel/stage5_kit_indel.sh`: the sbatch driver.
- `kit-indel/brc_debug.sh`, `kit-indel/brc_delq.py`: the bam-readcount deletion quality check.
- `kit-indel/giant_groups.py`: UMI group sizes at the site.
