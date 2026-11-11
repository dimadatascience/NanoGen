# Outcome-blind power check for the new diagnosis samples

Ticket: [Outcome-blind power check for the new diagnosis samples](https://github.com/dimadatascience/NanoGen/issues/76), map [Wayfinder: confirm genotype-v3 vs v3-ub on new diagnosis samples](https://github.com/dimadatascience/NanoGen/issues/72).
Aggregate numbers only; targets as T1…Tn (key in the HPC README). HPC folder `q76_power/`, jobs 9023095 (smoke), 9023102 (calibration), 9023167–9023171 (census), 9023175 (rescale + collect).

## Answer

**No target in the census is evaluable** under the carried-over rule: ≥ 100 covered negatives in B and ≥ 200 WT+MUT UMIs in A, in ≥ 10 of 20 seeds. **One target is borderline: PALM07E T1.** Its final status turns on the real v3 n.

Even an evaluable target would be **underpowered** for calibration. Detecting a 1.5× miscalibration (p_cal < 0.05, 80 % power) needs about **E\* = 31** expected null calls in B. Here the best single target has E_B ≈ 2.6 at c = 1 % and stays below E\* even at c = 4 %. Only a test **pooled across all targets and all five patients** reaches roughly the power of one earlier sample (see [Pooled power](#pooled-power)).

The cause is biological, not depth. The negatives (T/B/NK cells) seldom express the targeted transcripts. PALM21E has 4,431 negatives but at most 326 covered on any target. PALM14E T3–T5 reach only 2–3 covered negatives.

By the ticket's own conditions, therefore:
- **The census is underpowered overall.** The FASTQ-restart patch graduates from the fog.
- **The stop rule is triggered as things stand.** No sample has an evaluable target, clonal or not. PALM07E T1 is the only exception still open: borderline, with its clonality still to be fixed. Whether to stop, or to rescue power by pooling or by a FASTQ restart, is the user's call. It is not settled here.

## Method (outcome-blind)

No allele is ever read. The full v3 runs come only after pre-registration, so coverage is taken straight from the canonical tagged BAMs ([`power_check.py`](power-check/power_check.py) `coverage`):

- **UMI clusters:** reads are grouped into UMI clusters by genotype-v3's own `umi_grouping.py`, imported at `61d8809`. That code is allele-blind: primary reads, MAPQ ≥ 20, bundle per target, Levenshtein hybrid merge.
- **Usable reads:** at an SNV, a read is usable if it has an aligned base with Q ≥ 20 at the site, or a deletion over it. At an indel, every bundle read is usable. Only the base quality is checked, never the base itself.
- **n_proxy:** clusters with ≥ 3 usable reads (`min_reads`). This counts what v3's consensus labels WT, MUT or OTH, plus clusters it discards for no allele reaching 60 %. So it is an **upper bound on v3's n = WT + MUT**.
- **Evaluability:** as in the earlier [pre-registration](https://github.com/dimadatascience/NanoGen/blob/76c3f99/docs/research/negative-control-preregistration.md). Negatives are restricted to QC cells; the split sorts the barcodes and uses `default_rng(seed).permutation`, first half = A. The 20 seeds are the earlier map's, a stand-in until the new pre-registration lists its own. The split is random halves, so the choice of seeds hardly moves these counts.
- **Calibration expectation E_B:** Σ P(MUT | n_i) over covered B negatives, through the imported `call` at θ = P_s + c, m_min = 1. This is the nominal null of criterion 1. Contamination is unknown before scoring, so c runs over a **fixed scenario grid**, 0.25 / 0.5 / 1 / 2 / 4 %. That spans the c_s seen on the last map (0.09 %–4 %), with P_s at the 1e-4 floor.
- **Power:** O ~ Poisson(r·E_B) stands in for the Poisson-binomial tail.
  - r80 is the smallest O/E that p_cal < 0.05 detects with 80 % power.
  - fail_r1 and fail_r15 are the rates at which criterion 1 fails (p_cal < 0.05 or O/E > 1.5) at true O/E 1 and 1.5.
  - c_min is the smallest c at which E_B ≥ E\*.

### Proxy vs v3 (calibration on the already-scored PALM22E / PALM09E)

| | result |
|---|---|
| cells with v3 n > proxy n, any target | **0**: the proxy is a strict upper bound |
| v3 / proxy UMIs over negatives, targets with ≥ 100 v3-covered negatives (8) | pooled 0.907, per target 0.685–0.960 |
| v3 / proxy covered negatives, same targets | pooled 0.951, per target 0.709–0.990 |
| evaluable, proxy vs v3 | 10 / 11 agree; PALM22E T3 is proxy 20 seeds vs v3 6 seeds |

The worst ratio is PALM09E T3, the partial-molecule artefact site (see [PALM09E T3 cause](https://github.com/dimadatascience/NanoGen/issues/73)). There, v3 discards mixed clusters. A long insertion such as PALM22E T2 is far off: v3 calls it OTH. Full table: [`out/calib.md`](power-check/out/calib.md).

That gives a three-way **verdict rule**:
- **not evaluable** if the raw proxy fails. This is final, because v3 ≤ proxy.
- **evaluable** if it still passes at the worst measured scale (UMIs × 0.685, cells × 0.709).
- **borderline** otherwise.

Power is computed at the pooled central scale (× 0.907, × 0.951).

## Per-target table

Raw = proxy; cen = central-scale estimate of v3. E_B, r80 and the fail rates are at the c = 1 % scenario. The full grid is in [`out/census.tsv`](power-check/out/census.tsv).

| sample (negatives in QC / Leukemic) | target | covered neg. | B cells (raw) | A UMIs (raw) | seeds raw / worst | covered Leukemic (cen) | E_B c=1 % | r80 | fail at O/E 1 / 1.5 | verdict |
|---|---|---|---|---|---|---|---|---|---|---|
| PALM07E (3,513 / 2,334) | T1 | 270 | 134 | 318 | 20 / 7 | 302 | 2.6 | 3.1 | 0.26 / 0.54 | **borderline** |
| | T2 | 12 | 6 | 6 | 0 / 0 | 20 | 0.06 | – | | not evaluable |
| | T3 | 89 | 44 | 88 | 0 / 0 | 230 | 0.66 | 6.5 | | not evaluable |
| | T4 | 88 | 44 | 46 | 0 / 0 | 188 | 0.41 | 10.5 | | not evaluable |
| | T5 | 174 | 86 | 181 | 0 / 0 | 433 | 1.6 | 4.2 | | not evaluable |
| PALM08E (877 / 6,578) | T1 | 2 | 1 | 1 | 0 / 0 | 15 | 0.01 | – | | not evaluable |
| | T2 | 0 | 0 | 0 | 0 / 0 | 29 | 0 | – | | not evaluable |
| | T3 | 86 | 43 | 57 | 0 / 0 | 2,330 | 0.48 | 8.9 | | not evaluable |
| | T4 | 96 | 49 | 180 | 0 / 0 | 1,840 | 1.4 | 4.8 | | not evaluable |
| | T5 | 36 | 19 | 47 | 0 / 0 | 813 | 0.31 | 9.7 | | not evaluable |
| | T6 | 1 | 1 | 0 | 0 / 0 | 84 | 0.01 | – | | not evaluable |
| | T7 | 140 | 71 | 210 | 0 / 0 | 2,360 | 1.7 | 4.1 | | not evaluable |
| | T8–T11 | 0 | 0 | 0 | 0 / 0 | 0 | 0 | – | | not evaluable (uncovered) |
| PALM11E (293 / 4,664) | T1 | 156 | 81 | 153 | 0 / 0 | 3,150 | 1.4 | 4.8 | | not evaluable |
| | T2 | 0 | 0 | 0 | 0 / 0 | 0 | 0 | – | | not evaluable (uncovered) |
| | T3 | 13 | 7 | 6 | 0 / 0 | 1,010 | 0.06 | – | | not evaluable |
| | T4 | 20 | 10 | 27 | 0 / 0 | 204 | 0.27 | 11.1 | | not evaluable |
| | T5 | 93 | 47 | 104 | 0 / 0 | 1,120 | 0.68 | 6.3 | | not evaluable |
| PALM14E (1,064 / 3,698) | T1 | 10 | 5 | 10 | 0 / 0 | 88 | 0.07 | – | | not evaluable |
| | T2 | 5 | 3 | 3 | 0 / 0 | 72 | 0.04 | – | | not evaluable |
| | T3 | 3 | 2 | 2 | 0 / 0 | 103 | 0.03 | – | | not evaluable |
| | T4 | 2 | 1 | 1 | 0 / 0 | 81 | 0.02 | – | | not evaluable |
| | T5 | 2 | 1 | 1 | 0 / 0 | 83 | 0.02 | – | | not evaluable |
| | T6 | 87 | 42 | 110 | 0 / 0 | 725 | 0.92 | 6.0 | | not evaluable |
| | T7 | 85 | 43 | 83 | 0 / 0 | 647 | 0.75 | 5.7 | | not evaluable |
| PALM21E (4,431 / 4,285) | T1 | 89 | 45 | 98 | 0 / 0 | 243 | 0.82 | 6.7 | | not evaluable |
| | T2 | 78 | 40 | 110 | 0 / 0 | 192 | 0.68 | 6.3 | | not evaluable |
| | T3 | 3 | 2 | 1 | 0 / 0 | 7 | 0.02 | – | | not evaluable |
| | T4 | 19 | 9 | 10 | 0 / 0 | 37 | 0.08 | – | | not evaluable |
| | T5 | 19 | 9 | 19 | 0 / 0 | 90 | 0.11 | – | | not evaluable |
| | T6 | 0 | 0 | 0 | 0 / 0 | 0 | 0 | – | | not evaluable (uncovered) |
| | T7 | 182 | 91 | 208 | 0 / 0 | 254 | 1.8 | 3.8 | | not evaluable |
| | T8 | 267 | 134 | 190 | 8 / 0 | 722 | 1.6 | 4.2 | | not evaluable (A UMIs just short) |
| | T9 | 106 | 54 | 54 | 0 / 0 | 297 | 0.51 | 8.4 | | not evaluable |
| | T10 | 152 | 76 | 82 | 0 / 0 | 500 | 0.76 | 5.6 | | not evaluable |
| | T11 | 326 | 170 | 184 | 4 / 0 | 733 | 1.8 | 3.8 | | not evaluable (A UMIs just short) |
| | T12 | 100 | 52 | 136 | 0 / 0 | 184 | 0.85 | 6.5 | | not evaluable |

"–" means no O/E below 20 is detectable. **c_min (1.5×) is above the 4 % top of the grid for every target.**

**Few negatives:**
- **PALM11E** (293 negatives; at best 81 covered in B) cannot reach 100 covered B cells on any target, whatever its depth.
- **PALM08E** (877) reaches at most 71 covered in B.

**Tumour side:** covered Leukemic cells are plentiful (often > 1,000). Tumour-side criteria such as sensitivity would be measurable where a target is clonal. Calibration is the binding constraint.

## Pooled power

Here E_B is summed over every target, at central scale (B cells = covered negatives in B, summed over targets). This is a reading for the pass bar ([Pre-register the confirmation pass bar](https://github.com/dimadatascience/NanoGen/issues/78)), not a decision.

| pool | B cells | E_B c=1 % | E_B c=2 % | E_B c=4 % | r80 at 1 % / 2 % / 4 % |
|---|---|---|---|---|---|
| PALM07E | 298 | 5.3 | 8.9 | 11.2 | 2.37 / 2.05 / 1.92 |
| PALM08E | 174 | 3.9 | 5.6 | 6.4 | 2.65 / 2.45 / 2.33 |
| PALM11E | 137 | 2.4 | 4.0 | 5.4 | 3.25 / 2.83 / 2.33 |
| PALM14E | 92 | 1.9 | 2.8 | 3.8 | 3.65 / 3.31 / 2.70 |
| PALM21E | 646 | 9.0 | 15.8 | 25.7 | 2.02 / 1.77 / 1.56 |
| **all five** | 1,347 | 22.4 | 37.1 | 52.4 | **1.59 / 1.45 / 1.39** |
| *for scale: PALM22E / PALM09E (raw, last map)* | 1,922 / 1,135 | 55.5 / 32.8 | | 75.3 / 41.9 | |

E_B grows less than linearly in c, because m\*(n) rises to 2 for well-covered cells.

Only the five-patient pool approaches E\*. It does so at the cost the pass bar has to weigh: a single-site miss carrying a share s of the pooled E shows up as only 1 + 0.5·s. So pooling would have **hidden PALM09E T3**.

## Deviations and limits

- **The proxy, not v3's n, decides evaluability here.** It is a strict upper bound, so "not evaluable" is final. "Borderline" (PALM07E T1) needs the real v3 n.
- **E_B is scaled linearly** by the central UMI ratio. The Poisson stand-in is slightly conservative relative to the Poisson-binomial.
- **P_s is fixed at its 1e-4 floor.** A target with a high indel background (long-insertion-like) would have a larger E but a less informative test.
- **The seeds are a stand-in** (see Method).
- **Resources:** coverage peaked at 18.3 GB and 11.6 min (PALM08E, 54.6 M reads, 40 G requested). The other samples used 3.3–8.4 GB and 2.4–7.1 min.
