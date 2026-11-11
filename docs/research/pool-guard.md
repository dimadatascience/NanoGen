# Pooled calibration: guarding against a single-site miss

Ticket: [Underpowered census: stop, pool, or restart from FASTQ?](https://github.com/dimadatascience/NanoGen/issues/83), map [Wayfinder: confirm genotype-v3 vs v3-ub on new diagnosis samples](https://github.com/dimadatascience/NanoGen/issues/72).
Aggregate numbers only; targets as T1…Tn. HPC folder `q83_pool_guard/`, job 9063440 (`short`, 1 CPU, 2 GB, 1 min).

## Question

The primary calibration verdict pools all targets and patients (ΣO vs ΣE). A single-site miss with share s of the pooled E shows only as 1 + 0.5·s. Which guard should go next to the pooled test?

- **(a)** The pooled test alone.
- **(b)** The pooled test plus a per-target outlier screen that flags but doesn't veto.
- **(c)** A heterogeneity test that fails the verdict.

## Answer

**(b)**, with a one-sided exact Poisson tail per target and Holm at α = 0.10. A flagged target:
- stays in ΣO / ΣE;
- doesn't fail the verdict;
- is reported as a **site exception**.

On this census, though, **no guard can detect a T3-sized miss**. PALM09E T3 was O/E 1.70 on the held-out half (O = 12, E = 6.8, p = 0.048). Against a miss of that size, every option detects 3–19 %, compared with a 5 % false-alarm rate. The screen sees only r ≳ 5, and only on the larger-E targets.

## Method

[`pool_guard.py`](pool-guard/pool_guard.py) models O_i ~ Poisson(r_i·E_i) over the 34 targets with E > 0. E is the central-scale E_B from [`power-check/out/census.tsv`](power-check/out/census.tsv), extracted to [`eb.tsv`](pool-guard/eb.tsv). It runs 1e5 simulations per scenario.

- **Pooled verdict:** fails when p < 0.05 or O/E > 1.5.
- **(c):** G-tests across targets and across patients, with Monte-Carlo critical values conditional on ΣO. The asymptotic χ² is miscalibrated here: 0.2–1 % false alarm across targets, 6–7 % across patients.
- **Bad site:** placed on the largest-E target (2.6 at c = 1 %, about 11 % of pooled E) or on a median-E target (0.51).
- **Sensitivity row:** E × 1.27, for the regenerated BAMs.

Full results are in [`pool_guard_out.json`](pool-guard/pool_guard_out.json).

## Results, c = 1 %

Each r_site cell reads "largest-E target / median-E target".

| option | null false alarm | r=1.7 | r=3 | r=5 | r=10 | uniform 1.5× |
|---|---|---|---|---|---|---|
| (a) pooled alone | 0.049 | .10 / .06 | .28 / .08 | .64 / .11 | .99 / .25 | 0.70 |
| (b) pooled fail or site flagged (Holm α=.10) | 0.049 fail + 0.036 exception | .12 / .07 | .47 / .14 | .92 / .32 | 1.0 / .78 | 0.72 |
| (c) pooled or across-target G | 0.097 | .16 / .11 | .40 / .15 | .82 / .26 | 1.0 / .62 | 0.71 |

## Results, c = 2 %

| option | null false alarm | r=1.7 | r=3 | r=5 | r=10 | uniform 1.5× |
|---|---|---|---|---|---|---|
| (a) | 0.048 | .11 / .06 | .35 / .09 | .78 / .16 | 1.0 / .41 | 0.86 |
| (b) | 0.048 + 0.039 | .17 / .07 | .71 / .16 | .99 / .46 | 1.0 / .95 | 0.87 |
| (c) | 0.098 | .18 / .12 | .53 / .19 | .94 / .39 | 1.0 / .87 | 0.87 |

## Minimum E for the screen

This is the per-target E at which the screen (Holm α = 0.10, 34 targets) reaches 80 % power:

| r_site | 1.7 | 2 | 3 | 5 | 10 |
|---|---|---|---|---|---|
| min E | 34 | 18 | 5.7 | 1.9 | 0.7 |

No target here reaches E = 5.7. The maximum is 4.0 at c = 2 %, or 5.1 on the regenerated BAMs.

## Other findings

- Holm and BH agree to within 0.01, since there is only one bad site.
- An uncorrected screen flags at least one target in 50–64 % of null runs.

## Methods context

- **Pooled verdict plus a flag-don't-veto screen is standard practice.** It appears in provider SMR monitoring, where funnel plots use pre-set exact Poisson limits ([Spiegelhalter 2005](https://pubmed.ncbi.nlm.nih.gov/15568194/)), and in UK clinical-audit outlier policies, where a flagged unit is investigated rather than removed ([HQIP](https://www.hqip.org.uk/guidance/detection-and-management-of-outliers-for-national-clinical-audits/)).
- **No precedent was found for a heterogeneity test used as a veto.** Sum-type dispersion tests are weak against a single outlier and asymptotically invalid with most E < 1. The reference for the latter is the Cochran (1954) rule.
- **Choice of correction.** For 38 units, a 99.8 % funnel limit and Holm's first step are nearly the same. [Jones, Ohlssen & Spiegelhalter 2008](https://pubmed.ncbi.nlm.nih.gov/18226745/) argue for FDR in provider screening; here BH and Holm coincide.
