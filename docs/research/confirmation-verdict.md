# Verdict: genotype-v3 vs v3-ub on the new diagnosis pool

Wayfinder ticket [Verdict: does genotype-v3 hold, or does v3-ub replace it?](https://github.com/dimadatascience/NanoGen/issues/90), map [Wayfinder: confirm genotype-v3 vs v3-ub on new diagnosis samples](https://github.com/dimadatascience/NanoGen/issues/72).

- **Decision rule:** §7 of [`confirmation-preregistration.md`](https://github.com/dimadatascience/NanoGen/blob/e6f629f/docs/research/confirmation-preregistration.md) at `e6f629f`.
- **Evidence:** [`confirmation-score.md`](https://github.com/dimadatascience/NanoGen/blob/6a9a9fd/docs/research/confirmation-score.md) at `6a9a9fd` ([Score v3 vs v3-ub against the confirmation pre-registration](https://github.com/dimadatascience/NanoGen/issues/89)). Every number below is quoted from it.
- **Patient data:** aggregate numbers only; targets as T1…Tn.

The document has two parts that must not be confused:
1. **The pre-registered verdict** (§1–§3): the outcome of the §7 rule. Nothing here was decided after the scores were seen.
2. **A post-hoc interpretation** (§4–§5): a two-mode reading of the same evidence, decided after the scores were seen. It is the handoff to the later code ticket. It is **not** an outcome of the pre-registered rule and does not change it.

## 1. Pre-registered verdict: rule 4. Neither contender is confirmed; v3 stays provisional

| criterion | v3 | v3-ub |
|---|---|---|
| §5.1 pooled calibration | **fail**: O/E 1.66, p 0.0002 | **pass**: O/E 0.63, p 1.00 |
| §6.3 precision gate (FDR ≤ 0.2, 8 gated loci) | pass (max 0.17) | pass (max 0.09) |
| §6.5 tumour side evaluable | yes | yes |
| passes (§7) | **no** | **yes** |
| §6.4 tumour-call loss vs v3 ≤ 5 % (replacement only) | n/a | **fail: 11.2 %** |

- Rule 1 (v3 holds) does not apply: v3 fails the pooled calibration.
- Rule 2 (v3-ub replaces v3) does not apply: v3-ub passes, but its 11.2 % loss is above the 5 % cap.
- Rule 3 (both pass) does not apply.
- **Rule 4 applies: neither contender is confirmed as the single forward setting, and v3 stays the provisional setting.**

"Provisional" (see `CONTEXT.md`) means v3 remains the operating default, used with its measured limits stated, not as a confirmed setting.

## 2. Labels (§7), none of which changes the verdict

- **Fragile: no.** Both outcomes hold in the n_A ≥ 11 rerun (v3 O/E 1.61, fail; v3-ub 0.63, pass).
- **Site exceptions: AML5 T11, for v3 only.** O_B 1 against E_B 0.0009 (Holm-adjusted p 0.037). v3-ub: O_B 0 against E_B 0. Its follow-up is in §6.
- **Weak-null pass: PALM14E {T3–T5}.** v3's FDR there is 0.081–0.089 under the point-estimate null, and 1.06–1.16 with the null raised to CPU(x_s = 1, n_s = 6).
- **Blind spot (§1, verbatim):** "On this census, neither the pooled test nor the site screen can rule out a single-site miss the size of PALM09E T3's (O/E ≈ 1.7). A pooled pass therefore means 'no average miss ≥ 1.5× and no gross site miss on well-covered targets', **not** 'no T3-like site misses'."
- **Scope line (§6.6):** "Per-target precision confirmed on 8 covered, signal-bearing loci in 4 patients. 24 under-covered and 0 no-signal targets are listed. Thin-null loci: PALM14E {T3–T5}. Weak-null pass: PALM14E {T3–T5}."

## 3. Stated limits

- **Per-target calibration is not confirmable on this census** and stays provisional for both contenders (§1). The calibration half is pooled: it speaks to average miscalibration on new patients, not to any one target.
- **The outlier flag excluded 20 of 57 pool targets.** That includes **all five of sAML1's tumour-evaluable targets** (c_s 0.26–0.61), AML5 T1 and PALM14E T2, so 7 of 17 tumour-evaluable targets. The tumour-side claim does not reach sAML1. Whether this reflects impure negative cells or tumour labels in the old-panel samples, or germline/artefact targets, is unresolved. It stays in the map's fog ("Flag reach and negative-cell purity").
- **The tumour side covers 8 gated loci in 4 patients** (PALM08E, PALM11E, PALM14E, AML5). AML4 is calibration-only; PALM07E and PALM21E have no WES.

## 4. Finding (explanatory, not a recommendation): where v3's excess comes from

- v3's miss is **not driven by one patient**. Its point O/E is above 1 in 7 of 8 patients, and every leave-one-patient-out pool still fails (O/E 1.52–1.85).
- The 90 % interval excludes 1 only at **PALM14E** (2.90) and **AML5** (5.06).
- The excess concentrates at **thin nulls**: targets where the A-half contamination estimate c_A is near 0 but B still reads MUT (PALM14E T6, AML5 T3, AML5 T11). There v3's held-out null falls to the P_s floor, so a single ALT UMI in a negative cell becomes a MUT call.
- The same thin nulls are where v3-ub pays its loss (PALM14E 26.7 %, AML5 26.9 %; PALM14E T1 79 %, {T3–T5} 65–78 %, T7 58 %, AML5 T3 27 %). At well-sampled nulls the two contenders agree within 3 %.

A rule tuned between the two at thin nulls (a floor or shrinkage on c, or v3-ub only at thin nulls) would be chosen on these data and could not be confirmed here. It is **out of scope** for this map and would be a fresh effort on fresh samples.

## 5. Post-hoc interpretation: two modes (handoff to `/to-spec`)

The pre-registration asked for **one** forward setting, and §1 answers that. The scores also show that the two contenders give different, separately documented guarantees. The project wants both a stringent mode and a relaxed one that genotypes more cells. This interpretation offers each contender as a mode. **It was decided after the scores were seen**, and it is reported as such.

| | **relaxed mode = v3** (default) | **strict mode = v3-ub** (opt-in) |
|---|---|---|
| status | provisional setting | passes the pre-registered pooled calibration out of sample (v3-ub was chosen in-sample on PALM09E T3) |
| guarantee | per-target precision gate (FDR ≤ 0.2) at the 8 gated loci | pooled calibration + per-target precision gate at the 8 gated loci |
| measured cost | null under-estimated ≈ 1.66× on average (LOPO 1.52–1.85): more false MUT calls in non-tumour cells than the model expects | over-corrects (O/E 0.63); 11.2 % fewer gated tumour calls (≈ 27 % at PALM14E and AML5); −22 % tumour calls and 7 targets silenced outside the gated set |
| thin-null caveat | at PALM14E {T3–T5} precision holds only under the point-estimate null (weak-null pass) | calls at thin nulls are largely suppressed |

- **Default: relaxed (v3).** It keeps current behaviour and keeps the cells and targets that v3-ub drops.
- **No selection rule.** The verdict states what each mode guarantees and costs; choosing a mode is up to the pipeline's user, not an outcome of this map.
- **Handoff spec (for `/to-spec`):** expose a mode switch in the genotyping step: relaxed (v3, current behaviour, default) and strict (v3-ub: P_s with the upper 95 % bound of c). Document each mode's guarantee, cost and the limits in §3, and mark the relaxed mode as provisional. v3-ub is not in the merged code; implementing and merging it belongs to that later code ticket, not to this map.
- **Indel matching** is decided separately in [Indel arm: is position-tolerant matching worth a code ticket?](https://github.com/dimadatascience/NanoGen/issues/91) and feeds `/to-spec` on its own.

## 6. Site exception follow-up (§5.3): AML5 T11

Recorded in the map's fog under "Outlier-flag revision", not as a ticket:
- **Germline is implausible.** In-sample c_s is 0.077 (one MUT UMI over 11 covered negatives), far below the ≈ 0.5 a germline heterozygous site would show.
- **An artefact cannot be decided from one UMI.** The hit is the thin-null mechanism of §4: in most seeds A holds no MUT, so a single ALT UMI in B is called.
- It joins PALM09E T3 as a second example for any per-site guard a revised flag would add.
