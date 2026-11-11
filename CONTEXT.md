# NanoGen genotyping

Per-cell genotyping of targeted mutations from SCM-seq: 10x single-cell RNA-seq with targeted enrichment, sequenced on Oxford Nanopore.

## Language

### Comparing callers

**Call**:
The genotype assigned to one cell × target: MUT, WT or NA.
_Avoid_: genotype (when the per-cell × target result is meant)

**Hard flip**:
A call that changes between MUT and WT when two pipeline versions are compared.
_Avoid_: discordance, regression

**Soft flip**:
A call that changes to or from NA when two pipeline versions are compared.
_Avoid_: dropout, gain

**Attribution ladder**:
An ordered sequence of pipeline states, each adding one change on top of the previous one. Each step between states is a **rung**, and each rung's flips are attributed to the change it adds.
_Avoid_: ablation (an ablation removes changes one at a time from the final state)

### Calibration on pooled samples

**Pooled calibration test**:
The calibration test of negative-control calls summed over every target and patient in the census (total observed vs total expected), used when no single target has enough expected calls to be tested alone. A pass says nothing about any one target.
_Avoid_: per-target calibration, aggregate evaluability

**Site exception**:
A target flagged by the pre-registered per-target outlier screen that runs next to the pooled calibration test. It stays in the pool and doesn't fail the verdict; it is reported alongside it and investigated.
_Avoid_: outlier exclusion, dropped target
