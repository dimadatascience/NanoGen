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

### Negative controls

**Negative-control cell**:
A cell of a type that cannot carry the tumour mutation (T, B or NK), used to measure how often a wild-type cell reads MUT.
_Avoid_: normal cell, background cell

**Tumour compartment**:
The cells annotated as leukaemic in a sample. Ambiguous cells belong neither to it nor to the negative-control cells.
_Avoid_: blasts (when the annotated group is meant)

**Clonal target**:
A target whose mutation is carried by essentially every cell of the tumour compartment, so a sensitive caller should call most covered tumour cells MUT.
_Avoid_: driver, truncal

**Held-out split**:
A random split of a sample's negative-control cells into two halves: one half estimates the contamination rate, and the other half tests whether the caller's false MUT calls match what that rate predicts.

**Pass bar**:
The pre-registered criteria (calibration in negative-control cells, precision and sensitivity in the tumour compartment) a caller must meet on a sample. A target that lacks the data to test a criterion is *not evaluable* for it, which is distinct from failing it.

**Clone fraction**:
The fraction of covered cells in a compartment that carry a target's mutation, estimated per target from calls.
_Avoid_: VAF, allele frequency (those are per-molecule, not per-cell)
