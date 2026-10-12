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

### UMI consensus

**Usable read**:
A read in a UMI cluster that carries an allele at the target site. At an SNV: base quality ≥ 20, not N, and no deletion over the site. At an indel target: any read classified as REF, ALT (tolerant match) or OTH, with no base-quality filter. The usable reads of a cluster are r + a + o.
_Avoid_: cluster read, supporting read (when only reads with a site allele are meant)

**Consensus label**:
The allele (WT, MUT or OTH) assigned to a UMI cluster. A label is assigned when that allele has at least n_min usable reads and a fraction of at least f_min of r + a + o. The rule is the same for SNV and indel targets.
_Avoid_: UMI call

**Sub-n_min cluster**:
A UMI cluster covering the target with fewer than n_min usable reads. It is reported (`subnmin_umis`, `subnmin_reads`) but never labelled.
_Avoid_: small UMI (cluster size alone does not decide it)

**Discarded cluster**:
A UMI cluster with at least n_min usable reads where no allele reaches both n_min and f_min. Every cluster covering a target is exactly one of: sub-n_min, labelled or discarded.
_Avoid_: ambiguous UMI
