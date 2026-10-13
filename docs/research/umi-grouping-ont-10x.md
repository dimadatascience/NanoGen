# UMI grouping for ONT 10x reads: current tools vs `make_consensus.sh`

Researched 2026-10-01 for the ticket "UMI grouping approaches for Nanopore 10x long reads" (genotyping spec `261001_nanogen_genotype.pdf`, limitation 3).
Tags: **[src]** = read in the tool's source or docs at the pinned commit (links below); **[paper]** = stated in the cited paper; **[inf]** = my own inference, not yet checked on data.

## 1. What NanoGen does today

**Where UB comes from.** BLAZE assigns cell barcodes. It runs with `--max-edit-distance 1`, which is a *barcode* edit distance ([`blaze.nf`](../../subworkflows/long_reads/modules/blaze.nf), [`nextflow.config:11`](../../nextflow.config)). It writes `@<BC>_<UMI>#<readid>_<strand>` into the read name ([BLAZE `read_assignment.py:129`][blaze-ra]). `fix_tags.py` then copies `name[:29]` into CB and UB ([`fix_tags.py:21-24`](../../bin/long_reads/fix_tags.py)).
- BLAZE does **not** correct UMIs. Its UMI is the fixed 12 nt (10 nt for v2) window right after the barcode. The only adjustment is that the window's start moves to absorb indels in the *barcode* ([BLAZE README l.23, l.79][blaze-readme]; [`read_assignment.py:92-99`][blaze-ra]) [src]. The BLAZE paper describes barcode identification only ([doi:10.1186/s13059-023-02907-y][blaze-paper]) [paper].
- So UB is always exactly 12 nt (`fix_tags.py` hard-codes 16+1+12, so v2 kits would break it) [src]. An indel *inside* the UMI therefore shows up as a shifted tail: the bases after the indel slide by one, and a flanking or polyT base enters or leaves the window [inf, follows from fixed-window extraction].
- BLAZE re-orients reads to the transcript strand by default ([`parser.py:124-127`][blaze-parser]) [src].

**Pipeline vs script defaults.** The pipeline passes `window_width = 100` and `min_reads = 2` ([`nextflow.config:17,22`](../../nextflow.config)), not the script's own defaults of 1000 and 3. The genotyping spec says n_min = 3. That mismatch should be resolved whatever else is decided.

**What `umi_tools group … --method adjacency --edit-distance-threshold 1 --paired --umi-tag=UB --extract-umi-method=tag` actually does** ([`make_consensus.sh:114`](../../bin/long_reads/make_consensus.sh)):

1. **UMIs are only compared among reads with the same start position.** No `--per-gene` or `--per-contig` is given, so reads are bundled by `pos` together with `key = (strand, spliced, tlen*paired, read_length)` ([`sam_methods.py:478-510`][ut-sam]).
   - `pos` is the read's **5′ end, including soft clips**: alignment start minus the leading soft clip on the + strand, alignment end plus the trailing soft clip on the − strand ([`sam_methods.py:700-730`][ut-sam]) [src].
   - UMI clustering runs separately inside each bundle ([`group.py:213-235`][ut-group]) [src].
   - The 1000 bp in `check_output` only controls when buffered bundles are flushed. It is not a grouping window ([`sam_methods.py:301-310`][ut-sam]) [src].
   - **Effect [inf]:** after BLAZE re-stranding, the read's 5′ end is the variably truncated, TSO-side end of the cDNA, and ONT soft clips vary from read to read. Reads of one molecule therefore rarely share an exact `pos`. Most bundles would be singletons, the adjacency step would merge almost nothing, and BX would mostly equal UB. This is the key thing to measure (§4).
   - Exception [inf]: if SCM-seq enrichment uses a target-specific primer that fixes one read end, positions would be more stable. The protocol docs in the repo do not say.
2. **Distance is Hamming only.** `edit_distance` returns ∞ for unequal lengths and otherwise counts mismatches ([`_dedup_umi.pyx:2-19`][ut-pyx]). `UMIClusterer.__call__` asserts that all UMIs in a bundle have the same length ([`network.py:366-370`][ut-net]) [src]. umi_tools has no Levenshtein option for `group` or `dedup`. Its fuzzy indel matching exists only in `CellClusterer` for cell barcodes ([`network.py:448-500`][ut-net]) [src].
   - With fixed-window UMIs, one indel costs Hamming distance ≈ the number of bases after the indel, but Levenshtein distance ≤ 2 (one indel plus the swapped end base) [inf]. So even a Levenshtein-1 threshold would miss most indel siblings; **threshold 2 is the minimum that catches a single indel**.
3. **`--paired` on single-end reads is effectively a no-op, but it is wrong.** For a read without flag 0x1 it increments the "Unpaired reads" counter and continues, because `--unpaired-reads` defaults to `use` ([`sam_methods.py:331-350`][ut-sam]; [`Utilities.py:1008-1014`][ut-util]). In the bundle key, `tlen*paired` is 0 because TLEN = 0 for SE reads [src]. The mate-unmapped branch tests flag 0x8, which minimap2 does not set for SE reads [inf].
   - Result: grouping is the same with or without the flag, and only the log is misleading. Drop it.
4. **Splitting by `BX_chrom_int(POS/win)*win`** ([`make_consensus.sh:120-128`](../../bin/long_reads/make_consensus.sh)) uses BX, umi_tools' best UMI *per bundle*, not UG, its group ID. It also uses SAM POS, the leftmost aligned base.
   - So reads with identical BX at *different* start positions are merged again if they fall in the same window. The window effectively turns this into "exact UMI string within a window" grouping [inf].
   - Molecules are split when (a) the UMI carries any error, (b) POS crosses a window boundary, or (c) different bundles chose different representative UMIs [inf].
   - The split depends on strand [inf]. For genes on the + strand, POS is the variable 5′ (truncation) end, so 100 bp windows split often. For genes on the − strand, POS is near the polyA site, so it is fairly stable.
5. **Groups below the threshold are deleted** (`rm "$samfile"`, `make_consensus.sh:131-140`), and `count_consensus.py` applies `min_read` again per allele [src]. Secondary and supplementary alignments are filtered by neither umi_tools `group` nor the awk step, only by `-q` [src]. Supplementary records could add extra "reads" to a group [inf].

## 2. How other tools handle UMIs

| Tool | UMI error model / extraction | Distance | Clustering | Grouping key | Groups below threshold |
|---|---|---|---|---|---|
| **NanoGen now** | BLAZE fixed 12-nt window, no UMI correction | Hamming ≤1 | umi_tools adjacency | exact 5′ position + strand, then BX + 100 bp window on POS | groups with fewer than `min_reads` (2) reads deleted |
| **umi_tools** ([Smith 2017][ut-paper]) | substitutions only | Hamming, equal length required | adjacency / directional / cluster | position (default); `--per-gene --gene-tag=X` or `--per-contig` uses the tag or contig as the bundle key and ignores position ([`sam_methods.py:434-473`][ut-sam]) | no read threshold in `group` |
| **UMICollapse** ([Liu 2019][uc-paper]) | substitutions (`-k`, "substitution edits") | Hamming | dir / adj / cc | unclipped 5′ coordinate + strand ([`DeduplicateSAM.java:90-98`][uc-src]); UMI taken from the read name | none; collapses each group to one consensus read |
| **wf-single-cell** (successor of sockeye; [README l.10][wfsc-readme]) | UMI taken from a parasail alignment of an adapter + N + polyT probe, gaps removed, so length can vary ([`extract_barcode.py:199-214`][wfsc-ext]) | **Levenshtein ≤2** | umi_tools *directional*, monkey-patched for Levenshtein and to drop the equal-length assertion ([`create_matrix.py:65-111`][wfsc-cm]) | **gene × cell**; reads with no gene fall back to 1 kb windows on the alignment *midpoint* ([`create_matrix.py:148-185`][wfsc-cm]; [README l.337-347][wfsc-readme]) | singletons kept as UMIs; reads whose corrected UMI ≠ 12 nt are dropped. Its SNV branch then runs `umi_tools dedup --method unique --per-cell` *without* `--per-gene`, so it is position-based ([`snv.nf:79-85`][wfsc-snv]) |
| **sockeye** | upstream `nanoporetech/sockeye` returns 404 (checked 2026-10-01); wf-single-cell lists itself as its port ([CHANGELOG l.386][wfsc-cl]) | n/a | n/a | n/a | n/a |
| **BLAZE** | barcodes only; UMI is a fixed window | n/a for UMI | n/a | n/a | n/a |
| **FLAMES** ([Tian 2021][flames-paper]) | input from BLAZE or flexiplex | Levenshtein ≤1 (`fast_edit_distance`) | greedy: most abundant UMI absorbs neighbours, rarest first ([`count_gene.py:605-635`][flames-cg]) | **cell × gene**, optionally split by 3′ (polyT-side) position clusters 50 bp apart ([`count_gene.py:552-564, 589-603`][flames-cg]) | all UMIs counted; one read kept per UMI |
| **Sicelore 2.1** ([Lebrigand 2020][sic-paper]) | UMIs searched in the read | "ED ≤ 2" | cluster centre or highest-quality UMI | cell × genomic region (reads within 500 nt) ([README l.581-597][sic-readme]) | singletons kept; consensus is the read itself (n = 1), the best read (n = 2) or spoa (n > 2); `MINRN` default 0 ([README l.1138-1143, 1449][sic-readme]) |
| **scNanoGPS** ([Shiau 2023][gps-paper]) | UMI from the read | Levenshtein ≤1 (`--umi_ld`); the README says 2 | pairwise merge, no count rule ([`curator_io.py:141-203`][gps-cur]) | cell × reads whose `reference_start` is within 5 bp of the block's first read ([`curator_io.py:88-113`][gps-cur]); soft-clip filter applied first | singletons kept and merged into the final BAM; spoa consensus for ≥ 2 reads |
| **IsoQuant** | UMI tag | Levenshtein ≤3 (constructor default) | greedy by abundance | **gene × barcode** ([`umi_filtering.py:145-267, 586-593`][iq-umi]) | one representative read per molecule |

**Error rates for context.** In amplicon work, Karst et al. needed 15× (R10.3) to 25× (R9.4.1) reads per UMI for consensus error below 0.01% ([doi:10.1038/s41592-020-01041-y][karst]) [paper]. NanoGen's consensus with n_min = 2-3 only needs a per-site majority vote, so its per-UMI calls are much noisier than that [inf].

**The pattern [inf from the table].** Every ONT-specific single-cell tool groups by **cell × gene** (or a coarse region of ≥ 500 nt, or the 3′ end) rather than by exact 5′ position. Every one except BLAZE uses **Levenshtein** distance, and none of them discards single-read molecules at the grouping stage.

## 3. Options for NanoGen

| # | Option | Gain | Cost / risk |
|---|---|---|---|
| A | **Drop `--paired`** | correct semantics and logs | none (no change in grouping) |
| B | **Group per target instead of per window or position.** Tag each read with its BED target ID or gene (BED col 6), then run `umi_tools group --per-gene --gene-tag=XT` (or group by `(CB, target)` in Python) and split by UG, not BX | removes splitting from position and window boundaries; one group per molecule per target; per-site counting is naturally per target | a read covering two targets needs a rule (duplicate it per target, or use the gene); `--per-gene` holds a whole contig in memory (fine for per-cell BAMs) |
| C | **Levenshtein ≤2, directional**, as in wf-single-cell: a ~40-line Python step reusing the `UMIClusterer` patch, keyed on (CB, target) | merges indel siblings, which Hamming cannot | false merges. Rough estimate: about 10³ 12-mers lie within LD ≤ 2 of a UMI, so the chance per pair is ~10⁻⁴; with k molecules per cell × target that is about k²/2·10⁻⁴ (≈ 0.1 at k = 50) [inf]. A false merge of a WT and a MUT molecule produces a mixed group, which f_min then discards or which turns into OTH |
| D | Keep Hamming but go to `--edit-distance-threshold 2`+ | trivial change | still blind to indels (shifted tail), more false merges from substitutions; not recommended |
| E | Re-extract UMIs by alignment with flanks (wf-single-cell probe) or use flexiplex/BLAZE flanks, so a single indel becomes LD 1 | cleaner distances, tighter threshold | larger change; replaces part of BLAZE output |
| F | **Groups below n_min**: (i) keep discarding *after* B and C (most fragments should then be merged); (ii) keep 1-read UMIs as a separate low-confidence class, as Sicelore, scNanoGPS, wf-single-cell and FLAMES do, and report them; (iii) let them count as WT or MUT only when base quality is high | (i) cheapest; (ii)/(iii) recover WT sensitivity (spec limitation 2) | (ii)/(iii) raise MUT calls driven by single-read errors; the spec's error model (ε floor) would need re-checking |
| G | Filter secondary and supplementary alignments (`-F 0x900`) before grouping | stops one read being counted twice | a few chimeric reads lost |

**Recommendation [inf]:** A + G now (safe), then **B + C** as the main change, then decide F(i) vs F(ii) using the data below. Settle the n_min mismatch (2 in the config vs 3 in the spec) at the same time.

## 4. Measurements that decide (sibling ticket "Measure UMI splitting and window double-counting on real data")

1. **umi_tools bundle sizes**: run `--group-out` and check the share of bundles with a single read and the BX ≠ UB rate. A high singleton share means adjacency is doing nothing, which supports B.
2. **Splitting per (CB, target, exact UB)**: number of distinct windows (100 bp and 1 kb) and distinct umi_tools UG values. This gives the double-counting rate directly. Split it by gene strand to test the POS asymmetry.
3. **UMI distance spectrum within (CB, target)** vs a null of UMI pairs drawn from *different* cells for the same target: count pairs at Hamming 1, LD 1, LD 2 and LD 3. Pairs above the null are error siblings; the null rate is the false-merge rate. Together they fix the threshold for C.
4. **Allele concordance of candidate merges**: at the target site, how often do LD ≤ 2 siblings disagree (WT vs MUT) in cells with both alleles? Discordance means false merges.
5. **What n_min removes**: fraction of reads and UMIs dropped below n_min under current, B, and B + C grouping. Then compare WT/MUT/undetermined calls per cell between the methods.
6. Fraction of secondary and supplementary records in the per-cell BAMs (for G).

## Sources

All repos read at the commit in the link (cloned 2026-10-01).

- UMI-tools @c457eef: [sam_methods.py][ut-sam], [network.py][ut-net], [_dedup_umi.pyx][ut-pyx], [group.py][ut-group], [Utilities.py][ut-util]. Smith T et al. 2017, Genome Res, [doi:10.1101/gr.209601.116][ut-paper] (PMID 28100584).
- UMICollapse @efeab35: [README](https://github.com/Daniel-Liu-c0deb0t/UMICollapse/blob/efeab35f5d29dec1d496ade3f681eeb34d9c2057/README.md) (l.76, 107, 113), [DeduplicateSAM.java][uc-src]. Liu D 2019, PeerJ, [doi:10.7717/peerj.8275][uc-paper] (PMID 31871845).
- wf-single-cell @cb39932: [README][wfsc-readme], [CHANGELOG][wfsc-cl], [create_matrix.py][wfsc-cm], [extract_barcode.py][wfsc-ext], [snv.nf][wfsc-snv].
- BLAZE @af8d2cc: [README][blaze-readme], [read_assignment.py][blaze-ra], [parser.py][blaze-parser]. You Y et al. 2023, Genome Biol, [doi:10.1186/s13059-023-02907-y][blaze-paper] (PMID 37024980).
- FLAMES @164abac: [count_gene.py][flames-cg]. Tian L et al. 2021, Genome Biol, [doi:10.1186/s13059-021-02525-6][flames-paper] (PMID 34763716).
- Sicelore 2.1 @33209a0: [README][sic-readme]. Lebrigand K et al. 2020, Nat Commun, [doi:10.1038/s41467-020-17800-6][sic-paper] (PMID 32788667).
- scNanoGPS @8095d28: [curator_io.py][gps-cur], [preprocessing.py](https://github.com/gaolabtools/scNanoGPS/blob/8095d282ccb8250fb5d3ac9a05b7a4c57a3c2936/curator_core/preprocessing.py#L36). Shiau CK et al. 2023, Nat Commun, [doi:10.1038/s41467-023-39813-7][gps-paper] (PMID 37433798).
- IsoQuant @7dd4406: [umi_filtering.py][iq-umi].
- Karst SM et al. 2021, Nat Methods, [doi:10.1038/s41592-020-01041-y][karst] (PMID 33432244).
- Metadata for the papers above retrieved from PubMed.

[ut-sam]: https://github.com/CGATOxford/UMI-tools/blob/c457eefc5e0c47368cf3652b9b6d497b78c3050a/umi_tools/sam_methods.py
[ut-net]: https://github.com/CGATOxford/UMI-tools/blob/c457eefc5e0c47368cf3652b9b6d497b78c3050a/umi_tools/network.py
[ut-pyx]: https://github.com/CGATOxford/UMI-tools/blob/c457eefc5e0c47368cf3652b9b6d497b78c3050a/umi_tools/_dedup_umi.pyx
[ut-group]: https://github.com/CGATOxford/UMI-tools/blob/c457eefc5e0c47368cf3652b9b6d497b78c3050a/umi_tools/group.py
[ut-util]: https://github.com/CGATOxford/UMI-tools/blob/c457eefc5e0c47368cf3652b9b6d497b78c3050a/umi_tools/Utilities.py
[ut-paper]: https://doi.org/10.1101/gr.209601.116
[uc-src]: https://github.com/Daniel-Liu-c0deb0t/UMICollapse/blob/efeab35f5d29dec1d496ade3f681eeb34d9c2057/src/umicollapse/main/DeduplicateSAM.java#L90-L98
[uc-paper]: https://doi.org/10.7717/peerj.8275
[wfsc-readme]: https://github.com/epi2me-labs/wf-single-cell/blob/cb39932e3dfca08ad6dae1e6d08579c957874fcf/README.md
[wfsc-cl]: https://github.com/epi2me-labs/wf-single-cell/blob/cb39932e3dfca08ad6dae1e6d08579c957874fcf/CHANGELOG.md
[wfsc-cm]: https://github.com/epi2me-labs/wf-single-cell/blob/cb39932e3dfca08ad6dae1e6d08579c957874fcf/bin/workflow_glue/create_matrix.py
[wfsc-ext]: https://github.com/epi2me-labs/wf-single-cell/blob/cb39932e3dfca08ad6dae1e6d08579c957874fcf/bin/workflow_glue/extract_barcode.py#L167-L227
[wfsc-snv]: https://github.com/epi2me-labs/wf-single-cell/blob/cb39932e3dfca08ad6dae1e6d08579c957874fcf/subworkflows/snv.nf#L79-L85
[blaze-readme]: https://github.com/shimlab/BLAZE/blob/af8d2cc9f6374caf4d1006d5f8a0e99097211646/README.md
[blaze-ra]: https://github.com/shimlab/BLAZE/blob/af8d2cc9f6374caf4d1006d5f8a0e99097211646/blaze/read_assignment.py
[blaze-parser]: https://github.com/shimlab/BLAZE/blob/af8d2cc9f6374caf4d1006d5f8a0e99097211646/blaze/parser.py#L124-L127
[blaze-paper]: https://doi.org/10.1186/s13059-023-02907-y
[flames-cg]: https://github.com/mritchielab/FLAMES/blob/164abaced4e0838bf7cc167a45d08d819aa4907f/inst/python/count_gene.py
[flames-paper]: https://doi.org/10.1186/s13059-021-02525-6
[sic-readme]: https://github.com/ucagenomix/sicelore-2.1/blob/33209a00e859037af1baefe9b69047c1d4ca9b86/README.md
[sic-paper]: https://doi.org/10.1038/s41467-020-17800-6
[gps-cur]: https://github.com/gaolabtools/scNanoGPS/blob/8095d282ccb8250fb5d3ac9a05b7a4c57a3c2936/curator_core/curator_io.py
[gps-paper]: https://doi.org/10.1038/s41467-023-39813-7
[iq-umi]: https://github.com/ablab/IsoQuant/blob/7dd4406f1799573df1505a2b76ef2e48d9455ed5/isoquant_lib/barcode_calling/umi_filtering.py
[karst]: https://doi.org/10.1038/s41592-020-01041-y
