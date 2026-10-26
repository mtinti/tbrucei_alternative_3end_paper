# Methods

Each step of the analysis with its parameters, in the order the data flow.
The script of each step is in `scripts/`; `run_analysis.sh` runs them in this
order and `make all` wraps it. The run that produced the published tables is
recorded in `expected/run_manifest.json` (inputs with checksums, options,
software versions) and `expected/run_manifest_interpretation.json`.

Terms. A **poly(A) read** is a nanopore read with a measured poly(A) tail; its
3' end on the genome is the poly(A) site of the molecule it came from. A
**full-length read** also carries the spliced leader at its 5' end, so both
ends of the transcript are known. The **window** of a gene is the span in which
its read ends are counted. An **end cluster** is a range where ends pile up.
**Usage** of a cluster is the share of the window's ends that fall in it. A
**supported** gene or cluster passes stageR at a gene-level overall false
discovery rate (OFDR) of 5 percent.

## 0. Samples and inputs (`upstream/`)

Bloodstream-form (BSF) and procyclic-form (PCF) *T. brucei*, three biological
replicates each, sequenced by direct RNA nanopore sequencing (kit SQK-RNA002,
flow cell FLO-MIN106; ENA ERR10957780 to ERR10957785). Reads were basecalled
with Guppy and aligned with minimap2 (`-ax splice -uf -k14 -G 10000`) to the
TREU927 genome (TriTrypDB release 56). Nanopolish estimated the poly(A) tail
of every read; reads with `qc_tag` PASS and a tail longer than 25 nt are the
poly(A) reads. A fuzzy search for the spliced leader in the first 50 bases of
each read (up to 2 mismatches, soft-clipped match) marked the full-length
reads. The replicates were merged per stage into one poly(A) BAM and one
full-length BAM; the read names of each replicate are kept
(`data/read_ids/`). The 3' end of every primary alignment with at least 200
aligned bases was counted per position and strand, per replicate
(`data/end_sites/`). `upstream/README.md` lists the scripts.

The annotation (`data/annotation/merged_fixed_6.gtf.gz`) combines TriTrypDB
gene models with transcript models curated by hand in the Apollo genome
browser, which include the untranslated regions. The CDS of a gene is the
outermost span of its CDS records, as a half-open zero-based interval. The
annotated non-coding RNA genes come from Ensembl Protists release 116.

## 1. Windows and the gene family (`scripts/16`)

No stage label is read here or in section 2.

- **Nanopore window**: from the CDS 5' end, through the CDS, to the far end
  of the longest eligible full-length read plus 100 bases, in transcript
  direction. A full-length read is eligible when it lies on the gene's strand,
  spans the whole CDS and touches no other CDS. Ends inside the CDS are
  counted: some genes end a stage-specific share of their transcripts inside
  the annotated coding sequence.
- **Intergenic window**: from the CDS 5' end to the nearest CDS boundary
  downstream on either strand. A second definition of the same question,
  used as a sensitivity analysis of the window choice.
- **Family**: a gene enters the analysis when every replicate has at least 10
  poly(A) read ends in its window. The family is fixed by depth only, never
  by signal. Nanopore windows: 2891 genes; intergenic windows: 3182.

## 2. End clusters (`scripts/17`)

For each window, the end positions of all six replicates are pooled.

1. Per-base counts are smoothed with a centred 11-base running sum. A
   candidate peak is a local maximum holding at least 2 percent of the
   window's ends and at least 5 ends. On a plateau the centre base is the
   peak.
2. Peaks closer than 15 bases are one site. Peaks 15 to 30 bases apart merge
   when the smoothed count between them never drops below half the lower
   peak. Peaks further apart never merge. Merging chains.
3. Each cluster core (first to last merged peak) is extended by 15 bases on
   each side, cut at the midpoint between neighbouring cores where two
   extensions would overlap, and clipped at the window. Clusters are
   disjoint, and for every replicate the cluster counts plus the ends in no
   cluster (**other**) equal the window total. The script stops if this
   identity fails.
4. A cluster is **supported** when it holds at least 5 ends in at least 2
   replicates. Unsupported clusters fold into other before any fit.
5. A window is testable with two or more supported clusters.

Primary setting: merge 30 bases, peak 2 percent. Two sensitivity settings run
in full: merge 50 bases, peak 5 percent; merge 20 bases, peak 1 percent. The
script has a self-test on synthetic windows (two near peaks, a plateau, a
peak at the window edge, both strands).

## 3. Differential usage test (`scripts/18`, `scripts/09`, `scripts/19`, `scripts/20`)

- **DEXSeq.** For every feature (each supported cluster and other) the counts
  of this feature against the rest of its window; full design
  `~ sample + exon + stage:exon`, reduced design `~ sample + exon`.
  Dispersions per feature, shrunk toward the trend. DEXSeq's independent
  filtering is off, because sections 1 and 2 already filtered. A feature whose
  fit gives no p value is kept with p = 1 and flagged.
- **stageR.** DEXSeq's per-gene q value screens the genes
  (`pScreenAdjusted = TRUE`); the features' unadjusted p values are the
  confirmation stage. Overall false discovery rate 5 percent at the gene
  level. Nanopore windows and intergenic windows are two separate families.
- **Effect size.** For every cluster: the raw usage fraction per replicate,
  the mean per stage and their difference in percentage points (PCF minus
  BSF). DEXSeq's fold change is reported as a change in relative usage on the
  log scale.
- **Sensitivity.** DRIMSeq on the same observed counts (no pseudocount), and
  the two other cluster settings. A gene is "supported in all analyses" when
  the primary, both other settings and DRIMSeq support it and all four
  could test it; the table also records whether the same site changes in the
  same direction under the other settings.
- **Permutation diagnostic.** Of the 10 distinct three-against-three
  partitions of the six samples, one is the true stage split; the other nine
  mix the stages. All ten are run with fixed cluster definitions and the
  dispersions refitted. The nine mixed counts show how many genes a wrong
  labelling produces; they are not an estimate of the false discoveries in
  the true comparison, because the mixed splits still contain real stage
  effects.

Result: 230 supported genes with nanopore windows (388 supported clusters)
and 352 with intergenic windows, 227 in both.

## 4. Interpretation of every supported cluster (`scripts/21`)

One row per supported cluster of a supported gene, from the reads that end
in it (primary alignment, at least 200 aligned bases, the same reads the test
counted).

- An **informative** read is one whose alignment start is a transcript 5'
  end: it is in the full-length BAM, or it starts within 100 bases of the
  annotated transcript start. Every other read is incomplete.
- A **linked** read is an informative read that starts within 100 bases of
  the gene's transcript 5' end, or upstream of the CDS 5' end.
- **Transcript evidence** at the cluster: a record of the independent
  transcript catalogue (section 5) or an annotated non-coding RNA ending
  within 60 bases on the same strand, or a trans-splice acceptor cluster
  downstream of the CDS 3' end (a 100-base bin of read starts with at least
  10 reads and 15 percent of the reads ending in the cluster, with an AG
  preceded by a pyrimidine-rich stretch within 60 bases upstream).

Label, in this order: **unresolved, insufficient evidence** with fewer than 20
reads or fewer than 10 informative reads; **mixed** with transcript evidence
and at least a third of the informative reads linked; **transcript-associated**
with transcript evidence and fewer than a third linked; **gene-connected** with
no transcript evidence and at least two thirds linked; **unresolved,
incompatible evidence** otherwise. The thirds are documented heuristics fixed
before the labels were computed. A gene has no class; labels belong to
clusters. A sensitivity run also counts a read as informative when at least
5 other reads ending in the same cluster start within 10 bases of it
(`events_*_sharpstart_sensitivity.tsv`).

Nanopore windows: 291 gene-connected, 39 mixed, 5 transcript-associated, 53
unresolved clusters.

## 5. Independent transcripts and their abundance (`scripts/22`, `scripts/23`)

A **transcript record** is a strand, a start range and an end range. Records
come from the full-length reads that touch no CDS, grouped when their starts
lie within 100 bases and their ends within 60 bases, and from the annotated
non-coding RNAs (same tolerances to merge an annotated and a read-derived
record). Records that share any ambiguous read are joined into one locus, so
counting units are disjoint. Records with at least 5 full-length reads in at
least 2 replicates are tested with edgeR (quasi-likelihood F test, stage the
only factor, TMM library sizes computed from all protein-coding genes and all
records together), Benjamini-Hochberg over the records tested. 68 records
were tested and 16 change in relative abundance within the full-length read
population. This is the only test behind any statement that an independent
transcript changes abundance.

## 6. Ranking, summary and figures (`scripts/26`, `scripts/24`, `scripts/11`, `scripts/27`, `scripts/25`)

Supported genes are ranked by their strongest cluster: tier 1 when the
replicates are separated (every replicate of the stage with higher usage
above every replicate of the other), the change is 30 points or more, the
gene is supported in all analyses and the label is gene-connected or mixed;
tier 2 with 20 points; tier 3 with 20 points or separated; tier 4 the rest.
The summary (`RESULTS_v2.md`) is written from the saved tables of one run.
Two figures per gene: an overview of the locus with coverage, reads and
clusters, and a zoomed view of the tested region with the tested reads only,
reads of uncertain origin (section 7) in grey. Every figure reads 5' to 3'
from left to right; for minus-strand genes the x axis is reversed.

## 7. Gene-assigned end usage and the strong-evidence count (`scripts/31`, `scripts/28`, `scripts/29`, `scripts/32`)

The window test compares a cluster with the rest of its window, and a window
can hold ends of transcripts that are not the gene's. This section repeats the
test on the reads that carry the gene.

- **Assignment rule** (stage-blind, fixed). A poly(A) read is assigned to a
  gene when (a) at least 50 of its aligned bases lie inside the gene's CDS
  (aligned bases, so a deletion or a skipped stretch does not count); (b) its
  placement is resolved: mapping quality 20 or more (the Phred convention of
  the SAM specification, a 1 percent chance of a wrong placement) and no
  secondary alignment anywhere in the genome with an alignment score equal to
  or higher than the primary's (minimap2 chooses the primary by chain score,
  so such a secondary is a placement at least as good; `scripts/31` collects
  the scores in one pass over the BAMs); (c) no other CDS on the same strand
  receives 50 aligned bases from it. Every other read is of uncertain origin
  and is counted separately with its reason (no CDS overlap, low mapping
  quality, equal placement elsewhere, covers another CDS), never in the
  denominator. In the nanopore family 821678 reads are assigned and 379698
  set aside.
- **Same features, same model.** The clusters of section 2 are kept. Counts
  per cluster and replicate are the assigned reads ending in it; other is the
  assigned ends outside the clusters. A cluster is supported with 5 assigned
  ends in 2 replicates. A gene is **eligible** with 10 assigned ends in every
  replicate and two or more supported clusters (M = 742 nanopore, 762
  intergenic). DEXSeq with stageR as in section 3: 138 and 147 supported
  genes.
- **Strong evidence, one rule.** A gene is flagged when at least one of its
  supported clusters satisfies all of: gene and cluster supported; the
  gene-assigned usage changes by 30 points or more between stages (mean of
  the raw replicate fractions, computed from the count matrix so the
  threshold applies to the exact value); the replicates are separated; a
  supported cluster overlapping it changes in the same direction under both
  other cluster settings, with the gene supported there; at least 20 assigned
  ends in the cluster. The rule runs over every supported cluster and the
  gene is counted once. X = 52 of 742 (nanopore windows); 54 of 762 with
  intergenic windows, 52 in both.
- **Sensitivity to the placement rule** (`scripts/32`). The whole chain was
  refitted on the nanopore family under mapping quality 1, 10, 30 and 60 with
  the secondary rule, and under 20 without it; the thresholds were fixed
  before any was run. X is 51 to 52; 50 of the 52 genes are flagged under
  every variant (`mapq_sensitivity_nanopore.tsv`).

What the flag means. It is a fixed evidence filter on top of the stageR
result and carries no error rate of its own; the nominal 5 percent belongs
to the stageR procedure. The window and gene-assigned analyses use the same
reads, so their agreement is sensitivity evidence, not replication. An end
inside the CDS counts as an alternative 3' end, not as an alternative 3'
UTR. The result concerns these six libraries.

## Software

Docker image (`Dockerfile`): Python 3.11 with pandas, numpy, pysam and
matplotlib; R 4.5 with DEXSeq 1.56.0, stageR 1.32.0, DRIMSeq 1.38.0 and
edgeR 4.8.2. The exact versions of the published run are in
`expected/run_manifest.json`. Upstream, on the cluster: Guppy, minimap2,
nanopolish and samtools (`upstream/README.md`).
