# Data

Inputs of the analysis. `checksums.md5` lists the MD5 of each file; for a
`.gz` file it is the MD5 of the decompressed content. `make check-data`
verifies them. License: CC BY 4.0 (`LICENSE` in this folder).

`upstream/` describes how these files were made from the raw sequencing
runs (ENA ERR10957780 to ERR10957785).

## annotation/

| File | Content |
|---|---|
| `merged_fixed_6.gtf.gz` | Gene models of *T. brucei* TREU927 on the 11 megabase chromosomes: TriTrypDB gene models (source `VEuPathDB`) and transcript models curated by hand in the Apollo genome browser (source `.`), which include the untranslated regions. The analysis uses the CDS records (one outermost span per gene) and the transcript records (the annotated transcript 5' end). |
| `ensembl116_noncoding_genes.bed` | Annotated non-coding RNA genes from Ensembl Protists release 116 (BED, biotype in column 5). Used to interpret end clusters and to seed the catalogue of independent transcripts. |

## genome/

| File | Content |
|---|---|
| `tb927_56.fa.gz` | The reference the reads were aligned to: *T. brucei* TREU927 genome, TriTrypDB release 56 (11 chromosomes and unassigned contigs) plus Tb427 VSG sequences. 196 sequences. |

## end_sites/

| File | Content |
|---|---|
| `<replicate>.end_sites.bed.gz` | One line per genomic position and strand where poly(A) reads of the replicate end: chromosome, start, end, `polyA_site`, number of reads ending there, strand. Primary alignments with at least 200 aligned bases; the 3' end is the alignment end on the plus strand and the alignment start on the minus strand (`upstream/02_polya_end_sites.py`). |

## read_ids/

| File | Content |
|---|---|
| `<replicate>.polyA_read_ids.txt.gz` | The names of the reads of the replicate with a measured poly(A) tail (nanopolish `qc_tag` PASS, tail longer than 25 nt; `upstream/01_select_polya_reads.py`). The BAM files are merged per stage, and these lists give each read its replicate. |

Replicates: `bs_rep_1` to `bs_rep_3` (bloodstream form), `pc_rep_1` to
`pc_rep_3` (procyclic form).

## bam/ (from the Zenodo record, not in the GitHub repository)

| File | Size | Content |
|---|---|---|
| `bs.filter_polyA_tb927.sorted.bam`, `pc.filter_polyA_tb927.sorted.bam` | 1.5 GB, 1.0 GB | Poly(A) reads of the three replicates of each stage, aligned with minimap2 to `genome/tb927_56.fa`. Secondary alignments are kept (the analysis reads their scores). |
| `bs.filter_polyA_SL_tb927.sorted.bam`, `pc.filter_polyA_SL_tb927.sorted.bam` | 120 MB, 124 MB | The subset of those reads that also carry the spliced leader at their 5' end: full-length transcripts. |

Each BAM comes with its `.bai` index. `bam/checksums.md5` lists their MD5;
`make check-data` verifies them when they are present.
