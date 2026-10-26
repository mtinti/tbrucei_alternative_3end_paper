# From the raw reads to the inputs in `data/`

These steps ran on the University of Dundee compute cluster (Grid Engine
submission scripts, kept as they ran). They are not part of `make all`; their
outputs are the files in `data/` and in the Zenodo record.

| Step | Tool | Script | Output |
|---|---|---|---|
| 1. Basecalling | Guppy (`ont-guppy` 6.2.1 build in the script), flow cell FLO-MIN106, kit SQK-RNA002, `--trim_strategy none`, `--fast5_out` | `guppy_basecall.sh` | `<replicate>.pass.fastq.gz` and the basecalled fast5 files |
| 2. Alignment | minimap2 `-ax splice -uf -k14 -G 10000`, then `samtools sort` | `minimap_align.sh` | `<replicate>.bam` against `data/genome/tb927_56.fa` |
| 3. Poly(A) tail estimation | nanopolish polya on the basecalled fast5, fastq and BAM | `run_nanopolish.sh` | `nanopolish_polya_results.tsv` per replicate |
| 4. Spliced leader detection | a fuzzy search for the 3' end of the spliced leader in the first 50 bases of each read (up to 2 mismatches), keeping soft-clipped matches | notebook on the cluster (`find_sl.ipynb`) | a table of reads with the spliced leader |
| 5. Poly(A) read selection | nanopolish `qc_tag` PASS and `polya_length` > 25 | `01_select_polya_reads.py` | `<replicate>.filter_polyA_tb927.bam` and `data/read_ids/<replicate>.polyA_read_ids.txt` |
| 6. Full-length read selection | the poly(A) rule and a soft-clipped spliced leader match | notebook on the cluster | `<replicate>.filter_polyA_SL_tb927.bam` |
| 7. Merging per stage | `samtools merge`, then sort and index | `merge_bam.sh` | the four BAM files of the Zenodo record |
| 8. End sites | the 3' end of every primary alignment with at least 200 aligned bases, counted per position and strand | `02_polya_end_sites.py --bam <replicate>.filter_polyA_tb927.bam` | `data/end_sites/<replicate>.end_sites.bed` |

Raw runs: direct RNA nanopore sequencing, European Nucleotide Archive.

| Replicate | ENA run |
|---|---|
| bs_rep_1 | ERR10957783 |
| bs_rep_2 | ERR10957784 |
| bs_rep_3 | ERR10957785 |
| pc_rep_1 | ERR10957780 |
| pc_rep_2 | ERR10957781 |
| pc_rep_3 | ERR10957782 |
