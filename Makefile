# Reproduce the analysis. Run inside the Docker image (see README.md):
#
#   make all          check the inputs, run the whole analysis without figures, compare with expected/
#   make check-data   checksums of the input files (and of the BAM files when present)
#   make inputs       decompress the inputs into results/inputs and index the genome
#   make analysis     the analysis (run_analysis.sh), about 40 minutes with 4 threads
#   make verify       compare the tables in results/ with the published ones in expected/
#   make figures      the analysis again with the figures and the HTML indexes, about 2 hours
#   make zenodo       the data archive for Zenodo (BAM files, full results and figures)
#
# The four BAM files (data/bam/, from the Zenodo record) are needed by make analysis.

PYTHON ?= python
THREADS ?= 4
OUT ?= results
IN = $(OUT)/inputs
REPS = bs_rep_1 bs_rep_2 bs_rep_3 pc_rep_1 pc_rep_2 pc_rep_3

export GTF = $(IN)/merged_fixed_6.gtf
export GENOME = $(IN)/tb927_56.fa
export FAI = $(IN)/tb927_56.fa.fai
export NCRNA_BED = data/annotation/ensembl116_noncoding_genes.bed
export SITES_DIR = $(IN)/end_sites
export READ_IDS_DIR = $(IN)/read_ids
export BSF_POLYA_BAM = data/bam/bs.filter_polyA_tb927.sorted.bam
export PCF_POLYA_BAM = data/bam/pc.filter_polyA_tb927.sorted.bam
export BSF_FL_BAM = data/bam/bs.filter_polyA_SL_tb927.sorted.bam
export PCF_FL_BAM = data/bam/pc.filter_polyA_SL_tb927.sorted.bam
export OUT THREADS

.PHONY: all check-data inputs analysis verify figures zenodo clean

all: check-data inputs analysis verify

check-data:
	$(PYTHON) scripts/verify.py inputs

inputs:
	mkdir -p $(IN)/end_sites $(IN)/read_ids
	gzip -dc data/annotation/merged_fixed_6.gtf.gz > $(GTF)
	gzip -dc data/genome/tb927_56.fa.gz > $(GENOME)
	$(PYTHON) -c "import pysam; pysam.faidx('$(GENOME)')"
	for r in $(REPS); do \
		gzip -dc data/end_sites/$$r.end_sites.bed.gz > $(IN)/end_sites/$$r.end_sites.bed; \
		gzip -dc data/read_ids/$$r.polyA_read_ids.txt.gz > $(IN)/read_ids/$$r.polyA_read_ids.txt; \
	done
	@for f in $(BSF_POLYA_BAM) $(PCF_POLYA_BAM) $(BSF_FL_BAM) $(PCF_FL_BAM); do \
		[ -f $$f ] || { echo "missing $$f: download the BAM files from the Zenodo record into data/bam/ (see data/README.md)"; exit 1; }; \
	done

analysis:
	FIGURES=0 bash run_analysis.sh

figures:
	FIGURES=1 bash run_analysis.sh
	$(PYTHON) scripts/verify.py results $(OUT)

verify:
	$(PYTHON) scripts/verify.py results $(OUT)

zenodo:
	bash scripts/make_zenodo_archive.sh $(OUT) zenodo

clean:
	rm -rf $(OUT)
