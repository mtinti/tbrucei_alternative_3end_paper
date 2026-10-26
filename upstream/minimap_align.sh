#!/bin/bash
#$ -adds l_hard local_free 500G
#$ -mods l_hard m_mem_free 20G
#$ -adds l_hard avx 1 
#$ -cwd
#$ -V
#$ -j y
#$ -pe smp 24
#$ -o errors_minimap_$JOB_ID



# Exit immediately if a command exits with a non-zero status
set -e

# Check if TMPDIR is set
if [ -z "$TMPDIR" ]; then
    echo "Error: TMPDIR is not set"
    exit 1
fi

# Check if two arguments are provided
if [ $# -lt 2 ]; then
    echo "Error: Two arguments are required. Please provide the experiment name and genome name."
    exit 1
fi

experiment="$1"
echo "experiment: $experiment"
genome="$2"
echo "genome: $genome"

# Check if required files exist
[ -f "nanopolish_data_${experiment}/${experiment}.pass.fastq.gz" ] || { echo "Error: ${experiment}.pass.fastq.gz not found"; exit 1; }
[ -f "genomes/$genome/$genome.fa" ] || { echo "Error: genomes/$genome/$genome.fa not found"; exit 1; }

echo "# Copy the .fastq.gz file to $TMPDIR"
cp "nanopolish_data_${experiment}/${experiment}.pass.fastq.gz" "$TMPDIR/"

echo "# Copy the genome file to $TMPDIR"
cp "genomes/$genome/$genome.fa" "$TMPDIR/"

echo "# Running minimap2 for alignment and samtools for sorting"
minimap2 -ax splice -uf -k14 -G 10000 -t 24 \
"${TMPDIR}/$genome.fa" "${TMPDIR}/${experiment}.pass.fastq.gz" \
| samtools sort -@ 24 -o "${TMPDIR}/${experiment}.bam"

echo "# Indexing the BAM file"
samtools index "${TMPDIR}/${experiment}.bam"

mkdir -p "nanopolish_data_${experiment}/minmap_out_$genome"

echo "# Copying the BAM and BAI files back to the current directory"
cp "${TMPDIR}/${experiment}.bam" "nanopolish_data_${experiment}/minmap_out_$genome"
cp "${TMPDIR}/${experiment}.bam.bai" "nanopolish_data_${experiment}/minmap_out_$genome"

echo "minimap2 Script completed successfully"