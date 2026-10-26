#!/bin/bash
#$ -adds l_hard local_free 500G
#$ -mods l_hard m_mem_free 20G
#$ -adds l_hard avx 1 
#$ -cwd
#$ -V
#$ -j y
#$ -pe smp 24
#$ -o errors_$JOB_ID


# Exit immediately if a command exits with a non-zero status
set -e

# Check if TMPDIR is set
if [ -z "$TMPDIR" ]; then
    echo "Error: TMPDIR is not set"
    exit 1
fi

# Check if INPUT_FILES and OUTPUT_FILE are set
if [ -z "$INPUT_FILES" ] || [ -z "$OUTPUT_FILE" ]; then
    echo "Error: INPUT_FILES and OUTPUT_FILE must be set"
    echo "Usage: qsub -v INPUT_FILES='file1.bam file2.bam ...',OUTPUT_FILE='output.bam' $0"
    exit 1
fi

# Convert INPUT_FILES string to an array
read -ra input_files <<< "$INPUT_FILES"

output_file="$OUTPUT_FILE"

echo "Input files: ${input_files[@]}"
echo "Output file: $output_file"

# Copy input files to TMPDIR
echo "Copying input files to $TMPDIR"
for file in "${input_files[@]}"; do
    cp "$file" "$TMPDIR/"
    echo "Copied $file"
done

# Prepare input file list for samtools merge
input_files_tmp=()
for file in "${input_files[@]}"; do
    input_files_tmp+=("$TMPDIR/$(basename "$file")")
done

# Merge BAM files
echo "Merging BAM files"
samtools merge -@ 8 "$TMPDIR/$(basename "$output_file")" "${input_files_tmp[@]}"

# Index the merged BAM file
echo "Indexing merged BAM file"
samtools index "$TMPDIR/$(basename "$output_file")"

# Create output directory if it doesn't exist
output_dir=$(dirname "$output_file")
mkdir -p "$output_dir"

# Copy the merged BAM and BAI files to the output location
echo "Copying merged BAM and BAI files to output location"
cp "$TMPDIR/$(basename "$output_file")" "$output_file"
cp "$TMPDIR/$(basename "$output_file").bai" "${output_file}.bai"

echo "Script completed successfully"