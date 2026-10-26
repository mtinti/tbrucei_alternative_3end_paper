#!/bin/bash
#$ -adds l_hard local_free 500G
#$ -mods l_hard m_mem_free 20G
#$ -adds l_hard avx 1 
#$ -cwd
#$ -V
#$ -j y
#$ -mods l_hard gpu 1
#$ -adds l_hard cuda.0.name 'NVIDIA A40'
#$ -pe smp 4
#$ -o errors_guppy_$JOB_ID


# Exit immediately if a command exits with a non-zero status
set -e

# Check if TMPDIR is set
if [ -z "$TMPDIR" ]; then
    echo "Error: TMPDIR is not set"
    exit 1
fi


# Check if an argument is provided
if [ $# -eq 0 ]; then
    echo "Error: No argument provided. Please provide the experiment name."
    exit 1
fi

experiment=$1

echo 'experiment: '$experiment

echo 'CUDA_VISIBLE_DEVICES :'$CUDA_VISIBLE_DEVICES

echo '# Copy the tar.gz file to $TMPDIR'
cp "${experiment}.tar.gz" "$TMPDIR/"

echo '# Create the fast5_in directory'
mkdir -p "$TMPDIR/fast5_in"

echo '# Extract the contents'
tar -xzf "$TMPDIR/${experiment}.tar.gz" --transform 's:.*/::' -C "$TMPDIR/fast5_in"

echo '# Run guppy_basecaller'
##/cluster/majf_lab/mtinti/nanpore_data/ont-guppy/bin/guppy_basecaller -i "$TMPDIR/fast5_in" \
##--trim_strategy 'none' \


/cluster/majf_lab/mtinti/nanopore_software/ont-guppy-621/bin/guppy_basecaller \
-i "$TMPDIR/fast5_in" \
-s "$TMPDIR/fastq" \
--fast5_out \
--flowcell FLO-MIN106 \
--kit SQK-RNA002 \
-x "auto" \
--recursive \
--trim_strategy 'none' \
--disable_pings 1 \
--chunk_size 2000 --chunks_per_runner 550 --gpu_runners_per_device 12

ls -l $TMPDIR
ls -l "$TMPDIR/fastq/"


echo '# Run pycoQC'
pycoQC -f "$TMPDIR/fastq/sequencing_summary.txt" \
-o "$TMPDIR/pycoQC_${experiment}.html"

echo '# Concatenate and compress fastq files'
cat "$TMPDIR/fastq/pass/"*.fastq | gzip > "$TMPDIR/${experiment}.pass.fastq.gz"


out_dir="nanopolish_data_${experiment}/"
mkdir -p $out_dir


echo "# Copy results back to the current directory ${out_dir}"
cp "$TMPDIR/${experiment}.pass.fastq.gz"  $out_dir
#cp "$TMPDIR/pycoQC_${experiment}.html"  $out_dir

#cp "$TMPDIR/fastq/sequencing_summary.txt" "${out_dir}sequencing_summary_${experiment}.txt"

mkdir -p "${out_dir}guppy_out"
cp -r "$TMPDIR"/fastq/* "${out_dir}guppy_out/"
echo "# done"








