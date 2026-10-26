#!/bin/bash
#$ -adds l_hard local_free 500G
#$ -mods l_hard m_mem_free 20G
#$ -adds l_hard avx 1 
#$ -cwd
#$ -V
#$ -j y
#$ -pe smp 50
#$ -o errors_nanopolish_$JOB_ID

set -e  # Exit immediately if a command exits with a non-zero status

run_polya() {
    if [ $# -ne 2 ]; then
        echo "Error: Two arguments required. Usage: $0 <extension> <genome>"
        exit 1
    fi

    extension="$1"
    genome="$2"
    orig_dir="$(pwd)"
    
    echo "extension: $extension"
    echo "genome: $genome"    
    
    # Create temporary directory
    tmp_dir="$TMPDIR/nanopolish_data_$extension/"
    mkdir -p "$tmp_dir"
    
    echo "Copying input files to temporary directory"
    echo "Copying guppy_out"
    cp -r "nanopolish_data_$extension/guppy_out" "$tmp_dir/guppy_out"
    echo "Copying fastq"
    cp "nanopolish_data_$extension/$extension.pass.fastq.gz" "$tmp_dir/"
    echo "Copying bam"    
    cp "nanopolish_data_$extension/minmap_out_${genome}/${extension}.bam" "$tmp_dir/"
    cp "nanopolish_data_$extension/minmap_out_${genome}/${extension}.bam.bai" "$tmp_dir/"
    echo "Copying genomes"      
    cp "genomes/$genome/$genome.fa" "$tmp_dir/"
    
    echo "Deflating fastq"
    gzip -d "$tmp_dir/$extension.pass.fastq.gz"
    
    # Change to temporary directory
    cd "$tmp_dir" || exit
    
    echo "Running nanopolish index"
    nanopolish index \
    --directory="$tmp_dir/guppy_out/workspace/" \
    --sequencing-summary="$tmp_dir/guppy_out/sequencing_summary.txt" \
    "$tmp_dir/$extension.pass.fastq"
    
    echo "Running nanopolish polya"
    nanopolish polya --threads=50 \
    --reads="${tmp_dir}/${extension}.pass.fastq" \
    --bam="${tmp_dir}/${extension}.bam" \
    --genome="$tmp_dir/$genome.fa" > "$tmp_dir/nanopolish_${genome}_polya_results.tsv"
    
    echo "Copying results back to original directory"
    cp "$tmp_dir/nanopolish_${genome}_polya_results.tsv" "$orig_dir/nanopolish_data_$extension/"
    
    
    
    # Change back to original directory
    cd "$orig_dir" || exit
    
    # Uncomment the following line if you want to clean up the temporary directory
    # rm -rf "$tmp_dir"
}

run_polya "$1" "$2"
echo "# done"