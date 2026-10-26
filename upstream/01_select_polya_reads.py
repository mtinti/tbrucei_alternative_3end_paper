#!/usr/bin/env python3
"""Keep the reads with a measured poly(A) tail and write them to a new BAM.

Input: the nanopolish polya table and the minimap2 BAM of one replicate.
A read is kept when nanopolish gives qc_tag PASS and a tail longer than
--min-tail nucleotides. The kept read names are also written as a text file.
"""
import argparse
import pandas as pd
import pysam


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--polya-tsv", required=True, help="nanopolish polya output (TSV)")
    ap.add_argument("--bam", required=True, help="minimap2 BAM of the same replicate")
    ap.add_argument("--out-bam", required=True)
    ap.add_argument("--out-ids", required=True, help="text file with one kept read name per line")
    ap.add_argument("--min-tail", type=float, default=25)
    args = ap.parse_args()

    polya = pd.read_csv(args.polya_tsv, sep="\t")
    polya = polya[(polya["qc_tag"] == "PASS") & (polya["polya_length"] > args.min_tail)].dropna()
    keep = set(polya["readname"])
    pd.Series(sorted(keep)).to_csv(args.out_ids, index=False, header=False)

    written = 0
    with pysam.AlignmentFile(args.bam, "rb") as src, \
            pysam.AlignmentFile(args.out_bam, "wb", template=src) as dst:
        for read in src:
            if read.query_name in keep:
                dst.write(read)
                written += 1
    pysam.index(args.out_bam)
    print(f"{written} reads kept out of {len(keep)} poly(A) read names")


if __name__ == "__main__":
    main()
