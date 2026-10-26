#!/usr/bin/env python3
"""Count where poly(A) reads end on the genome and write the counts as BED.

For each primary alignment with at least --min-aligned-length aligned bases
and mapping quality at least --min-mapq, the 3' end is the reference end for
forward reads and the reference start for reverse reads. The output has one
line per position and strand: chrom, start, end, 'polyA_site', number of
reads ending there, strand. The original run wrote '.' for the strand and
kept secondary and supplementary alignments; pass --keep-all-alignments and
--no-strand to reproduce it.
"""
import argparse
from collections import defaultdict
import pysam


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--bam", required=True)
    ap.add_argument("--out-bed", required=True)
    ap.add_argument("--min-aligned-length", type=int, default=200)
    ap.add_argument("--min-mapq", type=int, default=0)
    ap.add_argument("--keep-all-alignments", action="store_true",
                    help="also count secondary and supplementary alignments (original behaviour)")
    ap.add_argument("--no-strand", action="store_true", help="write '.' as strand (original behaviour)")
    args = ap.parse_args()

    counts = defaultdict(int)
    total = used = 0
    with pysam.AlignmentFile(args.bam, "rb") as bam:
        for read in bam.fetch():
            total += 1
            if read.is_unmapped or read.query_alignment_length < args.min_aligned_length:
                continue
            if not args.keep_all_alignments and (read.is_secondary or read.is_supplementary):
                continue
            if read.mapping_quality < args.min_mapq:
                continue
            used += 1
            pos = read.reference_start if read.is_reverse else read.reference_end
            strand = "." if args.no_strand else ("-" if read.is_reverse else "+")
            counts[(read.reference_name, pos, strand)] += 1

    with open(args.out_bed, "w") as out:
        for (chrom, pos, strand) in sorted(counts):
            out.write(f"{chrom}\t{pos}\t{pos + 1}\tpolyA_site\t{counts[(chrom, pos, strand)]}\t{strand}\n")
    print(f"{used} of {total} reads had at least {args.min_aligned_length} aligned bases")


if __name__ == "__main__":
    main()
