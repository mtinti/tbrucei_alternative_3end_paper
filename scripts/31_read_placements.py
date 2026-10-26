#!/usr/bin/env python3
"""Collect the placements of every poly(A) read in one pass over the BAMs (docs/methods.md, section 7).

The assignment rule (assignlib.py) sets a read aside when a secondary
alignment of it, anywhere in the genome, has an alignment score (AS) equal to
or higher than its primary alignment's. The secondary records are not near
the primary in a sorted BAM, so this step reads each BAM once and writes one
row per read with a primary alignment: read, bam, mapq (of the primary),
primary_as, best_secondary_as (-1 when the read has no secondary record),
secondaries (their number). Supplementary records are ignored. The output is
gzipped; assignlib.Placements loads it.
"""
import argparse
import gzip
from collections import defaultdict

import pysam


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--bams", nargs="+", required=True, help="label=path, one per BAM")
    ap.add_argument("--out", required=True, help="read_placements.tsv.gz")
    args = ap.parse_args()
    with gzip.open(args.out, "wt") as out:
        out.write("read\tbam\tmapq\tprimary_as\tbest_secondary_as\tsecondaries\n")
        for spec in args.bams:
            label, path = spec.split("=", 1)
            primary, secondary = {}, defaultdict(list)
            with pysam.AlignmentFile(path) as bam:
                for r in bam.fetch(until_eof=True):
                    if r.is_unmapped or r.is_supplementary:
                        continue
                    if r.is_secondary:
                        secondary[r.query_name].append(int(r.get_tag("AS")))
                    else:
                        primary[r.query_name] = (r.mapping_quality, int(r.get_tag("AS")))
            n_eq = 0
            for name, (q, a) in primary.items():
                s = secondary.get(name, [])
                best = max(s) if s else -1
                n_eq += best >= a
                out.write(f"{name}\t{label}\t{q}\t{a}\t{best}\t{len(s)}\n")
            print(f"{label}: {len(primary)} reads with a primary alignment, {n_eq} with a secondary scoring at least as high")
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
