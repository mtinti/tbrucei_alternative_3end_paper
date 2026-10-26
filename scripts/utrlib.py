"""Shared helpers: the end-site BED format and window coordinates.

Terms used across the scripts:
  end site      A genomic position where a poly(A)-tailed nanopore read ends,
                with the number of reads ending there (the BED score).
  window        A genomic interval, written chrom:start..end:strand, in which
                end sites are collected for one gene.
"""
import pandas as pd

SITE_COLUMNS = ["chromosome", "start", "end", "name", "score", "strand"]


def parse_coord(coord):
    """Split 'chrom:start..end:strand' into its four parts."""
    chrom, span, strand = coord.split(":")
    start, end = span.split("..")
    return chrom, int(start), int(end), strand


def load_site_beds(paths_by_replicate):
    """Read one end-site BED per replicate. Returns {replicate: DataFrame}."""
    out = {}
    for rep, path in paths_by_replicate.items():
        df = pd.read_csv(path, sep="\t", header=None, names=SITE_COLUMNS)
        out[rep] = df
    return out
