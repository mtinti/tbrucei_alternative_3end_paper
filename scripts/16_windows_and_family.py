#!/usr/bin/env python3
"""Define the windows and the gene family for the revised analysis (docs/methods.md, section 1).

No stage label is read here. For every protein-coding gene in the GTF:

CDS: the outermost span of the gene's CDS records, written as a half-open
BED interval (GTF start minus one, GTF end). The transcript start from the
GTF transcript record, where one exists, is kept for the event rules later.

Nanopore window: from the CDS 5' end, through the CDS, to the far end of the
longest eligible full-length read plus --pad bases, in transcript direction.
Ends inside the CDS are counted, as in the original analysis, because some
genes end part of their transcripts there (Tb927.11.1710); their clusters get
a negative distance from the CDS 3' end. --utr-only starts the window at the
CDS 3' end instead. An eligible read is
a primary alignment on the gene's strand that spans the whole CDS and touches
no other CDS. On the minus strand the window runs toward lower coordinates and
its BED start is its biological end.

Intergenic window: from the CDS 5' end, through the CDS, to the nearest CDS
boundary of any gene on either strand downstream in transcript direction.

Depth: the window enters the family when every replicate has at least
--min-ends primary, strand-matched poly(A) read ends inside it.

Outputs in --outdir: cds_v2.bed, transcript_starts.tsv, family_nanopore.tsv
and family_intergenic.tsv (one row per gene with the window, the end count
per replicate, in_family and the exclusion reason where excluded).
"""
import argparse
import bisect
import os
import re
import sys
import numpy as np
import pandas as pd
import pysam

sys.path.insert(0, os.path.dirname(__file__))
import utrlib  # noqa: E402

GTF_COLUMNS = ["chromosome", "source", "annotation_type", "start", "end", "score", "strand", "frame", "attribute"]


def attr(attribute, key):
    m = re.search(key + r' "([^"]+)"', attribute)
    return m.group(1) if m else None


def read_annotation(gtf_path):
    """CDS table (half-open, zero-based) and transcript 5' ends per gene."""
    gtf = pd.read_csv(gtf_path, sep="\t", header=None, names=GTF_COLUMNS, comment="#")
    gtf["gene"] = gtf["attribute"].map(lambda a: attr(a, "gene_id"))
    cds = gtf[gtf["annotation_type"] == "CDS"].groupby("gene").agg(
        chromosome=("chromosome", "first"), strand=("strand", "first"), gtf_start=("start", "min"),
        gtf_end=("end", "max"), cds_records=("start", "size")).reset_index()
    cds["start"] = cds["gtf_start"] - 1          # GTF is one-based closed; BED is zero-based half-open
    cds["end"] = cds["gtf_end"]
    tr = gtf[gtf["annotation_type"] == "transcript"].groupby("gene").agg(
        strand=("strand", "first"), tr_start=("start", "min"), tr_end=("end", "max")).reset_index()
    # the 5' end in transcript direction, zero-based
    tr["transcript_5prime"] = np.where(tr["strand"] == "+", tr["tr_start"] - 1, tr["tr_end"] - 1)
    return cds.sort_values(["chromosome", "start"]).reset_index(drop=True), tr[["gene", "transcript_5prime"]]


class CdsIndex:
    """Per-chromosome sorted CDS intervals for overlap queries."""

    def __init__(self, cds):
        self.by_chrom = {}
        for chrom, grp in cds.groupby("chromosome"):
            g = grp.sort_values("start")
            self.by_chrom[chrom] = (g["start"].to_numpy(), g["end"].to_numpy(), g["strand"].to_numpy(),
                                    g["gene"].to_numpy(), np.maximum.accumulate(g["end"].to_numpy()))

    def overlapping(self, chrom, start, end):
        """Indices of CDS with start < end and end > start."""
        if chrom not in self.by_chrom:
            return []
        starts, ends, _, _, maxend = self.by_chrom[chrom]
        hi = bisect.bisect_left(starts, end)
        # walk back while some earlier CDS could still reach past start
        lo = bisect.bisect_right(maxend, start)
        idx = [i for i in range(lo, hi) if ends[i] > start]
        return idx


def longest_eligible_reads(bam_paths, cds, index, chrom_len):
    """For each gene, the eligible full-length read reaching furthest 3'."""
    best = {}   # gene -> (far_end, read_name, read_start, read_end)
    stats = dict(reads=0, no_cds=0, several_cds=0, strand_mismatch=0, not_spanning=0, eligible=0)
    gene_row = {g: i for i, g in enumerate(cds["gene"])}
    for path in bam_paths:
        with pysam.AlignmentFile(path) as bam:
            for r in bam.fetch():
                if r.is_unmapped or r.is_secondary or r.is_supplementary or r.reference_name not in chrom_len:
                    continue
                stats["reads"] += 1
                idx = index.overlapping(r.reference_name, r.reference_start, r.reference_end)
                if not idx:
                    stats["no_cds"] += 1
                    continue
                if len(idx) > 1:
                    stats["several_cds"] += 1
                    continue
                starts, ends, strands, genes, _ = index.by_chrom[r.reference_name]
                i = idx[0]
                strand = "-" if r.is_reverse else "+"
                if strand != strands[i]:
                    stats["strand_mismatch"] += 1
                    continue
                if not (r.reference_start <= starts[i] and r.reference_end >= ends[i]):
                    stats["not_spanning"] += 1
                    continue
                stats["eligible"] += 1
                far = r.reference_end if strand == "+" else r.reference_start
                g = genes[i]
                cur = best.get(g)
                if cur is None or (far > cur[0] if strand == "+" else far < cur[0]):
                    best[g] = (far, r.query_name, r.reference_start, r.reference_end)
    return best, stats


class EndCounter:
    """Strand-matched poly(A) end counts in any interval, per replicate, from the end-site BEDs."""

    def __init__(self, site_dfs):
        self.cum = {}
        for rep, df in site_dfs.items():
            for (chrom, strand), grp in df.groupby(["chromosome", "strand"]):
                g = grp.sort_values("start")
                pos = g["start"].to_numpy()
                cs = np.concatenate([[0], np.cumsum(g["score"].to_numpy())])
                self.cum[(rep, chrom, strand)] = (pos, cs)

    def count(self, rep, chrom, strand, start, end):
        key = (rep, chrom, strand)
        if key not in self.cum or end <= start:
            return 0
        pos, cs = self.cum[key]
        lo, hi = np.searchsorted(pos, start, "left"), np.searchsorted(pos, end, "left")
        return int(cs[hi] - cs[lo])


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--fai", required=True, help="FASTA index, for chromosome lengths")
    ap.add_argument("--full-length-bams", nargs="+", required=True, help="spliced leader plus poly(A) reads, any number of BAMs")
    ap.add_argument("--sites", nargs="+", required=True, help="replicate=path end-site BED files (stranded, primary only)")
    ap.add_argument("--pad", type=int, default=100)
    ap.add_argument("--min-ends", type=int, default=10)
    ap.add_argument("--utr-only", action="store_true", help="start the window at the CDS 3' end, not the 5' end")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    chrom_len = {l.split("\t")[0]: int(l.split("\t")[1]) for l in open(args.fai)}
    cds, tr = read_annotation(args.gtf)
    cds = cds[cds["chromosome"].isin(chrom_len)].reset_index(drop=True)
    cds[["chromosome", "start", "end", "gene", "cds_records", "strand"]].to_csv(
        os.path.join(args.outdir, "cds_v2.bed"), sep="\t", header=False, index=False)
    tr.to_csv(os.path.join(args.outdir, "transcript_starts.tsv"), sep="\t", index=False)
    print(f"{len(cds)} genes with a CDS, {int((cds['cds_records'] > 1).sum())} with several CDS records (outermost span used)")

    index = CdsIndex(cds)
    best, stats = longest_eligible_reads(args.full_length_bams, cds, index, chrom_len)
    print("full-length reads:", stats)
    site_dfs = utrlib.load_site_beds(dict(s.split("=", 1) for s in args.sites))
    reps = list(site_dfs)
    counter = EndCounter(site_dfs)

    # nearest CDS boundary in transcript direction, any strand, for the intergenic window
    def next_boundary(chrom, strand, cds_start, cds_end):
        starts, ends, _, _, _ = index.by_chrom[chrom]
        if strand == "+":
            i = bisect.bisect_right(starts, cds_end)
            return int(starts[i]) if i < len(starts) else None
        cands = ends[ends <= cds_start]
        return int(cands.max()) if len(cands) else None

    nano, inter = [], []
    for c in cds.itertuples():
        base = dict(gene=c.gene, chromosome=c.chromosome, strand=c.strand, cds_start=c.start, cds_end=c.end)
        # nanopore window
        row = dict(base, read_name="", read_start="", read_end="", window_start="", window_end="", exclusion_reason="")
        b = best.get(c.gene)
        if b is None:
            row["exclusion_reason"] = "no eligible full-length read (same strand, spanning the whole CDS, touching no other CDS)"
        else:
            far, name, rs, re_ = b
            row.update(read_name=name, read_start=rs, read_end=re_)
            if c.strand == "+":
                ws, we = (c.end if args.utr_only else c.start), min(far + args.pad, chrom_len[c.chromosome])
            else:
                ws, we = max(far - args.pad, 0), (c.start if args.utr_only else c.end)
            row.update(window_start=ws, window_end=we)
            if (far <= c.end if c.strand == "+" else far >= c.start):
                row["exclusion_reason"] = "invalid window (read 3' end at or before the CDS end)"
        nano.append(row)
        # intergenic window
        row = dict(base, window_start="", window_end="", exclusion_reason="")
        nb = next_boundary(c.chromosome, c.strand, c.start, c.end)
        if nb is None:
            row["exclusion_reason"] = "no CDS downstream on the chromosome"
        else:
            ws, we = (c.start, nb) if c.strand == "+" else (nb, c.end)
            if args.utr_only:
                ws, we = (c.end, nb) if c.strand == "+" else (nb, c.start)
            row.update(window_start=ws, window_end=we)
            if (nb <= c.end if c.strand == "+" else nb >= c.start):
                row["exclusion_reason"] = "invalid window (overlapping or abutting CDS)"
        inter.append(row)

    for name, rows in (("nanopore", nano), ("intergenic", inter)):
        df = pd.DataFrame(rows)
        for rep in reps:
            df[f"ends_{rep}"] = [counter.count(rep, r.chromosome, r.strand, r.window_start, r.window_end)
                                 if r.exclusion_reason == "" else 0 for r in df.itertuples()]
        depth = df[[f"ends_{r}" for r in reps]].min(axis=1)
        shallow = (df["exclusion_reason"] == "") & (depth < args.min_ends)
        df.loc[shallow, "exclusion_reason"] = f"fewer than {args.min_ends} poly(A) read ends in at least one replicate"
        df["min_ends_per_replicate"] = depth
        df["in_family"] = df["exclusion_reason"] == ""
        df.to_csv(os.path.join(args.outdir, f"family_{name}.tsv"), sep="\t", index=False)
        print(f"{name}: {int(df['in_family'].sum())} windows in the family of {len(df)} genes")
        print(df.loc[~df["in_family"], "exclusion_reason"].value_counts().to_string())


if __name__ == "__main__":
    main()
