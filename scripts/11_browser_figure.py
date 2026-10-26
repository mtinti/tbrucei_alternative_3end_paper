#!/usr/bin/env python3
"""Draw a genome-browser style figure for one gene from the BAM files.

Panels, top to bottom: gene models; read coverage per stage; a sample of
individual reads per stage with their 3' ends marked; the 3' end density of
each replicate; and the sites tested by DRIMSeq, shaded through every panel
with their fitted usage in each stage. With --transcripts, the gene panel
also draws the independent transcripts from step 13 (full-length reads that
overlap no CDS, the candidate non-coding RNAs) on the gene's strand, and with
--overrides the transcript recorded for a manual class B call. With --ncrna,
the annotated non-coding RNA genes (Ensembl BED with biotype in column 5)
are drawn in green on both strands.

Inputs: the per-stage poly(A) BAMs (merged replicates), the per-replicate
end-site BED files, the CDS BED, and the step 08/09 folder for the sites and
their usage.
"""
import argparse
import os
import re
import sys
import numpy as np
import pandas as pd
import pysam
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyArrow  # noqa: E402

sys.path.insert(0, os.path.dirname(__file__))
import utrlib  # noqa: E402

COL = {"BSF": "#1f77b4", "PCF": "#ff7f0e"}


def primary(read):
    return not (read.is_unmapped or read.is_secondary or read.is_supplementary)


def coverage(bam, chrom, start, end):
    cov = np.zeros(end - start)
    for r in bam.fetch(chrom, start, end):
        if not primary(r):
            continue
        for b_start, b_end in r.get_blocks():
            lo, hi = max(b_start, start), min(b_end, end)
            if hi > lo:
                cov[lo - start:hi - start] += 1
    return cov


def sample_reads(bam, chrom, start, end, n, seed=42):
    reads = [(r.reference_start, r.reference_end, r.is_reverse) for r in bam.fetch(chrom, start, end) if primary(r)]
    total = len(reads)
    rng = np.random.default_rng(seed)
    if total > n:
        reads = [reads[i] for i in sorted(rng.choice(total, n, replace=False))]
    # pack into rows so reads do not overlap
    rows, ends = [], []
    for s, e, rev in sorted(reads):
        for i, last in enumerate(ends):
            if s > last + 20:
                ends[i] = e
                rows.append((s, e, rev, i))
                break
        else:
            ends.append(e)
            rows.append((s, e, rev, len(ends) - 1))
    return rows, total


def independent_transcripts(args, chrom, strand, lo, hi):
    """Step 13 transcripts on the gene's strand inside the region, plus the override transcript.

    Returns a DataFrame with start, end and a label giving the read count per stage.
    """
    rows = []
    if args.transcripts:
        t = pd.read_csv(args.transcripts, sep="\t", header=None,
                        names=["chrom", "start", "end", "name", "reads", "strand", "reads_BSF", "reads_PCF"])
        t = t[(t["chrom"] == chrom) & (t["strand"] == strand) & (t["end"] > lo) & (t["start"] < hi)
              & (t["reads"] >= args.min_transcript_reads)]
        for _, r in t.iterrows():
            rows.append((r["start"], r["end"], f"independent transcript: {r['reads_BSF']} BSF, {r['reads_PCF']} PCF "
                                               f"full-length reads"))
    if args.catalogue:
        c = pd.read_csv(args.catalogue, sep="\t")
        reps = [x for x in c.columns if x.startswith(("bs_", "pc_"))]
        c = c[(c["chrom"] == chrom) & (c["strand"] == strand) & (c["end"] > lo) & (c["start"] < hi) & (c[reps].sum(axis=1) >= args.min_transcript_reads)]
        for _, r in c.iterrows():
            nb = int(r[[x for x in reps if x.startswith("bs_")]].sum()); npc = int(r[[x for x in reps if x.startswith("pc_")]].sum())
            rows.append((r["start"], r["end"], f"catalogue {r['kind']}: {nb} BSF, {npc} PCF full-length reads"))
    if args.overrides:
        ov = pd.read_csv(args.overrides, sep="\t")
        ov = ov[(ov["gene"] == args.gene) & ov["independent_transcript"].notna()]
        for coord in ov["independent_transcript"]:
            m = re.match(r"(.+):(\d+)-(\d+)\(([+-])\)", coord)
            if not m or m.group(1) != chrom:
                continue
            start, end = int(m.group(2)), int(m.group(3))
            if not any(abs(start - a) < 50 and abs(end - b) < 50 for a, b, _ in rows):
                rows.append((start, end, "independent transcript recorded in class_overrides.tsv"))
    return pd.DataFrame(rows, columns=["start", "end", "label"])


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gene", required=True)
    ap.add_argument("--bsf-bam", required=True)
    ap.add_argument("--pcf-bam", required=True)
    ap.add_argument("--cds-bed", required=True)
    ap.add_argument("--drimseq", required=True, help="step 08/09 folder with sites.tsv and drimseq_site.tsv")
    ap.add_argument("--sites", nargs="+", required=True, help="replicate=path end-site BED files")
    ap.add_argument("--bsf", nargs="+", default=["bs_rep_1", "bs_rep_2", "bs_rep_3"])
    ap.add_argument("--pcf", nargs="+", default=["pc_rep_1", "pc_rep_2", "pc_rep_3"])
    ap.add_argument("--flank", type=int, default=300)
    ap.add_argument("--reads-per-stage", type=int, default=60)
    ap.add_argument("--transcripts", help="independent_transcripts.bed from step 13 (earlier analysis)")
    ap.add_argument("--catalogue", help="catalogue.tsv from step 22; its units on the gene strand are drawn with their reads per stage")
    ap.add_argument("--overrides", help="class_overrides.tsv; its transcript for this gene is drawn too")
    ap.add_argument("--ncrna", help="BED of annotated non-coding genes: chrom, start, end, id, biotype, strand")
    ap.add_argument("--min-transcript-reads", type=int, default=5,
                    help="draw a step 13 transcript only when this many full-length reads support it")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    cds = pd.read_csv(args.cds_bed, sep="\t", names=["chrom", "start", "end", "gene", "score", "strand"])
    sites = pd.read_csv(os.path.join(args.drimseq, "sites.tsv"), sep="\t")
    usage_path = os.path.join(args.drimseq, "dexseq_site.tsv")
    if not os.path.exists(usage_path):
        usage_path = os.path.join(args.drimseq, "drimseq_site.tsv")
    if os.path.exists(usage_path):
        usage = pd.read_csv(usage_path, sep="\t")
        gsites = sites[sites["gene_id"] == args.gene].merge(usage, on=["gene_id", "feature_id"], how="left")
    else:  # clusters only, no test yet (the eye check of the cluster step)
        gsites = sites[sites["gene_id"] == args.gene].copy()
    if gsites.empty and os.path.exists(usage_path):
        sys.exit(f"{args.gene} has no sites in {args.drimseq}")
    if gsites.empty:  # untestable window: draw the clusters from clusters.bed instead
        cb = pd.read_csv(os.path.join(args.drimseq, "clusters.bed"), sep="\t", header=None,
                         names=["chromosome", "genomic_start", "genomic_end", "name", "pooled", "strand", "support"])
        cb = cb[cb["name"].str.startswith(args.gene + "|")]
        gsites = pd.DataFrame(dict(gene_id=args.gene, feature_id=cb["name"].str.split("|").str[1] + " " + cb["support"],
                                   chromosome=cb["chromosome"], genomic_start=cb["genomic_start"],
                                   genomic_end=cb["genomic_end"], strand=cb["strand"]))
    g = cds[cds["gene"] == args.gene].iloc[0]
    chrom, strand = g["chrom"], g["strand"]
    lo = min(g["start"], gsites["genomic_start"].min()) - args.flank
    hi = max(g["end"], gsites["genomic_end"].max()) + args.flank
    transcripts = independent_transcripts(args, chrom, strand, lo, hi)
    if not transcripts.empty:
        lo = min(lo, transcripts["start"].min() - args.flank)
        hi = max(hi, transcripts["end"].max() + args.flank)
    x = np.arange(lo, hi)

    bams = {"BSF": pysam.AlignmentFile(args.bsf_bam), "PCF": pysam.AlignmentFile(args.pcf_bam)}
    cov = {k: coverage(b, chrom, lo, hi) for k, b in bams.items()}
    reads = {k: sample_reads(b, chrom, lo, hi, args.reads_per_stage) for k, b in bams.items()}
    site_dfs = utrlib.load_site_beds(dict(s.split("=", 1) for s in args.sites))
    reps = args.bsf + args.pcf
    stage_of = {r: "BSF" for r in args.bsf} | {r: "PCF" for r in args.pcf}

    heights = [1.3, 1.4, 2.2, 2.2] + [0.5] * len(reps)
    fig, axes = plt.subplots(len(heights), 1, figsize=(14, sum(heights) + 0.8), sharex=True,
                             gridspec_kw={"height_ratios": heights, "hspace": 0.15, "top": 0.95})
    ax_gene, ax_cov, ax_bsf, ax_pcf = axes[:4]
    ax_reps = axes[4:]

    # gene models
    near = cds[(cds["chrom"] == chrom) & (cds["end"] > lo) & (cds["start"] < hi)]
    for _, c in near.iterrows():
        colour = "black" if c["gene"] == args.gene else "grey"
        s, e = max(c["start"], lo), min(c["end"], hi)
        dx = (e - s) * (1 if c["strand"] == "+" else -1)
        ax_gene.add_patch(FancyArrow(s if c["strand"] == "+" else e, 0.5, dx, 0, width=0.5,
                                     head_width=0.9, head_length=min(60, abs(dx) / 3), length_includes_head=True,
                                     color=colour))
        ax_gene.text((s + e) / 2, 1.05, c["gene"], ha="center", va="bottom", fontsize=8, color=colour)
    # independent transcripts (candidate non-coding RNAs) above the CDS models
    top = 1.6
    for _, t in transcripts.iterrows():
        dx = (t["end"] - t["start"]) * (1 if strand == "+" else -1)
        ax_gene.add_patch(FancyArrow(t["start"] if strand == "+" else t["end"], 2.0, dx, 0, width=0.4,
                                     head_width=0.8, head_length=min(60, abs(dx) / 3), length_includes_head=True,
                                     color="purple", alpha=0.7))
        ax_gene.text((t["start"] + t["end"]) / 2, 2.5, t["label"], ha="center", va="bottom", fontsize=7,
                     color="purple")
        top = 3.4
    if args.ncrna:
        nc = pd.read_csv(args.ncrna, sep="\t", header=None, usecols=range(6),
                         names=["chrom", "start", "end", "id", "biotype", "strand"])
        nc = nc[(nc["chrom"] == chrom) & (nc["end"] > lo) & (nc["start"] < hi) & (nc["biotype"] != "pseudogene")]
        for _, t in nc.iterrows():
            s_, e_ = max(t["start"], lo), min(t["end"], hi)
            dx = (e_ - s_) * (1 if t["strand"] == "+" else -1)
            ax_gene.add_patch(FancyArrow(s_ if t["strand"] == "+" else e_, 1.4, dx, 0, width=0.4,
                                         head_width=0.8, head_length=min(60, abs(dx) / 3), length_includes_head=True,
                                         color="seagreen", alpha=0.8))
            ax_gene.text((s_ + e_) / 2, 1.0, f"{t['id']} ({t['biotype']}, annotated)", ha="center", va="top",
                         fontsize=7, color="seagreen")
            top = max(top, 2.0)
    ax_gene.set_ylim(0, top)
    ax_gene.set_yticks([])
    ax_gene.set_ylabel("CDS and\nother\ntranscripts" if top > 1.6 else "CDS",
                       rotation=0, ha="right", va="center")

    # coverage
    for k in ("BSF", "PCF"):
        ax_cov.fill_between(x, cov[k], step="mid", alpha=0.5, color=COL[k], label=f"{k} poly(A) reads")
    ax_cov.set_ylabel("coverage\n(reads)")
    ax_cov.legend(loc="upper left", fontsize=8, frameon=False)

    # reads
    for ax, k in ((ax_bsf, "BSF"), (ax_pcf, "PCF")):
        rows, total = reads[k]
        for s, e, rev, row in rows:
            ax.plot([s, e], [row, row], color=COL[k], lw=1.2, alpha=0.8)
            end3 = s if rev else e
            ax.plot(end3, row, marker="|", color="black", ms=6, mew=1.2)
        ax.set_ylim(-1, max(r[3] for r in rows) + 1 if rows else 1)
        ax.set_yticks([])
        ax.set_ylabel(f"{k} reads\n({len(rows)} of {total} shown,\n3' end marked)", fontsize=8)

    # 3' end density per replicate, as fraction of the replicate's ends in the region
    for ax, r in zip(ax_reps, reps):
        d = site_dfs[r]
        d = d[(d["chromosome"] == chrom) & (d["start"] >= lo) & (d["start"] < hi) & (d["strand"].isin([strand, "."]))]
        dens = np.zeros(hi - lo)
        for p, sc in zip(d["start"], d["score"]):
            dens[p - lo] += sc
        total = dens.sum()
        ax.fill_between(x, dens / total if total else dens, step="mid", color=COL[stage_of[r]], alpha=0.9)
        ax.set_yticks([])
        ax.set_ylabel(f"{r}\n{int(total)} ends", rotation=0, ha="right", va="center", fontsize=8)

    # sites shaded through all panels, labels staggered so neighbours do not collide
    gsites = gsites.sort_values("genomic_start").reset_index(drop=True)
    for i, (_, s) in enumerate(gsites.iterrows()):
        label = s["feature_id"].replace("_", " ")
        if "usage_BSF" in s and not pd.isna(s.get("usage_BSF")):
            label += f": BSF {s['usage_BSF']:.2f}, PCF {s['usage_PCF']:.2f}"
        for ax in axes:
            ax.axvspan(s["genomic_start"], s["genomic_end"], color="gold", alpha=0.25, lw=0)
        y = -0.25 - 0.55 * (i % 3)
        centre = (s["genomic_start"] + s["genomic_end"]) / 2
        ax_gene.annotate(label, xy=(centre, 0.0), xytext=(centre, y), fontsize=7, ha="center", va="top",
                         arrowprops=dict(arrowstyle="-", color="goldenrod", lw=0.8))
    ax_gene.set_ylim(-2.0, top)

    # every figure reads 5' to 3' from left to right
    axes[-1].set_xlabel(f"{chrom} position (bp), 5' to 3'" + (" (minus strand, axis reversed)" if strand == "-" else ""))
    axes[-1].set_xlim((hi, lo) if strand == "-" else (lo, hi))
    fig.suptitle(f"{args.gene} ({strand} strand): 3' ends of poly(A) nanopore reads, "
                 f"bloodstream (BSF) and procyclic (PCF) forms. Shaded: tested end clusters, with their usage per stage.")
    fig.savefig(args.out, dpi=110, bbox_inches="tight")
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
