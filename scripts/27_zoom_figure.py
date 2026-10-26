#!/usr/bin/env python3
"""Draw the tested region of one gene with the tested reads only (docs/methods.md, section 6).

The frame runs from 100 bases before the first tested cluster to 100 bases
after the last one, with the CDS 3' end marked. Reads are selected with the
step 02 rules used for the counts (primary, same strand, at least 200 aligned
bases) and only reads whose 3' end lies inside the gene's window are drawn,
so every read in the figure is one the test counted. Panels: clusters with
their usage per stage; reads per stage sorted by 3' end, with the window
total in the label; end counts per replicate with the cluster boundaries.
Every figure reads 5' to 3' from left to right: for minus-strand genes the
x axis is reversed, so the CDS 3' end is on the left and the poly(A) ends on
the right in every figure. Reads of uncertain origin under the Part 7 rule
(assignlib.py: fewer than --min-cds-overlap aligned bases in the CDS, mapping
quality below --min-mapq, a secondary alignment scoring at least as high as the
primary (--placements, step 31), or another CDS covered) are drawn grey; with --assigned (a step 28 folder with its step 18 results) the
labels also give the usage among the gene's assigned reads.
"""
import argparse
import os
import numpy as np
import pandas as pd
import pysam
from assignlib import CdsIndex, Placements, assignment
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

COL = {"BSF": "#1f77b4", "PCF": "#ff7f0e"}


def tested_reads(bam, chrom, strand, win_start, win_end, min_aligned=200):
    out = []
    for r in bam.fetch(chrom, max(0, win_start - 1), win_end + 1):
        if r.is_unmapped or r.is_secondary or r.is_supplementary or r.query_alignment_length < min_aligned:
            continue
        if r.is_reverse != (strand == "-"):
            continue
        end3 = r.reference_start if strand == "-" else r.reference_end
        if win_start <= end3 < win_end:
            out.append((r.reference_start, r.reference_end, end3, r.query_name, r))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gene", required=True)
    ap.add_argument("--family", required=True, help="family_<mode>.tsv from step 16")
    ap.add_argument("--dexseq", required=True, help="step 17/18 folder (sites.tsv, dexseq_site.tsv, site_counts.tsv, samples.tsv)")
    ap.add_argument("--events", help="events_<mode>.tsv from step 21, for the labels")
    ap.add_argument("--bsf-bam", required=True)
    ap.add_argument("--pcf-bam", required=True)
    ap.add_argument("--read-ids", nargs="*", default=[], help="replicate=path, to split the end counts per replicate")
    ap.add_argument("--margin", type=int, default=100)
    ap.add_argument("--max-reads", type=int, default=150, help="reads drawn per stage (all are counted)")
    ap.add_argument("--min-cds-overlap", type=int, default=50, help="aligned bases inside the CDS a read needs to count as the gene's (Part 7)")
    ap.add_argument("--min-mapq", type=int, default=20)
    ap.add_argument("--placements", help="read_placements.tsv.gz from step 31")
    ap.add_argument("--cds-bed", help="cds_v2.bed from step 16; with it the other-gene check of Part 7 applies to the colouring")
    ap.add_argument("--assigned", help="step 28 folder with dexseq_site.tsv, to show the gene-assigned usage in the labels")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    fam = pd.read_csv(args.family, sep="\t").set_index("gene").loc[args.gene]
    chrom, strand = fam["chromosome"], fam["strand"]
    ws, we = int(fam["window_start"]), int(fam["window_end"])
    cds3 = int(fam["cds_end"]) if strand == "+" else int(fam["cds_start"])
    cs, ce = int(fam["cds_start"]), int(fam["cds_end"])
    index = CdsIndex(args.cds_bed) if args.cds_bed else None
    placements = Placements(args.placements) if args.placements else None
    assigned_usage = {}
    if args.assigned and os.path.exists(os.path.join(args.assigned, "dexseq_site.tsv")):
        au = pd.read_csv(os.path.join(args.assigned, "dexseq_site.tsv"), sep="\t")
        au = au[au["gene_id"] == args.gene].set_index("feature_id")
        assigned_usage = {k: (v["usage_BSF"], v["usage_PCF"], bool(v["site_confirmed"])) for k, v in au.iterrows()}
    sites = pd.read_csv(os.path.join(args.dexseq, "sites.tsv"), sep="\t")
    sites = sites[sites["gene_id"] == args.gene].sort_values("genomic_start")
    usage = pd.read_csv(os.path.join(args.dexseq, "dexseq_site.tsv"), sep="\t")
    usage = usage[usage["gene_id"] == args.gene].set_index("feature_id")
    counts = pd.read_csv(os.path.join(args.dexseq, "site_counts.tsv"), sep="\t")
    counts = counts[counts["gene_id"] == args.gene].set_index("feature_id")
    samples = pd.read_csv(os.path.join(args.dexseq, "samples.tsv"), sep="\t")
    reps = samples["sample_id"].tolist(); stage_of = dict(zip(samples["sample_id"], samples["stage"]))
    labels = {}
    if args.events and os.path.exists(args.events):
        ev = pd.read_csv(args.events, sep="\t")
        labels = dict(zip(ev.loc[ev["gene_id"] == args.gene, "feature_id"], ev.loc[ev["gene_id"] == args.gene, "label"]))
    rep_of = {}
    for spec in args.read_ids:
        rep, path = spec.split("=", 1)
        for line in open(path):
            rep_of[line.strip()] = rep

    lo = max(ws, int(sites["genomic_start"].min()) - args.margin)
    hi = min(we, int(sites["genomic_end"].max()) + args.margin)
    if strand == "+":
        lo = min(lo, cds3 - 20)
    else:
        hi = max(hi, cds3 + 20)
    bams = {"BSF": pysam.AlignmentFile(args.bsf_bam), "PCF": pysam.AlignmentFile(args.pcf_bam)}
    reads = {k: tested_reads(b, chrom, strand, ws, we) for k, b in bams.items()}

    heights = [2.2, 2.4, 2.4] + [0.6] * len(reps)
    fig, axes = plt.subplots(len(heights), 1, figsize=(12, sum(heights) + 1), sharex=True,
                             gridspec_kw={"height_ratios": heights, "hspace": 0.12, "top": 0.93})
    ax_model, ax_bsf, ax_pcf = axes[:3]; ax_reps = axes[3:]
    # clusters and labels
    sites = sites.sort_values("genomic_start", ascending=(strand == "+")).reset_index(drop=True)
    for i, s in enumerate(sites.itertuples()):
        for ax in axes:
            ax.axvspan(s.genomic_start, s.genomic_end, color="gold", alpha=0.25, lw=0)
            ax.axvline(s.genomic_start, color="goldenrod", lw=0.6); ax.axvline(s.genomic_end, color="goldenrod", lw=0.6)
        u = usage.loc[s.feature_id] if s.feature_id in usage.index else None
        text = s.feature_id.replace("_", " ")
        if u is not None:
            text += f"\nBSF {u['usage_BSF']:.2f}, PCF {u['usage_PCF']:.2f}"
            if u["site_confirmed"]:
                text += " (supported)"
        if s.feature_id in assigned_usage:
            ub, up_, conf = assigned_usage[s.feature_id]
            text += f"\nassigned reads: BSF {ub:.2f}, PCF {up_:.2f}" + (" (supported)" if conf else "")
        if s.feature_id in labels:
            text += f"\n{labels[s.feature_id]}"
        c = counts.loc[s.feature_id] if s.feature_id in counts.index else None
        if c is not None:
            text += "\nends " + "/".join(str(int(c[r])) for r in reps if stage_of[r] == "BSF") + " | " + "/".join(str(int(c[r])) for r in reps if stage_of[r] == "PCF")
        ax_model.text((s.genomic_start + s.genomic_end) / 2, 0.98 - 0.45 * (i % 3), text, ha="center", va="top", fontsize=7)
    ax_model.axvline(cds3, color="black", lw=1.5)
    ax_model.text(cds3, 1.0, " CDS 3' end", ha="left", va="top", fontsize=8)
    ax_model.set_ylim(-0.3, 1.05); ax_model.set_yticks([])
    ax_model.set_ylabel("tested\nclusters", rotation=0, ha="right", va="center")
    other = counts.loc["other"] if "other" in counts.index else None
    # reads per stage, sorted by 3' end, drawn as lines clipped to the frame
    for ax, k in ((ax_bsf, "BSF"), (ax_pcf, "PCF")):
        # longest reads at the top of the panel on both strands
        rr = sorted(reads[k], key=lambda x: x[2], reverse=(strand == "-"))
        total = len(rr)
        if total > args.max_reads:
            idx = np.linspace(0, total - 1, args.max_reads).astype(int)
            rr = [rr[i] for i in idx]
        n_unc = 0
        for row, (s, e, end3, _, r) in enumerate(rr):
            if index is not None:
                covers, _reason = assignment(r, args.gene, cs, ce, strand, index, args.min_cds_overlap, args.min_mapq, placements)
            else:
                covers = r.get_overlap(cs, ce) >= args.min_cds_overlap and r.mapping_quality >= args.min_mapq
            n_unc += not covers
            ax.plot([max(s, lo), min(e, hi)], [row, row], color=COL[k] if covers else "0.6", lw=0.9, alpha=0.8)
            ax.plot(end3, row, marker="|", color="black", ms=5, mew=1.0)
        ax.set_ylim(-1, max(len(rr), 1)); ax.set_yticks([])
        ax.set_ylabel(f"{k} reads ending\nin the window\n({len(rr)} of {total} drawn,\n{n_unc} grey: uncertain origin)", fontsize=8)
    # end counts per replicate (counts, not fractions), only tested reads
    for ax, rep in zip(ax_reps, reps):
        arr = np.zeros(hi - lo, dtype=int)
        for k in reads:
            for s, e, end3, name, _ in reads[k]:
                if rep_of.get(name) == rep and lo <= end3 < hi:
                    arr[end3 - lo] += 1
        ax.fill_between(np.arange(lo, hi), arr, step="mid", color=COL[stage_of[rep]], alpha=0.9)
        n_win = int(sum(1 for k in reads for *_, name, _r in reads[k] if rep_of.get(name) == rep))
        ax.set_ylabel(f"{rep}\n{n_win} ends", rotation=0, ha="right", va="center", fontsize=8)
        ax.set_yticks([])
    axes[-1].set_xlabel(f"{chrom} position (bp), 5' to 3'" + (" (minus strand, axis reversed)" if strand == "-" else ""))
    axes[-1].set_xlim((hi, lo) if strand == "-" else (lo, hi))
    n_bsf, n_pcf = len(reads["BSF"]), len(reads["PCF"])
    other_text = f"; other (ends outside the clusters) {int(other[[r for r in reps if stage_of[r]=='BSF']].sum())} BSF, {int(other[[r for r in reps if stage_of[r]=='PCF']].sum())} PCF" if other is not None else ""
    fig.suptitle(f"{args.gene} ({strand} strand), tested region. Window ends: {n_bsf} BSF, {n_pcf} PCF{other_text}.", fontsize=10)
    fig.savefig(args.out, dpi=110, bbox_inches="tight")
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
