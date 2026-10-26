#!/usr/bin/env python3
"""Count poly(A) ends per cluster among the reads that carry the gene (docs/methods.md, section 7).

The window test (Parts 2 and 3) compares a cluster with the rest of its
window, and the window can hold ends of other transcripts. This step keeps the
same windows and the same clusters but counts only the reads assigned to the
gene (assignlib.py): a read with the step 02 eligibility (primary alignment,
same strand, at least --min-aligned-length aligned bases, 3' end inside the
window) with at least --min-cds-overlap aligned bases inside the gene's CDS,
a resolved placement (mapping quality at least --min-mapq and, with
--placements from step 31, no secondary alignment scoring at least as high as
the primary), and no other same-strand CDS receiving --min-cds-overlap
aligned bases from it. Such a read contains coding sequence of this gene and
no other, at a placement no other placement matches, so its poly(A) end is an
end of the gene's transcript. Every other read is of uncertain origin (no CDS
overlap, low mapping quality, equal placement elsewhere, covers another CDS);
these are counted separately, never in the denominator. The rule does not
look at the stage.

Per window: counts per cluster and replicate over the assigned reads; other
is the assigned ends outside the clusters; a cluster is supported with
--support-reads assigned ends in --support-replicates replicates, otherwise
it folds into other; a gene is eligible with --min-ends assigned ends in
every replicate and two or more supported clusters. The identity clusters
plus other equal the assigned ends of the window is checked per replicate.

Output folder, in the layout step 18 and step 20 expect: site_counts.tsv,
samples.tsv, sites.tsv, windows.tsv (eligibility with the reason when a gene
is not eligible), plus assignment.tsv with the assigned and the uncertain
ends per cluster and replicate.
"""
import argparse
import os
import sys
import pandas as pd
import pysam
from assignlib import REASONS, CdsIndex, Placements, assignment


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--family", required=True, help="family_<mode>.tsv from step 16")
    ap.add_argument("--clusters", required=True, help="step 17 folder of one setting (sites.tsv, windows.tsv)")
    ap.add_argument("--polya-bams", nargs=2, required=True, metavar=("BSF", "PCF"))
    ap.add_argument("--read-ids", nargs="+", required=True, help="replicate=path, poly(A) read ids per replicate")
    ap.add_argument("--min-aligned-length", type=int, default=200, help="as in step 02")
    ap.add_argument("--cds-bed", required=True, help="cds_v2.bed from step 16, for the other-gene check")
    ap.add_argument("--min-cds-overlap", type=int, default=50, help="aligned bases inside the CDS")
    ap.add_argument("--min-mapq", type=int, default=20, help="mapping quality below this is a low mapping quality")
    ap.add_argument("--placements", help="read_placements.tsv.gz from step 31; without it the secondary-alignment rule is not applied")
    ap.add_argument("--min-ends", type=int, default=10, help="assigned ends per replicate for eligibility")
    ap.add_argument("--support-reads", type=int, default=5)
    ap.add_argument("--support-replicates", type=int, default=2)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    fam = pd.read_csv(args.family, sep="\t")
    fam = fam[fam["in_family"]].set_index("gene")
    sites = pd.read_csv(os.path.join(args.clusters, "sites.tsv"), sep="\t")
    win = pd.read_csv(os.path.join(args.clusters, "windows.tsv"), sep="\t").set_index("gene")
    rep_of = {}
    for spec in args.read_ids:
        rep, path = spec.split("=", 1)
        for line in open(path):
            rep_of[line.strip()] = rep
    reps = sorted(set(rep_of.values()))
    stage_of = {r: ("BSF" if r.startswith("bs") else "PCF") for r in reps}
    bams = [pysam.AlignmentFile(b) for b in args.polya_bams]
    index = CdsIndex(args.cds_bed)
    placements = Placements(args.placements) if args.placements else None

    count_rows, site_rows, win_rows, assign_rows = [], [], [], []
    n_assigned = n_uncertain = 0
    for gene, f in fam.iterrows():
        chrom, strand = f["chromosome"], f["strand"]
        cs, ce = int(f["cds_start"]), int(f["cds_end"])
        if gene in win.index:
            ws, we = int(win.loc[gene, "window_start"]), int(win.loc[gene, "window_end"])
        else:
            ws, we = int(f["window_start"]), int(f["window_end"])
        gs = sites[sites["gene_id"] == gene].sort_values("genomic_start")
        bounds = list(zip(gs["feature_id"], gs["genomic_start"].astype(int), gs["genomic_end"].astype(int)))
        assigned = {r: 0 for r in reps}
        unc_reason = {k: 0 for k in REASONS}
        per_cluster = {fid: {r: 0 for r in reps} for fid, _, _ in bounds}
        per_cluster_unc = {fid: {r: 0 for r in reps} for fid, _, _ in bounds}
        per_cluster_reason = {fid: {k: 0 for k in REASONS} for fid, _, _ in bounds}
        seen = set()
        for bam in bams:
            for r in bam.fetch(chrom, max(0, ws - 1), we + 1):
                if r.is_unmapped or r.is_secondary or r.is_supplementary or r.query_alignment_length < args.min_aligned_length:
                    continue
                if r.is_reverse != (strand == "-") or r.query_name in seen:
                    continue
                end3 = r.reference_start if strand == "-" else r.reference_end
                if not (ws <= end3 < we):
                    continue
                rep = rep_of.get(r.query_name)
                if rep is None:
                    continue
                seen.add(r.query_name)
                is_assigned, reason = assignment(r, gene, cs, ce, strand, index, args.min_cds_overlap, args.min_mapq, placements)
                if is_assigned:
                    assigned[rep] += 1
                    n_assigned += 1
                else:
                    n_uncertain += 1
                    unc_reason[reason] += 1
                for fid, a, b in bounds:
                    if a <= end3 < b:
                        (per_cluster if is_assigned else per_cluster_unc)[fid][rep] += 1
                        if not is_assigned:
                            per_cluster_reason[fid][reason] += 1
                        break
        supported = {fid: sum(per_cluster[fid][r] >= args.support_reads for r in reps) >= args.support_replicates for fid, _, _ in bounds}
        n_sup = sum(supported.values())
        other = {r: assigned[r] - sum(per_cluster[fid][r] for fid in per_cluster if supported[fid]) for r in reps}
        for r in reps:
            if other[r] < 0:
                sys.exit(f"count identity broken for {gene} {r}")
        depth_ok = all(assigned[r] >= args.min_ends for r in reps)
        reason = ""
        if not depth_ok:
            reason = f"fewer than {args.min_ends} assigned ends in a replicate"
        elif n_sup < 2:
            reason = "fewer than two supported clusters among the assigned ends"
        eligible = reason == ""
        win_rows.append(dict(gene=gene, chromosome=chrom, strand=strand, window_start=ws, window_end=we,
                             pooled_ends=sum(assigned.values()), clusters=len(bounds), supported_clusters=n_sup,
                             testable=eligible, exclusion_reason=reason,
                             **{f"assigned_{r}": assigned[r] for r in reps},
                             **{"uncertain_" + k.replace(" ", "_"): v for k, v in unc_reason.items()}))
        for fid, a, b in bounds:
            assign_rows.append(dict(gene_id=gene, feature_id=fid, supported_assigned=supported[fid],
                                    **{f"assigned_{r}": per_cluster[fid][r] for r in reps},
                                    **{f"uncertain_{r}": per_cluster_unc[fid][r] for r in reps},
                                    **{"uncertain_" + k.replace(" ", "_"): v for k, v in per_cluster_reason[fid].items()}))
        if not eligible:
            continue
        for fid, a, b in bounds:
            if supported[fid]:
                count_rows.append(dict(gene_id=gene, feature_id=fid, **per_cluster[fid]))
                site_rows.append(dict(gene_id=gene, feature_id=fid, chromosome=chrom, genomic_start=a, genomic_end=b, strand=strand))
        count_rows.append(dict(gene_id=gene, feature_id="other", **other))

    pd.DataFrame(count_rows).to_csv(os.path.join(args.outdir, "site_counts.tsv"), sep="\t", index=False)
    pd.DataFrame(site_rows).to_csv(os.path.join(args.outdir, "sites.tsv"), sep="\t", index=False)
    w = pd.DataFrame(win_rows)
    w.to_csv(os.path.join(args.outdir, "windows.tsv"), sep="\t", index=False)
    pd.DataFrame(assign_rows).to_csv(os.path.join(args.outdir, "assignment.tsv"), sep="\t", index=False)
    pd.DataFrame({"sample_id": reps, "stage": [stage_of[r] for r in reps]}).to_csv(os.path.join(args.outdir, "samples.tsv"), sep="\t", index=False)
    reasons_text = ", ".join(f"{int(w['uncertain_' + k.replace(' ', '_')].sum())} {k}" for k in REASONS)
    print(f"{len(fam)} family genes; {n_assigned} assigned and {n_uncertain} uncertain reads ({reasons_text}); "
          f"{int(w['testable'].sum())} eligible genes -> {args.outdir}")
    print(w.loc[~w["testable"], "exclusion_reason"].value_counts().to_string())


if __name__ == "__main__":
    main()
