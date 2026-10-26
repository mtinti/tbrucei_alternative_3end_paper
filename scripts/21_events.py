#!/usr/bin/env python3
"""Interpret every supported end cluster as an event (docs/methods.md, section 4).

One row per supported cluster of a supported gene (stageR decisions from step
18). The reads that end in the cluster are selected with the rules of step 02
(primary alignment, at least 200 aligned bases, 3' end = alignment end on the
plus strand and alignment start on the minus strand), so the evidence
population is the tested population.

Evidence columns:
- reads ending in the cluster, per stage and per replicate (--read-ids);
- informative reads: the alignment start is a transcript 5' end, because the
  read is in the full-length BAM (spliced leader detected) or starts within
  --link-distance bases of the annotated transcript start. Every other read is
  incomplete. With --sharp-start N (sensitivity analysis, off by default) a
  read is also informative when at least N other reads ending in the same
  cluster start within 10 bases of it;
- linked reads: informative reads that start within --link-distance bases of
  the gene's annotated transcript 5' end or upstream of the CDS 5' end (they
  span the CDS, so they carry the gene);
- catalogue record (step 22) whose 3' end lies within --tolerance bases on the
  same strand, with its full-length read counts per replicate and, when it
  was tested, the abundance result (step 23);
- annotated non-coding RNA (Ensembl BED) ending within --tolerance bases;
- trans-splice acceptor cluster downstream of the CDS 3' end: a 100-base bin
  of read starts holding at least --min-cluster reads and --min-fraction of
  the reads ending in the cluster, with an AG preceded by a pyrimidine-rich
  stretch within 60 bases upstream. A start cluster inside the gene (between
  the transcript 5' end and the CDS 3' end) is recorded separately as an
  intragenic alternative start and is not evidence for a separate transcript.

Transcript evidence = catalogue record, annotated non-coding RNA, or
downstream acceptor cluster. Label, in this order (shares are documented
heuristics):
  unresolved, insufficient evidence: fewer than --min-reads reads, or fewer
      than --min-informative informative reads;
  mixed: transcript evidence and at least one third of the informative reads
      linked to the gene;
  transcript-associated: transcript evidence, fewer than one third linked;
  gene-connected: no transcript evidence and at least two thirds linked;
  unresolved, incompatible evidence: otherwise.
A gene has no class; its most significant site is only a display choice.
"""
import argparse
import os
import re
import numpy as np
import pandas as pd
import pysam

PYR = set("CT")


def reads_ending_in(bam, chrom, start, end, strand, min_aligned):
    """Primary reads on the strand, with the step 02 eligibility, whose 3' end lies in [start, end)."""
    out = []
    for r in bam.fetch(chrom, max(0, start - 1), end + 1):
        if r.is_unmapped or r.is_secondary or r.is_supplementary or r.query_alignment_length < min_aligned:
            continue
        if r.is_reverse != (strand == "-"):
            continue
        end3 = r.reference_start if strand == "-" else r.reference_end
        if start <= end3 < end:
            out.append(r)
    return out


def acceptor_at(fasta, chrom, strand, mode):
    """AG with a pyrimidine-rich 15-base stretch before it, within 60 bases upstream of mode, in transcript direction."""
    if strand == "+":
        seq = fasta.fetch(chrom, max(0, mode - 60), mode + 2).upper()
        for i in range(len(seq) - 2, 0, -1):
            if seq[i:i + 2] == "AG":
                up = seq[max(0, i - 15):i]
                frac = sum(ch in PYR for ch in up) / max(1, len(up))
                if frac >= 0.6:
                    return max(0, mode - 60) + i, round(frac, 2)
    else:
        seq = fasta.fetch(chrom, max(0, mode - 1), mode + 61).upper()
        for i in range(0, len(seq) - 1):
            if seq[i:i + 2] == "CT":
                up = seq[i + 2:i + 17]
                frac = sum(ch in "AG" for ch in up) / max(1, len(up))
                if frac >= 0.6:
                    return max(0, mode - 1) + i, round(frac, 2)
    return "", ""


def start_clusters(starts, strand, gene5, cds3, fasta, chrom, n_reads, min_cluster, min_fraction, bin_size=100):
    """Strongest start bin downstream of the CDS 3' end (candidate separate transcript) and strongest bin inside the gene."""
    out = dict(downstream_start_cluster_pos="", downstream_start_cluster_reads=0, acceptor_AG="", pyrimidine_fraction="",
               acceptor_cluster=False, intragenic_start_cluster_pos="", intragenic_start_cluster_reads=0)
    if not starts:
        return out
    arr = np.array(starts)
    if strand == "+":
        down, inside = arr[arr >= cds3], arr[(arr > gene5 + 100) & (arr < cds3)]
    else:
        down, inside = arr[arr < cds3], arr[(arr < gene5 - 100) & (arr >= cds3)]
    for name, sub in (("downstream", down), ("intragenic", inside)):
        if len(sub) == 0:
            continue
        b, c = np.unique(sub // bin_size, return_counts=True)
        members = sub[(sub // bin_size) == b[c.argmax()]]
        mode = int(pd.Series(members).mode().iloc[0]); n = int(c.max())
        out[f"{name}_start_cluster_pos"], out[f"{name}_start_cluster_reads"] = mode, n
        if name == "downstream":
            ag, pyr = acceptor_at(fasta, chrom, strand, mode)
            out["acceptor_AG"], out["pyrimidine_fraction"] = ag, pyr
            out["acceptor_cluster"] = bool(ag != "" and n >= min_cluster and n >= min_fraction * max(1, n_reads))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dexseq", required=True, help="step 17/18 folder (sites.tsv, dexseq_site.tsv, dexseq_gene.tsv)")
    ap.add_argument("--family", required=True)
    ap.add_argument("--transcript-starts", required=True)
    ap.add_argument("--polya-bams", nargs=2, required=True, metavar=("BSF", "PCF"))
    ap.add_argument("--full-length-bams", nargs=2, required=True, metavar=("BSF", "PCF"))
    ap.add_argument("--read-ids", nargs="*", default=[], help="replicate=path, poly(A) read ids per replicate")
    ap.add_argument("--catalogue", help="catalogue.tsv from step 22")
    ap.add_argument("--abundance", help="transcript_abundance.tsv from step 23")
    ap.add_argument("--ncrna")
    ap.add_argument("--genome", required=True)
    ap.add_argument("--min-aligned-length", type=int, default=200, help="as in step 02")
    ap.add_argument("--tolerance", type=int, default=60)
    ap.add_argument("--link-distance", type=int, default=100)
    ap.add_argument("--min-reads", type=int, default=20)
    ap.add_argument("--min-informative", type=int, default=10)
    ap.add_argument("--min-cluster", type=int, default=10)
    ap.add_argument("--min-fraction", type=float, default=0.15)
    ap.add_argument("--sharp-start", type=int, default=0, help="sensitivity analysis: reads with this many companions within 10 bases count as informative (0 = off)")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    sites = pd.read_csv(os.path.join(args.dexseq, "sites.tsv"), sep="\t")
    st = pd.read_csv(os.path.join(args.dexseq, "dexseq_site.tsv"), sep="\t")
    g = pd.read_csv(os.path.join(args.dexseq, "dexseq_gene.tsv"), sep="\t")
    supported_genes = set(g.loc[g["gene_confirmed"], "gene_id"])
    ev = st[st["site_confirmed"] & (st["feature_id"] != "other") & st["gene_id"].isin(supported_genes)].merge(sites, on=["gene_id", "feature_id"])
    fam = pd.read_csv(args.family, sep="\t").set_index("gene")
    tstart = pd.read_csv(args.transcript_starts, sep="\t").set_index("gene")["transcript_5prime"]
    fasta = pysam.FastaFile(args.genome)
    polya = {"BSF": pysam.AlignmentFile(args.polya_bams[0]), "PCF": pysam.AlignmentFile(args.polya_bams[1])}
    full = {"BSF": pysam.AlignmentFile(args.full_length_bams[0]), "PCF": pysam.AlignmentFile(args.full_length_bams[1])}
    rep_of = {}
    for spec in args.read_ids:
        rep, path = spec.split("=", 1)
        for line in open(path):
            rep_of[line.strip()] = rep
    reps = sorted(set(rep_of.values()))
    cat = None
    if args.catalogue:
        cat = pd.read_csv(args.catalogue, sep="\t")
        if args.abundance:
            ab = pd.read_csv(args.abundance, sep="\t")[["record", "log2FC_PCF_vs_BSF", "FDR", "abundance_changed"]]
            cat = cat.merge(ab, on="record", how="left")
        # one row per member record, so that a combined locus is matched through any of its members
        mem = []
        for r in cat.itertuples():
            for m in str(r.members).split(";"):
                mm = re.match(r"(.+):(\d+)-(\d+)\(([+-])\)", m)
                if mm:
                    mem.append(dict(unit_index=r.Index, member=m, chrom=mm.group(1), strand=mm.group(4),
                                    end3=int(mm.group(2)) if mm.group(4) == "-" else int(mm.group(3))))
        members = pd.DataFrame(mem)
    nc = None
    if args.ncrna:
        nc = pd.read_csv(args.ncrna, sep="\t", header=None, usecols=range(6), names=["chrom", "start", "end", "id", "biotype", "strand"])
        nc = nc[nc["biotype"] != "pseudogene"]
        nc["end3"] = np.where(nc["strand"] == "-", nc["start"], nc["end"])

    rows = []
    for s in ev.itertuples():
        f = fam.loc[s.gene_id]
        chrom, strand = s.chromosome, s.strand
        cs, ce = int(f.cds_start), int(f.cds_end)
        gene5_cds = cs if strand == "+" else ce - 1
        cds3 = ce if strand == "+" else cs
        has_tr = s.gene_id in tstart.index
        gene5 = int(tstart[s.gene_id]) if has_tr else gene5_cds
        row = dict(gene_id=s.gene_id, feature_id=s.feature_id, site_location=f"{chrom}:{s.genomic_start}-{s.genomic_end}({strand})",
                   site_pvalue=s.pvalue, site_padj=s.site_padj, usage_BSF=s.usage_BSF, usage_PCF=s.usage_PCF,
                   usage_change_points=s.usage_change_points,
                   transcript_5prime_reference="annotated transcript" if has_tr else "none (CDS start recorded as a proxy; linkage from spliced-leader reads only)")
        sl_names = set()
        for bam in full.values():
            sl_names |= {r.query_name for r in reads_ending_in(bam, chrom, s.genomic_start, s.genomic_end, strand, args.min_aligned_length)}
        collected, n_stage = [], {}
        for stage, bam in polya.items():
            reads = reads_ending_in(bam, chrom, s.genomic_start, s.genomic_end, strand, args.min_aligned_length)
            n_stage[stage] = len(reads)
            for r in reads:
                start5 = r.reference_start if strand == "+" else r.reference_end - 1
                collected.append((r.query_name, start5, rep_of.get(r.query_name, "")))
        starts_arr = np.array(sorted(x[1] for x in collected))

        def sharp(pos):
            return args.sharp_start > 0 and (np.searchsorted(starts_arr, pos + 10, "right") - np.searchsorted(starts_arr, pos - 10, "left")) - 1 >= args.sharp_start

        per = {k: {r: 0 for r in reps} for k in ("reads", "informative", "linked")}
        informative = linked = incomplete = 0
        for name, start5, rep in collected:
            if rep:
                per["reads"][rep] += 1
            is_sl = name in sl_names
            near_tss = has_tr and abs(start5 - gene5) <= args.link_distance
            if is_sl or near_tss or sharp(start5):
                informative += 1
                if rep:
                    per["informative"][rep] += 1
                upstream_of_cds = (start5 <= gene5_cds) if strand == "+" else (start5 >= gene5_cds)
                if near_tss or upstream_of_cds:
                    linked += 1
                    if rep:
                        per["linked"][rep] += 1
            else:
                incomplete += 1
        n = len(collected)
        row.update(reads_BSF=n_stage["BSF"], reads_PCF=n_stage["PCF"], reads_total=n,
                   informative_reads=informative, incomplete_reads=incomplete, linked_reads=linked, spliced_leader_reads=len(sl_names),
                   linkage_share=round(linked / informative, 3) if informative >= args.min_informative else "")
        for k in ("reads", "informative", "linked"):
            row.update({f"{k}_{r}": per[k][r] for r in reps})
        # catalogue record at the cluster end, with its abundance result
        row.update(catalogue_record="", catalogue_member="", catalogue_kind="", catalogue_reads="", catalogue_abundance="")
        if cat is not None:
            hm = members[(members["chrom"] == chrom) & (members["strand"] == strand) & (members["end3"] >= s.genomic_start - args.tolerance)
                         & (members["end3"] <= s.genomic_end + args.tolerance)]
            if len(hm):
                units = cat.loc[hm["unit_index"].unique()]
                h = units.loc[units[reps].sum(axis=1).idxmax()] if reps and all(r in units.columns for r in reps) else units.iloc[0]
                row["catalogue_record"], row["catalogue_kind"] = h["record"], h["kind"]
                row["catalogue_member"] = ";".join(hm.loc[hm["unit_index"] == h.name, "member"])
                row["catalogue_reads"] = ",".join(str(int(h[r])) for r in reps if r in h)
                if "FDR" in h and not pd.isna(h["FDR"]):
                    row["catalogue_abundance"] = f"log2FC {h['log2FC_PCF_vs_BSF']:+.2f}, FDR {h['FDR']:.3g}" + (" (changed)" if h["abundance_changed"] else "")
                else:
                    row["catalogue_abundance"] = "not tested (below the expression filter)"
        hit_nc = ""
        if nc is not None:
            h = nc[(nc["chrom"] == chrom) & (nc["strand"] == strand) & (nc["end3"] >= s.genomic_start - args.tolerance) & (nc["end3"] <= s.genomic_end + args.tolerance)]
            hit_nc = "; ".join(f"{a} ({b})" for a, b in zip(h["id"], h["biotype"]))
        sc = start_clusters([x[1] for x in collected], strand, gene5, cds3, fasta, chrom, n, args.min_cluster, args.min_fraction)
        row.update(annotated_ncRNA=hit_nc, **sc)
        transcript_evidence = bool(row["catalogue_record"] or hit_nc or sc["acceptor_cluster"])
        share = linked / informative if informative else 0
        if n < args.min_reads or informative < args.min_informative:
            label = "unresolved, insufficient evidence"
        elif transcript_evidence and share >= 1 / 3:
            label = "mixed"
        elif transcript_evidence:
            label = "transcript-associated"
        elif share >= 2 / 3:
            label = "gene-connected"
        else:
            label = "unresolved, incompatible evidence"
        row["label"] = label
        rows.append(row)
    out = pd.DataFrame(rows).sort_values(["gene_id", "site_pvalue"])
    out.to_csv(args.out, sep="\t", index=False)
    print(f"wrote {args.out}: {len(out)} events in {out['gene_id'].nunique()} genes"
          + (f" (sensitivity: sharp-start {args.sharp_start})" if args.sharp_start else ""))
    print(out["label"].value_counts().to_string())


if __name__ == "__main__":
    main()
