#!/usr/bin/env python3
"""Build the catalogue of independent transcripts and count their full-length reads per replicate (docs/methods.md, section 5).

Without reading any usage result. A transcript record is a strand, a start
range and an end range. Records come from:

- the full-length reads (spliced leader and poly(A) tail; primary alignments)
  that touch no CDS, grouped by strand when their starts lie within
  --start-tolerance bases and their ends within --end-tolerance bases of a
  record's bounds (greedy, in start order);
- the annotated non-coding RNA genes (Ensembl BED, pseudogenes excluded).

An annotated record and a read-derived record are the same unit when their
starts and ends agree within the same tolerances on the same strand; the
annotated id is kept and the read-derived bounds are recorded. Records that
overlap but differ beyond the tolerances stay separate. Records that share
any ambiguous read (a read compatible with both) are joined, transitively,
into one combined locus, and every read of every member, unique or
ambiguous, counts to that locus. Counting units are therefore disjoint and
complete before the expression filter. Reads compatible with no record are
unassigned.

Protein-coding genes get their own full-length read counts (reads on the
gene's strand overlapping its CDS and no other); they serve only as the
normalisation reference in the abundance test.

Outputs in --outdir: catalogue.tsv (record or locus, members, bounds,
annotated ids, counts per replicate, in_family after the expression filter),
gene_full_length_counts.tsv, unassigned_reads.tsv (per replicate: reads
touching no CDS and no record).
"""
import argparse
import bisect
import os
import numpy as np
import pandas as pd
import pysam


def load_cds(path):
    cds = pd.read_csv(path, sep="\t", header=None, usecols=range(6), names=["chrom", "start", "end", "gene", "x", "strand"])
    idx = {}
    for chrom, g in cds.groupby("chrom"):
        g = g.sort_values("start")
        idx[chrom] = (g["start"].to_numpy(), g["end"].to_numpy(), g["strand"].to_numpy(), g["gene"].to_numpy(),
                      np.maximum.accumulate(g["end"].to_numpy()))
    return idx


def overlapping_cds(idx, chrom, start, end):
    if chrom not in idx:
        return []
    starts, ends, strands, genes, maxend = idx[chrom]
    hi = bisect.bisect_left(starts, end)
    lo = bisect.bisect_right(maxend, start)
    return [(genes[i], strands[i]) for i in range(lo, hi) if ends[i] > start]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--full-length-bams", nargs="+", required=True)
    ap.add_argument("--cds-bed", required=True)
    ap.add_argument("--ncrna", required=True)
    ap.add_argument("--read-ids", nargs="+", required=True, help="replicate=path")
    ap.add_argument("--start-tolerance", type=int, default=100)
    ap.add_argument("--end-tolerance", type=int, default=60)
    ap.add_argument("--min-reads", type=int, default=5)
    ap.add_argument("--min-replicates", type=int, default=2)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    rep_of = {}
    for spec in args.read_ids:
        rep, path = spec.split("=", 1)
        for line in open(path):
            rep_of[line.strip()] = rep
    reps = sorted(set(rep_of.values()))
    cds_idx = load_cds(args.cds_bed)

    free, gene_counts, unmapped_rep = [], {}, 0
    for path in args.full_length_bams:
        with pysam.AlignmentFile(path) as bam:
            for r in bam.fetch():
                if r.is_unmapped or r.is_secondary or r.is_supplementary or r.reference_name not in cds_idx:
                    continue
                rep = rep_of.get(r.query_name)
                if rep is None:
                    unmapped_rep += 1
                    continue
                strand = "-" if r.is_reverse else "+"
                hits = overlapping_cds(cds_idx, r.reference_name, r.reference_start, r.reference_end)
                if not hits:
                    free.append((r.reference_name, strand, r.reference_start, r.reference_end, rep))
                elif len(hits) == 1 and hits[0][1] == strand:
                    gene_counts.setdefault(hits[0][0], {x: 0 for x in reps})[rep] += 1
    print(f"{len(free)} full-length reads touch no CDS; {unmapped_rep} reads could not be mapped to a replicate")
    if unmapped_rep:
        raise SystemExit("every full-length read must map to exactly one replicate")
    reads = pd.DataFrame(free, columns=["chrom", "strand", "start", "end", "rep"])
    reads["five"] = np.where(reads["strand"] == "+", reads["start"], reads["end"])
    reads["three"] = np.where(reads["strand"] == "+", reads["end"], reads["start"])

    # read-derived records: greedy grouping in 5' order
    records = []   # dict(chrom, strand, five_lo, five_hi, three_lo, three_hi, n)
    reads["record"] = -1
    for (chrom, strand), grp in reads.groupby(["chrom", "strand"]):
        # full order (5' end, 3' end, replicate) with a stable sort, so the greedy grouping does not depend on the platform's tie order
        grp = grp.sort_values(["five", "three", "rep"], kind="mergesort")
        open_recs = []
        for i, r in grp.iterrows():
            placed = False
            for k in open_recs:
                rec = records[k]
                if abs(r["five"] - rec["five_mode"]) <= args.start_tolerance and abs(r["three"] - rec["three_mode"]) <= args.end_tolerance:
                    rec["n"] += 1
                    rec["five_lo"], rec["five_hi"] = min(rec["five_lo"], r["five"]), max(rec["five_lo"], r["five"])
                    rec["three_lo"], rec["three_hi"] = min(rec["three_lo"], r["three"]), max(rec["three_hi"], r["three"])
                    reads.at[i, "record"] = k
                    placed = True
                    break
            if not placed:
                records.append(dict(chrom=chrom, strand=strand, five_mode=r["five"], three_mode=r["three"], five_lo=r["five"],
                                    five_hi=r["five"], three_lo=r["three"], three_hi=r["three"], n=1, annotated=""))
                reads.at[i, "record"] = len(records) - 1
                open_recs.append(len(records) - 1)
            open_recs = [k for k in open_recs if r["five"] - records[k]["five_mode"] <= args.start_tolerance]
    rec = pd.DataFrame(records)
    print(f"{len(rec)} read-derived records")

    # annotated non-coding RNAs: same unit when both ends agree, else their own record
    nc = pd.read_csv(args.ncrna, sep="\t", header=None, usecols=range(6), names=["chrom", "start", "end", "id", "biotype", "strand"])
    nc = nc[(nc["biotype"] != "pseudogene") & nc["chrom"].isin(cds_idx)]
    nc["five"] = np.where(nc["strand"] == "+", nc["start"], nc["end"])
    nc["three"] = np.where(nc["strand"] == "+", nc["end"], nc["start"])
    new_rows = []
    for a in nc.itertuples():
        m = rec[(rec["chrom"] == a.chrom) & (rec["strand"] == a.strand) & ((rec["five_mode"] - a.five).abs() <= args.start_tolerance)
                & ((rec["three_mode"] - a.three).abs() <= args.end_tolerance)]
        if len(m):
            k = m.index[0]
            rec.at[k, "annotated"] = (rec.at[k, "annotated"] + ";" if rec.at[k, "annotated"] else "") + f"{a.id} ({a.biotype})"
        else:
            new_rows.append(dict(chrom=a.chrom, strand=a.strand, five_mode=a.five, three_mode=a.three, five_lo=a.five, five_hi=a.five,
                                 three_lo=a.three, three_hi=a.three, n=0, annotated=f"{a.id} ({a.biotype})"))
    rec = pd.concat([rec, pd.DataFrame(new_rows)], ignore_index=True)
    rec["start"] = np.where(rec["strand"] == "+", rec[["five_lo", "five_mode"]].min(axis=1), rec[["three_lo"]].min(axis=1))
    rec["end"] = np.where(rec["strand"] == "+", rec[["three_hi"]].max(axis=1), rec[["five_hi", "five_mode"]].max(axis=1))
    rec["record_id"] = [f"{c}:{int(s)}-{int(e)}({st})" for c, s, e, st in zip(rec["chrom"], rec["start"], rec["end"], rec["strand"])]

    # assign every free read to its compatible records; records sharing a read join one unit (union-find)
    parent = {rid: rid for rid in rec["record_id"]}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    assignments = []   # (list of compatible record ids, rep)
    n_unassigned = {r: 0 for r in reps}
    for (chrom, strand), grp in reads.groupby(["chrom", "strand"]):
        sub = rec[(rec["chrom"] == chrom) & (rec["strand"] == strand)]
        fives, threes, ids = sub["five_mode"].to_numpy(), sub["three_mode"].to_numpy(), sub["record_id"].to_numpy()
        for r in grp.itertuples():
            ok = np.where((np.abs(fives - r.five) <= args.start_tolerance) & (np.abs(threes - r.three) <= args.end_tolerance))[0]
            if len(ok) == 0:
                n_unassigned[r.rep] += 1
                continue
            for k in ok[1:]:
                union(ids[ok[0]], ids[k])
            assignments.append((ids[ok[0]], r.rep))
    units = {}
    for rid in rec["record_id"]:
        units.setdefault(find(rid), []).append(rid)
    counts = {u: {x: 0 for x in reps} for u in units}
    for rid, rep in assignments:
        counts[find(rid)][rep] += 1
    cat = rec.set_index("record_id")
    rows = []
    for u, members in units.items():
        c = counts[u]
        members = sorted(members)
        if sum(c.values()) == 0 and not any(cat.loc[members, "annotated"]):
            continue
        if len(members) == 1:
            key = members[0]
            rows.append(dict(record=key, kind="transcript", members=key, chrom=cat.loc[key, "chrom"], start=int(cat.loc[key, "start"]),
                             end=int(cat.loc[key, "end"]), strand=cat.loc[key, "strand"], annotated=cat.loc[key, "annotated"], **c))
        else:
            key = "locus:" + "+".join(members)
            rows.append(dict(record=key, kind="combined locus", members=";".join(members), chrom=cat.loc[members[0], "chrom"],
                             start=int(cat.loc[members, "start"].min()), end=int(cat.loc[members, "end"].max()),
                             strand=cat.loc[members[0], "strand"], annotated=";".join(a for a in cat.loc[members, "annotated"] if a), **c))
    out = pd.DataFrame(rows)
    out["in_family"] = (out[reps] >= args.min_reads).sum(axis=1) >= args.min_replicates
    out.sort_values(["chrom", "start"]).to_csv(os.path.join(args.outdir, "catalogue.tsv"), sep="\t", index=False)
    pd.DataFrame([dict(gene_id=g, **c) for g, c in gene_counts.items()]).to_csv(
        os.path.join(args.outdir, "gene_full_length_counts.tsv"), sep="\t", index=False)
    pd.DataFrame([dict(replicate=r, unassigned_reads=n) for r, n in n_unassigned.items()]).to_csv(
        os.path.join(args.outdir, "unassigned_reads.tsv"), sep="\t", index=False)
    print(f"catalogue: {len(out)} records ({int((out['kind'] == 'combined locus').sum())} combined loci, "
          f"{int((out['annotated'] != '').sum())} with an annotated non-coding RNA); {int(out['in_family'].sum())} pass the expression filter; "
          f"unassigned reads: {sum(n_unassigned.values())}; {len(gene_counts)} protein-coding genes counted")


if __name__ == "__main__":
    main()
