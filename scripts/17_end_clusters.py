#!/usr/bin/env python3
"""Define end clusters from the pooled read ends of all replicates (docs/methods.md, section 2).

No stage label is read. For each window in the family:

1. Pool the end positions of every replicate. The window total is the number
   of these ends; the peak percentage refers to it.
2. Smooth the per-base counts with a centred running sum of --smooth bases.
   A candidate peak is a local maximum of the smoothed count that holds at
   least --peak-fraction of the window total and at least --peak-min ends.
   On a plateau (a run of equal smoothed values, which a single sharp end
   position produces) the centre base is the peak.
3. Merge neighbouring peaks. Peaks closer than --always-merge bases are one
   site (nanopore 3' ends scatter over about ten bases around a site). Peaks
   between --always-merge and --merge bases apart merge only when the smoothed
   count between them never drops below --valley times the lower peak (one
   hump); a deeper dip keeps them apart. Peaks further than --merge apart never
   merge. Merging chains. The core of a cluster is the span from its first to
   its last peak.
4. Extend each core by --extend bases on each side, cut at the midpoint
   between neighbouring cores where two extensions would overlap, clip at the
   window bounds. Clusters are disjoint.
5. A cluster is supported when it holds at least --support-reads reads in at
   least --support-replicates replicates. Unsupported clusters are folded into
   "other" before anything is fitted. Nothing is removed after fitting.
6. A window is testable when it has at least two supported clusters;
   otherwise it is excluded with the reason "fewer than two resolved supported
   clusters".

The script checks, for every replicate of every window, that the cluster
counts plus other equal the window total, and stops if the identity fails.

Outputs in --outdir: clusters.bed (chrom, start, end, gene|cluster, pooled
ends, strand, supported), site_counts.tsv (gene_id, feature_id, one column per
replicate; only testable windows, supported clusters plus other), sites.tsv
(the cluster coordinates in the layout step 11 expects) and windows.tsv (per
window: total ends, clusters found, supported, testable, reason).

--self-test runs the synthetic checks (two near peaks, a plateau, a peak at
the window edge, both strands) and exits.
"""
import argparse
import os
import sys
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
import utrlib  # noqa: E402


def find_peaks(counts, smooth, peak_fraction, peak_min):
    """Positions (0-based within the window) of the candidate peaks."""
    total = counts.sum()
    if total == 0:
        return []
    k = np.ones(smooth)
    sm = np.convolve(counts, k, mode="same")
    threshold = max(peak_fraction * total, peak_min)
    peaks = []
    n = len(sm)
    i = 0
    while i < n:
        if sm[i] < threshold:
            i += 1
            continue
        # plateau: run of equal values; the centre base is the peak
        j = i
        while j + 1 < n and sm[j + 1] == sm[i]:
            j += 1
        left = sm[i - 1] if i > 0 else -1
        right = sm[j + 1] if j + 1 < n else -1
        if sm[i] > left and sm[i] > right:
            peaks.append((i + j) // 2)
        i = j + 1
    return peaks


def merge_peaks(peaks, merge, smoothed=None, valley=0.5, always_merge=15):
    """Groups of neighbouring peaks, chained. Returns (core_start, core_end) per group.

    Two neighbours merge when they are closer than merge bases and, if the
    smoothed profile is given, the dip between them stays above valley times
    the lower peak.
    """
    cores = []
    for p in peaks:
        if cores and p - cores[-1][1] < merge:
            q = cores[-1][1]
            close = p - q < always_merge
            one_hump = smoothed is None or smoothed[q:p + 1].min() >= valley * min(smoothed[q], smoothed[p])
            if close or one_hump:
                cores[-1][1] = p
                continue
        cores.append([p, p])
    return cores


def cluster_bounds(cores, extend, length):
    """Disjoint [start, end) per core: extend, cut at midpoints, clip."""
    bounds = []
    for i, (cs, ce) in enumerate(cores):
        s, e = cs - extend, ce + extend + 1
        if i > 0:
            prev_end = cores[i - 1][1]
            s = max(s, (prev_end + cs) // 2 + 1)
        if i + 1 < len(cores):
            next_start = cores[i + 1][0]
            e = min(e, (ce + next_start) // 2 + 1)
        bounds.append((max(s, 0), min(e, length)))
    for a, b in zip(bounds, bounds[1:]):
        assert a[1] <= b[0], f"overlapping clusters {a} {b}"
    return bounds


def cluster_window(per_rep_counts, args):
    """per_rep_counts: {rep: array of end counts per base}. Returns bounds, supported flags, counts."""
    reps = list(per_rep_counts)
    pooled = np.sum([per_rep_counts[r] for r in reps], axis=0)
    peaks = find_peaks(pooled, args.smooth, args.peak_fraction, args.peak_min)
    smoothed = np.convolve(pooled, np.ones(args.smooth), mode="same")
    cores = merge_peaks(peaks, args.merge, smoothed, args.valley, args.always_merge)
    bounds = cluster_bounds(cores, args.extend, len(pooled))
    counts = {r: [int(per_rep_counts[r][s:e].sum()) for s, e in bounds] for r in reps}
    supported = [sum(counts[r][i] >= args.support_reads for r in reps) >= args.support_replicates
                 for i in range(len(bounds))]
    return pooled, bounds, supported, counts


def self_test():
    class A:
        smooth, peak_fraction, peak_min, merge, extend, support_reads, support_replicates, valley, always_merge = 11, 0.02, 5, 30, 15, 5, 2, 0.5, 15
    L = 400
    def window(spikes):
        c = np.zeros(L, dtype=int)
        for pos, n in spikes:
            c[pos] += n
        return c
    # two sharp peaks 25 bases apart with nothing between: separate at merge 30 (deep dip), disjoint bounds
    w = {"a": window([(100, 50), (125, 40)]), "b": window([(100, 45), (125, 35)])}
    _, b, s, _ = cluster_window(w, A)
    assert len(b) == 2 and b[0][1] <= b[1][0] and b[0] == (85, 113) and b[1] == (113, 141), b
    # two sharp peaks 12 bases apart: always one site
    w = {"a": window([(100, 50), (112, 40)]), "b": window([(100, 45), (112, 35)])}
    _, b, s, _ = cluster_window(w, A)
    assert len(b) == 1 and s == [True], (b, s)
    # two peaks 25 apart joined by a shallow dip (ends spread between them): one hump, one cluster
    w = {"a": window([(100, 30), (112, 20), (125, 30)]), "b": window([(100, 30), (112, 20), (125, 30)])}
    _, b, s, _ = cluster_window(w, A)
    assert len(b) == 1, b
    # plateau: equal smoothed counts over a run; centre base is the peak, one cluster
    w = {"a": window([(200 + i, 10) for i in range(5)]), "b": window([(200 + i, 10) for i in range(5)])}
    _, b, s, _ = cluster_window(w, A)
    assert len(b) == 1 and b[0][0] <= 200 and b[0][1] > 204 and 200 <= (b[0][0] + b[0][1]) // 2 <= 204, b
    # peak at the window edge: clipped to the window
    w = {"a": window([(0, 30), (399, 30)]), "b": window([(0, 30), (399, 30)])}
    _, b, s, _ = cluster_window(w, A)
    assert b[0][0] == 0 and b[-1][1] == L, b
    # two strong peaks 50 bases apart with a small bump between them (the Tb927.2.4850 case):
    # the bump must not chain them into one cluster
    w = {"a": window([(100, 40), (125, 4), (150, 35)]), "b": window([(100, 42), (125, 4), (150, 33)])}
    _, b, s, _ = cluster_window(w, A)
    sup = [bb for bb, ss in zip(b, s) if ss]
    assert len(sup) == 2 and sup[0][1] <= 125 <= sup[1][0], (b, s)
    # the same two peaks with a shallow dip (one broad hump) merge into one cluster
    w = {"a": window([(100 + i, 10) for i in range(0, 51, 5)]), "b": window([(100 + i, 10) for i in range(0, 51, 5)])}
    _, b, s, _ = cluster_window(w, A)
    assert len(b) == 1, b
    # support rule: a 3-read peak in one replicate only is a cluster but unsupported
    w = {"a": window([(100, 60), (300, 3)]), "b": window([(100, 60)])}
    A.peak_fraction = 0.0; A.peak_min = 3
    _, b, s, _ = cluster_window(w, A)
    assert len(b) == 2 and s == [True, False], (b, s)
    print("self-test passed")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--family", help="family_<mode>.tsv from step 16")
    ap.add_argument("--sites", nargs="*", default=[], help="replicate=path end-site BED files")
    ap.add_argument("--smooth", type=int, default=11)
    ap.add_argument("--peak-fraction", type=float, default=0.02)
    ap.add_argument("--peak-min", type=int, default=5)
    ap.add_argument("--merge", type=int, default=30)
    ap.add_argument("--extend", type=int, default=15)
    ap.add_argument("--always-merge", type=int, default=15, help="peaks closer than this are always one site")
    ap.add_argument("--valley", type=float, default=0.5, help="dip between two peaks, as a fraction of the lower peak, below which they stay apart")
    ap.add_argument("--support-reads", type=int, default=5)
    ap.add_argument("--support-replicates", type=int, default=2)
    ap.add_argument("--outdir")
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()
    if args.self_test:
        self_test()
        return
    os.makedirs(args.outdir, exist_ok=True)
    fam = pd.read_csv(args.family, sep="\t")
    fam = fam[fam["in_family"]]
    site_dfs = utrlib.load_site_beds(dict(s.split("=", 1) for s in args.sites))
    reps = list(site_dfs)
    by_key = {}
    for rep, df in site_dfs.items():
        for (chrom, strand), grp in df.groupby(["chromosome", "strand"]):
            g = grp.sort_values("start")
            by_key[(rep, chrom, strand)] = (g["start"].to_numpy(), g["score"].to_numpy())

    def per_base(rep, chrom, strand, ws, we):
        arr = np.zeros(we - ws, dtype=int)
        if (rep, chrom, strand) in by_key:
            pos, sc = by_key[(rep, chrom, strand)]
            lo, hi = np.searchsorted(pos, ws, "left"), np.searchsorted(pos, we, "left")
            np.add.at(arr, pos[lo:hi] - ws, sc[lo:hi])
        return arr

    bed, counts_rows, sites_rows, win_rows = [], [], [], []
    for w in fam.itertuples():
        ws, we = int(w.window_start), int(w.window_end)
        per_rep = {r: per_base(r, w.chromosome, w.strand, ws, we) for r in reps}
        pooled, bounds, supported, counts = cluster_window(per_rep, args)
        total = {r: int(per_rep[r].sum()) for r in reps}
        # fold unsupported clusters into other; check the identity per replicate
        other = {r: total[r] - sum(c for c, s in zip(counts[r], supported) if s) for r in reps}
        for r in reps:
            kept = sum(c for c, s in zip(counts[r], supported) if s)
            if kept + other[r] != total[r] or other[r] < 0:
                sys.exit(f"count identity broken for {w.gene} {r}: {kept} + {other[r]} != {total[r]}")
        n_sup = sum(supported)
        testable = n_sup >= 2
        reason = "" if testable else "fewer than two resolved supported clusters"
        win_rows.append(dict(gene=w.gene, chromosome=w.chromosome, strand=w.strand, window_start=ws, window_end=we,
                             pooled_ends=int(pooled.sum()), clusters=len(bounds), supported_clusters=n_sup,
                             testable=testable, exclusion_reason=reason))
        # number clusters in transcript direction
        order = range(len(bounds)) if w.strand == "+" else range(len(bounds) - 1, -1, -1)
        k = 0
        for i in order:
            s, e = ws + bounds[i][0], ws + bounds[i][1]
            k += 1
            dist = (s - w.cds_end) if w.strand == "+" else (w.cds_start - e)
            fid = f"cluster_{k}_{int(dist)}nt"
            bed.append((w.chromosome, s, e, f"{w.gene}|{fid}", int(pooled[bounds[i][0]:bounds[i][1]].sum()), w.strand,
                        "supported" if supported[i] else "unsupported"))
            if testable and supported[i]:
                counts_rows.append(dict(gene_id=w.gene, feature_id=fid, **{r: counts[r][i] for r in reps}))
                sites_rows.append(dict(gene_id=w.gene, feature_id=fid, chromosome=w.chromosome, genomic_start=s,
                                       genomic_end=e, strand=w.strand))
        if testable:
            counts_rows.append(dict(gene_id=w.gene, feature_id="other", **other))
    pd.DataFrame(bed, columns=["chrom", "start", "end", "name", "pooled_ends", "strand", "support"]).to_csv(
        os.path.join(args.outdir, "clusters.bed"), sep="\t", header=False, index=False)
    pd.DataFrame(counts_rows).to_csv(os.path.join(args.outdir, "site_counts.tsv"), sep="\t", index=False)
    pd.DataFrame(sites_rows).to_csv(os.path.join(args.outdir, "sites.tsv"), sep="\t", index=False)
    win = pd.DataFrame(win_rows)
    win.to_csv(os.path.join(args.outdir, "windows.tsv"), sep="\t", index=False)
    pd.DataFrame({"sample_id": reps, "stage": ["BSF" if r.startswith("bs") else "PCF" for r in reps]}).to_csv(
        os.path.join(args.outdir, "samples.tsv"), sep="\t", index=False)
    print(f"{len(win)} windows: {int(win['testable'].sum())} testable (two or more supported clusters), "
          f"{len(bed)} clusters of which {sum(b[6] == 'supported' for b in bed)} supported; "
          f"settings smooth={args.smooth} peak={args.peak_fraction}/{args.peak_min} merge={args.merge} valley={args.valley} extend={args.extend} "
          f"support={args.support_reads}x{args.support_replicates}")
    print(win["supported_clusters"].value_counts().sort_index().to_string())


if __name__ == "__main__":
    main()
