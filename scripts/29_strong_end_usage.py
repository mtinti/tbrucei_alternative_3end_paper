#!/usr/bin/env python3
"""Flag the genes with strong evidence of alternative end usage (docs/methods.md, section 7).

One row per gene eligible for the gene-assigned analysis (step 28 and step 18
on its counts). A gene is flagged strong_end_usage when at least one of its
supported clusters (other excluded) satisfies every criterion of the plan:
supported in the primary gene-assigned analysis; gene-assigned usage change
of --min-change points or more; replicates separated; an overlapping
supported cluster changing in the same direction in both other cluster
settings of the gene-assigned analysis, with the gene supported there; at
least --min-cluster-ends assigned ends in the cluster. Fractions, stage
means, the change and the separation are computed from the count matrix
(site_counts.tsv), not from the rounded values of the model table, so the
threshold applies to the exact value; rounding is for display. The rule runs over
every supported cluster of the gene and the gene is counted once. The row
shows the qualifying cluster or, when none qualifies, the supported cluster
with the largest change, and the reason it fails.
"""
import argparse
import os
import pandas as pd


def load(folder):
    g = pd.read_csv(os.path.join(folder, "dexseq_gene.tsv"), sep="\t")
    s = pd.read_csv(os.path.join(folder, "dexseq_site.tsv"), sep="\t").merge(
        pd.read_csv(os.path.join(folder, "sites.tsv"), sep="\t"), on=["gene_id", "feature_id"], how="left")
    return g, s


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--primary", required=True, help="step 28 folder of the primary setting, with step 18 results")
    ap.add_argument("--settings", nargs=2, required=True, help="name=folder for the two other cluster settings")
    ap.add_argument("--window-dexseq", required=True, help="step 18 folder of the window analysis, primary setting")
    ap.add_argument("--events", help="events_<mode>.tsv from step 21, for the labels")
    ap.add_argument("--min-change", type=float, default=30.0)
    ap.add_argument("--min-cluster-ends", type=int, default=20)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    g, s = load(args.primary)
    win = pd.read_csv(os.path.join(args.primary, "windows.tsv"), sep="\t")
    cm = pd.read_csv(os.path.join(args.primary, "site_counts.tsv"), sep="\t")
    rep_cols = [c for c in cm.columns if c not in ("gene_id", "feature_id")]
    tot = cm.groupby("gene_id")[rep_cols].transform("sum").clip(lower=1)
    exact = cm[["gene_id", "feature_id"]].copy()
    for c in rep_cols:
        exact[c] = cm[c] / tot[c]
    exact = exact.set_index(["gene_id", "feature_id"])
    assign = pd.read_csv(os.path.join(args.primary, "assignment.tsv"), sep="\t")
    reps = [c[len("assigned_"):] for c in assign.columns if c.startswith("assigned_")]
    bs = [r for r in reps if r.startswith("bs")]; pc = [r for r in reps if r.startswith("pc")]
    assign = assign.set_index(["gene_id", "feature_id"])
    wg, ws = load(args.window_dexseq)
    wsup = ws.set_index(["gene_id", "feature_id"])
    wgene = wg.set_index("gene_id")["gene_confirmed"]
    others = {}
    for spec in args.settings:
        name, folder = spec.split("=", 1)
        og, os_ = load(folder)
        others[name] = (set(og.loc[og["gene_confirmed"], "gene_id"]), os_[os_["site_confirmed"] & (os_["feature_id"] != "other")])
    labels = {}
    if args.events:
        ev = pd.read_csv(args.events, sep="\t")
        labels = dict(zip(zip(ev["gene_id"], ev["feature_id"]), ev["label"]))

    gene_q = g.set_index("gene_id")
    rows = []
    for w in win[win["testable"]].itertuples():
        gid = w.gene
        gsup = bool(gene_q.loc[gid, "gene_confirmed"]) if gid in gene_q.index else False
        cands = s[(s["gene_id"] == gid) & (s["feature_id"] != "other")].copy()
        ex = exact.loc[gid]
        cands["usage_BSF"] = [ex.loc[f, bs].mean() for f in cands["feature_id"]]
        cands["usage_PCF"] = [ex.loc[f, pc].mean() for f in cands["feature_id"]]
        cands["usage_change_points"] = 100 * (cands["usage_PCF"] - cands["usage_BSF"])
        for r in reps:
            cands[f"frac_{r}"] = [ex.loc[f, r] for f in cands["feature_id"]]
        cands["abs_change"] = cands["usage_change_points"].abs()
        cands = cands.sort_values(["site_confirmed", "abs_change"], ascending=[False, False])
        best, best_reason = None, ""
        for c in cands.itertuples():
            reasons = []
            if not gsup:
                reasons.append("gene not supported in the gene-assigned analysis")
            elif not c.site_confirmed:
                reasons.append("cluster not supported")
            if abs(c.usage_change_points) < args.min_change:
                reasons.append(f"change below {args.min_change:.0f} points")
            fr = cands.loc[c.Index, [f"frac_{r}" for r in reps]]
            up = c.usage_change_points > 0
            sep = (fr[[f"frac_{r}" for r in pc]].min() > fr[[f"frac_{r}" for r in bs]].max()) if up else \
                  (fr[[f"frac_{r}" for r in bs]].min() > fr[[f"frac_{r}" for r in pc]].max())
            if not sep:
                reasons.append("replicates not separated")
            sens = {}
            for name, (ogenes, osites) in others.items():
                if gid not in ogenes:
                    sens[name] = "gene not supported"
                else:
                    o = osites[(osites["gene_id"] == gid) & (osites["genomic_start"] < c.genomic_end) & (osites["genomic_end"] > c.genomic_start)]
                    if o.empty:
                        sens[name] = "no overlapping supported cluster"
                    else:
                        sens[name] = "same direction" if (o["usage_change_points"] * c.usage_change_points > 0).any() else "opposite direction"
                if sens[name] != "same direction":
                    reasons.append(f"{name}: {sens[name]}")
            a = assign.loc[(gid, c.feature_id)]
            n_assigned = int(sum(a[f"assigned_{r}"] for r in reps)); n_unc = int(sum(a[f"uncertain_{r}"] for r in reps))
            if n_assigned < args.min_cluster_ends:
                reasons.append(f"fewer than {args.min_cluster_ends} assigned ends in the cluster")
            key = (gid, c.feature_id)
            wchange = wsup.loc[key, "usage_change_points"] if key in wsup.index else float("nan")
            wconf = bool(wsup.loc[key, "site_confirmed"]) if key in wsup.index else False
            row = dict(gene_id=gid, strong_end_usage=not reasons, feature_id=c.feature_id,
                       site_location=f"{c.chromosome}:{c.genomic_start}-{c.genomic_end}({c.strand})",
                       assigned_usage_BSF=round(c.usage_BSF, 4), assigned_usage_PCF=round(c.usage_PCF, 4),
                       assigned_change_points=round(c.usage_change_points, 1), replicates_separated=bool(sep),
                       cluster_supported=bool(c.site_confirmed), site_pvalue=c.pvalue, site_padj=c.site_padj,
                       **{f"{n}_same_site_direction": v for n, v in sens.items()},
                       assigned_ends_in_cluster=n_assigned, uncertain_ends_in_cluster=n_unc,
                       assigned_share=round(n_assigned / (n_assigned + n_unc), 3) if n_assigned + n_unc else "",
                       label_window_analysis=labels.get(key, ""),
                       window_change_points=round(wchange, 1) if wchange == wchange else "",
                       window_cluster_supported=wconf, window_gene_supported=bool(wgene.get(gid, False)),
                       reason_not_strong="; ".join(reasons),
                       **{f"frac_{r}": round(fr[f"frac_{r}"], 4) for r in reps})
            if not reasons:
                best, best_reason = row, ""
                break
            if best is None:
                best = row
        if best is None:
            best = dict(gene_id=gid, strong_end_usage=False, reason_not_strong="no supported cluster besides other")
        best.update(gene_supported=gsup, gene_qvalue=gene_q.loc[gid, "gene_qvalue"] if gid in gene_q.index else "",
                    supported_clusters=int(w.supported_clusters), assigned_ends_window=int(w.pooled_ends),
                    **{f"assigned_{r}": int(getattr(w, f"assigned_{r}")) for r in reps})
        rows.append(best)
    out = pd.DataFrame(rows)
    front = ["gene_id", "strong_end_usage", "gene_supported", "gene_qvalue", "supported_clusters", "assigned_ends_window"] + [f"assigned_{r}" for r in reps]
    out = out[front + [c for c in out.columns if c not in front]]
    out = out.sort_values(["strong_end_usage", "assigned_change_points"], ascending=[False, False], key=lambda c: c.abs() if c.name == "assigned_change_points" else c)
    out.to_csv(args.out, sep="\t", index=False)
    m, x = len(out), int(out["strong_end_usage"].sum())
    print(f"wrote {args.out}: M = {m} eligible genes, {int(out['gene_supported'].sum())} supported, X = {x} with strong end usage")
    print(out.loc[~out["strong_end_usage"] & out["gene_supported"], "reason_not_strong"].value_counts().head(12).to_string())


if __name__ == "__main__":
    main()
