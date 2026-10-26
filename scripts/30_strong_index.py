#!/usr/bin/env python3
"""Write an HTML index of the genes with strong evidence of alternative 3' end usage (docs/methods.md, section 7).

One section per window type, the primary (nanopore) first. The table lists
every flagged gene with its qualifying cluster, the gene-assigned usage per
stage, the change, the assigned ends in the cluster, the window-analysis
change of the same cluster and its label. Below the table each gene shows the
zoomed figure (step 27) and the overview (step 11).
"""
import argparse
import os
import pandas as pd


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--strong", nargs="+", required=True, help="label=strong_end_usage_<mode>.tsv=zoom_folder=overview_folder")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    base = os.path.dirname(os.path.abspath(args.out))
    html = ["<!doctype html><html><head><meta charset='utf-8'><title>Strong alternative 3' end usage</title>",
            "<style>body{font-family:sans-serif;margin:20px} table{border-collapse:collapse} td,th{border:1px solid #ccc;padding:3px 6px;font-size:13px}"
            " img{max-width:100%} .g{margin-top:30px;border-top:2px solid #ddd;padding-top:10px}</style></head><body>",
            "<h1>Genes with strong evidence of a different 3' end between BSF and PCF</h1>",
            "<p>A gene is here when, among the reads assigned to it, at least one end cluster changes its share by 30 points or "
            "more between stages, every replicate of one stage lies above every replicate of the other, the same cluster changes in "
            "the same direction under both other cluster settings, and the cluster holds 20 or more assigned reads; the gene and the "
            "cluster pass stageR at a nominal gene-level overall false discovery rate of 5 percent (docs/methods.md, section 7). A read "
            "is assigned to the gene when at least 50 of its aligned bases lie inside the gene's CDS, its mapping quality is 20 or "
            "more, no secondary alignment of it scores at least as high as the primary, and it covers no other CDS on the same "
            "strand. The rule resolves the read's origin as far as the alignments allow; it is not an experimental proof of origin. "
            "Usage is the share of the gene's assigned ends. The window change is the same cluster's change in the window analysis, "
            "whose denominator is every end in the window. In the zoomed figure, reads of uncertain origin under this rule (no CDS "
            "overlap, low mapping quality, equal placement elsewhere, another CDS covered) are grey. Every figure reads 5' to 3' "
            "from left to right; for minus-strand genes the x axis is reversed and the axis label says so.</p>"]
    for spec in args.strong:
        label, path, zoom, overview = spec.split("=")
        se = pd.read_csv(path, sep="\t")
        f = se[se["strong_end_usage"]]
        html.append(f"<h2>{label}: {len(f)} of {len(se)} eligible genes</h2>")
        html.append("<table><tr><th>Gene</th><th>Cluster</th><th>Location</th><th>Assigned usage BSF</th><th>PCF</th><th>Change (points)</th>"
                    "<th>Assigned ends in cluster (share of all ends there)</th><th>Window change</th><th>Label (window analysis)</th></tr>")
        for r in f.itertuples():
            wc = f"{float(r.window_change_points):+.0f}" if str(r.window_change_points) not in ("", "nan") else "not tested"
            lab = r.label_window_analysis if isinstance(r.label_window_analysis, str) else "none"
            html.append(f"<tr><td><a href='#{label}-{r.gene_id}'>{r.gene_id}</a></td><td>{r.feature_id}</td><td>{r.site_location}</td>"
                        f"<td>{r.assigned_usage_BSF:.2f}</td><td>{r.assigned_usage_PCF:.2f}</td><td>{r.assigned_change_points:+.0f}</td>"
                        f"<td>{r.assigned_ends_in_cluster} ({r.assigned_share})</td><td>{wc}</td><td>{lab}</td></tr>")
        html.append("</table>")
        for r in f.itertuples():
            parts = [f"<div class='g' id='{label}-{r.gene_id}'><h3>{r.gene_id}</h3><p>{r.feature_id}: assigned usage BSF {r.assigned_usage_BSF:.2f}, "
                     f"PCF {r.assigned_usage_PCF:.2f} ({r.assigned_change_points:+.0f} points)</p>"]
            for title, folder in (("Tested region, tested reads only.", zoom), ("Overview of the locus.", overview)):
                img = os.path.relpath(os.path.join(folder, f"{r.gene_id}.png"), base)
                if os.path.exists(os.path.join(base, img)):
                    parts.append(f"<p><b>{title}</b></p><img src='{img}'>")
            html.append("".join(parts) + "</div>")
    html.append("</body></html>")
    open(args.out, "w").write("\n".join(html))
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
