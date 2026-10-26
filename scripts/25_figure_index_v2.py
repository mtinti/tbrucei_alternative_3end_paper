#!/usr/bin/env python3
"""Write an HTML index of the browser figures of the supported genes of the revised analysis.

One section per window type. Each gene appears once with every supported
cluster and its event label (step 21), ordered by the gene's most
significant cluster. The figure comes from step 11 drawn on the step 17
clusters with the DEXSeq usage.
"""
import argparse
import os
import pandas as pd


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--events", nargs="+", required=True, help="label=events_<mode>.tsv=figure_folder")
    ap.add_argument("--ranked", nargs="*", default=[], help="label=ranked_<mode>.tsv from step 26; orders the genes by tier and size of change")
    ap.add_argument("--zoom", nargs="*", default=[], help="label=folder of the zoomed figures (step 27), shown above the overview")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    base = os.path.dirname(os.path.abspath(args.out))
    html = ["<!doctype html><html><head><meta charset='utf-8'><title>Supported genes, revised analysis</title>",
            "<style>body{font-family:sans-serif;margin:20px} table{border-collapse:collapse} td,th{border:1px solid #ccc;padding:3px 6px;font-size:13px}"
            " img{max-width:100%} .g{margin-top:30px;border-top:2px solid #ddd;padding-top:10px}</style></head><body>",
            "<h1>Genes with supported stage-specific 3' end usage, revised analysis</h1>",
            "<p>Supported: passes stageR at a nominal gene-level overall false discovery rate of 5 percent (DEXSeq on end "
            "clusters defined from the pooled ends of all replicates). Labels belong to clusters, not genes: gene-connected "
            "(the gene's own transcript ends there), transcript-associated (a separate downstream transcript ends there), "
            "mixed (both), unresolved (the reads do not settle it). Usage is the mean raw fraction of the window's ends per stage. "
            "Every figure reads 5' to 3' from left to right; for minus-strand genes the x axis is reversed and the axis label says so.</p>"]
    ranked = dict(s.split("=", 1) for s in args.ranked)
    zoom = dict(s.split("=", 1) for s in args.zoom)
    for spec in args.events:
        label, path, folder = spec.split("=")
        ev = pd.read_csv(path, sep="\t")
        tier = {}
        if label in ranked:
            rk = pd.read_csv(ranked[label], sep="\t")
            order = rk["gene_id"].tolist()
            tier = dict(zip(rk["gene_id"], rk["tier"]))
            note = ("<p>Ordered by priority tier, then by the size of the usage change of the strongest site (step 26). Tier 1: "
                    "replicates separated, change of 30 points or more, supported in all analyses, site gene-connected or mixed. "
                    "Tier 2: separated, 20 points or more, supported in all analyses. Tier 3: 20 points or more, or separated. "
                    "Tier 4: the rest.</p>")
        else:
            order = ev.sort_values("site_pvalue").drop_duplicates("gene_id")["gene_id"].tolist()
            note = "<p>Ordered by the p value of the strongest site.</p>"
        html.append(f"<h2>{label} ({len(order)} genes, {len(ev)} supported clusters)</h2>" + note)
        html.append("<table><tr><th>Tier</th><th>Gene</th><th>Cluster</th><th>Location</th><th>Usage BSF</th><th>Usage PCF</th><th>Change (points)</th>"
                    "<th>Reads</th><th>Label</th><th>Transcript evidence</th></tr>")
        for g in order:
            for r in ev[ev["gene_id"] == g].sort_values("site_pvalue").itertuples():
                evidence = "; ".join(x for x in [
                    (f"catalogue {r.catalogue_record} ({r.catalogue_abundance})" if isinstance(r.catalogue_record, str) and r.catalogue_record else ""),
                    str(r.annotated_ncRNA) if isinstance(r.annotated_ncRNA, str) else "",
                    "acceptor cluster downstream of the CDS" if r.acceptor_cluster else ""] if x)
                html.append(f"<tr><td>{tier.get(g, '')}</td><td><a href='#{label}-{g}'>{g}</a></td><td>{r.feature_id}</td><td>{r.site_location}</td>"
                            f"<td>{r.usage_BSF:.2f}</td><td>{r.usage_PCF:.2f}</td><td>{r.usage_change_points:+.0f}</td>"
                            f"<td>{r.reads_total}</td><td>{r.label}</td><td>{evidence}</td></tr>")
        html.append("</table>")
        for g in order:
            img = os.path.relpath(os.path.join(folder, f"{g}.png"), base)
            if not os.path.exists(os.path.join(base, img)):
                continue
            labels = ", ".join(f"{r.feature_id}: {r.label}" for r in ev[ev["gene_id"] == g].sort_values("site_pvalue").itertuples())
            zimg = ""
            if label in zoom:
                zpath = os.path.relpath(os.path.join(zoom[label], f"{g}.png"), base)
                if os.path.exists(os.path.join(base, zpath)):
                    zimg = f"<p><b>Tested region, tested reads only.</b></p><img src='{zpath}'><p><b>Overview of the locus.</b></p>"
            html.append(f"<div class='g' id='{label}-{g}'><h3>{g}" + (f" <small>(tier {tier[g]})</small>" if g in tier else "") + f"</h3><p>{labels}</p>{zimg}<img src='{img}'></div>")
    html.append("</body></html>")
    open(args.out, "w").write("\n".join(html))
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
