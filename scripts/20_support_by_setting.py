#!/usr/bin/env python3
"""Join the primary analysis with the sensitivity settings and the DRIMSeq fit.

One row per gene of the primary family. Columns per setting: testable and
supported; for the DRIMSeq sensitivity fit on the primary clusters: tested and
supported. Two summary columns: supported_in_all_analyses (supported by the
primary, by both other settings and by DRIMSeq, all of which tested the gene)
and supported_wherever_tested (supported by every analysis that could test
the gene, with missing_comparisons counting the analyses that could not).
Neither carries an error rate of its own. For genes supported in the primary and in another setting, the table
says whether a supported cluster of the other setting overlaps the primary's
most significant supported cluster and changes in the same direction.
"""
import argparse
import os
import pandas as pd


def load(folder):
    g = pd.read_csv(os.path.join(folder, "dexseq_gene.tsv"), sep="\t")
    s = pd.read_csv(os.path.join(folder, "dexseq_site.tsv"), sep="\t").merge(
        pd.read_csv(os.path.join(folder, "sites.tsv"), sep="\t"), on=["gene_id", "feature_id"], how="left")
    w = pd.read_csv(os.path.join(folder, "windows.tsv"), sep="\t")
    return g, s, w


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--primary", required=True, help="name=folder")
    ap.add_argument("--settings", nargs="*", default=[], help="name=folder of the other cluster settings")
    ap.add_argument("--drimseq", help="folder with drimseq_gene.tsv from step 09 on the primary clusters")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    pname, pfolder = args.primary.split("=", 1)
    pg, ps, pw = load(pfolder)
    out = pw[["gene", "testable"]].rename(columns={"gene": "gene_id", "testable": f"{pname}_testable"}).copy()
    out[f"{pname}_screened"] = out["gene_id"].isin(pg.loc[pg["gene_qvalue"].notna(), "gene_id"])
    out[f"{pname}_supported"] = out["gene_id"].isin(pg.loc[pg["gene_confirmed"], "gene_id"])
    best = ps[ps["site_confirmed"] & (ps["feature_id"] != "other")].sort_values("pvalue").drop_duplicates("gene_id").set_index("gene_id")
    out["primary_site"] = out["gene_id"].map(best["feature_id"])
    out["primary_site_change_points"] = out["gene_id"].map(best["usage_change_points"])
    stable = out[f"{pname}_supported"].copy()
    strict = out[f"{pname}_supported"].copy()
    missing = pd.Series(0, index=out.index)
    for spec in args.settings:
        name, folder = spec.split("=", 1)
        g, s, w = load(folder)
        out[f"{name}_testable"] = out["gene_id"].map(w.set_index("gene")["testable"]).fillna(False).astype(bool)
        out[f"{name}_supported"] = out["gene_id"].isin(g.loc[g["gene_confirmed"], "gene_id"])
        same = []
        sup = s[s["site_confirmed"] & (s["feature_id"] != "other")]
        for gid, row in zip(out["gene_id"], out.itertuples()):
            b = best.loc[gid] if gid in best.index else None
            if b is None or not out.loc[out["gene_id"] == gid, f"{name}_supported"].iloc[0]:
                same.append("")
                continue
            o = sup[(sup["gene_id"] == gid) & (sup["genomic_start"] < b["genomic_end"]) & (sup["genomic_end"] > b["genomic_start"])]
            if o.empty:
                same.append("no overlapping supported cluster")
            else:
                same.append("same direction" if (o["usage_change_points"] * b["usage_change_points"] > 0).any() else "opposite direction")
        out[f"{name}_same_site_direction"] = same
        stable &= ~out[f"{name}_testable"] | out[f"{name}_supported"]
        strict &= out[f"{name}_supported"]
        missing += (~out[f"{name}_testable"]).astype(int)
    if args.drimseq:
        d = pd.read_csv(os.path.join(args.drimseq, "drimseq_gene.tsv"), sep="\t")
        out["drimseq_tested"] = out["gene_id"].isin(d["gene_id"])
        out["drimseq_supported"] = out["gene_id"].isin(d.loc[d["gene_confirmed"], "gene_id"])
        stable &= ~out["drimseq_tested"] | out["drimseq_supported"]
        strict &= out["drimseq_supported"]
        missing += (~out["drimseq_tested"]).astype(int)
    out["supported_in_all_analyses"] = strict
    out["supported_wherever_tested"] = stable
    out["missing_comparisons"] = missing.where(out[f"{pname}_supported"], 0)
    out.to_csv(args.out, sep="\t", index=False)
    n = int(out[f"{pname}_supported"].sum())
    print(f"primary supported {n}; supported in all analyses {int(strict.sum())}; supported wherever tested {int(stable.sum())}; "
          + "; ".join(f"{c.replace('_supported', '')} supported {int(out[c].sum())}" for c in out.columns
                      if c.endswith("_supported") and not c.startswith(pname)))


if __name__ == "__main__":
    main()
