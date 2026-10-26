#!/usr/bin/env python3
"""Summarise the sensitivity of the Part 7 result to the mapping-confidence requirement (docs/methods.md, section 7).

One row per run of the gene-assigned chain (step 28, step 18, step 29) under
a different placement rule: the primary rule (mapping quality 20, secondary
rule on) and the variants fixed in the plan. Columns: the rule, the assigned
and uncertain reads of the family by reason (from windows.tsv), M (eligible
genes), the supported genes, X (flagged genes), and how the flagged set
differs from the primary one (genes lost and gained, with their names). A
last column, in_every_run, gives the number of primary genes flagged under
every run, and the primary row names the genes that are not.
"""
import argparse
import pandas as pd


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--runs", nargs="+", required=True,
                    help="label=min_mapq=secondary_rule=strong_end_usage.tsv=windows.tsv; the first is the primary rule")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    rows, primary = [], None
    for spec in args.runs:
        label, q, sec, strong_path, win_path = spec.split("=")
        se = pd.read_csv(strong_path, sep="\t")
        w = pd.read_csv(win_path, sep="\t")
        flagged = set(se.loc[se["strong_end_usage"], "gene_id"])
        if primary is None:
            primary = flagged
        unc = {c: int(w[c].sum()) for c in w.columns if c.startswith("uncertain_")}
        rows.append(dict(rule=label, min_mapq=int(q), secondary_rule=sec, assigned_reads=int(w["pooled_ends"].sum()), **unc,
                         eligible_M=len(se), supported=int(se["gene_supported"].sum()), strong_X=len(flagged),
                         in_primary_set=len(flagged & primary), lost_vs_primary=len(primary - flagged), gained_vs_primary=len(flagged - primary),
                         lost_genes=";".join(sorted(primary - flagged)), gained_genes=";".join(sorted(flagged - primary))))
    out = pd.DataFrame(rows)
    sets = [set(r["lost_genes"].split(";")) - {""} for r in rows]
    not_everywhere = set().union(*sets)
    out["in_every_run"] = len(primary - not_everywhere)
    out["primary_genes_not_in_every_run"] = ""
    out.loc[0, "primary_genes_not_in_every_run"] = ";".join(sorted(not_everywhere))
    out.to_csv(args.out, sep="\t", index=False)
    print(f"{len(primary - not_everywhere)} of {len(primary)} primary genes are flagged under every run; not everywhere: {sorted(not_everywhere)}")
    print(f"wrote {args.out}")
    print(out.drop(columns=["lost_genes", "gained_genes", "primary_genes_not_in_every_run"]).to_string(index=False))


if __name__ == "__main__":
    main()
