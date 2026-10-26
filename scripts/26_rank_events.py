#!/usr/bin/env python3
"""Rank the supported genes for follow-up (docs/methods.md, section 6).

Every gene in the table is already supported (stageR, gene-level OFDR 5
percent), so significance has done its job as a gate. The ranking then asks
how large and how consistent the change is, and how well it is understood.
Per gene, the strongest supported cluster (smallest p value) is the ranked
site. Criteria, all from saved tables:

- change: usage of the site in PCF minus BSF, in percentage points (mean of
  the raw replicate fractions);
- replicates separated: every replicate of the stage with higher usage has a
  higher raw fraction than every replicate of the other stage;
- supported in all analyses: both other cluster settings and DRIMSeq also
  support the gene (step 20);
- label of the site (step 21).

Tiers:
  1  separated, |change| >= 30, supported in all analyses, label gene-connected or mixed
  2  separated, |change| >= 20, supported in all analyses
  3  |change| >= 20 or separated
  4  the rest (supported, but a small or inconsistent shift)
Within a tier genes are ordered by |change|, the size of the usage change of
the strongest site. The ranking covers the supported genes that have a
supported cluster; genes whose only confirmed feature is "other" are not in
it. Support in all analyses is gene-level: the same site need not be the
supported one under every setting.
"""
import argparse
import pandas as pd


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--events", required=True)
    ap.add_argument("--dexseq-site", required=True)
    ap.add_argument("--support", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    ev = pd.read_csv(args.events, sep="\t")
    site = pd.read_csv(args.dexseq_site, sep="\t")
    sup = pd.read_csv(args.support, sep="\t").set_index("gene_id")
    frac_cols = [c for c in site.columns if c.startswith("frac_")]
    bs = [c for c in frac_cols if "bs_" in c]; pc = [c for c in frac_cols if "pc_" in c]
    best = ev.sort_values("site_pvalue").drop_duplicates("gene_id").merge(site[["gene_id", "feature_id"] + frac_cols], on=["gene_id", "feature_id"])
    up_pcf = best["usage_change_points"] > 0
    sep = (best[pc].min(axis=1) > best[bs].max(axis=1)) & up_pcf | (best[bs].min(axis=1) > best[pc].max(axis=1)) & ~up_pcf
    best["replicates_separated"] = sep
    best["supported_in_all_analyses"] = best["gene_id"].map(sup["supported_in_all_analyses"]).fillna(False).astype(bool)
    a = best["usage_change_points"].abs()
    good_label = best["label"].isin(["gene-connected", "mixed"])
    tier = pd.Series(4, index=best.index)
    tier[(a >= 20) | sep] = 3
    tier[sep & (a >= 20) & best["supported_in_all_analyses"]] = 2
    tier[sep & (a >= 30) & best["supported_in_all_analyses"] & good_label] = 1
    best["tier"] = tier
    best["abs_change_points"] = a.round(1)
    out = best.sort_values(["tier", "abs_change_points"], ascending=[True, False])
    cols = ["tier", "abs_change_points", "gene_id", "feature_id", "site_location", "usage_BSF", "usage_PCF", "usage_change_points",
            "replicates_separated", "supported_in_all_analyses", "label", "reads_total", "informative_reads", "site_pvalue"] + frac_cols
    out[cols].to_csv(args.out, sep="\t", index=False)
    print(f"wrote {args.out}: " + ", ".join(f"tier {t}: {n}" for t, n in out["tier"].value_counts().sort_index().items()))
    print(out[cols[:11]].head(15).to_string(index=False))


if __name__ == "__main__":
    main()
