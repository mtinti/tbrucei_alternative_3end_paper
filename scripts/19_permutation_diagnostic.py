#!/usr/bin/env python3
"""Run the DEXSeq and stageR step on every three-against-three split of the six samples.

Of the 20 ways to put six samples into two groups of three, 10 are distinct
partitions for a two-sided test. One is the true stage split; the other nine
mix the stages. The call counts of the nine are a diagnostic of how many genes
a wrong labelling produces with these cluster definitions. They are not an
estimate of the false discoveries in the true comparison, because the mixed
splits still contain real stage effects. Cluster definitions are fixed; the
dispersions and fits are redone for every split.

Output: permutation_calls.tsv in --folder (one row per partition: the samples
labelled as group A, whether it is the true split, genes and features
confirmed).
"""
import argparse
import itertools
import os
import shutil
import subprocess
import pandas as pd


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--folder", required=True, help="step 17 folder with site_counts.tsv and samples.tsv")
    ap.add_argument("--script", default=os.path.join(os.path.dirname(__file__), "18_dexseq_stager.R"))
    ap.add_argument("--fdr", default="0.05")
    ap.add_argument("--threads", default="4")
    args = ap.parse_args()
    samples = pd.read_csv(os.path.join(args.folder, "samples.tsv"), sep="\t")
    reps = samples["sample_id"].tolist()
    true_a = set(samples.loc[samples["stage"] == "BSF", "sample_id"])
    rows = []
    seen = set()
    for combo in itertools.combinations(reps, 3):
        a = frozenset(combo)
        b = frozenset(reps) - a
        if b in seen:
            continue
        seen.add(a)
        name = "perm_" + "_".join(sorted(a))
        sub = os.path.join(args.folder, "permutations", name)
        os.makedirs(sub, exist_ok=True)
        shutil.copy(os.path.join(args.folder, "site_counts.tsv"), sub)
        pd.DataFrame({"sample_id": reps, "stage": ["BSF" if r in a else "PCF" for r in reps]}).to_csv(
            os.path.join(sub, "samples.tsv"), sep="\t", index=False)
        r = subprocess.run(["Rscript", args.script, sub, args.fdr, args.threads], capture_output=True, text=True)
        if r.returncode != 0:
            print(r.stderr[-2000:])
            raise SystemExit(f"step 18 failed for {name}")
        g = pd.read_csv(os.path.join(sub, "dexseq_gene.tsv"), sep="\t")
        s = pd.read_csv(os.path.join(sub, "dexseq_site.tsv"), sep="\t")
        rows.append(dict(partition=name, group_A=",".join(sorted(a)), true_split=(a == true_a),
                         genes_screened=int(g["gene_qvalue"].notna().sum()), genes_confirmed=int(g["gene_confirmed"].sum()),
                         features_confirmed=int(s["site_confirmed"].sum())))
        print(rows[-1])
    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(args.folder, "permutation_calls.tsv"), sep="\t", index=False)
    mixed = out[~out["true_split"]]
    print(f"true split: {int(out.loc[out['true_split'], 'genes_confirmed'].iloc[0])} genes; "
          f"mixed splits: median {mixed['genes_confirmed'].median():.0f}, max {mixed['genes_confirmed'].max()} genes")


if __name__ == "__main__":
    main()
