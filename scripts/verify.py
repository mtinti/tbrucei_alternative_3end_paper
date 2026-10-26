"""Check the input data and compare a run with the published tables.

Usage:
    python scripts/verify.py inputs
        Every file in data/ against data/checksums.md5 (for a .gz file, the
        MD5 of the decompressed content) and, when the BAM files are present,
        data/bam/ against data/bam/checksums.md5.
    python scripts/verify.py results DIR
        The tables written by run_analysis.sh into DIR against expected/:
        the family, the supported genes of the window analysis, the eligible
        and supported genes of the gene-assigned analysis, the flagged genes
        with their qualifying cluster and change, and the mapping-rule
        sensitivity table. Expected: identical.

Each command prints OK or FAIL per check and exits with status 1 when a
check fails.
"""
import gzip
import hashlib
import os
import sys

import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
EXPECTED = os.path.join(ROOT, "expected")


def content_md5(path):
    opener = gzip.open if path.endswith(".gz") else open
    digest = hashlib.md5()
    with opener(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def report(ok, text):
    print(f'{"OK  " if ok else "FAIL"} {text}')
    return 0 if ok else 1


def check_list(listing, folder, optional=False):
    failures = 0
    for line in open(listing):
        if not line.strip() or line.startswith("#"):
            continue
        expected, name = line.split()
        path = os.path.join(folder, name)
        if not os.path.exists(path):
            if optional:
                print(f"skip {name} (not present)")
                continue
            failures += report(False, f"{name} missing")
            continue
        failures += report(content_md5(path) == expected, name)
    return failures


def check_inputs():
    failures = check_list(os.path.join(DATA, "checksums.md5"), DATA)
    failures += check_list(os.path.join(DATA, "bam", "checksums.md5"), os.path.join(DATA, "bam"), optional=True)
    return failures


def same_set(a, b, text):
    return report(set(a) == set(b), f"{text}: {len(set(a))} in the run, {len(set(b))} published, {len(set(a) ^ set(b))} differ")


def check_results(out):
    failures = 0
    for mode in ("nanopore", "intergenic"):
        fam_run = pd.read_csv(os.path.join(out, "family", f"family_{mode}.tsv"), sep="\t")
        fam_pub = pd.read_csv(os.path.join(EXPECTED, f"family_{mode}.tsv"), sep="\t")
        failures += same_set(fam_run.loc[fam_run["in_family"], "gene"], fam_pub.loc[fam_pub["in_family"], "gene"], f"{mode}: genes in the family")
        g_run = pd.read_csv(os.path.join(out, f"clusters_{mode}_m30p2", "dexseq_gene.tsv"), sep="\t")
        g_pub = pd.read_csv(os.path.join(EXPECTED, f"dexseq_gene_{mode}.tsv"), sep="\t")
        failures += same_set(g_run.loc[g_run["gene_confirmed"], "gene_id"], g_pub.loc[g_pub["gene_confirmed"], "gene_id"], f"{mode}: supported genes, window analysis")
        s_run = pd.read_csv(os.path.join(out, "docs", f"strong_end_usage_{mode}.tsv"), sep="\t")
        s_pub = pd.read_csv(os.path.join(EXPECTED, f"strong_end_usage_{mode}.tsv"), sep="\t")
        failures += same_set(s_run["gene_id"], s_pub["gene_id"], f"{mode}: eligible genes, gene-assigned analysis")
        failures += same_set(s_run.loc[s_run["gene_supported"], "gene_id"], s_pub.loc[s_pub["gene_supported"], "gene_id"], f"{mode}: supported genes, gene-assigned analysis")
        f_run = s_run[s_run["strong_end_usage"]].set_index("gene_id")
        f_pub = s_pub[s_pub["strong_end_usage"]].set_index("gene_id")
        failures += same_set(f_run.index, f_pub.index, f"{mode}: genes with strong end usage")
        shared = f_run.index.intersection(f_pub.index)
        same_cluster = (f_run.loc[shared, "feature_id"] == f_pub.loc[shared, "feature_id"]).all()
        same_change = (f_run.loc[shared, "assigned_change_points"] - f_pub.loc[shared, "assigned_change_points"]).abs().max() < 0.05
        failures += report(bool(same_cluster and same_change), f"{mode}: qualifying cluster and change of every flagged gene")
    ms_run = pd.read_csv(os.path.join(out, "docs", "mapq_sensitivity_nanopore.tsv"), sep="\t")
    ms_pub = pd.read_csv(os.path.join(EXPECTED, "mapq_sensitivity_nanopore.tsv"), sep="\t")
    cols = ["rule", "eligible_M", "supported", "strong_X", "lost_genes", "gained_genes", "in_every_run"]
    failures += report(ms_run[cols].fillna("").astype(str).equals(ms_pub[cols].fillna("").astype(str)), "mapping-rule sensitivity table")
    return failures


def main():
    if len(sys.argv) >= 2 and sys.argv[1] == "inputs":
        failures = check_inputs()
    elif len(sys.argv) == 3 and sys.argv[1] == "results":
        failures = check_results(sys.argv[2])
    else:
        sys.exit(__doc__)
    print("all checks passed" if failures == 0 else f"{failures} check(s) failed")
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
