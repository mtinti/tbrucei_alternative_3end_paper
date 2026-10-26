#!/usr/bin/env bash
# Build the data archive for Zenodo: the files too large for GitHub.
#
#   bash scripts/make_zenodo_archive.sh results zenodo
#
# Writes into the second argument:
#   bam/                      the four BAM files with their indexes and checksums (data/bam/)
#   results_full.tar.gz       every table of one complete run (the first argument), including the
#                             sensitivity fits and the read placements
#   figures_all.tar.gz        the figures of every supported gene and both HTML indexes
#   checksums.md5             MD5 of every file in the archive
set -euo pipefail
OUT="${1:-results}"; Z="${2:-zenodo}"
mkdir -p "$Z/bam"
cp data/bam/*.bam data/bam/*.bai data/bam/checksums.md5 "$Z/bam/"
tar -czf "$Z/results_full.tar.gz" --exclude='docs/figures' --exclude='inputs' -C "$(dirname "$OUT")" "$(basename "$OUT")"
if [[ -d "$OUT/docs/figures" ]]; then tar -czf "$Z/figures_all.tar.gz" -C "$OUT/docs" figures; fi
(cd "$Z" && find . -type f ! -name checksums.md5 -exec md5sum {} \; > checksums.md5 2>/dev/null || find . -type f ! -name checksums.md5 -exec md5 -r {} \; > checksums.md5)
du -sh "$Z"/*
