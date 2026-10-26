"""Read-to-gene assignment rule shared by steps 27 and 28 (docs/methods.md, section 7; corrected after the eighth and the ninth review).

A poly(A) read with the step 02 eligibility is assigned to a gene when

- at least `min_overlap` of its aligned bases fall inside the gene's CDS
  (aligned bases, so a deletion or a skipped reference stretch does not
  count; `pysam.AlignedSegment.get_overlap`);
- its placement is resolved: mapping quality at least `min_mapq` (20 by
  default, the mapper's estimate that the placement is wrong with probability
  1 percent or less), and no secondary alignment of the read, anywhere in the
  genome, with an alignment score equal to or higher than the primary's
  (minimap2 chooses the primary by chain score, so such a secondary is a
  placement at least as good; the scores come from step 31, `Placements`);
- no other CDS on the same strand receives `min_overlap` aligned bases from
  it, so the read does not carry two genes (a read-through or dicistronic
  read is not attributable to one gene).

Any other read is of uncertain origin, with the reason recorded: 'no CDS
overlap', 'low mapping quality', 'equal placement elsewhere', 'covers another
CDS'. The rule does not look at the stage.
"""
import bisect
import gzip
from collections import defaultdict

import pandas as pd

REASONS = ("no CDS overlap", "low mapping quality", "equal placement elsewhere", "covers another CDS")


class CdsIndex:
    """Same-strand CDS intervals per chromosome, for the other-gene check."""

    def __init__(self, cds_bed):
        cds = pd.read_csv(cds_bed, sep="\t", header=None, usecols=range(6), names=["chrom", "start", "end", "gene", "score", "strand"])
        self.by_key = defaultdict(list)
        for r in cds.itertuples():
            self.by_key[(r.chrom, r.strand)].append((int(r.start), int(r.end), r.gene))
        self.starts = {}
        for k, v in self.by_key.items():
            v.sort()
            self.by_key[k] = v
            self.starts[k] = [x[0] for x in v]
        self.max_len = {k: max((e - s for s, e, _ in v), default=0) for k, v in self.by_key.items()}

    def overlapping(self, chrom, strand, start, end):
        """CDS intervals on (chrom, strand) that overlap [start, end)."""
        k = (chrom, strand)
        if k not in self.by_key:
            return []
        v, st = self.by_key[k], self.starts[k]
        lo = bisect.bisect_left(st, start - self.max_len[k])
        out = []
        for s, e, g in v[lo:]:
            if s >= end:
                break
            if e > start:
                out.append((s, e, g))
        return out


class Placements:
    """Best secondary alignment score per read, from step 31 (read_placements.tsv.gz).

    `equal_elsewhere(read)` is True when a secondary alignment of the read
    scores at least as high as its primary. Reads absent from the table have
    no secondary record."""

    def __init__(self, path):
        self.best = {}
        with gzip.open(path, "rt") as f:
            head = f.readline().rstrip("\n").split("\t")
            i_read, i_best = head.index("read"), head.index("best_secondary_as")
            for line in f:
                x = line.rstrip("\n").split("\t")
                b = int(x[i_best])
                if b >= 0:
                    self.best[x[i_read]] = b

    def equal_elsewhere(self, read):
        b = self.best.get(read.query_name)
        return b is not None and b >= int(read.get_tag("AS"))


def assignment(read, gene, cds_start, cds_end, strand, index, min_overlap=50, min_mapq=20, placements=None):
    """Return (assigned, reason). reason is '' when assigned, else one of REASONS.
    With placements None the secondary-alignment rule is not applied."""
    if read.get_overlap(cds_start, cds_end) < min_overlap:
        return False, "no CDS overlap"
    if read.mapping_quality < min_mapq:
        return False, "low mapping quality"
    if placements is not None and placements.equal_elsewhere(read):
        return False, "equal placement elsewhere"
    for s, e, g in index.overlapping(read.reference_name, strand, read.reference_start, read.reference_end):
        if g != gene and read.get_overlap(s, e) >= min_overlap:
            return False, "covers another CDS"
    return True, ""
