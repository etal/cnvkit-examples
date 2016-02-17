#!/usr/bin/env python

"""Match up and aggregate gene coverages between 2 sets of samples."""
from __future__ import division, print_function

import sys

import numpy as np

import cnvlib
from cnvlib import params
from cnvlib.rary import RegionArray as RA


# --- by targeted interval ---

def read_paired_genes(cbs1, cbs2, interval):
    """Get the segment CN values for each targeted region.

    Get overlapping regions of two paired segment/gene sets.

    For genes with 2 or more segments, take the longest segment (or [weighted]
    average).

    Return a pandas.DataFrame with columns:
        chrom, start, end, label, value1, value2
    """
    segments1 = cnvlib.read(cbs1).autosomes()
    segments2 = cnvlib.read(cbs2).autosomes()
    non_overlapping = set(segments1.chromosome).symmetric_difference(
            set(segments2.chromosome))
    if non_overlapping:
        raise ValueError("Mismatched chromosomes: " +
                         ' '.join(sorted(non_overlapping)))
    segments1.sort()
    segments2.sort()

    genes = interval2genes(interval)
    print("#Genes tiled:", len(genes), file=sys.stderr)

    genes["value1"] = [segment_cn(sel)
                       for (_r, sel) in segments1.by_ranges(genes, mode="trim")]
    genes["value2"] = [segment_cn(sel)
                       for (_r, sel) in segments2.by_ranges(genes, mode="trim")]
    genes = genes.data.dropna()
    print("#Genes after dropna:", len(genes), file=sys.stderr)
    return genes


def segment_cn(segset):
    """Get or estimate the segment's copy number.

    If the segment set only contains one row (i.e. segment), just return that
    value. Otherwise, return the average of the segment log2 values.
    """
    if len(segset) == 0:
        return np.nan
    elif len(segset) == 1:
        return segset.log2[0]
    else:
        # Weight by segment sizes
        return np.average(segset.log2, weights=(segset.end - segset.start))
        # Median
        # return segset.log2.median()


# ENH - port to GA/CNA.by_genes, .squash_genes
def interval2genes(interval, min_gene_size=200):
    """Squash intervals into named genes."""
    rarr = RA.read(interval).autosomes()
    rarr = rarr[~rarr.data.name.isin(params.IGNORE_GENE_NAMES + ("Background",))]

    # XXX This will combine non-adjacent genes w/ same name
    # gb = rarr.data.groupby(["chromosome", "name"], as_index=False, sort=False)
    # ag = gb.aggregate({'start': np.min, 'end': np.max})
    # return rarr.as_dataframe(ag)

    curr_name = None
    curr_chrom = None
    curr_start = None
    curr_end = None
    curr_len = 0
    out_rows = []
    for chrom, start, end, name in rarr.coords(also=["name"]):
        if chrom != curr_chrom or name != curr_name:
            if (curr_name is not None
                # Skip CGH probes that are not real targeted genes
                and (curr_len > 1 or curr_end - curr_start >= min_gene_size)
               ):
                out_rows.append((curr_chrom, curr_start, curr_end, curr_name))
            # Reset
            curr_name = name
            curr_chrom = chrom
            curr_start = start
            curr_len = 0
        # Extend
        curr_end = end
        curr_len += 1
    if curr_name is not None and curr_len >= 1:
        out_rows.append((curr_chrom, curr_start, curr_end, curr_name))
    return RA.from_rows(out_rows,
                        columns=["chromosome", "start", "end", "label"])


def main(args):
    """Make and emit the table."""
    table = read_paired_genes(args.asegment, args.bsegment, args.interval)
    table.to_csv(args.output or sys.stdout, index=False)


if __name__ == '__main__':
    import argparse
    AP = argparse.ArgumentParser(description=__doc__)
    AP.add_argument("asegment", help="Segmentation calls")
    AP.add_argument("bsegment", help="Segmentation calls")
    AP.add_argument("-i", "--interval", help="Target intervals list")
    AP.add_argument("-o", "--output", help="Output CSV file name")
    main(AP.parse_args())
