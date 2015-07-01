#!/usr/bin/env python

"""Match up and aggregate gene coverages between 2 sets of samples."""
from __future__ import division, print_function

import sys

import pandas
import numpy as np

import cnvlib
from cnvlib import ngfrills


# --- by aCGH segment ---
MIN_ACGH_PROBES = 10

def read_paired_genes(cbs1, cbs2, interval):
    """Get the segment CN values for each targeted region.

    For genes with 2 or more segments, take the longest segment (or [weighted]
    average).
    """
    segments1 = cnvlib.read(cbs1)
    segments2 = cnvlib.read(cbs2)
    non_overlapping = set(segments1.chromosome).symmetric_difference(
            set(segments2.chromosome))
    non_overlapping = [chrom for chrom in non_overlapping
                       if not is_skipped_chromosome(chrom)]
    if non_overlapping:
        raise ValueError("Mismatched chromosomes: " +
                         ' '.join(sorted(non_overlapping)))
    segments1.sort()
    segments2.sort()

    for s1_chrom, s1_start, s1_end, s1_name, s1_value, s1_probes in segments1:
        if s1_probes < MIN_ACGH_PROBES or is_skipped_chromosome(s1_chrom):
            continue
        s1_name = "{}:{}-{}".format(s1_chrom, s1_start, s1_end)
        seglike2 = segments2.in_range(s1_chrom, s1_start, s1_end, trim=True)
        if len(seglike2) == 0:
            print("Skipping", s1_name, "-- covers no CNVkit segments")
            continue
        s2_value = segment_cn(seglike2)
        yield (s1_chrom, s1_value, s2_value, s1_start, s1_end, s1_name)


def is_skipped_chromosome(chrom):
    return (chrom in ('X', 'Y', 'chrX', 'chrY') or
            chrom.startswith(('chrUn_', 'Un_')) or
            chrom.endswith('_random'))


def segment_cn(segset):
    """Get or estimate the segment's copy number.

    If the segment set only contains one row (i.e. segment), just return that
    value. Otherwise, return the weighted mean of the segment log2 values.
    """
    if len(segset) == 0:
        raise ValueError("WTF: %s" % segset)
    elif len(segset) == 1:
        return segset.coverage[0]
    else:
        return np.average(segset.coverage, weights=(segset.end - segset.start))


def pairs_as_dframe(pairs):
    """Get overlapping regions of two paired segment/gene sets.

    Returns a pandas.DataFrame with columns: chrom, start, end, cns1, cns2
    """
    col_chrom, col_val1, col_val2, col_start, col_end, col_label = zip(*pairs)
    return pandas.DataFrame.from_items([
        ("chromosome", np.asarray(col_chrom)),
        ("start", np.asarray(col_start)),
        ("end", np.asarray(col_end)),
        ("label", np.asarray(col_label)),
        ("value1", np.asarray(col_val1)),
        ("value2", np.asarray(col_val2))])


def main(args):
    """Make and emit the table."""
    chrom_coords = read_paired_genes(args.asegment, args.bsegment,
                                     args.interval)
    table = pairs_as_dframe(list(chrom_coords))
    table.to_csv(args.output or sys.stdout, index=False)


if __name__ == '__main__':
    import argparse
    AP = argparse.ArgumentParser(description=__doc__)
    AP.add_argument("asegment", help="Segmentation calls")
    AP.add_argument("bsegment", help="Segmentation calls")
    AP.add_argument("-i", "--interval", help="Target intervals list")
    AP.add_argument("-o", "--output", help="Output PDF file name")
    main(AP.parse_args())
