#!/usr/bin/env python

"""Match up and aggregate gene coverages between 2 sets of samples."""
from __future__ import division, print_function

import sys

import pandas
import numpy as np

import cnvlib
from cnvlib import ngfrills


# --- by targeted interval ---

def is_skipped_chromosome(chrom):
    return (chrom in ('X', 'Y', 'chrX', 'chrY') or
            chrom.startswith(('chrUn_', 'Un_')) or
            chrom.endswith('_random'))

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

    genes = list(interval2genes(interval))
    print("#Genes tiled:", len(genes), file=sys.stderr)

    has_chr = segments1.chromosome[0].startswith('chr')
    for chrom, start, end, name in genes:
    # for chrom, start, end, name in interval2genes(interval):
        if is_skipped_chromosome(chrom):
            continue
        if not has_chr:
            # Remove the 'chr' prefix from target gene chromosome name
            chrom = chrom[3:]
        sel1 = segments1.in_range(chrom, start, end, trim=True)
        sel2 = segments2.in_range(chrom, start, end, trim=True)
        # if (sel1['probes'] < 3).all() or (sel2['probes'] < 3).all():
        if len(sel1) == 0 or len(sel2) == 0:
            print("Skipping", name, "-- not covered by a segment")
            continue
        val1 = segment_cn(sel1)
        val2 = segment_cn(sel2)
        yield (chrom, val1, val2, start, end, name)


def segment_cn(segset):
    """Get or estimate the segment's copy number.

    If the segment set only contains one probe (i.e. segment), just return that
    value. Otherwise, return the weighted mean of the segment CNs.
    """
    if len(segset) == 0:
        raise ValueError("WTF: %s" % segset)
    elif len(segset) == 1:
        return segset.coverage[0]
    else:
        ### Weighted mean
        coverages = segset.coverage
        sizes = segset.end - segset.start
        avg = sum(coverages * sizes) / sum(sizes)
        assert min(coverages) <= avg <= max(coverages)
        return avg
        # OR #
        ### Value of segment with the most probes
        # max_idx = segset['probes'].argmax()
        # if max_idx.shape:
        #     # Tie: average of the 2+ winners
        #     return segset.coverage[max_idx].mean()
        # else:
        #     return segset.coverage[max_idx]


def interval2genes(interval, #skip=('CGH', '-'),
                  ):
    """Squash intervals into named genes."""
    curr_name = None
    curr_chrom = None
    curr_start = None
    curr_end = None
    curr_len = 0
    for chrom, start, end, name in ngfrills.parse_regions(interval):
        if name == '-': #in skip:
            continue
        if chrom != curr_chrom or name != curr_name:
            if curr_name is not None:
                yield (curr_chrom, curr_start, curr_end, curr_name)
                # if curr_len > 1:
                #     # Emit
                #     yield (curr_chrom, curr_start, curr_end, curr_name)
                # else:
                #     print("Single-probe gene is probably CGH:", curr_name,
                #           file=sys.stderr)
            # Reset
            curr_name = name
            curr_chrom = chrom
            curr_start = start
            curr_len = 0
        # Extend
        curr_end = end
        curr_len += 1
    if curr_name is not None:
        if curr_len > 1:
            # Emit
            yield (curr_chrom, curr_start, curr_end, curr_name)
        else:
            print("Single-probe gene is probably CGH:", curr_name,
                    file=sys.stderr)


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


# def residual_stats(table):
#     # mids = .5 * (table['value1'] + table['value2'])
#     # resids = table['value1'] - mids
#     resids = .5*(table['value1'] - table['value2'])
#     mean = resids.mean()
#     two_sd = 2 * resids.std()
#     print("Mean:", mean)
#     print("2*SD:", two_sd)
#     return resids, mean, two_sd


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
