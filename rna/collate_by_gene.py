#!/usr/bin/env python
"""Load .cnr or .cns files into a table of genes vs. sample log2 ratios.

Output that combined table as TSV.
"""
from __future__ import absolute_import, division, print_function

import os
import sys

import pandas as pd

import cnvlib


def basename(path):
    fname = os.path.basename(path)
    return fname.split('.', 1)[0]


def cnx2series(fname, min_weight=0, is_segment=False):
    """Load .cnr or .cns file, extract 'log2' and 'gene' columns.

    With `is_segment`, unpack genes in each segment of the input .cns file.

    Returns: Series of log2 ratios indexed by gene names.
    """
    d = cnvlib.read(fname).autosomes().data
    if min_weight:
        ok_wt = d['weight'] >= min_weight
        d = d[ok_wt]
        print("Dropped", (~ok_wt).sum(), "rows with weight below", min_weight)
    if is_segment:
        # Skip segments that cover no genes
        d = d[d.gene != '-']
        # Create a row from each gene, building series by each segment
        chunks = []
        for gene, log2 in d.loc[:, ('gene', 'log2')].itertuples(index=False):
            genes = gene.split(',')
            chunks.append(pd.Series(log2, index=genes))
        ser = pd.concat(chunks)
    else:
        ser = pd.Series(d.log2.values, index=d.gene.values)
    ser = ser.sort_index().rename(basename(fname))
    # Drop any rows genes with duplicate gene names
    if not ser.index.is_unique:
        dup_idx = ser.index.duplicated(keep=False)
        print("Found", dup_idx.sum(), "duplicated gene names in",
              fname, file=sys.stderr)
        ser = ser[~dup_idx]
    return ser


if __name__ == '__main__':
    import argparse
    AP = argparse.ArgumentParser(description=__doc__)
    AP.add_argument('fnames', nargs='+')
    AP.add_argument('-s', '--segmented', action='store_true',
                    help="Input files are .cns, multiple genes per row.")
    AP.add_argument('-w', '--min-weight', type=int, default=0,
                    help="Minimum segment or bin weight to keep.")
    AP.add_argument('-o', '--output',
                    help="Output filename (*.tsv)")
    args = AP.parse_args()

    print("Expecting", ".cns" if args.segmented else ".cnr", "files",
          file=sys.stderr)
    sample_cols = [cnx2series(fname, args.min_weight, args.segmented)
                   for fname in args.fnames]
    print("Loaded", len(sample_cols), "samples", file=sys.stderr)

    data = pd.concat(sample_cols, axis=1)
    data.to_csv(args.output, sep='\t', index=True)
    print("Wrote", args.output, "with", len(data), "rows")
