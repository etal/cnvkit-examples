#!/usr/bin/env python
"""Load .cnr or .cns files into a table of genes vs. sample log2 ratios.

Output that combined table as TSV.
"""
from __future__ import absolute_import, division, print_function

import os
import sys

import pandas as pd

import cnvlib
from cnvlib.rna import load_gene_info
from skgenome.intersect import by_shared_chroms


def basename(path):
    fname = os.path.basename(path)
    return fname.split('.', 1)[0]


def load_cnx(fname, gene_info, min_weight=0, is_segment=False):
    """Load .cnr or .cns file, extract 'log2' and 'gene' columns.

    With `is_segment`, unpack genes in each segment of the input .cns file.

    Returns: Series of log2 ratios indexed by gene names.

    Example
    -------

    ::

    idx     Segments:     | Midpoints:
    0       0       90      0, 50
    X                       90, 99
    1       100     200     100, 150, 199
    2       200     2000    200, 201, 1000
    X                       3000
    3       5000    5050    5000, 5020
    4       5050    6000    5050, 5500
    X                       6600

    >>> starts.searchsorted(gene_mids, 'right')
    array([1, 1, 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 5, 5, 5])
    >>> ends.searchsorted(gene_mids, 'right')
    array([0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 5])
                 X  X                    X              X   <-gaps

    """
    d = cnvlib.read(fname).autosomes().data
    if min_weight:
        ok_wt = d['weight'] >= min_weight
        d = d[ok_wt]
        print("Dropped", (~ok_wt).sum(), "rows with weight below", min_weight)

    chunks = []
    for _chrom, info_rows, cnx_rows in by_shared_chroms(gene_info, d, False):
        info_midpoints = info_rows['midpoint']
        # Locate which segments/bins each gene midpoint falls within
        # - Compare both start and end to ensure (start <= midpoint < end)
        # - If not, then skip that gene
        starts_idx = cnx_rows['start'].searchsorted(info_midpoints, 'right')
        ends_idx = cnx_rows['end'].searchsorted(info_midpoints, 'right')
        ok_genes_idx = (starts_idx == ends_idx + 1)
        genes_in_cnx_idx = starts_idx.take(ok_genes_idx.nonzero()[0])
        gene_log2 = cnx_rows['log2'][genes_in_cnx_idx]
        gene_sizes = (cnx_rows['end'] - cnx_rows['start'])[genes_in_cnx_idx]
        # Stash 'em, including gene name
        chunk_df = pd.DataFrame({'gene': info_rows.loc[ok_genes_idx, 'gene'].values,
                                 'log2': gene_log2.values,
                                 'size': gene_sizes.values})
        chunks.append(chunk_df)

    df = pd.concat(chunks)
    # Drop any rows genes with duplicate gene names
    if not df['gene'].is_unique:
        dup_idx = df['gene'].duplicated(keep=False)
        print("Found", dup_idx.sum(), "duplicated gene names in",
              fname, file=sys.stderr)
        df = df[~dup_idx]
    # TODO choose log2 vs. weight column here?
    df = df.set_index('gene')
    ser = df['log2'].sort_index().rename(basename(fname))
    return ser


def load_gene_midpoints(gene_resource):
    """Load all genes' midpoint coordinates.

    Return a DataFrame including the columns 'start', 'end', and 'midpoint',
    where unique gene names are the index.
    """
    gene_info = load_gene_info(gene_resource, None, None)
    # For each gene, get the midpoint of start and end
    midpoints = 0.5 * (gene_info['start'] + gene_info['end'])
    return gene_info.assign(midpoint=midpoints)


if __name__ == '__main__':
    import argparse
    AP = argparse.ArgumentParser(description=__doc__)
    AP.add_argument('fnames', nargs='+')
    AP.add_argument('-g', '--gene-resource', metavar="FILE", required=True,
                    # default="data/ensembl-gene-info.hg38.tsv",
                    help="Ensembl BioMart-derived gene info table.")
    AP.add_argument('-s', '--segmented', action='store_true',
                    help="Input files are .cns, multiple genes per row.")
    AP.add_argument('-w', '--min-weight', type=int, default=0,
                    help="Minimum segment or bin weight to keep.")
    AP.add_argument('-o', '--output',
                    help="Output filename (*.tsv)")
    args = AP.parse_args()

    gene_info = load_gene_midpoints(args.gene_resource)

    print("Expecting", ".cns" if args.segmented else ".cnr", "files",
          file=sys.stderr)
    sample_cols = [load_cnx(fname, gene_info, args.min_weight, args.segmented)
                   for fname in args.fnames]
    print("Loaded", len(sample_cols), "samples", file=sys.stderr)

    # TODO write either log2 values or containing-segment sizes
    data = pd.concat(sample_cols, axis=1)
    data.to_csv(args.output, sep='\t', index=True)
    print("Wrote", args.output, "with", len(data), "rows")
