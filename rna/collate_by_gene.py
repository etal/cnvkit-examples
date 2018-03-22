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

    # Drop genes that aren't also listed in the .cnr/.cns file?
    print("Filtering out bad gene names from gene_info")
    if 'probes' in d.columns:
        # It's segments -- multiple genes
        ok_gene_names = set()
        for x in d['gene'].str.split(','):
            ok_gene_names.update(x)
    else:
        ok_gene_names = d['gene']
    mask_to_keep = gene_info['gene'].isin(ok_gene_names)
    print("Keeping", mask_to_keep.sum(), "/", len(mask_to_keep),
          "gene names in gene_info")
    gene_info = gene_info[mask_to_keep]

    chunks = []
    for _chrom, info_rows, cnx_rows in by_shared_chroms(gene_info, d, False):
        info_midpoints = info_rows['midpoint'].values
        info_genes = info_rows['gene'].values
        # Locate which segments/bins each gene midpoint falls within
        # - Compare both start and end to ensure (start <= midpoint < end)
        # - If not, then skip that gene
        cnx_starts = cnx_rows['start'].values
        starts_idx = cnx_starts.searchsorted(info_midpoints, 'right')
        cnx_ends = cnx_rows['end'].values
        ends_idx = cnx_ends.searchsorted(info_midpoints, 'right')
        ok_genes_mask = (starts_idx == ends_idx + 1)
        genes_in_cnx_idx = starts_idx.take(ok_genes_mask.nonzero()[0]) - 1
        gene_log2 = cnx_rows['log2'].values[genes_in_cnx_idx]
        gene_sizes = (cnx_ends - cnx_starts)[genes_in_cnx_idx]
        # Stash 'em, including gene name
        chunk_df = pd.DataFrame({'gene': info_genes[ok_genes_mask],
                                 'log2': gene_log2,
                                 'size': gene_sizes})
        chunks.append(chunk_df)

    df = pd.concat(chunks)
    # Drop any rows genes with duplicate gene names
    if not df['gene'].is_unique:
        dup_idx = df['gene'].duplicated(keep=False)
        print("Found", dup_idx.sum(), "duplicated gene names in",
              fname, file=sys.stderr)
        df = df[~dup_idx]
    df = df.set_index('gene').sort_index()
    return basename(fname), df


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
    AP.add_argument('-w', '--min-weight', type=int, default=0,
                    help="Minimum segment or bin weight to keep.")
    AP.add_argument('-o', '--output',
                    help="Output filename (*.tsv)")
    AP.add_argument('-s', '--sizes',
                    help="Output filename for containin-segment sizes (*.tsv).")
    args = AP.parse_args()

    gene_info = load_gene_midpoints(args.gene_resource)

    bnames, dframes = zip(*[load_cnx(fname, gene_info, args.min_weight)
                            for fname in args.fnames])
    print("Loaded", len(bnames), "samples", file=sys.stderr)

    # Write log2 values and containing-segment sizes to separate files
    all_log2 = pd.concat([df['log2'] for df in dframes], axis=1)
    all_log2.columns = bnames
    all_log2.to_csv(args.output, sep='\t', index=True)
    print("Wrote", args.output, "with", len(all_log2), "rows")

    if args.sizes:
        all_sizes = pd.concat([df['size'] for df in dframes], axis=1)
        all_sizes.columns = bnames
        all_sizes.to_csv(args.sizes, sep='\t', index=True)
        print("Wrote", args.sizes, "with", len(all_sizes), "rows")
