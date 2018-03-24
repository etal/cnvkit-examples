#!/usr/bin/env python
"""Plot CNVkit-RNA estimates as deviations from TCGA segmented log2 ratios.

Considering:

- CBS segment means, all
- CBS segment means, weight >= 20
- Arm-level segment means
- Smoothed with Kaiser window
- Smoothed with adaptive window

"""
from __future__ import absolute_import, division, print_function
from builtins import zip
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sn

sn.set(font='Sans', style='whitegrid')


def load_tables(acgh_fname, size_fname, rna_fnames):
    """Load and crunch the input datasets."""
    acgh_table = pd.read_table(acgh_fname, index_col=0)
    print("Loaded", acgh_table.shape, "genes x aCGH samples")
    if size_fname is not None:
        size_table = pd.read_table(size_fname, index_col=0)
        print("Loaded", size_table.shape,
              "genes x aCGH samples's segment sizes")
    else:
        size_table = None
    rna_tables, rna_labels = zip(*((pd.read_table(fname, index_col=0),
                                    basename(fname))
                                   for fname in rna_fnames))
    print("Loading RNA samples incrementally")
    return extract_residuals(acgh_table, size_table, rna_tables, rna_labels)


def basename(path):
    name = os.path.basename(path).split('.', 1)[0]
    if name.startswith("tcga-rna-"):
        name = name[len("tcga-rna-"):]
    if name.endswith("-genes"):
        name = name[:-len("-genes")]
    return name


def extract_residuals(acgh_table, size_table, rna_tables, rna_labels):
    """Concatenate paired samples as arrays of per-gene log2 values.

    Subtract aCGH value from RNA value for each gene.

    Emit 2 columns: differences, and table label (from filename)

    Returns a 2-column array/DataFrame: RNA estimate residuals, and labels
    indicated which input table the RNA estimates came from.

    Sample and gene name indices are lost here, and all sample-gene cells with
    any NA values (in either RNA or aCGH table) are dropped.
    """
    # ENH? calculate & show signal-to-noise ratio (SNR)
    dframes = []
    for rna_table, rna_label in zip(rna_tables, rna_labels):
        print("Processing table", rna_label, "with shape", rna_table.shape)
        # Match columns by sample and rows by gene, 1:1
        acgh_aligned, rna_aligned = acgh_table.align(rna_table, join='inner')
        if size_table is not None:
            ra, size_aligned = rna_aligned.align(size_table, join='inner')
            assert (ra.index == rna_aligned.index).all()
            assert (ra.index == acgh_aligned.index).all()
            assert (ra.index == size_aligned.index).all()
        print("Trimmed input tables to shape:", acgh_aligned.shape)
        # Calculate deviation of RNA gene-ratio estimates from aCGH segments
        resids = acgh_aligned.values.ravel() - rna_aligned.values.ravel()
        # Drop cells where either table is NaN/missing data
        nan_mask = ~np.isnan(resids)
        resids = np.abs(resids[nan_mask])
        dframe = pd.DataFrame({'Deviation': resids,
                               'Method': np.repeat(rna_label, len(resids))})
        if size_table is not None:
            # Handle size_aligned
            assert (len(size_aligned.values.ravel()) ==
                    len(acgh_aligned.values.ravel()))
            sizes = size_aligned.values.ravel()[nan_mask]
            size_labels = np.repeat('5-50MB', len(sizes))
            size_labels[sizes < 5e6] = '<5MB'
            size_labels[sizes > 5e7] = '>50MB'
            dframe['Size'] = size_labels
        dframes.append(dframe)
    result = pd.concat(dframes, axis=0)
    print("Final table shape:", result.shape)
    return result


def plot_residuals(table, output=None):
    """Violin(?) plot of residuals (y) from TCGA aCGH log2 (x)."""
    if 'Size' in table.columns:
        sn.factorplot(x='Method', y='Deviation', data=table, kind='box',
                      row='Size', row_order=['>50MB', '5-50MB', '<5MB'],
                      aspect=2.5, showfliers=False)
    else:
        sn.boxplot(x='Method', y='Deviation', data=table, showfliers=False)
    plt.xticks(rotation=30)
    if output:
        plt.savefig(output, format='pdf', bbox_inches='tight')
        print("Wrote", output, file=sys.stderr)
    else:
        plt.show()


if __name__ == '__main__':
    import argparse
    AP = argparse.ArgumentParser(description=__doc__)
    AP.add_argument('acgh_fname', help="TCGA aCGH-based table")
    AP.add_argument('rna_fnames', nargs='+', help="CNVkit RNA-based tables")
    AP.add_argument('-s', '--sizes', help="aCGH segment sizes for faceting")
    AP.add_argument('-o', '--output', help="Output filename (PDF)")
    args = AP.parse_args()

    table = load_tables(args.acgh_fname, args.sizes, args.rna_fnames)
    plot_residuals(table, args.output)
