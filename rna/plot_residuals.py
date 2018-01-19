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


def load_tables(acgh_fname, rna_fnames):
    """Load and crunch the input datasets."""
    acgh_table = pd.read_table(acgh_fname, index_col=0)
    print("Loaded", acgh_table.shape, "genes x aCGH samples")
    rna_tables, rna_labels = zip(*((pd.read_table(fname, index_col=0),
                                    basename(fname))
                                   for fname in rna_fnames))
    print("Loading RNA samples incrementally")
    return extract_residuals(acgh_table, rna_tables, rna_labels)


def basename(path):
    name = os.path.basename(path).split('.', 1)[0]
    if name.startswith("tcga-rna-"):
        name = name[len("tcga-rna-"):]
    if name.endswith("-genes"):
        name = name[:-len("-genes")]
    return name


def extract_residuals(acgh_table, rna_tables, rna_labels):
    """Concatenate paired samples as arrays of per-gene log2 values.

    Subtract aCGH value from RNA value for each gene.

    Emit 2 columns: differences, and table label (from filename)

    Returns a 2-column array/DataFrame: RNA estimate residuals, and labels
    indicated which input table the RNA estimates came from.

    Sample and gene name indices are lost here, and all sample-gene cells with
    any NA values (in either RNA or aCGH table) are dropped.
    """
    dframes = []
    for rna_table, rna_label in zip(rna_tables, rna_labels):
        print("Processing table", rna_label, "with shape", rna_table.shape)
        # Match columns by sample and rows by gene, 1:1
        acgh_aligned, rna_aligned = acgh_table.align(rna_table, join='inner')
        print("Trimmed input tables to shape:", acgh_aligned.shape)
        # Calculate residuals of RNA gene-ratio estimates from aCGH segments
        resids = acgh_aligned.values.ravel() - rna_aligned.values.ravel()
        # Drop cells where either table is NaN/missing data
        resids = resids[~np.isnan(resids)]
        dframe = pd.DataFrame({'residual': resids,
                               'estimator': np.repeat(rna_label, len(resids))})
        dframes.append(dframe)
    result = pd.concat(dframes, axis=0)
    print("Final table shape:", result.shape)
    return result


def plot_residuals(table, output=None):
    """Violin(?) plot of residuals (y) from TCGA aCGH log2 (x)."""
    # Fast
    sn.boxplot(x='estimator', y='residual', data=table,
               #  width=.8,
               # See matplotlib.pyplot.boxplot
               #  sym='',  # Hide "fliers", the outlier points
               showfliers=False,
               )
    # Slow
    # sn.violinplot(x='estimator', y='residual', data=table)
    if output:
        plt.savefig(output, format='png', bbox_inches=0)
        print("Wrote", output, file=sys.stderr)
    else:
        plt.show()


if __name__ == '__main__':
    import argparse
    AP = argparse.ArgumentParser(description=__doc__)
    AP.add_argument('acgh_fname', help="TCGA aCGH-based table")
    AP.add_argument('rna_fnames', nargs='+', help="CNVkit RNA-based tables")
    AP.add_argument('-o', '--output', help="Output filename (PDF)")
    args = AP.parse_args()

    table_resid = load_tables(args.acgh_fname, args.rna_fnames)
    plot_residuals(table_resid, args.output)
