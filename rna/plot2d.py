#!/usr/bin/env python
"""Plot CNVkit-RNA estimate vs. TCGA segmented log2 ratios.

Of interest:

- 2D plot with 1:1 line, regression coefficients & Pearson r
- TCGA segment values (x) vs. CNVkit-RNA deviations/residuals (y)
- CNVkit-RNA estimates:
    - CBS segment means, all
    - CBS segment means, weight >= 20
    - Smoothed with Kaiser window
    - Smoothed with adaptive window

Each gene = 1 datapoint --> ~3 million datapoints.

"""
from __future__ import absolute_import, division, print_function
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn

from scipy.stats import pearsonr

seaborn.set(font='Sans', style='darkgrid')


def plot_paired_genes(table, output=None):
    """Hexbin plot of RNA vs. aCGH log2 values."""
    xymin = min(min(table['RNA']), min(table['aCGH']))
    xymax = max(max(table['RNA']), max(table['aCGH']))
    pad = 0.3
    xy_limits = (xymin - pad, xymax + pad)
    nbins = 60

    grid = seaborn.jointplot('aCGH', 'RNA', data=table,
                             kind='hex',
                             space=.03,
                             xlim=xy_limits,
                             ylim=xy_limits,
                             stat_func=corr_stats,
                             annot_kws=dict(
                                 template="Pearson r = {val:.3f}\nN = {p}",
                                 # template="N = {p}",
                                 loc='upper left',
                             ),
                             joint_kws=dict(
                                 bins='log',
                                 gridsize=nbins,
                                 mincnt=1,
                             ),
                             marginal_kws=dict(
                                 bins=nbins,
                                 hist_kws=dict(
                                     # range=(xymin, xymax),
                                     histtype='stepfilled',
                                     alpha=None,
                                     color=color_alpha(0.4),
                                     edgecolor=color_alpha(1.0),
                                     linewidth=1,
                                 ),
                             ),
                            )
    # Plot a diagonal line
    grid.ax_joint.plot(xy_limits, xy_limits, color='white', linestyle='-',
                       linewidth=1, zorder=-1)

    if output:
        plt.savefig(output, format='pdf', bbox_inches=0)
        print("Wrote", output, file=sys.stderr)
    else:
        plt.show()


def color_alpha(alpha, core_color=(0.298, 0.447, 0.69)):
    return core_color + (alpha,)


def corr_stats(x, y):
    r_coef, _p_val = pearsonr(x, y)
    n = len(x)
    return r_coef, n


def extract_xy(rna_table, acgh_table):
    """Concatenate paired samples as arrays of per-gene log2 values.

    Returns a 2-column array/DataFrame of all samples' per-gene log2 values,
    concatenating all samples and losing the gene name index. All rows with any
    NA values are dropped.
    """
    # Match columns by sample and rows by gene, 1:1
    rna_table, acgh_table = rna_table.align(acgh_table, join='inner')
    print("Trimmed input tables to shape:", rna_table.shape)
    # Collect log2 values from paired samples
    chunks = []
    for rna_col, acgh_col in zip(rna_table.values.T, acgh_table.values.T):
        pairs = np.vstack([rna_col, acgh_col])
        # Drop rows where either value is missing/NaN
        row_nan = np.isnan(pairs).any(axis=0)
        pairs = pairs[:, ~row_nan]
        chunks.append(pairs)
    allpairs = np.hstack(chunks).T
    return pd.DataFrame(allpairs, columns=("RNA", "aCGH"))


if __name__ == '__main__':
    import argparse
    AP = argparse.ArgumentParser(description=__doc__)
    AP.add_argument('rna', help="CNVkit RNA-based table")
    AP.add_argument('acgh', help="TCGA aCGH-based table")
    AP.add_argument('-o', '--output',
                    help="Output filename (PDF).")
    args = AP.parse_args()

    all_rna = pd.read_table(args.rna, index_col=0)
    print("Loaded", all_rna.shape, "genes x RNA samples")

    all_acgh = pd.read_table(args.acgh, index_col=0)
    print("Loaded", all_acgh.shape, "genes x aCGH samples")

    table_xy = extract_xy(all_rna, all_acgh)
    plot_paired_genes(table_xy, args.output)
