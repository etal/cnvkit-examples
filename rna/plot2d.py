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

seaborn.set(font='Sans', style='whitegrid')


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


def extract_xys(rna_table, acgh_table, size_table):
    """Concatenate paired samples as arrays of per-gene log2 values.

    Returns a 3-column array/DataFrame of all samples' per-gene log2 values,
    concatenating all samples and losing the gene name index. All rows with any
    NA values are dropped.

    Segment sizes are replaced with the labels '<5MB', '5-50MB', '>50MB'.
    """
    tables = rna_table, acgh_table, size_table
    tables = align_indices(tables)
    print("Trimmed all 3 input tables to shape:", rna_table.shape)
    # Collect log2 and segment size values from matched samples
    chunks = []
    for rna_col, acgh_col, size_col in zip(*[t.values.T for t in tables]):
        trips = np.vstack([rna_col, acgh_col, size_col])
        # Drop rows where either value is missing/NaN
        row_nan = np.isnan(trips).any(axis=0)
        trips = trips[:, ~row_nan]
        chunks.append(trips)
    allpairs = np.hstack(chunks).T
    df = pd.DataFrame(allpairs, columns=('RNA', 'aCGH', 'Size'))
    size_labels = np.repeat('5-50MB', len(df))
    size_labels[df['Size'] < 5e6] = '<5MB'
    size_labels[df['Size'] > 5e7] = '>50MB'
    df['Size'] = size_labels
    return df


def align_indices(tables):
    """Reduce all tables to the shared index and column values.

    Return the same tables with intersected, sorted indices.
    """
    common_index = intersect_all([t.index.values for t in tables])
    common_columns = intersect_all([t.columns.values for t in tables])
    out_tables = [table.loc[common_index, common_columns]
                  for table in tables]
    return out_tables


def intersect_all(sers):
    common = sers[0]
    for other in sers[1:]:
        common = common[np.isin(common, other)]
    common.sort()
    return common


def color_alpha(alpha, core_color=(0.298, 0.447, 0.69)):
    return core_color + (alpha,)


def corr_stats(x, y):
    r_coef, _p_val = pearsonr(x, y)
    n = len(x)
    return r_coef, n


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


def plot_paired_genes_facet(table, output=None):
    """Hexbin plot of RNA vs. aCGH log2 values.

    Facet by segment size: <5MB<50MB<
    """
    xymin = min(min(table['RNA']), min(table['aCGH']))
    xymax = max(max(table['RNA']), max(table['aCGH']))
    pad = 0.3
    xy_limits = (xymin - pad, xymax + pad)
    nbins = 50
    size_labels = ['<5MB', '5-50MB', '>50MB']
    grid = seaborn.FacetGrid(table, col='Size',
                             col_order=size_labels,
                             xlim=xy_limits, ylim=xy_limits,
                             # subplot_kws={'aspect': 1},
                             margin_titles=True)

    # This shows blue lines around the hex bins :'(
    #  grid.map(plt.hexbin, 'aCGH', 'RNA',
    #           bins='log', gridsize=nbins, mincnt=1,
    #           edgecolors='red',
    #          )

    # Plot a diagonal line, summary stats, and hex bins
    for ax, size_label in zip(grid.axes.flat, size_labels):
        # Draw the 1:1 diagonal as an overlaid line plot
        # Use the legend to draw the annotation
        # (adapted from seaborn.axisgrid.JointGrid.annotate)
        data_subset = table[table['Size'] == size_label]
        pearson_r, _p = corr_stats(data_subset['aCGH'], data_subset['RNA'])
        N = len(data_subset)
        annotation = "Pearson r = {:.3f}\nN = {}".format(pearson_r, N)
        ax.text(xy_limits[0] + 1, xy_limits[1] - 1, annotation,
                fontsize='x-small', verticalalignment='top')
        ax.plot(xy_limits, xy_limits, color='lightgray',
                linestyle='-', linewidth=1, zorder=-1)
        ax.hexbin(data_subset['aCGH'], data_subset['RNA'],
                  bins='log', gridsize=nbins, mincnt=1,
                 )
        ax.set_xlabel("aCGH")
    grid.axes.flat[0].set_ylabel("RNA")

    if output:
        plt.savefig(output, format='pdf', bbox_inches='tight')
        print("Wrote", output, file=sys.stderr)
    else:
        plt.show()



if __name__ == '__main__':
    import argparse
    AP = argparse.ArgumentParser(description=__doc__)
    AP.add_argument('rna', help="CNVkit RNA-based table")
    AP.add_argument('acgh', help="TCGA aCGH-based table")
    AP.add_argument('-s', '--sizes', help="aCGH segment sizes for faceting")
    AP.add_argument('-o', '--output', help="Output filename (PDF).")
    args = AP.parse_args()

    all_rna = pd.read_table(args.rna, index_col=0)
    print("Loaded", all_rna.shape, "genes x RNA samples")

    all_acgh = pd.read_table(args.acgh, index_col=0)
    print("Loaded", all_acgh.shape, "genes x aCGH samples")

    if args.sizes:
        all_sizes = pd.read_table(args.sizes, index_col=0)
        print("Loaded", all_sizes.shape, "genes x aCGH sample segment sizes")
        table_xys = extract_xys(all_rna, all_acgh, all_sizes)
        plot_paired_genes_facet(table_xys, args.output)

    else:
        table_xy = extract_xy(all_rna, all_acgh)
        plot_paired_genes(table_xy, args.output)
