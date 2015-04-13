#!/usr/bin/env python

"""Bland & Altman's approach to evaluating alternative assays.

Summarize differences in assay values, and plot vs. the mean of both assays.
"""
from __future__ import division, print_function

import sys
# from os.path import basename

import numpy
import pandas
import seaborn
from matplotlib import pyplot

seaborn.set_style("ticks")


def plot_diffs_vs_means(name, diffs, means):
    grid = seaborn.jointplot(means, diffs, kind='scatter', stat_func=None,
                             space=0, size=6, ylim=(-2.0, 2.0),
                             marginal_kws={'bins': 60,
                                           'hist_kws': {
                                               'histtype': 'stepfilled',
                                               'range': (-2, 2),
                                           }},
                             joint_kws={'alpha': 0.2})
    grid.ax_joint.set_xlabel("Mean estimated log2 ratio")
    grid.ax_joint.set_ylabel("Difference from aCGH")
    grid.ax_joint.axhline(c='k', lw=1, zorder=-10)

    # Limits of reliability at 95%
    low, mid, hi = numpy.percentile(diffs, [2.5, 50.0, 97.5])
    grid.ax_joint.axhline(mid, c='k', ls=':', lw=1, zorder=-11)
    grid.ax_joint.axhline(low, c='k', ls=':', lw=1, zorder=-11)
    grid.ax_joint.axhline(hi, c='k', ls=':', lw=1, zorder=-11)

    fig_fname = name + ".alt.pdf"
    pyplot.savefig(fig_fname, format='pdf', bbox_inches="tight")
    print("Wrote", fig_fname, file=sys.stderr)
    pyplot.close()


if __name__ == '__main__':
    import argparse
    AP = argparse.ArgumentParser(description=__doc__)
    AP.add_argument("fnames", nargs='+', help="*.pair.csv")
    AP.add_argument("-n", "--name", required=True, help="Name")
    args = AP.parse_args()

    all_diffs = []
    all_means = []
    for fname in args.fnames:
        table = pandas.read_csv(fname)
        diffs = table['value2'] - table['value1']
        means = .5*(table['value1'] + table['value2'])
        # plot_diffs_vs_means(basename(fname).split('.', 1)[0], diffs, means)
        all_diffs.append(diffs)
        all_means.append(means)

    plot_diffs_vs_means(args.name,
                        numpy.concatenate(all_diffs),
                        numpy.concatenate(all_means))
    with open(args.name + ".diffs.dat", 'w') as handle:
        handle.writelines("%s\n" % d for d in numpy.concatenate(all_diffs))

