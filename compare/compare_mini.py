#!/usr/bin/env python

"""Plot each method's distribution of residuals as a boxplot."""
from __future__ import division, print_function

import os
import sys
import argparse

import numpy as np
import pandas as pd
import seaborn as sn
from matplotlib import pyplot as plt

sn.set_style("darkgrid")
Y_RANGE = 0.35

def get_argparser():
    AP = argparse.ArgumentParser(description=__doc__)
    AP.add_argument("cnvkit_pool_tr", help="CNVkit (TR cohort, pooled)")
    AP.add_argument("cnvkit_pair_tr", help="CNVkit (TR cohort, paired)")
    AP.add_argument("cnvkit_flat_tr", help="CNVkit (TR cohort, flat ref.)")
    AP.add_argument("cnvkit_pool_ex", help="CNVkit (EX cohort, pooled)")
    AP.add_argument("cnvkit_pair_ex", help="CNVkit (EX cohort, paired)")
    AP.add_argument("cnvkit_flat_ex", help="CNVkit (EX cohort, flat ref.)")
    AP.add_argument("cnvkit_pool_cl", help="CNVkit (CL cohort, pooled)")
    AP.add_argument("cnvkit_pair_cl", help="CNVkit (CL cohort, paired)")
    AP.add_argument("cnvkit_flat_cl", help="CNVkit (CL cohort, flat ref.)")
    AP.add_argument("-o", "--output", help="Output PDF filename.")
    return AP


def load_inputs(args):
    return pd.concat([
        as_dframe(args.cnvkit_pool_tr, 'CNVkit\npooled', 'TR'),
        as_dframe(args.cnvkit_pair_tr, 'CNVkit\npaired', 'TR'),
        as_dframe(args.cnvkit_flat_tr, 'CNVkit\nno ref.', 'TR'),
        as_dframe(args.cnvkit_pool_ex, 'CNVkit\npooled', 'EX'),
        as_dframe(args.cnvkit_pair_ex, 'CNVkit\npaired', 'EX'),
        as_dframe(args.cnvkit_flat_ex, 'CNVkit\nno ref.', 'EX'),
        as_dframe(args.cnvkit_pool_cl, 'CNVkit\npooled', 'CL'),
        as_dframe(args.cnvkit_pair_cl, 'CNVkit\npaired', 'CL'),
        as_dframe(args.cnvkit_flat_cl, 'CNVkit\nno ref.', 'CL'),
    ], ignore_index=True)


def as_dframe(fname, method, cohort):
    """Load an aCGH-vs-method file as a DataFrame. Print summary stats."""
    if fname and os.stat(fname).st_size > 1:
        arr = np.loadtxt(fname)
        print("Loaded", fname, file=sys.stderr)
    else:
        # Dummy data
        print("Dummy data for:", cohort, method, file=sys.stderr)
        arr = np.random.triangular(-1.0, 0, 1.0, 300)
    return pd.DataFrame({"Difference": arr, "Method": method, "Cohort": cohort})


def blank_dframe(method, cohort):
    arr = np.zeros(1, dtype=np.float_)
    return pd.DataFrame({"Difference": arr, "Method": method, "Cohort": cohort})


def make_plot(data):
    my_colors = sn.color_palette(sn.color_palette("Blues")[1:4], 3)

    # Cohorts split, methods in subplots
    grid = sn.FacetGrid(data, row="Cohort", row_order=['TR', 'EX', 'CL'],
                        # aspect=2.5, size=3,
                        aspect=1.8, size=3.4,
                        ylim=(-Y_RANGE, Y_RANGE),
                        margin_titles=True)
    return (grid.map(sn.boxplot, "Method", "Difference",
                     sym='', whis=[2.5, 97.5],
                     palette=my_colors, linewidth=1, width=0.7)
            .set_ylabels("Difference from aCGH"))



def label_plot(grid):
    """Calculate summary statistics and label the plot with them.

    Return a dataframe of the summary stats.
    """
    # columns = ["Cohort", "Method", "Mean", "Spread", "Limit", "Low95", "Median", "High95",
    #            "MaxAbs", "N"]
    dframes = []
    for ax, (_idx, data) in zip(grid.axes.flat, grid.facet_data()):
        # ax.grid(b=True, which='major', axis='y', color='w', linewidth=1.5)
        # ax.grid(b=True, which='minor', color='w', linewidth=0.5)
        cohort = data.iat[0, 0]
        for idx, (method, group) in enumerate(data.groupby("Method", sort=False)):
            # Calculate descriptives
            arr = np.asarray(group.Difference)
            low, mid, hi = np.percentile(arr, [2.5, 50.0, 97.5])
            spread = hi - low
            dframes.append(pd.DataFrame.from_items([
                ("Cohort", [cohort]),
                ("Method", method.replace('\n', '_')),
                ("Spread95", spread),
                ("TwoSD", 2 * arr.std()),
                ("Mean", arr.mean()),
                ("Median", mid),
                ("Low95", low),
                ("High95", hi),
                ("MaxAbs", np.absolute(arr).max()),
                ("N", len(arr)),
            ]))
            ax.get_ygridlines()[4].set_linewidth(3.0)
            if low > -Y_RANGE + .02:
                ax.text(idx,
                        low - .01,
                        "%.3g" % spread,
                        fontsize=8,
                        verticalalignment="top",
                        horizontalalignment="center")
            else:
                ax.text(idx,
                        -Y_RANGE + .01,
                        " %.3g" % spread,
                        fontsize=8,
                        horizontalalignment="left")
    return pd.concat(dframes, ignore_index=True)


def test_data():
    return load_inputs(get_argparser().parse_args([
        "tr-cnvkit-pool.diffs.dat", "tr-cnvkit-pair.diffs.dat",
        "tr-cnvkit-flat.diffs.dat",
        "ex-cnvkit-pool.diffs.dat", "ex-cnvkit-pair.diffs.dat",
        "ex-cnvkit-flat.diffs.dat",
        "cl-cnvkit-pool.diffs.dat", "cl-cnvkit-pair.diffs.dat",
        "cl-cnvkit-flat.diffs.dat",
    ]))



if __name__ == '__main__':
    args = get_argparser().parse_args()
    stats = label_plot(make_plot(load_inputs(args)))
    stats.to_csv(sys.stdout, sep='\t', float_format="%.5f", index=False)

    if args.output:
        plt.savefig(args.output, format='pdf', bbox_inches="tight")
        print("Wrote", args.output, file=sys.stderr)
        plt.close()
    else:
        plt.show()
