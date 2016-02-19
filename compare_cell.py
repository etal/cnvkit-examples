#!/usr/bin/env python

"""Plot each method's distribution of residuals as a boxplot."""
from __future__ import division, print_function

import os
import sys
import itertools

import numpy as np
import pandas as pd
import seaborn as sn
from matplotlib import pyplot as plt

# from compare_methods import as_dframe


sn.set(style="darkgrid", palette="Blues")
Y_RANGE = 0.5


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


def make_plot(data):
    """Lay out the plot."""
    _fig, ax = plt.subplots(figsize=(4.5, 5.0),
                            subplot_kw={"ylim": (-Y_RANGE, Y_RANGE)})
    sn.boxplot("Method", "Difference", data=data,
               sym='', whis=[2.5, 97.5],
               linewidth=1, width=0.7, ax=ax)
    ax.set_ylabel("Difference from aCGH")
    gridlines = ax.get_ygridlines()
    gridlines[len(gridlines) // 2].set_linewidth(3.0)
    return ax


def label_plot(ax, data):
    """Label the plot with interval sizes."""
    stats = []
    for idx, (method, group) in enumerate(data.groupby("Method", sort=False)):
        # Calculate descriptives
        arr = np.asarray(group.Difference)
        low, mid, hi = np.percentile(arr, [2.5, 50.0, 97.5])
        spread = hi - low
        ax.text(idx,
                low - .01,
                "%.3g" % spread,
                fontsize=9,
                verticalalignment="top",
                horizontalalignment="center")
        stats.append(pd.DataFrame.from_items([
            ("Cohort", ["CL"]),
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
    return pd.concat(stats, ignore_index=True)



if __name__ == '__main__':
    import argparse
    AP = argparse.ArgumentParser(description=__doc__)
    AP.add_argument("tables", nargs='+',
                    help="CL tables (CNVkit pooled, paired, flat)")
    AP.add_argument("-o", "--output", help="Output PDF filename.")
    args = AP.parse_args()

    # Load data
    data = pd.concat([as_dframe(tbl, title, "CL")
                      for tbl, title in itertools.izip(
                          args.tables,
                          ("CNVkit\npooled",
                           'CNVkit\npaired',
                           'CNVkit\nno ref.'))],
                     ignore_index=True)

    ax = make_plot(data)
    stats = label_plot(ax, data)
    # Stats table on stdout
    stats.to_csv(sys.stdout, sep='\t', float_format="%.5g", index=False)

    # Display the plot
    if args.output:
        plt.savefig(args.output, format='pdf', bbox_inches="tight")
        print("Wrote", args.output, file=sys.stderr)
        plt.close()
    else:
        plt.show()
