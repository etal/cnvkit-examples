#!/usr/bin/env python

"""Plot each method's distribution of residuals as a boxplot."""
from __future__ import division, print_function

import os
import sys

import numpy as np
import pandas as pd
import seaborn as sn
from matplotlib import pyplot as plt

sn.set(style="darkgrid", palette="Blues")


def as_dframe(fname, method, cohort, ymin=-1.0, ymax=1.0):
    """Load an aCGH-vs-method file as a DataFrame. Print summary stats."""
    if fname and os.stat(fname).st_size > 1:
        arr = np.loadtxt(fname)
        print("Loaded", fname, file=sys.stderr)
        # Stats
        mean = arr.mean()
        limit = 1.96 * arr.std()
        low, mid, hi = np.percentile(arr, [2.5, 50.0, 97.5])
        spread = hi - low
        absmax = np.absolute(arr).max()
        n = len(arr)
        print(fname, *(["%.5g" % val
                       for val in (mean, spread, limit, low, mid, hi, absmax)]
                      + [n]),
              sep='\t')
    else:
        # Dummy data
        arr = np.random.triangular(ymin, 0, ymax, 300)
    return pd.DataFrame({"Difference": arr, "Method": method, "Cohort": cohort})


if __name__ == '__main__':
    import argparse
    AP = argparse.ArgumentParser(description=__doc__)
    AP.add_argument("cnvkit_pool", help="CNVkit (CL cohort, pooled)")
    # AP.add_argument("cnvkit_pair", help="CNVkit (TR cohort, paired)")
    # AP.add_argument("cnvkit_flat", help="CNVkit (CL cohort, flat ref.)")
    # AP.add_argument("copywriter_pair", help="CopywriteR (TR cohort, paired)")
    # AP.add_argument("copywriter_noref", help="CopywriteR (TR cohort, no ref.)")
    # AP.add_argument("contra_pool", help="CONTRA (TR cohort, pooled)")
    # AP.add_argument("contra_pair", help="CONTRA (TR cohort, paired)")
    AP.add_argument("-o", "--output", help="Output PDF filename.")
    args = AP.parse_args()

    # Stats table on stdout
    print("File", "Mean", "Spread", "Limit", "Low95", "Median", "High95",
          "MaxAbs", "N", sep='\t')

    df = pd.concat([
        as_dframe(args.cnvkit_pool, 'CNVkit\npooled', 'CL'),
        # as_dframe(args.cnvkit_pair, 'CNVkit\npaired', 'CL'),
        # as_dframe(args.cnvkit_flat, 'CNVkit\nno ref.', 'CL'),
        # as_dframe(args.copywriter_pair, 'CopywriteR\npaired', 'CL'),
        # as_dframe(args.copywriter_noref, 'CopywriteR\nno ref.', 'CL'),
        # as_dframe(args.contra_pool, 'CONTRA\npooled', 'CL'),
        # as_dframe(args.contra_pair, 'CONTRA\npaired', 'CL'),
    ])

    # Cohorts split, methods in subplots
    _fig, ax = plt.subplots(figsize=(4.5, 5.0),
                            subplot_kw={"ylim": (-.25, .65)})
    sn.violinplot("Method", "Difference", inner="quartile",
                  cut=0, gridsize=1000, linewidth=1, width=0.95,
                  axis=ax, data=df)
    ax.set_ylabel("Difference from aCGH")

    if args.output:
        plt.savefig(args.output, format='pdf', bbox_inches="tight")
        print("Wrote", args.output, file=sys.stderr)
        plt.close()
    else:
        plt.show()
