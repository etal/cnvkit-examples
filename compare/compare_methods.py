#!/usr/bin/env python

"""Plot each method's distribution of residuals as a boxplot."""
from __future__ import division, print_function

import os
import sys

import numpy as np
import pandas as pd
import seaborn as sn
from matplotlib import pyplot as plt

sn.set_style("darkgrid")


def as_dframe(fname, method, cohort, ymin=-1.0, ymax=1.0):
    if fname and os.stat(fname).st_size > 1:
        arr = np.loadtxt(fname)
        print("Loaded", fname, file=sys.stderr)
        # Stats
        mean = arr.mean()
        limit = 1.96 * arr.std()
        low, mid, hi = np.percentile(arr, [2.5, 50.0, 97.5])
        absmax = np.absolute(arr).max()
        n = len(arr)
        print(fname, *(["%.5f" % val
                       for val in (mean, limit, low, mid, hi, absmax)]
                      + [n]),
              sep='\t')
    else:
        # Dummy data
        arr = np.random.standard_t(1, (300))/10
        arr = arr[(ymin <= arr) & (arr <= ymax)]
    return pd.DataFrame({"Difference": arr, "Method": method, "Cohort": cohort})


if __name__ == '__main__':
    import argparse
    AP = argparse.ArgumentParser(description=__doc__)
    AP.add_argument("cnvkit_pool_tr", help="CNVkit (TR cohort, pooled)")
    AP.add_argument("cnvkit_pool_ex", help="CNVkit (EX cohort, pooled)")
    AP.add_argument("cnvkit_pair_tr", help="CNVkit (TR cohort, paired)")
    AP.add_argument("cnvkit_pair_ex", help="CNVkit (EX cohort, paired)")
    AP.add_argument("cnvkit_flat_tr", help="CNVkit (TR cohort, flat ref.)")
    AP.add_argument("cnvkit_flat_ex", help="CNVkit (EX cohort, flat ref.)")
    AP.add_argument("copywriter_pair_tr", help="CopywriteR (TR cohort, paired)")
    AP.add_argument("copywriter_pair_ex", help="CopywriteR (EX cohort, paired)")
    AP.add_argument("copywriter_noref_tr", help="CopywriteR (TR cohort, no ref.)")
    AP.add_argument("copywriter_noref_ex", help="CopywriteR (EX cohort, no ref.)")
    AP.add_argument("contra_pool_tr", help="CONTRA (TR cohort, pooled)")
    AP.add_argument("contra_pool_ex", help="CONTRA (EX cohort, pooled)")
    AP.add_argument("contra_pair_tr", help="CONTRA (TR cohort, paired)")
    AP.add_argument("contra_pair_ex", help="CONTRA (EX cohort, paired)")
    AP.add_argument("-o", "--output", help="Output PDF filename.")
    args = AP.parse_args()

    # Stats table on stdout
    print("File", "Mean", "Limit", "Low95", "Median", "High95", "MaxAbs", "N",
          sep='\t')

    my_colors = sn.color_palette(
        sn.color_palette("Blues")[1:4] +
        ["khaki", "gold"] +
        sn.color_palette("Reds")[2:4],
        7)

    df = pd.concat([
        as_dframe(args.cnvkit_pool_tr, 'CNVkit\npooled', 'TR'),
        as_dframe(args.cnvkit_pool_ex, 'CNVkit\npooled', 'EX'),
        as_dframe(args.cnvkit_pair_tr, 'CNVkit\npaired', 'TR'),
        as_dframe(args.cnvkit_pair_ex, 'CNVkit\npaired', 'EX'),
        as_dframe(args.cnvkit_flat_tr, 'CNVkit\nno ref.', 'TR'),
        as_dframe(args.cnvkit_flat_ex, 'CNVkit\nno ref.', 'EX'),
        as_dframe(args.copywriter_pair_tr, 'CopywriteR\npaired', 'TR'),
        as_dframe(args.copywriter_pair_ex, 'CopywriteR\npaired', 'EX'),
        as_dframe(args.copywriter_noref_tr, 'CopywriteR\nno ref.', 'TR'),
        as_dframe(args.copywriter_noref_ex, 'CopywriteR\nno ref.', 'EX'),
        as_dframe(args.contra_pool_tr, 'CONTRA\npooled', 'TR'),
        as_dframe(args.contra_pool_ex, 'CONTRA\npooled', 'EX'),
        as_dframe(args.contra_pair_tr, 'CONTRA\npaired', 'TR'),
        as_dframe(args.contra_pair_ex, 'CONTRA\npaired', 'EX'),
    ])

    # Cohorts split, methods in subplots
    grid = sn.FacetGrid(df, row="Cohort", row_order=['TR', 'EX'],
                        aspect=2.0, size=3.5, ylim=(-.6, .6),
                        margin_titles=True)
    (grid.map(sn.violinplot, "Method", "Difference", inner='quartile',
              palette=my_colors,
              cut=0, gridsize=500, linewidth=1, width=0.95)
     .set_ylabels("Difference from aCGH")
    )

    if args.output:
        plt.savefig(args.output, format='pdf', bbox_inches="tight")
        print("Wrote", args.output, file=sys.stderr)
        plt.close()
    else:
        plt.show()
