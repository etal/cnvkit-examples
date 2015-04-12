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
    if fname and os.stat(fname).st_size > 0:
        arr = np.loadtxt(fname)
        print("Loaded", fname)
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
    AP.add_argument("contra_pool_tr", help="CONTRA (TR cohort, pooled)")
    AP.add_argument("contra_pool_ex", help="CONTRA (EX cohort, pooled)")
    AP.add_argument("contra_pair_tr", help="CONTRA (TR cohort, paired)")
    AP.add_argument("contra_pair_ex", help="CONTRA (EX cohort, paired)")
    AP.add_argument("copywriter_tr", help="CopywriteR (TR cohort)")
    AP.add_argument("copywriter_ex", help="CopywriteR (EX cohort)")
    AP.add_argument("-o", "--output", help="Output PDF filename.")
    args = AP.parse_args()
    df = pd.concat([
        as_dframe(args.cnvkit_pool_tr, 'CNVkit-pool', 'TR'),
        as_dframe(args.cnvkit_pool_ex, 'CNVkit-pool', 'EX'),
        as_dframe(args.cnvkit_pair_tr, 'CNVkit-pair', 'TR'),
        as_dframe(args.cnvkit_pair_ex, 'CNVkit-pair', 'EX'),
        as_dframe(args.cnvkit_flat_tr, 'CNVkit-flat', 'TR'),
        as_dframe(args.cnvkit_flat_ex, 'CNVkit-flat', 'EX'),
        as_dframe(args.contra_pool_tr, 'CONTRA-pool', 'TR'),
        as_dframe(args.contra_pool_ex, 'CONTRA-pool', 'EX'),
        as_dframe(args.contra_pair_tr, 'CONTRA-pair', 'TR'),
        as_dframe(args.contra_pair_ex, 'CONTRA-pair', 'EX'),
        as_dframe(args.copywriter_tr, 'CopywriteR', 'TR'),
        as_dframe(args.copywriter_ex, 'CopywriteR', 'EX'),
    ])

    sixcolors = (sn.color_palette("Blues")[1:4] +
                 sn.color_palette("Reds")[2:4] +
                 ["gold"])

    # Cohorts split, methods in subplots
    grid = sn.FacetGrid(df, row="Cohort", row_order=['TR', 'EX'],
                        aspect=2.0, size=3.5, ylim=(-.6, .6),
                        margin_titles=True)
    (grid.map(sn.violinplot, "Method", "Difference", inner='quartile',
              palette=sn.color_palette(sixcolors),
              cut=0, gridsize=500, linewidth=1, width=0.95)
     .set_ylabels("Difference from aCGH")
    )

    if args.output:
        plt.savefig(args.output, format='pdf', bbox_inches="tight")
        print("Wrote", args.output, file=sys.stderr)
        plt.close()
    else:
        plt.show()
