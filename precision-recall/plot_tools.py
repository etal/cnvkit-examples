#!/usr/bin/env python

"""Plot the benchmark data.

Facets:
    Size (Small/Large) -- separate plots? hue?
    Filename (tool x reference) -- X-axis categories
    CNum, IsGain                -- Y-axis categories

    values: y=precision, x=recall
        varying: size (2-3 points, hue)? pct overlap?

"""
from __future__ import division, print_function

import sys
import collections

import matplotlib.pyplot as plt
import pandas as pd
import seaborn


seaborn.set_style("white")
seaborn.set_context("talk")

def make_plot(data):
    filenames = collections.OrderedDict([
        ("CL_seq", "CNVkit (pooled)"),
        ("CL_pair", "CNVkit (paired)"),
        ("CL_flat", "CNVkit (no ref.)"),
        ("CL.cw-pair", "CopywriteR (paired)"),
        ("CL.cw-noref", "CopywriteR (no ref.)"),
        ("CL.contra-pool", "CONTRA (pooled)"),
        ("CL.contra-pair", "CONTRA (paired)"),
    ])
    data.Filename.replace(filenames, inplace=True)

    my_colors = seaborn.color_palette(
        seaborn.color_palette("Blues")[1:4] +
        ["khaki", "gold"] +
        seaborn.color_palette("Reds")[2:4],
        7)

    grid = seaborn.FacetGrid(data,
                             row="Size",
                             row_order=["All", "Large", "Small", "bp"],
                             col="CN",
                             hue="Filename",
                             hue_order=filenames.values(),
                             # ---
                             legend_out=True,
                             despine=True,
                             margin_titles=True,
                             xlim=(0, 1), ylim=(0, 1),
                             # sharex=True, sharey=True,
                             # size=2, aspect=1,
                             palette=my_colors,
                            )
    grid.map(plt.scatter, "Recall", "Precision",
             edgecolor='none', s=50, zorder=3)
    grid.add_legend()
    # grid.set_titles(col_template="{col_name}",
                    #row_template="{row_name}")
                   # )



if __name__ == "__main__":
    import argparse
    AP = argparse.ArgumentParser()
    AP.add_argument("table")
    AP.add_argument("-o", "--output")
    args = AP.parse_args()

    data = pd.read_table(args.table)
    make_plot(data)

    if args.output:
        plt.savefig(args.output, format='pdf', bbox_inches="tight")
        print("Wrote", args.output, file=sys.stderr)
        plt.close()
    else:
        plt.show()
