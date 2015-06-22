#!/usr/bin/env python

"""Scatter plot of gene coverages in 2 sample sets."""
from __future__ import division, print_function

import math
import sys

# import numpy as np
import pandas
import seaborn
from matplotlib import pyplot
from scipy import stats

seaborn.set(font='Sans', style='darkgrid')
# -- segment-wise approach --

# take two .cbs files
#   for each chromosome,
#       calculate overlapping regions of the segments
#       enh: drop overly short segments (do a hist and decide)
#   bubble plot: size ~ length of segment


def my_stats(x, y):
    val, p = stats.pearsonr(x, y)
    n = len(x)
    return val, n


def color_alpha(alpha, core_color=[0.298, 0.447, 0.69]):
    """Add an alpha channel (RBGA) to an RGB color triplet."""
    return tuple(core_color + [alpha])


def main(args):
    """."""
    table = pandas.read_csv(args.table)
    xvals = table['value1']
    yvals = table['value2']
    print(table[(xvals < -1) | (yvals < -1) | ((xvals - yvals).abs() > .5)],
          file=sys.stderr)
#     chromosome     start       end   label  value1    value2
# 4         chr1  44401651  44402437    ARTN -0.1335 -1.919546
# 102       chr6  32163209  32191718  NOTCH4 -0.0026 -2.128800
# 153       chr9  21968181  21994454  CDKN2A -3.8740 -4.511400
# 158       chr9  71043927  71142631    PGM5 -1.3695 -1.618300
# 329      chr19  33791879  33793321   CEBPA -1.5691 -0.133600

    # Set axes to be square (same max/min)
    xymin = min(min(xvals), min(yvals)) - 0.2
    xymax = max(max(xvals), max(yvals)) + 0.2
    xy_limits = (xymin, xymax)

    grid = seaborn.jointplot('value1', 'value2', data=table,
                             kind='scatter',
                             space=.03,
                             xlim=xy_limits,
                             ylim=xy_limits,
                             stat_func=my_stats,
                             annot_kws=dict(
                                 template="Pearson r = {val:.3f}\nN = {p}",
                                 loc='upper left',
                             ),
                             joint_kws=dict(
                                     s=50,
                                     color="#4C72B0",
                                     alpha=.2,
                             ),
                             marginal_kws=dict(
                                 bins=int(math.sqrt(len(xvals))),
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

    # --- muy viejo ---
    # fig = pyplot.figure()
    # axes = fig.add_subplot(1, 1, 1)
    # legend_plots = []
    # legend_labels = []
    # for chrom, coords in sorted(chrom_coords.iteritems()):
    #     print("Plotting", chrom, "with", len(coords), "points")
    #     x, y, starts, ends, names = zip(*coords)
    #     p = axes.scatter(x, y, 50, color='#009699', alpha=0.2, marker='o')
    #     # Label points that show copy number change
    #     for xval, yval, start, end, name in zip(x, y, starts, ends, names):
    #         if (abs(xval) > .8 or abs(yval) > .8 or
    #             (abs(xval) > .6 and abs(yval) > .6)):
    #             valign = "bottom" if yval > 0 else "top"
    #             halign = "left" if xval > 0 else "right"
    #             axes.text(xval, yval, name,
    #                       horizontalalignment=halign,
    #                       verticalalignment=valign)
    #     legend_plots.append(p)
    #     legend_labels.append(chrom)
    # axes.axhline(color='k', linestyle='-')
    # axes.axvline(color='k', linestyle='-')
    # Set axes to be square

    # Plot a diagonal line
    grid.ax_joint.plot(xy_limits, xy_limits, color='white', linestyle='-',
                       linewidth=1, zorder=-1)

    # axes.set_title("Copy number changes after TKI treatment")
    # pyplot.title("Copy ratios assessed by array CGH and CNVkit", zorder=10)
    # axes.legend(legend_plots, legend_labels)
    # TODO - fix legend aesthetics -- one circle each, reasonable size
    # TODO - axis labels as arguments
    grid.ax_joint.set_xlabel(args.x_label + " copy ratio (log2)")
    grid.ax_joint.set_ylabel(args.y_label + " copy ratio (log2)")
    # grid.ax_joint.set_xlabel("Lane 1 copy ratio (log2)")
    # grid.ax_joint.set_ylabel("Lane 2 copy ratio (log2)")

    if args.output:
        pyplot.savefig(args.output, format='pdf', bbox_inches=0)
        print("Wrote", args.output, file=sys.stderr)
    else:
        pyplot.show()


if __name__ == '__main__':
    import argparse
    AP = argparse.ArgumentParser(description=__doc__)
    AP.add_argument("table", help="Data table (.csv)")
    AP.add_argument("-o", "--output", help="Output PDF file name")
    AP.add_argument("-x", "--x-label", default="aCGH", help="x-axis label")
    AP.add_argument("-y", "--y-label", default="CNVkit", help="y-axis label")
    main(AP.parse_args())
