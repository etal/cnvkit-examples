#!/usr/bin/env python
"""Plot gene read depths vs. normalized log-ratios.

Show that the normalized log2 ratios are more variable at low transcript read
depths (i.e. low gene expression) than when read depths are deeper (higher
expression).

Columns used:

- depth -- The observed average read depth across the transcribed portion of the
  gene: 'expected_count' from the original RSEM or gene-level read counts (TCGA)
  input file, times read length (100bp), divided by transcript length.

- log2 -- The final bin/gene-level log2 ratio estimates (before smoothing or
  segmentation): depths normalized by sample and gene, log2-scaled; missing/zero
  values replaced by within-gene minima; with GC and transcript-length bias
  correction.

"""
from __future__ import division, print_function
import argparse
import os
import sys

import cnvlib
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn

AP = argparse.ArgumentParser(description=__doc__)
AP.add_argument('cnr_files', nargs='+', help="All sample .cnr files.")
AP.add_argument('-o', '--output', help="Output filename.")
args = AP.parse_args()

seaborn.plotting_context("poster")
seaborn.set(font="Sans", style="darkgrid")



def load_depths_logs(cnr_fnames):
    logs = []
    depths = []
    for fname in cnr_fnames:
        cnr = cnvlib.read(fname)
        logs.append(cnr['log2'])
        depths.append(cnr['depth'])
        # Ninja move
        if len(cnr_fnames) == 1:
            plt.title(cnr.sample_id)

    depths = np.concatenate(depths)
    keep_idx = (depths != 0)
    depths = depths[keep_idx]
    logs = np.concatenate(logs)[keep_idx]
    return pd.DataFrame({'depth': depths, 'log2': logs})



d = load_depths_logs(args.cnr_files)

max_x = d['depth'].max() + 1
nbins = 90
grid = seaborn.jointplot('depth', 'log2', data=d,
                         kind='hex',
                         space=0.03,
                         xlim=(1, max_x),
                         ylim=(-6, 6),
                         stat_func=None,
                         xscale='log',
                         joint_kws=dict(
                             bins='log',
                             #edgecolors='none',
                             gridsize=nbins,
                             mincnt=2,
                             #xscale='log',
                         ),
                         marginal_kws=dict(
                             bins=None,
                             hist_kws=dict(
                                 histtype='stepfilled',
                                 #alpha=None,
                                 #color=(0.298, 0.447, 0.69, 0.4),
                                 #edgecolor=(0.298, 0.447, 0.69, 1.0),
                                 linewidth=1,
                             ),
                         ),
                        )

grid.ax_joint.set_xlabel('Observed read depth')
grid.ax_joint.set_ylabel('Normalized log2 ratio')
grid.ax_joint.hlines(0, 0, max_x)


if args.output:
    oformat = os.path.splitext(args.output)[-1].replace(".", "")
    plt.savefig(args.output, format=oformat, bbox_inches="tight")
    print("Wrote", args.output, file=sys.stderr)
else:
    plt.show()
