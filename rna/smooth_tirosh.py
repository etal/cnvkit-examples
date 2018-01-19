#!/usr/bin/env python
"""Write a copy of the .cnr with smoothed log2 ratios.

Smoothing method is unweighted boxcar, with input signal windsorized at +/- 3.
"""
from __future__ import absolute_import, division, print_function
import argparse
import os
import sys

import numpy as np

import cnvlib
from skgenome import tabio


def clipped_rolling_mean(values, window):
    clipped = values.clip(-3, 3)
    smoothed = clipped.rolling(window, min_periods=1, center=True).mean()
    return smoothed.values


def smooth_by_arm(cnarr, window):
    logr_chunks = [clipped_rolling_mean(cnarm['log2'], window)
                   for _chrom, cnarm in cnarr.by_arm()]
    d = cnarr.data.assign(log2=np.concatenate(logr_chunks))
    return cnarr.as_dataframe(d)


AP = argparse.ArgumentParser(description=__doc__)
AP.add_argument('cnr_fnames', nargs='+')
AP.add_argument('-w', '--window', type=int, default=100,
                help="Window size for smoothing.")
AP.add_argument('-d', '--output-dir', default='.')
args = AP.parse_args()

for fname in args.cnr_fnames:
    cnr = cnvlib.read(fname)
    cnr = smooth_by_arm(cnr, args.window)
    base, ext = os.path.basename(fname).rsplit(".", 1)
    outfname = "{}/{}.tsmooth{}.{}".format(args.output_dir, base,
                                           args.window, ext)
    tabio.write(cnr, outfname)
    print("Wrote", outfname, file=sys.stderr)
