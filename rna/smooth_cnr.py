#!/usr/bin/env python
"""Write a copy of the .cnr with smoothed log2 ratios."""
from __future__ import absolute_import, division, print_function
import argparse
import os
import sys

import numpy as np
from numpy.lib.stride_tricks import as_strided

import cnvlib
from cnvlib import smoothing
from skgenome import tabio


def get_sliding_window(a, width):
    """Sliding window over a 2D array.

    Source: https://stackoverflow.com/questions/37447347/dataframe-representation-of-a-rolling-window/41406783#41406783
    """
    # NB: a = df.values or np.vstack([x, y]).T
    assert width >= 0, width
    s0, s1 = a.strides
    assert s0 >= 0, s0
    assert s1 >= 0, s1
    m, n = a.shape
    assert n >= 0, n
    assert m-width+1 >= 0, m-width+1
    return as_strided(a,
                      shape=(m-width+1, width, n),
                      strides=(s0, s0, s1))


def rolling_weighted_average(x, w, window_size):
    """Rolling weighted average with a uniform 'boxcar' window."""
    wing = window_size // 2
    window_size = 2 * wing + 1
    xp = np.r_[x[wing-1::-1], x, x[:-wing-1:-1]]
    wp = np.r_[w[wing-1::-1], w, w[:-wing-1:-1]]
    x_w = np.vstack([xp, wp]).T
    wins = get_sliding_window(x_w, window_size)
    result = np.average(wins[:,:,0], axis=1, weights=wins[:,:,1])
    return result


def clipped_smooth(cnarr, window_size):
    """Clip signal at +/- 3, then apply polynomial smoothing."""
    x = cnarr['log2'].clip(-3, 3).values
    w = cnarr['weight'].values
    wing = smoothing._width2wing(window_size, x)
    assert len(x) == len(w)
    #  signal = smoothing._pad_array(x, wing)
    #  weights = smoothing._pad_array(w, wing)
    return rolling_weighted_average(x, w, 2 * wing + 1)


def smooth_by_arm(cnarr, window):
    logr_chunks = [clipped_smooth(cnarm, window)
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
    outfname = "{}/{}.wsmooth{}.{}".format(args.output_dir, base,
                                           args.window, ext)
    tabio.write(cnr, outfname)
    print("Wrote", outfname, file=sys.stderr)
