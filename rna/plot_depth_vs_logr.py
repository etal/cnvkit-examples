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
from matplotlib import pyplot as plt
import numpy as np

AP = argparse.ArgumentParser(description=__doc__)
AP.add_argument('cnr_files', nargs='+', help="All sample .cnr files.")
AP.add_argument('-o', '--output', help="Output filename.")
args = AP.parse_args()

acc = []
logs = []
depths = []
for fname in args.cnr_files:
    cnr = cnvlib.read(fname)
    # XXX current .cnr files in out/ are clipped at +3
    logs.append(cnr['log2'])
    depths.append(cnr['depth'])
    if len(args.cnr_files) == 1:
        plt.title(cnr.sample_id)

logs = np.concatenate(logs)
depths = np.concatenate(depths)
plt.scatter(depths, logs, marker='o', alpha=.2, edgecolors='none', color='k')
max_x = depths.max() + 1
plt.hlines(0, 0, max_x)
plt.xlim(1, max_x)
plt.xscale('log')
plt.xlabel('Observed read depth')
plt.ylabel('Normalized log2 ratio')
if args.output:
    oformat = os.path.splitext(args.output)[-1].replace(".", "")
    plt.savefig(args.output, format=oformat, bbox_inches="tight")
    print("Wrote", args.output, file=sys.stderr)
else:
    plt.show()
