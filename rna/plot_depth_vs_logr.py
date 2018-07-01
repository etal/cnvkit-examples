#!/usr/bin/env python
"""
"""
from __future__ import division, print_function
import sys
import cnvlib
from matplotlib import pyplot as plt
import numpy as np

fnames = sys.argv[1:]
acc = []
logs = []
depths = []
for fname in fnames:
    cnr = cnvlib.read(fname)
    logs.append(cnr['log2'])
    # XXX what to show here? absolute expression, but normalized a bit?
    depths.append(cnr['depth'])

logs = np.concatenate(logs)
depths = np.concatenate(depths)
plt.scatter(depths, logs, marker='o', alpha=.2, edgecolors='none', color='k')
max_x = depths.max() + 1
plt.hlines(0, 0, max_x)
plt.xlim(1, max_x)
plt.xscale('log')
plt.xlabel('Depth')
plt.ylabel('Normalized log-ratio')
plt.show()
