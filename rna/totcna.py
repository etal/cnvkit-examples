#!/usr/bin/env python
"""Sum of arm-level absolute log2 values.

Input: *.cns, with arm-level segmentation (segment -m none)
Output: table
"""
from __future__ import division, print_function
import sys
import cnvlib

for fname in sys.argv[1:]:
    cna = cnvlib.read(fname)
    stat = (cna.autosomes()['log2'].abs() ** 2).sum()
    print("%.2f" % stat, cna.sample_id, sep='\t')
