#!/usr/bin/env python

"""Plot FISH observed counts as little heatmaps."""
from __future__ import print_function, division
import collections

import numpy
from matplotlib import pyplot, cm

data = {
    'ALK': [8, 8, 8, 8, 8, 8, 8, 8, 9, 3, 9, 8, 3, 3], # est. 8
    'BRAF': [6, 6, 6, 6, 6, 5], # est. 6
    'MET': [6, 6, 4, 6, 6, 5, 6, 5, 6, 6, 7, 6, 6, 6, 7], # est. 6
    'NTRK1': [8, 11, 12, 11, 9, 12], # est. 11-12
    'RET': [3, 3, 2, 6, 4, 2, 2, 4, 4, 4, 5, 2], # est. 4
    'ROS1': [6, 6, 3, 4, 3, 6, 3, 3, 6, 6, 5, 6], # est. 6
    'XXX': [1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6,
            7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 8,
            9, 9, 9, 9, 9, 9, 9, 9, 9,
            10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
            11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 
            12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12],
}


def as_rgba(val):
    """Val is between 0 and 1."""
    return numpy.asarray((250, 10, 10, int(val * 255)), dtype=numpy.uint8)


arr = numpy.zeros((13, len(data), 4), dtype=numpy.uint8)
for i, (gene, obs) in enumerate(sorted(data.items())):
    # print(gene)
    N = len(obs)
    cnt = collections.Counter(obs)
    max_cnt = cnt.most_common()[0][1]
    # most_common = cnt.most_common()
    for value, count in cnt.items():
        # frac = count / len(obs)
        frac = count / max_cnt
        arr[value, i, 0:4] = as_rgba(frac)
    print(gene, list(arr[:,i,3]))


pyplot.imshow(arr, #cmap=cm.gray,
              interpolation='none', origin='lower')

pyplot.show()

