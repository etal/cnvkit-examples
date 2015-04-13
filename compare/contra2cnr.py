#!/usr/bin/env python

"""Convert a CONTRA bin-level CNA table to CNVkit .cnr format.

Column names::

                    Targeted.Region.ID          1
                    Exon.Number                 0
    gene          = Gene.Sym                    "SSU72"
    chromosome    = Chr                         "chr1"
    start         = OriStCoordinate             1509023
    end           = OriEndCoordinate            1509125
                    Mean.of.LogRatio            0.153
    log2          = Adjusted.Mean.of.LogRatio   0.160315713205354
                    SD.of.LogRatio              0.083
                    Median.of.LogRatio          0.174
                    number.bases                102
                    P.Value                     0.726069918368404
                    Adjusted.P.Value            0.999933877657259
                    gain.loss                   "gain"
                    tumour.rd                   869.101
                    normal.rd                   779.272
                    tumour.rd.ori               892.716
                    normal.rd.ori               757.735
                    MinLogRatio                 -0.07
                    MaxLogRatio                 0.289
                    BinNumber                   2
"""

import argparse
import pandas as pd
from cnvlib.cnarray import CopyNumArray as CNA
from cnvlib.core import fbase


AP = argparse.ArgumentParser(description=__doc__)
AP.add_argument('contra_table')
AP.add_argument('-o', '--output')
args = AP.parse_args()


d = pd.read_table(args.contra_table)
cnarr = CNA(fbase(args.output),
            d['Chr'], d['OriStCoordinate'], d['OriEndCoordinate'],
            d['Gene.Sym'], d['Adjusted.Mean.of.LogRatio'])
cnarr.sort()
cnarr.write(args.output)
