#!/usr/bin/env python

"""Convert a CONTRA CBS-segmented table to CNVkit .cns format.

Column names::

    chromosome    = Chr                   chr1
                    Target.Start          1
                    Target.End            131
    probes        = NumberOfTargets       131
    start         = OriStCoordinate       1509023
    end           = OriEndCoordinate      39038936
                    CBS.Mean              -0.0291
    log2          = LogRatios             -0.0290663133430224
                    Above.PValues.Cutoff  0
    gene          = Calls                 No
"""

import argparse
import pandas as pd
from cnvlib.cnary import CopyNumArray as CNA


AP = argparse.ArgumentParser(description=__doc__)
AP.add_argument('contra_cbs1')
AP.add_argument('-o', '--output')
args = AP.parse_args()

d = pd.read_table(args.contra_cbs1)
cnarr = CNA.from_columns(dict(
    chromosome=d['Chr'],
    start=d['OriStCoordinate'],
    end=d['OriEndCoordinate'],
    gene=d['Calls'],
    log2=d['LogRatios'],
    probes=d['NumberOfTargets']))
cnarr.sort()
cnarr.write(args.output)
