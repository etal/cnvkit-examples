#!/usr/bin/env python

"""Consolidate CONTRA bin-level CNA tables into a single table.


Chr	Target.Start	Target.End	NumberOfTargets	OriStCoordinate	OriEndCoordinate	CBS.Mean	LogRatios	Above.PValues.Cutoff	Calls
chr1	1	131	131	1509023	39038936	-0.0291	-0.0290663133430224	0	No
chr1	132	170	39	40363043	45797982	-0.2451	-0.245057178206862	0	No
chr1	171	220	50	45798062	65349158	-0.0054	-0.00539518828914282	0	No
chr1	222	256	35	66049063	111034945	0.2002	0.200192308642405	0	No
chr1	257	448	192	112523599	206648338	-3e-04	-0.000349499565730335	0	No
chr1	449	462	14	206649520	206666453	-0.4329	-0.432890687111112	0	No
chr1	463	507	45	206666595	249066375	-0.0012	-0.00117366914479222	0	No
chr10	508	539	32	1500986	43600642	0.0671	0.0670689112093665	0	No
chr10	540	546	7	43601820	43613929	-0.4362	-0.436171857930254	0	No
"""

# Column names:
# chromosome    * Chr                   chr1
#                 Target.Start          1
#                 Target.End            131
# probes        * NumberOfTargets       131
# start         * OriStCoordinate       1509023
# end           * OriEndCoordinate      39038936
# log2?         * CBS.Mean              -0.0291
# log2?         * LogRatios             -0.0290663133430224
#                 Above.PValues.Cutoff  0
#                 Calls                 No

import sys
import pandas as pd
from cnvlib.cnarray import CopyNumArray as CNA
from cnvlib.core import fbase

for fname in sys.argv[1:]:
    sample_id = fbase(fname)
    d = pd.read_table(fname)
    cnarr = CNA(sample_id, d['Chr'], d['OriStCoordinate'],
                d['OriEndCoordinate'], d['Calls'], d['LogRatios'],
                probes=d['NumberOfTargets'])
    cnarr.sort()
    cnarr.write(sample_id + ".contra.cns")

