#!/usr/bin/env python

"""Consolidate CONTRA bin-level CNA tables into a single table.

Approaches:
    - convert to .cnr; replace missing rows with 0.0 or NA
        - be sure to drop missing vals when benchmarking
    - calculate per-gene averages directly
        - via mean/adjmean/median here

Checks to do first:
    - match up targeted exon index to .cnr/reference.cnn targets
        (which BED file was used as input? nimv2 targets?)
    - check gene names between targets, .cnr, CONTRA
        - did the gene name for MDM2 get magically dropped?

"""

# colnames = [
#     "Targeted.Region.ID", "Exon.Number", "Gene.Sym",
#     "Chr", "OriStCoordinate", "OriEndCoordinate",
#     "Mean.of.LogRatio", "Adjusted.Mean.of.LogRatio", "SD.of.LogRatio",
#     "Median.of.LogRatio", "number.bases", "P.Value", "Adjusted.P.Value",
#     "gain.loss", "tumour.rd", "normal.rd", "tumour.rd.ori", "normal.rd.ori",
#     "MinLogRatio", "MaxLogRatio", "BinNumber"]

# firstrow = [
#     1, 0, "SSU72",
#     "chr1", 1509023, 1509125,
#     0.153, 0.160315713205354, 0.083,
#     0.174, 102, 0.726069918368404, 0.999933877657259,
#     "gain", 869.101, 779.272, 892.716, 757.735,
#     -0.07, 0.289, 2]


# Keep: Chr, OriEndCoordinate, OriEndCoordinate, Gene.Sym, Adjusted.Mean.of.LogRatio
# -> chromosome, start, end, gene, log2, weight=1.0
# -> .cnr

import sys
import pandas as pd
from cnvlib.cnarray import CopyNumArray as CNA
from cnvlib.core import fbase

for fname in sys.argv[1:]:
    sample_id = fbase(fname)
    d = pd.read_table(fname)
    cnarr = CNA(sample_id, d['Chr'], d['OriStCoordinate'],
                d['OriEndCoordinate'], d['Gene.Sym'],
                d['Adjusted.Mean.of.LogRatio'])
    cnarr.sort()
    cnarr.write(sample_id + ".contra.cnr")

