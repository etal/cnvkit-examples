#!/usr/bin/env python
"""Plot the relationship (systematic bias) between counts and GC or tx length.

CNV-RNA bias plot:
- want plot of counts vs. GC for each gene
- fit a curve/trendline (e.g. rolling median)
- calculate % variance explained by the trendline: var(before)- var(after)


"""
from __future__ import absolute_import, division, print_function
import os
import sys
from glob import glob

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

from cnvlib import rna, import_rna as IR
from cnvlib.rna import load_gene_info
from cnvlib.smoothing import savgol, rolling_median

from skgenome import tabio

def do_import_rna(gene_count_fnames):
    """From cnvlib.import_rna.do_import_rna, minus the bias corrections."""
    # Deduplicate and ensure all normals are included in the analysis
    gene_count_fnames = sorted(set(gene_count_fnames))
    sample_counts = IR.aggregate_gene_counts(gene_count_fnames)
    sample_counts = rna.filter_probes(sample_counts)

    correlations_fname = "../../cnvkit/data/tcga-skcm.cnv-expr-corr.tsv"
    gene_resource_fname = "../../cnvkit/data/ensembl-gene-info.hg38.tsv"

    print("Loading gene metadata" +
          (" and TCGA gene expression/CNV profiles"
           if correlations_fname else ""))
    gene_info = rna.load_gene_info(gene_resource_fname, correlations_fname)

    print("Aligning gene info to sample gene counts")
    gene_info, sample_counts, sample_data_log2 = rna.align_gene_info_to_samples(
        gene_info, sample_counts, None, [])

    # Summary table has log2-normalized values, not raw counts
    # ENH show both, with column header suffixes to distinguish?
    all_data = pd.concat([gene_info, sample_data_log2], axis=1)
    # CNVkit files have both absolute and log2-normalized read counts
    cnrs = rna.attach_gene_info_to_cnr(sample_counts, sample_data_log2,
                                       gene_info)
    #  cnrs = (rna.correct_cnr(cnr) for cnr in cnrs) # XXX
    return all_data, cnrs


if __name__ == '__main__':

    fnames = sys.argv[1:] or glob("tcga-rna-counts/*.txt")
    _, cnrs = do_import_rna(fnames)

    for cnr in cnrs:
        # TODO calculate variance (or weighted SD) before & after
        d = (cnr.data.loc[:, ('gc', 'log2', 'weight')]
             .sort_values(by='gc')
             .reset_index(drop=True))

        sizes = 46 * d['weight'] ** 2 + 2
        plt.scatter(d['gc'], d['log2'], s=sizes,
                    color="#606060", edgecolor="none", alpha=.1)

        # Show the bias being corrected 
        trendline = rolling_median(d['log2'].values, .1)
        plt.plot(d['gc'], trendline,
                 color='darkorange', linewidth=2, zorder=2)

        plt.ylim(-2, 2)
        plt.xlim(.2, .8)
        out_fname = cnr.sample_id + ".gc-bias.png"
        plt.savefig(out_fname, format="png", bbox_inches="tight")
        print("Wrote", out_fname)
        plt.close()

        # plot vs. each column in d
