#!/usr/bin/env python

"""Calculate precision/recall for aCGH vs. HTS.

Split CNVs by:
    - tool
    - size
    - CNV type (gain/loss)

Choose 3 threshold for each CNV type, e.g. if neutral = 6 copies, calculate
gains >= 7, 8, 9; losses <= 5, 4, 3.

From the wao/wbo files, hits are >= 50% of bases in the reference region.


"""
from __future__ import division, print_function

# import sys

import pandas as pd
import numpy as np

NEUTRAL = 6
LARGE = 5e6


def load_wxo(bedfname):
    table = pd.read_table(bedfname,
                          names=["ref_chrom", "ref_start", "ref_end",
                                 "ref_sample", "ref_cn",
                                 "alt_chrom", "alt_start", "alt_end",
                                 "alt_sample", "alt_cn", "nbases"],
                          na_values=".")
    table["ref_size"] = table.ref_end - table.ref_start
    return table


def split_by_size(table):
    is_large = (table.ref_size >= LARGE)
    return table[is_large], table[~is_large]


def count_hits(table, copynum, is_gain, fraction=.5):
    """Count all positives and true positives."""
    if is_gain:
        table = table[table.ref_cn >= copynum]
    else:
        table = table[table.ref_cn <= copynum]
    denom = len(table)  # All positives

    if is_gain:
        is_alt_pos = (table.alt_cn >= copynum)
    else:
        is_alt_pos = (table.alt_cn <= copynum)
    is_overlap = ((table.nbases / table.ref_size) >= fraction)
    num = (is_overlap & is_alt_pos).sum()
    return num, denom


def count_bp_hits(table, copynum, is_gain):
    if is_gain:
        table = table[table.ref_cn >= copynum]
    else:
        table = table[table.ref_cn <= copynum]
    denom = table.ref_size.sum()  # All positives

    if is_gain:
        table = table[table.alt_cn >= copynum]
    else:
        table = table[table.alt_cn <= copynum]
    num = table.nbases.sum()
    return num, denom


def all_precision_recall(control_table, test_table):
    """Calculate precision and recall for a pair of tables.

    When a search engine returns 30 pages only 20 of which were relevant while
    failing to return 40 additional relevant pages, its precision is 20/30 = 2/3
    while its recall is 20/60 = 1/3.
            (at what threshold?)

    Total found = 30 (t_total)
    Total relevant = 60 (c_total)
    N true pos. = 20 (t_hits, or c_hits)
    N missed = 40 (c_total - t_hits or c_hits)

    Precision = n_true_pos / t_total
    Recall = n_true_pos / c_total

    """
    for cn in range(NEUTRAL - 2, NEUTRAL) + range(NEUTRAL + 1, NEUTRAL + 3):
        is_gain = (cn > NEUTRAL)
        # seen_pr = set()
        # for frac in np.arange(.02, 1.0, .02):
        frac = .5
        c_hits, c_total = count_hits(control_table, cn, is_gain, frac)
        t_hits, t_total = count_hits(test_table, cn, is_gain, frac)
        precision = t_hits / t_total if t_total else np.nan
        recall = c_hits / c_total if c_total else np.nan
        if precision or recall: # and (precision, recall) not in seen_pr:
            yield (cn, is_gain, c_hits, c_total, t_hits, t_total,
                    precision, recall, frac)
        # seen_pr.add((precision, recall))


def bp_precision_recall(control_table, test_table):
    """Calculate p/r by base pair."""
    for cn in range(NEUTRAL - 2, NEUTRAL) + range(NEUTRAL + 1, NEUTRAL + 3):
        is_gain = (cn > NEUTRAL)
        c_hits, c_total = count_bp_hits(control_table, cn, is_gain)
        t_hits, t_total = count_bp_hits(test_table, cn, is_gain)
        precision = t_hits / t_total if t_total else np.nan
        recall = c_hits / c_total if c_total else np.nan
        if precision or recall: # and (precision, recall) not in seen_pr:
            yield (cn, is_gain, c_hits, c_total, t_hits, t_total,
                    precision, recall, 1)


def enframe_pr(rows, size):
    df = pd.DataFrame.from_records(list(rows),
                                   columns=["CN", "IsGain",
                                            "CtrlHits", "CtrlTotal",
                                            "TestHits", "TestTotal",
                                            "Precision", "Recall",
                                            "FracOverlap"])
    df["Size"] = size
    return df


def main(args):
    wao, wbo = map(load_wxo, [args.wao, args.wbo])
    wao_large, wao_small = split_by_size(wao)
    wbo_large, wbo_small = split_by_size(wbo)
    results_all = enframe_pr(all_precision_recall(wao, wbo), "All")
    results_large = enframe_pr(all_precision_recall(wao_large, wbo_large),
                               "Large")
    results_small = enframe_pr(all_precision_recall(wao_small, wbo_small),
                               "Small")
    results_bp = enframe_pr(bp_precision_recall(wao, wbo), "bp")
    results = pd.concat([results_all, results_large, results_small, results_bp],
                        ignore_index=True)
    results["Filename"] = args.wao.rsplit(".", 2)[0]
    results.to_csv(args.output, sep='\t', index=False)



if __name__ == '__main__':
    import argparse
    AP = argparse.ArgumentParser(description=__doc__)
    AP.add_argument("wao", help="aCGH as reference BED")
    AP.add_argument("wbo", help="Seq. as reference BED")
    AP.add_argument("-o", "--output", help="Output filename")
    main(AP.parse_args())
