#!/usr/bin/env python

"""Plot raw coverages vs. various confounding properties for multiple files."""
from __future__ import division, print_function

import argparse
from itertools import izip
from os.path import basename

import numpy as np
import seaborn
from matplotlib import pyplot, cm

from cnvlib import read
from cnvlib import fix, params
from cnvlib.core import shift_xx
from cnvlib.ngfrills import echo
# from cnvlib.reference import mask_bad_probes
from cnvlib.smoothing import rolling_median, smoothed


seaborn.set(font='Sans', style="ticks")


def load_cna(fname, reference):
    """Read CNA, adjust gender. Subtract reference if given (for ratio)."""
    echo("Processing", fname)
    cnarr = read(fname)
    if reference:
        # Subtract the reference copy number values (to get the log2 ratio)
        cnarr = fix.load_adjust_coverages(cnarr, reference, False, False, False)
        cnarr = shift_xx(cnarr, male_normal=True)
    else:
        cnarr = shift_xx(cnarr, male_normal=True)
        # Drop low-coverage probes (otherwise done in load_adjust_coverages)
        # cnarr = cnarr.to_rows(
        #     cnarr[cnarr.coverage >= params.MIN_BIN_COVERAGE])
    return cnarr


def get_bias_func(mode, ref_pset, probes):
    if not ref_pset:
        raise ValueError("Must supply a reference for " + mode)
    ref_matched = fix.match_ref_to_probes(ref_pset, probes)

    if mode in ('gc', 'rmask'):
        return ref_matched[mode]
    elif mode == 'edge':
        return fix.make_edge_sorter(ref_matched, params.INSERT_SIZE)
    else:
        raise ValueError("Unknown mode: %s" % mode)


def sort_and_smooth(probes, biases):
    if callable(biases):
        biases = map(biases, probes)
    biases, coverages = zip(*sorted(
        ((bias, cvg) for bias, cvg in izip(biases, probes['coverage'])),
        key=lambda bc: bc[0]))
    # Smooth the biases
    cvg_fitted = rolling_median(coverages, .2)
    # Again! (for aesthetics)
    # cvg_fitted = smoothed(cvg_fitted, .05)

    # Print some stats
    coverages = np.asarray(coverages)
    orig_var = np.var(coverages)
    def improvement(fitvals):
        return 100 * (1 - (np.var(coverages - fitvals) / orig_var))

    # print("Sample \tRaw probes \tTrend line \tReduction")
    print(probes.sample_id,
          "\t %.5f    \t %.5f    \t %.4f"
          % (orig_var, np.var(cvg_fitted), improvement(cvg_fitted)),
          '%')
    return biases, coverages, cvg_fitted


def get_sort_and_smoother(cna_fname, ref_arr, mode):
    """Make a sort_and_smooth func from example CNA and reference."""
    ref_matched = fix.match_ref_to_probes(ref_arr, read(cna_fname))

    if mode in ('gc', 'rmask'):
        biases = ref_matched[mode]
    elif mode == 'edge':
        biases = map(fix.make_edge_sorter(ref_matched, params.INSERT_SIZE),
                     ref_arr)
    else:
        raise ValueError("Unknown mode: %s" % mode)


    def wrapped_sort_and_smooth(this_arr):
        """Sort and smooth."""
        assert len(this_arr) == len(biases)

        biases, coverages = zip(*sorted(
            ((bias, cvg) for bias, cvg in izip(biases, this_arr['coverage'])),
            key=lambda bc: bc[0]))
        # Smooth the biases
        cvg_fitted = rolling_median(coverages, .2)
        # Again! (for aesthetics)
        # cvg_fitted = smoothed(cvg_fitted, .05)

        # Print some stats
        coverages = np.asarray(coverages)
        orig_var = np.var(coverages)
        def improvement(fitvals):
            return 100 * (1 - (np.var(coverages - fitvals) / orig_var))

        # print("Sample \tRaw probes \tTrend line \tReduction")
        print(this_arr.sample_id,
            "\t %.5f    \t %.5f    \t %.4f"
            % (orig_var, np.var(cvg_fitted), improvement(cvg_fitted)))
        return biases, coverages, cvg_fitted

    return wrapped_sort_and_smooth


def plot_separate(filenames, ref_pset, bias_func, mode, do_ratio):
    """Plot coverages versus other factors to reveal systematic biases."""
    _fig, axes = pyplot.subplots(len(filenames), squeeze=False, sharex=True,
                                 figsize=(4, 4))

    for fname, ax in zip(filenames, axes[:, 0]):
        # Compute points to plot
        pset = load_cna(fname, ref_pset if do_ratio else None)
        bias, coverages, fitted = sort_and_smooth(pset, bias_func)
        ax.plot(bias, fitted, color='#F04040', alpha=0.7, lw=2, zorder=-.1)
        ax.scatter(bias, coverages, marker='.', color='#666666', zorder=-1, alpha=0.1)
        # Aesthetics
        if mode == 'edge':
            ax.set_xlim(xmin=-1, xmax=0)
        else:
            ax.set_xlim(xmin=0, xmax=1)
        # ax.set_ylim(ymin=min(bias) - .5, ymax=max(bias) + .5)
        ax.set_ylim(ymin=-1.1, ymax=1.1)
        ax.axhline(color='k', linestyle='-', zorder=-2)
        ax.set_title(basename(fname))
        ax.set_ylabel("Copy ratio (log2)" if do_ratio else "Copy number (log2)")

    pyplot.xlabel(mode)


def plot_overlaid(filenames, ref_pset, bias_func, mode, do_ratio, colorscheme):
    """Plot coverages versus other factors to reveal systematic biases."""
    _fig, ax = pyplot.subplots(figsize=(4, 4))
    ax.tick_params(labelsize='large')

    colors = map(getattr(cm, colorscheme), np.arange(.1, .9, .8/len(filenames)))
    for fname, color in zip(filenames, colors):
        # Compute points to plot
        pset = load_cna(fname, ref_pset if do_ratio else None)
        bias, _coverages, fitted = sort_and_smooth(pset, bias_func)
        ax.plot(bias, fitted, color=color, alpha=0.7, lw=2, zorder=-.1)

    # Aesthetics
    if mode == 'edge':
        ax.set_xlim(xmin=-1, xmax=0)
    else:
        ax.set_xlim(xmin=0, xmax=1)
    # ax.set_ylim(ymin=-1.5, ymax=1.5)
    ax.axhline(color='k', linestyle=':', zorder=-2)
    ax.set_xlabel(mode)
    ax.set_ylabel("Copy ratio (log2)" if do_ratio else "Copy number (log2)")


def main(args):
    """*"""
    do_ratio = bool(args.reference)
    ref_pset = read(args.reference or args.no_reference)
    bias_func = get_bias_func(args.mode, ref_pset, read(args.filenames[0]))

    print("Sample \tRaw probes \tTrend line \tReduction (%)")
    if args.batch:
        plot_overlaid(args.filenames, ref_pset, bias_func, args.mode, do_ratio, args.color)
    else:
        plot_separate(args.filenames, ref_pset, bias_func, args.mode, do_ratio)

    if args.output:
        pyplot.savefig(args.output, format='pdf', bbox_inches=0)
        echo("Wrote", args.output)
    else:
        pyplot.show()


if __name__ == '__main__':
    AP = argparse.ArgumentParser(
            description=__doc__,
            epilog="Contact Eric Talevich <etalevich@derm.ucsf.edu> for help.")
    AP.add_argument("filenames", nargs='+',
                    help="Raw coverage filenames (.cnn).")
    AP.add_argument('-m', "--mode",
                    default='gc', choices=['gc', 'edge', 'rmask'],
                    help="Type of bias to examine.")
    AP.add_argument('-b', '--batch', action='store_true',
                    help="Plot only the overlaid trendlines of the arguments.")
    AP.add_argument('-c', '--color', default="Reds",
                    help="matplotlib colormap name.")
    AP.add_argument('-r', '--reference',
                    help="Reference coverage table (to compute CN ratios).")
    AP.add_argument('-nr', '--no-reference',
                    help="""Reference coverage table for GC, but NOT to compute
                    CN ratios -- show raw copy numbers).""")
    AP.add_argument("-o", "--output", help="Output PDF file name")
    main(AP.parse_args())
