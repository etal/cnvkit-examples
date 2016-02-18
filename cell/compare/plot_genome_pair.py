#!/usr/bin/env python

"""Plot two whole-genome CNV profiles, vertically stacked."""
from __future__ import print_function

import sys

import numpy
from matplotlib import pyplot

import cnvlib
from cnvlib import commands, plots

PAD = 1e7

def limit(val, mini, maxi):
    return min(max(val, mini), maxi)


def main(args):
    """."""
    # Load data
    cnarr = cnvlib.read(args.cnr_fname)
    segarr = cnvlib.read(args.cns_fname)
    acgharr = cnvlib.read(args.cghr_fname)
    asegarr = cnvlib.read(args.cghs_fname)

    # Create a figure grid w/ 2 axes, vertically stacked, labels sandwiched
    _fig = pyplot.figure(figsize=(10, 3.5))
    axgrid = pyplot.GridSpec(2, 1, hspace=.37)
    topax = pyplot.subplot(axgrid[0])
    botax = pyplot.subplot(axgrid[1], sharex=topax, sharey=topax)
    botax.tick_params(labelbottom=False)
    topax.tick_params(labelbottom=True)
    # Twiddle y-axis limits
    all_y = numpy.concatenate((segarr.autosomes().log2, asegarr.autosomes().log2))
    topax.set_ylim(limit(min(all_y) - .2, -5.0, -.5),
                   limit(max(all_y) + .2, .5, 5.0))

    # Draw CNVkit and aCGH scatters
    plots.cnv_on_genome(botax, acgharr, asegarr, PAD)
    plots.cnv_on_genome(topax, cnarr, segarr, PAD)

    # Save it.
    if args.output:
        pyplot.savefig(args.output, format='pdf', bbox_inches=0)
        print("Wrote", args.output, file=sys.stderr)
    else:
        pyplot.show()


if __name__ == '__main__':
    import argparse
    AP = argparse.ArgumentParser(description=__doc__)
    AP.add_argument('cnr_fname', help="CNVkit .cnr filename")
    AP.add_argument('cns_fname', help="CNVkit .cns filename")
    AP.add_argument('cghr_fname', help="aCGH .cnr filename")
    AP.add_argument('cghs_fname', help="aCGH .cns filename")
    AP.add_argument('-o', '--output', help="Output PDF filename")
    main(AP.parse_args())

