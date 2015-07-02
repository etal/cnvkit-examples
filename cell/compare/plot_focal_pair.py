#!/usr/bin/env python

"""Plot two focal CNV profiles, horizontally stacked."""
from __future__ import print_function

import sys

from matplotlib import pyplot
import numpy

import cnvlib
from cnvlib import commands, plots


def get_plot_args(cnarr, segarr, chrom, window_coords):
    sel_probes = cnarr.in_range(chrom, *window_coords)
    print("Selected", len(sel_probes), "probes")
    sel_segs = segarr.in_range(chrom, *window_coords, trim=True)
    print("sel_segs:\n", sel_segs.data)
    return (sel_probes, sel_segs)


def get_xtick_values(probes):
    """FML"""
    min_x = probes.start.min() * 1e-6
    max_x = probes.end.max() * 1e-6
    min_tick = int(min_x)
    max_tick = int(max_x) + 1
    if max_x - min_x < 2.5:
        if min_tick < min_x:
            min_tick += .5
        if min_tick < min_x:
            min_tick += .5
        if max_tick > max_x:
            max_tick -= .5
        if max_tick > max_x:
            max_tick -= .5
        return numpy.arange(min_tick, max_tick + .5, .5)
    else:
        if min_tick < min_x:
            min_tick += 1
        if max_tick > max_x:
            max_tick -= 1
        return range(min_tick, max_tick + 1)


def main(args):
    """."""
    # Load data
    cnarr = cnvlib.read(args.cnr_fname)
    # cnarr['weight'] = numpy.repeat(.78, len(cnarr))
    segarr = cnvlib.read(args.cns_fname)
    acgharr = cnvlib.read(args.cghr_fname)
    asegarr = cnvlib.read(args.cghs_fname)

    # Find the genomic location matching the specified gene(s)
    gene_names = args.gene_name.split(',')
    gene_coords = plots.gene_coords_by_name(cnarr, gene_names)
    if not len(gene_coords) == 1:
        raise ValueError("Genes %s are split across chromosomes %s"
                         % (args.gene_name, gene_coords.keys()))
    chrom, genes = gene_coords.popitem()
    genes.sort()
    # Set the display window to the selected genes +/- a margin
    window_coords = (genes[0][0] - args.window_width,
                     genes[-1][1] + args.window_width)

    # Use plot_chromosome to draw CNVkit and aCGH scatters
    cnv_sel_probes, cnv_sel_segs = get_plot_args(cnarr, segarr,
                                                 chrom, window_coords)
    acgh_sel_probes, acgh_sel_segs = get_plot_args(acgharr, asegarr,
                                                   chrom, window_coords)

    # Create a figure grid w/ 2 side-by-side axes
    _fig = pyplot.figure(figsize=(3.5 * len(genes), 3.5))
    axgrid = pyplot.GridSpec(1, 2, wspace=0)
    leftax = pyplot.subplot(axgrid[0])
    rightax = pyplot.subplot(axgrid[1], sharex=leftax, sharey=leftax)

    plots.plot_chromosome(leftax, cnv_sel_probes, cnv_sel_segs,
                          chrom, cnarr.sample_id, genes)
    plots.plot_chromosome(rightax, acgh_sel_probes, acgh_sel_segs,
                          chrom, cnarr.sample_id, genes)

    # Tweak aesthetics
    rightax.tick_params(labelleft=False, left=False)
    leftax.tick_params(labelleft=True)
    leftax.set_xlabel("Position (Mb)")
    rightax.set_ylabel('')
    rightax.set_title('')
    # Rotate & cull x-axis (position) labels
    if len(genes) == 1:
        xlabels = get_xtick_values(acgh_sel_probes)
        leftax.set_xticks(xlabels)
        leftax.set_xticklabels(map(str, xlabels), rotation=60)
        rightax.set_xticks(xlabels)
        rightax.set_xticklabels(map(str, xlabels), rotation=60)

    # Set sensible y-axis limits
    # all_y = numpy.concatenate((cnv_sel_probes.coverage,
    #                            acgh_sel_probes.coverage))
    # leftax.set_ylim(plots.limit(min(all_y) - .1, -5.0, -.3),
    #                 plots.limit(max(all_y) + .25, .3, 5.0))
    if args.gene_name == 'CDKN2A':
        all_y = numpy.concatenate((cnv_sel_segs.coverage,
                                   acgh_sel_segs.coverage))
        print("all_y:", tuple(all_y))
        leftax.set_ylim(plots.limit(min(all_y) - .3, -5.0, -.5),
                        plots.limit(max(all_y) + .3, .5, 5.0))
    else:
        leftax.set_ylim(-2.1, 1.1)

    # Save it.
    if args.output:
        pyplot.savefig(args.output, format='pdf', bbox_inches='tight')
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
    AP.add_argument('gene_name', help="Name of the gene to show")
    AP.add_argument('-w', '--window-width', type=float, default=7e5,
                    help="Size of margin flanking the gene to show")
    AP.add_argument('-o', '--output', help="Output PDF filename")
    main(AP.parse_args())
