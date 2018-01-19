#!/usr/bin/env python
"""Convert TCGA SEG-formatted files to CNVkit's .cns format.

Get gene names from the GeneInfo table.
"""
from __future__ import absolute_import, division, print_function
import argparse
import os

from cnvlib.rna import load_gene_info
from skgenome import tabio, GenomicArray as GA


def basename(path):
    return os.path.basename(path).split('.', 1)[0]


def join_unique(vals):
    return ','.join(vals.unique())


def genes_in_segments(segarr, gene_info):
    return gene_info.into_ranges(segarr, 'gene', '-', join_unique)


if __name__ == '__main__':
    AP = argparse.ArgumentParser(description=__doc__)
    AP.add_argument('seg_files', nargs='+',
                    help="Segmented aCGH data in SEG format.")
    AP.add_argument('-g', '--gene-resource', metavar="FILE", required=True,
                    # default="data/ensembl-gene-info.hg38.tsv",
                    help="Ensembl BioMart-derived gene info table.")
    AP.add_argument('-d', '--output-dir', metavar='PATH', default='.',
                    help="Output directory.")

    args = AP.parse_args()
    gene_info = load_gene_info(args.gene_resource, None, None)
    bad_genes = ['Metazoa_SRP', '5S_rRNA', 'Y_RNA', 'U1', 'U2', 'U3', 'U4',
                 'U5', 'U6', 'U7', 'U8', 'uc_338', 'Clostridiales-1']
    gene_info = gene_info[~gene_info['gene'] .isin(bad_genes)]
    gene_info = GA(gene_info.loc[:, ('chromosome', 'start', 'end', 'gene')])

    for seg_fname in args.seg_files:
        seg = tabio.read(seg_fname, 'seg')
        # Assign gene names to segments using genomic coordinates from gene_info
        seg['gene'] = genes_in_segments(seg, gene_info)

        outfname = os.path.join(args.output_dir,
                                basename(seg_fname) + ".acgh.cns")
        tabio.write(seg, outfname, 'tab')
        print("Wrote", outfname)
