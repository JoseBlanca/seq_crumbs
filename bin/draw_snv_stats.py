#!/usr/bin/env python
'''
Created on 2013 urr 3

@author: peio
'''

import argparse
import os

from vcf_crumbs.vcf_stats import calc_density_per_chrom, get_data_from_vcf
from bam_crumbs.plot import draw_histogram, draw_scatter


def _setup_argparse():
    'It prepares the command line argument parsing.'
    description = 'Filter the snvs in an indexed vcf'
    parser = argparse.ArgumentParser(description=description)
    in_help = 'Indexed vcf file'
    parser.add_argument('input', help=in_help, type=argparse.FileType('rt'))
    parser.add_argument('-r', '--reference_file', type=argparse.FileType('rt'),
                        help="Reference file used to find the snvs",
                        required=True)
    parser.add_argument('-o', '--out_dir', default=os.getcwd(),
                        help='Directory to write all output files')

    return parser


def _parse_args(parser):
    '''It parses the command line and it returns a dict with the arguments.'''
    parsed_args = parser.parse_args()
    args = {}
    args['in_fhand'] = parsed_args.input
    args['out_dir'] = parsed_args.out_dir
    if not os.path.exists(args['out_dir']):
        os.mkdir(os.path.abspath(args['out_dir']))
    args['ref_fhand'] = parsed_args.reference_file
    return args


def main():
    parser = _setup_argparse()
    args = _parse_args(parser)
    out_dir = args['out_dir']
    ref_fhand = args['ref_fhand']

    data = get_data_from_vcf(args['in_fhand'].name)

    #densities
    densities = calc_density_per_chrom(data['snps_per_chromo'], ref_fhand,
                                       size=100)
    density_fhand = open(os.path.join(out_dir, 'distrib_snv_per_100_pb.png'),
                         'w')
    kwargs = {'title': 'Distribution of snps per 100 pb',
              'xlabel': 'Snvs per 100 pb', 'ylabel': 'Num. Seqs'}
    draw_histogram(densities.values(), density_fhand, bins=100, **kwargs)

    #mafs
    maf_fhand = open(os.path.join(out_dir, 'distrib_maf_per_snv.png'), 'w')
    kwargs = {'title': 'Maf per snv distribution', 'xlim': (0.5, 1),
              'xlabel': 'Maf', 'ylabel': 'Num. Snvs'}
    draw_histogram(data['maf_per_snp'], maf_fhand, bins=20, **kwargs)

    # call_scatter
    scatter_fhand = open(os.path.join(out_dir, 'call_scatter.png'), 'w')
    kwargs = {'title': 'Ref/alt allele count scatter',
              'ylabel': 'Alternative allele count',
              'xlabel': 'Reference_allele count', 'xlim': 0, 'ylim': 0}

    draw_scatter(data['call_data'].values(), scatter_fhand, **kwargs)


if __name__ == '__main__':
    main()
