'''
Created on 2012 eka 25

@author: peio
'''

import argparse

from crumbs.iterutils import sample, length
from crumbs.utils import (STDIN, STDOUT, get_inputs_from_args,
                                get_output_from_args, INFILES, OUTFILE)
from crumbs.seqio import (guess_format, write_seqrecords, read_seqrecords,
                          count_seqs_in_files)


def sample_setup_argparse():
    'It prepares the command line argument parsing.'
    description = 'Get first seqs from file[s]'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(INFILES, default=STDIN, nargs='*',
                        help="Sequence input file to process")
    parser.add_argument('-o', '--outfile', default=STDOUT, dest=OUTFILE,
                        help="Sequence output file to process")
    parser.add_argument('-n', '--num_seqs', default=10, type=int,
                        dest='num_seqs', help='Number of sequendes to print')

    return parser


def sample_parse_args(parser):
    '''It parses the command line and it returns a dict with the arguments.'''
    parsed_args = parser.parse_args()
    num_seqs = parsed_args.num_seqs
    in_fhands = get_inputs_from_args(parsed_args)
    out_fhand = get_output_from_args(parsed_args)

    return {'out_fhand': out_fhand, 'in_fhands': in_fhands,
            'num_seqs': num_seqs}


def get_random_seqs(fhands, num_seqs):
    'it gets the random seqs'
    total_seqs = count_seqs_in_files(fhands)
    seqs = read_seqrecords(fhands)
    return sample(seqs, total_seqs, num_seqs)

