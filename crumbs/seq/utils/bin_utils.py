# Copyright 2012 Jose Blanca, Peio Ziarsolo, COMAV-Univ. Politecnica Valencia
# This file is part of seq_crumbs.
# seq_crumbs is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.

# seq_crumbs is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR  PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with seq_crumbs. If not, see <http://www.gnu.org/licenses/>.

import sys
import argparse

from crumbs.utils.file_utils import (wrap_in_buffered_reader,
                                     uncompress_if_required, compress_fhand)
from crumbs.utils.tags import OUTFILE, GUESS_FORMAT
from crumbs.seq.utils.file_formats import get_format, set_format
from crumbs.utils.bin_utils import build_version_msg, get_requested_compression


def create_basic_argparse(**kwargs):
    'It returns a parser with several inputs and one output'
    parser = argparse.ArgumentParser(**kwargs)

    parser.add_argument('input', default=sys.stdin, nargs='*',
                        help='Sequence input files to process (default STDIN)',
                        type=argparse.FileType('rt'))

    hlp_fmt = 'Format of the input files (default: %(default)s)'
    parser.add_argument('-t', '--in_format', default=GUESS_FORMAT,
                        help=hlp_fmt)

    parser.add_argument('-o', '--outfile', default=sys.stdout, dest=OUTFILE,
                        help='Sequence output file (default: STDOUT)',
                        type=argparse.FileType('wt'))

    parser.add_argument('--version', action='version',
                        version=build_version_msg())

    group = parser.add_mutually_exclusive_group()
    group.add_argument('-z ', '--gzip', action='store_true',
                       help='Compress the output in gzip format')
    group.add_argument('-Z ', '--bgzf', action='store_true',
                       help='Compress the output in bgzf format')
    group.add_argument('-B ', '--bzip2', action='store_true',
                       help='Compress the output in bzip2 format')
    return parser


def create_basic_parallel_argparse(**kwargs):
    'It returns a cmd parser with inputs, output, format, num_processes'
    parser = create_basic_argparse(**kwargs)
    parser = argparse.ArgumentParser(parents=[parser], add_help=False)
    parser.add_argument('-p', '--processes', dest='processes', type=int,
                        help='Num. of processes to use (default: %(default)s)',
                        default=1)
    return parser


def _to_bool(string):
    if string.lower()[0] == 'f':
        return False
    elif string.lower()[0] == 't':
        return True
    elif string.isdigit():
        return bool(int(string))


def create_filter_argparse(add_reverse=True, **kwargs):
    'It returns a cmd parser for the filter executables'
    parser = create_basic_parallel_argparse(**kwargs)
    parser = argparse.ArgumentParser(parents=[parser], add_help=False)
    if add_reverse:
        parser.add_argument('-r', '--reverse', action='store_true',
                            help='Reverses the filtering')
    parser.add_argument('-e', '--filtered_file',
                        help='Filtered out sequences output file',
                        type=argparse.FileType('wt'))
    group = parser.add_argument_group('Pairing')
    group.add_argument('--paired_reads', action='store_true',
                       help='Filter considering interleaved pairs')
    help_msg = 'If one read fails the pair will be filtered out '
    help_msg += '(default: %(default)s)'
    group.add_argument('--fail_drags_pair', type=_to_bool, default='true',
                       choices=(True, False), help=help_msg)
    return parser


def create_trimmer_argparse(**kwargs):
    'It returns a cmd parser for the filter executables'
    parser = create_basic_parallel_argparse(**kwargs)
    parser = argparse.ArgumentParser(parents=[parser], add_help=False)
    parser.add_argument('-m', '--mask', dest='mask', action='store_true',
                        help='Do not trim, only mask by lowering the case')

    group = parser.add_argument_group('Pairing')
    group.add_argument('--paired_reads', action='store_true',
                       help='Trim considering interleaved pairs')
    group.add_argument('-e', '--orphan_file',
                       help='Orphan sequences output file',
                       type=argparse.FileType('wt'))
    return parser


def parse_basic_args(parser):
    'It parses the command line and it returns a dict with the arguments.'
    parsed_args = parser.parse_args()
    # we have to wrap the file in a BufferedReader to allow peeking into stdin
    wrapped_fhands = []
    # if input is stdin it will be a fhand not a list of fhands.
    # we have to convert to a list
    in_fhands = parsed_args.input
    if not isinstance(in_fhands, list):
        in_fhands = [in_fhands]
    for fhand in in_fhands:
        fhand = wrap_in_buffered_reader(fhand)
        fhand = uncompress_if_required(fhand)
        wrapped_fhands.append(fhand)

    # We have to add the one_line to the fastq files in order to get the
    # speed improvements of the seqitems
    in_format = parsed_args.in_format
    if in_format == GUESS_FORMAT:
        for wrapped_fhand in wrapped_fhands:
            get_format(wrapped_fhand)
    else:
        for wrapped_fhand in wrapped_fhands:
            set_format(wrapped_fhand, in_format)

    out_fhand = getattr(parsed_args, OUTFILE)

    comp_kind = get_requested_compression(parsed_args)
    if isinstance(out_fhand, list):
        new_out_fhands = []
        for out_f in out_fhand:
            try:
                out_f = compress_fhand(out_f, compression_kind=comp_kind)
            except RuntimeError, error:
                parser.error(error)

            new_out_fhands.append(out_f)
        out_fhand = new_out_fhands
    else:
        try:
            out_fhand = compress_fhand(out_fhand, compression_kind=comp_kind)
        except RuntimeError, error:
            parser.error(error)

    # The default output format is the same as the first file
    if 'fastq' in in_format or in_format == GUESS_FORMAT:
        out_format = get_format(wrapped_fhands[0])
    else:
        out_format = in_format

    # The original fhands should be stored, because otherwise they would be
    # closed
    args = {'out_fhand': out_fhand, 'in_fhands': wrapped_fhands,
            'out_format': out_format, 'original_in_fhands': in_fhands}
    return args, parsed_args


def parse_basic_parallel_args(parser):
    'It parses the command line and it returns a dict with the arguments.'
    args, parsed_args = parse_basic_args(parser)
    args['processes'] = parsed_args.processes
    return args, parsed_args


def parse_filter_args(parser, add_reverse=True):
    'It parses the command line and it returns a dict with the arguments.'
    args, parsed_args = parse_basic_parallel_args(parser)
    if add_reverse:
        args['reverse'] = parsed_args.reverse
    args['filtered_fhand'] = parsed_args.filtered_file
    paired_reads = parsed_args.paired_reads
    args['paired_reads'] = paired_reads
    if paired_reads:
        # in this case fail_drags_pair is required
        fail_drags_pair = parsed_args.fail_drags_pair
        if fail_drags_pair is None:
            msg = 'For pairs fail_drags_pair is required'
            parser.error(msg)
            # raise argparse.ArgumentError(parsed_args.fail_drags_pair, msg)
    else:
        fail_drags_pair = None
    args['fail_drags_pair'] = fail_drags_pair

    return args, parsed_args


def parse_trimmer_args(parser):
    'It parses the command line and it returns a dict with the arguments.'
    args, parsed_args = parse_basic_parallel_args(parser)
    args['mask'] = parsed_args.mask
    args['orphan_fhand'] = parsed_args.orphan_file
    paired_reads = parsed_args.paired_reads
    args['paired_reads'] = paired_reads
    return args, parsed_args
