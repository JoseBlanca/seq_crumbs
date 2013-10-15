#!/usr/bin/env python
'''
Created on 2013 uzt 16

@author: peio
'''

import argparse
import sys
import itertools

from vcf import Reader


def _setup_argparse():
    'It prepares the command line argument parsing.'
    description = 'Filter the snvs in an indexed vcf'
    parser = argparse.ArgumentParser(description=description)
    in_help = 'Indexed vcf file'
    parser.add_argument('input', help=in_help, type=argparse.FileType('rt'))
    parser.add_argument('-o', '--output',
                        help='output vcf file', type=argparse.FileType('wt'))
    parser.add_argument('-f', '--filters', action='append', required=True,
                        help='Filters you want to calculate')

    return parser


def _parse_args(parser):
    '''It parses the command line and it returns a dict with the arguments.'''
    parsed_args = parser.parse_args()
    args = {}
    args['in_fhand'] = parsed_args.input
    if parsed_args.output is not None:
        args['out_fhand'] = parsed_args.output
    else:
        args['out_fhand'] = sys.stdout
    args['filters'] = parsed_args.filters
    return args


def venndiv2(*sets):
    bins = tuple(set() for n in range(2 ** len(sets)))
    for i, si in enumerate(sets):  # Over every individual item of each set
        for item in si:
            binindex = 0
            for j, sj in enumerate(sets):  # Test for membership of each set
                if item in sj:
                    binindex += 2 ** j  # Binary count binning
            bins[binindex].add(item)
    return bins[1:]  # bins[0] - an item in no set is impossible.


def venndiv(*sets):
    def vd(x, *xs):
        bins = vd(*xs) if xs else (set.union(*sets),)
        return sum(tuple((bin - x, bin & x) for bin in bins), ())
    return vd(*sets)[1:]


def write_venn_headers(headers):
    #     'p u m\ne p u\np v 1\no 1 6\n-----\n')
    if not headers:
        raise RuntimeError('No headers')
    headers = headers[::-1]
    max_ = max([len(h) for h in headers])
    headers_str = []
    for index in range(0, max_)[::-1]:
        line = []
        for header in headers:
            try:
                letter = header[index]
            except IndexError:
                letter = ' '
            line.append(letter)
        line = " ".join(line)
        headers_str.append(line)
    headers_str = headers_str[::-1]
    headers_str.append('-' * len(headers_str[-1]))
    return '+' * len(headers_str[-1]) + '\n' + '\n'.join(headers_str) + '\n'


def main():
    parser = _setup_argparse()
    args = _parse_args(parser)
    w_filters = args['filters']
    out_fhand = args['out_fhand']
    reader = Reader(args['in_fhand'])
    for w_filter in w_filters:
        vcf_filter = reader.filters[w_filter]
        out_fhand.write('{} : {}\n'.format(vcf_filter.id, vcf_filter.desc))
    out_fhand.write(write_venn_headers(w_filters))
    codes = [" ".join(seq) for seq in itertools.product("01", repeat=len(w_filters))][1:]
    filter_stats = {f: set()for f in w_filters}
    count = 0
    for index, record in enumerate(reader):
#         if count == 1000:
#             break
        filters = record.FILTER
        record_id = index
        for w_filter in w_filters:
            if w_filter not in filters:
                filter_stats[w_filter].add(record_id)
        count += 1
    filt_result = [filter_stats[f] for f in w_filters]

    sets = venndiv(*filt_result)
    for code, set_ in zip(codes, sets):
        out_fhand.write('{}: {}\n'.format(code, len(set_)))


if __name__ == '__main__':
    main()
