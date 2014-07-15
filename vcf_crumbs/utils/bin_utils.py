'''
Created on 2014 uzt 15

@author: peio
'''
import argparse
import sys


def setup_basic_argparse(**kwargs):
    'It prepares the command line argument parsing.'

    parser = argparse.ArgumentParser(**kwargs)
    in_help = 'Input VCF file (default STDIN)'
    parser.add_argument('input', help=in_help)
    parser.add_argument('-o', '--output', default=sys.stdout,
                        help='Output VCF file (default STDOUT)',
                        type=argparse.FileType('w'))
    msg = 'Template VCF to get the header (default same as input)'
    parser.add_argument('-t', '--template', help=msg)
    msg = 'File to print some statistics (default STDERR)'
    parser.add_argument('-l', '--log', help=msg, type=argparse.FileType('r'),
                        default=sys.stderr)
    return parser


def setup_filter_argparse(**kwargs):
    'It prepares the command line argument parsing.'
    parser = setup_basic_argparse(**kwargs)
    parser.add_argument('-f', '--filtered',
                        help='Output for filtered SNVs',
                        type=argparse.FileType('w'))
    return parser


def parse_basic_args(parser):
    parsed_args = parser.parse_args()
    in_fpath = parsed_args.input
    template_fpath = parsed_args.template

    if in_fpath is None:
        in_fhand = sys.stdin
        if template_fpath is None:
            msg = 'The template file is required if the input is STDIN'
            parser.error(msg)
    else:
        in_fhand = open(in_fpath)

    if template_fpath is None:
        template_fhand = open(in_fpath)
    else:
        template_fhand = open(template_fhand)

    out_fhand = parsed_args.output
    log_fhand = parsed_args.log

    args = {'in_fhand': in_fhand, 'template_fhand': template_fhand,
            'log_fhand': log_fhand, 'out_fhand': out_fhand}

    return args, parsed_args


def parse_filter_args(parser):
    'It parses the command line and it returns a dict with the arguments.'
    args, parsed_args = parse_basic_args(parser)

    filtered_fhand = parsed_args.filtered
    args['filtered_fhand'] = filtered_fhand
    return args, parsed_args
