
import argparse
import sys

# Missing docstring
# pylint: disable=C0111


def setup_basic_argparse(**kwargs):
    'It prepares the command line argument parsing.'

    parser = argparse.ArgumentParser(**kwargs)
    in_help = 'Input VCF file (default STDIN)'
    parser.add_argument('input', help=in_help, nargs='?',
                        type=argparse.FileType('r'), default=sys.stdin)
    parser.add_argument('-o', '--output', default=sys.stdout,
                        help='Output VCF file (default STDOUT)',
                        type=argparse.FileType('w'))
    msg = 'Template VCF to get the header (default same as input)'
    parser.add_argument('-t', '--template', help=msg,
                        type=argparse.FileType('r'))
    msg = 'File to print some statistics (default STDERR)'
    parser.add_argument('-l', '--log', help=msg, type=argparse.FileType('w'),
                        default=sys.stderr)

    return parser


def setup_filter_argparse(**kwargs):
    'It prepares the command line argument parsing.'
    parser = setup_basic_argparse(**kwargs)
    parser.add_argument('-f', '--filtered',
                        help='Output for filtered SNVs',
                        type=argparse.FileType('w'))
    parser.add_argument('-s', '--samples', action='append',
                        help='samples to use')
    parser.add_argument('-p', '--samples_file',
                        help='File with samples to use. One per line',
                        type=argparse.FileType('r'))
    return parser


def parse_basic_args(parser):
    parsed_args = parser.parse_args()
    in_fhand = parsed_args.input
    template_fhand = parsed_args.template

    if template_fhand is None:
        in_fname = in_fhand.name
        if in_fname == '<stdin>':
            msg = 'A template file has to be provided when the input is given '
            msg += 'via STDIN'
            parser.error(msg)
        template_fhand = open(in_fname)

    out_fhand = parsed_args.output
    log_fhand = parsed_args.log

    args = {'in_fhand': in_fhand, 'template_fhand': template_fhand,
            'log_fhand': log_fhand, 'out_fhand': out_fhand}

    return args, parsed_args


def _parse_sample_file(fhand):
    samples = []
    for line in fhand:
        line = line.strip()
        if not line:
            continue
        samples.append(line)
    return samples


def parse_filter_args(parser):
    'It parses the command line and it returns a dict with the arguments.'
    filter_snvs_args, parsed_args = parse_basic_args(parser)

    filtered_fhand = parsed_args.filtered
    filter_snvs_args['filtered_fhand'] = filtered_fhand

    samples = set()
    if parsed_args.samples is not None:
        samples.update(parsed_args.samples)
    if parsed_args.samples_file is not None:
        samples.update(_parse_sample_file(parsed_args.samples_file))

    filter_class_kwargs = {'samples_to_consider': samples if samples else None}

    return filter_snvs_args, filter_class_kwargs, parsed_args 
