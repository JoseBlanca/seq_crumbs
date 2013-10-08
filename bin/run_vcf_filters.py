#!/usr/bin/env python
import sys
import argparse
import json
import os.path

from configobj import ConfigObj

from vcf import Reader, Writer
from vcf.parser import _Filter, _Info

FILTERS_USING_VCF = ('CloseToSnvFilter', 'HighVariableRegionFilter')


def _setup_argparse():
    'It prepares the command line argument parsing.'
    description = 'Filter the snvs in an indexed vcf'
    parser = argparse.ArgumentParser(description=description)
    in_help = 'Indexed vcf file'
    parser.add_argument('input', help=in_help, type=argparse.FileType('rt'))
    parser.add_argument('-o', '--output',
                        help='output vcf file', type=argparse.FileType('wt'))
    parser.add_argument('-f', '--filter_conf', required=True,
                        type=argparse.FileType('rt'),
                        help='File with the filter configuration')

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
    args['filter_fhand'] = parsed_args.filter_conf
    return args


def _parse_filters_conf(fpath, snv_in_fpath):
    config = ConfigObj(fpath, unrepr=True)
    filters = []
    for filter_confs in config.values():
        filter_class_name = filter_confs.keys()[0]
        args = filter_confs.values()[0]

        mod = __import__('vcf_crumbs.vcf_filters',
                         fromlist=[filter_class_name])
        filter_class = getattr(mod, filter_class_name)
        if filter_class_name in FILTERS_USING_VCF:
            args['vcf_fpath'] = os.path.abspath(snv_in_fpath)
        filter_ = filter_class(**args)
        filters.append(filter_)
    return filters


def _manage_mappers(reader, filters_in_conf):
    mappers_to_use = []
    for mapper in filters_in_conf:
        if mapper.name in reader.filters:
            new_conf = mapper.conf
            conf_in_vcf = reader.filters[mapper.name].desc.split('::')[1]
            conf_in_vcf = conf_in_vcf.replace("'", "\"")
            conf_in_vcf = json.loads(conf_in_vcf)
            if new_conf == conf_in_vcf:  # no need to run again
                continue
            else:
                mappers_to_use.append(mapper)
        else:
            mappers_to_use.append(mapper)

    # se anyaden la metadata de FILTERS e INFO al reader
    for mapper in mappers_to_use:
        reader.filters[mapper.name] = _Filter(mapper.name,
                                               mapper.description)
        if mapper.info:
            info = mapper.info
            reader.infos[info['id']] = _Info(**info)
    return mappers_to_use


def main():
    parser = _setup_argparse()
    args = _parse_args(parser)

    reader = Reader(args['in_fhand'])
    filters_in_conf = _parse_filters_conf(args['filter_fhand'],
                                            args['in_fhand'].name)
    mappers = _manage_mappers(reader, filters_in_conf)
    writer = Writer(args['out_fhand'], reader)

    for record in reader:
        for mapper in mappers:
            record = mapper(record)
        writer.write_record(record)
    writer.close()
if __name__ == '__main__':
    main()
