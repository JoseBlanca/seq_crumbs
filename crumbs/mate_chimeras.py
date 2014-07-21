'''
Created on 2014 uzt 10

@author: peio
'''
from __future__ import division
from tempfile import NamedTemporaryFile

try:
    from pysam import Samfile
except ImportError:
    # This is an optional requirement
    pass

from crumbs.utils.tags import (SEQS_PASSED, SEQS_FILTERED_OUT, CHIMERA,
                               NON_CHIMERIC, UNKNOWN)
from crumbs.statistics import IntCounter
from crumbs.settings import get_setting
from crumbs.mapping import (map_process_to_sortedbam, map_with_bwamem,
                            alignedread_to_seqitem)
from crumbs.seqio import write_seqs
from crumbs.pairs import group_pairs, group_pairs_by_name


def seq_to_filterpackets(seq_packets, group_paired_reads=False):
    'It yields packets suitable for the filters'

    for packet in seq_packets:
        if group_paired_reads:
            packet = list(group_pairs_by_name(packet))
        else:
            packet = list(group_pairs(packet, n_seqs_in_pair=1))
        yield {SEQS_PASSED: packet, SEQS_FILTERED_OUT: []}


def _group_alignments_reads_by_qname(sorted_alignments):
    prev_name = None
    group = None
    for alignment in sorted_alignments:
        name = alignment.qname
        if prev_name is None:
            group = [alignment]
        elif prev_name == name:
            group.append(alignment)
        else:
            yield group
            group = [alignment]
        prev_name = name
    yield group


def _split_mates(alignments_group):
    # We add a tag to differenciate between mates afterwards. This tag should
    # be probably taken from the file
    forwards = []
    reverses = []
    for alignment_read in alignments_group:
        if alignment_read.is_read1:
            alignment_read.qname += ' 1:N:0:GATCAG'
            forwards.append(alignment_read)
        elif alignment_read.is_read2:
            alignment_read.qname += ' 2:N:0:GATCAG'
            reverses.append(alignment_read)
    return [forwards, reverses]


def _get_primary_alignment(alignments_group):
    for alignment in alignments_group:
        if not alignment.is_secondary:
            return alignment


def _get_ref_lengths(samfile):
    lengths = {}
    for ref in samfile.header['SQ']:
        lengths[ref['SN']] = ref['LN']
    return lengths


def _count_cigar_char(cigar, numbers):
    counts = 0
    for element in cigar:
        if element[0] in numbers:
            counts += element[1]
    return counts


def _find_distances_to_ends(alignment_read, samfile, reference_lengths):

    start = alignment_read.pos
    end = alignment_read.aend
    name = samfile.getrname(alignment_read.rname)
    length = reference_lengths[name]
    if alignment_read.is_reverse:
        start = alignment_read.aend
        end = alignment_read.pos
    return [start, length - end]


def _read_is_totally_mapped(alignments_group, max_clipping):
    for alignment_read in alignments_group:
        if alignment_read.is_unmapped:
            return False
        cigar = alignment_read.cigar
        max_clipping_positions = max_clipping * alignment_read.alen
        if _count_cigar_char(cigar, [4, 5]) <= max_clipping_positions:
            return True
    return False


# when having hard clipping qend, qstart and qlen are always 0, so we need
# different ways to calculate them using cigar information
def _get_qstart(aligned_read):
    positions = 0
    for element in aligned_read.cigar:
        if element[0] == 0:
            return positions
        positions += element[1]


def _get_qend(aligned_read):
    qstart = _get_qstart(aligned_read)
    end = qstart + aligned_read.alen
    for element in aligned_read.cigar:
        if element[0] == 1:
            end += element[1]
    return end


def _get_qlen(aligned_read):
    positions = 0
    for element in aligned_read.cigar:
        positions += element[1]
    return positions


def _3end_mapped(alignment_read, max_clipping):
    # check differences between aend and qend, as well as alen, qlen, rlen,
    # tlen should we use a different max_clipping for this?
    if alignment_read.is_unmapped:
        return False
    max_clipping_positions = max_clipping * alignment_read.alen
    alignment_qstart = _get_qstart(alignment_read)
    alignment_qend = _get_qend(alignment_read)
    alignment_qlen = _get_qlen(alignment_read)
    if alignment_read.is_reverse:
        return alignment_qstart < max_clipping_positions
    else:
        return alignment_qend > alignment_qlen - max_clipping_positions


def _find_distance(aligned_reads):
    'It returns distance between two aligned_reads in the reference seq'
    aligned_reads.sort(key=lambda x: x.pos)
    return aligned_reads[1].pos - aligned_reads[0].aend


def _mates_are_outies(mates):
    mates.sort(key=lambda x: x.pos)
    return mates[0].is_reverse and (not mates[1].is_reverse)


def _mates_are_innies(mates):
    mates.sort(key=lambda x: x.pos)
    return (not mates[0].is_reverse) and mates[1].is_reverse


def _get_totally_mapped_alignments(reads, max_clipping):
    for aligned_read in reads:
        if _read_is_totally_mapped([aligned_read], max_clipping):
            yield aligned_read


def _5end_mapped(aligned_read, max_clipping):
    if aligned_read.is_unmapped:
        return False
    aligned_read_qend = _get_qend(aligned_read)
    aligned_read_qstart = _get_qstart(aligned_read)
    max_clipping_positions = max_clipping * aligned_read.alen
    if aligned_read.is_reverse:
        return aligned_read_qend > (_get_qlen(aligned_read) -
                                    max_clipping_positions)
    else:
        return aligned_read_qstart < max_clipping_positions


def _mates_are_not_chimeric(mates, max_clipping, mate_length_range, bamfile,
                            reference_lengths):
    for aligned_read1 in _get_totally_mapped_alignments(mates[0],
                                                        max_clipping):
        distances_to_end1 = _find_distances_to_ends(aligned_read1, bamfile,
                                                    reference_lengths)
        for aligned_read2 in _get_totally_mapped_alignments(mates[1],
                                                            max_clipping):
            if aligned_read1.rname == aligned_read2.rname:
                aligned_reads = [aligned_read1, aligned_read2]
                distance = _find_distance(aligned_reads)
                # print 'distance: ', distance
                if (_mates_are_outies(aligned_reads) and
                    distance > mate_length_range[0] and
                    distance < mate_length_range[1]):
                    return True
            else:
                distances_to_end2 = _find_distances_to_ends(aligned_read2,
                                                            bamfile,
                                                            reference_lengths)
                if aligned_read1.is_reverse:
                    distances_sum = distances_to_end1[1]
                    if aligned_read2.is_reverse:
                        distances_sum += distances_to_end2[1]
                    else:
                        distances_sum += distances_to_end2[0]
                else:
                    distances_sum = distances_to_end1[0]
                    if aligned_read2.is_reverse:
                        distances_sum += distances_to_end2[1]
                    else:
                        distances_sum += distances_to_end2[0]
                # print 'distances_sum: ', distances_sum
                if distances_sum < mate_length_range[1]:
                    return True
        return False


def _mates_are_chimeric(mates, samfile, max_clipping, max_insert_size,
                        reference_lengths):
    for aligned_read1 in mates[0]:
        if _3end_mapped(aligned_read1, max_clipping):
            distances_to_end1 = _find_distances_to_ends(aligned_read1, samfile,
                                                        reference_lengths)
            for aligned_read2 in mates[1]:
                if _3end_mapped(aligned_read2, max_clipping):
                    if aligned_read1.rname == aligned_read2.rname:
                        aligned_reads = [aligned_read1, aligned_read2]
                        distance = _find_distance(aligned_reads)
                        if (_mates_are_innies(aligned_reads) and
                            distance < max_insert_size):
                            return True
                    else:
                        distances_to_end2 = _find_distances_to_ends(aligned_read2,
                                                                    samfile,
                                                             reference_lengths)
                        if aligned_read1.is_reverse:
                            distances_sum = distances_to_end1[0]
                            if aligned_read2.is_reverse:
                                distances_sum += distances_to_end2[0]
                            else:
                                distances_sum += distances_to_end2[1]
                        else:
                            distances_sum = distances_to_end1[1]
                            if aligned_read2.is_reverse:
                                distances_sum += distances_to_end2[0]
                            else:
                                distances_sum += distances_to_end2[1]
                        if distances_sum < max_insert_size:
                            return True
    return False


def classify_mapped_reads(bam_fhand, mate_distance,
                          settings=get_setting('CHIMERAS_SETTINGS')):
    '''It classifies sequences from bam file in chimeric, unknown and
    non chimeric, according to its distance and orientation in the reference
    sequence'''
    bamfile = Samfile(bam_fhand.name)

    # settings. Include in function properties with default values
    max_clipping = settings['MAX_CLIPPING']
    max_pe_len = settings['MAX_PE_LEN']
    variation = settings['MATE_DISTANCE_VARIATION']
    mate_length_range = [mate_distance - variation, mate_distance + variation]
    reference_lengths = _get_ref_lengths(bamfile)
    # It tries to find out the kind of each pair of sequences
    for grouped_mates in _group_alignments_reads_by_qname(bamfile):
        mates_alignments = _split_mates(grouped_mates)
        if _mates_are_not_chimeric(mates_alignments, max_clipping,
                                   mate_length_range, bamfile,
                                   reference_lengths):
            kind = NON_CHIMERIC
        elif _mates_are_chimeric(mates_alignments, bamfile, max_clipping,
                                 max_pe_len, reference_lengths):
            kind = CHIMERA
        else:
            kind = UNKNOWN

        pair = [alignedread_to_seqitem(_get_primary_alignment(mates))
                for mates in mates_alignments]

        if None not in pair:
            yield pair, kind


def classify_chimeras(in_fhand, index_fpath, mate_distance, out_fhand,
                      chimeras_fhand=None, unknown_fhand=None, tempdir=None,
                      threads=None, settings=get_setting('CHIMERAS_SETTINGS')):

    '''It maps sequences from input files, sorts them and writes to output
    files according to its classification'''
    bam_fhand = NamedTemporaryFile(suffix='.bam')
    extra_params = ['-a', '-M']
    bwa = map_with_bwamem(index_fpath, interleave_fpath=in_fhand.name,
                          extra_params=extra_params)
    map_process_to_sortedbam(bwa, bam_fhand.name, key='queryname',
                             tempdir=tempdir)

    for pair, kind in classify_mapped_reads(bam_fhand, settings=settings,
                                            mate_distance=mate_distance):
        if kind is NON_CHIMERIC:
            write_seqs(pair, out_fhand)
        elif kind is CHIMERA and chimeras_fhand is not None:
            write_seqs(pair, chimeras_fhand)
        elif kind is UNKNOWN and unknown_fhand is not None:
            write_seqs(pair, unknown_fhand)


def calculate_distance_distribution(interleave_fhand, index_fpath,
                                    max_clipping, max_distance=None,
                                    tempdir=None, threads=None):
    bam_fhand = NamedTemporaryFile(suffix='.bam')
    extra_params = ['-a', '-M']
    bwa = map_with_bwamem(index_fpath, interleave_fpath=interleave_fhand.name,
                          extra_params=extra_params)
    map_process_to_sortedbam(bwa, bam_fhand.name, key='queryname',
                             temdir=tempdir)
    bamfile = Samfile(bam_fhand.name)
    stats = {'outies': IntCounter(), 'innies': IntCounter(),
             'others': IntCounter()}
    for grouped_mates in _group_alignments_reads_by_qname(bamfile):
        mates = _split_mates(grouped_mates)
        for aligned_read1 in _get_totally_mapped_alignments(mates[0],
                                                            max_clipping):
            for aligned_read2 in _get_totally_mapped_alignments(mates[1],
                                                                max_clipping):
                if aligned_read1.rname == aligned_read2.rname:
                    aligned_reads = [aligned_read1, aligned_read2]
                    distance = _find_distance(aligned_reads)
                    if _mates_are_outies(aligned_reads):
                        kind = 'outies'
                    elif _mates_are_innies(aligned_reads):
                        kind = 'innies'
                    else:
                        kind = 'others'
                    if max_distance is None or max_distance > distance:
                        stats[kind][distance] += 1
    return stats
