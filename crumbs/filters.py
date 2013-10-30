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

# pylint: disable=C0111

from __future__ import division
from tempfile import NamedTemporaryFile

try:
    import pysam
except ImportError:
    # This is an optional requirement
    pass

from crumbs.utils.tags import (SEQS_PASSED, SEQS_FILTERED_OUT, SEQITEM,
                               SEQRECORD)
from crumbs.utils.seq_utils import uppercase_length, get_uppercase_segments
from crumbs.seq import get_name, get_file_format, get_str_seq, get_length
from crumbs.exceptions import WrongFormatError
from crumbs.blast import Blaster
from crumbs.statistics import calculate_dust_score
from crumbs.settings import get_setting
from crumbs.mapping import  get_or_create_bowtie2_index, map_with_bowtie2
from crumbs.seqio import write_seqs
from crumbs.pairs import group_seqs_in_pairs


def seq_to_filterpackets(seq_packets, group_paired_reads=False):
    'It yields packets suitable for the filters'

    for packet in seq_packets:
        if group_paired_reads:
            packet = list(group_seqs_in_pairs(packet))
        else:
            packet = [[seq] for seq in packet]
        yield {SEQS_PASSED: packet, SEQS_FILTERED_OUT: []}


def _reverse_bools(bools):
    return (not bool_ for bool_ in bools)


class _BaseFilter(object):
    def __init__(self, failed_drags_pair=True, reverse=False):
        self.reverse = reverse
        self.failed_drags_pair = failed_drags_pair

    def _setup_checks(self, filterpacket):
        pass

    def _do_check(self, seq):
        raise NotImplementedError()

    def __call__(self, filterpacket):
        self._setup_checks(filterpacket)
        reverse = self.reverse
        failed_drags_pair = self.failed_drags_pair
        seqs_passed = []
        filtered_out = filterpacket[SEQS_FILTERED_OUT][:]
        for paired_seqs in filterpacket[SEQS_PASSED]:
            checks = (self._do_check(seq) for seq in paired_seqs)
            if reverse:
                checks = _reverse_bools(checks)

            if failed_drags_pair:
                pair_passed = True if all(checks) else False
            else:
                pair_passed = True if any(checks) else False

            if pair_passed:
                seqs_passed.append(paired_seqs)
            else:
                filtered_out.append(paired_seqs)

        return {SEQS_PASSED: seqs_passed, SEQS_FILTERED_OUT: filtered_out}


class FilterByFeatureTypes(_BaseFilter):
    'It filters out sequences not annotated with the given feature types'
    def __init__(self, feature_types, reverse=False, failed_drags_pair=True):
        '''The initiator

        feat_types is a list of types of features to use to filter'''

        self._feat_types = feature_types
        super(FilterByFeatureTypes, self).__init__(reverse=reverse,
                                           failed_drags_pair=failed_drags_pair)

    def _do_check(self, seq):
        seq = seq.object
        f_in_seq = [f.type for f in seq.features if f.type in self._feat_types]
        return True if f_in_seq else False


class FilterByRpkm(_BaseFilter):
    def __init__(self, read_counts, min_rpkm, reverse=False,
                 failed_drags_pair=True):
        '''The init

        read_counts its a dict:
            - keys are the sequence names
            - values should be dicts with length, mapped_reads and
            unmmapped_reads, counts
        '''
        self._read_counts = read_counts
        self._min_rpkm = min_rpkm
        total_reads = sum([v['mapped_reads'] + v['unmapped_reads'] for v in read_counts.values()])
        self._million_reads = total_reads / 1e6
        super(FilterByRpkm, self).__init__(reverse=reverse,
                                           failed_drags_pair=failed_drags_pair)

    def _do_check(self, seq):
        count = self._read_counts[get_name(seq)]

        kb_len = count['length'] / 1000
        rpk = count['mapped_reads'] / kb_len  # rpks
        rpkm = rpk / self._million_reads  # rpkms

        return True if rpkm >= self._min_rpkm else False


class FilterByLength(_BaseFilter):
    'It removes the sequences according to their length.'
    def __init__(self, minimum=None, maximum=None, ignore_masked=False,
                 failed_drags_pair=True):
        '''The initiator.

        threshold - minimum length to pass the filter (integer)
        reverse - if True keep the short sequences and discard the long ones
        ignore_masked - If True only uppercase letters will be counted.
        '''
        if min is None and max is None:
            raise ValueError('min or max threshold must be given')
        self.min = minimum
        self.max = maximum
        self.ignore_masked = ignore_masked
        super(FilterByLength, self).__init__(reverse=False,
                                           failed_drags_pair=failed_drags_pair)

    def _do_check(self, seq):
        min_ = self.min
        max_ = self.max
        length = uppercase_length(get_str_seq(seq)) if self.ignore_masked else get_length(seq)

        passed = True
        if min_ is not None and length < min_:
            passed = False
        if max_ is not None and length > max_:
            passed = False
        return passed


class FilterById(_BaseFilter):
    'It removes the sequences not found in the given set'
    def __init__(self, seq_ids, failed_drags_pair=True, reverse=False):
        '''The initiator.

        seq_ids - An iterator with the sequence ids to keep
        reverse - if True keep the sequences not found on the list
        '''
        if not isinstance(seq_ids, set):
            seq_ids = set(seq_ids)
        self.seq_ids = seq_ids
        super(FilterById, self).__init__(failed_drags_pair=failed_drags_pair,
                                              reverse=reverse)

    def _do_check(self, seq):
        return True if get_name(seq) in self.seq_ids else False


def _get_mapped_reads(bam_fpath, min_mapq=0):
    bam = pysam.Samfile(bam_fpath)
    return [read.qname for read in bam if not read.is_unmapped and (not min_mapq or read.mapq > min_mapq)]


class FilterByBam(FilterById):
    'It filters the reads not mapped in the given BAM files'
    def __init__(self, bam_fpaths, min_mapq=0, reverse=False):
        seq_ids = self._get_mapped_reads(bam_fpaths, min_mapq)
        super(FilterByBam, self).__init__(seq_ids, reverse=reverse)

    def _get_mapped_reads(self, bam_fpaths, min_mapq):
        mapped_reads = []
        for fpath in bam_fpaths:
            mapped_reads.extend(_get_mapped_reads(fpath, min_mapq=min_mapq))
        return mapped_reads


class FilterByQuality(_BaseFilter):
    'It removes the sequences according to its quality'
    def __init__(self, threshold, ignore_masked=False, failed_drags_pair=True,
                 reverse=False):
        '''The initiator.

        threshold - minimum quality to pass the filter (float)
        reverse - if True keep the sequences not found on the list
        '''
        self.threshold = float(threshold)
        self.ignore_masked = ignore_masked
        super(FilterByQuality, self).__init__(reverse=reverse,
                                           failed_drags_pair=failed_drags_pair)

    def _do_check(self, seq):
        seq_object = seq.object
        try:
            quals = seq_object.letter_annotations['phred_quality']
        except KeyError:
            msg = 'Some of the input sequences do not have qualities: {}'
            msg = msg.format(get_name(seq))
            raise WrongFormatError(msg)
        if self.ignore_masked:
            str_seq = str(seq_object.seq)
            seg_quals = [quals[segment[0]: segment[1] + 1]
                            for segment in get_uppercase_segments(str_seq)]
            qual = sum(sum(q) * len(q) for q in seg_quals) / len(quals)
        else:
            qual = sum(quals) / len(quals)
        return True if qual >= self.threshold else False


class FilterBlastMatch(_BaseFilter):
    'It filters a seq if there is a match against a blastdb'
    def __init__(self, database, program, filters, dbtype=None,
                 failed_drags_pair=True, reverse=False):
        '''The initiator
            database: path to a file with seqs or a blast database
            filter_params:
                expect_threshold
                similarty treshlod
                min_length_percentaje
        '''
        self._blast_db = database
        self._blast_program = program
        self._filters = filters
        self._dbtype = dbtype
        super(FilterBlastMatch, self).__init__(reverse=reverse,
                                          failed_drags_pair=failed_drags_pair)

    def _setup_checks(self, filterpacket):
        seqs = [s for seqs in filterpacket[SEQS_PASSED]for s in seqs]
        self._matcher = Blaster(seqs, self._blast_db, dbtype=self._dbtype,
                                program=self._blast_program,
                                filters=self._filters)

    def _do_check(self, seq):
        segments = self._matcher.get_matched_segments(get_name(seq))
        return True if segments is None else False


class FilterBowtie2Match(_BaseFilter):
    'It filters a seq if it maps against a bowtie2 index'
    def __init__(self, index_fpath, reverse=False, min_mapq=None,
                 failed_drags_pair=True):
        self._index_fpath = index_fpath
        self._reverse = reverse
        self.min_mapq = min_mapq
        super(FilterBowtie2Match, self).__init__(reverse=reverse,
                                          failed_drags_pair=failed_drags_pair)

    def _setup_checks(self, filterpacket):
        index_fpath = self._index_fpath
        get_or_create_bowtie2_index(index_fpath)
        seqs = [s for seqs in filterpacket[SEQS_PASSED]for s in seqs]
        seq_class = seqs[0].kind
        extra_params = []
        # Which format do we need for the bowtie2 input read file fasta or
        # fastq?
        if seq_class == SEQRECORD:
            if 'phred_quality' in seqs[0].object.letter_annotations.viewkeys():
                file_format = 'fastq'
            else:
                extra_params.append('-f')
                file_format = 'fasta'
        elif seq_class == SEQITEM:
            file_format = get_file_format(seqs[0])
            if 'illumina' in file_format:
                extra_params.append('--phred64')
            elif 'fasta' in file_format:
                extra_params.append('-f')
            elif 'fastq' in file_format:
                pass
            else:
                msg = 'For FilterBowtie2Match and SeqItems fastq or fasta '
                msg += 'files are required'
                raise RuntimeError(msg)
        else:
            raise NotImplementedError()

        reads_fhand = NamedTemporaryFile(suffix=file_format)
        write_seqs(seqs, reads_fhand, file_format=file_format)
        reads_fhand.flush()

        bam_fhand = NamedTemporaryFile(suffix='.bam')
        map_with_bowtie2(index_fpath, bam_fhand.name,
                         unpaired_fpaths=[reads_fhand.name],
                         extra_params=extra_params)

        self.mapped_reads = _get_mapped_reads(bam_fhand.name, self.min_mapq)

    def _do_check(self, seq):
        return False if get_name(seq) in self.mapped_reads else True


class FilterDustComplexity(_BaseFilter):
    'It filters a sequence according to its dust score'
    def __init__(self, threshold=get_setting('DEFATULT_DUST_THRESHOLD'),
                 reverse=False, failed_drags_pair=True):
        '''The initiator
        '''
        self._threshold = threshold
        super(FilterDustComplexity, self).__init__(reverse=reverse,
                                          failed_drags_pair=failed_drags_pair)

    def  _do_check(self, seq):
        threshold = self._threshold
        dustscore = calculate_dust_score(seq)
        return True if dustscore < threshold else False
