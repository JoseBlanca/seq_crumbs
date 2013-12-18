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
from crumbs.utils.file_formats import get_format, set_format
from bisect import bisect

try:
    import pysam
except ImportError:
    # This is an optional requirement
    pass

from crumbs.utils.tags import (SEQS_PASSED, SEQS_FILTERED_OUT, SEQITEM,
                               SEQRECORD, CHIMERA, NON_CHIMERIC, UNKNOWN)
from crumbs.utils.seq_utils import uppercase_length, get_uppercase_segments
from crumbs.seq import get_name, get_file_format, get_str_seq, get_length
from crumbs.exceptions import WrongFormatError
from crumbs.blast import Blaster
from crumbs.statistics import calculate_dust_score, IntCounter, draw_histogram
from crumbs.settings import get_setting
from crumbs.mapping import (get_or_create_bowtie2_index, map_with_bowtie2,
                            map_process_to_bam, sort_mapped_reads,
                            get_or_create_bwa_index, map_with_bwamem,
                            alignedread_to_seqitem)
from crumbs.seqio import write_seqs, read_seqs
from crumbs.pairs import group_pairs, group_pairs_by_name, deinterleave_pairs


def seq_to_filterpackets(seq_packets, group_paired_reads=False):
    'It yields packets suitable for the filters'

    for packet in seq_packets:
        if group_paired_reads:
            packet = list(group_pairs_by_name(packet))
        else:
            packet = list(group_pairs(packet, n_seqs_in_pair=1))
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
        map_process = map_with_bowtie2(index_fpath,
                                       unpaired_fpaths=[reads_fhand.name],
                                       extra_params=extra_params)
        map_process_to_bam(map_process, bam_fhand.name)

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


def _get_cigar_seq(cigar):
    cigar_seq = ''
    for element in cigar:
        cigar_seq += str(element[0]) * element[1]
    return cigar_seq


def _get_mapping_positions(alignment):
    cigar_seq = _get_cigar_seq(alignment.cigar)
    cigar_seq = list(cigar_seq)
    if alignment.is_reverse:
        cigar_seq.reverse()
    positions = set()
    for i in range(len(cigar_seq)):
        if cigar_seq[i] in ['0', '1', '2']:
            positions.add(i)
    return positions


def _group_alignments_by_reads(sorted_alignments):
    prev_name = None
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
    #We add a tag to differenciate between mates afterwards. This tag should
    #be probably changed .f, .r
    forward = []
    reverse = []
    for alignment_read in alignments_group:
        if alignment_read.is_read1:
            alignment_read.qname += '.f'
            forward.append(alignment_read)
        elif alignment_read.is_read2:
            alignment_read.qname += '.r'
            reverse.append(alignment_read)
        elif alignment_read.is_unmapped:
            print 'aaaa'
            print alignment_read
    return [forward, reverse]


def _pop_primary_alignment(alignments_group):
    primary_alignment = _get_primary_alignment(alignments_group)
    alignments_group.remove(primary_alignment)
    return primary_alignment


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


def _alignment_at_ends_of_reference(alignment_read, limit, samfile):
    lengths = _get_ref_lengths(samfile)
    start = alignment_read.pos
    end = alignment_read.aend
    name = samfile.getrname(alignment_read.rname)
    length = lengths[name]
    if alignment_read.is_reverse:
        start = alignment_read.aend
        end = alignment_read.pos
    if start > limit and end < length - limit:
        return False
    else:
        return True


def _read_is_totally_mapped(alignments_group, max_clipping):
    for alignment_read in alignments_group:
        if alignment_read.is_unmapped:
            return False
        cigar = alignment_read.cigar
        max_clipping_positions = max_clipping * alignment_read.alen
        if _count_cigar_char(cigar, [4, 5]) <= max_clipping_positions:
            return True
    return False


def _find_secondary_fragment(alignments_group, max_coincidences,
                        max_mapq_difference):
    if len(alignments_group) == 1:
        return None
    primary_alignment = _pop_primary_alignment(alignments_group)
    primary_mapping_positions = _get_mapping_positions(primary_alignment)
    for alignment_read in alignments_group:
        mapping_positions = _get_mapping_positions(alignment_read)
        coincidences = mapping_positions.intersection(primary_mapping_positions)
        n_coincidences = len(coincidences)
        difference = primary_alignment.mapq - alignment_read.mapq
        if difference < max_mapq_difference:
            if n_coincidences <= max_coincidences:
                alignments_group.append(primary_alignment)
                return alignment_read
    alignments_group.append(primary_alignment)
    return None


def _find_distances_with_mate(alignment_reads, mate_alignments):
    '''It calculates the distances between two alignments of the same sequence
    with the position to which the mate maps. It assumes that it is a typical
    chimera'''
    alignment_reads.sort(key=lambda x: x.pos)
    primary_mate = _get_primary_alignment(mate_alignments)
    distance1 = abs(primary_mate.pos - alignment_reads[0].aend)
    distance2 = abs(alignment_reads[1].pos - primary_mate.aend)
    return sorted([distance1, distance2])


def _guess_kind_fragments_at_ends(fragments, samfile, limit):
    fragments.sort(key=lambda x: x.qstart)
    lengths = _get_ref_lengths(samfile)
    names = [samfile.getrname(fragment.rname) for fragment in fragments]
    lengths = [lengths[name] for name in names]
    if  ((fragments[0].is_reverse and fragments[0].aend > lengths[0] - limit)
        or (not fragments[0].is_reverse and fragments[0].pos < limit)):
        return CHIMERA
    else:
        if ((not fragments[1].is_reverse and fragments[1].pos < limit) or (
        fragments[1].is_reverse and fragments[1].aend > lengths[1] - limit)):
            return NON_CHIMERIC
        else:
            return CHIMERA


def _3end_mapped(alignment_read, max_clipping):
    #check differences between aend and qend, as well as alen, qlen, rlen, tlen
    #should we use a different max_clipping for this?
    length = alignment_read.qend - alignment_read.qstart
    max_clipping_positions = max_clipping * length
    if alignment_read.is_reverse:
        return alignment_read.qstart < max_clipping_positions
    else:
        return alignment_read.qend > alignment_read.aend - max_clipping_positions


def _get_mate(i, mates_alignments):
    if i == 1:
        return mates_alignments[1]
    else:
        return mates_alignments[0]


def _find_distance(aligned_reads):
    'It returns distance between two aligned_reads in the reference seq'
    aligned_reads.sort(key=lambda x: x.pos)
    return aligned_reads[1].pos - aligned_reads[0].aend


def _get_distances_to_mate(alignment_read, mate_alignments,
                           max_mapq_difference):
    '''It returns the distances and alignment between an aligned_read and all
    the alignments of the mate. It is < 0 depending on positions in referece'''
    primary_mate = _get_primary_alignment(mate_alignments)
    for alignment_mate in mate_alignments:
        if primary_mate.mapq - alignment_mate.mapq <= max_mapq_difference:
            if alignment_read.rname == alignment_mate.rname:
                distance = _find_distance([alignment_mate, alignment_read])
                yield [alignment_mate, abs(distance)]


def _mapped_read_is_chimeric(alignment_read, mate_alignments,
                          max_mapq_difference, max_pe_len):
    for alignment_mate, distance in _get_distances_to_mate(alignment_read,
                                                        mate_alignments,
                                                        max_mapq_difference):
        if distance < max_pe_len and _mates_are_innies([alignment_read,
                                                       alignment_mate]):
            return True
    return False


def _mates_are_outies(mates):
    mates.sort(key=lambda x: x.pos)
    return mates[0].is_reverse and (not mates[1].is_reverse)


def _mates_are_innies(mates):
    mates.sort(key=lambda x: x.pos)
    return (not mates[0].is_reverse) and mates[1].is_reverse


def _alignments_in_same_ref(alignments):
    prev_ref = None
    for alignment in alignments:
        if prev_ref is None:
            prev_ref = alignment.rname
        elif alignment.rname != prev_ref:
            return False
    return True


def _read_is_typical_chimera(fragments, max_pe_len, min_mp_len, mate):
    positions = [fragment.pos for fragment in fragments]
    if (fragments[0].is_reverse == fragments[1].is_reverse and
        bisect(positions, fragments[0].pnext) == 1 and
        fragments[0].is_reverse != fragments[0].mate_is_reverse):
        len_pe, len_mp = sorted(_find_distances_with_mate(fragments, mate))
        return len_pe < max_pe_len and len_mp > min_mp_len
    return False


def classify_mapped_reads(ref_fpath, paired_fpaths=None,
                     unpaired_fpaths=None, directory=None,
                     settings=get_setting('CHIMERAS_SETTINGS'),
                     file_format=None, paired_result=True,
                     min_seed_len=None):
    fhand = open(paired_fpaths[0]) if paired_fpaths else open(unpaired_fpaths[0])
    if file_format is not None:
        set_format(fhand, file_format)
    else:
        file_format = get_format(fhand)
    index_fpath = get_or_create_bwa_index(ref_fpath, directory)
    extra_params = ['-a', '-M']
    if min_seed_len is not None:
        extra_params.extend(['-k', min_seed_len])
    bwa = map_with_bwamem(index_fpath, paired_fpaths=paired_fpaths,
                         unpaired_fpath=unpaired_fpaths,
                         extra_params=extra_params)
    bam_fhand = NamedTemporaryFile(dir='/home/carlos/tmp')
    sort_mapped_reads(bwa, bam_fhand.name, key='queryname')
    bamfile = pysam.Samfile(bam_fhand.name)

    #settings. Include in function properties with default values
    max_coincidences = settings['MAX_COINCIDENCES']
    max_mapq_difference = settings['MAX_MAPQ_DIFFERENCE']
    limit = settings['MAX_DISTANCE_TO_END']
    max_clipping = settings['MAX_CLIPPING']
    max_pe_len = settings['MAX_PE_LEN']
    min_mp_len = settings['MIN_MP_LEN']

    #parameters for statistics
    len_mp_non_chimeric = 0
    n_non_chimeric_mates = 0
    len_pe_sum = 0
    len_mp_sum = 0
    typic_chimeras = 0
    pe_like_chimeras = 0
    #It tries to find out the kind of each pair of sequences
    for grouped_mates in _group_alignments_by_reads(bamfile):
        mates_alignments = _split_mates(grouped_mates)
        i = 0
        pair = []
        for alignments_group in mates_alignments:
            i += 1
            mate = _get_mate(i, mates_alignments)
            primary_mate = _get_primary_alignment(mate)
            primary_alignment = _get_primary_alignment(alignments_group)
            if _read_is_totally_mapped(alignments_group, max_clipping):
                if  primary_alignment.mate_is_unmapped:
                    kind = UNKNOWN
                else:
                    if primary_alignment.rname != primary_alignment.rnext:
                        kind = NON_CHIMERIC
                    else:
                        mates = [primary_mate, primary_alignment]
                        distance = _find_distance(mates)
                        if _mapped_read_is_chimeric(primary_alignment, mate,
                                                    max_mapq_difference,
                                                    max_pe_len):
                            kind = CHIMERA
                        elif _mates_are_outies(mates) and distance > min_mp_len:
                            len_mp_non_chimeric += abs(distance)
                            n_non_chimeric_mates += 1
                            kind = NON_CHIMERIC
                        else:
                            kind = UNKNOWN
            else:
                fragment = _find_secondary_fragment(alignments_group,
                                                    max_coincidences,
                                                    max_mapq_difference=100)
                #For fragmented reads. likely to be chimeric
                if fragment is not None:
                    fragments = [primary_alignment, fragment]
                    if _alignments_in_same_ref(fragments):
                        if _read_is_typical_chimera(fragments, max_pe_len,
                                                    min_mp_len, mate):
                            len_pe, len_mp = _find_distances_with_mate(fragments, mate)
                            len_pe_sum += len_pe
                            len_mp_sum += len_mp
                            typic_chimeras += 1
                            #typic chimera
                        kind = CHIMERA
                    else:
                        if (_alignment_at_ends_of_reference(fragments[0],
                                                            limit, bamfile)
                            and _alignment_at_ends_of_reference(fragments[1],
                                                                limit, bamfile)):
                            kind = _guess_kind_fragments_at_ends(fragments,
                                                                 bamfile, limit)
                        else:
                            kind = CHIMERA
                else:
                    #Find PE-like chimeras in partially mapping reads
                    if primary_alignment.is_unmapped:
                        kind = UNKNOWN
                    elif not _3end_mapped(primary_alignment, max_clipping):
                        kind = UNKNOWN
                    else:
                        if primary_alignment.rname == primary_alignment.rnext:
                            if _mapped_read_is_chimeric(primary_alignment,
                                                    mate, max_mapq_difference,
                                                    max_pe_len):
                                pe_like_chimeras += 1
                                kind = CHIMERA
                            else:
                                kind = UNKNOWN
                        else:
                            kind = UNKNOWN
            read = [alignedread_to_seqitem(alignments_group[0], file_format),
                    kind]
            #read = [alignments_group, kind]
            if paired_result == False:
                yield [read[0]], read[1]
            else:
                pair.append(read)
        if paired_result:
            kinds = [read[1] for read in pair]
            reads = [read[0] for read in pair]
            if CHIMERA in kinds:
                yield [reads, CHIMERA]
            elif UNKNOWN in kinds:
                yield [reads, UNKNOWN]
            else:
                yield [reads, NON_CHIMERIC]
    try:
        mean_mp_non_chimeric = len_mp_non_chimeric / float(n_non_chimeric_mates)
        mean_pe_len = len_pe_sum / float(typic_chimeras)
        mean_mp_len = len_mp_sum / float(typic_chimeras)
        print 'MP mean length in non chimeric reads: ', mean_mp_non_chimeric
        print 'Typical chimera number: ', typic_chimeras
        print 'PE-like chimera number: ', pe_like_chimeras
        print 'PE mean length in chimeric reads: ', mean_pe_len
        print 'MP mean length in chimeric reads: ', mean_mp_len
    except ZeroDivisionError:
        print 'typic chimeras or non_chimeric_mates not found'


def filter_chimeras(ref_fpath, out_fhand, chimeras_fhand, in_fhands,
                    unknown_fhand, unpaired=False, paired_result=True,
                    settings=get_setting('CHIMERAS_SETTINGS'),
                    min_seed_len=None):
    file_format = get_format(in_fhands[0])
    if unpaired:
        unpaired_fpaths = [fhand.name for fhand in in_fhands]
        paired_fpaths = None
    else:
        f_fhand = NamedTemporaryFile()
        r_fhand = NamedTemporaryFile()
        seqs = read_seqs(in_fhands)
        deinterleave_pairs(seqs, f_fhand, r_fhand, file_format)
        paired_fpaths = [f_fhand.name, r_fhand.name]
        unpaired_fpaths = None

    total = 0
    chimeric = 0
    unknown = 0
    for pair, kind in classify_mapped_reads(ref_fpath=ref_fpath,
                                           unpaired_fpaths=unpaired_fpaths,
                                           paired_fpaths=paired_fpaths,
                                           settings=settings,
                                           file_format=file_format,
                                           paired_result=paired_result,
                                           min_seed_len=min_seed_len):
        if kind is NON_CHIMERIC:
            write_seqs(pair, out_fhand)
        elif kind is CHIMERA and chimeras_fhand is not None:
            write_seqs(pair, chimeras_fhand)
            chimeric += 1
        elif kind is UNKNOWN and unknown_fhand is not None:
            write_seqs(pair, unknown_fhand)
            unknown += 1
        total += 1
    mapped = total - chimeric - unknown
    print 'Total pairs analyzed: ', total
    print 'Chimeric pairs filtered: ', chimeric, '\t', chimeric / float(total)
    print 'Unknown pairs found: ', unknown, '\t', unknown / float(total)
    print 'Non-chimeric pairs: ', mapped, '\t', mapped / float(total)


def classify_mapped_reads_new(bam_fpath, n=None):
    bamfile = pysam.Samfile(bam_fpath)
    type1 = 0
    type2a = []
    type2b = []
    type3a = []
    type3b = []
    type4 = []
    type5 = []
    type6 = 0
    others = 0
    h = 0

    #It tries to find out the kind of each pair of sequences
    for grouped_mates in _group_alignments_by_reads(bamfile):
        if n is not None and h == n:
            break
        h += 1
        mates_alignments = _split_mates(grouped_mates)
        i = 0
        pair = []
        for alignments_group in mates_alignments:
            i += 1
            mate = _get_mate(i, mates_alignments)
            primary_mate = _get_primary_alignment(mate)
            primary_alignment = _get_primary_alignment(alignments_group)
            mates = [primary_alignment, primary_mate]
            if _read_is_totally_mapped(alignments_group, 0.05):
                if  primary_alignment.mate_is_unmapped:
                    kind = '1'
                else:
                    if primary_alignment.rname != primary_alignment.rnext:
                        kind = '6'
                    else:
                        if _mates_are_outies(mates):
                            kind = '3a'
                        elif _mates_are_innies(mates):
                            kind = '2a'
                        else:
                            kind = '5'
            else:
                fragment = _find_secondary_fragment(alignments_group, 5, 100)
                if fragment is not None:
                    fragments = [primary_alignment, fragment]
                    if (_alignments_in_same_ref([fragments[0], primary_mate])
                        or _alignments_in_same_ref([fragments[1], primary_mate])):
                        kind = '4'
                    else:
                        kind = 'other'
                else:
                    if primary_alignment.is_unmapped:
                        kind = '1'
                    else:
                        if primary_alignment.rname == primary_alignment.rnext:
                            if _mates_are_outies(mates):
                                kind = '3b'
                            elif _mates_are_innies(mates):
                                kind = '2b'
                            else:
                                kind = '5'
                        else:
                            kind = '6'
            pair.append(kind)
        if '1' in pair:
            type1 += 1
        elif '6' in pair:
            type6 += 1
        elif 'other' in pair:
            others += 1
        else:
            distance = _find_distance(mates)
            if '4' in pair:
                type4.append(distance)
            elif '2b' in pair:
                type2b.append(distance)
            elif '3b' in pair:
                type3b.append(distance)
            elif '2a' in pair:
                type2a.append(distance)
            elif '3a' in pair:
                type3a.append(distance)
            elif '5' in pair:
                type5.append(distance)
    stats1 = {'1': type1, '6': type6, 'other': others}
    stats2 = {'4': type4, '2b': type2b, '3b': type3b, '2a': type2a,
              '3a': type3a, '5': type5}
    for key in stats1.keys():
        print key.ljust(5), stats1[key]
    for key in stats2.keys():
        print key, len(stats2[key])
        if '2' in key:
            counter = IntCounter(iter(stats2[key]))
            distribution = counter.calculate_distribution(remove_outliers=True)
            counts = distribution['counts']
            bin_limits = distribution['bin_limits']
            print draw_histogram(bin_limits, counts)



