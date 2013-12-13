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

from operator import itemgetter

from crumbs.utils.optional_modules import Seq
from crumbs.utils.tags import (TRIMMING_RECOMMENDATIONS, QUALITY, OTHER,
                               VECTOR, TRIMMING_KINDS, SEQS_PASSED,
                               ORPHAN_SEQS)
from crumbs.utils.seq_utils import get_uppercase_segments
from crumbs.seq import (copy_seq, get_str_seq, get_annotations, get_length,
                        slice_seq, get_int_qualities, get_name)
from crumbs.utils.segments_utils import (get_longest_segment, get_all_segments,
                                         get_longest_complementary_segment,
                                         merge_overlaping_segments)
from crumbs.utils.tags import SEQRECORD
from crumbs.iterutils import rolling_window
from crumbs.blast import BlasterForFewSubjects
from crumbs.seqio import write_seqs
from crumbs.pairs import group_pairs_by_name, group_pairs

# pylint: disable=R0903


def seq_to_trim_packets(seq_packets, group_paired_reads=False):
    'It yields packets suitable for the filters'

    for packet in seq_packets:
        if group_paired_reads:
            packet = list(group_pairs_by_name(packet))
        else:
            packet = list(group_pairs(packet, n_seqs_in_pair=1))
        yield {SEQS_PASSED: packet, ORPHAN_SEQS: []}


class _BaseTrim(object):
    'Base Trim class'
    def __call__(self, trim_packet):
        'It trims the seqs'
        self._pre_trim(trim_packet)
        trimmed_seqs = []
        for paired_seqs in trim_packet[SEQS_PASSED]:
            trimmed_seqs.append([self._do_trim(s) for s in paired_seqs])
        return {SEQS_PASSED: trimmed_seqs,
                ORPHAN_SEQS: trim_packet[ORPHAN_SEQS]}

    def _do_trim(self, seq):
        raise NotImplementedError()

    def _pre_trim(self, trim_packet):
        pass


class TrimLowercasedLetters(_BaseTrim):
    'It trims the masked segments of the seqrecords.'

    def _do_trim(self, seq):
        str_seq = get_str_seq(seq)
        unmasked_segments = get_uppercase_segments(str_seq)
        segment = get_longest_segment(unmasked_segments)
        if segment is not None:
            segments = []
            if segment[0] != 0:
                segments.append((0, segment[0] - 1))
            len_seq = len(str_seq)
            if segment[1] != len_seq - 1:
                segments.append((segment[1] + 1, len_seq - 1))

            _add_trim_segments(segments, seq, kind=OTHER)

        else:
            segments = [(0, len(seq))]
            _add_trim_segments(segments, seq, kind=OTHER)
        return seq


def _add_trim_segments(segments, sequence, kind):
    'It adds segments to the trimming recommendation in the annotation'
    assert kind in TRIMMING_KINDS
    if not segments:
        return
    annotations = sequence.object.annotations
    if TRIMMING_RECOMMENDATIONS not in annotations:
        annotations[TRIMMING_RECOMMENDATIONS] = {}
        for trim_kind in TRIMMING_KINDS:
            annotations[TRIMMING_RECOMMENDATIONS][trim_kind] = []

    trim_rec = annotations[TRIMMING_RECOMMENDATIONS]
    trim_rec[kind].extend(segments)


class TrimEdges(_BaseTrim):
    'It adds a trimming recommendation a fixed number of bases from the seqs.'
    def __init__(self, left=0, right=0):
        '''The initiator.

        left - number of bases to trim from the left side
        right - number of bases to trim from the right side
        mask - If True the edges will be masked instead of trimmed
        '''
        self.left = left
        self.right = right
        super(TrimEdges, self).__init__()

    def _do_trim(self, seq):
        'It trims the edges of the given seqs.'
        left = self.left
        right = self.right
        segments = [(0, left - 1)] if left else []
        if right:
            seq_len = get_length(seq)
            segments.append((seq_len - right, seq_len - 1))
        _add_trim_segments(segments, seq, kind=OTHER)
        return seq


def _mask_sequence(seq, segments):
    'It masks the given segments of the sequence'

    if not segments:
        return seq
    segments = merge_overlaping_segments(segments)
    segments = get_all_segments(segments, get_length(seq))
    str_seq = get_str_seq(seq)
    new_seq = ''
    for segment in segments:
        start = segment[0][0]
        end = segment[0][1] + 1
        str_seq_ = str_seq[start:end]

        if segment[1]:
            str_seq_ = str_seq_.lower()
        new_seq += str_seq_
    if seq.kind == SEQRECORD:
        new_seq = Seq(new_seq, alphabet=seq.object.seq.alphabet)
    return copy_seq(seq, seq=new_seq)


class TrimOrMask(object):
    'It trims and masks the Seq following the trimming recommendations.'
    def __init__(self, mask=False):
        '''The initiator.'''
        self.mask = mask

    def __call__(self, trim_packet):
        'It trims the seqs'
        trimmed_seqs = []
        orphan_seqs = trim_packet[ORPHAN_SEQS]
        for paired_seqs in trim_packet[SEQS_PASSED]:
            trimmed_paired_seqs = [self._do_trim(s) for s in paired_seqs]
            # all sequences are trimed, no lost
            if None not in trimmed_paired_seqs:
                trimmed_seqs.append(trimmed_paired_seqs)
            # all secuences are lost because of trimming
            elif (len(trimmed_paired_seqs) == 1 or
                  trimmed_paired_seqs == (None, None)):
                continue
            # one of the pairs is lost in trimming
            else:
                orphans = [s for s in trimmed_paired_seqs if s is not None]
                orphan_seqs.extend(orphans)
        orphan_seqs = self._trim_orphans(orphan_seqs)
        return {SEQS_PASSED: trimmed_seqs, ORPHAN_SEQS: orphan_seqs}

    def _do_trim(self, seq):
        'It trims the edges of the given seqs.'
        annots = get_annotations(seq)
        if not TRIMMING_RECOMMENDATIONS in annots:
            return seq

        trim_rec = annots[TRIMMING_RECOMMENDATIONS]
        # fixing the trimming recommendations
        if TRIMMING_RECOMMENDATIONS in annots:
            del annots[TRIMMING_RECOMMENDATIONS]

        trim_segments = []
        for trim_kind in TRIMMING_KINDS:
            trim_segments.extend(trim_rec.get(trim_kind, []))

        # masking
        if self.mask:
            seq = _mask_sequence(seq, trim_segments)
        else:
            # trimming
            if trim_segments:
                trim_limits = get_longest_complementary_segment(
                                            trim_segments, get_length(seq))
                if trim_limits is None:
                    # there's no sequence left
                    return None
            else:
                trim_limits = []

            if trim_limits:
                seq = slice_seq(seq, trim_limits[0], trim_limits[1] + 1)

        return seq

    def _trim_orphans(self, seqs):
        new_seqs = []
        for seq in seqs:
            seq = self._do_trim(seq)
            if seq is not None:
                new_seqs.append(seq)
        return new_seqs


def _get_bad_quality_segments(quals, window, threshold, trim_left=True,
                              trim_right=True):
    '''It returns the regions with quality above the threshold.

    The algorithm is similar to the one used by qclip in Staden.
    '''
    # do window quality means
    mean = lambda l: float(sum(l)) / len(l) if len(l) > 0 else float('nan')

    wquals = [mean(win_quals) for win_quals in rolling_window(quals, window)]

    if not wquals:
        return [(0, len(quals) - 1)]

    index_max, max_val = max(enumerate(wquals), key=itemgetter(1))

    if max_val < threshold:
        return [(0, len(quals) - 1)]

    if trim_left:
        wleft_index = 0
        for wleft_index in range(index_max - 1, -1, -1):
            if wquals[wleft_index] < threshold:
                wleft_index += 1
                break
    else:
        wleft_index = 0
    if trim_right:
        wright_index = index_max
        for wright_index in range(index_max, len(wquals)):
            if wquals[wright_index] < threshold:
                wright_index -= 1
                break
    else:
        wright_index = len(wquals) - 1
    left = wleft_index
    right = wright_index + window - 1
    segments = []
    if left:
        segments.append((0, left - 1))
    if right < len(quals) - 1:
        segments.append((right + 1, len(quals) - 1))
    if not segments:
        return None
    return segments


class TrimByQuality(_BaseTrim):
    'It trims the low quality regions of the SeqRecords.'

    def __init__(self, window, threshold, trim_left=True, trim_right=True):
        'The initiator'
        self.window = int(window)
        self.threshold = threshold
        self.trim_left = trim_left
        self.trim_right = trim_right
        super(TrimByQuality, self).__init__()

    def _do_trim(self, seq):
        'It trims the masked segments of the seqrecords.'
        window = self.window
        threshold = self.threshold
        trim_left = self.trim_left
        trim_right = self.trim_right
        try:
            quals = list(get_int_qualities(seq))
        except KeyError:
            msg = 'Some of the input sequences do not have qualities: {}'
            msg = msg.format(get_name(seq))
        segments = _get_bad_quality_segments(quals, window, threshold,
                                            trim_left, trim_right)
        if segments is not None:
            _add_trim_segments(segments, seq, kind=QUALITY)

        return seq


class TrimWithBlastShort(_BaseTrim):
    'It trims adaptors with the blast short algorithm'
    def __init__(self, oligos):
        'The initiator'
        self.oligos = oligos
        super(TrimWithBlastShort, self).__init__()

    def _pre_trim(self, trim_packet):
        seqs = [s for seqs in trim_packet[SEQS_PASSED]for s in seqs]
        db_fhand = write_seqs(seqs, file_format='fasta')
        db_fhand.flush()
        params = {'task': 'blastn-short', 'expect': '0.0001'}
        filters = [{'kind': 'score_threshold', 'score_key': 'identity',
                    'min_score': 87},
                   {'kind': 'min_length', 'min_num_residues': 13,
                    'length_in_query': False}]
        self._matcher = BlasterForFewSubjects(db_fhand.name, self.oligos,
                                             program='blastn', filters=filters,
                                             params=params,
                                             elongate_for_global=True)

    def _do_trim(self, seq):
        'It trims the masked segments of the SeqWrappers.'
        segments = self._matcher.get_matched_segments_for_read(get_name(seq))
        if segments is not None:
            _add_trim_segments(segments[0], seq, kind=VECTOR)
        return seq


# can not put in seqio due to cross imports

