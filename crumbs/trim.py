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

from Bio.Seq import Seq

from crumbs.utils.tags import (PROCESSED_PACKETS, PROCESSED_SEQS, YIELDED_SEQS,
                               TRIMMING_RECOMMENDATIONS, QUALITY, OTHER,
                               VECTOR, TRIMMING_KINDS)
from crumbs.utils.seq_utils import copy_seqrecord, get_uppercase_segments
from crumbs.utils.segments_utils import (get_longest_segment, get_all_segments,
                                         get_longest_complementary_segment,
                                         merge_overlaping_segments)
from crumbs.iterutils import rolling_window
from crumbs.blast import BlastMatcher
from crumbs.seqio import write_seqrecords

# pylint: disable=R0903


class TrimLowercasedLetters(object):
    'It trims the masked segments of the seqrecords.'

    def __init__(self):
        'The initiator'
        self._stats = {PROCESSED_SEQS: 0,
                       PROCESSED_PACKETS: 0,
                       YIELDED_SEQS: 0}

    @property
    def stats(self):
        'The process stats'
        return self._stats

    def __call__(self, seqrecords):
        'It trims the masked segments of the seqrecords.'
        stats = self._stats
        stats[PROCESSED_PACKETS] += 1
        trimmed_seqs = []
        for seqrecord in seqrecords:
            stats[PROCESSED_SEQS] += 1
            seq = str(seqrecord.seq)
            unmasked_segments = get_uppercase_segments(seq)
            segment = get_longest_segment(unmasked_segments)
            if segment is not None:
                stats[YIELDED_SEQS] += 1
                segments = []
                if segment[0] != 0:
                    segments.append((0, segment[0] - 1))
                len_seq = len(seqrecord)
                if segment[1] != len_seq - 1:
                    segments.append((segment[1] + 1, len_seq - 1))

                _add_trim_segments(segments, seqrecord, kind=OTHER)
                trimmed_seqs.append(seqrecord)
            #trimmed_seqs.append(seqrecord[segment[0]:segment[1] + 1])
        return trimmed_seqs


def _add_trim_segments(segments, sequence, kind):
    'It adds segments to the trimming recommendation in the annotation'
    assert kind in TRIMMING_KINDS
    if not segments:
        return
    if TRIMMING_RECOMMENDATIONS not in sequence.annotations:
        sequence.annotations[TRIMMING_RECOMMENDATIONS] = {}
        for trim_kind in TRIMMING_KINDS:
            sequence.annotations[TRIMMING_RECOMMENDATIONS][trim_kind] = []

    trim_rec = sequence.annotations[TRIMMING_RECOMMENDATIONS]
    trim_rec[kind].extend(segments)


class TrimEdges(object):
    'It adds a trimming recommendation a fixed number of bases from the seqs.'
    def __init__(self, left=0, right=0):
        '''The initiator.

        left - number of bases to trim from the left side
        right - number of bases to trim from the right side
        mask - If True the edges will be masked instead of trimmed
        '''
        self.left = left
        self.right = right
        self._stats = {PROCESSED_SEQS: 0,
                       PROCESSED_PACKETS: 0,
                       YIELDED_SEQS: 0}

    @property
    def stats(self):
        'The process stats'
        return self._stats

    def __call__(self, seqrecords):
        'It trims the edges of the given seqrecords.'
        stats = self._stats
        left = self.left
        right = self.right
        stats[PROCESSED_PACKETS] += 1
        for seqrecord in seqrecords:
            stats[PROCESSED_SEQS] += 1

            segments = [(0, left - 1)] if left else []
            if right:
                seq_len = len(seqrecord)
                segments.append((seq_len - right, seq_len - 1))
            _add_trim_segments(segments, seqrecord, kind=OTHER)
            stats[YIELDED_SEQS] += 1
        return seqrecords


def _mask_sequence(seqrecord, segments):
    'It masks the given segments of the sequence'

    if not segments:
        return seqrecord
    segments = merge_overlaping_segments(segments)
    segments = get_all_segments(segments, len(seqrecord))
    seq = str(seqrecord.seq)
    new_seq = ''
    for segment in segments:
        start = segment[0][0]
        end = segment[0][1] + 1
        seq_ = seq[start:end]

        if segment[1]:
            seq_ = seq_.lower()
        new_seq += seq_
    return copy_seqrecord(seqrecord, seq=Seq(new_seq,
                                             alphabet=seqrecord.seq.alphabet))


class TrimOrMask(object):
    'It trims and masks the SeqRecords following the trimming recommendations.'

    def __init__(self, mask=False):
        '''The initiator.'''
        self.mask = mask
        self._stats = {PROCESSED_SEQS: 0,
                       PROCESSED_PACKETS: 0,
                       YIELDED_SEQS: 0}

    @property
    def stats(self):
        'The process stats'
        return self._stats

    def __call__(self, seqrecords):
        'It trims the edges of the given seqrecords.'
        stats = self._stats
        stats[PROCESSED_PACKETS] += 1
        mask = self.mask
        processed_seqs = []
        for seqrecord in seqrecords:
            stats[PROCESSED_SEQS] += 1

            if not TRIMMING_RECOMMENDATIONS in seqrecord.annotations:
                processed_seqs.append(copy_seqrecord(seqrecord))
                continue

            trim_rec = seqrecord.annotations[TRIMMING_RECOMMENDATIONS]
            #fixing the trimming recommendations
            if TRIMMING_RECOMMENDATIONS in seqrecord.annotations:
                del seqrecord.annotations[TRIMMING_RECOMMENDATIONS]

            trim_segments = []
            for trim_kind in TRIMMING_KINDS:
                trim_segments.extend(trim_rec.get(trim_kind, []))

            #masking
            if mask:
                seqrecord = _mask_sequence(seqrecord, trim_segments)
            else:
                #trimming
                if trim_segments:
                    trim_limits = get_longest_complementary_segment(
                                                 trim_segments, len(seqrecord))
                    if trim_limits is None:
                        # there's no sequence left
                        continue
                else:
                    trim_limits = []

                if trim_limits:
                    seqrecord = seqrecord[trim_limits[0]:trim_limits[1] + 1]

            processed_seqs.append(seqrecord)

            stats[YIELDED_SEQS] += 1
        return processed_seqs


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


class TrimByQuality(object):
    'It trims the low quality regions of the SeqRecords.'

    def __init__(self, window, threshold, trim_left=True, trim_right=True):
        'The initiator'
        self.window = int(window)
        self.threshold = threshold
        self.trim_left = trim_left
        self.trim_right = trim_right
        self._stats = {PROCESSED_SEQS: 0,
                       PROCESSED_PACKETS: 0,
                       YIELDED_SEQS: 0}

    @property
    def stats(self):
        'The process stats'
        return self._stats

    def __call__(self, seqrecords):
        'It trims the masked segments of the seqrecords.'
        window = self.window
        threshold = self.threshold
        trim_left = self.trim_left
        trim_right = self.trim_right
        stats = self._stats
        stats[PROCESSED_PACKETS] += 1
        trimmed_seqs = []
        for seqrecord in seqrecords:
            stats[PROCESSED_SEQS] += 1
            try:
                quals = seqrecord.letter_annotations['phred_quality']
            except KeyError:
                msg = 'Some of the input sequences do not have qualities: {}'
                msg = msg.format(seqrecord.id)
            segments = _get_bad_quality_segments(quals, window, threshold,
                                                trim_left, trim_right)
            if segments is not None:
                _add_trim_segments(segments, seqrecord, kind=QUALITY)
            stats[YIELDED_SEQS] += 1
            trimmed_seqs.append(seqrecord)
        return trimmed_seqs


class TrimWithBlastShort(object):
    'It trims adaptors with the blast short algorithm'
    def __init__(self, oligos):
        'The initiator'
        self.oligos = oligos
        self._stats = {PROCESSED_SEQS: 0,
                       PROCESSED_PACKETS: 0,
                       YIELDED_SEQS: 0}

    @property
    def stats(self):
        'The process stats'
        return self._stats

    def __call__(self, seqrecords):
        'It trims the masked segments of the seqrecords.'
        stats = self.stats
        db_fhand = write_seqrecords(seqrecords, file_format='fasta')
        db_fhand.flush()
        params = {'task': 'blastn-short', 'expect': '0.0001'}
        filters = [{'kind': 'score_threshold', 'score_key': 'identity',
                    'min_score': 89},
                   {'kind': 'min_length', 'min_num_residues': 13,
                    'length_in_query': False}]
        matcher = BlastMatcher(open(db_fhand.name), self.oligos,
                               program='blastn', filters=filters,
                               params=params, elongate_for_global=True)
        for seqrec in seqrecords:
            stats[PROCESSED_SEQS] += 1
            segments = matcher.get_matched_segments_for_read(seqrec.id)
            if segments is not None:
                _add_trim_segments(segments[0], seqrec, kind=VECTOR)
            stats[YIELDED_SEQS] += 1
        return seqrecords
