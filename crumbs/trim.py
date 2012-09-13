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

from crumbs.utils.tags import (PROCESSED_PACKETS, PROCESSED_SEQS, YIELDED_SEQS,
                               TRIMMING_RECOMMENDATIONS)
from crumbs.utils.seq_utils import (replace_seq_same_length,
                                    get_uppercase_segments)
from crumbs.utils.segments_utils import (get_longest_segment, get_all_segments,
                                         get_longest_complementary_segment,
                                         merge_overlaping_segments)

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
                trimmed_seqs.append(seqrecord[segment[0]:segment[1] + 1])
        return trimmed_seqs


def _add_trim_segments(segments, sequence, vector=True, trim=True):
    'It adds segments to the trimming recommendation in the annotation'
    if not segments:
        return
    if TRIMMING_RECOMMENDATIONS not in sequence.annotations:
        sequence.annotations[TRIMMING_RECOMMENDATIONS] = {'vector': [],
                                                           'quality': [],
                                                           'mask': []}
    trim_rec = sequence.annotations[TRIMMING_RECOMMENDATIONS]
    if vector and trim:
        trim_rec['vector'].extend(segments)
    elif not vector and trim:
        trim_rec['quality'].extend(segments)
    elif not trim:
        trim_rec['mask'].extend(segments)


class TrimEdges(object):
    'It adds a trimming recommendation a fixed number of bases from the seqs.'
    def __init__(self, left=0, right=0, mask=False):
        '''The initiator.

        left - number of bases to trim from the left side
        right - number of bases to trim from the right side
        mask - If True the edges will be masked instead of trimmed
        '''
        self.left = left
        self.right = right
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
        left = self.left
        right = self.right
        trim = not(self.mask)
        stats[PROCESSED_PACKETS] += 1
        for seqrecord in seqrecords:
            stats[PROCESSED_SEQS] += 1

            segments = [(0, left - 1)] if left else []
            if right:
                seq_len = len(seqrecord)
                segments.append((seq_len - right, seq_len - 1))
            _add_trim_segments(segments, seqrecord, vector=False, trim=trim)
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
    return replace_seq_same_length(seqrecord, new_seq)


class TrimAndMask(object):
    'It trims and masks the seqrecords following the trimming recommendations.'

    def __init__(self):
        '''The initiator.'''
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
        processed_seqs = []
        for seqrecord in seqrecords:
            stats[PROCESSED_SEQS] += 1

            if not TRIMMING_RECOMMENDATIONS in seqrecord.annotations:
                processed_seqs.append(seqrecord)
                continue

            trim_rec = seqrecord.annotations[TRIMMING_RECOMMENDATIONS]

            #masking
            segments = trim_rec.get('mask', [])
            seqrecord = _mask_sequence(seqrecord, segments)

            #trimming
            trim_segments = trim_rec.get('vector', [])
            trim_segments.extend(trim_rec.get('quality', []))
            if trim_segments:
                trim_limits = get_longest_complementary_segment(trim_segments,
                                                                len(seqrecord))
                if trim_limits is None:
                    # there's no sequence left
                    continue
            else:
                trim_limits = []

            if trim_limits:
                seqrecord = seqrecord[trim_limits[0]:trim_limits[1] + 1]

            #fixing the trimming recommendations
            seqrecord.annotations[TRIMMING_RECOMMENDATIONS] = {}
            processed_seqs.append(seqrecord)

            stats[YIELDED_SEQS] += 1
        return processed_seqs
