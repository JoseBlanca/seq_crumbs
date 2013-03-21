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

from crumbs.utils.tags import (TRIMMING_RECOMMENDATIONS, QUALITY, OTHER,
                               VECTOR, TRIMMING_KINDS)
from crumbs.utils.seq_utils import get_uppercase_segments
from crumbs.seq import (copy_seq, get_str_seq, get_annotations, get_length,
                        slice_seq, get_qualities, get_name)
from crumbs.utils.segments_utils import (get_longest_segment, get_all_segments,
                                         get_longest_complementary_segment,
                                         merge_overlaping_segments)
from crumbs.utils.tags import SEQRECORD
from crumbs.iterutils import rolling_window
from crumbs.blast import BlasterForFewSubjects
from crumbs.seqio import write_seqs

# pylint: disable=R0903


class TrimLowercasedLetters(object):
    'It trims the masked segments of the seqrecords.'

    def __call__(self, seqs):
        'It trims the masked segments of the seqrecords.'
        trimmed_seqs = []
        for seq in seqs:
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
                trimmed_seqs.append(seq)
        return trimmed_seqs


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

    def __call__(self, seqs):
        'It trims the edges of the given seqs.'
        left = self.left
        right = self.right
        for seq in seqs:
            segments = [(0, left - 1)] if left else []
            if right:
                seq_len = get_length(seq)
                segments.append((seq_len - right, seq_len - 1))
            _add_trim_segments(segments, seq, kind=OTHER)
        return seqs


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

    def __call__(self, seqs):
        'It trims the edges of the given seqs.'
        mask = self.mask
        processed_seqs = []
        for seq in seqs:
            annots = get_annotations(seq)
            if not TRIMMING_RECOMMENDATIONS in annots:
                processed_seqs.append(copy_seq(seq))
                continue

            trim_rec = annots[TRIMMING_RECOMMENDATIONS]
            # fixing the trimming recommendations
            if TRIMMING_RECOMMENDATIONS in annots:
                del annots[TRIMMING_RECOMMENDATIONS]

            trim_segments = []
            for trim_kind in TRIMMING_KINDS:
                trim_segments.extend(trim_rec.get(trim_kind, []))

            # masking
            if mask:
                seq = _mask_sequence(seq, trim_segments)
            else:
                # trimming
                if trim_segments:
                    trim_limits = get_longest_complementary_segment(
                                                trim_segments, get_length(seq))
                    if trim_limits is None:
                        # there's no sequence left
                        continue
                else:
                    trim_limits = []

                if trim_limits:
                    seq = slice_seq(seq, trim_limits[0], trim_limits[1] + 1)

            processed_seqs.append(seq)

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

    def __call__(self, seqs):
        'It trims the masked segments of the seqrecords.'
        window = self.window
        threshold = self.threshold
        trim_left = self.trim_left
        trim_right = self.trim_right
        trimmed_seqs = []
        for seq in seqs:
            try:
                quals = list(get_qualities(seq))
            except KeyError:
                msg = 'Some of the input sequences do not have qualities: {}'
                msg = msg.format(get_name(seq))
            segments = _get_bad_quality_segments(quals, window, threshold,
                                                trim_left, trim_right)
            if segments is not None:
                _add_trim_segments(segments, seq, kind=QUALITY)
            trimmed_seqs.append(seq)
        return trimmed_seqs


class TrimWithBlastShort(object):
    'It trims adaptors with the blast short algorithm'
    def __init__(self, oligos):
        'The initiator'
        self.oligos = oligos

    def __call__(self, seqs):
        'It trims the masked segments of the SeqWrappers.'
        db_fhand = write_seqs(seqs, file_format='fasta')
        db_fhand.flush()
        params = {'task': 'blastn-short', 'expect': '0.0001'}
        filters = [{'kind': 'score_threshold', 'score_key': 'identity',
                    'min_score': 89},
                   {'kind': 'min_length', 'min_num_residues': 13,
                    'length_in_query': False}]
        matcher = BlasterForFewSubjects(db_fhand.name, self.oligos,
                                        program='blastn', filters=filters,
                                        params=params,
                                        elongate_for_global=True)
        for seq in seqs:
            segments = matcher.get_matched_segments_for_read(get_name(seq))
            if segments is not None:
                _add_trim_segments(segments[0], seq, kind=VECTOR)
        return seqs
