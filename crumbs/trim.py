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

'''
Created on 2012 eka 26

@author: peio
'''
import random
import itertools

from crumbs.settings import PROCESSED_PACKETS, PROCESSED_SEQS, YIELDED_SEQS

# pylint: disable=R0903


def _get_uppercase_segments(string):
    '''It detects the unmasked regions of a sequence

    It returns a list of (start, end) tuples'''
    start = 0
    for is_upper, group in itertools.groupby(string, lambda x: x.isupper()):
        group = list(group)
        end = start + len(group) - 1
        if is_upper:
            yield start, end
        start = end + 1


def _get_longest_segment(segments):
    'it returns the longest segment'
    longest = None
    longest_size = None
    for segment in segments:
        size = segment[1] - segment[0]
        if longest is None:
            longest = [segment]
            longest_size = size
        elif size > longest_size:
            longest = [segment]
            longest_size = size
        elif size == longest_size:
            longest.append(segment)

    if longest is None:
        return None
    elif len(longest) == 1:
        return longest[0]
    else:
        return random.choice(longest)


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
            unmasked_segments = _get_uppercase_segments(seq)
            segment = _get_longest_segment(unmasked_segments)
            if segment is not None:
                stats[YIELDED_SEQS] += 1
                trimmed_seqs.append(seqrecord[segment[0]:segment[1] + 1])
        return trimmed_seqs


class TrimEdges(object):
    'It trims a fixed number of bases from the seqrecords.'
    def __init__(self, left=0, right=0, mask=False):
        '''The initiator.

        left - number of bases to trim from the left side
        right - number of bases to trim from the right side
        mask - if True it will mask by lowercasing instead of trimming.
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
        mask = self.mask
        stats[PROCESSED_PACKETS] += 1
        processed_seqs = []
        for seqrecord in seqrecords:
            stats[PROCESSED_SEQS] += 1
            right = self.right
            if mask:
                seq = str(seqrecord.seq)
                new_seq = []
                if left:
                    new_seq.append(seq[:left].lower())
                max_right = len(seq) - left
                if right > max_right:
                    right = max_right
                if right:
                    new_seq.append(seq[left:-right])
                    new_seq.append(seq[-right:].lower())
                else:
                    new_seq.append(seq[left:])
                annots = seqrecord.letter_annotations
                seqrecord.letter_annotations = {}
                seqrecord.seq = ''.join(new_seq)
                seqrecord.letter_annotations = annots
            else:
                if right:
                    seqrecord = seqrecord[left:-right]
                else:
                    seqrecord = seqrecord[left:]
            if len(seqrecord.seq):
                processed_seqs.append(seqrecord)
                stats[YIELDED_SEQS] += 1
        return processed_seqs
