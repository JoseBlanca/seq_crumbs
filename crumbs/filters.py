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

from crumbs.utils.tags import PROCESSED_PACKETS, PROCESSED_SEQS, YIELDED_SEQS
from crumbs.utils.seq_utils import uppercase_length


class FilterByLength(object):
    'It removes the sequences according to their length.'
    def __init__(self, threshold, reverse=False, ignore_masked=False):
        '''The initiator.

        threshold - minimum length to pass the filter (integer)
        reverse - if True keep the short sequences and discard the long ones
        ignore_masked - If True only uppercase letters will be counted.
        '''
        self.threshold = threshold
        self.reverse = reverse
        self.ignore_masked = ignore_masked
        self._stats = {PROCESSED_SEQS: 0,
                       PROCESSED_PACKETS: 0,
                       YIELDED_SEQS: 0}

    @property
    def stats(self):
        'The process stats'
        return self._stats

    def __call__(self, seqrecords):
        'It filters out the short seqrecords.'
        stats = self._stats
        threshold = self.threshold
        reverse = self.reverse
        ignore_masked = self.ignore_masked
        stats[PROCESSED_PACKETS] += 1
        processed_seqs = []
        for seqrecord in seqrecords:
            stats[PROCESSED_SEQS] += 1
            seq = str(seqrecord.seq)
            length = uppercase_length(seq) if ignore_masked else len(seq)
            passed = True if length >= threshold else False
            if reverse:
                passed = not(passed)
            if passed:
                processed_seqs.append(seqrecord)
                stats[YIELDED_SEQS] += 1
        return processed_seqs
