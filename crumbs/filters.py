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

from __future__ import division

from crumbs.utils.tags import PROCESSED_PACKETS, PROCESSED_SEQS, YIELDED_SEQS
from crumbs.utils.seq_utils import uppercase_length, get_uppercase_segments
from crumbs.exceptions import WrongFormatError
from crumbs.blast import Blaster
# pylint: disable=R0903


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


class FilterById(object):
    'It removes the sequences not found in the given set'
    def __init__(self, seq_ids, reverse=False):
        '''The initiator.

        seq_ids - An iterator with the sequence ids to keep
        reverse - if True keep the sequences not found on the list
        '''
        if not isinstance(seq_ids, set):
            seq_ids = set(seq_ids)
        self.seq_ids = seq_ids
        self.reverse = reverse
        self._stats = {PROCESSED_SEQS: 0,
                       PROCESSED_PACKETS: 0,
                       YIELDED_SEQS: 0}

    @property
    def stats(self):
        'The process stats'
        return self._stats

    def __call__(self, seqrecords):
        'It filters out the seqrecords not found in the list.'
        stats = self._stats
        seq_ids = self.seq_ids
        reverse = self.reverse
        stats[PROCESSED_PACKETS] += 1
        processed_seqs = []
        for seqrecord in seqrecords:
            stats[PROCESSED_SEQS] += 1
            passed = True if seqrecord.id in seq_ids else False
            if reverse:
                passed = not(passed)
            if passed:
                processed_seqs.append(seqrecord)
                stats[YIELDED_SEQS] += 1
        return processed_seqs


class FilterByQuality(object):
    'It removes the sequences according to its quality'
    def __init__(self, threshold, reverse=False, ignore_masked=False):
        '''The initiator.

        threshold - minimum quality to pass the filter (float)
        reverse - if True keep the sequences not found on the list
        '''
        self.threshold = float(threshold)
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
        'It filters out the seqrecords not found in the list.'
        stats = self._stats
        threshold = self.threshold
        reverse = self.reverse
        ignore_masked = self.ignore_masked

        stats[PROCESSED_PACKETS] += 1
        processed_seqs = []
        for seqrecord in seqrecords:
            stats[PROCESSED_SEQS] += 1
            try:
                quals = seqrecord.letter_annotations['phred_quality']
            except KeyError:
                msg = 'Some of the input sequences do not have qualities: {}'
                msg = msg.format(seqrecord.id)
                raise WrongFormatError(msg)
            if ignore_masked:
                seq = str(seqrecord.seq)
                seg_quals = [quals[segment[0]: segment[1] + 1]
                                    for segment in get_uppercase_segments(seq)]
                qual = sum(sum(q) * len(q) for q in seg_quals) / len(quals)
            else:
                qual = sum(quals) / len(quals)
            passed = True if qual >= threshold else False
            if reverse:
                passed = not(passed)
            if passed:
                processed_seqs.append(seqrecord)
                stats[YIELDED_SEQS] += 1
        return processed_seqs


class FilterBlastMatch(object):
    'It filters a seq if there is a match against a blastdb'
    def __init__(self, database, program, filters, dbtype=None, reverse=False):
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
        self._reverse = reverse
        self._stats = {PROCESSED_SEQS: 0,
                       PROCESSED_PACKETS: 0,
                       YIELDED_SEQS: 0}

    @property
    def stats(self):
        'The process stats'
        return self._stats

    def __call__(self, seqrecords):
        'It filters the seq by blast match'
        filtered_seqrecords = []
        stats = self._stats

        matcher = Blaster(seqrecords, self._blast_db,
                               program=self._blast_program,
                               filters=self._filters)
        for seqrec in seqrecords:
            stats[PROCESSED_SEQS] += 1
            segments = matcher.get_matched_segments(seqrec.id)
            if ((not self._reverse and segments is None) or
                (self._reverse and segments)):
                filtered_seqrecords.append(seqrec)
                stats[YIELDED_SEQS] += 1

        return filtered_seqrecords
