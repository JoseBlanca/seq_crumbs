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

from collections import Counter

from crumbs.utils.tags import PROCESSED_PACKETS, PROCESSED_SEQS, YIELDED_SEQS
from crumbs.utils.seq_utils import uppercase_length, get_uppercase_segments
from crumbs.exceptions import WrongFormatError
from crumbs.blast import Blaster
from crumbs.iterutils import rolling_window
from crumbs.settings import DUST_WINDOWSIZE as WINDOWSIZE
from crumbs.settings import DUST_WINDOWSTEP as WINDOWSTEP
from crumbs.settings import DEFATULT_DUST_THRESHOLD

# pylint: disable=R0903


class FilterByLength(object):
    'It removes the sequences according to their length.'
    def __init__(self, minimum=None, maximum=None, ignore_masked=False):
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
        min_ = self.min
        max_ = self.max
        ignore_masked = self.ignore_masked
        stats[PROCESSED_PACKETS] += 1
        processed_seqs = []
        for seqrecord in seqrecords:
            stats[PROCESSED_SEQS] += 1
            seq = str(seqrecord.seq)
            length = uppercase_length(seq) if ignore_masked else len(seq)
            passed = True
            if min_ is not None and length < min_:
                passed = False
            if max_ is not None and length > max_:
                passed = False

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
        self._dbtype = dbtype
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
        matcher = Blaster(seqrecords, self._blast_db, dbtype=self._dbtype,
                          program=self._blast_program, filters=self._filters)
        for seqrec in seqrecords:
            stats[PROCESSED_SEQS] += 1
            segments = matcher.get_matched_segments(seqrec.id)
            if ((not self._reverse and segments is None) or
                (self._reverse and segments)):
                filtered_seqrecords.append(seqrec)
                stats[YIELDED_SEQS] += 1

        return filtered_seqrecords


def _calculate_rawscore(string):
    'It returns a non-normalized dustscore'
    triplet_counts = Counter()
    for triplet in rolling_window(string, 3):
        # It should do something with non ATCG, but we sacrifice purity for
        # speed. Maybe we should reconsider this
        triplet_counts[triplet.upper()] += 1

    return sum(tc * (tc - 1) * 0.5 for tc in triplet_counts.viewvalues())


def _calculate_dust_score(seqrecord):
    '''It returns the dust score.

    From: "A Fast and Symmetric DUST Implementation to Mask Low-Complexity DNA
    Sequences"
    doi:10.1089/cmb.2006.13.1028

    and re-implemented from PRINSEQ
    '''
    seq = str(seqrecord.seq)
    length = len(seq)
    if length == 3:
        return 0
    if length <= 5:
        return None

    dustscores = []
    if length > WINDOWSIZE:
        windows = 0
        for seq_in_win in rolling_window(seq, WINDOWSIZE, WINDOWSTEP):
            score = _calculate_rawscore(seq_in_win)
            dustscores.append(score / (WINDOWSIZE - 2))
            windows += 1
        #remaining_seq = seq[length - windows * WINDOWSTEP:]
        #remaining = len(seq) - windows * WINDOWSTEP
        remaining_seq = seq[windows * WINDOWSTEP:]
        #print 'remaining', length - windows * WINDOWSTEP, remaining_seq
    else:
        remaining_seq = seq

    if remaining_seq > 5:
        length = len(remaining_seq)
        score = _calculate_rawscore(remaining_seq)
        dustscore = score / (length - 3) * (WINDOWSIZE - 2) / (length - 2)
        dustscores.append(dustscore)

    # max score should be 100 not 31
    dustscore = sum(dustscores) / len(dustscores) * 100 / 31
    return dustscore


class FilterDustComplexity(object):
    'It filters a sequence according to its dust score'
    def __init__(self, threshold=DEFATULT_DUST_THRESHOLD, reverse=False):
        '''The initiator
        '''
        self._threshold = threshold
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
        filtered_seqs = []
        stats = self._stats
        threshold = self._threshold
        reverse = self._reverse
        for seqrec in seqrecords:
            stats[PROCESSED_SEQS] += 1
            dustscore = _calculate_dust_score(seqrec)
            passed = True if dustscore < threshold else False
            if reverse:
                passed = not(passed)
            if passed:
                filtered_seqs.append(seqrec)
                stats[YIELDED_SEQS] += 1
        return filtered_seqs
