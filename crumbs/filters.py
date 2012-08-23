'''
Created on 23/08/2012

@author: jose
'''

from crumbs.settings import PROCESSED_PACKETS, PROCESSED_SEQS, YIELDED_SEQS
from crumbs.seq_utils import uppercase_length


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
