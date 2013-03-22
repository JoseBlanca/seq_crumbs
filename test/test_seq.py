# Copyright 2013 Jose Blanca, Peio Ziarsolo, COMAV-Univ. Politecnica Valencia
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

# pylint: disable=R0201
# pylint: disable=R0904
# pylint: disable=C0111

import unittest

from crumbs.seq import (get_length, get_str_seq, get_qualities, slice_seq,
                        copy_seq, SeqItem, SeqWrapper)
from crumbs.utils.tags import SEQITEM


class SeqMethodsTest(unittest.TestCase):
    def test_len(self):
        # with fasta
        seq = SeqItem(name='s1', lines=['>s1\n', 'ACTG\n', 'GTAC\n'])
        seq = SeqWrapper(SEQITEM, seq, 'fasta')
        assert get_length(seq) == 8

        # with fastq
        seq = SeqItem(name='seq',
                      lines=['@seq\n', 'aaaa\n', '+\n', '????\n'])
        seq = SeqWrapper(SEQITEM, seq, 'fastq')
        assert get_length(seq) == 4

    def test_str_seq(self):
        # with fasta
        seq = SeqItem(name='s1', lines=['>s1\n', 'ACTG\n', 'GTAC\n'])
        seq = SeqWrapper(SEQITEM, seq, 'fasta')
        assert get_str_seq(seq) == 'ACTGGTAC'

        # with fastq
        seq = SeqItem(name='seq',
                      lines=['@seq\n', 'aaaa\n', '+\n', '????\n'])
        seq = SeqWrapper(SEQITEM, seq, 'fastq')
        assert get_str_seq(seq) == 'aaaa'

        # with fastq-multiline
        seq = SeqItem(name='seq', lines=['@seq\n', 'aaaa\n', 'tttt\n',
                                         '+\n', '????\n', '????\n'])
        seq = SeqWrapper(SEQITEM, seq, 'fastq-multiline')
        assert get_str_seq(seq) == 'aaaatttt'

    def test_qualities(self):
        # with fasta
        seq = SeqItem(name='s1', lines=['>s1\n', 'ACTG\n', 'GTAC\n'])
        seq = SeqWrapper(SEQITEM, seq, 'fasta')
        try:
            assert get_qualities(seq)
            self.fail('AttributeError expected')
        except AttributeError:
            pass

        # with fastq
        seq = SeqItem(name='seq',
                      lines=['@seq\n', 'aaaa\n', '+\n', '!???\n'])
        seq = SeqWrapper(SEQITEM, seq, 'fastq')
        assert list(get_qualities(seq)) == [0, 30, 30, 30]

        # with multiline fastq
        seq = SeqItem(name='seq', lines=['@seq\n', 'aaaa\n', 'aaaa\n', '+\n',
                                         '@AAA\n', 'BBBB\n'])
        seq = SeqWrapper(SEQITEM, seq, 'fastq-illumina-multiline')
        assert list(get_qualities(seq)) == [0, 1, 1, 1, 2, 2, 2, 2]

    def test_slice(self):
        # with fasta
        seq = SeqItem(name='s1', lines=['>s1\n', 'ACTG\n', 'GTAC\n'])
        seq = SeqWrapper(SEQITEM, seq, 'fasta')
        expected_seq = SeqItem(name='s1', lines=['>s1\n', 'CTGG\n'])
        expected_seq = SeqWrapper(SEQITEM, expected_seq, 'fasta')
        assert slice_seq(seq, 1, 5) == expected_seq

        # with fastq
        seq = SeqItem(name='seq',
                      lines=['@seq\n', 'aata\n', '+\n', '!?!?\n'])
        seq = SeqWrapper(SEQITEM, seq, 'fastq')
        seq = slice_seq(seq, 1, 3)
        assert list(get_qualities(seq)) == [30, 0]
        assert get_str_seq(seq) == 'at'
        assert seq.object.lines == ['@seq\n', 'at\n', '+\n', '?!\n']

        # with multiline fastq
        seq = SeqItem(name='seq', lines=['@seq\n', 'aaat\n', 'caaa\n', '+\n',
                                         '@AAA\n', 'BBBB\n'])
        seq = SeqWrapper(SEQITEM, seq, 'fastq-illumina-multiline')
        seq_ = slice_seq(seq, 1, 5)
        assert list(get_qualities(seq_)) == [1, 1, 1, 2]
        assert get_str_seq(seq_) == get_str_seq(seq)[1: 5]

        # It tests the stop is None
        seq = SeqItem('seq', ['>seq\n', 'aCTG'])
        seq = SeqWrapper(SEQITEM, seq, 'fasta')
        assert get_str_seq(slice_seq(seq, 1, None)) == 'aCTG'[1:]

        assert get_str_seq(slice_seq(seq, None, 1)) == 'aCTG'[:1]

    def test_copy(self):
        # with fasta
        seq = SeqItem(name='s1', lines=['>s1\n', 'ACTG\n', 'GTAC\n'],
                      annotations={'a': 'b'})
        seq = SeqWrapper(SEQITEM, seq, 'fasta')
        seq2 = copy_seq(seq, seq='ACTG')
        assert seq2.object == SeqItem(name='s1', lines=['>s1\n', 'ACTG\n'],
                                      annotations={'a': 'b'})
        assert seq.object is not seq2.object
        assert seq.object.lines is not seq2.object.lines

        # with fastq
        seq = SeqItem(name='seq',
                      lines=['@seq\n', 'aaaa\n', '+\n', '!???\n'])
        seq = SeqWrapper(SEQITEM, seq, 'fastq')
        seq2 = copy_seq(seq, seq='ACTG')
        assert seq2.object == SeqItem(name='seq',
                               lines=['@seq\n', 'ACTG\n', '+\n', '!???\n'])

        # with multiline fastq
        seq = SeqItem(name='seq', lines=['@seq\n', 'aaaa\n', 'aaaa\n', '+\n',
                                         '@AAA\n', 'BBBB\n'])
        seq = SeqWrapper(SEQITEM, seq, 'fastq-illumina-multiline')
        seq2 = copy_seq(seq, seq='ACTGactg')
        assert seq2.object == SeqItem(name='seq',
                                      lines=['@seq\n', 'ACTGactg\n', '+\n',
                                             '@AAABBBB\n'])

    def test_change_name(self):
        seq = SeqItem(name='seq',
                      lines=['@seq\n', 'aaaa\n', '+seq\n', '!???\n'])
        seq = SeqWrapper(SEQITEM, seq, 'fastq')
        seq = copy_seq(seq, name='seq2')
        assert seq.object == ('seq2', ['@seq2\n', 'aaaa\n', '+\n', '!???\n'],
                              {})

        seq = SeqItem(name='seq', lines=['>seq\n', 'aaaa\n'])
        seq = SeqWrapper(SEQITEM, seq, 'fasta')
        seq = copy_seq(seq, name='seq2')
        assert seq.object == ('seq2', ['>seq2\n', 'aaaa\n'],
                              {})

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'ChangeCaseTest.test_bin']
    unittest.main()
