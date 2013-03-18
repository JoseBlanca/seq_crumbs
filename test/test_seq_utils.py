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

# pylint: disable=R0201
# pylint: disable=R0904
# pylint: disable=C0111

import unittest
from cStringIO import StringIO
from tempfile import NamedTemporaryFile
import os.path
from subprocess import check_output

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from crumbs.utils.file_utils import fhand_is_seekable, wrap_in_buffered_reader
from crumbs.utils.seq_utils import (uppercase_length, ChangeCase,
                                    get_uppercase_segments, get_length,
                                    get_str_seq, get_qualities, slice_seq,
                                    copy_seq)
from crumbs.utils.tags import SWAPCASE, UPPERCASE, LOWERCASE, SEQITEM
from crumbs.utils.bin_utils import BIN_DIR
from crumbs.seqio import SeqItem, SeqWrapper


class SeekableFileTest(unittest.TestCase):

    def test_is_seekable(self):
        'It tests wether the fhands are seekable or not'

        # StringIO
        fhand = StringIO('hola')
        assert fhand_is_seekable(fhand)

        # standard file
        fhand = NamedTemporaryFile()
        fhand.seek(0)
        assert fhand_is_seekable(fhand)

        # a wrapped BufferedReader
        fhand2 = wrap_in_buffered_reader(fhand)
        assert fhand_is_seekable(fhand2)

        # pylint: disable=R0903
        # pylint: disable=C0111
        class NonSeekable(object):
            'Just for testing'
            pass

        assert not fhand_is_seekable(NonSeekable())

        class NonSeekable2(object):
            def seek(self):
                pass

            def seekable(self):
                return False

        assert not fhand_is_seekable(NonSeekable2())


class UppercaseLengthTest(unittest.TestCase):
    'It tests the uppercase character count'
    def test_uppercase_length(self):
        'It counts the number of uppercase letters in a string'
        assert uppercase_length('aCTaGGt') == 4
        assert uppercase_length('acagt') == 0


def _make_fhand(content=''):
    'It makes temporary fhands'
    fhand = NamedTemporaryFile()
    fhand.write(content)
    fhand.flush()
    return fhand


class MaskedSegmentsTest(unittest.TestCase):
    'It tests the lower case segments location functions'

    @staticmethod
    def test_masked_locations():
        'It test the masked locations function'

        assert list(get_uppercase_segments('aaATTTTTTaa')) == [(2, 8)]

        assert list(get_uppercase_segments('aaATTTaTTaa')) == [(2, 5), (7, 8)]

        assert list(get_uppercase_segments('AAATaaa')) == [(0, 3)]

        assert list(get_uppercase_segments('aaaaAAAA')) == [(4, 7)]

        seq = 'AATTaaTTaaTTT'
        assert list(get_uppercase_segments(seq)) == [(0, 3), (6, 7), (10, 12)]

        assert list(get_uppercase_segments('AATT')) == [(0, 3)]
        assert not list(get_uppercase_segments('aatt'))


class ChangeCaseTest(unittest.TestCase):
    'It tests the case change'
    def test_case_change(self):
        'It changes the case of the sequences'
        seqs = [SeqRecord(Seq('aCCg'), letter_annotations={'dummy': 'dddd'})]
        change_case = ChangeCase(action=UPPERCASE)
        strs = [str(s.seq) for s in change_case(seqs)]
        assert strs == ['ACCG']

        seqs = [SeqRecord(Seq('aCCg'))]
        change_case = ChangeCase(action=LOWERCASE)
        strs = [str(s.seq) for s in change_case(seqs)]
        assert strs == ['accg']

        seqs = [SeqRecord(Seq('aCCg'))]
        change_case = ChangeCase(action=SWAPCASE)
        strs = [str(s.seq) for s in change_case(seqs)]
        assert strs == ['AccG']

    def test_bin(self):
        'It tests the trim seqs binary'
        change_bin = os.path.join(BIN_DIR, 'change_case')
        assert 'usage' in check_output([change_bin, '-h'])

        fastq = '@seq1\naTCgt\n+\n?????\n@seq2\natcGT\n+\n?????\n'
        fastq_fhand = _make_fhand(fastq)

        result = check_output([change_bin, '-a', 'upper', fastq_fhand.name])
        assert '@seq1\nATCGT\n+' in result


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
            self.fail('ValueError expected')
        except ValueError:
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

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'ChangeCaseTest.test_bin']
    unittest.main()
