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

import unittest
from cStringIO import StringIO
from tempfile import NamedTemporaryFile
import os.path
from subprocess import check_output

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from crumbs.utils.file_utils import fhand_is_seekable, wrap_in_buffered_reader
from crumbs.utils.seq_utils import (uppercase_length, guess_format, ChangeCase,
                                    _guess_format, get_uppercase_segments)
from crumbs.utils.tags import SWAPCASE, UPPERCASE, LOWERCASE
from crumbs.exceptions import UnknownFormatError, UndecidedFastqVersionError
from crumbs.utils.bin_utils import BIN_DIR
from crumbs.seqio import guess_seq_type


class GuessFormatTest(unittest.TestCase):
    'It tests the function that guess the sequence format'

    def test_fasta(self):
        'It guess fasta formats'
        fhand = StringIO('>seq\nACTC\n')
        assert guess_format(fhand) == 'fasta'

        # qual
        fhand = StringIO('>seq\n10 20\n')
        assert guess_format(fhand) == 'qual'

        # qual
        qual = ">seq1\n30 30 30 30 30 30 30 30\n>seq2\n30 30 30 30 30 30 30"
        qual += " 30\n>seq3\n30 30 30 30 30 30 30 30\n"

        fhand = StringIO(qual)
        assert guess_format(fhand) == 'qual'

    def test_unkown(self):
        'It tests unkown formats'
        fhand = StringIO('xseq\nACTC\n')
        try:
            guess_format(fhand)
            self.fail('UnknownFormatError expected')
        except UnknownFormatError:
            pass

    def test_empty_file(self):
        'It guesses the format of an empty file'
        fhand = StringIO()
        try:
            guess_format(fhand)
            self.fail('UnknownFormatError expected')
        except UnknownFormatError:
            pass

    def test_fastq(self):
        'It guesses the format for the solexa and illumina fastq'

        txt = '@HWI-EAS209_0006_FC706VJ:5:58:5894:21141#ATCACG/1\n'
        txt += 'TTAATTGGTAAATAAATCTCCTAATAGCTTAGATNTTACCTTNNNNNNNNNNTAGTTTCT\n'
        txt += '+HWI-EAS209_0006_FC706VJ:5:58:5894:21141#ATCACG/1\n'
        txt += 'efcfffffcfeefffcffffffddf`feed]`]_Ba_^__[YBBBBBBBBBBRTT\]][]\n'
        fhand = StringIO(txt)
        assert guess_format(fhand) == 'fastq-illumina'

        fhand = StringIO('@HWI-EAS209\n@')
        try:
            assert guess_format(fhand) == 'fasta'
            self.fail('UnknownFormatError expected')
        except UnknownFormatError:
            pass

        # sanger
        txt = '@HWI-EAS209_0006_FC706VJ:5:58:5894:21141#ATCACG/1\n'
        txt += 'TTAATTGGTAAATAAATCTCCTAATAGCTTAGATNTTACCTTNNNNNNNNNNTAGTTTCT\n'
        txt += '+HWI-EAS209_0006_FC706VJ:5:58:5894:21141#ATCACG/1\n'
        txt += '000000000000000000000000000000000000000000000000000000000000\n'
        fhand = StringIO(txt)
        assert guess_format(fhand) == 'fastq'

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

    def test_long_illumina(self):
        'The qualities seem illumina, but the reads are too lengthly'
        txt = '@read\n'
        txt += 'T' * 400 + '\n'
        txt += '+\n'
        txt += '@' * 400 + '\n'
        fhand = StringIO(txt)
        try:
            guess_format(fhand)
            self.fail('UndecidedFastqVersionError expected')
        except UndecidedFastqVersionError:
            pass

    def test_non_seekable(self):
        'Fastq version guessing using the non-seekable route'
        txt = '@HWI-EAS209_0006_FC706VJ:5:58:5894:21141#ATCACG/1\n'
        txt += 'TTAATTGGTAAATAAATCTCCTAATAGCTTAGATNTTACCTTNNNNNNNNNNTAGTTTCT\n'
        txt += '+HWI-EAS209_0006_FC706VJ:5:58:5894:21141#ATCACG/1\n'
        txt += 'efcfffffcfeefffcffffffddf`feed]`]_Ba_^__[YBBBBBBBBBBRTT\]][]\n'
        fhand = StringIO(txt)
        assert _guess_format(fhand, True) == 'fastq-illumina'

        fhand = StringIO('@HWI-EAS209\n@')
        try:
            assert _guess_format(fhand, True) == 'fasta'
            self.fail('UnknownFormatError expected')
        except UnknownFormatError:
            pass

        # sanger
        txt = '@HWI-EAS209_0006_FC706VJ:5:58:5894:21141#ATCACG/1\n'
        txt += 'TTAATTGGTAAATAAATCTCCTAATAGCTTAGATNTTACCTTNNNNNNNNNNTAGTTTCT\n'
        txt += '+HWI-EAS209_0006_FC706VJ:5:58:5894:21141#ATCACG/1\n'
        txt += '000000000000000000000000000000000000000000000000000000000000\n'
        fhand = StringIO(txt)
        assert _guess_format(fhand, True) == 'fastq'


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

if __name__ == "__main__":
#    import sys;sys.argv = ['', 'TestPool']
    unittest.main()
