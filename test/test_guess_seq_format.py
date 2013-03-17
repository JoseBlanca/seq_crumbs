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

import unittest
import os.path
from subprocess import check_output, CalledProcessError
from tempfile import NamedTemporaryFile
from StringIO import StringIO

from crumbs.utils.bin_utils import BIN_DIR
from crumbs.utils.seq_utils import guess_format, _guess_format
from crumbs.exceptions import UnknownFormatError, UndecidedFastqVersionError

# pylint: disable=R0201
# pylint: disable=R0904


class GuessFormatBinTest(unittest.TestCase):
    'It tests the guess_seq_format binary'

    def test_guess_format(self):
        'It tests guess_seq_format'

        guess_bin = os.path.join(BIN_DIR, 'guess_seq_format')

        assert 'usage' in check_output([guess_bin, '-h'])

        # a fasta file
        fasta_fhand = NamedTemporaryFile()
        fasta_fhand.write('>seq\nACTA\n')
        fasta_fhand.flush()
        assert check_output([guess_bin, fasta_fhand.name]) == 'fasta\n'

        # Unknown_format
        bad_fhand = NamedTemporaryFile()
        bad_fhand.write('bad file')
        bad_fhand.flush()
        stderr = NamedTemporaryFile()
        try:
            check_output([guess_bin, bad_fhand.name], stderr=stderr)
            self.fail('Error expected')
        except CalledProcessError:
            assert open(stderr.name).read().startswith('Sequence file of unkn')

    def test_stdin(self):
        'It works with stdin'
        guess_bin = os.path.join(BIN_DIR, 'guess_seq_format')
        fasta_fhand = NamedTemporaryFile()
        fasta_fhand.write('>seq\nACTA\n')
        fasta_fhand.flush()
        fmt = check_output([guess_bin], stdin=open(fasta_fhand.name))
        assert fmt == 'fasta\n'

    def test_version(self):
        'It can return its version number'
        guess_bin = os.path.join(BIN_DIR, 'guess_seq_format')

        stderr = NamedTemporaryFile()
        check_output([guess_bin, '--version'], stderr=stderr)
        assert 'from seq_crumbs version:' in open(stderr.name).read()


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

        fhand = StringIO(txt + txt)
        assert guess_format(fhand) == 'fastq-illumina-one_line'

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

if __name__ == '__main__':
    #import sys;sys.argv = ['', 'SffExtractTest.test_items_in_gff']
    unittest.main()
