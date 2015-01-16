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
from crumbs.seq.utils.file_formats import get_format, _guess_format
from crumbs.exceptions import (UnknownFormatError, FileIsEmptyError,
                               UndecidedFastqVersionError)

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
        assert get_format(fhand) == 'fasta'

        # multiline fasta
        fhand = StringIO('>seq\nACTC\nACTG\n>seq2\nACTG\n')
        assert get_format(fhand) == 'fasta'

        # qual
        fhand = StringIO('>seq\n10 20\n')
        assert get_format(fhand) == 'qual'

        # qual
        qual = ">seq1\n30 30 30 30 30 30 30 30\n>seq2\n30 30 30 30 30 30 30"
        qual += " 30\n>seq3\n30 30 30 30 30 30 30 30\n"

        fhand = StringIO(qual)
        assert get_format(fhand) == 'qual'

    def test_with_long_desc(self):
        fhand = StringIO('''>comp27222_c1_seq1 len=4926 path=[89166356:0-46 89167522:47-85 89315292:86-121 89170132:122-176 89377211:177-217 89377235:218-244 89172846:245-247 89172856:248-251 89173028:252-276 89174386:277-292 89174684:293-506 89377352:507-582 89183669:583-587 89183821:588-613 89184868:614-644 89185624:645-719 89187914:720-723 89187935:724-870 89191280:871-887 89377494:888-907 89191517:908-927 89193046:928-1071 89198507:1072-1109 89199632:1110-1170 89201544:1171-1194 89202607:1195-1247 89377606:1248-1252 89377611:1253-1591 89215759:1592-1606 89215815:1607-1636 89216359:1637-1664 89377693:1665-1678 88727916:1679-2152 88743802:2153-2171 88744738:2172-2623 88759485:2624-2648 88759762:2649-2953 88769199:2954-2971 88769596:2972-3657 88791809:3658-3665 88792014:3666-3723 88793720:3724-3731 88794381:3732-3812 88799277:3813-3813 88799328:3814-3996 88807093:3997-3999 88807177:4000-4215 88813164:4216-4246 88814188:4247-4287 88815355:4288-4308 88816198:4309-4352 88817845:4353-4369 88818294:4370-4403 88818879:4404-4465 88821150:4466-4469 88821188:4470-4925]
GAAGGATCGATCGGCCTCGGCGGTGTTCCCAAAAATCTAAGAGCGTTTACTCCAAGCTTC''')
        get_format(fhand)

    def test_unkown(self):
        'It tests unkown formats'
        fhand = StringIO('xseq\nACTC\n')
        try:
            get_format(fhand)
            self.fail('UnknownFormatError expected')
        except UnknownFormatError:
            pass

    def test_empty_file(self):
        'It guesses the format of an empty file'
        fhand = StringIO()
        try:
            get_format(fhand)
            self.fail('FileIsEmptyError expected')
        except FileIsEmptyError:
            pass

    def test_fastq(self):
        'It guesses the format for the solexa and illumina fastq'

        txt = '@HWI-EAS209_0006_FC706VJ:5:58:5894:21141#ATCACG/1\n'
        txt += 'TTAATTGGTAAATAAATCTCCTAATAGCTTAGATNTTACCTTNNNNNNNNNNTAGTTTCT\n'
        txt += '+HWI-EAS209_0006_FC706VJ:5:58:5894:21141#ATCACG/1\n'
        txt += 'efcfffffcfeefffcffffffddf`feed]`]_Ba_^__[YBBBBBBBBBBRTT\]][]\n'
        fhand = StringIO(txt)
        assert get_format(fhand) == 'fastq-illumina'

        txt = '@HWI-EAS209_0006_FC706VJ:5:58:5894:21141#ATCACG/1\n'
        txt += 'TTAATTGGTAAATAAATCTCCTAATAGCTTAGATNTTACCTTNNNNNNNNNNTAGTTTCT\n'
        txt += 'TTAATTGGTAAATAAATCTCCTAATAGCTTAGATNTTACCTTNNNNNNNNNNTAGTTTCT\n'
        txt += '+HWI-EAS209_0006_FC706VJ:5:58:5894:21141#ATCACG/1\n'
        txt += 'efcfffffcfeefffcffffffddf`feed]`]_Ba_^__[YBBBBBBBBBBRTT\]][]\n'
        txt += 'efcfffffcfeefffcffffffddf`feed]`]_Ba_^__[YBBBBBBBBBBRTT\]][]\n'

        fhand = StringIO(txt + txt)
        assert get_format(fhand) == 'fastq-illumina'

        fhand = StringIO('@HWI-EAS209\n@')
        try:
            assert get_format(fhand) == 'fasta'
            self.fail('UndecidedFastqVersionError expected')
        except UndecidedFastqVersionError:
            pass

        # sanger
        txt = '@HWI-EAS209_0006_FC706VJ:5:58:5894:21141#ATCACG/1\n'
        txt += 'TTAATTGGTAAATAAATCTCCTAATAGCTTAGATNTTACCTTNNNNNNNNNNTAGTTTCT\n'
        txt += '+HWI-EAS209_0006_FC706VJ:5:58:5894:21141#ATCACG/1\n'
        txt += '000000000000000000000000000000000000000000000000000000000000\n'
        fhand = StringIO(txt)
        assert get_format(fhand) == 'fastq'

    def test_long_illumina(self):
        'The qualities seem illumina, but the reads are too lengthly'
        txt = '@read\n'
        txt += 'T' * 400 + '\n'
        txt += '+\n'
        txt += '@' * 400 + '\n'
        fhand = StringIO(txt)
        try:
            get_format(fhand)
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
            self.fail('UndecidedFastqVersionError expected')
        except UndecidedFastqVersionError:
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
