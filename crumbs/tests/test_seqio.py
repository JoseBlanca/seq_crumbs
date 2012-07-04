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
from tempfile import NamedTemporaryFile
from StringIO import StringIO
from subprocess import check_output, CalledProcessError

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from crumbs.seqio import guess_format, seqio, write_seqrecords, read_seqrecords
from crumbs.exceptions import UnknownFormatError
from crumbs.tests.utils import BIN_DIR

# pylint: disable=R0201
# pylint: disable=R0904

FASTA = ">seq1\natctagtc\n>seq2\natctagtc\n>seq3\natctagtc\n"
QUAL = ">seq1\n30 30 30 30 30 30 30 30\n>seq2\n30 30 30 30 30 30 30 30\n"
QUAL += ">seq3\n30 30 30 30 30 30 30 30\n"
FASTQ = '@seq1\natcgt\n+\n?????\n@seq2\natcgt\n+\n?????\n@seq3\natcgt\n+\n?????\n'


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
        fhand = StringIO(QUAL)
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


class SeqIOTest(unittest.TestCase):
    'It tests the seqio functions'

    @staticmethod
    def _make_fhand(content=None):
        'It makes temporary fhands'
        if content is None:
            content = ''
        fhand = NamedTemporaryFile()
        fhand.write(content)
        fhand.flush()
        return fhand

    def test_seqio(self):
        'It tets the seqio function'

        #fasta-qual to fastq
        in_fhands = (self._make_fhand(FASTA), self._make_fhand(QUAL))
        out_fhands = (self._make_fhand(),)
        out_format = 'fastq'
        seqio(in_fhands, out_fhands, out_format)
        assert "@seq1\natctagtc\n+\n???????" in open(out_fhands[0].name).read()

        #fastq to fasta-qual
        out_fhands = [self._make_fhand(), self._make_fhand()]
        seqio([self._make_fhand(FASTQ)], out_fhands, 'fasta')
        assert ">seq1\natcgt" in open(out_fhands[0].name).read()
        assert ">seq1\n30 30 30" in open(out_fhands[1].name).read()

        #fastq to fasta
        out_fhands = [self._make_fhand()]
        seqio([self._make_fhand(FASTQ)], out_fhands, 'fasta')
        assert ">seq1\natcgt" in open(out_fhands[0].name).read()

        #fastq to fastq-illumina
        out_fhands = [self._make_fhand()]
        seqio([self._make_fhand(FASTQ)], out_fhands, 'fastq-illumina')
        assert "@seq1\natcgt\n+\n^^^^" in open(out_fhands[0].name).read()

        #fasta-qual to fasta-qual
        in_fhands = (self._make_fhand(FASTA), self._make_fhand(QUAL))
        out_fhands = (self._make_fhand(), self._make_fhand())
        out_format = 'fasta'
        seqio(in_fhands, out_fhands, out_format)
        assert FASTA == open(out_fhands[0].name).read()
        assert QUAL == open(out_fhands[1].name).read()


class ReadWriteSeqsTest(unittest.TestCase):
    'It writes seqrecords in a file'
    def test_write_empy_seq(self):
        'It does not write an empty sequence'
        seq1 = SeqRecord(Seq('ACTG'), id='seq1')
        fhand = StringIO()
        write_seqrecords(fhand, [seq1, None, SeqRecord(Seq(''), id='seq2')],
                         file_format='fasta')
        assert fhand.getvalue() == '>seq1\nACTG\n'

    def test_read_fasta(self):
        'It tests the reading of a fasta file'
        fhand = StringIO('>seq1\nACTG\n')
        assert not list(read_seqrecords([fhand]))[0].description


class SeqioBinTest(unittest.TestCase):
    'It test the seqio binary'

    @staticmethod
    def _make_fhand(content=None):
        'It makes temporary fhands'
        if content is None:
            content = ''
        fhand = NamedTemporaryFile()
        fhand.write(content)
        fhand.flush()
        return fhand

    def test_seqio_bin(self):
        'It test the seqio binary'
        seqio_bin = os.path.join(BIN_DIR, 'seqio')
        assert check_output([seqio_bin]).startswith('usage')

        #get one se
        fasta_fhand = self._make_fhand(FASTA)
        qual_fhand = self._make_fhand(QUAL)
        fastq_fhand = self._make_fhand(FASTQ)

        # fasta-qual to fastq
        out_fhand = NamedTemporaryFile()
        check_output([seqio_bin, '-o', out_fhand.name, '-f', 'fastq',
                      fasta_fhand.name, qual_fhand.name])
        assert "@seq1\natctagtc\n+" in  open(out_fhand.name).read()

        #fastq to fast-qual
        fasta_out_fhand = NamedTemporaryFile()
        qual_out_fhand = NamedTemporaryFile()
        check_output([seqio_bin, '-o', fasta_out_fhand.name,
                      qual_out_fhand.name, '-f', 'fasta', fastq_fhand.name])
        assert ">seq1\natcgt" in  open(fasta_out_fhand.name).read()
        assert ">seq1\n30 30 30 30 30" in  open(qual_out_fhand.name).read()

        #bad_format_fasta
        bad_fasta_fhand = self._make_fhand(FASTA + 'asdsa')
        out_fhand = NamedTemporaryFile()
        stderr = NamedTemporaryFile()
        try:
            check_output([seqio_bin, '-o', out_fhand.name, '-f', 'fastq',
                         bad_fasta_fhand.name, qual_fhand.name], stderr=stderr)
            self.fail('error expected')
        except CalledProcessError:
            assert 'Sequence length and number' in open(stderr.name).read()

        #bad_format_fastq
        bad_fastq_fhand = self._make_fhand(FASTQ + 'aklsjhdas')
        fasta_out_fhand = NamedTemporaryFile()
        qual_out_fhand = NamedTemporaryFile()
        stderr = NamedTemporaryFile()
        try:
            check_output([seqio_bin, '-o', fasta_out_fhand.name,
                          qual_out_fhand.name, '-f', 'fasta',
                           bad_fastq_fhand.name], stderr=stderr)
            self.fail('error expected')
        except CalledProcessError:
            assert 'Lengths of sequence and qualit' in open(stderr.name).read()

        #malformed fastq to fastq
        bad_fastq_fhand = self._make_fhand(FASTQ + 'aklsjhdas')
        fastq_out_fhand = NamedTemporaryFile()
        stderr = NamedTemporaryFile()
        try:
            check_output([seqio_bin, '-o', fastq_out_fhand.name, '-f',
                          'fastq-illumina', bad_fastq_fhand.name],
                          stderr=stderr)
            self.fail('error expected')
        except CalledProcessError:
            assert 'Lengths of sequence and qualit' in open(stderr.name).read()

if __name__ == '__main__':
    #import sys;sys.argv = ['', 'class.fuanction']
    unittest.main()
