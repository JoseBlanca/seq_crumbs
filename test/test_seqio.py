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

import os
import unittest
from  cStringIO import StringIO
from tempfile import NamedTemporaryFile
from subprocess import Popen, PIPE

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from crumbs.utils.test_utils import TEST_DATA_DIR
from crumbs.utils.bin_utils import BIN_DIR

from crumbs.seqio import (guess_seq_type, fastaqual_to_fasta, seqio,
                          write_seqrecords, read_seqrecords, _itemize_fasta,
                          _itemize_fastq, read_seqs, write_seqs)

from crumbs.utils.tags import SEQITEM, SEQRECORD
from crumbs.exceptions import IncompatibleFormatError, MalformedFile


FASTA = ">seq1\natctagtc\n>seq2\natctagtc\n>seq3\natctagtc\n"
QUAL = ">seq1\n30 30 30 30 30 30 30 30\n>seq2\n30 30 30 30 30 30 30 30\n"
QUAL += ">seq3\n30 30 30 30 30 30 30 30\n"
FASTQ = '@seq1\natcgt\n+\n?????\n@seq2\natcgt\n+\n?????\n@seq3\natcgt\n+\n'
FASTQ += '?????\n'


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

    def test_guess_seq_type(self):
        'It guesses if the sequence is nucleotide or protein'
        fpath = os.path.join(TEST_DATA_DIR, 'arabidopsis_genes')
        assert guess_seq_type(open(fpath)) == 'nucl'

        fpath = os.path.join(TEST_DATA_DIR, 'pairend2.sfastq')
        assert guess_seq_type(open(fpath)) == 'nucl'

    @staticmethod
    def test_fastaqual_to_fasta():
        seq_fhand = StringIO('>seq1\nattct\n>seq2\natc\n')
        qual_fhand = StringIO('>seq1\n2 2 2 2 2\n>seq2\n2 2 2\n')
        out_fhand = NamedTemporaryFile()
        fastaqual_to_fasta(seq_fhand, qual_fhand, out_fhand)
        fastq = open(out_fhand.name).read()
        assert fastq == "@seq1\nattct\n+\n#####\n@seq2\natc\n+\n###\n"

    def test_seqio(self):
        'It tets the seqio function'

        # fastq to fasta
        out_fhand = NamedTemporaryFile()
        seqio([self._make_fhand(FASTQ)], out_fhand, 'fasta')
        assert ">seq1\natcgt" in open(out_fhand.name).read()

        # fastq to fastq-illumina
        out_fhand = NamedTemporaryFile()
        seqio([self._make_fhand(FASTQ)], out_fhand, 'fastq-illumina')
        assert "@seq1\natcgt\n+\n^^^^" in open(out_fhand.name).read()

        out_fhand = NamedTemporaryFile()
        seqio([self._make_fhand(FASTQ), self._make_fhand(FASTQ)],
              out_fhand, 'fastq-illumina')

        assert "@seq3\natcgt\n+\n^^^^^\n@seq1" in open(out_fhand.name).read()

        # fasta to fastq
        out_fhand = NamedTemporaryFile()
        try:
            seqio([self._make_fhand(FASTA)], out_fhand, 'fastq')
            self.fail("error previously expected")
        except IncompatibleFormatError as error:
            assert 'No qualities available' in str(error)

        # bad_format fastq
        bad_fastq_fhand = self._make_fhand(FASTQ + 'aklsjhdas')
        try:
            seqio([bad_fastq_fhand], out_fhand, 'fasta')
            self.fail("error previously expected")
        except MalformedFile as error:
            assert 'Lengths of sequence and quality'  in str(error)

        # genbank to fasta
        out_fhand = NamedTemporaryFile()
        genbank_fhand = open(os.path.join(TEST_DATA_DIR, 'sequence.gb'))
        seqio([genbank_fhand], out_fhand, 'fasta')
        result = open(out_fhand.name).read()
        assert '>NM_019354.2' in result


class ReadWriteSeqsTest(unittest.TestCase):
    'It writes seqrecords in a file'
    def test_write_empy_seq(self):
        'It does not write an empty sequence'
        seq1 = SeqRecord(Seq('ACTG'), id='seq1')
        fhand = StringIO()
        write_seqrecords([seq1, None, SeqRecord(Seq(''), id='seq2')], fhand,
                         file_format='fasta')
        fhand.flush()
        assert fhand.getvalue() == '>seq1\nACTG\n'

    def test_read_fasta(self):
        'It tests the reading of a fasta file'
        fhand = StringIO('>seq1\nACTG\n')
        assert not list(read_seqrecords([fhand]))[0].description


class SimpleIOTest(unittest.TestCase):
    'It tests the simple input and output read'
    def test_fasta_itemizer(self):
        'It tests the fasta itemizer'
        fhand = StringIO('>s1\nACTG\n>s2 desc\nACTG\n')
        seqs = list(_itemize_fasta(fhand))
        assert seqs == [('s1', ['>s1\n', 'ACTG\n'], {}),
                        ('s2', ['>s2 desc\n', 'ACTG\n'], {})]

        # with several lines
        fhand = StringIO('>s1\nACTG\nGTAC\n>s2 desc\nACTG\n')
        seqs = list(_itemize_fasta(fhand))
        assert seqs == [('s1', ['>s1\n', 'ACTG\n', 'GTAC\n'], {}),
                        ('s2', ['>s2 desc\n', 'ACTG\n'], {})]

    def test_fastq_itemizer(self):
        'It tests the fasta itemizer'
        fhand = StringIO('@s1\nACTG\n+\n1234\n@s2 desc\nACTG\n+\n4321\n')
        seqs = list(_itemize_fastq(fhand))
        assert seqs == [('s1', ['@s1\n', 'ACTG\n', '+\n', '1234\n'], {}),
                        ('s2', ['@s2 desc\n', 'ACTG\n', '+\n', '4321\n'], {})]

    def test_seqitems_io(self):
        'It checks the different seq class streams IO'
        fhand = StringIO('>s1\nACTG\n>s2 desc\nACTG\n')
        seqs = list(read_seqs([fhand], 'fasta',
                              prefered_seq_classes=[SEQITEM]))
        assert seqs[0].kind == SEQITEM
        fhand = StringIO()
        write_seqs(seqs, fhand)
        assert fhand.getvalue() == '>s1\nACTG\n>s2 desc\nACTG\n'
        assert seqs[0].object.name == 's1'

        # SeqRecord
        fhand = StringIO('>s1\nACTG\n>s2 desc\nACTG\n')
        seqs = list(read_seqs([fhand], 'fasta',
                              prefered_seq_classes=[SEQRECORD]))
        assert seqs[0].kind == SEQRECORD
        fhand = StringIO()
        write_seqs(seqs, fhand, 'fasta')
        assert fhand.getvalue() == '>s1\nACTG\n>s2 desc\nACTG\n'

        # seqitem not possible with different input and output formats
        fhand = StringIO('>s1\nACTG\n>s2 desc\nACTG\n')
        try:
            seqs = list(read_seqs([fhand], 'fasta', out_format='fastq',
                        prefered_seq_classes=[SEQITEM]))
            self.fail('ValueError expected')
        except ValueError:
            pass

        fhand = StringIO('>s1\nACTG\n>s2 desc\nACTG\n')
        seqs = list(read_seqs([fhand], 'fasta', out_format='fasta',
                        prefered_seq_classes=[SEQITEM]))
        fhand = StringIO()
        write_seqs(seqs, fhand)
        assert fhand.getvalue() == '>s1\nACTG\n>s2 desc\nACTG\n'


class PipingTest(unittest.TestCase):
    'It tests that we get no error when trying to write in a closed pipe'
    def test_write_closed_pipe(self):
        seq_fhand = NamedTemporaryFile(suffix='.fasta')
        n_seqs = 1000
        for i in range(n_seqs):
            seq_fhand.write('>s\nACTG\n')
        seq_fhand.flush()
        in_fpath = seq_fhand.name
        seq_head = os.path.join(BIN_DIR, 'seq_head')

        process_seq = Popen([seq_head, '-n', str(n_seqs), in_fpath],
                            stdout=PIPE)
        stdout = NamedTemporaryFile(suffix='.stdout')
        process_head = Popen(['head', '-n', '1'], stdin=process_seq.stdout,
                             stdout=stdout)
        process_seq.stdout.close()  # Allow seq_head to receive a SIGPIPE if
                                    # head exits.
        process_head.communicate()

        assert open(stdout.name).read() == '>s\n'
        seq_fhand.close()
        stdout.close()

        # With SeqRecords
        gb_fpath = os.path.join(TEST_DATA_DIR, 'sequence.gb')
        gb_content = open(gb_fpath).read()
        seq_fhand = NamedTemporaryFile(suffix='.gb')
        n_seqs = 100
        for i in range(n_seqs):
            seq_fhand.write(gb_content)
        seq_fhand.flush()
        in_fpath = seq_fhand.name

        process_seq = Popen([seq_head, '-n', str(n_seqs), in_fpath],
                            stdout=PIPE)
        stdout = NamedTemporaryFile(suffix='.stdout')
        process_head = Popen(['head', '-n', '1'], stdin=process_seq.stdout,
                             stdout=stdout)
        process_seq.stdout.close()  # Allow seq_head to receive a SIGPIPE if
                                    # head exits.
        process_head.communicate()

        seq_fhand.close()
        assert 'LOCUS' in open(stdout.name).read()
        stdout.close()

if __name__ == '__main__':
    # import sys;sys.argv = ['', 'SffExtractTest.test_items_in_gff']
    unittest.main()
