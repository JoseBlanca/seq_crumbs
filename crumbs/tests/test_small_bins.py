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
import os
from subprocess import check_output, CalledProcessError
from tempfile import NamedTemporaryFile
from StringIO import StringIO

from crumbs.tests.utils import BIN_DIR
from crumbs.seqio import count_seqs_in_files

# pylint: disable=R0201
# pylint: disable=R0904


class CatHeadTest(unittest.TestCase):
    'It tests the cat seqsbinary'
    counter = 1

    def make_fasta(self):
        'it returns a fasta fhand'
        fhand = NamedTemporaryFile()
        fhand.write('>seq{0:d}\nACTATCATGGCAGATA\n'.format(self.counter))
        fhand.flush()
        self.counter += 1
        return fhand

    def test_cat_seqs(self):
        'It test the cat seqs'
        cat_bin = os.path.join(BIN_DIR, 'cat_seqs')
        assert check_output([cat_bin]).startswith('usage')

        #fasta to fasta
        in_fhand1 = self.make_fasta()
        in_fhand2 = self.make_fasta()
        result = check_output([cat_bin, '-f', 'fasta', in_fhand1.name,
                               in_fhand2.name])
        assert '>seq1\nACTATCATGGCAGATA\n>seq2\nACTATCATGGCAGATA' in result

        #from fasta to fastq
        try:
            stderr = NamedTemporaryFile()
            result = check_output([cat_bin, '-f', 'fastq', in_fhand1.name,
                               in_fhand2.name], stderr=stderr)
            self.fail()
        except CalledProcessError:
            assert 'output format is incompatible'  in open(stderr.name).read()


class SeqHeadTest(unittest.TestCase):
    'It tests the seq head binary'

    def test_seq_head(self):
        'It tests the seq head'
        head_bin = os.path.join(BIN_DIR, 'seq_head')
        assert check_output([head_bin]).startswith('usage')

        #get one seq
        fasta_fhand = NamedTemporaryFile()
        fasta_fhand.write('>seq\nACTA\n>seq2\nACTA\n>seq3\nACTA\n')
        fasta_fhand.flush()

        result = check_output([head_bin, '-n', '1', fasta_fhand.name])
        assert result == '>seq\nACTA\n'


class SampleSeqTest(unittest.TestCase):
    'It tests the seq head binary'

    def test_sample_seq(self):
        'It tests the seq head'
        sample_seq = os.path.join(BIN_DIR, 'sample_seqs')
        assert check_output([sample_seq]).startswith('usage')

        fasta_fhand = NamedTemporaryFile()
        fasta_fhand.write('>seq\nACTA\n>seq2\nACTA\n>seq3\nACTA\n')
        fasta_fhand.flush()

        #random sample
        result = check_output([sample_seq, '-n', '1', fasta_fhand.name])
        assert count_seqs_in_files([StringIO(result)]) == 1

        #random sample
        result = check_output([sample_seq, '-n', '2', fasta_fhand.name])
        assert count_seqs_in_files([StringIO(result)]) == 2

if __name__ == '__main__':
    #import sys;sys.argv = ['', 'SffExtractTest.test_items_in_gff']
    unittest.main()
