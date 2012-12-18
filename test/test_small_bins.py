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
from gzip import GzipFile

from crumbs.utils.bin_utils import BIN_DIR
from crumbs.seqio import count_seqs_in_files

# pylint: disable=R0201
# pylint: disable=R0904


class CatTest(unittest.TestCase):
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

        # help
        assert 'usage' in check_output([cat_bin, '-h'])

        # fasta to fasta
        in_fhand1 = self.make_fasta()
        in_fhand2 = self.make_fasta()
        result = check_output([cat_bin, in_fhand1.name, in_fhand2.name])
        assert '>seq1\nACTATCATGGCAGATA\n>seq2\nACTATCATGGCAGATA' in result

        # from fastq to fastq
        fhand = NamedTemporaryFile()
        fhand.write('@seq1\nACTA\n+\nqqqq\n')
        fhand.flush()
        result = check_output([cat_bin, fhand.name])
        assert result == '@seq1\nACTA\n+\nqqqq\n'

        # No input
        fhand = NamedTemporaryFile()
        fhand.write('')
        fhand.flush()
        try:
            stderr = NamedTemporaryFile()
            result = check_output([cat_bin, fhand.name], stderr=stderr)
            self.fail()
        except CalledProcessError:
            assert 'The file is empty'  in open(stderr.name).read()

        # No format
        in_fhand1 = self.make_fasta()
        in_fhand2 = self.make_fasta()
        result = check_output([cat_bin, in_fhand1.name, in_fhand2.name])
        assert '>seq3\nACTATCATGGCAGATA\n>seq4\nACTATCATGGCAGATA' in result

        # bad input format
        in_fhand1 = self.make_fasta()
        try:
            stderr = NamedTemporaryFile()
            result = check_output([cat_bin, '-t', 'fastq', in_fhand1.name],
                                  stderr=stderr)
            self.fail()
        except CalledProcessError:
            stderr_str = open(stderr.name).read()
            assert 'output format is incompatible with input' in stderr_str

    def test_gziped_output(self):
        'It writes a gziped file'
        in_fhand = self.make_fasta()
        cat_bin = os.path.join(BIN_DIR, 'cat_seqs')

        # bgzf does not work with STDOUT
        try:
            stderr = NamedTemporaryFile()
            check_output([cat_bin, '-Z', in_fhand.name], stderr=stderr)
            self.fail('CalledProcessError expected')
        except CalledProcessError:
            msg = 'bgzf is only available to seekable files'
            assert msg in open(stderr.name).read()

        # bgzf
        out_fhand = NamedTemporaryFile()
        check_output([cat_bin, '-Z', '-o', out_fhand.name, in_fhand.name])
        result = GzipFile(out_fhand.name).read()
        assert '\nACTATCATGGCAGATA\n' in  result

        # gzip
        result = check_output([cat_bin, '-z', in_fhand.name])
        result = GzipFile(fileobj=StringIO(result)).read()
        assert '\nACTATCATGGCAGATA\n' in  result

    def test_compressed_input(self):
        'It can read compressed files'
        # we need a compressed file
        in_fhand = self.make_fasta()
        cat_bin = os.path.join(BIN_DIR, 'cat_seqs')
        out_fhand = NamedTemporaryFile()
        check_output([cat_bin, '-Z', '-o', out_fhand.name, in_fhand.name])
        result = GzipFile(out_fhand.name).read()
        assert '\nACTATCATGGCAGATA\n' in  result

        result = check_output([cat_bin, out_fhand.name])
        assert '>seq1\nACTATCATGGCAGATA\n' in  result

    def test_version(self):
        'It can return its version number'
        guess_bin = os.path.join(BIN_DIR, 'cat_seqs')

        stderr = NamedTemporaryFile()
        check_output([guess_bin, '--version'], stderr=stderr)
        assert 'from seq_crumbs version:' in open(stderr.name).read()


class SeqHeadTest(unittest.TestCase):
    'It tests the seq head binary'

    def test_seq_head(self):
        'It tests the seq head'
        head_bin = os.path.join(BIN_DIR, 'seq_head')
        #assert check_output([head_bin, '-h']).startswith('usage')

        #get one seq
        fasta_fhand = NamedTemporaryFile()
        fasta_fhand.write('>seq\nACTA\n>seq2\nACTA\n>seq3\nACTA\n')
        fasta_fhand.flush()

        result = check_output([head_bin, '-n', '1', fasta_fhand.name])
        assert result == '>seq\nACTA\n'


class SampleSeqTest(unittest.TestCase):
    'It tests the seq sample binary'

    def test_sample_seq(self):
        'It tests the seq head'
        sample_seq = os.path.join(BIN_DIR, 'sample_seqs')
        assert 'usage' in check_output([sample_seq, '-h'])

        fasta_fhand = NamedTemporaryFile()
        fasta_fhand.write('>seq\nACTA\n>seq2\nACTA\n>seq3\nACTA\n')
        fasta_fhand.flush()

        #random sample
        result = check_output([sample_seq, '-n', '1', fasta_fhand.name])
        assert count_seqs_in_files([StringIO(result)], ['fasta']) == 1

        #random sample
        result = check_output([sample_seq, '-n', '2', fasta_fhand.name])
        assert count_seqs_in_files([StringIO(result)], ['fasta']) == 2

        #random sample with stdin
        result = check_output([sample_seq, '-n', '2'],
                              stdin=open(fasta_fhand.name))
        assert count_seqs_in_files([StringIO(result)], ['fasta']) == 2

if __name__ == '__main__':
#    import sys;sys.argv = ['', 'SampleSeqTest']
    unittest.main()
