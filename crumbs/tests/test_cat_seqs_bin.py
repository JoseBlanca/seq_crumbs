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

from crumbs.tests.utils import BIN_DIR
from crumbs.exceptions import IncompatibleFormatError


class CatHead(unittest.TestCase):
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


if __name__ == '__main__':
    #import sys;sys.argv = ['', 'SffExtractTest.test_items_in_gff']
    unittest.main()
