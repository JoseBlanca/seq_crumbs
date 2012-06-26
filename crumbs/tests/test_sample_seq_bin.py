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

'''
Created on 2012 eka 22

@author: peio
'''
import unittest
import os
from StringIO import StringIO
from tempfile import NamedTemporaryFile
from subprocess import check_output, CalledProcessError, check_call

from crumbs.tests.utils import BIN_DIR
from crumbs.seqio import count_seqs_in_files


# pylint: disable=R0201
# pylint: disable=R0904


class SampleSeq(unittest.TestCase):
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
