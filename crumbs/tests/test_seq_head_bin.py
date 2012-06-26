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
Created on 2012 eka 21

@author: peio
'''
import unittest
import os
from tempfile import NamedTemporaryFile
from subprocess import check_output

from crumbs.tests.utils import BIN_DIR

# pylint: disable=R0201
# pylint: disable=R0904


class SeqHead(unittest.TestCase):
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

if __name__ == '__main__':
    #import sys;sys.argv = ['', 'SffExtractTest.test_items_in_gff']
    unittest.main()
