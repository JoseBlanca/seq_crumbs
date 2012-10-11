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

import os
import unittest
from crumbs.utils.test_utils import TEST_DATA_DIR
from crumbs.seqio import guess_seq_type


class SeqioTest(unittest.TestCase):
    'It test seqIO functions'

    @staticmethod
    def test_guess_seq_type():
        fpath = os.path.join(TEST_DATA_DIR, 'arabidopsis_genes')
        assert guess_seq_type(open(fpath)) == 'nucl'

        fpath = os.path.join(TEST_DATA_DIR, 'pairend2.sfastq')
        assert guess_seq_type(open(fpath)) == 'nucl'

if __name__ == '__main__':
    #import sys;sys.argv = ['', 'SffExtractTest.test_items_in_gff']
    unittest.main()
