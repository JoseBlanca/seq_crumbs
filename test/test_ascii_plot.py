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

import unittest

from crumbs.ascii_plot import draw_histogram


class HistogramTest(unittest.TestCase):
    'It tests the ASCII histogram plotting'
    def test_ascii_histogram(self):
        'It plots an ASCII histogram'
        hist = draw_histogram(bin_limits=[-2, -1, 0, 1, 2],
                              counts=[9, 20, 30, 40])
        expected = '[-2 , -1[ ( 9): ******************\n'
        expected += '[-1 ,  0[ (20): ****************************************'
        expected += '\n'
        expected += '[ 0 ,  1[ (30): *****************************************'
        expected += '*******************\n'
        expected += '[ 1 ,  2[ (40): *****************************************'
        expected += '***************************************\n'
        assert hist == expected

if __name__ == '__main__':
    #import sys;sys.argv = ['', 'SffExtractTest.test_items_in_gff']
    unittest.main()
