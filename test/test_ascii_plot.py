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
from crumbs.statistics import IntsStats


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


class IntsStatsTest(unittest.TestCase):
    'It test the extensible array class'
    @staticmethod
    def create_test_array():
        d = {'9': '5', '10': '288', '11': '002556688', '12': '00012355555',
             '13': '0000013555688', '14': '00002555558',
             '15': '0000000000355555555557', '16': '000045', '17': '000055',
             '18': '0005', '19': '00005', '21': '5'}
        ext_array = IntsStats()
        for key, values in d.items():
            for num in values:
                ext_array.append(int(key + num))
        return ext_array

    @staticmethod
    def test_array():
        'Create an extensible array'
        ext_array = IntsStats(init_len=5)
        ext_array.append(6)
        ext_array.append(2)
        assert  ext_array.min == 2
        assert  ext_array.max == 6
        ext_array.append(200)
        assert ext_array.max == 200

        assert list(ext_array.flat) == [2, 6, 200]

        input_ = (3, 5, 7, 7, 38)
        ext_array = IntsStats(input_)
        assert ext_array.median == 7
        assert list(ext_array.flat) == [3, 5, 7, 7, 38]

    def test_distribution(self):
        'It tests the histogram function'

        ints_array = self.create_test_array()
        distrib = ints_array.calculate_distribution(bins=10,
                                                        remove_outliers=5)

        assert distrib['distrib'] == [7L, 13L, 7L, 10L, 7L, 22L, 6L, 4L, 5L,
                                      5L]
        assert distrib['bin_edges'] == [110, 118, 126, 134, 142, 150, 158, 166,
                                        174, 182, 190]

        assert 'average' in str(ints_array)

        ints_array = IntsStats([0, 0, 1, 3])
        assert ints_array.calculate_distribution(bins=3)['distrib'] == [2, 1,
                                                                        1]

    def test_stats_functs(self):
        'It test the statistical functions of the class'
        ext_array = IntsStats()
        ext_array.append(3)
        ext_array.append(5)
        ext_array.append(7)
        ext_array.append(7)
        ext_array.append(38)
        assert ext_array.median == 7

        # median with two middle numbers
        ext_array = IntsStats([3, 5, 7, 7])
        assert ext_array.median == 6

        ints = IntsStats([34, 43, 81, 106, 106, 115])
        assert ints.average - 80.83 < 0.01

        ext_array = self.create_test_array()
        assert ext_array.median == 145
        assert round(ext_array.average, 2) == 145.15

        assert ext_array.sum == 13354
        assert ext_array.count == 92
        assert round(ext_array.variance, 2) == 557.43

if __name__ == '__main__':
    #import sys;sys.argv = ['', 'SffExtractTest.test_items_in_gff']
    unittest.main()
