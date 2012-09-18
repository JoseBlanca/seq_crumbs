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

from crumbs.statistics import IntSumarizedArray, draw_histogram, IntBoxplot


class HistogramTest(unittest.TestCase):
    'It tests the ASCII histogram plotting'
    def test_ascii_histogram(self):
        'It plots an ASCII histogram'
        hist = draw_histogram(bin_limits=[-2, -1, 0, 1, 2],
                              counts=[9, 20, 30, 40])
        assert '[-2 , -1[ ( 9): ****************' in hist


class IntsStatsTest(unittest.TestCase):
    'It test the extensible array class'
    # pylint: disable=C0111
    # pylint: disable=C0103
    @staticmethod
    def create_test_array():
        d = {'9': '5', '10': '288', '11': '002556688', '12': '00012355555',
             '13': '0000013555688', '14': '00002555558',
             '15': '0000000000355555555557', '16': '000045', '17': '000055',
             '18': '0005', '19': '00005', '21': '5'}
        ext_array = IntSumarizedArray()
        for key, values in d.items():
            for num in values:
                ext_array.append(int(key + num))
        return ext_array

    @staticmethod
    def test_array():
        'Create an extensible array'
        ext_array = IntSumarizedArray(init_len=5)
        ext_array.append(6)
        ext_array.append(2)
        assert  ext_array.min == 2
        assert  ext_array.max == 6
        ext_array.append(200)
        assert ext_array.max == 200

        assert list(ext_array.flat) == [2, 6, 200]

        input_ = (3, 5, 7, 7, 38)
        ext_array = IntSumarizedArray(input_)
        assert ext_array.median == 7
        assert list(ext_array.flat) == [3, 5, 7, 7, 38]

    @staticmethod
    def test__add__():
        ext_array = IntSumarizedArray(init_len=5)
        ext_array.append(6)
        ext_array.append(2)

        ext_array2 = IntSumarizedArray(init_len=7)
        ext_array2.append(7)
        ext_array2.append(2)

        new_array = ext_array + ext_array2
        assert list(new_array.flat) == [2, 2, 6, 7]

    def test_distribution(self):
        'It tests the histogram function'

        ints_array = self.create_test_array()
        distrib = ints_array.calculate_distribution(bins=10,
                                                        remove_outliers=5)

        assert distrib['counts'] == [7L, 13L, 7L, 10L, 7L, 22L, 6L, 4L, 5L,
                                      5L]
        assert distrib['bin_limits'] == [110, 118, 126, 134, 142, 150, 158,
                                         166, 174, 182, 190]

        assert 'average' in str(ints_array)

        ints_array = IntSumarizedArray([0, 0, 1, 3])
        assert ints_array.calculate_distribution(bins=3)['counts'] == [2, 1, 1]

    def test_stats_functs(self):
        'It test the statistical functions of the class'
        ints = IntSumarizedArray([3, 5, 7, 7, 38])
        assert ints.median == 7

        ints = IntSumarizedArray([3, 5, 7, 38])
        assert ints.median == 6

        # median with two middle numbers
        ext_array = IntSumarizedArray([3, 5, 7, 7])
        assert ext_array.median == 6

        ints = IntSumarizedArray([34, 43, 81, 106, 106, 115])
        assert ints.average - 80.83 < 0.01

        ext_array = self.create_test_array()
        assert ext_array.median == 145
        assert round(ext_array.average, 2) == 145.15

        assert ext_array.sum == 13354
        assert ext_array.count == 92
        assert round(ext_array.variance, 2) == 557.43

        ints = IntSumarizedArray([3, 4, 4, 5, 6, 8, 8])
        assert ints.median == 5
        assert ints.quartiles == (4, 5, 8)

        ints = IntSumarizedArray([1, 2, 3, 4, 5])
        assert ints.quartiles == (1.5, 3, 4.5)

        ints = IntSumarizedArray([1, 2, 3, 4, 5, 6])
        assert ints.quartiles == (1.5, 3.5, 5.5)

        ints = IntSumarizedArray([1, 2, 3, 4, 5, 6, 7])
        assert ints.quartiles == (2, 4, 6)

        ints = IntSumarizedArray([1, 2, 3, 4, 5, 6, 7, 8])
        assert ints.quartiles == (2.5, 4.5, 6.5)

        assert ints.irq == 4.0
        assert ints.outlier_limits == (-3, 12)

        try:
            ints = IntSumarizedArray([0, 1, 2])
            assert ints.quartiles
            self.fail('RuntimeError')
        except RuntimeError:
            pass

    def test_value_for_index_test(self):
        'We can get the integer for a given index'
        # pylint: disable=W0212
        ints = IntSumarizedArray([3, 5, 7, 7, 38])
        assert ints._get_value_for_index(0) == 3
        assert ints._get_value_for_index(1) == 5
        assert ints._get_value_for_index(2) == 7
        assert ints._get_value_for_index(3) == 7
        assert ints._get_value_for_index(4) == 38
        try:
            assert ints._get_value_for_index(5) == 38
            self.fail('IndexError expected')
        except IndexError:
            pass


class IntsBoxplot(unittest.TestCase):
    'It tests the boxplot for integers'
    def test_boxplot(self):
        'It does a bloxplot for integers'
        box = IntBoxplot()
        box.append(1, 50)
        box.append(1, 40)
        box.append(1, 30)
        box.append(1, 40)
        box.append(2, 30)
        box.append(2, 10)
        box.append(2, 20)
        box.append(2, 40)
        box.append('no distrib', 40)
        counts = box.aggregated_array
        assert len(list(counts.flat)) == 9

        plot = box.ascii_plot
        assert '2 <---------' in plot

if __name__ == '__main__':
    #import sys;sys.argv = ['', 'IntsBoxplot.test_boxplot']
    unittest.main()
