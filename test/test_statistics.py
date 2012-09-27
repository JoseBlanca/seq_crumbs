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
from os.path import join
import unittest
from subprocess import check_output
from tempfile import NamedTemporaryFile
import operator

from crumbs.statistics import (IntCounter, draw_histogram, IntBoxplot,
                               calculate_sequence_stats, NuclFreqsPlot)
from crumbs.utils.test_utils import TEST_DATA_DIR
from crumbs.utils.bin_utils import BIN_DIR
from crumbs.seqio import read_seqrecords


class HistogramTest(unittest.TestCase):
    'It tests the ASCII histogram plotting'
    def test_ascii_histogram(self):
        'It plots an ASCII histogram'
        hist = draw_histogram(bin_limits=[-2, -1, 0, 1, 2],
                              counts=[9, 20, 30, 40])
        assert '[-2 , -1[ ( 9): ****************' in hist


class CounterTest(unittest.TestCase):
    'It tests that we can use a Counter'
    @staticmethod
    def create_test_counter():
        counter = IntCounter()
        d = {'9': '5', '10': '288', '11': '002556688', '12': '00012355555',
             '13': '0000013555688', '14': '00002555558',
             '15': '0000000000355555555557', '16': '000045', '17': '000055',
             '18': '0005', '19': '00005', '21': '5'}

        for key, values in d.items():
            for num in values:
                counter[int(key + num)] += 1
        return counter

    @staticmethod
    def test_counter():
        'create a counter'
        # initialize with values
        counter = IntCounter({2: 2})
        counter[2] += 1
        counter[6] += 1
        assert counter.min == 2
        assert counter.max == 6

        counter = IntCounter({3: 1, 5: 1, 7: 2, 38: 1})
        assert counter.min == 3
        assert counter.max == 38
        assert counter.sum == 60
        assert counter.count == 5
        assert counter.median == 7

    @staticmethod
    def test__add__():
        ext_counter = IntCounter({6: 1, 2: 1})
        ext_counter2 = IntCounter({7: 1, 2: 1})

        new_array = ext_counter + ext_counter2
        assert new_array[6] == 1
        assert new_array[7] == 1
        assert new_array[2] == 2

        #assert list(new_array.flat) == [2, 2, 6, 7]
    def test_stats_functs(self):
        'It test the statistical functions of the class'
        ints = IntCounter({3: 1, 5: 1, 7: 2, 38: 1})
        assert ints.median == 7

        ints = IntCounter({3: 1, 5: 1, 7: 1, 38: 1})
        assert ints.median == 6

        # median with two middle numbers
        ext_array = IntCounter({3: 1, 5: 1, 7: 2})
        assert ext_array.median == 6

        ints = IntCounter({34: 1, 43: 1, 81: 1, 106: 2, 115: 1})
        assert ints.average - 80.83 < 0.01

        ext_counter = self.create_test_counter()
        assert ext_counter.median == 145
        assert round(ext_counter.average, 2) == 145.15

        assert ext_counter.sum == 13354
        assert ext_counter.count == 92
        assert round(ext_counter.variance, 2) == 557.43

        ints = IntCounter({3: 1, 4: 2, 5: 1, 6: 1, 8: 2})
        assert ints.median == 5
        assert ints.quartiles == (4, 5, 8)

        ints = IntCounter({1: 1, 2: 1, 3: 1, 4: 1, 5: 1})
        assert ints.quartiles == (1.5, 3, 4.5)

        ints = IntCounter({1: 1, 2: 1, 3: 1, 4: 1, 5: 1, 6: 1})
        assert ints.quartiles == (1.5, 3.5, 5.5)

        ints = IntCounter({1: 1, 2: 1, 3: 1, 4: 1, 5: 1, 6: 1, 7: 1})
        assert ints.quartiles == (2, 4, 6)

        ints = IntCounter({1: 1, 2: 1, 3: 1, 4: 1, 5: 1, 6: 1, 7: 1, 8: 1})
        assert ints.quartiles == (2.5, 4.5, 6.5)
        assert ints.irq == 4.0
        assert ints.outlier_limits == (-3, 12)

        try:
            ints = IntCounter({0: 1, 1: 1, 2: 1})
            assert ints.quartiles
            self.fail('RuntimeError')
        except RuntimeError:
            pass

    @staticmethod
    def test_sum_with_treshold_function():
        'It tests the function that calculates Q30 and Q20'

        quals = IntCounter({15: 10, 21: 13, 30: 12})
        assert quals.count_relative_to_value(20, operator.ge) == 25
        assert quals.count_relative_to_value(30, operator.ge) == 12

    def test_value_for_index_test(self):
        'We can get the integer for a given index'
        # pylint: disable=W0212
        ints = IntCounter({3: 1, 5: 1, 7: 2, 38: 1})
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

    def test_distribution(self):
        'It tests the histogram function'

        ints_counter = self.create_test_counter()
        distrib = ints_counter.calculate_distribution(bins=10,
                                                        remove_outliers=5)

        assert distrib['counts'] == [7L, 13L, 7L, 10L, 7L, 22L, 6L, 4L, 5L,
                                      5L]
        assert distrib['bin_limits'] == [110, 118, 126, 134, 142, 150, 158,
                                         166, 174, 182, 190]
        assert 'average' in str(ints_counter)

        ints_counter = IntCounter({0: 2, 1: 1, 3: 1})
        result = [2, 1, 1]
        assert ints_counter.calculate_distribution(bins=3)['counts'] == result


class IntsBoxplotTest(unittest.TestCase):
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
        assert sum(counts.values()) == 9

        plot = box.ascii_plot
        assert '2:10.0,15.0,25.0,35.0,40.0 <----------[=========' in plot


class CalculateStatsTest(unittest.TestCase):
    'It tests the calculate stats functions'

    def make_fasta(self):
        'it returns a fasta fhand'
        fhand = NamedTemporaryFile()
        fhand.write('>seq\nACTATCATGGCAGATA\n')
        fhand.flush()
        return fhand

    @staticmethod
    def test_calculate_stats():
        'It tests the calculate stat function'
        in_fhands = []
        for val in range(1, 6):
            fhand = open(join(TEST_DATA_DIR, 'pairend{0}.sfastq'.format(val)))
            in_fhands.append(fhand)
        seqs = read_seqrecords(in_fhands, file_format='fastq')
        (lengths_srt, qual_str, freq_str,
                                 qual_boxplot) = calculate_sequence_stats(seqs)
        assert 'maximum: 4' in lengths_srt
        assert 'Q30: 100.0' in qual_str
        assert '1:30.0,30.0,30.0,30.0,30.0 <[|]>' in qual_boxplot
        assert '[30 , 31[ (96): **********' in qual_str
        assert '0 (A: 1.00, C: 0.00, G: 0.00, T: 0.00, N: 0.00) |' in  freq_str

    def test_stats_bin(self):
        'It tests the statistics binary'

        bin_ = join(BIN_DIR, 'calculate_stats')

        # help
        assert 'usage' in check_output([bin_, '-h'])

        # fasta
        in_fhand1 = self.make_fasta()
        in_fhand2 = self.make_fasta()
        result = check_output([bin_, in_fhand1.name, in_fhand2.name])
        assert 'Length stats and distribution.\n-------------------' in result

        # fastq
        cmd = [bin_]
        for val in range(1, 6):
            cmd.append(join(TEST_DATA_DIR, 'pairend{0}.sfastq'.format(val)))
        assert 'Quality stats and distribution' in check_output(cmd)


class BaseFreqPlotTest(unittest.TestCase):
    'It tests the nucleotide frequency plot'
    def test_base_freq_plot(self):
        'It does a NuclBaseFreqPlot'
        plot = NuclFreqsPlot()
        plot.append(0, 'A')
        plot.append(0, 'T')
        plot.append(0, 'A')
        plot.append(1, 'N')
        plot.append(1, 'Y')
        plot.append(1, 'C')
        ascii = plot.ascii_plot
        expected = '0 (A: 0.67, C: 0.00, G: 0.00, T: 0.33, N: 0.00) |aaaaaaaaa'
        assert expected in ascii

if __name__ == '__main__':
    #import sys;sys.argv = ['', 'CounterTest']
    unittest.main()
