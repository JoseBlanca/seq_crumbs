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
# pylint: disable=C0111

from os.path import join
import unittest
from subprocess import check_output
from tempfile import NamedTemporaryFile
import operator

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from crumbs.statistics import (IntCounter, draw_histogram, IntBoxplot,
                               calculate_sequence_stats, NuclFreqsPlot,
                               KmerCounter, calculate_dust_score,
                               calculate_nx, BestItemsKeeper)
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

        # assert list(new_array.flat) == [2, 2, 6, 7]
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
        assert '2:10.0,15.0,25.0,35.0,40.0 <-----[============|=======' in plot


class KmerCounterTest(unittest.TestCase):
    'It tests the kmer counter test'
    @staticmethod
    def test_kmer_counter():
        kmers = KmerCounter(3)
        kmers.count_seq('ATCATGGCTACGACT')
        assert list(kmers.values) == [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
        kmers.count_seq('ATCATGGCTACGACT')
        assert list(kmers.values) == [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]


class DustCalculationTest(unittest.TestCase):
    'It calculates dust scores'
    @staticmethod
    def test_dustscore_calculation():
        'It calculates the dust score'
        seqs = ['TTTTTTTTTTTTTTTTTTTTTTTTTTTT', 'TATATATATATATATATATATATATATA',
                'GAAGAAGAAGAAGAAGAAGAAGAAGAAG', 'AACTGCAGTCGATGCTGATTCGATCGAT',
                'AACTGAAAAAAAATTTTTTTAAAAAAAA']

        # short sequences
        scores = [100, 48, 30.76, 4.31, 23.38]
        scoresx3 = [100, 48.68, 28.65, 5.62, 27.53]
        scoresx4 = [100, 48.55, 28.25, 5.79, 28.00]
        for seq, score, scorex3, scorex4 in zip(seqs, scores, scoresx3,
                                                scoresx4):
            seqrec = SeqRecord(Seq(seq))
            assert calculate_dust_score(seqrec) - score < 0.01
            seqrec = SeqRecord(Seq(seq * 3))
            assert calculate_dust_score(seqrec) - scorex3 < 0.01
            seqrec = SeqRecord(Seq(seq * 4))
            assert calculate_dust_score(seqrec) - scorex4 < 0.01


class NxCalculationTest(unittest.TestCase):
    'It calculates N50 and N95'
    @staticmethod
    def test_n50_calculation():
        'It calculates N50.'
        assert calculate_nx(IntCounter([2, 2, 2, 3, 3, 4, 8, 8]), 50) == 8
        assert calculate_nx(IntCounter([2, 2, 2, 3, 3, 4, 8, 8]), 95) == 2
        assert calculate_nx(IntCounter([8, 8, 8, 8, 8, 8, 8, 8]), 50) == 8
        assert calculate_nx(IntCounter(), 50) is None


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
        results = calculate_sequence_stats(seqs, nxs=[50])
        assert 'maximum: 4' in results['length']
        assert 'N50' in results['length']
        assert '1:30.0,30.0,30.0,30.0,30.0 <[|]>' in results['qual_boxplot']
        assert '[30 , 31[ (96): **********' in results['quality']
        assert 'Q30: 100.0' in results['quality']
        assert '0 (A: 1.00, C: 0.00, G: 0.00, T: 0.00' in  results['nucl_freq']
        assert results['kmer'] == ''

        infhands = [open(join(TEST_DATA_DIR, 'arabidopsis_genes'))]
        seqs = list(read_seqrecords(infhands, file_format='fasta'))
        kmers = calculate_sequence_stats(seqs)['kmer']
        assert not 'Kmer distribution' in kmers

        kmers = calculate_sequence_stats(seqs, kmer_size=3)['kmer']
        assert 'Kmer distribution' in kmers
        assert 'TCT: 167' in kmers

        # dust
        dust = calculate_sequence_stats(seqs)['dustscore']
        assert not dust
        dust = calculate_sequence_stats(seqs, do_dust_stats=True)['dustscore']
        assert 'average: 1.83\nvariance: 0.14\nnum. seqs.: 6\n' in dust
        assert '% above 7 (low complexity): 0.00' in dust

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
        assert 'N95:' in result

        # fastq
        cmd = [bin_]
        for val in range(1, 6):
            cmd.append(join(TEST_DATA_DIR, 'pairend{0}.sfastq'.format(val)))
        result = check_output(cmd)
        assert 'Quality stats and distribution' in result
        assert 'Kmer distribution' not in result

        # kmer distribution
        cmd = [bin_, '-c', '-k', '3']
        for val in range(1, 6):
            cmd.append(join(TEST_DATA_DIR, 'pairend{0}.sfastq'.format(val)))
        result = check_output(cmd)
        assert 'Quality stats and distribution' in result
        assert 'Kmer distribution' in result
        assert 'aaa: 48' in result

        # kmer distribution
        cmd = [bin_, '-k', '3']
        for val in range(1, 6):
            cmd.append(join(TEST_DATA_DIR, 'pairend{0}.sfastq'.format(val)))
        result = check_output(cmd)
        assert 'Kmer distribution' in result

        # dustscore distribution
        cmd = [bin_, '-d']
        cmd.append(join(TEST_DATA_DIR, 'arabidopsis_genes'))
        result = check_output(cmd)
        assert 'Dustscores' in result


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
        expected = '0 (A: 0.67, C: 0.00, G: 0.00, T: 0.33, N: 0.00) | aaaaaaaa'
        assert expected in ascii

        plot = NuclFreqsPlot()
        plot.append(0, 'A')
        plot.append(0, 'A')
        plot.append(0, 'A')
        plot.append(0, 'C')
        plot.append(0, 'T')
        plot.append(0, 'T')
        plot.append(0, 'T')

        ascii = plot.ascii_plot
        assert '0 (A: 0.43, C: 0.14, G: 0.00, T: 0.43, N: 0.00) | aaa' in ascii


class BestItemsKeeperTest(unittest.TestCase):
    def test_best_items(self):
        items = reversed(range(1, 10))
        best_items = BestItemsKeeper(5, items)
        assert best_items == [5, 6, 7, 8, 9]
        best_items.add(10)
        assert best_items == [6, 7, 8, 9, 10]
        assert [6, 7, 8, 9, 10] == best_items
        assert best_items[0] == 6
        assert best_items[:1] == [6]
        best_items.add(6)
        assert best_items == [6, 7, 8, 9, 10]

        best_items.add(6)

        items = [[i] for i in range(1, 10)]
        key = lambda x: -(x[0])
        best_items = BestItemsKeeper(5, key=key)
        best_items.update(items)
        assert best_items == [[5], [4], [3], [2], [1]]

        items = [[i] for i in range(1, 10)]
        key = lambda x: x[0]
        best_items = BestItemsKeeper(5, key=key, reverse=True)
        best_items.update(items)
        assert best_items == [[5], [4], [3], [2], [1]]



if __name__ == '__main__':
    # import sys;sys.argv = ['', 'CalculateStatsTest.test_stats_bin']
    unittest.main()
