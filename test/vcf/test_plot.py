'''
Created on 2014 uzt 24

@author: peio
'''
from os.path import join
import itertools
import random
import unittest
from collections import OrderedDict

from crumbs.utils.file_utils import TemporaryDir
from crumbs.vcf.plot import plot_in_genome, plot_haplotypes
from crumbs.utils.test_utils import TEST_DATA_DIR
from tempfile import NamedTemporaryFile


class GenomePlot(unittest.TestCase):
    def gen_windows(self, chrom='ch1'):
        step = 10
        for x in range(0, 1000, step):
            start = x
            end = x + step
            value = random.random()
            value2 = random.random()
            yield {'start': start, 'end': end, 'chrom': chrom,
                   'values': {'val1': value, 'val2': value2}}

    def test_plot_window(self):
        iterator = itertools.chain(self.gen_windows(),
                                   self.gen_windows('ch2'))
        tempdir = TemporaryDir()
        out_base = join(tempdir.name, 'out')
        labels = OrderedDict({'val1': {'title': 'val1 title',
                                       'ylabel': 'val1 ylabel'},
                              'val2': {'title': 'val2 title',
                                       'ylabel': 'val2 ylabel'}})

        plot_in_genome(iterator, out_base=out_base, labels=labels)
        # raw_input(tempdir.name)
        tempdir.close()


class HaplotypePlot(unittest.TestCase):
    def test_plot_haplo(self):
        vcf_fhand = open(join(TEST_DATA_DIR, 'freebayes_multisample.vcf.gz'))
        plot_fhand = NamedTemporaryFile(suffix='.png')
        plot_haplotypes(vcf_fhand, plot_fhand)
        # raw_input(plot_fhand.name)

if __name__ == "__main__":
    #import sys; sys.argv = ['', 'HaplotypePlot']
    unittest.main()
