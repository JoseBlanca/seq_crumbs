
import os.path
import unittest
from subprocess import check_output

import pysam

from bam_crumbs.utils.test import TEST_DATA_DIR
from bam_crumbs.utils.bin import BIN_DIR
from bam_crumbs.statistics import (count_reads, RpkmCounter, MapqCounter,
                                   CoverageCounter)

# pylint: disable=R0201
# pylint: disable=R0904
# pylint: disable=C0111


class CountTest(unittest.TestCase):
    def test_count(self):
        bam_fpath = os.path.join(TEST_DATA_DIR, 'seqs.bam')
        bam = pysam.Samfile(bam_fpath)
        assert count_reads('reference1', bam) == 9
        assert count_reads('reference2', bam) == 9
        assert count_reads('reference1', bam, start=1, end=10) == 0
        assert count_reads('reference1', bam, start=0, end=500) == 9

    def test_rpkm_distrib(self):
        bam_fpath = os.path.join(TEST_DATA_DIR, 'seqs.bam')
        bam = pysam.Samfile(bam_fpath)
        rpkms = RpkmCounter(bam)
        assert '(1)' in rpkms.ascii_histogram()

    def test_rpkm_bin(self):
        bam_fpath = os.path.join(TEST_DATA_DIR, 'seqs.bam')

        bin_ = os.path.join(BIN_DIR, 'calculate_rpkm_distrib')
        # help
        assert 'usage' in check_output([bin_, '-h'])

        assert '(1)' in check_output([bin_, bam_fpath])

    def test_mapq_distrib(self):
        bam_fpath = os.path.join(TEST_DATA_DIR, 'seqs.bam')
        bam = pysam.Samfile(bam_fpath)
        mapqs = MapqCounter(bam)
        assert mapqs.count == 18
        assert mapqs.min == 28
        assert mapqs.max == 149

    def test_coverage_distrib(self):
        bam_fpath = os.path.join(TEST_DATA_DIR, 'seqs.bam')
        bam = pysam.Samfile(bam_fpath)
        cov = CoverageCounter(bam)
        assert cov.count == 147
        assert cov.min == 6
        assert cov.max == 9


if __name__ == "__main__":
    # import sys;sys.argv = ['', 'PlotTest.test_geographic_map']
    unittest.main()
