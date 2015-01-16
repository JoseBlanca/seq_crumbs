'''
Created on 2014 api 17

@author: peio
'''
import unittest
import os

from subprocess import check_call, check_output
from tempfile import NamedTemporaryFile

from crumbs.utils.bin_utils import BIN_DIR
from crumbs.utils.test_utils import TEST_DATA_DIR


class BinTests(unittest.TestCase):
    def test_draw_coverage(self):
        bin_ = os.path.join(BIN_DIR, 'draw_coverage_hist')
        bam_fpath = os.path.join(TEST_DATA_DIR, 'seqs.bam')
        fhand = NamedTemporaryFile(suffix='.png')
        out = check_output([bin_, bam_fpath, '-o', fhand.name])
        assert '147' in out

    def test_mapq_hist(self):
        bin_ = os.path.join(BIN_DIR, 'draw_mapq_hist')
        bam_fpath = os.path.join(TEST_DATA_DIR, 'seqs.bam')
        fhand = NamedTemporaryFile(suffix='.png')
        null = open(os.devnull, 'w')
        check_call([bin_, bam_fpath, '-o', fhand.name], stdout=null)
        # raw_input(fhand.name)

        fhand = NamedTemporaryFile(suffix='.png')
        check_call([bin_, bam_fpath, '-o', fhand.name, '-t'], stdout=null)
        res = open(fhand.name).read()
        assert "[147 , 154[ (3):" in res

if __name__ == "__main__":
    # import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
