
import os.path
import unittest
from subprocess import check_output, check_call
from tempfile import NamedTemporaryFile

from bam_crumbs.utils.test import TEST_DATA_DIR
from bam_crumbs.utils.bin import BIN_DIR
from bam_crumbs.bam_tools import filter_bam

# pylint: disable=C0111


class SortTest(unittest.TestCase):
    def test_sort_bam_bin(self):
        bin_ = os.path.join(BIN_DIR, 'sort_bam')
        assert 'usage' in check_output([bin_, '-h'])

        bam_fpath = os.path.join(TEST_DATA_DIR, 'seqs.bam')
        sorted_fhand = NamedTemporaryFile(suffix='.sorted.bam')
        check_call([bin_, bam_fpath, '-o', sorted_fhand.name])
        assert "@HD\tVN:1.4" in check_output(['samtools', 'view', '-h',
                                              sorted_fhand.name])
        assert os.path.exists(sorted_fhand.name + '.bai')
        os.remove(sorted_fhand.name + '.bai')
        sorted_fhand.close()

        # no index
        sorted_fhand = NamedTemporaryFile()
        check_call([bin_, bam_fpath, '-o', sorted_fhand.name, '--no-index'])
        assert not os.path.exists(sorted_fhand.name + '.bai')

        # sort the sam file
        fhand = NamedTemporaryFile()
        fhand.write(open(bam_fpath).read())
        fhand.flush()
        check_call([bin_, fhand.name])
        assert "@HD\tVN:1.4" in check_output(['samtools', 'view', '-h',
                                              bam_fpath])
        assert os.path.exists(fhand.name + '.bai')
        os.remove(fhand.name + '.bai')


class FilterTest(unittest.TestCase):
    def test_filter_mapq(self):
        bam_fpath = os.path.join(TEST_DATA_DIR, 'seqs.bam')
        out_fhand = NamedTemporaryFile()
        filter_bam(bam_fpath, out_fhand.name, min_mapq=100)
        assert len(open(out_fhand.name).read(20)) == 20


if __name__ == "__main__":
    # import sys;sys.argv = ['', 'FilterTest.test_filter_mapq']
    unittest.main()
