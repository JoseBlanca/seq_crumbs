
import os.path
import unittest
from subprocess import check_output, check_call

from bam_crumbs.utils.test import TEST_DATA_DIR
from bam_crumbs.utils.bin import BIN_DIR
from tempfile import NamedTemporaryFile


class SortTest(unittest.TestCase):
    def test_sort_bam_bin(self):
        bin_ = bin_ = os.path.join(BIN_DIR, 'sort_bam')
        assert 'usage' in check_output([bin_, '-h'])

        bam_fpath = os.path.join(TEST_DATA_DIR, 'seqs.bam')
        sorted_fhand = NamedTemporaryFile()
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

if __name__ == "__main__":
    # import sys;sys.argv = ['', 'PlotTest.test_geographic_map']
    unittest.main()
