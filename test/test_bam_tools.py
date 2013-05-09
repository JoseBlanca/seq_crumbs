
import os.path
import unittest
from subprocess import check_output, check_call
from tempfile import NamedTemporaryFile

from bam_crumbs.utils.test import TEST_DATA_DIR
from bam_crumbs.utils.bin import BIN_DIR
from bam_crumbs.bam_tools import filter_bam, realign_bam, calmd_bam

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


class RealignTest(unittest.TestCase):
    def test_realign_bamself(self):
        ref_fpath = os.path.join(TEST_DATA_DIR, 'CUUC00007_TC01.fasta')
        bam_fpath = os.path.join(TEST_DATA_DIR, 'sample.bam')
        out_bam = NamedTemporaryFile()
        realign_bam(bam_fpath, ref_fpath, out_bam.name)

    def xtest_realign_bin(self):
        bin_ = os.path.join(BIN_DIR, 'realign_bam')
        assert 'usage' in check_output([bin_, '-h'])

        bam_fpath = os.path.join(TEST_DATA_DIR, 'sample.bam')
        ref_fpath = os.path.join(TEST_DATA_DIR, 'CUUC00007_TC01.fasta')
        realigned_fhand = NamedTemporaryFile(suffix='.realigned.bam')
        check_call([bin_, bam_fpath, '-o', realigned_fhand.name, '-f',
                    ref_fpath])
        assert open(realigned_fhand.name).read()

        # in parallel
        realigned_fhand = NamedTemporaryFile(suffix='.realigned.bam')
        check_call([bin_, bam_fpath, '-o', realigned_fhand.name, '-f',
                    ref_fpath, '-t', '2'])
        assert open(realigned_fhand.name).read()


class CalmdTest(unittest.TestCase):
    def test_calmd_bam(self):
        ref_fpath = os.path.join(TEST_DATA_DIR, 'CUUC00007_TC01.fasta')
        bam_fpath = os.path.join(TEST_DATA_DIR, 'sample.bam')
        out_bam = NamedTemporaryFile()
        calmd_bam(bam_fpath, ref_fpath, out_bam.name)
        assert  open(out_bam.name).read()

    def test_calmd_bin(self):
        bin_ = os.path.join(BIN_DIR, 'calmd_bam')
        assert 'usage' in check_output([bin_, '-h'])

        bam_fpath = os.path.join(TEST_DATA_DIR, 'sample.bam')
        ref_fpath = os.path.join(TEST_DATA_DIR, 'CUUC00007_TC01.fasta')
        calmd_fhand = NamedTemporaryFile(suffix='.calmd.bam')
        check_call([bin_, bam_fpath, '-o', calmd_fhand.name, '-f',
                    ref_fpath])
        assert open(calmd_fhand.name).read()

if __name__ == "__main__":
    import sys;sys.argv = ['', 'CalmdTest']
    unittest.main()
