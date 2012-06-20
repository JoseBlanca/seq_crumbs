'''
Created on 18/06/2012

@author: jose
'''

import unittest
import os.path
from subprocess import check_output, CalledProcessError
from tempfile import NamedTemporaryFile

from crumbs.tests.utils import TEST_DATA_DIR

# pylint: disable=R0201
# pylint: disable=R0904


BIN_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))


class SffExtractBinTest(unittest.TestCase):
    'It tests the sff_extract binary'

    def test_extract_sff(self):
        'It tests the sff_extract binary'
        sff_bin = os.path.join(BIN_DIR, 'sff_extract')
        assert check_output([sff_bin]).startswith('usage')

        # clipping warning
        sff_fpath = os.path.join(TEST_DATA_DIR, '10_454_reads.sff')
        cmd = [sff_bin, sff_fpath]
        stderr = NamedTemporaryFile()
        try:
            check_output(cmd, stderr=stderr)
            self.fail()
        except CalledProcessError:
            assert 'Countermeasures' in open(stderr.name).read()

        # min left clip
        cmd = [sff_bin, '--min_left_clip', '5', sff_fpath]
        check_output(cmd)

        # clip
        cmd = [sff_bin, '--min_left_clip', '5', '--clip', sff_fpath]
        stdout = check_output(cmd)
        assert stdout.startswith('@E3MFGYR02JWQ7T\nGTCTACATGTTG')

        # file does not exist
        cmd = [sff_bin, 'no_file']
        stderr = NamedTemporaryFile()
        try:
            check_output(cmd, stderr=stderr)
            self.fail('Error expected')
        except CalledProcessError:
            assert 'A file was not found: no_file' in open(stderr.name).read()

        #version

if __name__ == '__main__':
    #import sys;sys.argv = ['', 'SffExtractTest.test_items_in_gff']
    unittest.main()
