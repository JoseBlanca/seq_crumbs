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

'''
Created on 18/06/2012

@author: jose
'''

import unittest
import os.path
from subprocess import check_output, CalledProcessError
from tempfile import NamedTemporaryFile

from crumbs.tests.utils import TEST_DATA_DIR, BIN_DIR

# pylint: disable=R0201
# pylint: disable=R0904


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
