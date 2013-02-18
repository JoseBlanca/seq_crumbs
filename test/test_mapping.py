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

import unittest
import subprocess
import os.path
from tempfile import NamedTemporaryFile

from crumbs.utils.test_utils import TEST_DATA_DIR
from crumbs.mapping import (get_or_create_bowtie2_index, _bowtie2_index_exists,
                            map_with_bowtie2, get_or_create_bwa_index,
                            _bwa_index_exists, map_with_bwa_sw)
from crumbs.utils.file_utils import TemporaryDir
from crumbs.utils.bin_utils import get_binary_path


# pylint: disable=R0201
# pylint: disable=R0904
# pylint: disable=W0402
# pylint: disable=C0111


class Bowtie2Test(unittest.TestCase):
    def test_get_or_create_index(self):
        db_name = 'arabidopsis_genes'
        seq_fpath = os.path.join(TEST_DATA_DIR, db_name)
        assert not _bowtie2_index_exists(seq_fpath)

        directory = TemporaryDir()
        index_fpath = get_or_create_bowtie2_index(seq_fpath, directory.name)
        expected_index = os.path.join(directory.name,
                                      os.path.basename(db_name))
        assert index_fpath == expected_index
        assert _bowtie2_index_exists(index_fpath)

        # already exists
        index_fpath = get_or_create_bowtie2_index(seq_fpath, directory.name)
        assert index_fpath == expected_index
        assert _bowtie2_index_exists(index_fpath)
        directory.close()

    def test_map_with_bowtie2(self):
        reference_fpath = os.path.join(TEST_DATA_DIR, 'arabidopsis_genes')
        reads_fpath = os.path.join(TEST_DATA_DIR, 'arabidopsis_reads.fastq')
        directory = TemporaryDir()
        index_fpath = get_or_create_bowtie2_index(reference_fpath,
                                                  directory.name)
        bam_fhand = NamedTemporaryFile(suffix='.bam')
        map_with_bowtie2(index_fpath, bam_fhand.name,
                         unpaired_fpaths=[reads_fpath])

        directory.close()


class Bwa2Test(unittest.TestCase):
    def test_get_or_create_index(self):
        db_name = 'arabidopsis_genes'
        seq_fpath = os.path.join(TEST_DATA_DIR, db_name)
        assert not _bwa_index_exists(seq_fpath)

        directory = TemporaryDir()
        index_fpath = get_or_create_bwa_index(seq_fpath, directory.name)
        expected_index = os.path.join(directory.name,
                                      os.path.basename(db_name))
        assert index_fpath == expected_index
        assert _bwa_index_exists(index_fpath)

        # already exists
        index_fpath = get_or_create_bwa_index(seq_fpath, directory.name)
        assert index_fpath == expected_index
        assert _bwa_index_exists(index_fpath)
        directory.close()

    def test_map_with_bwa(self):
        reference_fpath = os.path.join(TEST_DATA_DIR, 'arabidopsis_genes')
        reads_fpath = os.path.join(TEST_DATA_DIR, 'arabidopsis_reads.fastq')
        directory = TemporaryDir()
        index_fpath = get_or_create_bwa_index(reference_fpath, directory.name)
        bam_fhand = NamedTemporaryFile(suffix='.bam')
        map_with_bwa_sw(index_fpath, bam_fhand.name, unpaired_fpath=reads_fpath)
        out = subprocess.check_output([get_binary_path('samtools'), 'view',
                                       bam_fhand.name])
        assert  'TTCTGATTCAATCTACTTCAAAGTTGGCTTTATCAATAAG' in out

        directory.close()
if __name__ == '__main__':
    # import sys;sys.argv = ['', 'BlastTest.test_blastdb']
    unittest.main()
