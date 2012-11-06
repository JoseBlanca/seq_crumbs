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
import os.path
from tempfile import NamedTemporaryFile

from crumbs.blast import (do_blast, BlastMatcherForFewSubjects, get_or_create_blastdb,
                          _blastdb_exists)
from crumbs.utils.file_utils import TemporaryDir
from crumbs.settings import LINKERS, TITANIUM_LINKER
from crumbs.utils.test_utils import TEST_DATA_DIR
from crumbs.utils.tags import NUCL

# pylint: disable=R0201
# pylint: disable=R0904


class BlastTest(unittest.TestCase):
    'It tests the blast infraestructure'

    def test_blastdb(self):
        'It creates a blast database.'
        db_name = 'arabidopsis_genes'
        seq_fpath = os.path.join(TEST_DATA_DIR, db_name)
        db_dir = TemporaryDir(prefix='blast_dbs_')
        try:
            db_path1 = get_or_create_blastdb(seq_fpath, directory=db_dir.name,
                                    dbtype='nucl')
            db_path = os.path.join(db_dir.name, db_name)
            assert 'CATAGGGTCACCAATGGC' in open(db_path1).read(100)
            assert db_path1 == db_path
            assert os.path.exists(db_path)
            index_fpath = os.path.join(db_dir.name, db_name + '.nsq')
            assert os.path.exists(index_fpath)

        finally:
            db_dir.close()

    def test_blast_search(self):
        'It does a blast search'
        db_name = 'arabidopsis_genes'
        seq_fpath = os.path.join(TEST_DATA_DIR, db_name)
        db_dir = TemporaryDir(prefix='blast_dbs_')
        try:
            db_fpath = get_or_create_blastdb(seq_fpath, directory=db_dir.name,
                                              dbtype='nucl')
            query_fhand = NamedTemporaryFile()
            query_fhand.write(open(seq_fpath).read(200))
            query_fhand.flush()
            out_fhand = NamedTemporaryFile()
            do_blast(seq_fpath, db_fpath, program='blastn',
                     out_fpath=out_fhand.name)
            assert '</BlastOutput>' in open(out_fhand.name).read()
        finally:
            db_dir.close()

    @staticmethod
    def test_get_or_create_blastdb():
        'It test the blastdb kind'
        blastdb = os.path.join(TEST_DATA_DIR, 'arabidopsis_genes')

        directory = TemporaryDir()
        assert not _blastdb_exists(blastdb, NUCL)
        get_or_create_blastdb(blastdb, NUCL, directory.name)
        new_blast_path = os.path.join(directory.name,
                                      os.path.basename(blastdb))
        assert _blastdb_exists(new_blast_path, NUCL)
        get_or_create_blastdb(blastdb, NUCL, directory.name)
        assert _blastdb_exists(new_blast_path, NUCL)
        directory.close()

        # already exists
        blastdb = os.path.join(TEST_DATA_DIR, 'blastdbs', 'arabidopsis_genes')
        assert _blastdb_exists(blastdb, NUCL)
        get_or_create_blastdb(blastdb, NUCL)
        assert _blastdb_exists(blastdb, NUCL)


def create_a_matepair_file():
    'It creates a matepair fasta file'

    seq_5 = 'CTAGTCTAGTCGTAGTCATGGCTGTAGTCTAGTCTACGATTCGTATCAGTTGTGTGAC'
    seq_3 = 'ATCGATCATGTTGTATTGTGTACTATACACACACGTAGGTCGACTATCGTAGCTAGT'
    mate_seq = seq_5 + TITANIUM_LINKER + seq_3
    mate_fhand = NamedTemporaryFile(suffix='.fasta')
    mate_fhand.write('>seq1\n' + mate_seq + '\n')
    mate_fhand.flush()
    return mate_fhand


class BlastMater(unittest.TestCase):
    'It tests the splitting of mate pairs'
    def test_matching_segments(self):
        'It tests the detection of oligos in sequence files'
        seq_5 = 'CTAGTCTAGTCGTAGTCATGGCTGTAGTCTAGTCTACGATTCGTATCAGTTGTGTGAC'
        mate_fhand = create_a_matepair_file()

        expected_region = (len(seq_5), len(seq_5 + TITANIUM_LINKER) - 1)
        matcher = BlastMatcherForFewSubjects(mate_fhand.name, LINKERS,
                                             program='blastn',
                                             elongate_for_global=True)
        linker_region = matcher.get_matched_segments_for_read('seq1')[0]
        assert [expected_region] == linker_region

if __name__ == '__main__':
    # import sys;sys.argv = ['', 'BlastTest.test_blastdb']
    unittest.main()
