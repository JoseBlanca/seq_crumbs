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
from StringIO import StringIO
from tempfile import NamedTemporaryFile
from subprocess import check_output

from crumbs.trim import (_get_uppercase_segments, _get_longest_segment,
                         trim_lowercased_seqs)
from crumbs.seqio import read_seqrecords
from crumbs.tests.utils import BIN_DIR

FASTQ = '@seq1\naTCgt\n+\n?????\n@seq2\natcGT\n+\n?????\n'


class TrimTest(unittest.TestCase):
    'It tests the rim functions'

    @staticmethod
    def test_masked_locations():
        'It test the masked loacations function'

        seq = 'aaATTTTTTaa'
        assert list(_get_uppercase_segments(seq)) == [(2, 8)]

        seq = 'aaATTTaTTaa'
        assert list(_get_uppercase_segments(seq)) == [(2, 5), (7, 8)]

        seq = 'AAATaaa'
        assert list(_get_uppercase_segments(seq)) == [(0, 3)]

        seq = 'aaaaAAAA'
        assert list(_get_uppercase_segments(seq)) == [(4, 7)]

    @staticmethod
    def test_get_longest_section():
        'It gets the longest section from a list of sections'

        segments = [(0, 3), (10, 34)]
        assert (10, 34) == _get_longest_segment(segments)

        segments = [(0, 3), (10, 13)]
        segment = _get_longest_segment(segments)
        assert segment == (0, 3) or segment == (10, 13)

    @staticmethod
    def test_trim_seqs():
        'it tests the trim seq function'
        fasta = '>seq1\naaCTTTC\n>seq2\nCTTCaa\n>seq3\naaCTCaa\n>s\nactg\n'
        seqs = read_seqrecords([StringIO(fasta)])
        seqs = [str(seq.seq) for seq in trim_lowercased_seqs(seqs)]
        assert seqs == ['CTTTC', 'CTTC', 'CTC']


class TrimByCaseTest(unittest.TestCase):
    'It tests the trim_by_case binary'
    @staticmethod
    def _make_fhand(content=None):
        'It makes temporary fhands'
        if content is None:
            content = ''
        fhand = NamedTemporaryFile()
        fhand.write(content)
        fhand.flush()
        return fhand

    def test_trim_case_bin(self):
        'It tests the trim seqs binary'
        trim_bin = os.path.join(BIN_DIR, 'trim_by_case')
        assert check_output([trim_bin]).startswith('usage')

        fastq_fhand = self._make_fhand(FASTQ)

        result = check_output([trim_bin, fastq_fhand.name])
        assert '@seq1\nTC\n+' in result


if __name__ == '__main__':
    #import sys;sys.argv = ['', 'SffExtractTest.test_items_in_gff']
    unittest.main()
