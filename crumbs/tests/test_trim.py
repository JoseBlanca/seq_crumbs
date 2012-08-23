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
from itertools import chain

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from crumbs.trim import (_get_uppercase_segments, _get_longest_segment,
                         TrimLowercasedLetters, TrimEdges)
from crumbs.seqio import read_seq_packets
from crumbs.tests.utils import BIN_DIR

FASTQ = '@seq1\naTCgt\n+\n?????\n@seq2\natcGT\n+\n?????\n'
FASTQ2 = '@seq1\nATCGT\n+\nA???A\n@seq2\nATCGT\n+\n?????\n'

# pylint: disable=R0201
# pylint: disable=R0904


def _make_fhand(content=''):
    'It makes temporary fhands'
    fhand = NamedTemporaryFile()
    fhand.write(content)
    fhand.flush()
    return fhand


class TrimTest(unittest.TestCase):
    'It tests the trim functions'

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
        'It tests the trim seq function'
        fasta = '>seq1\naaCTTTC\n>seq2\nCTTCaa\n>seq3\naaCTCaa\n>s\nactg\n'
        seq_packets = read_seq_packets([StringIO(fasta)])
        trim_lowercased_seqs = TrimLowercasedLetters()
        # pylint: disable=W0141
        seq_packets = map(trim_lowercased_seqs, seq_packets)
        seqs = [str(s.seq) for s in chain.from_iterable(seq_packets)]
        assert seqs == ['CTTTC', 'CTTC', 'CTC']


class TrimByCaseTest(unittest.TestCase):
    'It tests the trim_by_case binary'

    def test_trim_case_bin(self):
        'It tests the trim seqs binary'
        trim_bin = os.path.join(BIN_DIR, 'trim_by_case')
        assert 'usage' in check_output([trim_bin, '-h'])

        fastq_fhand = _make_fhand(FASTQ)

        result = check_output([trim_bin, fastq_fhand.name])
        assert '@seq1\nTC\n+' in result

    def test_trim_in_parallel(self):
        'It trims sequences in parallel'
        trim_bin = os.path.join(BIN_DIR, 'trim_by_case')
        fastq_fhand = _make_fhand(FASTQ)

        result = check_output([trim_bin, '-p', '2', fastq_fhand.name])
        assert '@seq1\nTC\n+' in result


class TrimEdgesTest(unittest.TestCase):
    'It test the fixed number of bases trimming'

    def _some_seqs(self):
        'It returns some seqrecords.'
        seqs = []
        seqs.append(SeqRecord(Seq('ACCG'),
                              letter_annotations={'dummy': 'dddd'}))
        seqs.append(SeqRecord(Seq('AAACCCGGG')))
        return seqs

    def test_edge_trimming(self):
        'It trims the edges'
        trim_edges = TrimEdges(left=1)
        res = [str(s.seq) for s in trim_edges(self._some_seqs())]
        assert res == ['CCG', 'AACCCGGG']

        trim_edges = TrimEdges(right=1)
        res = [str(s.seq) for s in trim_edges(self._some_seqs())]
        assert res == ['ACC', 'AAACCCGG']

        trim_edges = TrimEdges(left=1, right=1)
        res = [str(s.seq) for s in trim_edges(self._some_seqs())]
        assert res == ['CC', 'AACCCGG']

        trim_edges = TrimEdges(left=2, right=2)
        res = [str(s.seq) for s in trim_edges(self._some_seqs())]
        assert res == ['ACCCG']

        trim_edges = TrimEdges(left=3, right=3)
        res = [str(s.seq) for s in trim_edges(self._some_seqs())]
        assert res == ['CCC']

        trim_edges = TrimEdges(left=1, mask=True)
        res = [str(s.seq) for s in trim_edges(self._some_seqs())]
        assert res == ['aCCG', 'aAACCCGGG']

        trim_edges = TrimEdges(right=1, mask=True)
        res = [str(s.seq) for s in trim_edges(self._some_seqs())]
        assert res == ['ACCg', 'AAACCCGGg']

        trim_edges = TrimEdges(left=1, right=1, mask=True)
        res = [str(s.seq) for s in trim_edges(self._some_seqs())]
        assert res == ['aCCg', 'aAACCCGGg']

        trim_edges = TrimEdges(left=2, right=2, mask=True)
        res = [str(s.seq) for s in trim_edges(self._some_seqs())]
        assert res == ['accg', 'aaACCCGgg']

        trim_edges = TrimEdges(left=3, right=3, mask=True)
        res = [str(s.seq) for s in trim_edges(self._some_seqs())]
        assert res == ['accg', 'aaaCCCggg']

    def test_trim_edges_bin(self):
        'It tests the trim_edges binary'
        trim_bin = os.path.join(BIN_DIR, 'trim_edges')
        assert 'usage' in check_output([trim_bin, '-h'])

        fastq_fhand = _make_fhand(FASTQ2)
        result = check_output([trim_bin, fastq_fhand.name])
        assert '@seq1\nATCGT\n+' in result

        result = check_output([trim_bin, '-l', '1', '-r', '1',
                               fastq_fhand.name])
        assert '@seq1\nTCG\n+\n???\n' in result
        result = check_output([trim_bin, '-l', '1', '-r', '1', '-m',
                               fastq_fhand.name])
        assert '@seq1\naTCGt\n+\nA???A\n' in result


if __name__ == '__main__':
    #import sys;sys.argv = ['', 'SffExtractTest.test_items_in_gff']
    unittest.main()
