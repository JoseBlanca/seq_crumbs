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

# pylint: disable=R0201
# pylint: disable=R0904

import unittest
# pylint: disable=W0402
from string import ascii_lowercase
from random import choice
from subprocess import check_output, call
import os.path
from tempfile import NamedTemporaryFile

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from crumbs.filters import FilterByLength, FilterById, FilterByQuality
from crumbs.utils.bin_utils import BIN_DIR


def _create_seqrecord(string):
    'Given an string it returns a SeqRecord'
    # pylint: disable=W0612
    return SeqRecord(Seq(string),
                     id=''.join([choice(ascii_lowercase) for i in range(6)]))


def _make_fhand(content=''):
    'It makes temporary fhands'
    fhand = NamedTemporaryFile()
    fhand.write(content)
    fhand.flush()
    return fhand


class LengthFilterTest(unittest.TestCase):
    'It tests the filtering according to the sequence length'
    def test_length_filter(self):
        'It filters the seqs according to its length'
        seq1 = _create_seqrecord('aCTg')
        seq2 = _create_seqrecord('AC')
        seqs = [seq1, seq2]

        filter_by_length = FilterByLength(threshold=4)
        assert [str(s.seq) for s in filter_by_length(seqs)] == ['aCTg']

        filter_by_length = FilterByLength(threshold=4, reverse=True)
        assert [str(s.seq) for s in filter_by_length(seqs)] == ['AC']

        filter_by_length = FilterByLength(threshold=5)
        assert [str(s.seq) for s in filter_by_length(seqs)] == []

        filter_by_length = FilterByLength(threshold=5, reverse=True)
        assert [str(s.seq) for s in filter_by_length(seqs)] == ['aCTg', 'AC']

        filter_by_length = FilterByLength(threshold=2, ignore_masked=True)
        assert [str(s.seq) for s in filter_by_length(seqs)] == ['aCTg', 'AC']

        filter_by_length = FilterByLength(threshold=3, ignore_masked=True)
        assert [str(s.seq) for s in filter_by_length(seqs)] == []

        filter_by_length = FilterByLength(threshold=3, ignore_masked=True,
                                          reverse=True)
        assert [str(s.seq) for s in filter_by_length(seqs)] == ['aCTg', 'AC']

    def test_filter_by_length_bin(self):
        'It uses the filter_by_length binary'
        filter_bin = os.path.join(BIN_DIR, 'filter_by_length')
        assert 'usage' in check_output([filter_bin, '-h'])

        fasta = '>s1\naCTg\n>s2\nAC\n'
        fasta_fhand = _make_fhand(fasta)
        result = check_output([filter_bin, '-l', '4', fasta_fhand.name])
        assert '>s1\naCTg\n' in result

        result = check_output([filter_bin, '-l', '4', '-r', fasta_fhand.name])
        assert '>s2\nAC\n' in result

        result = check_output([filter_bin, '-rml', '4', fasta_fhand.name])
        assert '>s1\naCTg\n>s2\nAC\n' in result


class SeqListFilterTest(unittest.TestCase):
    'It tests the filtering using a list of sequences'
    def test_seq_list_filter(self):
        'It filters the reads given a list of ids'
        seq1 = SeqRecord(Seq('ACTG'), id='seq1')
        seq2 = SeqRecord(Seq('ACTG'), id='seq2')
        seqs = [seq1, seq2]

        ids = ['seq1']
        filter_by_id = FilterById(ids)
        assert [s.id for s in filter_by_id(seqs)] == ['seq1']

        filter_by_id = FilterById(set(ids), reverse=True)
        assert [s.id for s in filter_by_id(seqs)] == ['seq2']

    def test_filter_by_name_bin(self):
        'It uses the filter_by_name binary'
        filter_bin = os.path.join(BIN_DIR, 'filter_by_name')
        assert 'usage' in check_output([filter_bin, '-h'])

        fasta = '>s1\naCTg\n>s2\nAC\n'
        fasta_fhand = _make_fhand(fasta)
        list_fhand = _make_fhand('s1\n')
        result = check_output([filter_bin, '-rl', list_fhand.name,
                               fasta_fhand.name])
        assert '>s2\n' in result
        assert '>s1\n' not in result


class QualityFilterTest(unittest.TestCase):
    'It tests the filtering by a quality threshold'
    def test_quality_filter(self):
        'It filters the reads given a quality threshold'
        seq1 = SeqRecord(Seq('AAcTg'), id='seq1',
                    letter_annotations={'phred_quality': [42, 42, 40, 42, 40]})
        seq2 = SeqRecord(Seq('AAcTg'), id='seq2',
                    letter_annotations={'phred_quality': [40, 40, 42, 40, 42]})
        seqs = [seq1, seq2]

        filter_by_id = FilterByQuality(threshold=41)
        assert [s.id for s in filter_by_id(seqs)] == ['seq1']

        filter_by_id = FilterByQuality(threshold=41, reverse=True)
        assert [s.id for s in filter_by_id(seqs)] == ['seq2']

        filter_by_id = FilterByQuality(threshold=41.5, ignore_masked=True)
        assert [s.id for s in filter_by_id(seqs)] == ['seq1']

    def test_filter_by_qual_bin(self):
        'It uses the filter_by_quality binary'
        filter_bin = os.path.join(BIN_DIR, 'filter_by_quality')
        assert 'usage' in check_output([filter_bin, '-h'])

        fastq = '@s1\naCTg\n+\n"DD"\n@s2\nAC\n+\n""\n'
        fastq_fhand = _make_fhand(fastq)
        result = check_output([filter_bin, '-mrq', '34', fastq_fhand.name])
        assert '@s1\n' not in result
        assert '@s2\n' in result

        # Using a fasta file it will fail
        fasta = '>s1\naCTg\n>s2\nAC\n'
        fasta_fhand = _make_fhand(fasta)
        stderr = NamedTemporaryFile()
        result = call([filter_bin, '-q', '35', fasta_fhand.name],
                      stderr=stderr)
        assert result
        assert 'sequences do not have qualities' in open(stderr.name).read()


if __name__ == "__main__":
#    import sys;sys.argv = ['', 'TestPool']
    unittest.main()
