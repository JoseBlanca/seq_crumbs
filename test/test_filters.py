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
from subprocess import check_output, call, CalledProcessError
import os.path
from tempfile import NamedTemporaryFile

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from crumbs.filters import (FilterByLength, FilterById, FilterByQuality,
                            FilterBlastMatch, _calculate_dust_score,
                            FilterDustComplexity)
from crumbs.utils.bin_utils import BIN_DIR
from crumbs.utils.test_utils import TEST_DATA_DIR
from crumbs.utils.tags import NUCL


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

        filter_by_length = FilterByLength(minimum=4)
        assert [str(s.seq) for s in filter_by_length(seqs)] == ['aCTg']

        filter_by_length = FilterByLength(minimum=5)
        assert [str(s.seq) for s in filter_by_length(seqs)] == []

        filter_by_length = FilterByLength(minimum=2, ignore_masked=True)
        assert [str(s.seq) for s in filter_by_length(seqs)] == ['aCTg', 'AC']

        filter_by_length = FilterByLength(minimum=3, ignore_masked=True)
        assert [str(s.seq) for s in filter_by_length(seqs)] == []

        filter_by_length = FilterByLength(maximum=3, ignore_masked=True)
        assert [str(s.seq) for s in filter_by_length(seqs)] == ['aCTg', 'AC']

        filter_by_length = FilterByLength(maximum=3)
        assert [str(s.seq) for s in filter_by_length(seqs)] == ['AC']

        seq1 = _create_seqrecord('aCTTATg')
        seq2 = _create_seqrecord('ACCGCGC')
        seqs = [seq1, seq2]
        filter_by_length = FilterByLength(minimum=7, maximum=8,
                                          ignore_masked=True)
        assert len([str(s.seq) for s in filter_by_length(seqs)]) == 1

        filter_by_length = FilterByLength(minimum=7, maximum=8)
        assert len([str(s.seq) for s in filter_by_length(seqs)]) == 2

    def test_filter_by_length_bin(self):
        'It uses the filter_by_length binary'
        filter_bin = os.path.join(BIN_DIR, 'filter_by_length')
        assert 'usage' in check_output([filter_bin, '-h'])

        fasta = '>s1\naCTg\n>s2\nAC\n'
        fasta_fhand = _make_fhand(fasta)
        result = check_output([filter_bin, '-n', '4', fasta_fhand.name])
        assert '>s1\naCTg\n' in result

        result = check_output([filter_bin, '-x', '4', '-r', fasta_fhand.name])
        assert '>s2\nAC\n' in result

        result = check_output([filter_bin, '-rmx', '4', fasta_fhand.name])
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


class BlastMatchFilterTest(unittest.TestCase):
    'It tests the filtering by blast match'
    @staticmethod
    def test_blastmatch_filter():
        'it test filter by blast'
        blastdb = os.path.join(TEST_DATA_DIR, 'blastdbs', 'arabidopsis_genes')

        match = 'CCAAAGTACGGTCTCCCAAGCGGTCTCTTACCGGACACCGTCACCGATTTCACCCTCT'
        seq = 'ATCATGTAGTTACACATGAACACACACATG'
        seq += match

        seqs = [SeqRecord(Seq(seq), id='seq')]
        filters = [{'kind': 'score_threshold', 'score_key': 'expect',
                    'max_score': 0.001},
                   {'kind': 'score_threshold', 'score_key': 'identity',
                    'min_score': 80},
                   {'kind': 'min_length', 'min_percentage': 60,
                    'length_in_query': True}]

        filter_ = FilterBlastMatch(blastdb, 'blastn', filters=filters,
                                   dbtype=NUCL)
        new_seqs = filter_(seqs)
        assert  new_seqs == []

        filters = [{'kind': 'score_threshold', 'score_key': 'expect',
                    'max_score': 1e-28}]
        filter_ = FilterBlastMatch(blastdb, 'blastn', filters)
        new_seqs = filter_(seqs)
        assert len(new_seqs) == 1

        filters = [{'kind': 'score_threshold', 'score_key': 'expect',
                    'max_score': 1e-28}]
        filter_ = FilterBlastMatch(blastdb, 'blastn', filters, reverse=True)
        new_seqs = filter_(seqs)
        assert new_seqs == []

    def test_filter_blast_bin(self):
        'It test the binary of the filter_by_blast'
        filter_bin = os.path.join(BIN_DIR, 'filter_by_blast')
        assert 'usage' in check_output([filter_bin, '-h'])
        blastdb = os.path.join(TEST_DATA_DIR, 'blastdbs', 'arabidopsis_genes')

        match = 'CCAAAGTACGGTCTCACAAGCGGTCTCTTACCGGACACCGTCACCGATTTCACCCTCT'
        seq = 'ATCATGTAGTTACACATGAACACACACATG'
        seq += match
        seq_fhand = _make_fhand('>seq\n' + seq + '\n')
        # fail if no filters
        stderr = NamedTemporaryFile()
        try:
            check_output([filter_bin, '-b', blastdb, seq_fhand.name],
                         stderr=stderr)
            self.fail()
        except CalledProcessError:
            assert 'At least a filter must be Used' in open(stderr.name).read()

        # expected_filter
        result = check_output([filter_bin, '-b', blastdb, '-e', '1e-27',
                               seq_fhand.name])

        assert 'CATGAACACACACAT' in result

        # similarity_filter
        result = check_output([filter_bin, '-b', blastdb, '-s', '99',
                               seq_fhand.name])
        assert 'CATGAACACACACAT' in result

        # fail if -a an -l given i in the command line
        stderr = NamedTemporaryFile()
        try:
            check_output([filter_bin, '-b', blastdb, seq_fhand.name, '-a', '3',
                          '-l', '4'], stderr=stderr)
            self.fail()
        except CalledProcessError:
            assert 'not allowed with argument'  in open(stderr.name).read()

        # minlen percentaje_filter
        result = check_output([filter_bin, '-b', blastdb, '-l', '40',
                               seq_fhand.name])
        assert result == ''

        result = check_output([filter_bin, '-b', blastdb, '-l', '80',
                               seq_fhand.name])
        assert 'CATGAACACACACAT' in result

        # min nucleotides in query filter
        result = check_output([filter_bin, '-b', blastdb, '-a', '40',
                               seq_fhand.name])
        assert result == ''

        result = check_output([filter_bin, '-b', blastdb, '-a', '80',
                               seq_fhand.name])
        assert 'CATGAACACACACAT' in result

        # reverse
        result = check_output([filter_bin, '-b', blastdb, '-a', '80', '-r',
                               seq_fhand.name])
        assert result == ''


class ComplexityFilterTest(unittest.TestCase):
    'It tests the filtering by complexity'
    @staticmethod
    def test_dustscore_calculation():
        'It calculates the dust score'
        seqs = ['TTTTTTTTTTTTTTTTTTTTTTTTTTTT', 'TATATATATATATATATATATATATATA',
                'GAAGAAGAAGAAGAAGAAGAAGAAGAAG', 'AACTGCAGTCGATGCTGATTCGATCGAT',
                'AACTGAAAAAAAATTTTTTTAAAAAAAA']

        # short sequences
        scores = [100, 48, 30.76, 4.31, 23.38]
        scoresx3 = [100, 48.68, 28.65, 5.62, 27.53]
        scoresx4 = [100, 48.55, 28.25, 5.79, 28.00]
        for seq, score, scorex3, scorex4 in zip(seqs, scores, scoresx3,
                                                scoresx4):
            seqrec = SeqRecord(Seq(seq))
            assert _calculate_dust_score(seqrec) - score < 0.01
            seqrec = SeqRecord(Seq(seq * 3))
            assert _calculate_dust_score(seqrec) - scorex3 < 0.01
            seqrec = SeqRecord(Seq(seq * 4))
            assert _calculate_dust_score(seqrec) - scorex4 < 0.01

    @staticmethod
    def test_dust_filter():
        'It tests the complexity filter'
        seq1 = 'TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTAAAAAAAAAAAAAAAAAAAAAAAAA'
        seq2 = 'CATCGATTGCGATCGATCTTGTTGCACGACTAGCTATCGATTGCTAGCTTAGCTAGCTAGTT'
        seqs = [SeqRecord(Seq(seq1), id='seq1'),
                SeqRecord(Seq(seq2), id='seq2')]
        filter_dust = FilterDustComplexity()
        new_seqs = filter_dust(seqs)
        assert len(new_seqs) == 1
        assert new_seqs[0].id == 'seq2'

    @staticmethod
    def test_filter_by_dust_bin():
        'It uses the filter_by_complexity binary'
        filter_bin = os.path.join(BIN_DIR, 'filter_by_complexity')
        assert 'usage' in check_output([filter_bin, '-h'])

        fasta = '>s1\nTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n'
        fasta += '>s2\nCATCGATTGCGATCGATCTTGTTGCACGACTAGCTATCGATTGCTAGCTTAGT\n'
        fasta_fhand = _make_fhand(fasta)
        result = check_output([filter_bin, fasta_fhand.name])
        assert '>s1\n' not in result
        assert '>s2\n' in result

        result = check_output([filter_bin, '-r', fasta_fhand.name])
        assert '>s1\n' in result
        assert '>s2\n' not in result

        result = check_output([filter_bin, '-c', '1', fasta_fhand.name])
        assert result == ''


if __name__ == "__main__":
    # import sys;sys.argv = ['', 'BlastMatchFilterTest.test_blastmatch_filter']
    # import sys;sys.argv = ['', 'BlastMatchFilterTest.test_blastmatch_bin']
    unittest.main()
