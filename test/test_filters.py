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
                            FilterBlastMatch, FilterDustComplexity,
                            seq_to_filterpackets, FilterByReadCount)
from crumbs.utils.bin_utils import BIN_DIR
from crumbs.utils.test_utils import TEST_DATA_DIR
from crumbs.utils.tags import NUCL, SEQS_FILTERED_OUT, SEQS_PASSED


class PacketConversionTest(unittest.TestCase):
    'It tests the seqs and filter packet conversion'
    def test_seqs_to_filter_packets(self):
        'It converts seq packets into filter packets'
        seqpackets = [['ACT'], ['CTG', 'TTT']]
        filter_packets = list(seq_to_filterpackets(iter(seqpackets)))
        assert [p[SEQS_PASSED] for p in filter_packets] == seqpackets
        assert [p[SEQS_FILTERED_OUT] for p in filter_packets] == [[], []]


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
        seqs = {SEQS_PASSED: [seq1, seq2], SEQS_FILTERED_OUT: []}

        filter_by_length = FilterByLength(minimum=4)
        passed = [str(s.seq) for s in filter_by_length(seqs)[SEQS_PASSED]]
        assert passed == ['aCTg']
        filtr = [str(s.seq) for s in filter_by_length(seqs)[SEQS_FILTERED_OUT]]
        assert filtr == ['AC']

        filter_by_length = FilterByLength(minimum=5)
        passed = [str(s.seq) for s in filter_by_length(seqs)[SEQS_PASSED]]
        assert passed == []
        filtr = [str(s.seq) for s in filter_by_length(seqs)[SEQS_FILTERED_OUT]]
        assert filtr == ['aCTg', 'AC']

        filter_by_length = FilterByLength(minimum=2, ignore_masked=True)
        passed = [str(s.seq) for s in filter_by_length(seqs)[SEQS_PASSED]]
        assert passed == ['aCTg', 'AC']
        filtr = [str(s.seq) for s in filter_by_length(seqs)[SEQS_FILTERED_OUT]]
        assert filtr == []

        filter_by_length = FilterByLength(minimum=3, ignore_masked=True)
        passed = [str(s.seq) for s in filter_by_length(seqs)[SEQS_PASSED]]
        assert passed == []
        filtr = [str(s.seq) for s in filter_by_length(seqs)[SEQS_FILTERED_OUT]]
        assert filtr == ['aCTg', 'AC']

        filter_by_length = FilterByLength(maximum=3, ignore_masked=True)
        passed = [str(s.seq) for s in filter_by_length(seqs)[SEQS_PASSED]]
        assert passed == ['aCTg', 'AC']
        filtr = [str(s.seq) for s in filter_by_length(seqs)[SEQS_FILTERED_OUT]]
        assert filtr == []

        filter_by_length = FilterByLength(maximum=3)
        passed = [str(s.seq) for s in filter_by_length(seqs)[SEQS_PASSED]]
        assert passed == ['AC']
        filtr = [str(s.seq) for s in filter_by_length(seqs)[SEQS_FILTERED_OUT]]
        assert filtr == ['aCTg']

        seq1 = _create_seqrecord('aCTTATg')
        seq2 = _create_seqrecord('ACCGCGC')
        filter_packet = {SEQS_PASSED: [seq1, seq2], SEQS_FILTERED_OUT: []}
        filter_by_length = FilterByLength(minimum=7, maximum=8,
                                          ignore_masked=True)
        assert len(filter_by_length(filter_packet)[SEQS_PASSED]) == 1

        filter_by_length = FilterByLength(minimum=7, maximum=8)
        assert len(filter_by_length(filter_packet)[SEQS_PASSED]) == 2

    def test_filter_by_length_bin(self):
        'It uses the filter_by_length binary'
        filter_bin = os.path.join(BIN_DIR, 'filter_by_length')
        assert 'usage' in check_output([filter_bin, '-h'])

        fasta = '>s1\naCTg\n>s2\nAC\n'
        fasta_fhand = _make_fhand(fasta)
        result = check_output([filter_bin, '-n', '4', fasta_fhand.name])
        assert '>s1\naCTg\n' in result

        result = check_output([filter_bin, '-x', '4', fasta_fhand.name])
        assert '>s2\nAC\n' in result

        result = check_output([filter_bin, '-mx', '4', fasta_fhand.name])
        assert '>s1\naCTg\n>s2\nAC\n' in result

        filtered_fhand = NamedTemporaryFile()
        result = check_output([filter_bin, '-mx', '1',
                               '-e', filtered_fhand.name, fasta_fhand.name])
        assert not result
        assert '>s1\naCTg\n>s2\nAC\n' in open(filtered_fhand.name).read()


class SeqListFilterTest(unittest.TestCase):
    'It tests the filtering using a list of sequences'
    def test_seq_list_filter(self):
        'It filters the reads given a list of ids'
        seq1 = SeqRecord(Seq('ACTG'), id='seq1')
        seq2 = SeqRecord(Seq('ACTG'), id='seq2')
        seqs = {SEQS_PASSED: [seq1, seq2], SEQS_FILTERED_OUT: []}
        ids = ['seq1']
        filter_by_id = FilterById(ids)
        passed = [str(s.id) for s in filter_by_id(seqs)[SEQS_PASSED]]
        assert passed == ['seq1']

        filter_by_id = FilterById(set(ids), reverse=True)
        filtered = [str(s.id) for s in filter_by_id(seqs)[SEQS_PASSED]]
        assert filtered == ['seq2']

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
        seqs = {SEQS_PASSED: [seq1, seq2], SEQS_FILTERED_OUT: []}

        filter_by_id = FilterByQuality(threshold=41)
        passed = [str(s.id) for s in filter_by_id(seqs)[SEQS_PASSED]]
        assert passed == ['seq1']

        filter_by_id = FilterByQuality(threshold=41, reverse=True)
        filtered = [str(s.id) for s in filter_by_id(seqs)[SEQS_PASSED]]
        assert filtered == ['seq2']

        filter_by_id = FilterByQuality(threshold=41.5, ignore_masked=True)
        passed = [str(s.id) for s in filter_by_id(seqs)[SEQS_PASSED]]
        assert passed == ['seq1']

    def test_filter_by_qual_bin(self):
        'It uses the filter_by_quality binary'
        filter_bin = os.path.join(BIN_DIR, 'filter_by_quality')
        assert 'usage' in check_output([filter_bin, '-h'])

        fastq = '@s1\naCTg\n+\n"DD"\n@s2\nAC\n+\n""\n'
        fastq_fhand = _make_fhand(fastq)
        filtered_fhand = NamedTemporaryFile()
        result = check_output([filter_bin, '-mq', '34', fastq_fhand.name,
                               '-e', filtered_fhand.name])
        assert '@s1\n' in result
        filtered = open(filtered_fhand.name).read()
        assert '@s2\n' in filtered

        # reverse
        filtered_fhand = NamedTemporaryFile()
        result = check_output([filter_bin, '-rmq', '34', fastq_fhand.name,
                               '-e', filtered_fhand.name])
        assert '@s2\n' in result
        filtered = open(filtered_fhand.name).read()
        assert '@s1\n' in filtered

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
        seq1 = SeqRecord(Seq(seq), id='seq')
        seqs = {SEQS_PASSED: [seq1], SEQS_FILTERED_OUT: []}
        filters = [{'kind': 'score_threshold', 'score_key': 'expect',
                    'max_score': 0.001},
                   {'kind': 'score_threshold', 'score_key': 'identity',
                    'min_score': 80},
                   {'kind': 'min_length', 'min_percentage': 60,
                    'length_in_query': True}]

        filter_ = FilterBlastMatch(blastdb, 'blastn', filters=filters,
                                   dbtype=NUCL)
        new_seqs = filter_(seqs)[SEQS_PASSED]
        assert  new_seqs == []

        filters = [{'kind': 'score_threshold', 'score_key': 'expect',
                    'max_score': 1e-28}]
        filter_ = FilterBlastMatch(blastdb, 'blastn', filters)
        new_seqs = filter_(seqs)[SEQS_PASSED]
        assert len(new_seqs) == 1

        filters = [{'kind': 'score_threshold', 'score_key': 'expect',
                    'max_score': 1e-28}]
        filter_ = FilterBlastMatch(blastdb, 'blastn', filters, reverse=True)
        filter_packets = filter_(seqs)
        assert  filter_packets[SEQS_PASSED] == []
        assert len(filter_packets[SEQS_FILTERED_OUT]) == 1

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
        result = check_output([filter_bin, '-b', blastdb, '-x', '1e-27',
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
        filtered_fhand = NamedTemporaryFile()
        result = check_output([filter_bin, '-b', blastdb, '-l', '40',
                               seq_fhand.name, '-e', filtered_fhand.name])
        assert ">seq\nATCATGTAGTTACACATG" in open(filtered_fhand.name).read()
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
        filtered_fhand = NamedTemporaryFile()
        result = check_output([filter_bin, '-b', blastdb, '-a', '80', '-r',
                               seq_fhand.name, '-e', filtered_fhand.name])
        assert result == ''
        assert ">seq\nATCATGTAGTTACACATG" in open(filtered_fhand.name).read()


class ComplexityFilterTest(unittest.TestCase):
    'It tests the filtering by complexity'
    @staticmethod
    def test_dust_filter():
        'It tests the complexity filter'
        seq1 = 'TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTAAAAAAAAAAAAAAAAAAAAAAAAA'
        seq2 = 'CATCGATTGCGATCGATCTTGTTGCACGACTAGCTATCGATTGCTAGCTTAGCTAGCTAGTT'
        seqs = {SEQS_PASSED: [SeqRecord(Seq(seq1), id='seq1'),
                              SeqRecord(Seq(seq2), id='seq2')],
                SEQS_FILTERED_OUT: []}

        filter_dust = FilterDustComplexity()
        filter_packet = filter_dust(seqs)
        assert len(filter_packet[SEQS_PASSED]) == 1
        assert len(filter_packet[SEQS_FILTERED_OUT]) == 1

        assert filter_packet[SEQS_PASSED][0].id == 'seq2'
        assert filter_packet[SEQS_FILTERED_OUT][0].id == 'seq1'

        # reverse
        filter_dust = FilterDustComplexity(reverse=True)
        filter_packet = filter_dust(seqs)
        assert len(filter_packet[SEQS_PASSED]) == 1
        assert len(filter_packet[SEQS_FILTERED_OUT]) == 1

        assert filter_packet[SEQS_PASSED][0].id == 'seq1'
        assert filter_packet[SEQS_FILTERED_OUT][0].id == 'seq2'

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

        filtered_fhand = NamedTemporaryFile()
        result = check_output([filter_bin, fasta_fhand.name,
                               '-e', filtered_fhand.name])
        assert '>s1\n' in open(filtered_fhand.name).read()
        assert '>s2\n' in result


class ReadCountFIlterTest(unittest.TestCase):
    def test_filter_by_read_count(self):
        seq1 = 'TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTAAAAAAAAAAAAAAAAAAAAAAAAA'
        seq2 = 'CATCGATTGCGATCGATCTTGTTGCACGACTAGCTATCGATTGCTAGCTTAGCTAGCTAGTT'
        seqs = {SEQS_PASSED: [SeqRecord(Seq(seq1), id='seq1'),
                              SeqRecord(Seq(seq2), id='seq2')],
                SEQS_FILTERED_OUT: []}
        read_counts = {'seq1': {'mapped_reads': 1000000000,
                                'unmapped_reads': 1000000000,
                                'length': len(seq1)},
                       'seq2': {'mapped_reads': 100, 'unmapped_reads': 100,
                                'length': len(seq1)}}
        filter_ = FilterByReadCount(read_counts, 1000)
        seqs2 = filter_(seqs)
        assert seqs2[SEQS_FILTERED_OUT][0].id == 'seq2'
        assert seqs2[SEQS_PASSED][0].id == 'seq1'

        filter_ = FilterByReadCount(read_counts, 1000, reverse=True)
        seqs2 = filter_(seqs)
        assert seqs2[SEQS_FILTERED_OUT][0].id == 'seq1'
        assert seqs2[SEQS_PASSED][0].id == 'seq2'


if __name__ == "__main__":
    # import sys;sys.argv = ['', 'ComplexityFilterTest']
    unittest.main()
