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
# pylint: disable=W0402
# pylint: disable=C0111

import unittest

from cStringIO import StringIO
from string import ascii_lowercase
from random import choice
from subprocess import check_output, call, CalledProcessError
import os.path
from tempfile import NamedTemporaryFile

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation

from crumbs.filters import (FilterByLength, FilterById, FilterByQuality,
                            FilterBlastMatch, FilterDustComplexity,
                            seq_to_filterpackets, FilterByRpkm, FilterByBam,
                            FilterBowtie2Match, FilterByFeatureTypes,
    FilterDuplicates)
from crumbs.utils.bin_utils import BIN_DIR
from crumbs.utils.test_utils import TEST_DATA_DIR
from crumbs.utils.tags import (NUCL, SEQS_FILTERED_OUT, SEQS_PASSED, SEQITEM,
                               SEQRECORD)
from crumbs.utils.file_utils import TemporaryDir
from crumbs.seq import get_name, get_str_seq, SeqWrapper
from crumbs.mapping import get_or_create_bowtie2_index
from crumbs.seqio import read_seq_packets


_seqs_to_names = lambda seqs: [get_name(s) for pair in seqs for s in pair]
_seqs_to_str_seqs = lambda seqs: [get_str_seq(s) for pai in seqs for s in pai]


class PacketConversionTest(unittest.TestCase):
    'It tests the seqs and filter packet conversion'
    def test_seqs_to_filter_packets(self):
        'It converts seq packets into filter packets'
        seqpackets = [['ACT'], ['CTG', 'TTT']]
        filter_packets = list(seq_to_filterpackets(iter(seqpackets)))
        expected = [[['ACT']], [['CTG'], ['TTT']]]
        assert [p[SEQS_PASSED] for p in filter_packets] == expected
        assert [p[SEQS_FILTERED_OUT] for p in filter_packets] == [[], []]


def _create_seqrecord(string):
    'Given an string it returns a SeqRecord'
    # pylint: disable=W0612
    seq = SeqRecord(Seq(string),
                     id=''.join([choice(ascii_lowercase) for i in range(6)]))
    return SeqWrapper(kind=SEQRECORD, object=seq, file_format=None)


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
        seqs = {SEQS_PASSED: [[seq1], [seq2]], SEQS_FILTERED_OUT: []}

        filter_by_length = FilterByLength(minimum=4)
        passed = _seqs_to_str_seqs(filter_by_length(seqs)[SEQS_PASSED])
        assert passed == ['aCTg']
        filtr = _seqs_to_str_seqs(filter_by_length(seqs)[SEQS_FILTERED_OUT])
        assert filtr == ['AC']

        filter_by_length = FilterByLength(minimum=5)
        passed = _seqs_to_str_seqs(filter_by_length(seqs)[SEQS_PASSED])
        assert passed == []
        filtr = _seqs_to_str_seqs(filter_by_length(seqs)[SEQS_FILTERED_OUT])
        assert filtr == ['aCTg', 'AC']

        filter_by_length = FilterByLength(minimum=2, ignore_masked=True)
        passed = _seqs_to_str_seqs(filter_by_length(seqs)[SEQS_PASSED])
        assert passed == ['aCTg', 'AC']
        filtr = _seqs_to_str_seqs(filter_by_length(seqs)[SEQS_FILTERED_OUT])
        assert filtr == []

        filter_by_length = FilterByLength(minimum=3, ignore_masked=True)
        passed = _seqs_to_str_seqs(filter_by_length(seqs)[SEQS_PASSED])
        assert passed == []
        filtr = _seqs_to_str_seqs(filter_by_length(seqs)[SEQS_FILTERED_OUT])
        assert filtr == ['aCTg', 'AC']

        filter_by_length = FilterByLength(maximum=3, ignore_masked=True)
        passed = _seqs_to_str_seqs(filter_by_length(seqs)[SEQS_PASSED])
        assert passed == ['aCTg', 'AC']
        filtr = _seqs_to_str_seqs(filter_by_length(seqs)[SEQS_FILTERED_OUT])
        assert filtr == []

        filter_by_length = FilterByLength(maximum=3)
        passed = _seqs_to_str_seqs(filter_by_length(seqs)[SEQS_PASSED])
        assert passed == ['AC']
        filtr = _seqs_to_str_seqs(filter_by_length(seqs)[SEQS_FILTERED_OUT])
        assert filtr == ['aCTg']

        seq1 = _create_seqrecord('aCTTATg')
        seq2 = _create_seqrecord('ACCGCGC')
        filter_packet = {SEQS_PASSED: [[seq1], [seq2]], SEQS_FILTERED_OUT: []}
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

        # with pairs
        fasta = '>s1.f\naCTg\n>s1.r\nAC\n>s2.f\naTg\n>s2.r\nAC\n'
        fasta += '>s3.f\naCTg\n>s3.r\nACTG\n'
        fasta_fhand = _make_fhand(fasta)
        stderr = NamedTemporaryFile(suffix='.stderr')
        result = check_output([filter_bin, '-n', '4', '--paired_reads',
                               fasta_fhand.name], stderr=stderr)
        # By default fail drags pair
        assert '>s3.f\naCTg\n>s3.r\nACTG\n' in result

        result = check_output([filter_bin, '-n', '4', '--paired_reads',
                               '--fail_drags_pair', 'true', fasta_fhand.name])
        assert '>s3.f\naCTg\n>s3.r\nACTG\n' in result

        result = check_output([filter_bin, '-n', '4', '--paired_reads',
                               '--fail_drags_pair', 'false', fasta_fhand.name])
        assert '>s1.f\naCTg\n>s1.r\nAC\n>s3.f\naCTg\n>s3.r\nACTG\n' in result


class SeqListFilterTest(unittest.TestCase):
    'It tests the filtering using a list of sequences'
    def test_seq_list_filter(self):
        'It filters the reads given a list of ids'
        seq1 = SeqRecord(Seq('ACTG'), id='seq1')
        seq1 = SeqWrapper(object=seq1, kind=SEQRECORD, file_format=None)
        seq2 = SeqRecord(Seq('ACTG'), id='seq2')
        seq2 = SeqWrapper(object=seq2, kind=SEQRECORD, file_format=None)
        seqs = {SEQS_PASSED: [[seq1], [seq2]], SEQS_FILTERED_OUT: []}
        ids = ['seq1']
        filter_by_id = FilterById(ids)

        passed = _seqs_to_names(filter_by_id(seqs)[SEQS_PASSED])
        assert passed == ['seq1']

        filter_by_id = FilterById(set(ids), reverse=True)
        passed = _seqs_to_names(filter_by_id(seqs)[SEQS_PASSED])
        assert passed == ['seq2']

    def test_with_pairs(self):
        seq1 = SeqRecord(Seq('ACTG'), id='seq1')
        seq1 = SeqWrapper(object=seq1, kind=SEQRECORD, file_format=None)
        seq2 = SeqRecord(Seq('ACTG'), id='seq2')
        seq2 = SeqWrapper(object=seq2, kind=SEQRECORD, file_format=None)
        seqs = {SEQS_PASSED: [[seq1, seq2]], SEQS_FILTERED_OUT: []}
        ids = ['seq1']

        filter_by_id = FilterById(ids, failed_drags_pair=True)
        passed = _seqs_to_names(filter_by_id(seqs)[SEQS_PASSED])
        assert not passed

        filter_by_id = FilterById(ids, failed_drags_pair=False)
        passed = _seqs_to_names(filter_by_id(seqs)[SEQS_PASSED])
        assert passed == ['seq1', 'seq2']

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
        seq1 = SeqWrapper(object=seq1, kind=SEQRECORD, file_format=None)
        seq2 = SeqRecord(Seq('AAcTg'), id='seq2',
                    letter_annotations={'phred_quality': [40, 40, 42, 40, 42]})
        seq2 = SeqWrapper(object=seq2, kind=SEQRECORD, file_format=None)
        seqs = {SEQS_PASSED: [[seq1], [seq2]], SEQS_FILTERED_OUT: []}

        filter_ = FilterByQuality(threshold=41)
        passed = _seqs_to_names(filter_(seqs)[SEQS_PASSED])
        assert passed == ['seq1']

        filter_ = FilterByQuality(threshold=41, reverse=True)
        passed = _seqs_to_names(filter_(seqs)[SEQS_PASSED])
        assert passed == ['seq2']

        filter_ = FilterByQuality(threshold=41.5, ignore_masked=True)
        passed = _seqs_to_names(filter_(seqs)[SEQS_PASSED])
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
        seq1 = SeqWrapper(object=seq1, kind=SEQRECORD, file_format=None)
        seqs = {SEQS_PASSED: [[seq1]], SEQS_FILTERED_OUT: []}
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

    def test_filter_blast_bin_fastq(self):
        filter_bin = os.path.join(BIN_DIR, 'filter_by_blast')
        blastdb = os.path.join(TEST_DATA_DIR, 'blastdbs', 'arabidopsis_genes')
        # With fastq
        match = 'CCAAAGTACGGTCTCACAAGCGGTCTCTTACCGGACACCGTCACCGATTTCACCCTCT'
        seq = 'ATCATGTAGTTACACATGAACACACACATG'
        seq += match

        fastq = '@seq\n' + seq + '\n+\n' + 'a' * len(seq) + '\n'
        seq_fhand = _make_fhand(fastq + fastq)

        result = check_output([filter_bin, '-b', blastdb, '-x', '1e-27',
                               seq_fhand.name])

        assert 'CATGAACACACACAT' in result


class ComplexityFilterTest(unittest.TestCase):
    'It tests the filtering by complexity'
    @staticmethod
    def test_dust_filter():
        'It tests the complexity filter'
        seq1 = 'TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTAAAAAAAAAAAAAAAAAAAAAAAAA'
        seq2 = 'CATCGATTGCGATCGATCTTGTTGCACGACTAGCTATCGATTGCTAGCTTAGCTAGCTAGTT'
        seq1 = SeqRecord(Seq(seq1), id='seq1')
        seq2 = SeqRecord(Seq(seq2), id='seq2')
        seq1 = SeqWrapper(SEQRECORD, seq1, None)
        seq2 = SeqWrapper(SEQRECORD, seq2, None)
        seqs = {SEQS_PASSED: [[seq1], [seq2]], SEQS_FILTERED_OUT: []}

        filter_dust = FilterDustComplexity()
        filter_packet = filter_dust(seqs)
        assert len(filter_packet[SEQS_PASSED]) == 1
        assert len(filter_packet[SEQS_FILTERED_OUT]) == 1

        assert _seqs_to_names(filter_packet[SEQS_PASSED])[0] == 'seq2'
        assert _seqs_to_names(filter_packet[SEQS_FILTERED_OUT])[0] == 'seq1'

        # reverse
        filter_dust = FilterDustComplexity(reverse=True)
        filter_packet = filter_dust(seqs)
        assert len(filter_packet[SEQS_PASSED]) == 1
        assert len(filter_packet[SEQS_FILTERED_OUT]) == 1

        assert _seqs_to_names(filter_packet[SEQS_PASSED])[0] == 'seq1'
        assert _seqs_to_names(filter_packet[SEQS_FILTERED_OUT])[0] == 'seq2'

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


class RpkmFilterTest(unittest.TestCase):
    def test_filter_by_read_count(self):
        seq1 = 'T' * 1000
        seq2 = 'A' * 1000
        seq1 = SeqRecord(Seq(seq1), id='seq1')
        seq2 = SeqRecord(Seq(seq2), id='seq2')
        seq1 = SeqWrapper(SEQRECORD, seq1, None)
        seq2 = SeqWrapper(SEQRECORD, seq2, None)
        seqs = {SEQS_PASSED: [[seq1], [seq2]], SEQS_FILTERED_OUT: []}

        read_counts = {'seq1': {'mapped_reads': 10,
                                'unmapped_reads': 999989,
                                'length': len(seq1.object)},
                       'seq2': {'mapped_reads': 1, 'unmapped_reads': 0,
                                'length': len(seq2.object)}}
        filter_ = FilterByRpkm(read_counts, 2)
        seqs2 = filter_(seqs)
        assert _seqs_to_names(seqs2[SEQS_FILTERED_OUT]) == ['seq2']
        assert _seqs_to_names(seqs2[SEQS_PASSED]) == ['seq1']

        filter_ = FilterByRpkm(read_counts, 1)
        seqs2 = filter_(seqs)
        assert not seqs2[SEQS_FILTERED_OUT]

        filter_ = FilterByRpkm(read_counts, 2, reverse=True)
        seqs2 = filter_(seqs)
        assert _seqs_to_names(seqs2[SEQS_FILTERED_OUT]) == ['seq1']
        assert _seqs_to_names(seqs2[SEQS_PASSED]) == ['seq2']


class BamFilterTest(unittest.TestCase):
    'It tests the filtering of the mapped reads in a BAM file'
    @staticmethod
    def test_bam_filter():
        'it test filter by being mapped in a BAM file'
        reads = [SeqRecord(seq=Seq('aaa'), id='seq{}'.format(n)) for n in range(16, 23)]
        reads = [[SeqWrapper(SEQRECORD, r, None)] for r in reads]
        bam_fpath = os.path.join(TEST_DATA_DIR, 'seqs.bam')
        filter_ = FilterByBam([bam_fpath])
        filterpacket = {SEQS_PASSED: reads, SEQS_FILTERED_OUT: []}
        new_filterpacket = filter_(filterpacket)
        passed = _seqs_to_names(new_filterpacket[SEQS_PASSED])
        assert  passed == ['seq16', 'seq17', 'seq18']
        filtered_out = _seqs_to_names(new_filterpacket[SEQS_FILTERED_OUT])
        assert filtered_out == ['seq19', 'seq20', 'seq21', 'seq22']


class FilterBowtie2Test(unittest.TestCase):
    @staticmethod
    def test_filter_by_bowtie2():
        directory = TemporaryDir()
        index_fpath = get_or_create_bowtie2_index(os.path.join(TEST_DATA_DIR,
                                                          'arabidopsis_genes'),
                                                  directory=directory.name)
        fastq_fpath = os.path.join(TEST_DATA_DIR, 'arabidopsis_reads.fastq')
        fasta_fpath = os.path.join(TEST_DATA_DIR, 'arabidopsis_reads.fasta')

        passed = ['no_arabi']
        for preffered_classes in [[SEQITEM], [SEQRECORD]]:
            for reads_fpath in [fastq_fpath, fasta_fpath]:
                seq_packets = read_seq_packets([open(reads_fpath)],
                                        prefered_seq_classes=preffered_classes)
                filter_packets = seq_to_filterpackets(seq_packets)
                filter_ = FilterBowtie2Match(index_fpath)
                filter_packet = list(filter_packets)[0]
                filter_packets = filter_(filter_packet)
                assert _seqs_to_names(filter_packets[SEQS_PASSED]) == passed
                assert _seqs_to_names(filter_packets[SEQS_FILTERED_OUT]) == [
                                                     'read1', 'read2', 'read3']
        directory.close()

    @staticmethod
    def test_filter_by_bowtie2_bin():
        filter_bin = os.path.join(BIN_DIR, 'filter_by_bowtie2')
        assert 'usage' in check_output([filter_bin, '-h'])
        directory = TemporaryDir()
        index_fpath = get_or_create_bowtie2_index(os.path.join(TEST_DATA_DIR,
                                                          'arabidopsis_genes'),
                                                  directory=directory.name)

        fastq_fpath = os.path.join(TEST_DATA_DIR, 'arabidopsis_reads.fastq')
        fasta_fpath = os.path.join(TEST_DATA_DIR, 'arabidopsis_reads.fasta')
        for reads_fpath in [fastq_fpath, fasta_fpath]:
            out_fhand = NamedTemporaryFile(suffix='.seqs')
            filtered_fhand = NamedTemporaryFile(suffix='.seqs')
            cmd = [filter_bin, '-i', index_fpath, '-o', out_fhand.name,
                   '-e', filtered_fhand.name, reads_fpath]
            check_output(cmd)
            assert 'no_arabi' in open(out_fhand.name).read()
            assert 'read1' in open(filtered_fhand.name).read()
        directory.close()


class FilterByFeatureTypeTest(unittest.TestCase):
    def test_filter_by_feat_type(self):
        orf = SeqFeature(FeatureLocation(3, 4), type='ORF')
        seq1 = SeqRecord(Seq('aaaa'), id='seq1', features=[orf])
        seq2 = SeqRecord(Seq('aaaa'), id='seq2')
        seq1 = SeqWrapper(SEQRECORD, seq1, None)
        seq2 = SeqWrapper(SEQRECORD, seq2, None)
        seqs = {SEQS_PASSED: [[seq1], [seq2]], SEQS_FILTERED_OUT: []}

        filter_ = FilterByFeatureTypes(['ORF'])
        seqs = filter_(seqs)
        assert len(seqs[SEQS_FILTERED_OUT]) == 1
        assert len(seqs[SEQS_PASSED]) == 1


class FilterByDuplicity(unittest.TestCase):
    def test_filter_by_duplicity(self):

        seq1 = SeqRecord(Seq('aaaa'), id='seq1')
        seq2 = SeqRecord(Seq('aaaa'), id='seq2')
        seq3 = SeqRecord(Seq('aaab'), id='seq2')
        seq1 = SeqWrapper(SEQRECORD, seq1, None)
        seq2 = SeqWrapper(SEQRECORD, seq2, None)
        seq3 = SeqWrapper(SEQRECORD, seq3, None)
        seqs = {SEQS_PASSED: [[seq1], [seq2], [seq3]], SEQS_FILTERED_OUT: []}

        filter_ = FilterDuplicates()
        seqs = filter_(seqs)
        assert len(seqs[SEQS_FILTERED_OUT]) == 1
        assert len(seqs[SEQS_PASSED]) == 2

        seq1 = SeqRecord(Seq('aaaa'), id='seq1.f')
        seq2 = SeqRecord(Seq('aaaa'), id='seq1.r')
        seq3 = SeqRecord(Seq('aaaa'), id='seq2.f')
        seq4 = SeqRecord(Seq('aaab'), id='seq2.r')
        seq5 = SeqRecord(Seq('aaaa'), id='seq3.f')
        seq6 = SeqRecord(Seq('aaab'), id='seq3.r')

        seq1 = SeqWrapper(SEQRECORD, seq1, None)
        seq2 = SeqWrapper(SEQRECORD, seq2, None)
        seq3 = SeqWrapper(SEQRECORD, seq3, None)
        seq4 = SeqWrapper(SEQRECORD, seq4, None)
        seq5 = SeqWrapper(SEQRECORD, seq5, None)
        seq6 = SeqWrapper(SEQRECORD, seq6, None)
        seqs = {SEQS_PASSED: [[seq1, seq2], [seq3, seq4], [seq5, seq6]],
                SEQS_FILTERED_OUT: []}

        seqs = filter_(seqs)
        assert len(seqs[SEQS_FILTERED_OUT]) == 1
        assert len(seqs[SEQS_PASSED]) == 2

    def test_dup_bin(self):
        seqs = '>seq1.f\naaa\n>seq1.r\naaa\n>seq2.f\naab\n>seq2.r\naaa\n'
        in_fhand = NamedTemporaryFile()
        in_fhand.write(seqs)
        in_fhand.flush()

        filter_bin = os.path.join(BIN_DIR, 'filter_duplicates')
        assert 'usage' in check_output([filter_bin, '-h'])
        result = check_output([filter_bin, in_fhand.name])
        assert'>seq1.f\naaa\n>seq2.f\naab' in result

        assert 'usage' in check_output([filter_bin, '-h'])
        result = check_output([filter_bin, in_fhand.name, '--paired_reads'])
        assert seqs in result

if __name__ == "__main__":
    # import sys;sys.argv = ['', 'FilterBowtie2Test.test_filter_by_bowtie2']
    unittest.main()
