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
from tempfile import NamedTemporaryFile
from StringIO import StringIO
from subprocess import check_output, CalledProcessError
import os.path
from random import choice

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from crumbs.split_mates import MatePairSplitter
from crumbs.settings import get_setting
from crumbs.seqio import (read_seqrecord_packets, write_seqrecord_packets,
                          read_seqrecords)
from crumbs.utils.bin_utils import BIN_DIR
from crumbs.utils.test_utils import TEST_DATA_DIR
from crumbs.utils.seq_utils import process_seq_packets

TITANIUM_LINKER = get_setting('TITANIUM_LINKER')
FLX_LINKER = get_setting('FLX_LINKER')

# pylint: disable=R0201
# pylint: disable=R0904


def create_a_matepair_file():
    'It creates a matepair fasta file'

    seq_5 = 'CTAGTCTAGTCGTAGTCATGGCTGTAGTCTAGTCTACGATTCGTATCAGTTGTGTGAC'
    seq_3 = 'ATCGATCATGTTGTATTGTGTACTATACACACACGTAGGTCGACTATCGTAGCTAGT'
    mate_seq = seq_5 + TITANIUM_LINKER + seq_3
    mate_fhand = NamedTemporaryFile(suffix='.fasta')
    mate_fhand.write('>seq1\n' + mate_seq + '\n')
    mate_fhand.flush()
    return mate_fhand


class MateSplitterTest(unittest.TestCase):
    'It tests the splitting of mate pairs'
    def test_split_mate(self):
        'It tests the function that splits seqs using segments'
        # pylint: disable=W0212
        seq = 'aaatttccctt'
        seqrecord = SeqRecord(Seq(seq), id='seq')
        # fake class to test
        splitter = MatePairSplitter([seqrecord])
        # segment beginning
        seqs = splitter._split_by_mate_linker(seqrecord, ([(0, 3)], False))
        assert str(seqs[0].seq) == 'ttccctt'
        assert seqs[0].id == 'seq'

        # segment at end
        seqs = splitter._split_by_mate_linker(seqrecord, ([(7, 10)], False))
        assert  str(seqs[0].seq) == 'aaatttc'
        assert seqs[0].id == 'seq'

        # segmnent in the middle
        seqs = splitter._split_by_mate_linker(seqrecord, ([(4, 7)], True))
        assert str(seqs[0].seq) == 'aaat'
        assert str(seqs[1].seq) == 'ctt'
        assert seqs[0].id == 'seq_pl.part1'
        assert seqs[1].id == 'seq_pl.part2'

        seqs = splitter._split_by_mate_linker(seqrecord, ([(4, 7)], False))
        assert seqs[0].id == r'seq\1'
        assert seqs[1].id == r'seq\2'

        seqs = splitter._split_by_mate_linker(seqrecord, ([(4, 6), (8, 9)],
                                                          False))
        assert str(seqs[0].seq) == 'aaat'
        assert str(seqs[1].seq) == 'c'
        assert str(seqs[2].seq) == 't'
        assert seqs[0].id == 'seq_mlc.part1'

        # all sequence is linker
        seqs = splitter._split_by_mate_linker(seqrecord, ([(0, 10)], False))
        assert not str(seqs[0].seq)

        # there's no segments
        seqs = splitter._split_by_mate_linker(seqrecord, ([], False))
        assert seqrecord.id == seqs[0].id
        assert str(seqrecord.seq) == str(seqs[0].seq)

    def test_many_reads(self):
        'It splits lots of reads to check that blast finds everything'

        linker = TITANIUM_LINKER

        def create_seq(index):
            'It creates a random seq with a linker'
            seq1 = ''.join(choice('ACTG') for i in range(100))
            seq2 = ''.join(choice('ACTG') for i in range(100))
            seq = seq1 + linker + seq2
            return SeqRecord(id='seq_' + str(index), seq=Seq(seq))

        # We want to test that blast reports all reads
        packet_size = get_setting('PACKET_SIZE')
        default_blast_max_target_size = 500
        assert packet_size > default_blast_max_target_size
        seqs = [create_seq(i) for i in range(1000)]
        splitter = MatePairSplitter()

        for index, seq in enumerate(splitter(seqs)):
            seq_index = index // 2
            pair_index = (index % 2) + 1
            expected_id = 'seq_' + str(seq_index) + '\\' + str(pair_index)
            assert  seq.id == expected_id

    def test_split_mates(self):
        'It tests the detection of oligos in sequence files'

        mate_fhand = NamedTemporaryFile(suffix='.fasta')
        linker = TITANIUM_LINKER

        # a complete linker
        seq5 = 'CTAGTCTAGTCGTAGTCATGGCTGTAGTCTAGTCTACGATTCGTATCAGTTGTGTGAC'
        seq3 = 'ATCGATCATGTTGTATTGTGTACTATACACACACGTAGGTCGACTATCGTAGCTAGT'

        mate_fhand.write('>seq1\n' + seq5 + linker + seq3 + '\n')
        # no linker
        mate_fhand.write('>seq2\n' + seq5 + '\n')
        # a partial linker
        mate_fhand.write('>seq3\n' + seq5 + linker[2:25] + seq3 + '\n')
        # the linker is 5 prima
        mate_fhand.write('>seq4\n' + linker[10:] + seq3 + '\n')
        # two linkers
        mate_fhand.write('>seq5\n' + linker + seq3 + FLX_LINKER + seq5 + '\n')
        # reverse linker
        rev_linker = str(get_setting('LINKERS')[1].seq.reverse_complement())
        mate_fhand.write('>seq6\n' + seq5 + rev_linker + seq3 + '\n')
        mate_fhand.flush()

        splitter = MatePairSplitter()
        new_seqs = []
        for packet in read_seqrecord_packets([mate_fhand], 2):
            new_seqs.append(splitter(packet))

        out_fhand = StringIO()
        write_seqrecord_packets(out_fhand, new_seqs, file_format='fasta')

        result = out_fhand.getvalue()
        xpect = r'>seq1\1'
        xpect += '\n'
        xpect += 'CTAGTCTAGTCGTAGTCATGGCTGTAGTCTAGTCTACGATTCGTATCAGTTGTGTGAC\n'
        xpect += r'>seq1\2'
        xpect += '\n'
        xpect += 'ATCGATCATGTTGTATTGTGTACTATACACACACGTAGGTCGACTATCGTAGCTAGT\n'
        xpect += '>seq2\n'
        xpect += 'CTAGTCTAGTCGTAGTCATGGCTGTAGTCTAGTCTACGATTCGTATCAGTTGTGTGAC\n'
        xpect += '>seq3_pl.part1\n'
        xpect += 'CTAGTCTAGTCGTAGTCATGGCTGTAGTCTAGTCTACGATTCGTATCAGTTGTGTG\n'
        xpect += '>seq3_pl.part2\n'
        xpect += 'GTGTACTATACACACACGTAGGTCGACTATCGTAGCTAGT\n'
        xpect += '>seq4\n'
        xpect += 'ATCGATCATGTTGTATTGTGTACTATACACACACGTAGGTCGACTATCGTAGCTAGT\n'
        xpect += '>seq5_mlc.part1\n'
        xpect += 'TCGTATAACTTCGTATAATGTATGCTATACGAAGTTATTACGATCGATCATGTTGTAT'
        xpect += 'TG\n'
        xpect += 'TGTACTATACACACACGTAGGTCGACTATCGTAGCTAGT\n'
        xpect += '>seq5_mlc.part2\n'
        xpect += 'ACCTAGTCTAGTCGTAGTCATGGCTGTAGTCTAGTCTACGATTCGTATCAGTTGTGTGAC'
        xpect += '\n'
        xpect += r'>seq6\1'
        xpect += '\n'
        xpect += 'CTAGTCTAGTCGTAGTCATGGCTGTAGTCTAGTCTACGATTCGTATCAGTTGTGTGAC\n'
        xpect += r'>seq6\2'
        xpect += '\n'
        xpect += 'ATCGATCATGTTGTATTGTGTACTATACACACACGTAGGTCGACTATCGTAGCTAGT\n'
        assert xpect == result

    @staticmethod
    def test_giuseppe_reads():
        'It splits some real reads'
        seq_fpath = os.path.join(TEST_DATA_DIR, '454_reads.fastq')
        linker_fpath = os.path.join(TEST_DATA_DIR, 'linkers.fasta')
        linkers = list(read_seqrecords([open(linker_fpath)]))

        splitter = MatePairSplitter(linkers=linkers)
        new_seqs = []
        for packet in read_seqrecord_packets([open(seq_fpath)], 2):
            new_seqs.extend(splitter(packet))
        seq_names = [seq.name for seq in new_seqs]
        assert 'G109AZL01BJHT8\\1' in seq_names
        assert 'G109AZL01BJHT8\\2' in seq_names
        assert len(new_seqs) == 20

        # test with process_seq_packet
        seq_fpath = os.path.join(TEST_DATA_DIR, '454_reads.fastq')
        linker_fpath = os.path.join(TEST_DATA_DIR, 'linkers.fasta')
        linkers = list(read_seqrecords([open(linker_fpath)]))

        splitter = MatePairSplitter(linkers=linkers)
        seq_packets = read_seqrecord_packets([open(seq_fpath)], 2)
        seq_packets = process_seq_packets(seq_packets, [splitter])[0]

        new_seqs = [seq for l in list(seq_packets) for seq in l]
        seq_names = [seq.name for seq in new_seqs]
        assert 'G109AZL01BJHT8\\1' in seq_names
        assert 'G109AZL01BJHT8\\2' in seq_names
        assert len(new_seqs) == 20

        # reads 2
        seq_fpath = os.path.join(TEST_DATA_DIR, '454_reads2.fastq')
        linker_fpath = os.path.join(TEST_DATA_DIR, 'linkers.fasta')
        linkers = list(read_seqrecords([open(linker_fpath)]))

        splitter = MatePairSplitter(linkers=linkers)
        new_seqs = []
        for packet in read_seqrecord_packets([open(seq_fpath)], 2):
            new_seqs.extend(splitter(packet))
        seq_names = [seq.name for seq in new_seqs]
        assert 'G109AZL01D8U3X\\1' in seq_names
        assert 'G109AZL01D8U3X\\2' in seq_names
        assert len(new_seqs) == 20

        # test with process_seq_packet
        seq_fpath = os.path.join(TEST_DATA_DIR, '454_reads2.fastq')
        linker_fpath = os.path.join(TEST_DATA_DIR, 'linkers.fasta')
        linkers = list(read_seqrecords([open(linker_fpath)]))

        splitter = MatePairSplitter(linkers=linkers)
        seq_packets = read_seqrecord_packets([open(seq_fpath)], 2)
        seq_packets = process_seq_packets(seq_packets, [splitter])[0]

        new_seqs = [seq for l in list(seq_packets) for seq in l]
        seq_names = [seq.name for seq in new_seqs]
        assert 'G109AZL01D8U3X\\1' in seq_names
        assert 'G109AZL01D8U3X\\2' in seq_names
        assert len(new_seqs) == 20


class SplitMatesBinTest(unittest.TestCase):
    'It tests the sff_extract binary'

    def test_matepair_bin(self):
        'It tests the split mate pairs binary'
        mate_bin = os.path.join(BIN_DIR, 'split_matepairs')
        stdout = check_output([mate_bin, '-h'])
        assert 'usage' in stdout

        out_fhand = NamedTemporaryFile(suffix='.fasta')
        seq_fpath = os.path.join(TEST_DATA_DIR, '454_reads2.fastq')
        linkers = '454'
        cmd = [mate_bin, '-o', out_fhand.name, '-l', linkers,
               seq_fpath]
        check_output(cmd)
        result = open(out_fhand.name).read()

        assert r'G109AZL01D8U3X\1' in  result
        assert r'G109AZL01D8U3X\2' in  result

        out_fhand = NamedTemporaryFile(suffix='.fasta')
        seq_fpath = os.path.join(TEST_DATA_DIR, '454_reads.fastq')
        linkers = '454'
        cmd = [mate_bin, '-o', out_fhand.name, '-l', linkers,
               seq_fpath]
        check_output(cmd)
        result = open(out_fhand.name).read()
        assert r'@G109AZL01BJHT8\1' in  result
        assert r'@G109AZL01BJHT8\2' in  result

        mate_fhand = create_a_matepair_file()
        out_fhand = NamedTemporaryFile(suffix='.fasta')

        cmd = [mate_bin, '-o', out_fhand.name, '-l', TITANIUM_LINKER,
               mate_fhand.name]
        from subprocess import call
        call(cmd)
        result = open(out_fhand.name).read()
        assert result.startswith(r'>seq1\1')

        cmd = [mate_bin, '-o', out_fhand.name, '-l', '454',
               mate_fhand.name]
        check_output(cmd)
        result = open(out_fhand.name).read()
        assert result.startswith(r'>seq1\1')

        # Error if a file path is given
        cmd = [mate_bin, '-o', out_fhand.name, '-l', 'adaptors.fasta',
               mate_fhand.name]
        stderr = NamedTemporaryFile(suffix='.stderr')
        try:
            check_output(cmd, stderr=stderr)
            self.fail('Error expected')
        except CalledProcessError:
            pass

        # Error if not DNA is given
        cmd = [mate_bin, '-o', out_fhand.name, '-l', 'iamnotasequence',
               mate_fhand.name]
        stderr = NamedTemporaryFile(suffix='.stderr')
        try:
            check_output(cmd, stderr=stderr)
            self.fail('Error expected')
        except CalledProcessError:
            pass

        # ion_torrent option is recognized
        cmd = [mate_bin, '-o', out_fhand.name, '-l', 'ion_torrent',
               seq_fpath]
        check_output(cmd)

    def test_parallel_bin(self):
        'The mate pairs binary runs in parallel'
        mate_bin = os.path.join(BIN_DIR, 'split_matepairs')
        mate_fhand = create_a_matepair_file()
        mate_fhand.seek(0)

        out_fhand = NamedTemporaryFile(suffix='.fasta')

        cmd = [mate_bin, '-o', out_fhand.name, '-l', TITANIUM_LINKER,
               mate_fhand.name, '-p', '2']
        check_output(cmd)
        result = open(out_fhand.name).read()
        assert result.startswith(r'>seq1\1')


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'SplitMatesBinTest.test_matepair_bin']
    unittest.main()
