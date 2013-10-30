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
import os
from subprocess import check_output
from tempfile import NamedTemporaryFile
from StringIO import StringIO

from crumbs.seq import (SeqWrapper, SeqItem, get_title,
                        get_str_seq, get_str_qualities)
from crumbs.utils.tags import SEQITEM
from crumbs.bulk_filters import (filter_duplicates, _tabbed_pairs_equal,
                                 _read_pairs, _convert_fastq_to_tabbed_pairs,
                                 _tabbed_pair_to_seqs_seqitem,
                                 _pair_to_tabbed_str, _seqitem_pairs_equal,
                                 _unique_and_to_pairs)
from crumbs.utils.bin_utils import BIN_DIR

FASTQ_NO_DUPS1 = '''@CUESXEL822 1:Y:18:ATCACG
TAATACACCCAGTCTCAATTCCATCCTGGGAACTAAGT
+
AEGDFG5GGEGGF;EGD=D@>GCCGFFGGGCECFE:D@
@CUESXEL822 2:Y:18:ATCACG
TCATTACGTAGCTCCGGCTCCGCCATGTCTGTTCCTTC
+
CG?BEGGGGFGGGGGGGGGGGGGGGGBGGGA<EE=515
@CUESXEL824 1:Y:18:ATCACG
GATTGAAGCTCCAAACCGCCATGTTCACCACCGCAAGC
+
HHGEHD8EEHHHDGHHHHHHHHHCEHHHHDHHHHEHHH
@CUESXEL824 2:Y:18:ATCACG
TGCTTGCTGCACTTTGATGTTATTATCTGTGTTGTGTT
+
AA=AF7CDEDAFFDF@5D>D;FCF;GGGDGGEGGGFGE
'''
FASTQ_NO_DUPS2 = '''@CUESXEL830 1:Y:18:ATCACG
CGCTCTTCTTCGCCACCGTGTTCTTGATATCGGCCTTC
+
HHHHHHHHHHHHHHIHHHHHHGGHHHFHHHHHHHEHHH
@CUESXEL830 2:Y:18:ATCACG
TAATACACCCAGTCTCAATTCCATCCTGGGAACTAAGT
+
AEGDFG5GGEGGF;EGD=D@>GCCGFFGGGCECFE:D@
@CUESXEL832 1:Y:18:ATCACG
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
+
CG?BEGGGGFGGGGGGGGGGGGGGGGBGGGA<EE=515
@CUESXEL832 2:Y:18:ATCACG
TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
+
CG?BEGGGGFGGGGGGGGGGGGGGGGBGGGA<EE=515
@CUESXEL834 1:Y:18:ATCACG
TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
+
CG?BEGGGGFGGGGGGGGGGGGGGGGBGGGA<EE=515
@CUESXEL834 2:Y:18:ATCACG
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
+
CG?BEGGGGFGGGGGGGGGGGGGGGGBGGGA<EE=515
'''
FASTQ_DUPS = '''@CUESXEL836 1:Y:18:ATCACG
TAATACACCCAGTCTCAATTCCATCCTGGGAACTAAGT
+
TTTDFG5GGEGGF;EGD=D@>GCCGFFGGGCECFE:D@
@CUESXEL836 2:Y:18:ATCACG
TCATTACGTAGCTCCGGCTCCGCCATGTCTGTTCCTTC
+
ABCBEGGGGFGGGGGGGGGGGGGGGGBGGGA<EE=515
'''
TABBED_PAIR1 = '@CUESXEL836 1:Y:18:ATCACG\tTAATACACCCAGTCTCAATTCCATCCTGGGAACTA'
TABBED_PAIR1 += 'AGT\tTTTDFG5GGEGGF;EGD=D@>GCCGFFGGGCECFE:D@\t@CUESXEL836 2:Y:'
TABBED_PAIR1 += '18:ATCACG\tTCATTACGTAGCTCCGGCTCCGCCATGTCTGTTCCTTC\tABCBEGGGGF'
TABBED_PAIR1 += 'GGGGGGGGGGGGGGGGBGGGA<EE=515\n'
TABBED_PAIR2 = '@CUESXEL834 1:Y:18:ATCACG\tTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT'
TABBED_PAIR2 += 'TTT\tCG?BEGGGGFGGGGGGGGGGGGGGGGBGGGA<EE=515\t@CUESXEL834 2:Y:'
TABBED_PAIR2 += '18:ATCACG\tAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\tCG?BEGGGGF'
TABBED_PAIR2 += 'GGGGGGGGGGGGGGGGBGGGA<EE=515\n'

TABBED_SEQ1 = '@CUESXEL836 1:Y:18:ATCACG\tTAATACACCCAGTCTCAATTCCATCCTGGGAACTA'
TABBED_SEQ1 += 'AGT\tTTTDFG5GGEGGF;EGD=D@>GCCGFFGGGCECFE:D@\n'
TABBED_SEQ2 = '@CUESXEL836 2:Y:18:ATCACG\tTCATTACGTAGCTCCGGCTCCGCCATGTCTGTTCC'
TABBED_SEQ2 += 'TTC\tABCBEGGGGFGGGGGGGGGGGGGGGGBGGGA<EE=515\n'
TABBED_SEQ3 = '@CUESXEL834 1:Y:18:ATCACG\tTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT'
TABBED_SEQ3 += 'TTT\tCG?BEGGGGFGGGGGGGGGGGGGGGGBGGGA<EE=515\n'
TABBED_SEQ4 = '@CUESXEL834 2:Y:18:ATCACG\tAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA'
TABBED_SEQ4 += 'AAA\tCG?BEGGGGFGGGGGGGGGGGGGGGGBGGGA<EE=515\n'

FASTQ_PAIR = '''@CUESXEL834 1:Y:18:ATCACG
TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
+
CG?BEGGGGFGGGGGGGGGGGGGGGGBGGGA<EE=515
@CUESXEL834 2:Y:18:ATCACG
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
+
CG?BEGGGGFGGGGGGGGGGGGGGGGBGGGA<EE=515
'''
TABBED_PAIR1_FASTA = '>CUESXEL836 1:Y:18:ATCACG\tTAATACACCCAGTCTCAATTCCATCCTGG'
TABBED_PAIR1_FASTA += 'GAACTAAGT\t>CUESXEL836 2:Y:18:ATCACG\tTCATTACGTAGCTCCGG'
TABBED_PAIR1_FASTA += 'CTCCGCCATGTCTGTTCCTTC\n'
TABBED_SEQ1_FASTA = '>CUESXEL836 1:Y:18:ATCACG\tTAATACACCCAGTCTCAATTCCATCCTGGG'
TABBED_SEQ1_FASTA += 'AACTAAGT\n'


def _test_filter_duplicates(paired_reads):
    in_fhand = NamedTemporaryFile()
    fastq_with_dups = FASTQ_NO_DUPS1 + FASTQ_DUPS + FASTQ_NO_DUPS2
    in_fhand.write(fastq_with_dups)
    in_fhand.flush()
    out_fhand = NamedTemporaryFile()
    filter_duplicates([in_fhand], out_fhand, paired_reads)
    filtered_pairs = list(_read_pairs([open(out_fhand.name)],
                                      paired_reads))
    fastq_no_dups = FASTQ_NO_DUPS1 + FASTQ_NO_DUPS2
    expected_pairs = list(_read_pairs([StringIO(fastq_no_dups)],
                                        paired_reads))
    #print 'filtered_pairs ->', filtered_pairs
    #print 'expected_pairs ->', expected_pairs
    #assert len(filtered_pairs) == len(expected_pairs)
    for pair1 in expected_pairs:
        counts = 0
        for pair2 in filtered_pairs:
            if _seqitem_pairs_equal(pair1, pair2):
                counts += 1
        assert counts == 1
    in_fhand.close()


def _test_fastq_to_tabbed_pairs(paired_reads):
    fastq1 = FASTQ_DUPS
    fastq2 = FASTQ_PAIR
    in_fhand1 = StringIO(fastq1)
    in_fhand2 = StringIO(fastq2)
    lines = list(_convert_fastq_to_tabbed_pairs([in_fhand1, in_fhand2],
                                                    paired_reads))
    in_fhand1.seek(0)
    lines2 = list(_convert_fastq_to_tabbed_pairs([in_fhand1],
                                                     paired_reads))

    expected_paired_lines = [TABBED_PAIR1, TABBED_PAIR2]
    expected_single_lines = [TABBED_SEQ1, TABBED_SEQ2, TABBED_SEQ3,
                             TABBED_SEQ4]
    expected_lines = {True: expected_paired_lines,
                      False: expected_single_lines}

    assert len(lines) == len(expected_lines[paired_reads])
    for i in range(len(lines)):
        assert lines[i].strip() == expected_lines[paired_reads][i].strip()
    for i in range(len(lines2)):
        assert lines2[i].strip() == expected_lines[paired_reads][i].strip()


def _test_pair_to_tabbed_str(paired_reads):
    expected_tabbed_paired_pair = TABBED_PAIR1.rstrip()
    expected_tabbed_single_pair = TABBED_SEQ1.rstrip()
    expected_tabbed_pair = {True: expected_tabbed_paired_pair,
                            False: expected_tabbed_single_pair}
    seq1 = SeqWrapper(SEQITEM, SeqItem('CUESXEL836',
                        ['@CUESXEL836 1:Y:18:ATCACG\n',
                        'TAATACACCCAGTCTCAATTCCATCCTGGGAACTAAGT\n', '+\n',
                        'TTTDFG5GGEGGF;EGD=D@>GCCGFFGGGCECFE:D@\n']),
                        'fastq')
    seq2 = SeqWrapper(SEQITEM, SeqItem('CUESXEL836',
                        ['@CUESXEL836 2:Y:18:ATCACG\n',
                        'TCATTACGTAGCTCCGGCTCCGCCATGTCTGTTCCTTC\n', '+\n',
                        'ABCBEGGGGFGGGGGGGGGGGGGGGGBGGGA<EE=515\n']),
                        'fastq')
    pair = {True: (seq1, seq2), False: [seq1]}
    tabbed_pair = _pair_to_tabbed_str(pair[paired_reads]).rstrip()
    #print tabbed_pair
    #print expected_tabbed_pair
    assert expected_tabbed_pair[paired_reads] == tabbed_pair


def _test_tabbed_pair_to_seqitem(paired_reads, file_format):
    tabbed_paired_pair = {'fastq': TABBED_PAIR1.rstrip(),
                                   'fasta': TABBED_PAIR1_FASTA.rstrip()}
    tabbed_single_pair = {'fastq': TABBED_SEQ1.rstrip(),
                                   'fasta': TABBED_SEQ1_FASTA.strip()}
    tabbed_pair = {True: tabbed_paired_pair, False: tabbed_single_pair}
    seq1_fastq = SeqWrapper(SEQITEM, SeqItem('CUESXEL836',
                        ['@CUESXEL836 1:Y:18:ATCACG\n',
                        'TAATACACCCAGTCTCAATTCCATCCTGGGAACTAAGT\n', '+\n',
                        'TTTDFG5GGEGGF;EGD=D@>GCCGFFGGGCECFE:D@\n']),
                        'fastq')
    seq2_fastq = SeqWrapper(SEQITEM, SeqItem('CUESXEL836',
                        ['@CUESXEL836 2:Y:18:ATCACG\n',
                        'TCATTACGTAGCTCCGGCTCCGCCATGTCTGTTCCTTC\n', '+\n',
                        'ABCBEGGGGFGGGGGGGGGGGGGGGGBGGGA<EE=515\n']),
                        'fastq')
    seq1_fasta = SeqWrapper(SEQITEM, SeqItem('CUESXEL836',
                        ['@CUESXEL836 1:Y:18:ATCACG\n',
                        'TAATACACCCAGTCTCAATTCCATCCTGGGAACTAAGT\n']),
                        'fasta')
    seq2_fasta = SeqWrapper(SEQITEM, SeqItem('CUESXEL836',
                        ['@CUESXEL836 2:Y:18:ATCACG\n',
                        'TCATTACGTAGCTCCGGCTCCGCCATGTCTGTTCCTTC\n']),
                        'fasta')
    expected_seqs = {True: {'fastq': [seq1_fastq, seq2_fastq],
                     'fasta': [seq1_fasta, seq2_fasta]},
                     False: {'fastq': [seq1_fastq],
                     'fasta': [seq1_fasta]}}

    seqs = _tabbed_pair_to_seqs_seqitem(tabbed_pair[paired_reads][file_format],
                                        file_format)
    for read1, read2 in zip(seqs, expected_seqs[paired_reads][file_format]):
        assert get_title(read1) == get_title(read2)
        assert get_str_seq(read1) == get_str_seq(read2)
        if 'fastq' in file_format:
            assert get_str_qualities(read1) == get_str_qualities(read2)


class FilterDuplicatesTest(unittest.TestCase):

    def test_tabbed_pairs_equal(self):
        seq1 = '@seq1\tTAATAC\tTTTDFG'
        seq2 = '@seq2\tTCATTA\tABCBEG'
        seq3 = '@seq3\tTAATAC\tAGKDFG'
        seq4 = '@seq4\tACGCGT\tASDJFG'

        pair1 = '\t'.join([seq1, seq2])
        pair2 = '\t'.join([seq2, seq4])
        pair3 = '\t'.join([seq3, seq2])
        pair4 = '\t'.join([seq2, seq1])

        assert _tabbed_pairs_equal(pair1, pair3)
        assert _tabbed_pairs_equal(pair1, pair2) == False
        assert _tabbed_pairs_equal(pair1, pair4) == False
        assert _tabbed_pairs_equal(seq1, seq3)
        assert _tabbed_pairs_equal(seq1, seq2) == False
        assert _tabbed_pairs_equal(seq1, pair1) == False
        assert _tabbed_pairs_equal(pair1, seq2) == False

    def test_seqitem_pairs_equal(self):
        seq1 = SeqWrapper(SEQITEM, SeqItem('seq1',
                            ['@seq1\n', 'TAATAC\n', '+\n',
                             'TTTDFG\n']), 'fastq')
        seq2 = SeqWrapper(SEQITEM, SeqItem('seq2',
                            ['@seq2\n', 'TCATTA\n', '+\n',
                             'ABCBEG\n']), 'fastq')
        seq3 = SeqWrapper(SEQITEM, SeqItem('seq3',
                            ['@seq3\n', 'TAATAC\n', '+\n',
                             'TTTDFG\n']), 'fastq')
        seq4 = SeqWrapper(SEQITEM, SeqItem('seq4',
                            ['@seq4\n', 'ACGCGT\n', '+\n',
                             'ABCBEG\n']), 'fastq')
        pair1 = (seq1, seq2)
        pair2 = (seq2, seq4)
        pair3 = (seq3, seq2)
        pair4 = (seq2, seq1)

        assert _seqitem_pairs_equal(pair1, pair3)
        assert _seqitem_pairs_equal(pair1, pair2) == False
        assert _seqitem_pairs_equal(pair1, pair4) == False
        assert _seqitem_pairs_equal([seq1], [seq3])
        assert _seqitem_pairs_equal([seq1], [seq2]) == False
        assert _seqitem_pairs_equal([seq1], pair1) == False
        assert _seqitem_pairs_equal(pair1, seq2) == False

    def test_pair_to_tabbed_str(self):
        for option in [True, False]:
            _test_pair_to_tabbed_str(option)

    def test_fastq_to_tabbed_pairs(self):
        for option in [True, False]:
            _test_fastq_to_tabbed_pairs(option)

    def test_tabbed_pair_to_seqitem(self):
        for option in [True, False]:
            for file_format in ['fastq', 'fasta']:
                _test_tabbed_pair_to_seqitem(option, file_format)

    def test_unique_and_to_pairs(self):
        seq1 = '@seq1\tTAATAC\tTTTDFG'
        seq2 = '@seq2\tTCATTA\tABCBEG'
        seq3 = '@seq3\tTAATAC\tAGKDFG'
        seq4 = '@seq4\tACGCGT\tASDJFG'

        pair1 = '\t'.join([seq1, seq2])
        pair2 = '\t'.join([seq2, seq4])
        pair3 = '\t'.join([seq3, seq2])
        pair4 = '\t'.join([seq2, seq1])
        in_pairs = [pair1, pair3, pair2, pair4]
        sorted_pairs = '\n'.join(in_pairs)
        sorted_pairs2 = '\n'.join([seq4, seq3, seq1, seq2])
        in_fhand = StringIO(sorted_pairs)
        in_fhand2 = StringIO(sorted_pairs2)
        filtered_pairs = list(_unique_and_to_pairs(in_fhand, 'fastq'))
        filtered_pairs2 = list(_unique_and_to_pairs(in_fhand2, 'fastq'))
        seq1 = SeqWrapper(SEQITEM, SeqItem('seq1',
                            ['@seq1\n', 'TAATAC\n', '+\n',
                             'TTTDFG\n']), 'fastq')
        seq2 = SeqWrapper(SEQITEM, SeqItem('seq2',
                            ['@seq2\n', 'TCATTA\n', '+\n',
                             'ABCBEG\n']), 'fastq')
        seq3 = SeqWrapper(SEQITEM, SeqItem('seq3',
                            ['@seq3\n', 'TAATAC\n', '+\n',
                             'TTTDFG\n']), 'fastq')
        seq4 = SeqWrapper(SEQITEM, SeqItem('seq4',
                            ['@seq4\n', 'ACGCGT\n', '+\n',
                             'ABCBEG\n']), 'fastq')
        pair1 = (seq1, seq2)
        pair2 = (seq2, seq4)
        pair3 = (seq3, seq2)
        pair4 = (seq2, seq1)
        expected_pairs = [pair1, pair2, pair4]
        expected_pairs2 = [[seq3], [seq2], [seq4]]

        assert len(expected_pairs2) == len(filtered_pairs2)
        for pair1 in expected_pairs2:
            counts = 0
            for pair2 in filtered_pairs2:
                if _seqitem_pairs_equal(pair1, pair2):
                    counts += 1
            assert counts == 1
        in_fhand2.close()

        assert len(expected_pairs) == len(filtered_pairs)
        for pair1 in expected_pairs:
            counts = 0
            for pair2 in filtered_pairs:
                if _seqitem_pairs_equal(pair1, pair2):
                    counts += 1
            assert counts == 1
        in_fhand.close()

    def test_filter_duplicates(self):
        _test_filter_duplicates(paired_reads=True)
        _test_filter_duplicates(paired_reads=False)

    def test_dup_bin(self):
        seqs = '>seq1.f\naaa\n>seq1.r\naaa\n>seq2.f\naab\n>seq2.r\naaa\n'
        in_fhand = NamedTemporaryFile()
        in_fhand.write(seqs)
        in_fhand.flush()

        filter_bin = os.path.join(BIN_DIR, 'filter_duplicates')
        assert 'usage' in check_output([filter_bin, '-h'])
        result = check_output([filter_bin, in_fhand.name])
        assert'>seq1.f\naaa\n>seq2.f\naab' in result
        result = check_output([filter_bin, in_fhand.name, '--paired_reads'])
        assert seqs in result


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'FilterDuplicatesTest.test_filter_duplicates']
    unittest.main()
