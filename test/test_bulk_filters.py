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
from tempfile import NamedTemporaryFile
from StringIO import StringIO

from crumbs.seq import (SeqWrapper, SeqItem, get_title,
                        get_str_seq, get_str_qualities)
from crumbs.utils.tags import SEQITEM
from crumbs.bulk_filters import (filter_duplicates, _tabbed_pairs_equal,
                                 _read_pairs, _convert_fastq_to_tabbed_pairs,
                                 _tabbed_pair_to_seqs_seqitem,
                                 _pair_to_tabbed_str, _seqitem_pairs_equal)

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
TABBED_PAIR1 = '''@CUESXEL836 1:Y:18:ATCACG\tTAATACACCCAGTCTCAATTCCATCCTGGGAACTAAGT\tTTTDFG5GGEGGF;EGD=D@>GCCGFFGGGCECFE:D@\t@CUESXEL836 2:Y:18:ATCACG\tTCATTACGTAGCTCCGGCTCCGCCATGTCTGTTCCTTC\tABCBEGGGGFGGGGGGGGGGGGGGGGBGGGA<EE=515\n'''
TABBED_PAIR2 = '''@CUESXEL834 1:Y:18:ATCACG\tTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\tCG?BEGGGGFGGGGGGGGGGGGGGGGBGGGA<EE=515\t@CUESXEL834 2:Y:18:ATCACG\tAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\tCG?BEGGGGFGGGGGGGGGGGGGGGGBGGGA<EE=515\n'''

FASTQ_PAIR = '''@CUESXEL834 1:Y:18:ATCACG
TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
+
CG?BEGGGGFGGGGGGGGGGGGGGGGBGGGA<EE=515
@CUESXEL834 2:Y:18:ATCACG
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
+
CG?BEGGGGFGGGGGGGGGGGGGGGGBGGGA<EE=515
'''


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
        assert _seqitem_pairs_equal(seq1, seq3)
        assert _seqitem_pairs_equal(seq1, seq2) == False
        assert _seqitem_pairs_equal(seq1, pair1) == False
        assert _seqitem_pairs_equal(pair1, seq2) == False

    def test_pair_to_tabbed_str(self):
        expected_tabbed_pair = TABBED_PAIR1.rstrip()
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
        pair = (seq1, seq2)
        tabbed_pair = _pair_to_tabbed_str(pair).rstrip()
        #print tabbed_pair
        #print expected_tabbed_pair
        assert expected_tabbed_pair == tabbed_pair

    def test_fastq_to_tabbed_pairs(self):
        fastq1 = FASTQ_DUPS
        fastq2 = FASTQ_PAIR
        in_fhand1 = StringIO(fastq1)
        in_fhand2 = StringIO(fastq2)
        out_fhand = StringIO()
        out_fhand2 = StringIO()
        _convert_fastq_to_tabbed_pairs([in_fhand1, in_fhand2], out_fhand)
        _convert_fastq_to_tabbed_pairs([in_fhand1], out_fhand2)
        lines2 = out_fhand.getvalue().splitlines()
        lines = out_fhand.getvalue().splitlines()
        expected_lines = [TABBED_PAIR1, TABBED_PAIR2]
        assert len(lines) == 2
        assert lines[0].strip() == expected_lines[0].strip()
        assert lines[1].strip() == expected_lines[1].strip()
        assert lines[0].strip() == lines2[0].strip()

    def test_line_to_seqitem(self):
        tabbed_pair = TABBED_PAIR1
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
        expected_seqs = [seq1, seq2]
        seqs = _tabbed_pair_to_seqs_seqitem(tabbed_pair, 'fastq')
        for read1, read2 in zip(seqs, expected_seqs):
            assert get_title(read1) == get_title(read2)
            assert get_str_seq(read1) == get_str_seq(read2)
            assert get_str_qualities(read1) == get_str_qualities(read2)

    def test_filter_duplicates(self):
        in_fhand = NamedTemporaryFile()
        fastq_with_dups = FASTQ_NO_DUPS1 + FASTQ_DUPS + FASTQ_NO_DUPS2
        #print fastq_with_dups
        in_fhand.write(fastq_with_dups)
        in_fhand.flush()
        filtered_pairs = list(filter_duplicates([in_fhand.name]))

        fastq_no_dups = FASTQ_NO_DUPS1 + FASTQ_NO_DUPS2
        expected_pairs = list(_read_pairs([StringIO(fastq_no_dups)]))
        #print 'filtered_pairs ->', filtered_pairs
        #print 'expected_pairs ->', expected_pairs
        assert len(filtered_pairs) == len(expected_pairs)
        for pair1 in expected_pairs:
            counts = 0
            for pair2 in filtered_pairs:
                if _seqitem_pairs_equal(pair1, pair2):
                    counts += 1
            assert counts == 1
        in_fhand.close()

    def test_no_pairs(self):
        in_fhand1 = NamedTemporaryFile()
        in_fhand2 = NamedTemporaryFile()
        in_fpaths = [in_fhand1.name, in_fhand2.name]
        filtered_pairs = list(filter_duplicates(in_fpaths))
        in_fhand1.close()
        in_fhand2.close()
        assert not filtered_pairs


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'FilterDuplicatesTest.test_no_pairs']
    unittest.main()
