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

from crumbs.seq import SeqWrapper, SeqItem
from crumbs.utils.tags import SEQITEM
from crumbs.bulk_filters import (filter_duplicates, _read_pairs,
                                 _seqitem_pairs_equal)
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
FASTQ_NO_DUPS3 = '''@CUESXEL837 1:Y:18:ATCACG
TAATACACCCAGTCTCAATTCCATCCTGGGAACTAAGT
+
AEGDFG5GGEGGF;EGD=D@>GCCGFFGGGCECFE:D@
@CUESXEL837 2:Y:18:ATCACG
CGCTCTTCTTCGCCACCGTGTTCTTGATATCGGCCTTC
+
ABCBEGGGGFGGGGGGGGGGGGGGGGBGGGA<EE=515
'''
FASTQ_PAIR = '''@CUESXEL834 1:Y:18:ATCACG
TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
+
CG?BEGGGGFGGGGGGGGGGGGGGGGBGGGA<EE=515
@CUESXEL834 2:Y:18:ATCACG
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
+
CG?BEGGGGFGGGGGGGGGGGGGGGGBGGGA<EE=515
'''


def _test_filter_duplicates(paired_reads, n_seqs_packet):
    assert isinstance(n_seqs_packet, int) or n_seqs_packet == None
    in_fhand = NamedTemporaryFile()
    fastq_with_dups = (FASTQ_NO_DUPS1 + FASTQ_DUPS + FASTQ_NO_DUPS2
                       + FASTQ_NO_DUPS3)
    in_fhand.write(fastq_with_dups)
    in_fhand.flush()
    in_fhand = open(in_fhand.name)
    out_fhand = NamedTemporaryFile()
    filter_duplicates([in_fhand], out_fhand, paired_reads, n_seqs_packet)
    filtered_pairs = list(_read_pairs([open(out_fhand.name)],
                                      paired_reads))
    fastq_no_dups = FASTQ_NO_DUPS1 + FASTQ_NO_DUPS2 + FASTQ_NO_DUPS3
    expected_pairs = list(_read_pairs([StringIO(fastq_no_dups)],
                                        paired_reads))
    #print 'filtered_pairs ->', filtered_pairs
    #print 'expected_pairs ->', expected_pairs
    #print len(filtered_pairs), len(expected_pairs)
    #assert len(filtered_pairs) == len(expected_pairs)
    for pair1 in expected_pairs:
        counts = 0
        for pair2 in filtered_pairs:
            if _seqitem_pairs_equal(pair1, pair2):
                counts += 1
        assert counts == 1
    in_fhand.close()


class FilterDuplicatesTest(unittest.TestCase):

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

    def test_filter_duplicates(self):
        options1 = [True, False]
        options2 = [None, 3]
        for option1 in options1:
            for option2 in options2:
                _test_filter_duplicates(paired_reads=option1,
                                        n_seqs_packet=option2)

    def test_dup_bin(self):
        seqs = '>seq1.f\naaa\n+\nHHHH\n>seq1.r\naaa\n+\nHHHH\n'
        seqs += '>seq2.f\naab\n+\nHHHH\n>seq2.r\naaa\n+\nHHHH\n'
        in_fhand = NamedTemporaryFile()
        in_fhand.write(seqs)
        in_fhand.flush()

        filter_bin = os.path.join(BIN_DIR, 'filter_duplicates')
        assert 'usage' in check_output([filter_bin, '-h'])
        result = check_output([filter_bin, in_fhand.name])
        assert'>seq1.f\naaa\n+\nHHHH\n>seq2.f\naab\n+\nHHHH\n' in result
        result = check_output([filter_bin], stdin=in_fhand)
        assert'>seq1.f\naaa\n+\nHHHH\n>seq2.f\naab\n+\nHHHH\n' in result
        result = check_output([filter_bin, in_fhand.name, '-d',
                               '/home/carlos/devel/tmp'])
        assert'>seq1.f\naaa\n+\nHHHH\n>seq2.f\naab\n+\nHHHH\n' in result
        result = check_output([filter_bin, in_fhand.name, '-m', '3'])
        assert'>seq1.f\naaa\n+\nHHHH\n>seq2.f\naab\n+\nHHHH\n' in result
        result = check_output([filter_bin, in_fhand.name, '--paired_reads'])
        assert seqs in result


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'FilterDuplicatesTest.test_dup_bin']
    unittest.main()
