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
from tempfile import NamedTemporaryFile
from subprocess import check_output
from cStringIO import StringIO

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from crumbs.trim import (TrimLowercasedLetters, TrimEdges, TrimOrMask,
                         TrimByQuality, TrimWithBlastShort,
    seq_to_trim_packets)
from crumbs.utils.bin_utils import BIN_DIR
from crumbs.utils.tags import (SEQRECORD, SEQITEM, TRIMMING_RECOMMENDATIONS,
                               VECTOR, ORPHAN_SEQS, SEQS_PASSED)
from crumbs.seq import get_str_seq, get_annotations, get_int_qualities, get_name
from crumbs.seqio import read_seq_packets
from crumbs.seq import SeqWrapper, SeqItem

FASTQ = '@seq1\naTCgt\n+\n?????\n@seq2\natcGT\n+\n?????\n'
FASTQ2 = '@seq1\nATCGT\n+\nA???A\n@seq2\nATCGT\n+\n?????\n'
FASTQ3 = '@seq1\nAAAAAATCGTTTTTTT\n+\n00000A???A000000\n'

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
    def test_trim_seqs():
        'It tests the trim seq function'
        seqs = []
        seqs.append([SeqWrapper(SEQRECORD, SeqRecord(Seq('aaCTTTC')), None)])
        seqs.append([SeqWrapper(SEQRECORD, SeqRecord(Seq('CTTCaa')), None)])
        seqs.append([SeqWrapper(SEQRECORD, SeqRecord(Seq('aaCTCaa')), None)])
        seqs.append([SeqWrapper(SEQRECORD, SeqRecord(Seq('actg')), None)])
        seqs.append([SeqWrapper(SEQRECORD, SeqRecord(Seq('AC')), None)])

        trim_packet = {SEQS_PASSED: seqs, ORPHAN_SEQS: []}
        trim_lowercased_seqs = TrimLowercasedLetters()
        trim = TrimOrMask()
        # pylint: disable=W0141
        trim_packet = trim(trim_lowercased_seqs(trim_packet))
        res = [get_str_seq(s) for l in trim_packet[SEQS_PASSED] for s in l]

        assert res == ['CTTTC', 'CTTC', 'CTC', 'AC']

        seqs = []
        seq = SeqItem('s', ['>s\n', 'aaCTTTC\n'])
        seqs.append([SeqWrapper(SEQITEM, seq, 'fasta')])
        trim_packet = {SEQS_PASSED: seqs, ORPHAN_SEQS: []}
        trim_packet = trim(trim_lowercased_seqs(trim_packet))
        res = [get_str_seq(s) for l in trim_packet[SEQS_PASSED] for s in l]

        assert res == ['CTTTC']

        # with pairs
        seq = SeqItem('s.f', ['>s.f\n', 'aaCTTTC\n'])
        seq1 = SeqItem('s.r', ['>s.r\n', 'aaCTTTC\n'])
        seq2 = SeqItem('s1.f', ['>s1.f\n', 'aa\n'])
        seq3 = SeqItem('s1.r', ['>s1.r\n', 'aaCTTTC\n'])
        seqs = []
        seqs.append([SeqWrapper(SEQITEM, seq, 'fasta'),
                     SeqWrapper(SEQITEM, seq1, 'fasta')])
        seqs.append([SeqWrapper(SEQITEM, seq2, 'fasta'),
                     SeqWrapper(SEQITEM, seq3, 'fasta')])
        trim_packet = {SEQS_PASSED: seqs, ORPHAN_SEQS: []}
        trim_lowercased_seqs = TrimLowercasedLetters()
        trim = TrimOrMask()
        # pylint: disable=W0141
        trim_packet = trim(trim_lowercased_seqs(trim_packet))
        res = [get_str_seq(s) for l in trim_packet[SEQS_PASSED] for s in l]
        orphan_res = [get_str_seq(s) for s in trim_packet[ORPHAN_SEQS]]
        assert orphan_res == ['CTTTC']
        assert ['CTTTC', 'CTTTC'] == res

        # no drag
        trim_packet = {SEQS_PASSED: seqs, ORPHAN_SEQS: []}
        trim_lowercased_seqs = TrimLowercasedLetters()
        trim = TrimOrMask()
        # pylint: disable=W0141
        trim_packet = trim(trim_lowercased_seqs(trim_packet))
        res = [get_str_seq(s) for l in trim_packet[SEQS_PASSED] for s in l]
        orphan_res = [get_name(s) for s in trim_packet[ORPHAN_SEQS]]
        assert orphan_res == ['s1.r']
        assert ['CTTTC', 'CTTTC'] == res


class TrimByCaseBinTest(unittest.TestCase):
    'It tests the trim_by_case binary'

    def test_trim_case_bin(self):
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

    def test_trim_in_pairs(self):
        trim_bin = os.path.join(BIN_DIR, 'trim_by_case')
        content = '>s.f\naattACT\n>s.r\naattACT\n>s1.f\naattACT\n>s1.r\naatt\n'
        fastq_fhand = _make_fhand(content)
        orphan_fhand = NamedTemporaryFile()
        result = check_output([trim_bin, fastq_fhand.name, '--paired_reads',
                               '--orphan_file', orphan_fhand.name])
        assert ">s.f\nACT\n>s.r\nACT" in  result
        assert  open(orphan_fhand.name).read() == '>s1.f\nACT\n'

        # no orphan file
        fastq_fhand = _make_fhand(content)
        result = check_output([trim_bin, fastq_fhand.name, '--paired_reads'])
        assert ">s.f\nACT\n>s.r\nACT" in  result


class TrimEdgesTest(unittest.TestCase):
    'It test the fixed number of bases trimming'

    def _some_seqs(self):
        'It returns some seqrecords.'
        seqs = []
        seq = SeqRecord(Seq('ACCG'), letter_annotations={'dummy': 'dddd'})
        seq = SeqWrapper(SEQRECORD, seq, None)
        seqs.append([seq])
        seq = SeqRecord(Seq('AAACCCGGG'))
        seq = SeqWrapper(SEQRECORD, seq, None)
        seqs.append([seq])
        trim_packet = {SEQS_PASSED: seqs, ORPHAN_SEQS: []}
        return trim_packet

    def test_edge_trimming(self):
        'It trims the edges'
        trim = TrimOrMask()
        trim_edges = TrimEdges(left=1)
        trim_packet = trim(trim_edges(self._some_seqs()))
        res = [get_str_seq(s) for l in trim_packet[SEQS_PASSED] for s in l]
        assert res == ['CCG', 'AACCCGGG']

        trim_edges = TrimEdges(right=1)
        trim_packet = trim(trim_edges(self._some_seqs()))
        res = [get_str_seq(s) for l in trim_packet[SEQS_PASSED] for s in l]
        assert res == ['ACC', 'AAACCCGG']

        trim_edges = TrimEdges(left=1, right=1)
        trim_packet = trim(trim_edges(self._some_seqs()))
        res = [get_str_seq(s) for l in trim_packet[SEQS_PASSED] for s in l]
        assert res == ['CC', 'AACCCGG']

        trim_edges = TrimEdges(left=2, right=2)
        trim_packet = trim(trim_edges(self._some_seqs()))
        res = [get_str_seq(s) for l in trim_packet[SEQS_PASSED] for s in l]
        assert res == ['ACCCG']

        trim_edges = TrimEdges(left=3, right=3)
        trim_packet = trim(trim_edges(self._some_seqs()))
        res = [get_str_seq(s) for l in trim_packet[SEQS_PASSED] for s in l]
        assert res == ['CCC']

        trim = TrimOrMask(mask=True)
        trim_edges = TrimEdges(left=1)
        trim_packet = trim(trim_edges(self._some_seqs()))
        res = [get_str_seq(s) for l in trim_packet[SEQS_PASSED] for s in l]
        assert res == ['aCCG', 'aAACCCGGG']

        trim_edges = TrimEdges(right=1)
        trim_packet = trim(trim_edges(self._some_seqs()))
        res = [get_str_seq(s) for l in trim_packet[SEQS_PASSED] for s in l]
        assert res == ['ACCg', 'AAACCCGGg']

        trim_edges = TrimEdges(left=1, right=1)
        trim_packet = trim(trim_edges(self._some_seqs()))
        res = [get_str_seq(s) for l in trim_packet[SEQS_PASSED] for s in l]
        assert res == ['aCCg', 'aAACCCGGg']

        trim_edges = TrimEdges(left=2, right=2)
        trim_packet = trim(trim_edges(self._some_seqs()))
        res = [get_str_seq(s) for l in trim_packet[SEQS_PASSED] for s in l]
        assert res == ['accg', 'aaACCCGgg']

        trim_edges = TrimEdges(left=3, right=3)
        trim_packet = trim(trim_edges(self._some_seqs()))
        res = [get_str_seq(s) for l in trim_packet[SEQS_PASSED] for s in l]
        assert res == ['accg', 'aaaCCCggg']

        # test overlapping mask
        trim1 = TrimEdges(left=3, right=3)
        trim2 = TrimEdges(left=4, right=4)
        trim_packet = trim(trim2(trim1(self._some_seqs())))
        res = [get_str_seq(s) for l in trim_packet[SEQS_PASSED] for s in l]
        assert res == ['accg', 'aaacCcggg']

        # With a SeqItem
        trim = TrimOrMask(mask=False)
        trim_edges = TrimEdges(left=1, right=1)

        seq = SeqItem('s', ['>s\n', 'ACTTTC\n'])
        seqs = [[SeqWrapper(SEQITEM, seq, 'fasta')]]
        trim_packet = {SEQS_PASSED: seqs, ORPHAN_SEQS: []}

        trim_packet = trim(trim_edges(trim_packet))
        res = [get_str_seq(s) for l in trim_packet[SEQS_PASSED] for s in l]
        assert res == ['CTTT']

        trim = TrimOrMask(mask=True)
        seq = SeqItem('s', ['>s\n', 'ACTTTC\n'])
        seqs = [[SeqWrapper(SEQITEM, seq, 'fasta')]]
        trim_packet = {SEQS_PASSED: seqs, ORPHAN_SEQS: []}
        trim_packet = trim(trim_edges(trim_packet))
        res = [get_str_seq(s) for l in trim_packet[SEQS_PASSED] for s in l]
        assert res == ['aCTTTc']

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


class TrimAndMaskTest(unittest.TestCase):
    'It tests the trimming and masking according to the recommendations.'
    def test_trimming(self):
        'The sequences are trimmed according to the recommendations.'
        seq1 = 'gggtctcatcatcaggg'.upper()
        seq = SeqRecord(Seq(seq1), annotations={TRIMMING_RECOMMENDATIONS: {}})
        seq = SeqWrapper(SEQRECORD, seq, None)
        seqs = [seq]
        trim_packet = {SEQS_PASSED: [seqs], ORPHAN_SEQS: []}

        trim_rec = get_annotations(seq)[TRIMMING_RECOMMENDATIONS]
        seq_trimmer = TrimOrMask()

        trim_rec['vector'] = [(0, 3), (8, 13)]
        get_annotations(seq)[TRIMMING_RECOMMENDATIONS] = trim_rec
        trim_packet2 = seq_trimmer(trim_packet)
        res = [get_str_seq(s) for l in trim_packet2[SEQS_PASSED] for s in l]
        assert res == ['CTCA']

        trim_rec['vector'] = [(0, 0), (8, 13)]
        get_annotations(seq)[TRIMMING_RECOMMENDATIONS] = trim_rec
        trim_packet2 = seq_trimmer(trim_packet)
        res = [get_str_seq(s) for l in trim_packet2[SEQS_PASSED] for s in l]
        assert res == ['GGTCTCA']

        trim_rec['vector'] = [(0, 1), (8, 12)]
        trim_rec['quality'] = [(1, 8), (13, 17)]
        get_annotations(seq)[TRIMMING_RECOMMENDATIONS] = trim_rec
        trim_packet2 = seq_trimmer(trim_packet)
        assert not trim_packet2[SEQS_PASSED]

        trim_rec['vector'] = [(0, 0), (8, 13)]
        trim_rec['quality'] = []
        get_annotations(seq)[TRIMMING_RECOMMENDATIONS] = trim_rec
        trim_packet2 = seq_trimmer(trim_packet)
        res = [get_str_seq(s) for l in trim_packet2[SEQS_PASSED] for s in l]
        assert res == ['GGTCTCA']
        trim_packet2[SEQS_PASSED][0][0]
        assert TRIMMING_RECOMMENDATIONS not in get_annotations(trim_packet2[SEQS_PASSED][0][0])


class TrimByQualityTest(unittest.TestCase):
    'It test the quality trimming'

    def test_quality_trimming(self):
        'It trims the edges'
        trim = TrimOrMask()

        trim_quality = TrimByQuality(window=5, threshold=30)

        seq = SeqRecord(Seq('ACTGCTGCATAAAA'))
        quals = [10, 10, 20, 30, 30, 30, 40, 40, 30, 30, 20, 20, 10, 10]
        seq.letter_annotations['phred_quality'] = quals
        seq = SeqWrapper(SEQRECORD, seq, None)
        trim_packet = {SEQS_PASSED: [[seq]], ORPHAN_SEQS: []}
        trim_packet2 = trim(trim_quality(trim_packet))
        seq2 = trim_packet2[SEQS_PASSED][0][0]
        assert get_int_qualities(seq2) == [20, 30, 30, 30, 40, 40, 30, 30, 20]

        # all bad
        trim_quality = TrimByQuality(window=5, threshold=60)
        trim_packet = {SEQS_PASSED: [[seq]], ORPHAN_SEQS: []}
        trim_packet2 = trim(trim_quality(trim_packet))
        assert not trim_packet2[SEQS_PASSED]

        # all OK
        trim_quality = TrimByQuality(window=5, threshold=5)
        trim_packet = {SEQS_PASSED: [[seq]], ORPHAN_SEQS: []}
        trim_packet2 = trim(trim_quality(trim_packet))
        seq2 = trim_packet2[SEQS_PASSED][0][0]
        assert get_int_qualities(seq2) == quals

        seq = SeqRecord(Seq('ACTGCTGCATAA'))
        quals = [20, 20, 20, 60, 60, 60, 60, 60, 20, 20, 20, 20]

        trim_quality = TrimByQuality(window=5, threshold=50)
        seq.letter_annotations['phred_quality'] = quals
        seq = SeqWrapper(SEQRECORD, seq, None)
        trim_packet = {SEQS_PASSED: [[seq]], ORPHAN_SEQS: []}
        trim_packet2 = trim(trim_quality(trim_packet))
        seq2 = trim_packet2[SEQS_PASSED][0][0]
        assert get_int_qualities(seq2) == [20, 60, 60, 60, 60, 60, 20]

        quals = [40, 18, 10, 40, 40, 5, 8, 30, 14, 3, 40, 40, 40, 11, 6, 5, 3,
                 20, 10, 12, 8, 5, 4, 7, 1]
        seq = SeqRecord(Seq('atatatatagatagatagatagatg'))
        seq.letter_annotations['phred_quality'] = quals
        seq = SeqWrapper(SEQRECORD, seq, None)
        trim_packet = {SEQS_PASSED: [[seq]], ORPHAN_SEQS: []}
        trim_quality = TrimByQuality(window=5, threshold=25)
        trim_packet2 = trim(trim_quality(trim_packet))
        seq2 = trim_packet2[SEQS_PASSED][0][0]
        assert get_int_qualities(seq2) == [40, 18, 10, 40, 40]

        quals = [40, 40, 13, 11, 40, 9, 40, 4, 27, 38, 40, 4, 11, 40, 40, 10,
                 10, 21, 3, 40, 9, 9, 12, 10, 9]
        seq = SeqRecord(Seq('atatatatatatatatatatatata'))
        seq.letter_annotations['phred_quality'] = quals
        seq = SeqWrapper(SEQRECORD, seq, None)
        trim_quality = TrimByQuality(window=5, threshold=25)
        trim_packet = {SEQS_PASSED: [[seq]], ORPHAN_SEQS: []}

        trim_packet2 = trim(trim_quality(trim_packet))
        seq2 = trim_packet2[SEQS_PASSED][0][0]
        expected = [40, 4, 27, 38, 40]
        assert get_int_qualities(seq2) == expected

        quals = [40, 40, 13, 11, 40, 9, 40, 4, 27, 38, 40, 4, 11, 40, 40, 10,
                 10, 21, 3, 40, 9, 9, 12, 10, 9]
        seq = SeqRecord(Seq('atatatatatatatatatatatata'))
        seq.letter_annotations['phred_quality'] = quals
        seq = SeqWrapper(SEQRECORD, seq, None)
        trim_packet = {SEQS_PASSED: [[seq]], ORPHAN_SEQS: []}
        trim_quality = TrimByQuality(window=5, threshold=25, trim_left=False)
        trim_packet2 = trim(trim_quality(trim_packet))
        seq2 = trim_packet2[SEQS_PASSED][0][0]
        assert get_int_qualities(seq2) == [40, 40, 13, 11, 40, 9, 40, 4, 27,
                                           38, 40]

        quals = [40, 40, 13, 11, 40, 9, 40, 4, 27, 38, 40, 4, 11, 40, 40, 10,
                 10, 21, 3, 40, 9, 9, 12, 10, 9]
        seq = SeqRecord(Seq('atatatatatatatatatatatata'))
        seq.letter_annotations['phred_quality'] = quals
        seq = SeqWrapper(SEQRECORD, seq, None)
        trim_packet = {SEQS_PASSED: [[seq]], ORPHAN_SEQS: []}
        trim_quality = TrimByQuality(window=5, threshold=25, trim_right=False)
        trim_packet2 = trim(trim_quality(trim_packet))
        seq2 = trim_packet2[SEQS_PASSED][0][0]
        assert get_int_qualities(seq2) == [40, 4, 27, 38, 40, 4, 11, 40, 40, 10,
                                       10, 21, 3, 40, 9, 9, 12, 10, 9]

        quals = [40, 40, 13, 11, 40, 9, 40, 4, 27, 38, 40, 4, 11, 40, 40, 10,
                 10, 21, 3, 40, 9, 9, 12, 10, 9]
        seq = SeqRecord(Seq('atatatatatatatatatatatata'))
        seq.letter_annotations['phred_quality'] = quals
        seq = SeqWrapper(SEQRECORD, seq, None)
        trim_quality = TrimByQuality(window=5, threshold=25, trim_right=False,
                                     trim_left=False)
        trim_packet = {SEQS_PASSED: [[seq]], ORPHAN_SEQS: []}
        trim_packet2 = trim(trim_quality(trim_packet))
        seq2 = trim_packet2[SEQS_PASSED][0][0]
        assert get_int_qualities(seq2) == quals

        # With SeqItems
        seq = SeqItem('s', ['@s\n', 'atatatatatatatatatatatata\n', '\n',
                            'II.,I*I%<GI%,II++6$I**-+*\n'])
        seq = SeqWrapper(SEQITEM, seq, 'fastq')
        trim_quality = TrimByQuality(window=5, threshold=25, trim_right=True,
                                     trim_left=False)

        trim_packet = {SEQS_PASSED: [[seq]], ORPHAN_SEQS: []}
        trim_packet2 = trim(trim_quality(trim_packet))
        seq2 = trim_packet2[SEQS_PASSED][0][0]
        assert seq2.object.lines[3] == 'II.,I*I%<GI\n'

    def test_trim_quality_bin(self):
        'It tests the trim_edges binary'
        trim_bin = os.path.join(BIN_DIR, 'trim_quality')
        assert 'usage' in check_output([trim_bin, '-h'])

        fastq_fhand = _make_fhand(FASTQ2)
        result = check_output([trim_bin, fastq_fhand.name])
        assert '@seq1\nATCGT\n+' in result

        fastq_fhand = _make_fhand(FASTQ3)
        result = check_output([trim_bin, fastq_fhand.name])
        assert result == '@seq1\nAATCGTT\n+\n0A???A0\n'

        fastq_fhand = _make_fhand(FASTQ3)
        result = check_output([trim_bin, fastq_fhand.name, '-r'])
        assert result == '@seq1\nAATCGTTTTTTT\n+\n0A???A000000\n'

        fastq_fhand = _make_fhand(FASTQ3)
        result = check_output([trim_bin, fastq_fhand.name, '-l'])
        assert result == '@seq1\nAAAAAATCGTT\n+\n00000A???A0\n'

# pylint: disable=C0301

FASTQ4 = '''@HWI-ST1203:122:C130PACXX:4:1101:13499:4144 1:N:0:CAGATC
AAGCAGTGGTATCAAAGCAGAGTACTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTCCAACCCTTTGCTTTTTTTTTTTTTCGAGGAGGAGGGT
+
@@@DBDDDDHFBHGGHFEA@GG<?FHHIIIIIIIIIGCCCCCCCCCCCCCCCCBCCBC###########################################
@HWI-ST1203:122:C130PACXX:4:1101:13623:4101 1:N:0:CAGATC
AAGCAGTGGTATCAACGCAGAGTACATGGGCGAGAAGAAGGATCCAAGTGGTGCCAAGGTTACCAAATCTGCAGCCAAGAAGGCTGGAAAGTGAACCGTGC
+
CCCFFFFFHHHHHJJIJJJHIIHGIIIIIIIJGGHIIJGHJIJIGHFG@FHGGHIJHHHHHFFFFFDEEEEEEDDBDCBDDDDDDDBADCD>C@DCDD<<<
@HWI-ST1203:122:C130PACXX:4:1101:13615:4163 1:N:0:CAGATC
GGAAGAGGAACAAGTGAGCAGCAGGACTGTATGATATTCTCATCTGAAGACAGGGACCATCATATTCCCCGGGAAACTCCGATGCCAGAGTATTAGCATGC
+
@1?DFFFFGHHHHIBGGHGHGEICCAGHHCFGHHIGGHFIHIIIJJJIJIJJJIJIIIJJJJJJJICEEHHFDADBCCDDDDDBBDDCAB@CCDEEDDEDC
'''


class TrimBlastShortTest(unittest.TestCase):
    'It tests the blast short adaptor trimming'
    def test_blast_short_trimming(self):
        'It trims oligos using blast-short'

        oligo1 = SeqRecord(Seq('AAGCAGTGGTATCAACGCAGAGTACATGGG'))
        oligo2 = SeqRecord(Seq('AAGCAGTGGTATCAACGCAGAGTACTTTTT'))
        oligo1 = SeqWrapper(SEQRECORD, oligo1, None)
        oligo2 = SeqWrapper(SEQRECORD, oligo2, None)

        adaptors = [oligo1, oligo2]

        blast_trim = TrimWithBlastShort(oligos=adaptors)
        fhand = StringIO(FASTQ4)
        seq_packets = read_seq_packets([fhand],
                                            prefered_seq_classes=[SEQRECORD])
        trim_packets = list(seq_to_trim_packets(seq_packets))
        trim_packets2 = blast_trim(trim_packets[0])
        # It should trim the first and the second reads.
        res = [get_annotations(s).get(TRIMMING_RECOMMENDATIONS, {}).get(VECTOR,
                                                                        [])
                            for l in trim_packets2[SEQS_PASSED] for s in l]
        assert res == [[(0, 29)], [(0, 29)], []]

        # With SeqItems
        oligo1 = SeqItem('oligo1', ['>oligo1\n',
                                    'AAGCAGTGGTATCAACGCAGAGTACATGGG\n'])
        oligo2 = SeqItem('oligo2', ['>oligo2\n',
                                    'AAGCAGTGGTATCAACGCAGAGTACTTTTT\n'])
        oligo1 = SeqWrapper(SEQITEM, oligo1, 'fasta')
        oligo2 = SeqWrapper(SEQITEM, oligo2, 'fasta')

        adaptors = [oligo1, oligo2]

        blast_trim = TrimWithBlastShort(oligos=adaptors)
        fhand = StringIO(FASTQ4)
        seq_packets = list(read_seq_packets([fhand],
                                            prefered_seq_classes=[SEQITEM]))
        trim_packets = list(seq_to_trim_packets(seq_packets))
        trim_packets2 = blast_trim(trim_packets[0])
        # It should trim the first and the second reads.
        res = [get_annotations(s).get(TRIMMING_RECOMMENDATIONS, {}).get(VECTOR,
                                                                        [])
                            for l in trim_packets2[SEQS_PASSED] for s in l]
        assert res == [[(0, 29)], [(0, 29)], []]

    def test_trim_oligos_bin(self):
        'It tests the trim_blast_short binary'
        trim_bin = os.path.join(BIN_DIR, 'trim_blast_short')
        assert 'usage' in check_output([trim_bin, '-h'])

        fastq_fhand = _make_fhand(FASTQ4)
        result = check_output([trim_bin,
                               '-l', 'AAGCAGTGGTATCAACGCAGAGTACATGGG',
                               '-l', 'AAGCAGTGGTATCAACGCAGAGTACTTTTT',
                               fastq_fhand.name])
        assert '\nTTTTTTTTTTTTTTTTTTTT' in result
        assert '\nCGAGAAGAAGGATCCAAGT' in result

if __name__ == '__main__':
#     import sys; sys.argv = ['', 'TrimTest']
    unittest.main()
