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
                         TrimByQuality, TrimWithBlastShort)
from crumbs.utils.bin_utils import BIN_DIR
from crumbs.utils.tags import TRIMMING_RECOMMENDATIONS, VECTOR
from crumbs.seqio import read_seqrecord_packets

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
        seqs.append(SeqRecord(Seq('aaCTTTC')))
        seqs.append(SeqRecord(Seq('CTTCaa')))
        seqs.append(SeqRecord(Seq('aaCTCaa')))
        seqs.append(SeqRecord(Seq('actg')))
        seqs.append(SeqRecord(Seq('AC')))

        trim_lowercased_seqs = TrimLowercasedLetters()
        trim = TrimOrMask()
        # pylint: disable=W0141

        res = [str(s.seq) for s in trim(trim_lowercased_seqs(seqs))]
        assert res == ['CTTTC', 'CTTC', 'CTC', 'AC']


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
        trim = TrimOrMask()

        trim_edges = TrimEdges(left=1)
        res = [str(s.seq) for s in trim(trim_edges(self._some_seqs()))]
        assert res == ['CCG', 'AACCCGGG']

        trim_edges = TrimEdges(right=1)
        res = [str(s.seq) for s in trim(trim_edges(self._some_seqs()))]
        assert res == ['ACC', 'AAACCCGG']

        trim_edges = TrimEdges(left=1, right=1)
        res = [str(s.seq) for s in trim(trim_edges(self._some_seqs()))]
        assert res == ['CC', 'AACCCGG']

        trim_edges = TrimEdges(left=2, right=2)
        res = [str(s.seq) for s in trim(trim_edges(self._some_seqs()))]
        assert res == ['ACCCG']

        trim_edges = TrimEdges(left=3, right=3)
        res = [str(s.seq) for s in trim(trim_edges(self._some_seqs()))]
        assert res == ['CCC']

        trim = TrimOrMask(mask=True)
        trim_edges = TrimEdges(left=1)
        res = [str(s.seq) for s in trim(trim_edges(self._some_seqs()))]
        assert res == ['aCCG', 'aAACCCGGG']

        trim_edges = TrimEdges(right=1)
        res = [str(s.seq) for s in trim(trim_edges(self._some_seqs()))]
        assert res == ['ACCg', 'AAACCCGGg']

        trim_edges = TrimEdges(left=1, right=1)
        res = [str(s.seq) for s in trim(trim_edges(self._some_seqs()))]
        assert res == ['aCCg', 'aAACCCGGg']

        trim_edges = TrimEdges(left=2, right=2)
        res = [str(s.seq) for s in trim(trim_edges(self._some_seqs()))]
        assert res == ['accg', 'aaACCCGgg']

        trim_edges = TrimEdges(left=3, right=3)
        res = [str(s.seq) for s in trim(trim_edges(self._some_seqs()))]
        assert res == ['accg', 'aaaCCCggg']

        # test overlapping mask
        trim1 = TrimEdges(left=3, right=3)
        trim2 = TrimEdges(left=4, right=4)
        res = [str(s.seq) for s in trim(trim2(trim1(self._some_seqs())))]
        assert res == ['accg', 'aaacCcggg']

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

        trim_rec = seq.annotations[TRIMMING_RECOMMENDATIONS]
        seq_trimmer = TrimOrMask()

        trim_rec['vector'] = [(0, 3), (8, 13)]
        seq.annotations[TRIMMING_RECOMMENDATIONS] = trim_rec
        seqs2 = seq_trimmer([seq])
        assert str(seqs2[0].seq) == 'CTCA'

        trim_rec['vector'] = [(0, 0), (8, 13)]
        seq.annotations[TRIMMING_RECOMMENDATIONS] = trim_rec
        seqs2 = seq_trimmer([seq])
        assert str(seqs2[0].seq) == 'GGTCTCA'

        trim_rec['vector'] = [(0, 1), (8, 12)]
        trim_rec['quality'] = [(1, 8), (13, 17)]
        seq.annotations[TRIMMING_RECOMMENDATIONS] = trim_rec
        seqs2 = seq_trimmer([seq])
        assert not seqs2

        trim_rec['vector'] = [(0, 0), (8, 13)]
        trim_rec['quality'] = []
        seq.annotations[TRIMMING_RECOMMENDATIONS] = trim_rec
        seqs2 = seq_trimmer([seq])
        assert str(seqs2[0].seq) == 'GGTCTCA'
        assert TRIMMING_RECOMMENDATIONS not in seqs2[0].annotations


class TrimByQualityTest(unittest.TestCase):
    'It test the quality trimming'

    def test_quality_trimming(self):
        'It trims the edges'
        trim = TrimOrMask()

        trim_quality = TrimByQuality(window=5, threshold=30)

        seq = SeqRecord(Seq('ACTGCTGCATAAAA'))
        quals = [10, 10, 20, 30, 30, 30, 40, 40, 30, 30, 20, 20, 10, 10]
        seq.letter_annotations['phred_quality'] = quals
        seqs = trim(trim_quality([seq]))
        assert seqs[0].letter_annotations['phred_quality'] == [20, 30, 30, 30,
                                                            40, 40, 30, 30, 20]

        # all bad
        trim_quality = TrimByQuality(window=5, threshold=60)
        seqs = trim(trim_quality([seq]))
        assert not seqs

        # all OK
        trim_quality = TrimByQuality(window=5, threshold=5)
        seqs = trim(trim_quality([seq]))
        assert seqs[0].letter_annotations['phred_quality'] == quals

        quals = [20, 20, 20, 60, 60, 60, 60, 60, 20, 20, 20, 20]
        trim_quality = TrimByQuality(window=5, threshold=50)
        seq = SeqRecord(Seq('ataataataata'))
        seq.letter_annotations['phred_quality'] = quals
        seqs = trim(trim_quality([seq]))
        expected = [20, 60, 60, 60, 60, 60, 20]
        assert seqs[0].letter_annotations['phred_quality'] == expected

        quals = [40, 18, 10, 40, 40, 5, 8, 30, 14, 3, 40, 40, 40, 11, 6, 5, 3,
                 20, 10, 12, 8, 5, 4, 7, 1]
        seq = SeqRecord(Seq('atatatatagatagatagatagatg'))
        seq.letter_annotations['phred_quality'] = quals
        trim_quality = TrimByQuality(window=5, threshold=25)
        seqs = trim(trim_quality([seq]))
        assert seqs[0].letter_annotations['phred_quality'] == [40, 18, 10, 40,
                                                               40]

        quals = [40, 40, 13, 11, 40, 9, 40, 4, 27, 38, 40, 4, 11, 40, 40, 10,
                 10, 21, 3, 40, 9, 9, 12, 10, 9]
        seq = SeqRecord(Seq('atatatatatatatatatatatata'))
        seq.letter_annotations['phred_quality'] = quals
        trim_quality = TrimByQuality(window=5, threshold=25)
        seqs = trim(trim_quality([seq]))
        expected = [40, 4, 27, 38, 40]
        assert seqs[0].letter_annotations['phred_quality'] == expected

        quals = [40, 40, 13, 11, 40, 9, 40, 4, 27, 38, 40, 4, 11, 40, 40, 10,
                 10, 21, 3, 40, 9, 9, 12, 10, 9]
        seq = SeqRecord(Seq('atatatatatatatatatatatata'))
        seq.letter_annotations['phred_quality'] = quals
        trim_quality = TrimByQuality(window=5, threshold=25, trim_left=False)
        seqs = trim(trim_quality([seq]))
        expected = [40, 40, 13, 11, 40, 9, 40, 4, 27, 38, 40]
        assert seqs[0].letter_annotations['phred_quality'] == expected

        quals = [40, 40, 13, 11, 40, 9, 40, 4, 27, 38, 40, 4, 11, 40, 40, 10,
                 10, 21, 3, 40, 9, 9, 12, 10, 9]
        seq = SeqRecord(Seq('atatatatatatatatatatatata'))
        seq.letter_annotations['phred_quality'] = quals
        trim_quality = TrimByQuality(window=5, threshold=25, trim_right=False)
        seqs = trim(trim_quality([seq]))
        expected = [40, 4, 27, 38, 40, 4, 11, 40, 40, 10, 10, 21, 3, 40, 9, 9,
                    12, 10, 9]
        assert seqs[0].letter_annotations['phred_quality'] == expected

        quals = [40, 40, 13, 11, 40, 9, 40, 4, 27, 38, 40, 4, 11, 40, 40, 10,
                 10, 21, 3, 40, 9, 9, 12, 10, 9]
        seq = SeqRecord(Seq('atatatatatatatatatatatata'))
        seq.letter_annotations['phred_quality'] = quals
        trim_quality = TrimByQuality(window=5, threshold=25, trim_right=False,
                                     trim_left=False)
        seqs = trim(trim_quality([seq]))
        expected = quals
        assert seqs[0].letter_annotations['phred_quality'] == expected

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

        adaptors = [SeqRecord(Seq('AAGCAGTGGTATCAACGCAGAGTACATGGG')),
                    SeqRecord(Seq('AAGCAGTGGTATCAACGCAGAGTACTTTTT'))]

        blast_trim = TrimWithBlastShort(oligos=adaptors)
        fhand = StringIO(FASTQ4)
        seq_packets = list(read_seqrecord_packets([fhand]))
        # It should trim the first and the second reads.
        res = [s.annotations.get(TRIMMING_RECOMMENDATIONS, {}).get(VECTOR, [])
                                           for s in blast_trim(seq_packets[0])]
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
#    import sys;sys.argv = ['', 'TrimTest.test_trim_seqs']
    unittest.main()
