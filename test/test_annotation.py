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
import os
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from crumbs.annotation import (EstscanOrfAnnotator, _detect_polya_tail,
                               PolyaAnnotator)
from crumbs.utils.test_utils import TEST_DATA_DIR
from crumbs.seqio import read_seqrecords
from crumbs.utils.tags import FIVE_PRIME, THREE_PRIME


class AnnotationTest(unittest.TestCase):
    'Test annotator classes'

    def test_orf_annotator(self):
        'It tests orf annotator'
        fpath = os.path.join(TEST_DATA_DIR, 'orf_test.fasta')
        estscan_matrix = os.path.join(TEST_DATA_DIR,
                                      'Arabidopsis_thaliana.smat')
        seq_records = list(read_seqrecords([open(fpath)]))
        orf_annotator = EstscanOrfAnnotator(estscan_matrix)
        seq_records = orf_annotator(seq_records)
        orf1 = seq_records[0].features[0]
        orf2 = seq_records[1].features[0]
        assert orf1.strand == 1
        assert orf1.location.start.position == 0
        assert orf1.location.end.position == 541
        assert orf2.strand == -1
        assert orf2.location.start.position == 0
        assert orf2.location.end.position == 541
        assert not seq_records[2].features

    def test_polya_annotator(self):
        'It annotates poly-A or poly-T regions'
        seq1 = SeqRecord(seq=Seq('atccgtcagcatcCAATAAAAA'), id='seq1')
        seq2 = SeqRecord(seq=Seq('TTTTcTTcatccgtcag'), id='seq2')
        seq3 = SeqRecord(seq=Seq('TTTTcTTatccgtcagcatcCAATAAAAA'), id='seq3')
        seqs = [seq1, seq2, seq3]
        annotator = PolyaAnnotator(min_len=5, max_cont_mismatches=0)
        seqs = annotator(seqs)
        polya1 = seqs[0].features[0]
        assert not seqs[1].features
        polya3 = seqs[2].features[0]

        assert polya1.type == 'polyA_sequence'
        assert polya1.location.start.position == 17
        assert polya1.location.end.position == 22
        assert polya1.location.strand == 1

        assert polya3.type == 'polyA_sequence'
        assert polya3.location.start.position == 24
        assert polya3.location.end.position == 29
        assert polya3.location.strand == 1

        seq1 = SeqRecord(seq=Seq('atccgtcagcatcCAATAAAAA'), id='seq1')
        seq2 = SeqRecord(seq=Seq('TTTTcTTcatccgtcag'), id='seq2')
        seq3 = SeqRecord(seq=Seq('TTTTcTTatccgtcagcatcCAATAAAAA'), id='seq3')
        seqs = [seq1, seq2, seq3]
        annotator = PolyaAnnotator(min_len=4, max_cont_mismatches=1)
        seqs = annotator(seqs)
        polya1 = seqs[0].features[0]
        polya2 = seqs[1].features[0]
        polya3 = seqs[2].features[0]

        assert polya1.type == 'polyA_sequence'
        assert polya1.location.start.position == 17
        assert polya1.location.end.position == 22
        assert polya1.location.strand == 1

        assert polya2.type == 'polyA_sequence'
        assert polya2.location.start.position == 0
        assert polya2.location.end.position == 4
        assert polya2.location.strand == -1

        assert polya2.location.start.position == 0
        assert polya2.location.end.position == 4
        assert polya2.location.strand == -1

        assert polya3.location.start.position == 24
        assert polya3.location.end.position == 29
        assert polya3.location.strand == 1

    def test_polya_detection(self):
        'It detects poly-A regions'
        seq = 'CAATAAAAA'
        assert _detect_polya_tail(seq, THREE_PRIME, 5, 0) == (4, 9)
        assert _detect_polya_tail(seq, THREE_PRIME, 2, 1) == (1, 9)
        assert not _detect_polya_tail(seq, THREE_PRIME, 6, 0)
        assert _detect_polya_tail(seq, THREE_PRIME, 2, 2) == (1, 9)

        seq = 'TTTTcTTc'
        assert _detect_polya_tail(seq, FIVE_PRIME, 2, 0) == (0, 4)
        assert _detect_polya_tail(seq, FIVE_PRIME, 2, 1) == (0, 7)
        assert _detect_polya_tail(seq, FIVE_PRIME, 3, 3) == (0, 4)

        seq = 'AAAAA'
        assert _detect_polya_tail(seq, THREE_PRIME, 4, 0) == (0, 5)
        seq = 'AAAAAC'
        assert _detect_polya_tail(seq, THREE_PRIME, 4, 1) == (0, 6)
        seq = 'TTTTT'
        assert _detect_polya_tail(seq, FIVE_PRIME, 2, 0) == (0, 5)


if __name__ == "__main__":
    # import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
