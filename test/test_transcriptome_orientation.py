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
import os
from subprocess import check_output, CalledProcessError

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from crumbs.utils.test_utils import TEST_DATA_DIR
from crumbs.transcript_orientations import TranscriptOrientator
from crumbs.settings import POLYA_ANNOTATOR_MISMATCHES
from crumbs.utils.bin_utils import BIN_DIR
from tempfile import NamedTemporaryFile
from crumbs.seqio import read_seqrecords

# pylint: disable=R0201
# pylint: disable=R0904


class TestTranscriptomeOrientator(unittest.TestCase):
    'the class'

    def test_transcriptome_orientator(self):
        '''tests the orientator class'''
        estscan_matrix = os.path.join(TEST_DATA_DIR,
                                      'Arabidopsis_thaliana.smat')
        seq1 = SeqRecord(seq=Seq('atccgtcagcatcCAATAAAAA'), id='seq1_polia+')
        seq2 = SeqRecord(seq=Seq('TTTTcTTcatccgtcag'), id='seq2_polia-')
        seq3 = SeqRecord(seq=Seq('cTTcatccgtcag'), id='seq3')
        seq_forward = 'CATAGGGTCACCAATGGCTTCTTCTTTGCTTGCACTCTTCTCCTGTCTCTTCCTC'
        seq_forward += 'TCTCTCTTATCTCTCTCCTCCTCCCTAAATCTCCGCCGTCCGATCTTCTCTCAA'
        seq_forward += 'TCCAACGACCTCGATCTCTTCTCTTCTCTAAATCTCGACCGTCCATCTCTCGCC'
        seq_forward += 'GCCGATGACATCCACGATCTTCTCCCACGCTACGGATTCCCGAAAGGTCTTCTT'
        seq_forward += 'CCCAACAACGTCAAATCGTACACTATCTCCGACGACGGCGATTTCACCGTTGAC'
        seq_forward += 'CTGATTTCCAGTTGCTACGTCAAGTTCTCCGATCAACTCGTTTTCTACGGCAAG'
        seq_forward += 'AATATCGCCGGAAAACTCAGTTACGGATCTGTTAAAGACGTCCGTGGAATCCAA'
        seq_forward += 'GCTAAAGAAGCTTTCCTTTGGCTACCAATCACCGCCATGGAATCGGATCCAAGC'
        seq_forward += 'TCTGCCACGGTTGTGTTCTCCGTCGGATTTGTGTCCAAGACTTTACCTGCTTCC'
        seq_forward += 'ATGTTCGAAAATGTTCCTTCTTGCTCAAGAAACCTAAATCTTCAAGACTCTTGA'
        seq_forward += 'ATCCACCTGAAACGATCTCAAGATTCAACATTCCCTCCACCCTTTATAGTTTTG'
        seq_forward += 'TATTTCAGAAGTATTTTGCTTGGTTTCGTAGATATAGGTTCGAATTGGAAAAGA'
        seq_forward += 'TACTATCTTAATTATTCGAATCAGATTATGTTATACTGCCCAAA'

        seq_reverse = 'TTTGGGCAGTATAACATAATCTGATTCGAATAATTAAGATAGTATCTTTTCCAAT'
        seq_reverse += 'TCGAACCTATATCTACGAAACCAAGCAAAATACTTCTGAAATACAAAACTATAA'
        seq_reverse += 'AGGGTGGAGGGAATGTTGAATCTTGAGATCGTTTCAGGTGGATTCAAGAGTCTT'
        seq_reverse += 'GAAGATTTAGGTTTCTTGAGCAAGAAGGAACATTTTCGAACATGGAAGCAGGTA'
        seq_reverse += 'AAGTCTTGGACACAAATCCGACGGAGAACACAACCGTGGCAGAGCTTGGATCCG'
        seq_reverse += 'ATTCCATGGCGGTGATTGGTAGCCAAAGGAAAGCTTCTTTAGCTTGGATTCCAC'
        seq_reverse += 'GGACGTCTTTAACAGATCCGTAACTGAGTTTTCCGGCGATATTCTTGCCGTAGA'
        seq_reverse += 'AAACGAGTTGATCGGAGAACTTGACGTAGCAACTGGAAATCAGGTCAACGGTGA'
        seq_reverse += 'AATCGCCGTCGTCGGAGATAGTGTACGATTTGACGTTGTTGGGAAGAAGACCTT'
        seq_reverse += 'TCGGGAATCCGTAGCGTGGGAGAAGATCGTGGATGTCATCGGCGGCGAGAGATG'
        seq_reverse += 'GACGGTCGAGATTTAGAGAAGAGAAGAGATCGAGGTCGTTGGATTGAGAGAAGA'
        seq_reverse += 'TCGGACGGCGGAGATTTAGGGAGGAGGAGAGAGATAAGAGAGAGAGGAAGAGAC'
        seq_reverse += 'AGGAGAAGAGTGCAAGCAAAGAAGAAGCCATTGGTGACCCTATG'

        seq4 = SeqRecord(seq=Seq(seq_forward), id='seq_orf_forward')
        seq5 = SeqRecord(seq=Seq(seq_reverse), id='seq_orf_reverse')

        seq_forward = 'CTAAATCTCCGCCGTCCGATCTTCTCTCAATCCAACGACCTCGATCTCTTCTCTT'
        seq_forward += 'TCTCCGATCAACTCGTTTTCTACGGCAAGAATATCGCCGGAAAACTCAGTTACG'

        seq_reverse = 'TTTAACAGATCCGTAACTGAGTTTTCCGGCGATATTCTTGCCGTAGAAAACGAGT'
        seq_reverse += 'CGGAGATTTAG'

        seq6 = SeqRecord(seq=Seq(seq_forward), id='seq_blast_forward')
        seq7 = SeqRecord(seq=Seq(seq_reverse), id='seq_blast_reverse')

        seq_forward = 'GTTCGTTTCTCTTCTGAATTTCTGTAATCTGTAACGATGTCTCAGACTACTG'
        seq_forward += 'TCCTCAAGGTTGCTATGTCATGTCAG'

        seq_reverse = 'AGGCAGTCTTCTTCCCAGTTTTCGAAACGGTTTGGAAAACTACATCGC'

        seq8 = SeqRecord(seq=Seq(seq_forward), id='seq_blast2_forward')
        seq9 = SeqRecord(seq=Seq(seq_reverse), id='seq_blast2_reverse')

        seqrecords = [seq1, seq2, seq3, seq4, seq5, seq6, seq7, seq8, seq9]
        estscan_params = {'usage_matrix': estscan_matrix}
        polya_params = {'min_len': 4,
                        'max_cont_mismatches': POLYA_ANNOTATOR_MISMATCHES}
        ara_blastdb = os.path.join(TEST_DATA_DIR, 'blastdbs',
                                   'arabidopsis_genes')
        cala_blastdb = os.path.join(TEST_DATA_DIR, 'blastdbs', 'calabaza')
        filters = [{'kind': 'score_threshold', 'score_key': 'expect',
                    'max_score': 1e-10}]
        blast_params = [{'blastdb': ara_blastdb, 'program': 'blastn',
                         'filters':filters},
                        {'blastdb': cala_blastdb, 'program': 'blastn'}]

        orientator = TranscriptOrientator(polya_params, estscan_params,
                                          blast_params)
        seqrecords = orientator(seqrecords)

        assert str(seq1.seq) == str(seqrecords[0].seq)
        assert str(seq2.seq) == str(seqrecords[1].seq.reverse_complement())
        assert str(seq4.seq) == str(seqrecords[3].seq)
        assert str(seq5.seq) == str(seqrecords[4].seq.reverse_complement())
        assert str(seq6.seq) == str(seqrecords[5].seq)
        assert str(seq7.seq) == str(seqrecords[6].seq.reverse_complement())

    def test_bin_transcrip_orientator(self):
        'it tests the transcript orientator binary'
        orientate_bin = os.path.join(BIN_DIR, 'orientate_transcripts')
        assert 'usage' in check_output([orientate_bin, '-h'])

        in_fpath = os.path.join(TEST_DATA_DIR, 'seqs_to_orientate.fasta')
        estscan_matrix = os.path.join(TEST_DATA_DIR,
                                      'Arabidopsis_thaliana.smat')
        blastdb1 = os.path.join(TEST_DATA_DIR, 'blastdbs', 'arabidopsis_genes')
        blastdb2 = os.path.join(TEST_DATA_DIR, 'blastdbs', 'calabaza')

        out_fhand = NamedTemporaryFile()
        cmd = [orientate_bin, '-u', estscan_matrix, '-d', blastdb1, '-d',
               blastdb2, '-g', 'blastn', '-g', 'blastn', '-v', '0.0001',
               '-v', '0.0001', in_fpath, '-o', out_fhand.name,
               '--polya_min_len', '4']
        check_output(cmd)

        out_seqs = list(read_seqrecords([open(out_fhand.name)]))
        init_seqs = list(read_seqrecords([open(in_fpath)]))

        assert str(init_seqs[0].seq) == str(out_seqs[0].seq)
        out_seq1 = str(out_seqs[1].seq.reverse_complement())
        assert str(init_seqs[1].seq) == out_seq1
        assert 'polyA' in  out_seqs[1].description
        assert str(init_seqs[3].seq) == str(out_seqs[3].seq)
        out_seq4 = str(out_seqs[4].seq.reverse_complement())
        assert str(init_seqs[4].seq) == out_seq4
        assert 'estscan_orf' in  out_seqs[4].description
        assert str(init_seqs[5].seq) == str(out_seqs[5].seq)
        out_seq6 = str(out_seqs[6].seq.reverse_complement())
        assert str(init_seqs[6].seq) == out_seq6
        assert 'blast arabidopsis_genes' in  out_seqs[6].description
        cmd = [orientate_bin, '-u', estscan_matrix, '-d', blastdb1, '-d',
               blastdb2, '-g', 'blastn', '-g', 'blastn', '-v', '0.0001',
               in_fpath]
        stderr = NamedTemporaryFile()
        try:
            check_output(cmd, stderr=stderr)
            self.fail()
        except CalledProcessError:
            stde = open(stderr.name).read()
            assert 'Blast parameters are not well defined' in stde

        # witouth parameters
        out_fhand = NamedTemporaryFile()
        check_output([orientate_bin, in_fpath, '-o', out_fhand.name,
                      '--polya_min_len', '4'])

        out_seqs = list(read_seqrecords([open(out_fhand.name)]))
        init_seqs = list(read_seqrecords([open(in_fpath)]))

        assert str(init_seqs[0].seq) == str(out_seqs[0].seq)
        out_seq1 = str(out_seqs[1].seq.reverse_complement())
        assert str(init_seqs[1].seq) == out_seq1
        assert str(init_seqs[3].seq) == str(out_seqs[3].seq)
        assert str(init_seqs[4].seq) == str(out_seqs[4].seq)
        assert str(init_seqs[5].seq) == str(out_seqs[5].seq)
        assert str(init_seqs[6].seq) == str(out_seqs[6].seq)

        # only with orf annotator
        check_output([orientate_bin, in_fpath, '-o', out_fhand.name, '-u',
                      estscan_matrix, '--polya_min_len', '4'])

        out_seqs = list(read_seqrecords([open(out_fhand.name)]))
        init_seqs = list(read_seqrecords([open(in_fpath)]))

        assert str(init_seqs[0].seq) == str(out_seqs[0].seq)
        out_seq1 = str(out_seqs[1].seq.reverse_complement())
        assert str(init_seqs[1].seq) == out_seq1
        assert str(init_seqs[3].seq) == str(out_seqs[3].seq)
        out_seq4 = str(out_seqs[4].seq.reverse_complement())
        assert str(init_seqs[4].seq) == out_seq4
        assert str(init_seqs[5].seq) == str(out_seqs[5].seq)
        assert str(init_seqs[6].seq) == str(out_seqs[6].seq)

        # multiprocessor
        out_fhand = NamedTemporaryFile()
        cmd = [orientate_bin, '-u', estscan_matrix, '-d', blastdb1, '-d',
               blastdb2, '-g', 'blastn', '-g', 'blastn', '-v', '0.0001',
               '-v', '0.0001', in_fpath, '-o', out_fhand.name, '-p', '2',
               '--polya_min_len', '4']
        check_output(cmd)
        out_seqs = list(read_seqrecords([open(out_fhand.name)]))
        init_seqs = list(read_seqrecords([open(in_fpath)]))

        assert str(init_seqs[0].seq) == str(out_seqs[0].seq)
        out_seq1 = str(out_seqs[1].seq.reverse_complement())
        assert str(init_seqs[1].seq) == out_seq1
        assert 'polyA' in  out_seqs[1].description
        assert str(init_seqs[3].seq) == str(out_seqs[3].seq)
        out_seq4 = str(out_seqs[4].seq.reverse_complement())
        assert str(init_seqs[4].seq) == out_seq4
        assert 'estscan_orf' in  out_seqs[4].description
        assert str(init_seqs[5].seq) == str(out_seqs[5].seq)
        out_seq6 = str(out_seqs[6].seq.reverse_complement())
        assert str(init_seqs[6].seq) == out_seq6
        assert 'blast arabidopsis_genes' in  out_seqs[6].description


if __name__ == "__main__":
#    import sys;sys.argv = ['', 'TestTranscriptomeOrientator.test_bin_transcrip_orientator']
    unittest.main()
