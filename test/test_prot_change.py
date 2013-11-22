import unittest
from cStringIO import StringIO

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from vcf_crumbs.prot_change import (SeqCoords, IsIndelError, OutsideAlignment,
                                    BetweenSegments, get_amino_change)


class FakeClass():
    pass


class ProteinChangeTest(unittest.TestCase):

    def test_seq_coord(self):
        seq1 = 'ATCTAGGCTGCTACGATTAGCTGATCGATGTTATCGTAGATCTAGCTGATCATGCTAGCTGATCG'
        seq2 = 'ATCTAGGCTGCTACGATAGCTGATCGATGTTATCGTAGATCTAGCTGATCATGCAGCTGATCG'
        seq1 = SeqRecord(id='seq1', seq=Seq(seq1))
        seq2 = SeqRecord(id='seq2', seq=Seq(seq2))

        seq_coord = SeqCoords(seq1, seq2)
        assert seq_coord.coord_system['seq1'] == [(0, 15), (17, 54), (56, 64)]
        assert seq_coord.coord_system['seq2'] == [(0, 15), (16, 53), (54, 62)]

        seq1 = 'ATCTAGGCTGCTACGATTAGCTGACGATGTTATCGTAGATCTAGCTGATCATCTAGCTGATCG'
        seq2 = 'ATCTAGGCTCTACGATTAGCTGATCGATGTTATCGTAGATCAGCTGATCATGCTAGCTGATCG'
        seq1 = SeqRecord(id='seq1', seq=Seq(seq1))
        seq2 = SeqRecord(id='seq2', seq=Seq(seq2))
        seq_coord = SeqCoords(seq1, seq2)
        assert seq_coord.coord_system['seq1'] == [(0, 8), (10, 23), (24, 40),
                                                  (42, 51), (52, 62)]
        assert seq_coord.coord_system['seq2'] == [(0, 8), (9, 22), (24, 40),
                                                  (41, 50), (52, 62)]
        try:
            seq_coord.to_seq1_pos(23)
            self.fail('Should not reach here')
        except BetweenSegments:
            pass
        assert seq_coord.to_seq1_pos(62) == 62
        assert seq_coord.to_seq1_pos(10) == 11
        assert seq_coord.to_seq1_pos(47) == 48

        try:
            print seq_coord.to_seq2_pos(9)
            self.fail('Should not reach here')
        except BetweenSegments:
            pass
        assert seq_coord.to_seq2_pos(24) == 24
        assert seq_coord.to_seq2_pos(47) == 46
        try:
            seq_coord.to_seq1_pos(65)
            self.fail('Should not reach here')
        except OutsideAlignment:
            pass

        #slices
        assert seq_coord.to_seq1_slice(2, 4) == (2, 4)

        try:
            seq_coord.to_seq1_slice(5, 10)
            self.fail('Should not reach here')
        except BetweenSegments:
            pass
        try:
            seq_coord.to_seq1_slice(61, 72)
            self.fail()
        except OutsideAlignment:
            pass

        assert seq_coord.to_seq1_slice(42, 45) == (43, 46)

        assert seq_coord.to_seq2_slice(2, 4) == (2, 4)
        try:
            print seq_coord.to_seq2_slice(5, 10)
            self.fail('Should not reach here')
        except BetweenSegments:
            pass
        try:
            seq_coord.to_seq2_slice(61, 72)
            self.fail()
        except OutsideAlignment:
            pass
        assert seq_coord.to_seq2_slice(43, 46) == (42, 45)

        seq1 = 'ATCTAGGCTGCTACGATTAGCTGACGATGTTATCGTAGATCTAGCTGATCATCTAGCTGATC'
        seq2 = 'GAGGCTCTACGATTAGCTGATCGATGTTATC'
        seq1 = SeqRecord(id='seq1', seq=Seq(seq1))
        seq2 = SeqRecord(id='seq2', seq=Seq(seq2))
        seq_coord = SeqCoords(seq1, seq2)
        assert seq_coord.coord_system['seq1'] == [(4, 8), (10, 23), (24, 33)]
        assert seq_coord.coord_system['seq2'] == [(0, 4), (5, 18), (20, 29)]
        assert seq_coord.to_seq1_slice(22, 29) == (26, 33)
        #reverse alignment
        seq1 = 'ATCTAGGCTGCTACGATTAGCTGACGATGTTATCGTAGATCTAGCTGATCATCTAGCTGATC'
        seq2 = 'GGATAACATCGATCAGCTAATCGTAGAGCCT'
        seq1 = SeqRecord(id='seq1', seq=Seq(seq1))
        seq2 = SeqRecord(id='seq2', seq=Seq(seq2))
        seq_coord = SeqCoords(seq1, seq2)
        assert seq_coord.reverse

        assert seq_coord.to_seq2_pos(14) == 20
        assert seq_coord.to_seq2_pos(4) == 29
        assert seq_coord.to_seq2_pos(31) == 2
        try:
            seq_coord.to_seq2_pos(9)
            self.fail()
        except BetweenSegments:
            pass

        assert seq_coord.to_seq1_pos(3) == 30
        try:
            seq_coord.to_seq1_pos(10)
            self.fail()
        except BetweenSegments:
            pass
        assert seq_coord.to_seq1_pos(16) == 18

        assert seq_coord.to_seq1_slice(0, 3) == (30, 33)
        assert seq_coord.to_seq1_slice(0, 7) == (26, 33)

        try:
            seq_coord.to_seq1_slice(5, 12)
            self.fail()
        except BetweenSegments:
            pass
        try:
            seq_coord.to_seq1_slice(24, 26)
            self.fail()
        except BetweenSegments:
            pass
        try:
            seq_coord.to_seq1_slice(35, 45)
            self.fail()
        except OutsideAlignment:
            pass

        assert seq_coord.to_seq2_slice(4, 7) == (26, 29)
        try:
            seq_coord.to_seq2_slice(5, 10)
            self.fail()
        except BetweenSegments:
            pass
        assert seq_coord.to_seq2_slice(25, 30) == (3, 8)
        assert seq_coord.to_seq2_slice(20, 23) == (11, 14)

    def test_get_aminos(self):
        # test using
        seq_ref = """>SEUC00016_TC01
    CACGCTAAACAACGATCATTGTCATCGGTACCGATTGTTACAAGTTGTGTGCAGTGTCGT
    GCTATTTGTGTGTACATTCCTTCTAAGATGTCGTCAACAAAGTGGTTGGTGTGTGCGCTA
    GTGGTGGTGTGCGTGAGCGTAAGGCAAGCAACATCTGCGCCGGCGCCGCAGGAACAAGAA
    TACCCGCCTATGCCCTACGAGTACAAATATGACGTTGAAGATCAAGAGCTTGAAGAGAAA
    GCTCTCTACTTCGGAGCCAACGAAGCAGGAGATGCCCAGGGCAAGGTCATCGGAGGATAC
    CGAGTTCTCCTCCCCGATGGTCGTCTTATGACCGTCGAGTACAGTGTGGAGGGAGAAAGC
    GGTTTCGTTCCCAAAATCACCTTCGAAGACAACGCCAGCCCCTTCGGCAAAGGAAAGTAG
    ACCTTATAACGACGCCTACAAGACTGGTACCGCGATCAATTGATACTAGTTCAATTTGAT
    TTCTGAATTCTATGCCGTAAAACATTTTCTTTTATTAATTATACCGATTTCGATAAATAG
    ACATCTTTACCTACTTAACGAATTTCTCATAGGATTCAGAAGTCGAAACCGAAAAAAGTT
    ACTTCAGTTTTCATTAGATTGTAAATGTGTGTAAATTATTATTATTATTATATCAGGGAT
    CCTTAAGTTGATATTAGTGGTGATATAAACGATATTTATGAACGACAATCAGGTATCGTC
    ACTGGCTTGAGTAATGTTAGAAAAAATATAATTTTACCGAAAGCATTAGTAACTTTTTTC
    ACGATTATAATCTCCCATACATACTGTATACTTACGTTACGTATAATAATTTTGATTGTC
    TTCATAGTGTACTCTATAATATATGTAGGTGTAGGCAAAACTCATTCGCCAATAAGATAA
    TATGTACAGTCAGCGATTTCTAAGATAAATTTGTACCGCAAATATCGAGTTACCGATACT
    GTGATCAATTAGAACG"""
        orf_seq = '''>SEUC00016_TC01_orf_seq start=89 end=421 strand=forward
    ATGTCGTCAACAAAGTGGTTGGTGTGTGCGCTAGTGGTGGTGTGCGTGAGCGTAAGGCAAGCAACATCTGCGCC
    GGCGCCGCAGGAACAAGAATACCCGCCTATGCCCTACGAGTACAAATATGACGTTGAAGATCAAGAGCTTGAAG
    AGAAAGCTCTCTACTTCGGAGCCAACGAAGCAGGAGATGCCCAGGGCAAGGTCATCGGAGGATACCGAGTTCTC
    CTCCCCGATGGTCGTCTTATGACCGTCGAGTACAGTGTGGAGGGAGAAAGCGGTTTCGTTCCCAAAATCACCTT
    CGAAGACAACGCCAGCCCCTTCGGCAAAGGAAAGTAG'''
        seq_ref = StringIO(seq_ref)
        orf_seq = StringIO(orf_seq)
        seq_ref = SeqIO.read(seq_ref, 'fasta')
        orf_seq = SeqIO.read(orf_seq, 'fasta')

        record = FakeClass()
        alt_allele = FakeClass()
        alt_allele.sequence = 'C'

        record.is_indel = False
        record.POS = 112
        record.alleles = ['T', alt_allele]
        aminos = get_amino_change(seq_ref, orf_seq, record)
        assert aminos == {'ref_amino': 'C', 'alt_amino': ['R']}

        record = FakeClass()
        record.is_indel = True
        record.POS = 111
        try:
            aminos = get_amino_change(seq_ref, orf_seq, record)
            raise RuntimeError('We should not reach here')
        except IsIndelError:
            pass

        record = FakeClass()
        alt_allele = FakeClass()
        alt_allele.sequence = 'C'
        alt_allele2 = FakeClass()
        alt_allele2.sequence = 'A'

        record.is_indel = False
        record.POS = 112
        record.alleles = ['T', alt_allele, alt_allele2]
        aminos = get_amino_change(seq_ref, orf_seq, record)
        assert aminos == {'ref_amino': 'C', 'alt_amino': ['R', 'S']}

        ## Outside orf
        seq_ref = 'ATCTAGGCTGCTACGATTAGCTGACGATGTTATCGTAGATCTAGCTGATCATCTAGCT'
        seq_ref += 'GATCG'
        orf_seq = 'AGGCTCTACGATTAGCTGATCGATGTTATC'
        seq_ref = SeqRecord(seq=Seq(seq_ref), id='ref')
        orf_seq = SeqRecord(seq=Seq(orf_seq), id='orf')
        record = FakeClass()
        alt_allele = FakeClass()
        alt_allele.sequence = 'C'
        record.is_indel = False
        record.POS = 35
        record.alleles = ['G', alt_allele]

        try:
            aminos = get_amino_change(seq_ref, orf_seq, record)
            self.fail()
        except OutsideAlignment:
            pass

        ## frameshift
        record = FakeClass()
        alt_allele = FakeClass()
        alt_allele.sequence = 'C'
        record.is_indel = False
        record.POS = 10
        record.alleles = ['G', alt_allele]
        try:
            aminos = get_amino_change(seq_ref, orf_seq, record)
            self.fail()
        except BetweenSegments:
            pass

if __name__ == '__main__':
    unittest.main()
