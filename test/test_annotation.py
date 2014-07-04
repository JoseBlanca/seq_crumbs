'''
Created on 2014 mai 28

@author: peio
'''
import math
import unittest
from os.path import join, dirname
from tempfile import NamedTemporaryFile
from subprocess import check_output

from vcf import Reader

from vcf_crumbs.utils import TEST_DATA_DIR
from vcf_crumbs.annotation import (calculate_maf, count_alleles, CloseToSnv,
                                   HighVariableRegion, CloseToLimit, MafLimit,
                                   CapEnzyme, AminoChangeAnnotator,
                                   AminoSeverityChangeAnnotator,
                                   HeterozigoteInSamples, IsVariableAnnotator)
from vcf_crumbs.statistics import VARSCAN, FREEBAYES, GATK


VCF_PATH = join(TEST_DATA_DIR, 'sample.vcf.gz')
VCF_INDEL_PATH = join(TEST_DATA_DIR, 'sample_indel.vcf.gz')
REF_PATH = join(TEST_DATA_DIR, 'sample_ref.fasta')
VARI_VCF_PATH = join(TEST_DATA_DIR, 'vari_filter.vcf')
GATK_VCF_PATH = join(TEST_DATA_DIR, 'gatk_sample.vcf.gz')
FREEBAYES_VCF_PATH = join(TEST_DATA_DIR, 'freebayes_sample.vcf.gz')
FREEBAYES2_VCF_PATH = join(TEST_DATA_DIR, 'freebayes_sample2.vcf.gz')
FREEBAYES_MULTI_VCF_PATH = join(TEST_DATA_DIR, 'freebayes_multisample2.vcf.gz')
FREEBAYES4_VCF_PATH = join(TEST_DATA_DIR, 'variable_in_sample_5.vcf.gz')
FREEBAYES3_VCF_PATH = join(TEST_DATA_DIR, 'freebayes_sample3.vcf.gz')
FREEBAYES5_VCF_PATH = join(TEST_DATA_DIR, 'freebayes5.vcf.gz')
FREEBAYES6_VCF_PATH = join(TEST_DATA_DIR, 'freebayes6.vcf.gz')
REF_FREEBAYES = join(TEST_DATA_DIR, 'calabaza_selection.fasta')


def floats_are_equal(num1, num2):
    'Given two numbers it returns True if they are similar'
    if num1 == 0.0:
        if num2 == 0.0:
            return True
        else:
            return False
    log1 = math.log(float(num1))
    log2 = math.log(float(num2))
    return abs(log1 - log2) < 0.01


class AnnotateRecordTests(unittest.TestCase):

    def test_mac_calcule(self):
        records = Reader(filename=FREEBAYES_VCF_PATH)
        records.next()
        records.next()
        rec1 = records.next()
        assert floats_are_equal(calculate_maf(rec1, vcf_variant=FREEBAYES),
                                0.75)

    def test_calculate_alleles(self):
        records = Reader(filename=GATK_VCF_PATH)
        rec1 = records.next()
        counts = count_alleles(rec1, vcf_variant=GATK)
        assert counts == {'all': {'T': 5, 'G': 2}}
        for record in Reader(filename=GATK_VCF_PATH):
            if record.CHROM == 'CAUC00157_TC01' and record.POS == 198:
                counts = count_alleles(record, vcf_variant=GATK)
                assert counts == {'all': {'C': 10, 'G': 41}}

        for record in Reader(filename=FREEBAYES_VCF_PATH):
            if record.CHROM == 'CUUC60606_TC01' and record.POS == 341:
                counts = count_alleles(record, vcf_variant=FREEBAYES)
                assert counts == {'all': {'TT': 2, 'TCT': 20, 'CT': 15}}


class FakeClass(object):
    FILTER = []
    INFO = {}

    def add_filter(self, name):
        self.FILTER.append(name)

    def add_info(self, info, value):
        self.INFO[info] = value


class AnnotatorsTest(unittest.TestCase):

    def test_close_to_filter(self):
        records = Reader(filename=FREEBAYES_VCF_PATH)
        records.next()
        rec1 = records.next()
        filter_ = CloseToSnv(distance=300, max_maf=None,
                             vcf_fpath=FREEBAYES_VCF_PATH)
        filter_.vcf_variant = FREEBAYES
        rec1 = filter_(rec1)
        assert filter_.name in rec1.FILTER

        records = Reader(filename=FREEBAYES_VCF_PATH)
        records.next()
        rec1 = records.next()
        filter_ = CloseToSnv(distance=300, max_maf=0.5,
                             vcf_fpath=FREEBAYES_VCF_PATH)
        filter_.vcf_variant = FREEBAYES
        rec1 = filter_(rec1)
        assert rec1.FILTER is None

        records = Reader(filename=FREEBAYES_VCF_PATH)
        records.next()
        rec1 = records.next()
        filter_ = CloseToSnv(distance=300, max_maf=0.8,
                             vcf_fpath=FREEBAYES_VCF_PATH)
        filter_.vcf_variant = FREEBAYES
        rec1 = filter_(rec1)
        assert filter_.name in rec1.FILTER

        assert filter_.name == 'cs300_0.80'
        desc = 'The snv is closer than 300 nucleotides to another snv, '
        desc += 'with maf:0.80'
        assert desc in filter_.description

    def test_high_variable_region_filter(self):
        records = Reader(filename=VCF_PATH)
        rec1 = records.next()
        filter_ = HighVariableRegion(max_variability=0.002,
                                           window=None,
                                           vcf_fpath=VCF_PATH,
                                           ref_fpath=REF_PATH)
        rec1 = filter_(rec1)
        assert filter_.name in rec1.FILTER

        filter_ = HighVariableRegion(max_variability=0.005,
                                           window=None,
                                           vcf_fpath=VCF_PATH,
                                           ref_fpath=REF_PATH)
        rec1 = filter_(rec1)
        assert not filter_.name in rec1.FILTER

        assert  filter_.name == 'hv0.005'
        desc = 'The region has more than 0.5 snvs per 100 bases'
        assert desc in filter_.description

        records = Reader(filename=VCF_PATH)

        records = Reader(filename=VCF_PATH)
        filter_ = HighVariableRegion(max_variability=0.003,
                                           vcf_fpath=VCF_PATH,
                                           ref_fpath=REF_PATH)
        rec1 = filter_(rec1)
        assert filter_.name in rec1.FILTER

        filter_ = HighVariableRegion(max_variability=0.004,
                                           vcf_fpath=VCF_PATH,
                                           ref_fpath=REF_PATH)
        rec1 = filter_(rec1)
        assert not filter_.name in rec1.FILTER

        filter_ = HighVariableRegion(max_variability=0.003,
                                           window=10,
                                           vcf_fpath=VCF_PATH,
                                           ref_fpath=REF_PATH)
        rec1 = filter_(rec1)
        assert filter_.name in rec1.FILTER

        filter_ = HighVariableRegion(max_variability=0.003,
                                           window=100,
                                           vcf_fpath=VCF_PATH,
                                           ref_fpath=REF_PATH)

        rec1 = filter_(rec1)
        assert filter_.name in rec1.FILTER

    def test_close_to_limit_filter(self):
        records = Reader(filename=VCF_PATH)
        rec1 = records.next()
        filter_ = CloseToLimit(distance=60, ref_fpath=REF_PATH)
        rec1 = filter_(rec1)
        assert not filter_.name in rec1.FILTER

        filter_ = CloseToLimit(distance=800, ref_fpath=REF_PATH)
        rec1 = filter_(rec1)
        assert filter_.name in rec1.FILTER

        assert filter_.name == 'cl800'
        desc = 'The snv is closer than 800 nucleotides to the reference edge'
        assert  desc in filter_.description

    def test_maf_limit(self):
        records = Reader(filename=FREEBAYES_VCF_PATH)
        records.next()
        records.next()
        rec1 = records.next()
        filter_ = MafLimit(max_maf=0.8)
        filter_.vcf_variant = FREEBAYES
        rec1 = filter_(rec1)
        assert rec1.FILTER is None
        #assert not filter_.name in rec1.FILTER

        filter_ = MafLimit(max_maf=0.5)
        filter_.vcf_variant = FREEBAYES
        rec1 = filter_(rec1)
        assert filter_.name in rec1.FILTER

        assert  filter_.name == 'maf0.5'
        desc = 'The most frequent allele in samples: all.'
        desc += 'frequency greater than 0.5'
        assert  desc in filter_.description

    def test_cap_enzyme_filter(self):
        seq_str = '>seq1\nATGATGATGgaaattcATGATGATGTGGGAT\n'
        seq_str += '>seq2\nATGATGATGATGATGATGTGGGAT\n'
        fhand = NamedTemporaryFile()
        fhand.write(seq_str)
        fhand.flush()
        filter_ = CapEnzyme(all_enzymes=True, ref_fpath=fhand.name)
        assert filter_.name == 'cet'
        desc = 'SNV is not a CAP detectable by the enzymes: all'
        assert desc in filter_.description

        record = FakeClass()
        record.POS = 11
        record.end = 12
        record.CHROM = 'seq1'
        record.FILTER = []

        alt_allele = FakeClass()
        alt_allele.sequence = 'A'
        record.alleles = ['AA', alt_allele]
        record.FILTER = []
        record.add_filter = lambda x: self.FILTER.append(x)

        rec1 = filter_(record)
        assert not filter_.name in rec1.FILTER

        record = FakeClass()
        record.POS = 12
        record.end = 12
        record.CHROM = 'seq2'
        record.FILTER = []
        record.alleles = ['A', alt_allele]
        rec1 = filter_(record)
        assert filter_.name in rec1.FILTER

    def test_amino_change_filter(self):
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
    ATGTCGTCAACAAAGTGGTTGGTGTGTGCGCTAGTGGTGGTGTGCGTGAGCGTAAGGCAAGCAACATCTGCGC
    CGGCGCCGCAGGAACAAGAATACCCGCCTATGCCCTACGAGTACAAATATGACGTTGAAGATCAAGAGCTTGAA
    GAGAAAGCTCTCTACTTCGGAGCCAACGAAGCAGGAGATGCCCAGGGCAAGGTCATCGGAGGATACCGAGTTCT
    CCTCCCCGATGGTCGTCTTATGACCGTCGAGTACAGTGTGGAGGGAGAAAGCGGTTTCGTTCCCAAAATCACCT
    TCGAAGACAACGCCAGCCCCTTCGGCAAAGGAAAGTAG'''
        ref_fhand = NamedTemporaryFile()
        ref_fhand.write(seq_ref)
        ref_fhand.flush()

        orf_fhand = NamedTemporaryFile()
        orf_fhand.write(orf_seq)
        orf_fhand.flush()

        f = AminoChangeAnnotator(ref_fpath=ref_fhand.name,
                                 orf_seq_fpath=orf_fhand.name)

        record = FakeClass()
        alt_allele = FakeClass()
        alt_allele.sequence = 'C'
        record.CHROM = 'SEUC00016_TC01'
        record.is_indel = False
        record.POS = 112
        record.alleles = ['T', alt_allele]
        record = f(record)
        assert f.name in record.FILTER
        assert record.INFO == {'AAC': 'C->R'}

        f = AminoSeverityChangeAnnotator(ref_fpath=ref_fhand.name,
                                         orf_seq_fpath=orf_fhand.name)
        record = FakeClass()
        alt_allele = FakeClass()
        alt_allele.sequence = 'C'
        record.CHROM = 'SEUC00016_TC01'
        record.is_indel = False
        record.POS = 112
        record.alleles = ['T', alt_allele]
        record = f(record)
        assert f.name in record.FILTER

    def test_is_variable_annotator(self):
        snv = Reader(filename=FREEBAYES6_VCF_PATH).next()
        annotator = IsVariableAnnotator(max_maf=None, samples=['rg1'],
                                        in_union=False,
                                        in_all_groups=True,
                                        reference_free=True,
                                        min_reads=None,
                                        min_reads_per_allele=1,
                                        filter_id=1)
        annotator.vcf_variant = FREEBAYES
        snv = annotator(snv)
        info_id = annotator.info_id
        assert snv.INFO[info_id] == 'True'

        annotator = IsVariableAnnotator(max_maf=None, samples=['rg1'],
                                        in_union=True, in_all_groups=True,
                                        reference_free=True,
                                        min_reads=None, min_reads_per_allele=1,
                                        filter_id=1)

        annotator.vcf_variant = FREEBAYES
        snv = Reader(filename=FREEBAYES6_VCF_PATH).next()
        snv = annotator(snv)
        assert snv.INFO[info_id] == 'True'

        annotator = IsVariableAnnotator(max_maf=None, samples=['rg1'],
                                        in_union=True, in_all_groups=True,
                                        reference_free=True, min_reads=None,
                                        min_reads_per_allele=2, filter_id=1)
        annotator.vcf_variant = FREEBAYES
        snv = Reader(filename=FREEBAYES6_VCF_PATH).next()
        snv = annotator(snv)
        assert snv.INFO[info_id] == 'None'

        annotator = IsVariableAnnotator(max_maf=None, samples=['rg2', 'rg4'],
                             in_union=True, in_all_groups=True,
                             reference_free=True, min_reads=None,
                             min_reads_per_allele=1, filter_id=1)
        annotator.vcf_variant = FREEBAYES
        snv = Reader(filename=FREEBAYES6_VCF_PATH).next()
        snv = annotator(snv)
        assert snv.INFO[info_id] == 'False'

        annotator = IsVariableAnnotator(max_maf=None, samples=['fake'],
                             in_union=True, in_all_groups=True,
                             reference_free=True, min_reads=None,
                             min_reads_per_allele=1, filter_id=1)
        annotator.vcf_variant = FREEBAYES
        snv = Reader(filename=FREEBAYES6_VCF_PATH).next()
        snv = annotator(snv)
        assert snv.INFO[info_id] == 'None'

        annotator = IsVariableAnnotator(max_maf=None, samples=['rg2', 'rg3'],
                             in_union=True, in_all_groups=True,
                             reference_free=True, min_reads=None,
                             min_reads_per_allele=1, filter_id=1)
        annotator.vcf_variant = FREEBAYES
        snv = Reader(filename=FREEBAYES6_VCF_PATH).next()
        snv = annotator(snv)
        assert snv.INFO[info_id] == 'True'

        annotator = IsVariableAnnotator(max_maf=None, samples=['rg5'],
                             in_union=True, in_all_groups=True,
                             reference_free=False, min_reads=None,
                             min_reads_per_allele=1, filter_id=1)
        annotator.vcf_variant = FREEBAYES
        snv = Reader(filename=FREEBAYES6_VCF_PATH).next()
        snv = annotator(snv)
        assert snv.INFO[info_id] == 'None'

        annotator = IsVariableAnnotator(max_maf=None, samples=['rg2'],
                             in_union=True, in_all_groups=True,
                             reference_free=False, min_reads=None,
                             min_reads_per_allele=1, filter_id=1)
        annotator.vcf_variant = FREEBAYES
        snv = Reader(filename=FREEBAYES6_VCF_PATH).next()
        snv = annotator(snv)
        assert snv.INFO[info_id] == 'True'

        annotator = IsVariableAnnotator(max_maf=None, samples=['rg2'],
                             in_union=True, in_all_groups=True,
                             reference_free=True, min_reads=None,
                             min_reads_per_allele=1, filter_id=1)
        annotator.vcf_variant = FREEBAYES
        snv = Reader(filename=FREEBAYES6_VCF_PATH).next()
        snv = annotator(snv)
        assert snv.INFO[info_id] == 'False'

        annotator = IsVariableAnnotator(filter_id=1, samples=['rg1'],
                                        in_union=True, reference_free=True,
                                        max_maf=None, min_reads=1,
                                        min_reads_per_allele=1)
        annotator.vcf_variant = FREEBAYES
        snv = Reader(filename=FREEBAYES6_VCF_PATH).next()
        snv = annotator(snv)
        assert snv.INFO[info_id] == 'True'

        annotator = IsVariableAnnotator(filter_id=1, samples=['rg1'],
                                        in_union=True, reference_free=False,
                                        max_maf=None, min_reads=1,
                                        min_reads_per_allele=1)
        annotator.vcf_variant = FREEBAYES
        snv = Reader(filename=FREEBAYES6_VCF_PATH).next()
        snv = annotator(snv)
        assert snv.INFO[info_id] == 'True'

        annotator = IsVariableAnnotator(filter_id=1, samples=['rg1'],
                                        in_union=False,
                                        reference_free=False, max_maf=None,
                                        min_reads=1, min_reads_per_allele=1)
        annotator.vcf_variant = FREEBAYES
        snv = Reader(filename=FREEBAYES6_VCF_PATH).next()
        snv = annotator(snv)
        assert snv.INFO[info_id] == 'True'

        annotator = IsVariableAnnotator(filter_id=1, samples=['rg2', 'rg4'],
                                        in_union=True, reference_free=True,
                                        max_maf=None, min_reads=1,
                                        min_reads_per_allele=1)
        annotator.vcf_variant = FREEBAYES
        snv = Reader(filename=FREEBAYES6_VCF_PATH).next()
        snv = annotator(snv)
        assert snv.INFO[info_id] == 'False'

        annotator = IsVariableAnnotator(filter_id=1, samples=['rg3', 'rg4'],
                                in_union=True, reference_free=True,
                                max_maf=None, min_reads=1,
                                min_reads_per_allele=1)
        annotator.vcf_variant = FREEBAYES
        snv = Reader(filename=FREEBAYES6_VCF_PATH).next()
        snv = annotator(snv)
        assert snv.INFO[info_id] == 'True'

        annotator = IsVariableAnnotator(filter_id=1, samples=['rg3', 'rg4'],
                                in_union=False, reference_free=True,
                                max_maf=None, min_reads=1,
                                min_reads_per_allele=1)
        annotator.vcf_variant = FREEBAYES
        snv = Reader(filename=FREEBAYES6_VCF_PATH).next()
        snv = annotator(snv)
        assert snv.INFO[info_id] == 'False'

        annotator = IsVariableAnnotator(filter_id=1, samples=['rg5'],
                                in_union=True, reference_free=True,
                                max_maf=None, min_reads=1,
                                min_reads_per_allele=1)
        annotator.vcf_variant = FREEBAYES
        snvs = Reader(filename=FREEBAYES6_VCF_PATH)
        snv = annotator(snvs.next())
        assert snv.INFO[info_id] == 'None'

        annotator = IsVariableAnnotator(filter_id=1, samples=['rg2'],
                                in_union=True, reference_free=True,
                                max_maf=None, min_reads=1,
                                min_reads_per_allele=1)
        annotator.vcf_variant = FREEBAYES
        snvs = Reader(filename=FREEBAYES6_VCF_PATH)
        snvs.next()
        snv = annotator(snvs.next())
        assert snv.INFO[info_id] == 'None'

        annotator = IsVariableAnnotator(filter_id=1, samples=['rg2', 'rg4'],
                                in_union=True, reference_free=True,
                                max_maf=None, min_reads=1,
                                min_reads_per_allele=1)
        annotator.vcf_variant = FREEBAYES
        snvs = Reader(filename=FREEBAYES6_VCF_PATH)
        snvs.next()
        snv = annotator(snvs.next())
        assert snv.INFO[info_id] == 'None'

        snvs = Reader(filename=FREEBAYES5_VCF_PATH)
        snv = snvs.next()
        annotator = IsVariableAnnotator(max_maf=0.7,
                                        samples=['mu16', 'upv196'],
                                        filter_id=1)
        annotator.vcf_variant = FREEBAYES
        snv = annotator(snv)
        info_id = annotator.info_id
        assert snv.INFO[info_id] == 'True'

        # in union False
        snvs = Reader(filename=FREEBAYES5_VCF_PATH)
        snv = snvs.next()
        annotator = IsVariableAnnotator(max_maf=0.7,
                                        samples=['mu16', 'upv196'],
                                        in_union=False, filter_id=1)
        annotator.vcf_variant = FREEBAYES
        snv = annotator(snv)
        info_id = annotator.info_id
        assert snv.INFO[info_id] == 'False'

        # in all_groups False
        snvs = Reader(filename=FREEBAYES5_VCF_PATH)
        snv = snvs.next()
        annotator = IsVariableAnnotator(max_maf=0.7, filter_id=1,
                                        samples=['mu16', 'upv196'],
                                        in_union=False, in_all_groups=False)
        annotator.vcf_variant = FREEBAYES
        snv = annotator(snv)
        info_id = annotator.info_id
        assert snv.INFO[info_id] == 'False'

        annotator = IsVariableAnnotator(filter_id=1, samples=['rg1'],
                                in_union=True, reference_free=True,
                                max_maf=0.95, min_reads=50,
                                min_reads_per_allele=1)
        annotator.vcf_variant = FREEBAYES
        snvs = Reader(filename=FREEBAYES6_VCF_PATH)
        snvs.next()
        snvs.next()
        snv = annotator(snvs.next())
        assert snv.INFO[info_id] == 'True'

        annotator = IsVariableAnnotator(filter_id=1, samples=['rg1'],
                                in_union=True, reference_free=True,
                                max_maf=0.6, min_reads=50,
                                min_reads_per_allele=1)
        annotator.vcf_variant = FREEBAYES
        snvs = Reader(filename=FREEBAYES6_VCF_PATH)
        snvs.next()
        snvs.next()
        snv = annotator(snvs.next())
        assert snv.INFO[info_id] == 'False'

        annotator = IsVariableAnnotator(filter_id=1, samples=['rg1'],
                                in_union=True, reference_free=True,
                                max_maf=0.95, min_reads=200,
                                min_reads_per_allele=1)
        annotator.vcf_variant = FREEBAYES
        snvs = Reader(filename=FREEBAYES6_VCF_PATH)
        snvs.next()
        snvs.next()
        snv = annotator(snvs.next())
        assert snv.INFO[info_id] == 'None'

        annotator = IsVariableAnnotator(filter_id=1, samples=['rg1'],
                                in_union=True, reference_free=True,
                                max_maf=0.6, min_reads=200,
                                min_reads_per_allele=1)
        annotator.vcf_variant = FREEBAYES
        snvs = Reader(filename=FREEBAYES6_VCF_PATH)
        snvs.next()
        snvs.next()
        snv = annotator(snvs.next())
        assert snv.INFO[info_id] == 'None'

        #freebayes
        snvs = Reader(filename=FREEBAYES2_VCF_PATH)
        snv = snvs.next()
        annotator = IsVariableAnnotator(filter_id=1, samples=['sample01_gbs'],
                                in_union=True, reference_free=True,
                                max_maf=0.6, min_reads=3,
                                min_reads_per_allele=1)
        annotator.vcf_variant = FREEBAYES
        snv = annotator(snv)
        assert snv.INFO[info_id] == 'True'

        snvs = Reader(filename=FREEBAYES_MULTI_VCF_PATH)
        snv = snvs.next()
        annotator = IsVariableAnnotator(filter_id=2, samples=['sample06_gbs'])
        annotator.vcf_variant = FREEBAYES
        snv = annotator(snv)
        info_id = annotator.info_id
        assert snv.INFO[info_id] == 'False'


class TestInfoMappers(unittest.TestCase):

    def test_hetegorigot_percent(self):
        reader = Reader(open(FREEBAYES4_VCF_PATH))
        snp = reader.next()
        het_in_samples = HeterozigoteInSamples(filter_id=1)
        het_in_samples.vcf_variant = FREEBAYES
        het_in_samples(snp)
        info_id = het_in_samples.info_id
        assert snp.INFO[info_id] == 'True'

        snp = reader.next()
        het_in_samples = HeterozigoteInSamples(filter_id=1, gq_threshold=30,
                                               min_num_called=3)
        het_in_samples.vcf_variant = FREEBAYES
        het_in_samples(snp)
        info_id = het_in_samples.info_id
        assert snp.INFO[info_id] == 'None'

        reader = Reader(open(FREEBAYES4_VCF_PATH))
        snp = reader.next()
        het_in_samples = HeterozigoteInSamples(filter_id=1, gq_threshold=50,
                                               min_num_called=8)
        het_in_samples.vcf_variant = FREEBAYES
        het_in_samples(snp)
        info_id = het_in_samples.info_id
        assert snp.INFO[info_id] == 'None'

        reader = Reader(open(FREEBAYES4_VCF_PATH))
        het_in_samples = HeterozigoteInSamples(filter_id=1, gq_threshold=30,
                                               min_num_called=3,
                                               min_percent_het_gt=30)
        het_in_samples.vcf_variant = FREEBAYES
        for snp in reader:
            if snp.POS == 272668159 and snp.CHROM == 'Pepper.v.1.55.chr01':
                het_in_samples(snp)
                info_id = het_in_samples.info_id
                assert snp.INFO[info_id] == 'False'
                break

        reader = Reader(open(FREEBAYES4_VCF_PATH))
        het_in_samples = HeterozigoteInSamples(filter_id=1, gq_threshold=30,
                                               min_num_called=3,
                                               min_percent_het_gt=30,
                                               samples=['sample05_gbs',
                                                        'sample06_gbs',
                                                        'sample07_gbs'])
        het_in_samples.vcf_variant = FREEBAYES
        for snp in reader:
            if snp.POS == 228123401 and snp.CHROM == 'Pepper.v.1.55.chr10':
                het_in_samples(snp)
                info_id = het_in_samples.info_id
                assert snp.INFO[info_id] == 'None'
                break


class BinaryTest(unittest.TestCase):
    def test_run_binary(self):
        binary = join(dirname(__file__), '..', 'bin', 'annotate_snvs')
        assert 'usage' in check_output([binary, '-h'])

        config = '''[1]
    [[CloseToSnv]]
        distance = 60
        max_maf = 0.7
[2]
    [[HighVariableRegion]]
        max_variability = 0.05
        ref_fpath = '{sample_fasta}'
[3]
    [[CapEnzyme]]
        all_enzymes = True
        ref_fpath = '{sample_fasta}'

[4]
    [[HeterozigoteInSamples]]
        filter_id = 1
[5]
    [[IsVariableAnnotator]]
        filter_id = 1
        samples = ['mu16']
'''
        config = config.format(sample_fasta=REF_FREEBAYES)

        config_fhand = NamedTemporaryFile(suffix='.config')
        config_fhand.write(config)
        config_fhand.flush()
        cmd = [binary, FREEBAYES3_VCF_PATH, '-f', config_fhand.name]
        #raw_input(' '.join(cmd))
        result = check_output(cmd)
        assert 'cs60_0.70\t' in result
        assert 'CAP=MmeI' in result
        assert 'HIS1=True' in result

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'AnnotatorsTest.test_is_variable_annotator']
    unittest.main()
