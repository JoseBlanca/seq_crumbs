'''
Created on 2014 mai 28

@author: peio
'''
import math
import unittest
from os.path import join, dirname
from tempfile import NamedTemporaryFile
from subprocess import check_output
from StringIO import StringIO

from vcf import Reader

from vcf_crumbs.utils import TEST_DATA_DIR
from vcf_crumbs.annotation import (CloseToSnv, HighVariableRegion,
                                   CloseToLimit, MafDepthLimit, CapEnzyme,
                                   AminoChangeAnnotator, IsVariableAnnotator,
                                   AminoSeverityChangeAnnotator,
                                   HeterozigoteInSamples)
from vcf_crumbs.snv import VCFReader, SNV
from test.test_snv import VCF_HEADER


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


class FakeClass(object):
    filters = []
    infos = {}

    def add_filter(self, name):
        self.filters.append(name)

    def add_info(self, info, value):
        self.infos[info] = value


class AnnotatorsTest(unittest.TestCase):

    def test_close_to_filter(self):
        records = list(VCFReader(open(FREEBAYES_VCF_PATH),
                                 min_calls_for_pop_stats=1).parse_snvs())
        rec1 = records[1]
        filter_ = CloseToSnv(distance=300, max_maf_depth=None)
        rec1 = filter_(rec1)
        assert filter_.name in rec1.filters

        records = list(VCFReader(open(FREEBAYES_VCF_PATH),
                                 min_calls_for_pop_stats=1).parse_snvs())
        rec1 = records[1]
        filter_ = CloseToSnv(distance=300, max_maf_depth=0.5)
        rec1 = filter_(rec1)
        assert rec1.filters is None

        records = list(VCFReader(open(FREEBAYES_VCF_PATH),
                                 min_calls_for_pop_stats=1).parse_snvs())
        rec1 = records[1]

        filter_ = CloseToSnv(distance=300, max_maf_depth=0.8)
        rec1 = filter_(rec1)
        assert filter_.name in rec1.filters

        assert filter_.name == 'cs300_0.80'
        desc = 'The snv is closer than 300 nucleotides to another snv, '
        desc += 'with maf:0.80'
        assert desc in filter_.description

        records = list(VCFReader(open(FREEBAYES_VCF_PATH),
                                 min_calls_for_pop_stats=1).parse_snvs())
        rec1 = records[1]

        filter_ = CloseToSnv(distance=300, max_maf_depth=0.8, snv_type='snp')
        rec1 = filter_(rec1)
        assert filter_.name in rec1.filters

    def test_high_variable_region_filter(self):
        records = list(VCFReader(open(VCF_PATH),
                                 min_calls_for_pop_stats=1).parse_snvs())
        rec1 = records[0]
        filter_ = HighVariableRegion(max_variability=0.002,
                                     window=None, ref_fpath=REF_PATH)
        rec1 = filter_(rec1)
        assert filter_.name in rec1.filters

        filter_ = HighVariableRegion(max_variability=0.005,
                                     window=None, ref_fpath=REF_PATH)
        rec1 = filter_(rec1)
        assert filter_.name not in rec1.filters

        assert filter_.name == 'hv0.005'
        desc = 'The region has more than 0.5 snvs per 100 bases'
        assert desc in filter_.description
        records = list(VCFReader(open(VCF_PATH),
                                 min_calls_for_pop_stats=1).parse_snvs())
        filter_ = HighVariableRegion(max_variability=0.003, ref_fpath=REF_PATH)
        rec1 = filter_(rec1)
        assert filter_.name in rec1.filters

        filter_ = HighVariableRegion(max_variability=0.004, ref_fpath=REF_PATH)
        rec1 = filter_(rec1)
        assert filter_.name not in rec1.filters

        filter_ = HighVariableRegion(max_variability=0.003, window=10,
                                     ref_fpath=REF_PATH)
        rec1 = filter_(rec1)
        assert filter_.name in rec1.filters

        filter_ = HighVariableRegion(max_variability=0.003, window=100,
                                     ref_fpath=REF_PATH)

        rec1 = filter_(rec1)
        assert filter_.name in rec1.filters

    def test_close_to_limit_filter(self):
        records = list(VCFReader(open(VCF_PATH),
                                 min_calls_for_pop_stats=1).parse_snvs())
        rec1 = records[0]
        filter_ = CloseToLimit(distance=60, ref_fpath=REF_PATH)
        rec1 = filter_(rec1)
        assert filter_.name not in rec1.filters

        filter_ = CloseToLimit(distance=800, ref_fpath=REF_PATH)
        rec1 = filter_(rec1)
        assert filter_.name in rec1.filters

        assert filter_.name == 'cl800'
        desc = 'The snv is closer than 800 nucleotides to the reference edge'
        assert desc in filter_.description

    def test_maf_limit(self):
        records = list(VCFReader(open(FREEBAYES_VCF_PATH),
                                 min_calls_for_pop_stats=1).parse_snvs())
        rec1 = records[2]
        filter_ = MafDepthLimit(max_maf=0.8)
        rec1 = filter_(rec1)
        assert rec1.filters is None

        # assert not filter_.name in rec1.FILTER

        filter_ = MafDepthLimit(max_maf=0.5)
        rec1 = filter_(rec1)
        assert filter_.name in rec1.filters

        assert filter_.name == 'maf0.5'
        desc = 'The most frequent allele in samples: all.'
        desc += 'frequency greater than 0.5'
        assert desc in filter_.description

    def test_cap_enzyme_filter(self):
        seq_str = '>seq1\nATGATGATGgaaattcATGATGATGTGGGAT\n'
        seq_str += '>seq2\nATGATGATGATGATGATGTGGGAT\n'
        fhand = NamedTemporaryFile()
        fhand.write(seq_str)
        fhand.flush()
        vcf = '''#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT NA00001 NA00002 NA00003
seq1\t11\trs6054257\tAA\tA\t29\tPASS\tNS=3;DP=14;AF=0.5;DB;H2\tGT:GQ:DP:HQ\t0|0:48:1:51,51\t1|0:48:8:51,51\t1/1:43:5:.,.
seq2\t12\t.\tA\tAA\t3\tq10\tNS=3;DP=11;AF=0.017\tGT:GQ:DP:HQ\t0|0:49:3:58,50\t0|1:3:5:65,3\t0/0:41:3
20\t1110696\trs6040355\tA\tG,T\t67\tPASS\tNS=2;DP=10;AF=0.333,0.667;AA=T;DB\tGT:GQ:DP:HQ\t1|2:21:6:23,27\t2|1:2:0:18,2\t2/2:35:4
20\t1230237\t.\tT\t.\t47\tPASS\tNS=3;DP=13;AA=T\tGT:GQ:DP:HQ\t0|0:54:7:56,60\t0|0:48:4:51,51\t0/0:61:2
20\t1234567\tmicrosat1\tGTC\tG,GTCT\t50\tPASS\tNS=3;DP=9;AA=G\tGT:GQ:DP\t0/1:35:4\t0/2:17:2\t1/1:40:3
20\t1234567\tmicrosat1\tGTC\tG,GTCT\t50\tPASS\tNS=3;DP=9;AA=G\tGT:GQ:DP\t./.:35:4\t0/2:17:2\t1/1:40:3
'''
        vcf = StringIO(VCF_HEADER + vcf)
        snps = list(VCFReader(vcf).parse_snvs())

        filter_ = CapEnzyme(all_enzymes=True, ref_fpath=fhand.name)
        assert filter_.name == 'cet'
        desc = 'SNV is not a CAP detectable by the enzymes: all'
        assert desc in filter_.description

        record = snps[0]
        rec1 = filter_(record)
        assert filter_.name not in rec1.filters

        record = snps[1]
        rec1 = filter_(record)
        assert filter_.name in rec1.filters

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
        vcf = '''#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT NA00001 NA00002 NA00003
SEUC00016_TC01\t112\trs6054257\tT\tC\t29\tPASS\tNS=3;DP=14;AF=0.5;DB;H2\tGT:GQ:DP:HQ\t0|0:48:1:51,51\t1|0:48:8:51,51\t1/1:43:5:.,.
'''
        vcf = StringIO(VCF_HEADER + vcf)
        snps = list(VCFReader(vcf).parse_snvs())
        f = AminoChangeAnnotator(ref_fpath=ref_fhand.name,
                                 orf_seq_fpath=orf_fhand.name)

        snv = snps[0]

        snv = f(snv)
        assert f.name in snv.filters
        assert snv.infos['AAC'] == 'C->R'

        f = AminoSeverityChangeAnnotator(ref_fpath=ref_fhand.name,
                                         orf_seq_fpath=orf_fhand.name)
        record = snps[0]
        record = f(record)
        assert f.name in record.filters

    def test_is_variable_annotator(self):
        records = list(VCFReader(open(FREEBAYES6_VCF_PATH),
                                 min_calls_for_pop_stats=1).parse_snvs())
        snv = records[0]
        annotator = IsVariableAnnotator(max_maf_depth=None, samples=['rg1'],
                                        in_union=False,
                                        in_all_groups=True,
                                        reference_free=True,
                                        min_reads=None,
                                        min_reads_per_allele=1,
                                        filter_id=1)
        snv = annotator(snv)
        info_id = annotator.info_id
        assert snv.infos[info_id] == 'True'

        annotator = IsVariableAnnotator(max_maf_depth=None, samples=['rg1'],
                                        in_union=True, in_all_groups=True,
                                        reference_free=True,
                                        min_reads=None, min_reads_per_allele=1,
                                        filter_id=1)

        records = list(VCFReader(open(FREEBAYES6_VCF_PATH),
                                 min_calls_for_pop_stats=1).parse_snvs())
        snv = records[0]
        snv = annotator(snv)
        assert snv.infos[info_id] == 'True'

        annotator = IsVariableAnnotator(max_maf_depth=None, samples=['rg1'],
                                        in_union=True, in_all_groups=True,
                                        reference_free=True, min_reads=None,
                                        min_reads_per_allele=2, filter_id=1)
        records = list(VCFReader(open(FREEBAYES6_VCF_PATH),
                                 min_calls_for_pop_stats=1).parse_snvs())
        snv = records[0]
        snv = annotator(snv)
        assert snv.infos[info_id] == 'None'

        annotator = IsVariableAnnotator(max_maf_depth=None,
                                        samples=['rg2', 'rg4'],
                                        in_union=True, in_all_groups=True,
                                        reference_free=True, min_reads=None,
                                        min_reads_per_allele=1, filter_id=1)
        records = list(VCFReader(open(FREEBAYES6_VCF_PATH),
                                 min_calls_for_pop_stats=1).parse_snvs())
        snv = records[0]
        snv = annotator(snv)
        assert snv.infos[info_id] == 'False'

        annotator = IsVariableAnnotator(max_maf_depth=None,
                                        samples=['rg2', 'rg3'],
                                        in_union=True, in_all_groups=True,
                                        reference_free=True, min_reads=None,
                                        min_reads_per_allele=1, filter_id=1)
        records = list(VCFReader(open(FREEBAYES6_VCF_PATH),
                                 min_calls_for_pop_stats=1).parse_snvs())
        snv = annotator(records[0])
        assert snv.infos[info_id] == 'True'

        annotator = IsVariableAnnotator(max_maf_depth=None, samples=['rg2'],
                                        in_union=True, in_all_groups=True,
                                        reference_free=False, min_reads=None,
                                        min_reads_per_allele=1, filter_id=1)
        records = list(VCFReader(open(FREEBAYES6_VCF_PATH),
                                 min_calls_for_pop_stats=1).parse_snvs())
        snv = annotator(records[0])
        assert snv.infos[info_id] == 'True'

        annotator = IsVariableAnnotator(max_maf_depth=None, samples=['rg2'],
                                        in_union=True, in_all_groups=True,
                                        reference_free=True, min_reads=None,
                                        min_reads_per_allele=1, filter_id=1)
        records = list(VCFReader(open(FREEBAYES6_VCF_PATH),
                                 min_calls_for_pop_stats=1).parse_snvs())
        snv = annotator(records[0])
        assert snv.infos[info_id] == 'False'

        annotator = IsVariableAnnotator(filter_id=1, samples=['rg1'],
                                        in_union=True, reference_free=True,
                                        max_maf_depth=None, min_reads=1,
                                        min_reads_per_allele=1)
        records = list(VCFReader(open(FREEBAYES6_VCF_PATH),
                                 min_calls_for_pop_stats=1).parse_snvs())
        snv = annotator(records[0])
        assert snv.infos[info_id] == 'True'

        annotator = IsVariableAnnotator(filter_id=1, samples=['rg1'],
                                        in_union=True, reference_free=False,
                                        max_maf_depth=None, min_reads=1,
                                        min_reads_per_allele=1)
        records = list(VCFReader(open(FREEBAYES6_VCF_PATH),
                                 min_calls_for_pop_stats=1).parse_snvs())
        snv = annotator(records[0])
        assert snv.infos[info_id] == 'True'

        annotator = IsVariableAnnotator(filter_id=1, samples=['rg1'],
                                        in_union=False, max_maf_depth=None,
                                        reference_free=False,
                                        min_reads=1, min_reads_per_allele=1)
        records = list(VCFReader(open(FREEBAYES6_VCF_PATH),
                                 min_calls_for_pop_stats=1).parse_snvs())
        snv = annotator(records[0])
        assert snv.infos[info_id] == 'True'

        annotator = IsVariableAnnotator(filter_id=1, samples=['rg2', 'rg4'],
                                        in_union=True, reference_free=True,
                                        max_maf_depth=None, min_reads=1,
                                        min_reads_per_allele=1)
        records = list(VCFReader(open(FREEBAYES6_VCF_PATH),
                                 min_calls_for_pop_stats=1).parse_snvs())
        snv = annotator(records[0])
        assert snv.infos[info_id] == 'False'

        annotator = IsVariableAnnotator(filter_id=1, samples=['rg3', 'rg4'],
                                        in_union=True, reference_free=True,
                                        max_maf_depth=None, min_reads=1,
                                        min_reads_per_allele=1)
        records = list(VCFReader(open(FREEBAYES6_VCF_PATH),
                                 min_calls_for_pop_stats=1).parse_snvs())
        snv = annotator(records[0])
        assert snv.infos[info_id] == 'True'

        annotator = IsVariableAnnotator(filter_id=1, samples=['rg3', 'rg4'],
                                        in_union=False, reference_free=True,
                                        max_maf_depth=None, min_reads=1,
                                        min_reads_per_allele=1)
        records = list(VCFReader(open(FREEBAYES6_VCF_PATH),
                                 min_calls_for_pop_stats=1).parse_snvs())
        snv = annotator(records[0])
        assert snv.infos[info_id] == 'False'

        annotator = IsVariableAnnotator(filter_id=1, samples=['rg2'],
                                        in_union=True, reference_free=True,
                                        max_maf_depth=None, min_reads=1,
                                        min_reads_per_allele=1)
        records = list(VCFReader(open(FREEBAYES6_VCF_PATH),
                                 min_calls_for_pop_stats=1).parse_snvs())
        snv = annotator(records[1])
        assert snv.infos[info_id] == 'None'

        annotator = IsVariableAnnotator(filter_id=1, samples=['rg2', 'rg4'],
                                        in_union=True, reference_free=True,
                                        max_maf_depth=None, min_reads=1,
                                        min_reads_per_allele=1)
        records = list(VCFReader(open(FREEBAYES6_VCF_PATH),
                                 min_calls_for_pop_stats=1).parse_snvs())
        snv = annotator(records[1])
        assert snv.infos[info_id] == 'None'

        annotator = IsVariableAnnotator(max_maf_depth=0.7,
                                        samples=['mu16', 'upv196'],
                                        filter_id=1)
        records = list(VCFReader(open(FREEBAYES5_VCF_PATH),
                                 min_calls_for_pop_stats=1).parse_snvs())
        snv = annotator(records[0])
        info_id = annotator.info_id
        assert snv.infos[info_id] == 'True'

        # in union False
        annotator = IsVariableAnnotator(max_maf_depth=0.7,
                                        samples=['mu16', 'upv196'],
                                        in_union=False, filter_id=1)
        records = list(VCFReader(open(FREEBAYES5_VCF_PATH),
                                 min_calls_for_pop_stats=1).parse_snvs())
        snv = annotator(records[0])
        info_id = annotator.info_id
        assert snv.infos[info_id] == 'False'

        # in all_groups False
        annotator = IsVariableAnnotator(max_maf_depth=0.7, filter_id=1,
                                        samples=['mu16', 'upv196'],
                                        in_union=False, in_all_groups=False)
        records = list(VCFReader(open(FREEBAYES5_VCF_PATH),
                                 min_calls_for_pop_stats=1).parse_snvs())
        snv = annotator(records[0])
        info_id = annotator.info_id
        assert snv.infos[info_id] == 'False'

        annotator = IsVariableAnnotator(filter_id=1, samples=['rg1'],
                                        in_union=True, reference_free=True,
                                        max_maf_depth=0.95, min_reads=50,
                                        min_reads_per_allele=1)
        records = list(VCFReader(open(FREEBAYES6_VCF_PATH),
                                 min_calls_for_pop_stats=1).parse_snvs())
        snv = annotator(records[2])
        assert snv.infos[info_id] == 'True'

        annotator = IsVariableAnnotator(filter_id=1, samples=['rg1'],
                                        in_union=True, reference_free=True,
                                        max_maf_depth=0.6, min_reads=50,
                                        min_reads_per_allele=1)
        records = list(VCFReader(open(FREEBAYES6_VCF_PATH),
                                 min_calls_for_pop_stats=1).parse_snvs())
        snv = annotator(records[2])
        assert snv.infos[info_id] == 'False'

        annotator = IsVariableAnnotator(filter_id=1, samples=['rg1'],
                                        in_union=True, reference_free=True,
                                        max_maf_depth=0.95, min_reads=200,
                                        min_reads_per_allele=1)
        records = list(VCFReader(open(FREEBAYES6_VCF_PATH),
                                 min_calls_for_pop_stats=1).parse_snvs())
        snv = annotator(records[2])
        assert snv.infos[info_id] == 'None'

        annotator = IsVariableAnnotator(filter_id=1, samples=['rg1'],
                                        in_union=True, reference_free=True,
                                        max_maf_depth=0.6, min_reads=200,
                                        min_reads_per_allele=1)
        records = list(VCFReader(open(FREEBAYES6_VCF_PATH),
                                 min_calls_for_pop_stats=1).parse_snvs())
        snv = annotator(records[2])
        assert snv.infos[info_id] == 'None'

        annotator = IsVariableAnnotator(filter_id=1, samples=['sample01_gbs'],
                                        in_union=True, reference_free=True,
                                        max_maf_depth=0.6, min_reads=3,
                                        min_reads_per_allele=1)
        records = list(VCFReader(open(FREEBAYES2_VCF_PATH),
                                 min_calls_for_pop_stats=1).parse_snvs())
        snv = annotator(records[0])
        assert snv.infos[info_id] == 'True'

        annotator = IsVariableAnnotator(filter_id=2, samples=['sample06_gbs'])
        records = list(VCFReader(open(FREEBAYES_MULTI_VCF_PATH),
                                 min_calls_for_pop_stats=1).parse_snvs())
        snv = annotator(records[0])
        info_id = annotator.info_id
        assert snv.infos[info_id] == 'False'


class TestInfoMappers(unittest.TestCase):

    def test_hetegorigot_percent(self):
        het_in_samples = HeterozigoteInSamples(filter_id=1)
        records = list(VCFReader(open(FREEBAYES4_VCF_PATH),
                                 min_calls_for_pop_stats=1).parse_snvs())
        snp = het_in_samples(records[0])
        info_id = het_in_samples.info_id
        assert snp.infos[info_id] == 'True'

        snp = records[1]
        het_in_samples = HeterozigoteInSamples(filter_id=1, gq_threshold=30,
                                               min_num_called=3)
        het_in_samples(snp)
        info_id = het_in_samples.info_id
        assert snp.infos[info_id] == 'None'

        het_in_samples = HeterozigoteInSamples(filter_id=1, gq_threshold=50,
                                               min_num_called=8)
        records = list(VCFReader(open(FREEBAYES4_VCF_PATH),
                                 min_calls_for_pop_stats=1).parse_snvs())
        snp = het_in_samples(records[0])
        info_id = het_in_samples.info_id
        assert snp.infos[info_id] == 'None'

        het_in_samples = HeterozigoteInSamples(filter_id=1, gq_threshold=30,
                                               min_num_called=3,
                                               min_percent_het_gt=30)
        records = list(VCFReader(open(FREEBAYES4_VCF_PATH),
                                 min_calls_for_pop_stats=1).parse_snvs())
        snp = het_in_samples(records[0])
        for snp in records:
            if snp.pos == 272668159 and snp.chrom == 'Pepper.v.1.55.chr01':
                het_in_samples(snp)
                info_id = het_in_samples.info_id
                assert snp.infos[info_id] == 'False'
                break

        records = list(VCFReader(open(FREEBAYES4_VCF_PATH),
                                 min_calls_for_pop_stats=1).parse_snvs())
        het_in_samples = HeterozigoteInSamples(filter_id=1, gq_threshold=30,
                                               min_num_called=3,
                                               min_percent_het_gt=30,
                                               samples=['sample05_gbs',
                                                        'sample06_gbs',
                                                        'sample07_gbs'])
        for snp in records:
            if snp.pos == 228123401 and snp.chrom == 'Pepper.v.1.55.chr10':
                het_in_samples(snp)
                info_id = het_in_samples.info_id
                assert snp.infos[info_id] == 'None'
                break


class BinaryTest(unittest.TestCase):
    def test_run_binary(self):
        binary = join(dirname(__file__), '..', 'bin', 'annotate_snvs')
        assert 'usage' in check_output([binary, '-h'])

        config = '''[1]
    [[CloseToSnv]]
        distance = 60
        max_maf_depth = 0.7
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
        samples = ['pep']
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
        assert '\tPASS\t' in result

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'TestInfoMappers']
    unittest.main()
