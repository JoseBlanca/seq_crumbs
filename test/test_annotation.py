
import unittest
from os.path import join
from tempfile import NamedTemporaryFile
from subprocess import check_output
from StringIO import StringIO


from vcf_crumbs.utils.file_utils import TEST_DATA_DIR, BIN_DIR
from vcf_crumbs.annotation import (CloseToSnv, HighVariableRegion,
                                   CloseToLimit, MafDepthLimit, CapEnzyme,
                                   AminoChangeAnnotator, IsVariableAnnotator,
                                   IsVariableDepthAnnotator,
                                   AminoSeverityChangeAnnotator,
                                   HeterozigoteInSamples)
from vcf_crumbs.snv import VCFReader
from test.test_snv import VCF_HEADER

# Method could be a function
# pylint: disable=R0201
# Too many public methods
# pylint: disable=R0904
# Missing docstring
# pylint: disable=C0111


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
VCF_HEADER = '''##fileformat=VCFv4.1
##fileDate=20090805
##source=myImputationProgramV3.1
##reference=file:///seq/references/1000GenomesPilot-NCBI36.fasta
##contig=<ID=20,length=62435964,assembly=B36>
##phasing=partial
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">
##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership, build 129">
##INFO=<ID=H2,Number=0,Type=Flag,Description="HapMap2 membership">
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples have data">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype Quality">
'''


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
        rec1 = records[1].copy()
        filter_ = CloseToSnv(distance=300, max_maf_depth=None)
        filter_(rec1)
        assert filter_.name in rec1.filters

        rec1 = records[1].copy()
        filter_ = CloseToSnv(distance=300, max_maf_depth=0.5)
        filter_(rec1)
        assert rec1.filters is None

        rec1 = records[1].copy()
        filter_ = CloseToSnv(distance=300, max_maf_depth=0.8)
        filter_(rec1)
        assert filter_.name in rec1.filters

        assert filter_.name == 'cs300_0.80'
        desc = 'The snv is closer than 300 nucleotides to another snv, '
        desc += 'with maf:0.80'
        assert desc in filter_.description

        rec1 = records[1].copy()
        filter_ = CloseToSnv(distance=300, max_maf_depth=0.8, snv_type='snp')
        filter_(rec1)
        assert filter_.name in rec1.filters

    def test_high_variable_region_filter(self):
        records = list(VCFReader(open(VCF_PATH),
                                 min_calls_for_pop_stats=1).parse_snvs())
        rec1 = records[0].copy()
        filter_ = HighVariableRegion(max_variability=0.002,
                                     window=None, ref_fpath=REF_PATH)
        filter_(rec1)
        assert filter_.name in rec1.filters

        filter_ = HighVariableRegion(max_variability=0.005,
                                     window=None, ref_fpath=REF_PATH)
        rec1 = records[0].copy()
        filter_(rec1)
        assert filter_.name not in rec1.filters

        assert filter_.name == 'hv0.005'
        desc = 'The region has more than 0.5 snvs per 100 bases'
        assert desc in filter_.description
        records = list(VCFReader(open(VCF_PATH),
                                 min_calls_for_pop_stats=1).parse_snvs())
        filter_ = HighVariableRegion(max_variability=0.003, ref_fpath=REF_PATH)
        rec1 = records[0].copy()
        filter_(rec1)
        assert filter_.name in rec1.filters

        filter_ = HighVariableRegion(max_variability=0.004, ref_fpath=REF_PATH)
        rec1 = records[0].copy()
        filter_(rec1)
        assert filter_.name not in rec1.filters

        filter_ = HighVariableRegion(max_variability=0.003, window=10,
                                     ref_fpath=REF_PATH)
        rec1 = records[0].copy()
        filter_(rec1)
        assert filter_.name in rec1.filters

        filter_ = HighVariableRegion(max_variability=0.003, window=100,
                                     ref_fpath=REF_PATH)
        rec1 = records[0].copy()
        filter_(rec1)
        assert filter_.name in rec1.filters

    def test_close_to_limit_filter(self):
        records = list(VCFReader(open(VCF_PATH),
                                 min_calls_for_pop_stats=1).parse_snvs())
        rec1 = records[0].copy()
        filter_ = CloseToLimit(distance=60, ref_fpath=REF_PATH)
        filter_(rec1)
        assert filter_.name not in rec1.filters

        filter_ = CloseToLimit(distance=800, ref_fpath=REF_PATH)
        rec1 = records[0].copy()
        filter_(rec1)
        assert filter_.name in rec1.filters

        assert filter_.name == 'cl800'
        desc = 'The snv is closer than 800 nucleotides to the reference edge'
        assert desc in filter_.description

    def test_maf_limit(self):
        records = list(VCFReader(open(FREEBAYES_VCF_PATH),
                                 min_calls_for_pop_stats=1).parse_snvs())
        rec1 = records[2].copy()
        filter_ = MafDepthLimit(max_maf=0.8)
        filter_(rec1)
        assert rec1.filters is None

        # assert not filter_.name in rec1.FILTER

        filter_ = MafDepthLimit(max_maf=0.5)
        rec1 = records[2].copy()
        filter_(rec1)
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

        rec1 = snps[0].copy()
        filter_(rec1)
        assert filter_.name not in rec1.filters

        rec1 = snps[1].copy()
        filter_(rec1)
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

        snv = snps[0].copy()
        f(snv)
        assert f.name in snv.filters
        assert snv.infos['AAC'] == 'C->R'

        f = AminoSeverityChangeAnnotator(ref_fpath=ref_fhand.name,
                                         orf_seq_fpath=orf_fhand.name)
        record = snps[0].copy()
        f(record)
        assert f.name in record.filters

    def test_is_variable_dp_annotator(self):
        records = list(VCFReader(open(FREEBAYES6_VCF_PATH),
                                 min_calls_for_pop_stats=1).parse_snvs())
        snv = records[0].copy()
        annotator = IsVariableDepthAnnotator(max_maf_depth=None,
                                             samples=['rg1'],
                                             in_union=False,
                                             in_all_groups=True,
                                             reference_free=True,
                                             min_reads=None,
                                             min_reads_per_allele=1,
                                             filter_id=1)
        annotator(snv)
        info_id = annotator.info_id
        assert snv.infos[info_id] == 'True'

        annotator = IsVariableDepthAnnotator(max_maf_depth=None,
                                             samples=['rg1'],
                                             in_union=True,
                                             in_all_groups=True,
                                             reference_free=True,
                                             min_reads=None,
                                             min_reads_per_allele=1,
                                             filter_id=1)

        snv = records[0].copy()
        annotator(snv)
        assert snv.infos[info_id] == 'True'

        annotator = IsVariableDepthAnnotator(max_maf_depth=None,
                                             samples=['rg1'],
                                             in_union=True,
                                             in_all_groups=True,
                                             reference_free=True,
                                             min_reads=None,
                                             min_reads_per_allele=2,
                                             filter_id=1)
        snv = records[0].copy()
        annotator(snv)
        assert snv.infos[info_id] == 'None'

        annotator = IsVariableDepthAnnotator(max_maf_depth=None,
                                             samples=['rg2', 'rg4'],
                                             in_union=True,
                                             in_all_groups=True,
                                             reference_free=True,
                                             min_reads=None,
                                             min_reads_per_allele=1,
                                             filter_id=1)
        snv = records[0].copy()
        annotator(snv)
        assert snv.infos[info_id] == 'False'

        annotator = IsVariableDepthAnnotator(max_maf_depth=None,
                                             samples=['rg2', 'rg3'],
                                             in_union=True,
                                             in_all_groups=True,
                                             reference_free=True,
                                             min_reads=None,
                                             min_reads_per_allele=1,
                                             filter_id=1)
        snv = records[0].copy()
        annotator(snv)
        assert snv.infos[info_id] == 'True'

        annotator = IsVariableDepthAnnotator(max_maf_depth=None,
                                             samples=['rg2'],
                                             in_union=True,
                                             in_all_groups=True,
                                             reference_free=False,
                                             min_reads=None,
                                             min_reads_per_allele=1,
                                             filter_id=1)
        snv = records[0].copy()
        annotator(snv)
        assert snv.infos[info_id] == 'True'

        annotator = IsVariableDepthAnnotator(max_maf_depth=None,
                                             samples=['rg2'],
                                             in_union=True,
                                             in_all_groups=True,
                                             reference_free=True,
                                             min_reads=None,
                                             min_reads_per_allele=1,
                                             filter_id=1)
        snv = records[0].copy()
        annotator(snv)
        assert snv.infos[info_id] == 'False'

        annotator = IsVariableDepthAnnotator(filter_id=1, samples=['rg1'],
                                             in_union=True,
                                             reference_free=True,
                                             max_maf_depth=None, min_reads=1,
                                             min_reads_per_allele=1)
        snv = records[0].copy()
        annotator(snv)
        assert snv.infos[info_id] == 'True'

        annotator = IsVariableDepthAnnotator(filter_id=1, samples=['rg1'],
                                             in_union=True,
                                             reference_free=False,
                                             max_maf_depth=None, min_reads=1,
                                             min_reads_per_allele=1)
        snv = records[0].copy()
        annotator(snv)
        assert snv.infos[info_id] == 'True'

        annotator = IsVariableDepthAnnotator(filter_id=1, samples=['rg1'],
                                             in_union=False,
                                             max_maf_depth=None,
                                             reference_free=False,
                                             min_reads=1,
                                             min_reads_per_allele=1)
        snv = records[0].copy()
        annotator(snv)
        assert snv.infos[info_id] == 'True'

        annotator = IsVariableDepthAnnotator(filter_id=1,
                                             samples=['rg2', 'rg4'],
                                             in_union=True,
                                             reference_free=True,
                                             max_maf_depth=None, min_reads=1,
                                             min_reads_per_allele=1)
        snv = records[0].copy()
        annotator(snv)
        assert snv.infos[info_id] == 'False'

        annotator = IsVariableDepthAnnotator(filter_id=1,
                                             samples=['rg3', 'rg4'],
                                             in_union=True,
                                             reference_free=True,
                                             max_maf_depth=None, min_reads=1,
                                             min_reads_per_allele=1)
        snv = records[0].copy()
        annotator(snv)
        assert snv.infos[info_id] == 'True'

        annotator = IsVariableDepthAnnotator(filter_id=1,
                                             samples=['rg3', 'rg4'],
                                             in_union=False,
                                             reference_free=True,
                                             max_maf_depth=None, min_reads=1,
                                             min_reads_per_allele=1)
        snv = records[0].copy()
        annotator(snv)
        assert snv.infos[info_id] == 'False'

        annotator = IsVariableDepthAnnotator(filter_id=1, samples=['rg2'],
                                             in_union=True,
                                             reference_free=True,
                                             max_maf_depth=None, min_reads=1,
                                             min_reads_per_allele=1)
        snv = records[1].copy()
        annotator(snv)
        assert snv.infos[info_id] == 'None'

        annotator = IsVariableDepthAnnotator(filter_id=1,
                                             samples=['rg2', 'rg4'],
                                             in_union=True,
                                             reference_free=True,
                                             max_maf_depth=None, min_reads=1,
                                             min_reads_per_allele=1)
        snv = records[1].copy()
        annotator(snv)
        assert snv.infos[info_id] == 'None'

        annotator = IsVariableDepthAnnotator(max_maf_depth=0.7,
                                             samples=['mu16', 'upv196'],
                                             filter_id=1)
        records = list(VCFReader(open(FREEBAYES5_VCF_PATH),
                                 min_calls_for_pop_stats=1).parse_snvs())
        snv = records[0].copy()
        annotator(snv)
        info_id = annotator.info_id
        assert snv.infos[info_id] == 'True'

        # in union False
        annotator = IsVariableDepthAnnotator(max_maf_depth=0.7,
                                             samples=['mu16', 'upv196'],
                                             in_union=False, filter_id=1)
        snv = records[0].copy()
        annotator(snv)
        info_id = annotator.info_id
        assert snv.infos[info_id] == 'False'

        # in all_groups False
        annotator = IsVariableDepthAnnotator(max_maf_depth=0.7, filter_id=1,
                                             samples=['mu16', 'upv196'],
                                             in_union=False,
                                             in_all_groups=False)
        snv = records[0].copy()
        annotator(snv)
        info_id = annotator.info_id
        assert snv.infos[info_id] == 'False'

        annotator = IsVariableDepthAnnotator(filter_id=1, samples=['rg1'],
                                             in_union=True,
                                             reference_free=True,
                                             max_maf_depth=0.95, min_reads=50,
                                             min_reads_per_allele=1)
        records = list(VCFReader(open(FREEBAYES6_VCF_PATH),
                                 min_calls_for_pop_stats=1).parse_snvs())
        snv = records[2].copy()
        annotator(snv)
        assert snv.infos[info_id] == 'True'

        annotator = IsVariableDepthAnnotator(filter_id=1, samples=['rg1'],
                                             in_union=True,
                                             reference_free=True,
                                             max_maf_depth=0.6, min_reads=50,
                                             min_reads_per_allele=1)
        snv = records[2].copy()
        annotator(snv)
        assert snv.infos[info_id] == 'False'

        annotator = IsVariableDepthAnnotator(filter_id=1, samples=['rg1'],
                                             in_union=True,
                                             reference_free=True,
                                             max_maf_depth=0.95, min_reads=200,
                                             min_reads_per_allele=1)
        snv = records[2].copy()
        annotator(snv)
        assert snv.infos[info_id] == 'None'

        annotator = IsVariableDepthAnnotator(filter_id=1, samples=['rg1'],
                                             in_union=True,
                                             reference_free=True,
                                             max_maf_depth=0.6, min_reads=200,
                                             min_reads_per_allele=1)
        snv = records[2].copy()
        annotator(snv)
        assert snv.infos[info_id] == 'None'

        annotator = IsVariableDepthAnnotator(filter_id=1,
                                             samples=['sample01_gbs'],
                                             in_union=True,
                                             reference_free=True,
                                             max_maf_depth=0.6, min_reads=3,
                                             min_reads_per_allele=1)
        records = list(VCFReader(open(FREEBAYES2_VCF_PATH),
                                 min_calls_for_pop_stats=1).parse_snvs())
        snv = records[0].copy()
        annotator(snv)
        assert snv.infos[info_id] == 'True'

        annotator = IsVariableDepthAnnotator(filter_id=2,
                                             samples=['sample06_gbs'])
        records = list(VCFReader(open(FREEBAYES_MULTI_VCF_PATH),
                                 min_calls_for_pop_stats=1).parse_snvs())
        snv = records[0].copy()
        annotator(snv)
        info_id = annotator.info_id
        assert snv.infos[info_id] == 'False'

    def test_is_variable_annotator(self):
        
        vcf = '''#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT 1 2 3 4 5 6 7
20\t14\t.\tA\tA\t29\tPASS\t.\tGT\t./.\t./.\t1/1\t0/0\t1/.\t1/1\t1/1'''
        vcf = VCF_HEADER + vcf
        fhand = StringIO(vcf)
        snvs = list(VCFReader(fhand).parse_snvs())
        snv = snvs[0]

        annotator = IsVariableAnnotator(samples=['1', '2'])

        snvc = snv.copy()
        annotator(snvc)
        assert snvc.infos['IV0'] == 'None'

        snvc = snv.copy()
        annotator = IsVariableAnnotator(filter_id=0)
        annotator(snvc)
        assert snvc.infos['IV0'] == 'True'

        snvc = snv.copy()
        annotator = IsVariableAnnotator(filter_id=0, samples=['6', '7'])
        annotator(snvc)
        assert snvc.infos['IV0'] == 'False'

        # min num samples for not variable
        snvc = snv.copy()
        annotator = IsVariableAnnotator(filter_id=0, samples=['6', '7'],
                                        min_samples_for_non_var=3)
        annotator(snvc)
        assert snvc.infos['IV0'] == 'None'

        # with reference as one allele
        snvc = snv.copy()
        annotator = IsVariableAnnotator(filter_id=0, samples=['6', '7'],
                                        consider_reference=True)
        annotator(snvc)
        assert snvc.infos['IV0'] == 'True'

        snvc = snv.copy()
        annotator = IsVariableAnnotator(filter_id=0, samples=['1', '2'],
                                        consider_reference=True)
        annotator(snvc)
        assert snvc.infos['IV0'] == 'None'


class TestInfoMappers(unittest.TestCase):

    def test_hetegorigot_percent(self):
        het_in_samples = HeterozigoteInSamples(filter_id=1)
        records = list(VCFReader(open(FREEBAYES4_VCF_PATH),
                                 min_calls_for_pop_stats=1).parse_snvs())
        snp = records[0].copy()
        het_in_samples(snp)
        info_id = het_in_samples.info_id
        assert snp.infos[info_id] == 'True'

        snp = records[1].copy()
        het_in_samples = HeterozigoteInSamples(filter_id=1, gq_threshold=30,
                                               min_num_called=3)
        het_in_samples(snp)
        info_id = het_in_samples.info_id
        assert snp.infos[info_id] == 'None'

        het_in_samples = HeterozigoteInSamples(filter_id=1, gq_threshold=50,
                                               min_num_called=8)
        snp = records[0].copy()
        het_in_samples(snp)
        info_id = het_in_samples.info_id
        assert snp.infos[info_id] == 'None'

        het_in_samples = HeterozigoteInSamples(filter_id=1, gq_threshold=30,
                                               min_num_called=3,
                                               min_percent_het_gt=30)
        snp = records[0].copy()
        het_in_samples(snp)
        for snp in records:
            if snp.pos == 272668159 and snp.chrom == 'Pepper.v.1.55.chr01':
                snp = snp.copy()
                het_in_samples(snp)
                info_id = het_in_samples.info_id
                assert snp.infos[info_id] == 'False'
                break

        het_in_samples = HeterozigoteInSamples(filter_id=1, gq_threshold=30,
                                               min_num_called=3,
                                               min_percent_het_gt=30,
                                               samples=['sample05_gbs',
                                                        'sample06_gbs',
                                                        'sample07_gbs'])
        for snp in records:
            if snp.pos == 228123401 and snp.chrom == 'Pepper.v.1.55.chr10':
                snp = snp.copy()
                het_in_samples(snp)
                info_id = het_in_samples.info_id
                assert snp.infos[info_id] == 'None'
                break


class BinaryTest(unittest.TestCase):
    def test_run_binary(self):
        binary = join(BIN_DIR, 'annotate_snvs')
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
    [[IsVariableDepthAnnotator]]
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
    # import sys;sys.argv = ['', 'AnnotatorsTest.test_is_variable_annotator']
    unittest.main()
