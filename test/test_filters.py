from os.path import join, dirname
import unittest
from tempfile import NamedTemporaryFile
from subprocess import check_output, Popen, PIPE
from StringIO import StringIO

from vcf_crumbs.snv import VCFReader

from vcf_crumbs.filters import (PASSED, FILTERED_OUT, group_in_filter_packets,
                                CallRateFilter, BiallelicFilter, IsSNPFilter,
                                GenotypeQualFilter, ObsHetFilter, MafFilter,
                                filter_snvs)
from vcf_crumbs.utils.file_utils import TEST_DATA_DIR, BIN_DIR

# Method could be a function
# pylint: disable=R0201
# Too many public methods
# pylint: disable=R0904
# Missing docstring
# pylint: disable=C0111

VCF_PATH = join(TEST_DATA_DIR, 'sample.vcf.gz')
VCF_INDEL_PATH = join(TEST_DATA_DIR, 'sample_indel.vcf.gz')
VARI_VCF_PATH = join(TEST_DATA_DIR, 'vari_filter.vcf')


class FiltersTest(unittest.TestCase):
    def test_group_in_filter_packets(self):
        items = list(range(10))
        packets = list(group_in_filter_packets(items, 4))
        res = [{FILTERED_OUT: [], PASSED: (0, 1, 2, 3)},
               {FILTERED_OUT: [], PASSED: (4, 5, 6, 7)},
               {FILTERED_OUT: [], PASSED: (8, 9)}]
        assert res == packets

    @staticmethod
    def eval_prop_in_packet(packet, prop):
        eval_prop = lambda snps: [getattr(snp, prop)for snp in snps]
        packet = {PASSED: eval_prop(packet[PASSED]),
                  FILTERED_OUT: eval_prop(packet[FILTERED_OUT])}
        return packet

    @staticmethod
    def filter_vcf(vcf_fhand, filter_):
        snps = VCFReader(vcf_fhand, min_calls_for_pop_stats=1).parse_snvs()
        packet = list(group_in_filter_packets(snps, 10))[0]
        filtered_packet = filter_(packet)
        return filtered_packet

    def test_missing_genotypes(self):
        packet = self.filter_vcf(open(VCF_PATH),
                                 filter_=CallRateFilter(min_calls=2))
        res = self.eval_prop_in_packet(packet, 'num_called')
        expected = {FILTERED_OUT: [0, 1, 1, 1, 1], PASSED: [2, 2, 2, 2, 2]}
        assert res == expected

        packet = self.filter_vcf(open(VCF_PATH),
                                 filter_=CallRateFilter(min_calls=1))
        res = self.eval_prop_in_packet(packet, 'num_called')
        expected = {FILTERED_OUT: [0], PASSED: [2, 2, 2, 2, 1, 2, 1, 1, 1]}
        assert res == expected

        packet = self.filter_vcf(open(VCF_PATH),
                                 filter_=CallRateFilter(min_calls=1,
                                                        reverse=True))
        res = self.eval_prop_in_packet(packet, 'num_called')
        expected = {PASSED: [0], FILTERED_OUT: [2, 2, 2, 2, 1, 2, 1, 1, 1]}
        assert res == expected

        packet = self.filter_vcf(open(VCF_PATH),
                                 filter_=CallRateFilter(min_calls=3))
        res = self.eval_prop_in_packet(packet, 'num_called')
        expected = {PASSED: [], FILTERED_OUT: [2, 2, 2, 2, 0, 1, 2, 1, 1, 1]}
        assert res == expected

    def test_biallelic(self):
        packet = self.filter_vcf(open(VCF_PATH),
                                 filter_=BiallelicFilter())
        res = self.eval_prop_in_packet(packet, 'alleles')
        assert not(res[FILTERED_OUT])
        assert len(res[PASSED]) == 10

        packet = self.filter_vcf(open(VCF_INDEL_PATH),
                                 filter_=BiallelicFilter())
        res = self.eval_prop_in_packet(packet, 'alleles')
        assert len(res[FILTERED_OUT]) == 1
        assert len(res[PASSED]) == 6

    def test_is_snp(self):
        packet = self.filter_vcf(open(VCF_PATH),
                                 filter_=IsSNPFilter())
        res = self.eval_prop_in_packet(packet, 'is_snp')
        assert not(res[FILTERED_OUT])
        assert res[PASSED] == [True] * 10

        packet = self.filter_vcf(open(VCF_INDEL_PATH),
                                 filter_=IsSNPFilter())
        res = self.eval_prop_in_packet(packet, 'is_snp')
        assert res[FILTERED_OUT] == [False] * 7
        assert not res[PASSED]

    def test_gt_qual(self):
        packet = self.filter_vcf(open(VCF_PATH),
                                 filter_=GenotypeQualFilter(20))
        res = self.eval_prop_in_packet(packet, 'qual')
        assert res[FILTERED_OUT] == [None] * 10
        assert not res[PASSED]

        fpath = join(TEST_DATA_DIR, 'freebayes_al_depth.vcf')
        packet = self.filter_vcf(open(fpath),
                                 filter_=GenotypeQualFilter(1))
        res = self.eval_prop_in_packet(packet, 'qual')
        assert len(res[FILTERED_OUT]) == 1
        assert len(res[PASSED]) == 4

    def test_obs_het(self):
        packet = self.filter_vcf(open(VCF_PATH),
                                 filter_=ObsHetFilter(min_het=0.5))
        res = self.eval_prop_in_packet(packet, 'obs_het')
        assert res[FILTERED_OUT] == [0.0, 0.0, 0.0, None, 0.0, 0.0, 0.0]
        assert res[PASSED] == [0.5, 1.0, 0.5]

        packet = self.filter_vcf(open(VCF_PATH),
                                 filter_=ObsHetFilter(max_het=0.5))
        res = self.eval_prop_in_packet(packet, 'obs_het')
        assert res[FILTERED_OUT] == [None, 1.0]
        assert res[PASSED] == [0.0, 0.0, 0.0, 0.5, 0.5, 0.0, 0.0, 0.0]

        packet = self.filter_vcf(open(VCF_PATH),
                                 filter_=ObsHetFilter(min_het=0.1,
                                                      max_het=0.9))
        res = self.eval_prop_in_packet(packet, 'obs_het')
        assert res[FILTERED_OUT] == [0.0, 0.0, 0.0, None, 1.0, 0.0, 0.0, 0.0]
        assert res[PASSED] == [0.5, 0.5]

    def test_maf(self):
        packet = self.filter_vcf(open(VCF_PATH),
                                 filter_=MafFilter(min_maf=0.6))
        res = self.eval_prop_in_packet(packet, 'maf')
        assert res[FILTERED_OUT] == [0.5, 0.5, 0.5, None, 0.5]
        assert res[PASSED] == [0.75, 0.75, 1.0, 1.0, 1.0]

        packet = self.filter_vcf(open(VCF_PATH),
                                 filter_=MafFilter(max_maf=0.6))
        res = self.eval_prop_in_packet(packet, 'maf')
        assert res[FILTERED_OUT] == [0.75, None, 0.75, 1.0, 1.0, 1.0]
        assert res[PASSED] == [0.5, 0.5, 0.5, 0.5]

        packet = self.filter_vcf(open(VCF_PATH),
                                 filter_=MafFilter(min_maf=0.6, max_maf=0.8))
        res = self.eval_prop_in_packet(packet, 'maf')
        assert res[FILTERED_OUT] == [0.5, 0.5, 0.5, None, 0.5, 1.0, 1.0, 1.0]
        assert res[PASSED] == [0.75, 0.75]

VCF_HEADER = '''##fileformat=VCFv4.1
##fileDate=20090805
##source=myImputationProgramV3.1
##reference=file:///seq/references/1000GenomesPilot-NCBI36.fasta
##contig=<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species="Homo sapiens",taxonomy=x>
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
VCF = '''#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT NA00001 NA00002 NA00003
20\t14370\trs6054257\tG\tA\t29\tPASS\tNS=3;DP=14;AF=0.5;DB;H2\tGT:GQ:DP:HQ\t0|0:48:1:51,51\t1|0:48:8:51,51\t1/1:43:5:.,.
20\t17330\t.\tT\tA\t3\tq10\tNS=3;DP=11;AF=0.017\tGT:GQ:DP:HQ\t0|0:49:3:58,50\t0|1:3:5:65,3\t0/0:41:3
20\t1110696\trs6040355\tA\tG,T\t67\tPASS\tNS=2;DP=10;AF=0.333,0.667;AA=T;DB\tGT:GQ:DP:HQ\t1|2:21:6:23,27\t2|1:2:0:18,2\t2/2:35:4
20\t1230237\t.\tT\t.\t47\tPASS\tNS=3;DP=13;AA=T\tGT:GQ:DP:HQ\t0|0:54:7:56,60\t0|0:48:4:51,51\t0/0:61:2
20\t1234567\tmicrosat1\tGTC\tG,GTCT\t50\tPASS\tNS=3;DP=9;AA=G\tGT:GQ:DP\t0/1:35:4\t0/2:17:2\t1/1:40:3
20\t1234567\tmicrosat1\tGTC\tG,GTCT\t50\tPASS\tNS=3;DP=9;AA=G\tGT:GQ:DP\t./.:35:4\t0/2:17:2\t1/1:40:3
'''


class BinaryFilterTest(unittest.TestCase):

    def get_snv_pos(self, vcf_fhand):
        pos = []
        for line in vcf_fhand:
            if line.startswith('#'):
                continue
            pos.append(int(line.split()[1]))
        return pos

    def test_filter_fhand(self):
        in_fhand = StringIO(VCF_HEADER + VCF)
        out_fhand = StringIO()
        template_fhand = StringIO(VCF_HEADER + VCF)
        filter_snvs(in_fhand, out_fhand, filters=[],
                    template_fhand=template_fhand)
        res = self.get_snv_pos(StringIO(out_fhand.getvalue()))
        in_pos = self.get_snv_pos(StringIO(in_fhand.getvalue()))
        assert in_pos == res

        in_fhand = StringIO(VCF_HEADER + VCF)
        out_fhand = StringIO()
        filtered_fhand = StringIO()
        log_fhand = StringIO()
        template_fhand = StringIO(VCF_HEADER + VCF)
        filter_snvs(in_fhand, out_fhand, filters=[BiallelicFilter()],
                    template_fhand=template_fhand,
                    filtered_fhand=filtered_fhand, log_fhand=log_fhand)
        res = self.get_snv_pos(StringIO(out_fhand.getvalue()))
        filtered = self.get_snv_pos(StringIO(filtered_fhand.getvalue()))
        in_pos = self.get_snv_pos(StringIO(in_fhand.getvalue()))
        assert res == [14370, 17330, 1230237]
        assert filtered == [1110696, 1234567, 1234567]
        assert 'SNVs passsed: 3' in log_fhand.getvalue()

    def test_biallelic_binary(self):
        binary = join(BIN_DIR, 'filter_vcf_by_biallelic')

        assert 'positional' in check_output([binary, '-h'])

        in_fhand = NamedTemporaryFile()
        in_fhand.write(VCF_HEADER + VCF)
        in_fhand.flush()
        out_fhand = NamedTemporaryFile()
        filtered_fhand = NamedTemporaryFile()
        cmd = [binary, '-o', out_fhand.name, '-f', filtered_fhand.name,
               in_fhand.name]
        process = Popen(cmd, stderr=PIPE, stdout=PIPE)
        stderr = process.communicate()[1]
        assert "passsed: 3" in stderr

        res = self.get_snv_pos(open(out_fhand.name))
        filtered = self.get_snv_pos(open(filtered_fhand.name))

        assert res == [14370, 17330, 1230237]
        assert filtered == [1110696, 1234567, 1234567]

        # with stdout
        cmd = [binary, '-f', filtered_fhand.name, in_fhand.name]
        process = Popen(cmd, stderr=PIPE, stdout=PIPE)
        stdout, stderr = process.communicate()
        assert "passsed: 3" in stderr
        res = self.get_snv_pos(StringIO(stdout))
        assert res == [14370, 17330, 1230237]

if __name__ == "__main__":
    # import sys;sys.argv = ['', 'BinaryFilterTest.test_biallelic_binary']
    unittest.main()
