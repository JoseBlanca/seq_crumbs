
import unittest
from StringIO import StringIO
from os.path import join

from vcf_crumbs.snv import SNV, VCFReader, FREEBAYES, VARSCAN, GATK
from vcf_crumbs.utils import TEST_DATA_DIR

# Method could be a function
# pylint: disable=R0201
# Too many public methods
# pylint: disable=R0904
# Missing docstring
# pylint: disable=C0111


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


class SNVTests(unittest.TestCase):
    def test_init(self):
        vcf = '''#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT NA00001 NA00002 NA00003
20\t14370\trs6054257\tG\tA\t29\tPASS\tNS=3;DP=14;AF=0.5;DB;H2\tGT:GQ:DP:HQ\t0|0:48:1:51,51\t1|0:48:8:51,51\t1/1:43:5:.,.
20\t17330\t.\tT\tA\t3\tq10\tNS=3;DP=11;AF=0.017\tGT:GQ:DP:HQ\t0|0:49:3:58,50\t0|1:3:5:65,3\t0/0:41:3
20\t1110696\trs6040355\tA\tG,T\t67\tPASS\tNS=2;DP=10;AF=0.333,0.667;AA=T;DB\tGT:GQ:DP:HQ\t1|2:21:6:23,27\t2|1:2:0:18,2\t2/2:35:4
20\t1230237\t.\tT\t.\t47\tPASS\tNS=3;DP=13;AA=T\tGT:GQ:DP:HQ\t0|0:54:7:56,60\t0|0:48:4:51,51\t0/0:61:2
20\t1234567\tmicrosat1\tGTC\tG,GTCT\t50\tPASS\tNS=3;DP=9;AA=G\tGT:GQ:DP\t0/1:35:4\t0/2:17:2\t1/1:40:3
20\t1234567\tmicrosat1\tGTC\tG,GTCT\t50\tPASS\tNS=3;DP=9;AA=G\tGT:GQ:DP\t./.:35:4\t0/2:17:2\t1/1:40:3
'''
        vcf = StringIO(VCF_HEADER + vcf)
        snps = list(VCFReader(vcf).parse_snps())
        assert len(snps) == 6
        assert snps[0].pos == 14370
        assert snps[1].is_snp
        assert snps[1].num_called == 3
        assert [call.depth for call in snps[2].calls] == [6, 0, 4]
        assert snps[5].call_rate - 0.6666 < 0.0001
        assert [snp.num_called for snp in snps] == [3, 3, 3, 3, 3, 2]
        assert snp.alleles == ['GTC', 'G', 'GTCT']

    def test_heterozygosity(self):
        # 0/0 1/0 0/0
        vcf = '''#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT NA00001 NA00002 NA00003
20\t14370\t.\tG\tA\t29\tPASS\tNS=3\tGT:GQ:DP:HQ\t0|0:48:1:51,51\t1|0:48:8:51,51\t1/1:43:5:.,.
20\t14370\t.\tG\tA\t29\tPASS\tNS=3\tGT:GQ:DP:HQ\t0|0:48:1:51,51\t0|0:48:8:51,51\t1/1:43:5:.,.'''
        vcf = StringIO(VCF_HEADER + vcf)
        snps = list(VCFReader(vcf).parse_snps())
        snp = snps[1]
        assert snp.obs_het is None
        assert snp.exp_het is None

        snp.min_calls_for_pop_stats = 3
        assert abs(snp.obs_het) < 0.001
        assert abs(snp.exp_het - 0.44444) < 0.001
        assert abs(snp.inbreed_coef - 1.0000) < 0.001

    def test_allele_depths(self):
        vcf = open(join(TEST_DATA_DIR, 'freebayes_al_depth.vcf'))
        snps = list(VCFReader(vcf).parse_snps())
        snp = snps[0]
        result = [None, None, (1, 0), None, None, (0, 1)]
        for sample, res in zip(snp.calls, result):
            if res is None:
                assert sample.ref_depth is None
                assert not sample.allele_depths
            else:
                assert sample.ref_depth == res[0]
                assert sample.allele_depths[1] == res[1]

    def test_mafs(self):
        vcf = open(join(TEST_DATA_DIR, 'freebayes_al_depth.vcf'))
        snps = list(VCFReader(vcf).parse_snps())
        assert snps[0].maf_depth - 0.5 < 0.001
        assert snps[0].allele_depths == {0: 1, 1: 1}
        assert snps[0].depth == 2
        assert snps[1].maf_depth - 1.0 < 0.001
        assert snps[1].allele_depths == {0: 2, 1: 0}
        assert snps[4].maf_depth - 0.9890 < 0.001
        assert snps[4].allele_depths == {0: 90, 1: 1}
        assert snps[4].depth == 91

        result = [1, 1, 1, 1, 1, 0.944444]
        for call, res in zip(snps[4].calls, result):
            assert call.maf_depth - res < 0.001
        assert snps[0].mac

        snps[0].min_calls_for_pop_stats = 3
        assert snps[0].maf is None
        snps[3].min_calls_for_pop_stats = 3
        assert snps[3].maf - 0.75 < 0.0001
        snps[4].min_calls_for_pop_stats = 3
        assert snps[4].maf - 1.0 < 0.0001
        assert snps[0].mac == 2

        # varscan
        varscan_fhand = open(join(TEST_DATA_DIR, 'sample.vcf.gz'))
        reader = VCFReader(fhand=varscan_fhand)
        snp = list(reader.parse_snps())[0]
        snp.min_calls_for_pop_stats = 1
        assert snp.maf_depth is None

        # gatk
        fhand = open(join(TEST_DATA_DIR, 'gatk_sample.vcf.gz'))
        reader = VCFReader(fhand=fhand)
        snp = list(reader.parse_snps())[0]
        assert 0.7 < snp.maf_depth < 0.72
        assert 0.7 < snp.get_call('hib_amarillo').maf_depth < 0.72

        # freebayes
        fhand = open(join(TEST_DATA_DIR, 'freebayes_sample.vcf.gz'))
        reader = VCFReader(fhand=fhand)
        snp = list(reader.parse_snps())[0]
        assert 0.99 < snp.maf_depth < 1.01
        assert 0.99 < snp.get_call('pep').maf_depth < 1.01

    def test_modify_call(self):
        vcf = '''#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT NA00001 NA00002 NA00003
20\t14370\trs6054257\tG\tA\t29\tPASS\tNS=3;DP=14;AF=0.5;DB;H2\tGT:GQ:DP:HQ\t0|0:48:1:51,51\t1|0:48:8:51,51\t1/1:43:5:.,.
20\t17330\t.\tT\tA\t3\tq10\tNS=3;DP=11;AF=0.017\tGT:GQ:DP:HQ\t0|0:49:3:58,50\t0|1:3:5:65,3\t0/0:41:3
20\t1110696\trs6040355\tA\tG,T\t67\tPASS\tNS=2;DP=10;AF=0.333,0.667;AA=T;DB\tGT:GQ:DP:HQ\t1|2:21:6:23,27\t2|1:2:0:18,2\t2/2:35:4
20\t1230237\t.\tT\t.\t47\tPASS\tNS=3;DP=13;AA=T\tGT:GQ:DP:HQ\t0|0:54:7:56,60\t0|0:48:4:51,51\t0/0:61:2
20\t1234567\tmicrosat1\tGTC\tG,GTCT\t50\tPASS\tNS=3;DP=9;AA=G\tGT:GQ:DP\t0/1:35:4\t0/2:17:2\t1/1:40:3
20\t1234567\tmicrosat1\tGTC\tG,GTCT\t50\tPASS\tNS=3;DP=9;AA=G\tGT:GQ:DP\t./.:35:4\t0/2:17:2\t1/1:40:3
'''
        vcf = StringIO(VCF_HEADER + vcf)
        snps = list(VCFReader(vcf).parse_snps())
        call0 = snps[0].calls[0]

        assert call0.called
        call_mod = call0.copy_setting_gt_to_none()
        assert not call_mod.called


class ReaderTest(unittest.TestCase):
    def test_get_snpcaller(self):
        varscan = open(join(TEST_DATA_DIR, 'sample.vcf.gz'))
        gatk = open(join(TEST_DATA_DIR, 'gatk_sample.vcf.gz'))
        freebayes = open(join(TEST_DATA_DIR, 'freebayes_sample.vcf.gz'))
        assert VCFReader(fhand=varscan).snpcaller == VARSCAN
        assert VCFReader(fhand=gatk).snpcaller == GATK
        assert VCFReader(fhand=freebayes).snpcaller == FREEBAYES

    def test_samples(self):
        freebayes = open(join(TEST_DATA_DIR, 'freebayes_sample.vcf.gz'))
        assert VCFReader(fhand=freebayes).samples == ['pep']

if __name__ == "__main__":
    # import sys;sys.argv = ['', 'SNVTests.test_allele_depths']
    unittest.main()
