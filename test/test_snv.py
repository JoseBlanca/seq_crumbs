
import unittest
from StringIO import StringIO
from os.path import join

from vcf_crumbs.snv import (VCFReader, FREEBAYES, VARSCAN, GATK, VCFWriter,
                            GENERIC)
from vcf_crumbs.utils.file_utils import TEST_DATA_DIR
from tempfile import NamedTemporaryFile

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
        snps = list(VCFReader(vcf).parse_snvs())
        assert len(snps) == 6
        assert snps[0].pos == 14369
        assert snps[1].is_snp
        assert snps[1].num_called == 3
        assert [call.depth for call in snps[2].calls] == [6, 0, 4]
        assert snps[5].call_rate - 0.6666 < 0.0001
        assert [snp.num_called for snp in snps] == [3, 3, 3, 3, 3, 2]
        assert snps[5].alleles == ['GTC', 'G', 'GTCT']
        assert snps[5].filters == []
        snps[5].filters = None
        assert snps[5].filters is None

    def test_heterozygosity(self):
        # 0/0 1/0 0/0
        vcf = '''#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT NA00001 NA00002 NA00003
20\t14370\t.\tG\tA\t29\tPASS\tNS=3\tGT:GQ:DP:HQ\t0|0:48:1:51,51\t1|0:48:8:51,51\t1/1:43:5:.,.
20\t14370\t.\tG\tA\t29\tPASS\tNS=3\tGT:GQ:DP:HQ\t0|0:48:1:51,51\t0|0:48:8:51,51\t1/1:43:5:.,.'''
        vcf = StringIO(VCF_HEADER + vcf)
        snps = list(VCFReader(vcf).parse_snvs())
        snp = snps[1]
        assert snp.obs_het is None
        assert snp.exp_het is None

        snp.min_calls_for_pop_stats = 3
        assert abs(snp.obs_het) < 0.001
        assert abs(snp.exp_het - 0.44444) < 0.001
        assert abs(snp.inbreed_coef - 1.0000) < 0.001

    def test_allele_depths(self):
        vcf = open(join(TEST_DATA_DIR, 'freebayes_al_depth.vcf'))
        snps = list(VCFReader(vcf).parse_snvs())
        snp = snps[0]
        result = [None, None, (1, 0), None, None, (0, 1)]
        for sample, res in zip(snp.calls, result):
            if res is None:
                assert sample.ref_depth is None
                assert not sample.allele_depths
            else:
                assert sample.ref_depth == res[0]
                assert sample.allele_depths[1] == res[1]

    def test_vcf_only_with_gt(self):
        vcf = open(join(TEST_DATA_DIR, 'generic.vcf.gz'))
        snvs = list(VCFReader(vcf).parse_snvs())
        snv = snvs[0]
        for call in snv.calls:
            assert call.depth is None
            assert call.gt_qual is None
            assert call.ref_depth is None
            assert call.alt_sum_depths is None
            assert call.allele_depths == {}
            assert call.has_alternative_counts is None
            assert call.maf_depth is None
        assert snv.depth is None
        assert snv.maf_depth is None
        assert snv.allele_depths is None

    def test_mafs(self):
        vcf = open(join(TEST_DATA_DIR, 'freebayes_al_depth.vcf'))
        snps = list(VCFReader(vcf).parse_snvs())
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
        snp = list(reader.parse_snvs())[0]
        snp.min_calls_for_pop_stats = 1
        assert snp.maf_depth is None

        # gatk
        fhand = open(join(TEST_DATA_DIR, 'gatk_sample.vcf.gz'))
        reader = VCFReader(fhand=fhand)
        snp = list(reader.parse_snvs())[0]
        assert 0.7 < snp.maf_depth < 0.72
        assert 0.7 < snp.get_call('hib_amarillo').maf_depth < 0.72

        # freebayes
        fhand = open(join(TEST_DATA_DIR, 'freebayes_sample.vcf.gz'))
        reader = VCFReader(fhand=fhand)
        snp = list(reader.parse_snvs())[0]
        assert 0.99 < snp.maf_depth < 1.01
        assert 0.99 < snp.get_call('pep').maf_depth < 1.01

    def test_modify_calls(self):
        vcf = '''#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT NA00001 NA00002 NA00003
20\t14370\trs6054257\tG\tA\t29\tPASS\tNS=3;DP=14;AF=0.5;DB;H2\tGT:GQ:DP:HQ\t0|0:48:1:51,51\t1|0:48:8:51,51\t1/1:43:5:.,.
20\t17330\t.\tT\tA\t3\tq10\tNS=3;DP=11;AF=0.017\tGT:GQ:DP:HQ\t0|0:49:3:58,50\t0|1:3:5:65,3\t0/0:41:3
20\t1110696\trs6040355\tA\tG,T\t67\tPASS\tNS=2;DP=10;AF=0.333,0.667;AA=T;DB\tGT:GQ:DP:HQ\t1|2:21:6:23,27\t2|1:2:0:18,2\t2/2:35:4
20\t1230237\t.\tT\t.\t47\tPASS\tNS=3;DP=13;AA=T\tGT:GQ:DP:HQ\t0|0:54:7:56,60\t0|0:48:4:51,51\t0/0:61:2
20\t1234567\tmicrosat1\tGTC\tG,GTCT\t50\tPASS\tNS=3;DP=9;AA=G\tGT:GQ:DP\t0/1:35:4\t0/2:17:2\t1/1:40:3
20\t1234567\tmicrosat1\tGTC\tG,GTCT\t50\tPASS\tNS=3;DP=9;AA=G\tGT:GQ:DP\t./.:35:4\t0/2:17:2\t1/1:40:3
'''
        vcf = StringIO(VCF_HEADER + vcf)
        snps = list(VCFReader(vcf).parse_snvs())
        call0 = snps[0].calls[0]

        assert call0.called
        mod_call = call0.copy_setting_gt_to_none()
        assert not mod_call.called

        assert snps[0].get_call('NA00003').called
        mod_snp = snps[0].remove_gt_from_low_qual_calls(45)
        assert not mod_snp.get_call('NA00003').called

    def test_filter_calls(self):
        vcf = '''#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT NA00001 NA00002 NA00003
20\t14370\trs6054257\tG\tA\t29\tPASS\tNS=3;DP=14;AF=0.5;DB;H2\tGT:GQ:DP:HQ\t0|0:48:1:51,51\t1|0:48:8:51,51\t1/1:43:5:.,.
20\t17330\t.\tT\tA\t3\tq10\tNS=3;DP=11;AF=0.017\tGT:GQ:DP:HQ\t0|0:49:3:58,50\t0|1:3:5:65,3\t0/0:41:3
20\t1110696\trs6040355\tA\tG,T\t67\tPASS\tNS=2;DP=10;AF=0.333,0.667;AA=T;DB\tGT:GQ:DP:HQ\t1|2:21:6:23,27\t2|1:2:0:18,2\t2/2:35:4
20\t1230237\t.\tT\t.\t47\tPASS\tNS=3;DP=13;AA=T\tGT:GQ:DP:HQ\t0|0:54:7:56,60\t0|0:48:4:51,51\t0/0:61:2
20\t1234567\tmicrosat1\tGTC\tG,GTCT\t50\tPASS\tNS=3;DP=9;AA=G\tGT:GQ:DP\t0/1:35:4\t0/2:17:2\t1/1:40:3
20\t1234567\tmicrosat1\tGTC\tG,GTCT\t50\tPASS\tNS=3;DP=9;AA=G\tGT:GQ:DP\t./.:35:4\t0/2:17:2\t1/1:40:3
'''
        vcf = StringIO(VCF_HEADER + vcf)
        snps = list(VCFReader(vcf).parse_snvs())
        snp = snps[4]
        assert len(snp.alleles) == 3
        snp_filtered = snp.filter_calls_by_sample(samples=('NA00003',))
        assert len(snp_filtered.alleles) == 2

        snp = snps[1]
        snp_filtered = snp.filter_calls_by_sample(samples=('NA00003',))
        assert len(snp_filtered.calls) == 1

        snp = snps[1]
        snp_filtered = snp.filter_calls_by_sample(samples=('NA00003',),
                                                  reverse=True)
        assert len(snp_filtered.calls) == 2

        try:
            snp_filtered = snp.filter_calls_by_sample(samples=('NA0003',),
                                                      reverse=True)
            self.fail("KeyError Expected")
        except KeyError:
            pass
        assert len(snp_filtered.calls) == 2


class ReaderTest(unittest.TestCase):
    def get_snv_pos(self, snps):
        pos = []
        for snp in snps:
            pos.append(snp.pos)
        return pos

    def test_get_snpcaller(self):
        varscan = open(join(TEST_DATA_DIR, 'sample.vcf.gz'))
        gatk = open(join(TEST_DATA_DIR, 'gatk_sample.vcf.gz'))
        freebayes = open(join(TEST_DATA_DIR, 'freebayes_sample.vcf.gz'))
        assert VCFReader(fhand=varscan).snpcaller == VARSCAN
        assert VCFReader(fhand=gatk).snpcaller == GATK
        assert VCFReader(fhand=freebayes).snpcaller == FREEBAYES
        tassel = open(join(TEST_DATA_DIR, 'generic.vcf.gz'))
        assert VCFReader(fhand=tassel).snpcaller == GENERIC

    def test_samples(self):
        freebayes = open(join(TEST_DATA_DIR, 'freebayes_sample.vcf.gz'))
        assert VCFReader(fhand=freebayes).samples == ['pep']

    def test_vcf_writer(self):
        varscan = open(join(TEST_DATA_DIR, 'vari_filter.vcf'))
        reader = VCFReader(fhand=varscan)
        out_fhand = NamedTemporaryFile()
        writer = VCFWriter(out_fhand, reader)
        for snv in reader.parse_snvs():
            writer.write_snv(snv)
        writer.flush()
        assert 'CUUC00027_TC01' in open(out_fhand.name).read()
        writer.close()

    def test_sliding_window(self):
        fhand = open(join(TEST_DATA_DIR, 'sample_to_window.vcf.gz'))
        reader = VCFReader(fhand=fhand)
        # snps in this vcf [9, 19, 29, 0, 11, 20]
        windows = list(reader.sliding_windows(size=10, min_num_snps=1))
        assert [snp.pos for snp in windows[0]['snps']] == [9]
        assert [snp.pos for snp in windows[1]['snps']] == [19]
        assert [snp.pos for snp in windows[2]['snps']] == [29]
        assert [snp.pos for snp in windows[3]['snps']] == [0]
        assert [snp.pos for snp in windows[4]['snps']] == [11]

        windows = list(reader.sliding_windows(size=20, min_num_snps=1))
        assert [snp.pos for snp in windows[0]['snps']] == [9, 19]
        assert [snp.pos for snp in windows[1]['snps']] == [0, 11]

        ref = '>CUUC00007_TC01\nCTGATGCTGATCGTGATCGAGTCGTAGTCTAGTCGATGTCGACG\n'
        ref += '>CUUC00029_TC01\nCTGATGCTGATCGTGATCGAGTCGTAGTCTAGTCGATGTCGAA\n'
        fhand = open(join(TEST_DATA_DIR, 'sample_to_window.vcf.gz'))
        reader = VCFReader(fhand=fhand)
        windows = list(reader.sliding_windows(size=10, min_num_snps=1,
                                              ref_fhand=StringIO(ref)))
        assert [snp.pos for snp in windows[0]['snps']] == [9]
        assert [snp.pos for snp in windows[1]['snps']] == [19]
        assert [snp.pos for snp in windows[2]['snps']] == [29]
        assert [snp.pos for snp in windows[3]['snps']] == [0]
        assert [snp.pos for snp in windows[4]['snps']] == [11]
        assert [snp.pos for snp in windows[5]['snps']] == [20]

        # with fasta
        fhand = open(join(TEST_DATA_DIR, 'sample_to_window.vcf.gz'))
        reader = VCFReader(fhand=fhand)
        windows = list(reader.sliding_windows(size=20, min_num_snps=1,
                                              ref_fhand=StringIO(ref)))
        assert [snp.pos for snp in windows[0]['snps']] == [9, 19]
        assert [snp.pos for snp in windows[1]['snps']] == [29]
        assert [snp.pos for snp in windows[2]['snps']] == [0, 11]
        assert [snp.pos for snp in windows[3]['snps']] == [20]

        # we skip windows that have no snps
        fhand = open(join(TEST_DATA_DIR, 'sample_to_window.vcf.gz'))
        reader = VCFReader(fhand=fhand)
        windows = list(reader.sliding_windows(size=5, min_num_snps=1,
                                              ref_fhand=StringIO(ref)))
        assert [snp.pos for snp in windows[0]['snps']] == [9]
        assert [snp.pos for snp in windows[1]['snps']] == [19]
        assert [snp.pos for snp in windows[2]['snps']] == [29]
        assert [snp.pos for snp in windows[3]['snps']] == [0]
        assert [snp.pos for snp in windows[4]['snps']] == [11]
        assert [snp.pos for snp in windows[5]['snps']] == [20]

        # we skip no window
        fhand = open(join(TEST_DATA_DIR, 'sample_to_window.vcf.gz'))
        reader = VCFReader(fhand=fhand)
        windows = list(reader.sliding_windows(size=5, min_num_snps=0,
                                              ref_fhand=StringIO(ref)))
        assert [snp.pos for snp in windows[0]['snps']] == []
        assert [snp.pos for snp in windows[1]['snps']] == [9]
        assert [snp.pos for snp in windows[2]['snps']] == []
        assert [snp.pos for snp in windows[3]['snps']] == [19]

        fhand = open(join(TEST_DATA_DIR, 'sample_to_window.vcf.gz'))
        reader = VCFReader(fhand=fhand)
        windows = list(reader.sliding_windows(size=10, min_num_snps=0,
                                              ref_fhand=StringIO(ref),
                                              step=5))
        assert [snp.pos for snp in windows[0]['snps']] == [9]
        assert [snp.pos for snp in windows[1]['snps']] == [9]
        assert [snp.pos for snp in windows[2]['snps']] == [19]
        assert [snp.pos for snp in windows[3]['snps']] == [19]

    def test_het_unknown(self):
        vcf = '''#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT 1 2 3 4 5 6 7 8
20\t2\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/0\t0/0\t0/0\t0/0\t1/1\t1/1\t./.\t1/.\t
'''
        vcf = StringIO(VCF_HEADER + vcf)
        reader = VCFReader(vcf)
        snps = list(reader.parse_snvs())
        snp = snps[0]
        expected = [[0, 0], [0, 0], [0, 0], [0, 0], [1, 1], [1, 1], [],
                    [1, None]]
        assert [call.int_alleles for call in snps[0].calls] == expected
        assert snp.num_called == 7
        out_fhand = StringIO()
        writer = VCFWriter(out_fhand, reader)
        for snv in snps:
            writer.write_snv(snv)
        assert '1/1\t./.\t1/.' in out_fhand.getvalue()

if __name__ == "__main__":
    # import sys;sys.argv = ['', 'SNVTests.test_allele_depths']
    unittest.main()
