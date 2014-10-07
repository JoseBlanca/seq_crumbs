
import os
import unittest
from tempfile import NamedTemporaryFile
from StringIO import StringIO
from subprocess import check_call, CalledProcessError, check_output

from vcf_crumbs.utils.file_utils import BIN_DIR
from vcf_crumbs.snv import VCFReader
from vcf_crumbs.genotype_filters import (LowEvidenceAlleleFilter, RIL_SELF,
                                         prob_aa_given_n_a_reads_hw)

# Method could be a function
# pylint: disable=R0201
# Too many public methods
# pylint: disable=R0904
# Missing docstring
# pylint: disable=C0111

VCF_HEADER = '''##fileformat=VCFv4.1
##fileDate=20090805
##source=freebayes
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
##FORMAT=<ID=AO,Number=A,Type=Integer,Description="Read Depth">
##FORMAT=<ID=RO,Number=1,Type=Integer,Description="Read Depth">
'''
VCF_HEADER2 = '''##fileformat=VCFv4.1
##fileDate=20090805
##source=mysnpprogram
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
##FORMAT=<ID=AO,Number=A,Type=Integer,Description="Read Depth">
##FORMAT=<ID=RO,Number=1,Type=Integer,Description="Read Depth">
'''

class GenotypeFilterTests(unittest.TestCase):
    def test_het_filter(self):
        vcf = '''#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT NA00001 NA00002 NA00003
20\t14370\trs6054257\tG\tA\t29\tPASS\tNS=3;DP=14;AF=0.5;DB;H2\tGT:GQ:DP:HQ\t0|0:48:1:51,51\t1|0:48:8:51,51\t1/1:43:5:.,.
'''
        in_fhand = StringIO(VCF_HEADER + vcf)
        snps = list(VCFReader(in_fhand).parse_snvs())
        exp = [[0, 0], [1, 0], [1, 1]]
        assert [call.int_alleles for call in snps[0].calls] == exp
        res = [call.int_alleles for call in snps[0].remove_gt_from_het_calls().calls]
        assert res == [[0, 0], [], [1, 1]]
        

    def test_het_filter_binary(self):
        vcf = '''#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT NA00001 NA00002 NA00003
20\t14370\trs6054257\tG\tA\t29\tPASS\tNS=3;DP=14;AF=0.5;DB;H2\tGT:GQ:DP:HQ\t0|0:48:1:51,51\t1|0:48:8:51,51\t1/1:43:5:.,.
'''
        in_fhand = NamedTemporaryFile()
        in_fhand.write(VCF_HEADER + vcf)
        in_fhand.flush()
        out_fhand = NamedTemporaryFile()
        binary = os.path.join(BIN_DIR, 'filter_het_genotypes')
        cmd = [binary, in_fhand.name, '-o', out_fhand.name]
        check_call(cmd)
        assert './.:48:8:51,51' in open(out_fhand.name).read()

    def test_low_qual_gt_filter_binary(self):
        vcf = '''#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT NA00001 NA00002 NA00003
20\t14370\trs6054257\tG\tA\t29\tPASS\tNS=3;DP=14;AF=0.5;DB;H2\tGT:GQ:DP:HQ\t0|0:48:1:51,51\t1|0:48:8:51,51\t1/1:43:5:.,.
20\t17330\t.\tT\tA\t3\tq10\tNS=3;DP=11;AF=0.017\tGT:GQ:DP:HQ\t0|0:49:3:58,50\t0|1:3:5:65,3\t0/0:41:3
20\t1110696\trs6040355\tA\tG,T\t67\tPASS\tNS=2;DP=10;AF=0.333,0.667;AA=T;DB\tGT:GQ:DP:HQ\t1|2:21:6:23,27\t2|1:2:0:18,2\t2/2:35:4
20\t1230237\t.\tT\t.\t47\tPASS\tNS=3;DP=13;AA=T\tGT:GQ:DP:HQ\t0|0:54:7:56,60\t0|0:48:4:51,51\t0/0:61:2
20\t1234567\tmicrosat1\tGTC\tG,GTCT\t50\tPASS\tNS=3;DP=9;AA=G\tGT:GQ:DP\t0/1:35:4\t0/2:17:2\t1/1:40:3
20\t1234567\tmicrosat1\tGTC\tG,GTCT\t50\tPASS\tNS=3;DP=9;AA=G\tGT:GQ:DP\t./.:35:4\t0/2:17:2\t1/1:40:3
'''
        in_fhand = NamedTemporaryFile()
        in_fhand.write(VCF_HEADER + vcf)
        in_fhand.flush()
        out_fhand = NamedTemporaryFile()
        binary = os.path.join(BIN_DIR, 'filter_low_qual_genotypes')
        cmd = [binary, in_fhand.name, '-o', out_fhand.name, '-m', '20']
        check_call(cmd)
        assert './.:17:2' in open(out_fhand.name).read()


class LowQualAlleleTest(unittest.TestCase):
    def test_no_geno_no_alle_freq(self):
        vcf = '''#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT 1 2 3 4 5 6
20\t14\t.\tG\tA\t29\tPASS\tNS=3\tGT\t./.\t./.\t./.\t./.\t./.\t./.
20\t14\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/0\t1/1\t1/1\t0/0\t0/0\t0/1'''

        vcf = StringIO(VCF_HEADER + vcf)
        snps = list(VCFReader(vcf).parse_snvs())
        filter_ = LowEvidenceAlleleFilter(0.99)
        snps = [filter_(snp) for snp in snps]
        expected = [False] * 12
        res = [call.called for snp in snps for call in snp.calls]
        assert filter_.log == {'not_enough_individuals': 12, 'tot': 12}
        assert expected == res

    def test_no_allele_depths(self):
        vcf = '''#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT 1 2 3 4 5 6
20\t14\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/0\t1/1\t1/1\t0/0\t0/0\t0/1'''

        vcf = StringIO(VCF_HEADER2 + vcf)
        snps = list(VCFReader(vcf, min_calls_for_pop_stats=4).parse_snvs())
        filter_ = LowEvidenceAlleleFilter(0.99)
        try:
            snps = [filter_(snp) for snp in snps]
            self.fail('RuntimeError expected')
        except RuntimeError:
            pass

    def test_filter_low_alle_evidence_hw(self):
        vcf = '''#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT 1 2 3 4 5 6
20\t14\t.\tG\tA\t29\tPASS\tNS=3\tGT:RO:AO\t0/0:14:0\t1/1:0:15\t1/1:0:1\t0/0:1:0\t0/0:9:0\t0/1:1:1'''

        vcf_f = StringIO(VCF_HEADER + vcf)
        snps = list(VCFReader(vcf_f, min_calls_for_pop_stats=4).parse_snvs())
        filter_ = LowEvidenceAlleleFilter()
        snps = [filter_(snp) for snp in snps]
        assert filter_.log == {'tot': 6, 'not_enough_evidence': 3,
                               'enough_evidence': 2, 'was_het': 1}
        res = [call.call.data.GT for snp in snps for call in snp.calls]
        assert res == ['0/0', '1/1', '1/.', '0/.', '0/.', '0/1']

    def test_filter_low_alle_evidence_ril(self):
        vcf = '''#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT 1 2 3 4 5 6
20\t14\t.\tG\tA\t29\tPASS\tNS=3\tGT:RO:AO\t0/0:14:0\t1/1:0:15\t1/1:0:1\t0/0:1:0\t0/0:9:0\t0/1:1:1'''
        
        # ril n 7
        vcf_f = StringIO(VCF_HEADER + vcf)
        snps = list(VCFReader(vcf_f).parse_snvs())
        kwargs = {'n_generation': 7}
        filter_ = LowEvidenceAlleleFilter(genotypic_freqs_method=RIL_SELF,
                                          genotypic_freqs_kwargs=kwargs)
        snps = [filter_(snp) for snp in snps]
        assert filter_.log == {'tot': 6, 'enough_evidence': 3,
                               'not_enough_evidence': 2, 'was_het': 1}
        res = [call.call.data.GT for snp in snps for call in snp.calls]
        assert res == ['0/0', '1/1', '1/.', '0/.', '0/0', '0/1']

    def test_prob_aa_given_a_reads(self):
        res = prob_aa_given_n_a_reads_hw(30, freq_a_in_pop=0.1)
        self.assertAlmostEqual(res, 0.999999983236, 4)
        res = prob_aa_given_n_a_reads_hw(3, freq_a_in_pop=0.1)
        self.assertAlmostEqual(res, 0.307692307692, 4)
        res = prob_aa_given_n_a_reads_hw(3, freq_a_in_pop=0.99)
        self.assertAlmostEqual(res, 0.997481108312, 4)
        res = prob_aa_given_n_a_reads_hw(1, freq_a_in_pop=0.99)
        self.assertAlmostEqual(res, 0.99, 4)

    def test_low_qual_gt_filter_binary(self):
        vcf = '''#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT 1 2 3 4 5 6
20\t11\t.\tG\tA\t29\tPASS\tNS=3\tGT:RO:AO\t./.\t./.\t./.\t./.\t./.\t./.
20\t14\t.\tG\tA\t29\tPASS\tNS=3\tGT:RO:AO\t0/0:14:0\t1/1:0:15\t1/1:0:1\t0/0:1:0\t0/0:9:0\t0/1:1:1'''

        in_fhand = NamedTemporaryFile()
        in_fhand.write(VCF_HEADER + vcf)
        in_fhand.flush()
        out_fhand = NamedTemporaryFile()
        binary = os.path.join(BIN_DIR, 'filter_low_evidence_alleles')
        cmd = [binary, in_fhand.name, '-o', out_fhand.name, '-c', '2']
        stdout = check_output(cmd)
        assert 'Tot. SNVs' in stdout
        exp = '0/0:14:0\t1/1:0:15\t1/.:0:1\t0/.:1:0\t0/.:9:0\t0/1:1:1'
        assert exp in open(out_fhand.name).read()

    def test_low_qual_gt_filter_binary_ril(self):
        vcf = '''#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT 1 2 3 4 5 6
20\t11\t.\tG\tA\t29\tPASS\tNS=3\tGT:RO:AO\t./.\t./.\t./.\t./.\t./.\t./.
20\t14\t.\tG\tA\t29\tPASS\tNS=3\tGT:RO:AO\t0/0:14:0\t1/1:0:15\t1/1:0:1\t0/0:1:0\t0/0:9:0\t0/1:1:1'''

        in_fhand = NamedTemporaryFile()
        in_fhand.write(VCF_HEADER + vcf)
        in_fhand.flush()
        out_fhand = NamedTemporaryFile()
        binary = os.path.join(BIN_DIR, 'filter_low_evidence_alleles')
        cmd = [binary, in_fhand.name, '-o', out_fhand.name, '-n', '7',
               '-g', 'ril_self']
        stdout = check_output(cmd)
        assert 'Tot. SNVs' in stdout
        exp = '0/0:14:0\t1/1:0:15\t1/.:0:1\t0/.:1:0\t0/0:9:0\t0/1:1:1'
        assert exp in open(out_fhand.name).read()

        binary = os.path.join(BIN_DIR, 'filter_low_evidence_alleles')
        cmd = [binary, in_fhand.name, '-o', out_fhand.name, '-g', 'ril_self']
        err_fhand = NamedTemporaryFile()
        try:
            check_call(cmd, stderr=err_fhand)
            self.fail('Error expected')
        except CalledProcessError:
            pass

if __name__ == "__main__":
    # import sys;sys.argv = ['', 'test_filter_low_alle_evidence_ril']
    unittest.main()
