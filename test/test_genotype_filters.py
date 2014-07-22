
import os
import unittest
from tempfile import NamedTemporaryFile
from StringIO import StringIO

from vcf_crumbs.utils.file_utils import BIN_DIR
from subprocess import check_call
from vcf_crumbs.snv import VCFReader

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

if __name__ == "__main__":
    # import sys;sys.argv = ['', 'SNVTests.test_allele_depths']
    unittest.main()
