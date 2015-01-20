
import unittest

from StringIO import StringIO

from crumbs.vcf.ab_coding import ABCoder, ENOUGH_SUPPORT, NOT_ENOUGH_SUPPORT

# Method could be a function
# pylint: disable=R0201
# Too many public methods
# pylint: disable=R0904
# Missing docstring
# pylint: disable=C0111


class ABCodingTest(unittest.TestCase):
    VCF_HEADER = '''##fileformat=VCFv4.1
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
    vcf = '''#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4\tS5\tS6
20\t11\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/.\t1/.\t0/0\t0/0\t1/1\t1/1
20\t14\t.\tG\tA\t29\tPASS\tNS=3\tGT\t./.\t./.\t1/1\t0/1\t0/1\t0/1
20\t15\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/0\t1/1\t0/0\t./.\t1/1\t1/1
20\t16\t.\tG\tA\t29\tPASS\tNS=3\tGT\t1/1\t0/0\t1/1\t1/1\t0/0\t0/0
20\t17\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/0\t1/1\t2/2\t2/2\t1/1\t1/1
20\t18\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/0\t1/1\t1/1\t0/0\t0/0\t1/1
'''

    def test_ab_coding(self):
        fhand = StringIO(self.VCF_HEADER + self.vcf)

        coder = ABCoder(fhand, parents_a=['S1'], parents_b=['S2'],
                        threshold=0.9)
        assert coder.offspring == ['S3', 'S4', 'S5', 'S6']
        try:
            list(coder.recode_genotypes())
            self.fail('RuntimeError expected')
        except RuntimeError:
            pass

        fhand = StringIO(self.VCF_HEADER + self.vcf)

        coder = ABCoder(fhand, parents_a=['S1'], parents_b=['S2'],
                        threshold=0.9)
        result = coder.recode_genotypes(samples=coder.offspring)
        string = ''
        for snp, geno in result:
            string += str(snp.POS) + ' '
            string += ','.join(''.join(geno) for geno in geno.values())
            string += '\n'
        assert sum(coder.log.values()) == 6
        assert coder.log[NOT_ENOUGH_SUPPORT] == 2
        assert coder.log[ENOUGH_SUPPORT] == 3
        fhand = StringIO()
        coder.write_log(fhand)
        assert '6 SNPs ' in fhand.getvalue()

if __name__ == "__main__":
#     import sys;sys.argv = ['', 'FilterTest.test_close_to_filter']
    unittest.main()
