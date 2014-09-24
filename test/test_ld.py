'''
Created on 23/09/2014

@author: jose
'''
import unittest
from StringIO import StringIO

from vcf_crumbs.snv import VCFReader 
from vcf_crumbs.ld import (_count_biallelic_haplotypes, calculate_r_sqr,
                           HaploCount, _calculate_r_sqr)

from test_snv import VCF_HEADER

# Method could be a function
# pylint: disable=R0201
# Too many public methods
# pylint: disable=R0904
# Missing docstring
# pylint: disable=C0111


class LDTests(unittest.TestCase):
    def test_count_homo_haplotypes(self):
        vcf = '''#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT S1 S2 S3 S4 S5 S6
20\t14\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/0\t1/1\t1/1\t0/0\t0/0\t0/1
20\t14\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/0\t1/1\t1/1\t0/0\t0/0\t0/1'''

        vcf = StringIO(VCF_HEADER + vcf)
        snps = list(VCFReader(vcf).parse_snvs())
        call1 =  snps[0].record.samples
        call2 =  snps[1].record.samples
        counts = _count_biallelic_haplotypes(call1, call2)
        assert counts.AB == 3
        assert counts.ab == 2
        
        vcf = '''#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT S1 S2 S3 S4 S5 S6 S7 S8
20\t14\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/0\t1/1\t1/1\t0/0\t0/0\t0/1\t1/1\t0/0
20\t14\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/0\t1/1\t1/1\t0/0\t0/0\t0/1\t1/1\t0/0'''

        vcf = StringIO(VCF_HEADER + vcf)
        snps = list(VCFReader(vcf).parse_snvs())
        call1 =  snps[0].record.samples
        call2 =  snps[1].record.samples
        counts = _count_biallelic_haplotypes(call1, call2)
        assert counts.AB == 4
        assert counts.ab == 3

        vcf = '''#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT S1 S2 S3 S4 S5 S6 S7 S8
20\t14\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/0\t1/1\t1/1\t0/0\t0/0\t0/1\t1/1\t0/0
20\t14\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/0\t2/2\t1/1\t0/0\t0/0\t0/1\t1/1\t0/0'''

        vcf = StringIO(VCF_HEADER + vcf)
        snps = list(VCFReader(vcf).parse_snvs())
        call1 =  snps[0].record.samples
        call2 =  snps[1].record.samples
        counts = _count_biallelic_haplotypes(call1, call2)
        assert counts.AB == 4
        assert counts.ab == 3

        r_sqr = calculate_r_sqr(snps[0], snps[1])
        self.assertAlmostEqual(r_sqr,  1)

        vcf = '''#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT S1 S2 S3 S4 S5 S6 S7 S8
20\t14\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/0\t1/1\t1/1\t0/0\t0/0\t0/1\t1/1\t0/0
20\t14\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/0\t0/0\t1/1\t1/1\t0/0\t1/1\t0/0\t1/1'''

        vcf = StringIO(VCF_HEADER + vcf)
        snps = list(VCFReader(vcf).parse_snvs())
        call1 =  snps[0].record.samples
        call2 =  snps[1].record.samples
        counts = _count_biallelic_haplotypes(call1, call2)

        r_sqr = calculate_r_sqr(snps[0], snps[1])
        assert r_sqr - 1.0 < 0.0001
        
        vcf = '''#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT S1 S2 S3 S4 S5 S6 S7 S8
20\t14\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/0\t1/1\t1/1\t0/0\t0/0\t0/1\t1/1\t0/0
20\t14\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/0\t1/1\t1/1\t0/0\t0/0\t0/1\t1/1\t1/1'''

        vcf = StringIO(VCF_HEADER + vcf)
        snps = list(VCFReader(vcf).parse_snvs())
        call1 =  snps[0].record.samples
        call2 =  snps[1].record.samples
        counts = _count_biallelic_haplotypes(call1, call2)

        r_sqr = calculate_r_sqr(snps[0], snps[1])
        assert r_sqr - 1.0 < 0.0001

        # monomorphic
        vcf = '''#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT S1 S2 S3 S4 S5 S6 S7 S8
20\t2\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0
20\t3\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0'''

        vcf = StringIO(VCF_HEADER + vcf)
        snps = list(VCFReader(vcf).parse_snvs())
        call1 =  snps[0].record.samples
        call2 =  snps[1].record.samples
        assert _count_biallelic_haplotypes(call1, call2) is None
        r_sqr = calculate_r_sqr(snps[0], snps[1])
        assert r_sqr is None

        # Ab and aB
        vcf = '''#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT S1 S2 S3 S4 S5 S6 S7
20\t2\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/0\t0/0\t0/0\t1/1\t1/1\t1/1\t0/0\t
20\t3\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/0\t0/0\t0/0\t1/1\t1/1\t0/0\t1/1\t'''

        vcf = StringIO(VCF_HEADER + vcf)
        snps = list(VCFReader(vcf).parse_snvs())
        call1 =  snps[0].record.samples
        call2 =  snps[1].record.samples
        counts = _count_biallelic_haplotypes(call1, call2)
        assert counts.AB == 3
        assert counts.ab == 2
        assert counts.aB == 1
        assert counts.Ab == 1

        # different major allele names in snp1 (1, 2) and snp2 (2,3)
        vcf = '''#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT S1 S2 S3 S4 S5 S6 S7
20\t2\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/0\t0/0\t0/0\t1/1\t1/1\t1/1\t0/0\t
20\t3\t.\tG\tA\t29\tPASS\tNS=3\tGT\t3/3\t3/3\t3/3\t2/2\t2/2\t3/3\t3/3\t'''

        vcf = StringIO(VCF_HEADER + vcf)
        snps = list(VCFReader(vcf).parse_snvs())
        call1 =  snps[0].record.samples
        call2 =  snps[1].record.samples
        counts = _count_biallelic_haplotypes(call1, call2)
        assert counts.AB == 4
        assert counts.ab == 2
        assert counts.aB == 1
        assert counts.Ab == 0

        # missing data
        vcf = '''#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT S1 S2 S3 S4 S5 S6 S7
20\t2\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/0\t0/0\t0/0\t1/1\t1/1\t1/1\t./.\t
20\t3\t.\tG\tA\t29\tPASS\tNS=3\tGT\t3/3\t3/3\t3/3\t2/2\t2/2\t3/3\t3/3\t'''

        vcf = StringIO(VCF_HEADER + vcf)
        snps = list(VCFReader(vcf).parse_snvs())
        call1 =  snps[0].record.samples
        call2 =  snps[1].record.samples
        counts = _count_biallelic_haplotypes(call1, call2)
        assert counts.AB == 3
        assert counts.ab == 2
        assert counts.aB == 1
        assert counts.Ab == 0

    def test_r_example(self):
        # r examples
        assert _calculate_r_sqr(HaploCount(441, 435, 13, 111)) - 0.592 < 0.001

if __name__ == "__main__":
    # import sys;sys.argv = ['', 'ImputationTest.test_stars']
    # import sys;sys.argv = ['', 'VCFReadingPerChr']
    unittest.main()
