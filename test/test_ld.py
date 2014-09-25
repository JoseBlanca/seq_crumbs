'''
Created on 23/09/2014

@author: jose
'''
import unittest
from StringIO import StringIO

from vcf_crumbs.snv import VCFReader
from vcf_crumbs.ld import (_count_biallelic_haplotypes, calculate_r_sqr,
                           HaploCount, _calculate_r_sqr, _fisher_exact,
                           calculate_ld_stats, filter_snvs_by_ld)

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


class LDTests(unittest.TestCase):
    def test_count_homo_haplotypes(self):
        vcf = '''#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT 1 2 3 4 5 6
20\t14\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/0\t1/1\t1/1\t0/0\t0/0\t0/1
20\t14\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/0\t1/1\t1/1\t0/0\t0/0\t0/1'''

        vcf = StringIO(VCF_HEADER + vcf)
        snps = list(VCFReader(vcf).parse_snvs())
        call1 =  snps[0].record.samples
        call2 =  snps[1].record.samples
        counts = _count_biallelic_haplotypes(call1, call2)
        assert counts.AB == 3
        assert counts.ab == 2

        vcf = '''#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT 1 2 3 4 5 6 7 8
20\t14\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/0\t1/1\t1/1\t0/0\t0/0\t0/1\t1/1\t0/0
20\t14\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/0\t1/1\t1/1\t0/0\t0/0\t0/1\t1/1\t0/0'''

        vcf = StringIO(VCF_HEADER + vcf)
        snps = list(VCFReader(vcf).parse_snvs())
        call1 =  snps[0].record.samples
        call2 =  snps[1].record.samples
        counts = _count_biallelic_haplotypes(call1, call2)
        assert counts.AB == 4
        assert counts.ab == 3

        vcf = '''#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT 1 2 3 4 5 6 7 8
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

        vcf = '''#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT 1 2 3 4 5 6 7 8
20\t14\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/0\t1/1\t1/1\t0/0\t0/0\t0/1\t1/1\t0/0
20\t14\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/0\t0/0\t1/1\t1/1\t0/0\t1/1\t0/0\t1/1'''

        vcf = StringIO(VCF_HEADER + vcf)
        snps = list(VCFReader(vcf).parse_snvs())
        call1 =  snps[0].record.samples
        call2 =  snps[1].record.samples
        counts = _count_biallelic_haplotypes(call1, call2)

        r_sqr = calculate_r_sqr(snps[0], snps[1])
        assert r_sqr - 1.0 < 0.0001

        vcf = '''#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT 1 2 3 4 5 6 7 8
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
        vcf = '''#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT 1 2 3 4 5 6 7 8
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
        vcf = '''#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT 1 2 3 4 5 6 7
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
        vcf = '''#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT 1 2 3 4 5 6 7
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
        vcf = '''#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT 1 2 3 4 5 6 7
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
        self.assertAlmostEqual(_calculate_r_sqr(HaploCount(10, 10, 10, 10)),
                               0)
        self.assertAlmostEqual(_fisher_exact(HaploCount(10, 10, 10, 10)), 1)
        self.assertAlmostEqual(_calculate_r_sqr(HaploCount(10, 0, 0, 10)), 1)

        self.assertAlmostEqual(_calculate_r_sqr(HaploCount(441, 13, 111, 435)),
                               0.591332576)
        self.assertAlmostEqual(_fisher_exact(HaploCount(6, 6, 2, 6)),
                               0.3728506787)
        self.assertAlmostEqual(_fisher_exact(HaploCount(1, 0, 5, 7)),
                               0.4615385)
        self.assertAlmostEqual(_fisher_exact(HaploCount(5, 23, 1, 20)),
                               0.219157345)

        vcf = '''#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT 1 2 3 4 5 6 7
20\t2\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/0\t0/0\t0/0\t1/1\t1/1\t1/1\t./.\t
20\t3\t.\tG\tA\t29\tPASS\tNS=3\tGT\t3/3\t3/3\t3/3\t2/2\t2/2\t3/3\t3/3\t'''

        vcf = StringIO(VCF_HEADER + vcf)
        snps = list(VCFReader(vcf).parse_snvs())
        ld_stats = calculate_ld_stats(snps[0], snps[1])
        self.assertAlmostEqual(ld_stats.fisher, 0.39999999999)
        self.assertAlmostEqual(ld_stats.r_sqr, 0.49999999)


class FilterTest(unittest.TestCase):
    def test_filter(self):
        vcf = '''#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT 1 2 3 4 5 6 7 8
20\t2\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/0\t0/0\t0/0\t0/0\t1/1\t1/1\t1/1\t1/1\t
21\t3\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/0\t0/0\t0/0\t0/0\t1/1\t1/1\t1/1\t1/1\t'''
        vcf = StringIO(VCF_HEADER + vcf)
        snps = VCFReader(vcf).parse_snvs()
        snvs = filter_snvs_by_ld(snps, p_bonferroni=0.03)
        assert not list(snvs)

        vcf = '''#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT 1 2 3 4 5 6 7 8
20\t2\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/0\t0/0\t0/0\t0/0\t1/1\t1/1\t1/1\t1/1\t
20\t3\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/0\t0/0\t0/0\t0/0\t1/1\t1/1\t1/1\t1/1\t'''
        vcf = StringIO(VCF_HEADER + vcf)
        snps = VCFReader(vcf).parse_snvs()
        snvs = filter_snvs_by_ld(snps, p_bonferroni=0.03)
        assert not list(snvs)

        vcf = '''#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT 1 2 3 4 5 6 7 8
20\t2\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/0\t0/0\t0/0\t0/0\t1/1\t1/1\t1/1\t1/1\t
20\t703\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/0\t0/0\t0/0\t0/0\t1/1\t1/1\t1/1\t1/1\t'''
        vcf = StringIO(VCF_HEADER + vcf)
        snps = VCFReader(vcf).parse_snvs()
        snvs = filter_snvs_by_ld(snps, p_bonferroni=0.03)
        #print [s.pos for s in snvs]
        assert len(list(snvs)) == 2

        vcf = '''#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT 1 2 3 4 5 6 7 8
20\t2\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/0\t0/0\t0/0\t0/0\t1/1\t1/1\t1/1\t1/1\t
20\t703\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/0\t0/0\t0/0\t0/0\t1/1\t1/1\t1/1\t1/1\t
20\t2003\t.\tG\tA\t29\tPASS\tNS=3\tGT\t1/1\t0/0\t1/1\t0/0\t1/1\t0/0\t1/1\t0/0\t'''
        vcf = StringIO(VCF_HEADER + vcf)
        snps = VCFReader(vcf).parse_snvs()
        snvs = filter_snvs_by_ld(snps, p_bonferroni=0.03)
        assert [s.pos for s in snvs] == [1, 702]

    def test_check_backwards(self):
        vcf = '''#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT 1 2 3 4 5 6 7 8
20\t2\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/0\t0/0\t0/0\t0/0\t1/1\t1/1\t1/1\t1/1\t
20\t703\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/0\t0/0\t0/0\t0/0\t1/1\t1/1\t1/1\t1/1\t
20\t2003\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/0\t0/0\t0/0\t0/0\t1/1\t1/1\t1/1\t1/1\t'''
        vcf = StringIO(VCF_HEADER + vcf)
        snps = VCFReader(vcf).parse_snvs()
        snvs = filter_snvs_by_ld(snps, p_bonferroni=0.03)
        assert [s.pos for s in snvs] == [1, 702, 2002]

if __name__ == "__main__":
    # import sys;sys.argv = ['', 'ImputationTest.test_stars']
    # import sys;sys.argv = ['', 'VCFReadingPerChr']
    unittest.main()
