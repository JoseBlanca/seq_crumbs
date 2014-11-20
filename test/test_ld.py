import unittest
from os.path import join
from StringIO import StringIO
from tempfile import NamedTemporaryFile

from vcf_crumbs.snv import VCFReader
from vcf_crumbs.ld import (_count_biallelic_haplotypes, calculate_r_sqr,
                           HaploCount, _calculate_r_sqr, _fisher_exact,
                           calculate_ld_stats, filter_snvs_by_ld, fisher_exact,
                           _LDStatsCache, _calc_recomb_rate,
                           calc_recomb_rates_along_chroms)
from vcf_crumbs.utils.file_utils import BIN_DIR, TEST_DATA_DIR
from subprocess import check_call

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
##contig=<ID=20,length=62435964,assembly=B36>
##phasing=partial
##INFO=<ID=IV0,Number=1,Type=String,Description="True">
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

FREEBAYES_VCF_PATH = join(TEST_DATA_DIR, 'freebayes_multisample.vcf.gz')

class LDTests(unittest.TestCase):
    def test_empy_snv(self):
        vcf = '''#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT 1 2 3 4 5 6
20\t14\t.\tG\tA\t29\tPASS\tNS=3\tGT\t./.\t./.\t./.\t./.\t./.\t./.
20\t14\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/0\t1/1\t1/1\t0/0\t0/0\t0/1'''

        vcf = StringIO(VCF_HEADER + vcf)
        snps = list(VCFReader(vcf).parse_snvs())
        call1 = snps[0].record.samples
        call2 = snps[1].record.samples
        counts = _count_biallelic_haplotypes(call1, call2)
        assert counts is None

    def test_count_homo_haplotypes(self):
        vcf = '''#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT 1 2 3 4 5 6
20\t14\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/0\t1/1\t1/1\t0/0\t0/0\t0/1
20\t14\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/0\t1/1\t1/1\t0/0\t0/0\t0/1'''

        vcf = StringIO(VCF_HEADER + vcf)
        snps = list(VCFReader(vcf).parse_snvs())
        call1 = snps[0].record.samples
        call2 = snps[1].record.samples
        counts = _count_biallelic_haplotypes(call1, call2)
        assert counts.AB == 3
        assert counts.ab == 2

        vcf = '''#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT 1 2 3 4 5 6 7 8
20\t14\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/0\t1/1\t1/1\t0/0\t0/0\t0/1\t1/1\t0/0
20\t14\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/0\t1/1\t1/1\t0/0\t0/0\t0/1\t1/1\t0/0'''

        vcf = StringIO(VCF_HEADER + vcf)
        snps = list(VCFReader(vcf).parse_snvs())
        call1 = snps[0].record.samples
        call2 = snps[1].record.samples
        counts = _count_biallelic_haplotypes(call1, call2)
        assert counts.AB == 4
        assert counts.ab == 3

        vcf = '''#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT 1 2 3 4 5 6 7 8
20\t14\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/0\t1/1\t1/1\t0/0\t0/0\t0/1\t1/1\t0/0
20\t14\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/0\t2/2\t1/1\t0/0\t0/0\t0/1\t1/1\t0/0'''

        vcf = StringIO(VCF_HEADER + vcf)
        snps = list(VCFReader(vcf).parse_snvs())
        call1 = snps[0].record.samples
        call2 = snps[1].record.samples
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
        call1 = snps[0].record.samples
        call2 = snps[1].record.samples
        counts = _count_biallelic_haplotypes(call1, call2)

        r_sqr = calculate_r_sqr(snps[0], snps[1])
        assert r_sqr - 1.0 < 0.0001

        vcf = '''#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT 1 2 3 4 5 6 7 8
20\t14\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/0\t1/1\t1/1\t0/0\t0/0\t0/1\t1/1\t0/0
20\t14\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/0\t1/1\t1/1\t0/0\t0/0\t0/1\t1/1\t1/1'''

        vcf = StringIO(VCF_HEADER + vcf)
        snps = list(VCFReader(vcf).parse_snvs())
        call1 = snps[0].record.samples
        call2 = snps[1].record.samples
        counts = _count_biallelic_haplotypes(call1, call2)

        r_sqr = calculate_r_sqr(snps[0], snps[1])
        assert r_sqr - 1.0 < 0.0001

        # monomorphic
        vcf = '''#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT 1 2 3 4 5 6 7 8
20\t2\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0
20\t3\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0'''

        vcf = StringIO(VCF_HEADER + vcf)
        snps = list(VCFReader(vcf).parse_snvs())
        call1 = snps[0].record.samples
        call2 = snps[1].record.samples
        assert _count_biallelic_haplotypes(call1, call2) is None
        r_sqr = calculate_r_sqr(snps[0], snps[1])
        assert r_sqr is None
        assert fisher_exact(snps[0], snps[1]) is None

        # Ab and aB
        vcf = '''#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT 1 2 3 4 5 6 7
20\t2\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/0\t0/0\t0/0\t1/1\t1/1\t1/1\t0/0\t
20\t3\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/0\t0/0\t0/0\t1/1\t1/1\t0/0\t1/1\t'''

        vcf = StringIO(VCF_HEADER + vcf)
        snps = list(VCFReader(vcf).parse_snvs())
        call1 = snps[0].record.samples
        call2 = snps[1].record.samples
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
        call1 = snps[0].record.samples
        call2 = snps[1].record.samples
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
        call1 = snps[0].record.samples
        call2 = snps[1].record.samples
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
        snvs = filter_snvs_by_ld(snps, p_val=0.03)
        assert not list(snvs)

        vcf = '''#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT 1 2 3 4 5 6 7 8
20\t2\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/0\t0/0\t0/0\t0/0\t1/1\t1/1\t1/1\t1/1\t
20\t3\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/0\t0/0\t0/0\t0/0\t1/1\t1/1\t1/1\t1/1\t'''
        vcf = StringIO(VCF_HEADER + vcf)
        snps = VCFReader(vcf).parse_snvs()
        snvs = filter_snvs_by_ld(snps, p_val=0.03)
        assert not list(snvs)

        vcf = '''#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT 1 2 3 4 5 6 7 8
20\t2\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/0\t0/0\t0/0\t0/0\t1/1\t1/1\t1/1\t1/1\t
20\t703\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/0\t0/0\t0/0\t0/0\t1/1\t1/1\t1/1\t1/1\t
'''
        vcf = StringIO(VCF_HEADER + vcf)
        snps = VCFReader(vcf).parse_snvs()
        snvs = filter_snvs_by_ld(snps, p_val=0.03, bonferroni=False)
        assert len(list(snvs)) == 2

        vcf = '''#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT 1 2 3 4 5 6 7 8
20\t2\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/0\t0/0\t0/0\t0/0\t1/1\t1/1\t1/1\t1/1\t
20\t703\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/0\t0/0\t0/0\t0/0\t1/1\t1/1\t1/1\t1/1\t
'''
        vcf = StringIO(VCF_HEADER + vcf)
        snps = VCFReader(vcf).parse_snvs()
        snvs = filter_snvs_by_ld(snps, p_val=0.03)
        assert not list(snvs)

        vcf = '''#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT 1 2 3 4 5 6 7 8
20\t2\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/0\t0/0\t0/0\t0/0\t1/1\t1/1\t1/1\t1/1\t
20\t703\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/0\t0/0\t0/0\t0/0\t1/1\t1/1\t1/1\t1/1\t
20\t2003\t.\tG\tA\t29\tPASS\tNS=3\tGT\t1/1\t0/0\t1/1\t0/0\t1/1\t0/0\t1/1\t0/0\t
'''
        vcf = StringIO(VCF_HEADER + vcf)
        snps = VCFReader(vcf).parse_snvs()
        snvs = filter_snvs_by_ld(snps, p_val=0.03, bonferroni=False)
        assert [s.pos for s in snvs] == [1, 702]

        vcf = '''#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT 1 2 3 4 5 6 7 8
20\t2\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0
20\t3\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0'''

        vcf = StringIO(VCF_HEADER + vcf)
        snps = VCFReader(vcf).parse_snvs()
        snvs = filter_snvs_by_ld(snps, p_val=0.03, bonferroni=False)
        assert not list(snvs)

        vcf = '''#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT 1 2 3 4 5 6 7 8
20\t2\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0
20\t3\t.\tG\tA\t29\tPASS\tNS=3\tGT\t1/1\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0'''

        vcf = StringIO(VCF_HEADER + vcf)
        snps = VCFReader(vcf).parse_snvs()
        snvs = filter_snvs_by_ld(snps, p_val=0.03, bonferroni=False)
        assert not list(snvs)

    def test_check_backwards(self):
        vcf = '''#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT 1 2 3 4 5 6 7 8
20\t2\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/0\t0/0\t0/0\t0/0\t1/1\t1/1\t1/1\t1/1\t
20\t703\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/0\t0/0\t0/0\t0/0\t1/1\t1/1\t1/1\t1/1\t
20\t2003\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/0\t0/0\t0/0\t0/0\t1/1\t1/1\t1/1\t1/1\t
20\t2403\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/0\t0/0\t0/0\t0/0\t1/1\t1/1\t1/1\t1/1\t
'''
        vcf = StringIO(VCF_HEADER + vcf)
        snps = VCFReader(vcf).parse_snvs()
        snvs = filter_snvs_by_ld(snps, p_val=0.03, bonferroni=False)
        assert [s.pos for s in snvs] == [1, 702, 2002, 2402]

    def test_nobreak_generator(self):
        vcf = '''#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT 1 2 3 4 5 6 7 8
20\t2\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/0\t0/0\t0/0\t0/0\t1/1\t1/1\t1/1\t1/1\t
20\t703\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/0\t0/0\t0/0\t0/0\t1/1\t1/1\t1/1\t1/1\t
20\t2003\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/0\t0/0\t0/0\t0/0\t1/1\t1/1\t1/1\t1/1\t
20\t3003\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/0\t0/0\t0/0\t0/0\t1/1\t1/1\t1/1\t1/1\t
20\t3403\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/0\t0/0\t0/0\t0/0\t1/1\t1/1\t1/1\t1/1\t
'''
        vcf = StringIO(VCF_HEADER + vcf)
        snps = VCFReader(vcf).parse_snvs()
        snvs = filter_snvs_by_ld(snps, p_val=0.03, bonferroni=False,
                                 snv_win=3)
        assert [s.pos for s in snvs] == [1, 702, 2002, 3002]

    def test_cache(self):
        vcf = '''#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT 1 2 3 4 5 6 7 8
20\t2\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/0\t0/0\t0/0\t0/0\t1/1\t1/1\t1/1\t1/1\t
20\t703\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/0\t0/0\t0/0\t0/0\t1/1\t1/1\t1/1\t1/1\t
20\t2003\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/0\t0/0\t0/0\t0/0\t1/1\t1/1\t1/1\t1/1\t
20\t3003\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/0\t0/0\t0/0\t0/0\t1/1\t1/1\t1/1\t1/1\t
20\t3403\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/0\t0/0\t0/0\t0/0\t1/1\t1/1\t1/1\t1/1\t
'''
        vcf = StringIO(VCF_HEADER + vcf)
        snps = VCFReader(vcf).parse_snvs()
        snvs = filter_snvs_by_ld(snps, p_val=0.001, bonferroni=False,
                                 snv_win=3)
        assert not list(snvs)

    def test_binary(self):
        vcf = '''#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT 1 2 3 4 5 6 7 8
20\t2\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/0\t0/0\t0/0\t0/0\t1/1\t1/1\t1/1\t1/1\t
20\t703\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/0\t0/0\t0/0\t0/0\t1/1\t1/1\t1/1\t1/1\t
20\t2003\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/0\t0/0\t0/0\t0/0\t1/1\t1/1\t1/1\t1/1\t
'''
        fhand = NamedTemporaryFile()
        fhand.write(VCF_HEADER + vcf)
        fhand.flush()
        stderr = NamedTemporaryFile()
        out_fhand = NamedTemporaryFile()

        binary = join(BIN_DIR, 'filter_vcf_by_ld')
        cmd = [binary, '-o', out_fhand.name, fhand.name,
               '--no_bonferroni_correction', '--p_val', '0.03']
        check_call(cmd, stderr=stderr)
        assert len(list(VCFReader(open(out_fhand.name)).parse_snvs())) == 3

        log_fhand = NamedTemporaryFile()
        binary = join(BIN_DIR, 'filter_vcf_by_ld')
        cmd = [binary, '-o', out_fhand.name, fhand.name,
               '--no_bonferroni_correction', '--p_val', '0.03',
               '-l', log_fhand.name]
        check_call(cmd, stderr=stderr)
        assert 'filtered' in open(log_fhand.name).read()


class RecombRateTest(unittest.TestCase):
    def test_recomb_rate(self):
        # samples
        vcf = '''#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT 1 2 3 4 5 6 7 8
20\t2\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/0\t0/0\t0/0\t0/0\t1/1\t1/1\t1/1\t1/1\t
20\t3\t.\tG\tA\t29\tPASS\tNS=3\tGT\t0/0\t0/0\t0/0\t0/0\t1/1\t1/1\t1/1\t1/1\t
20\t4\t.\tG\tA\t29\tPASS\tNS=3\tGT\t1/1\t0/0\t1/1\t0/0\t0/0\t1/1\t0/0\t1/1\t
20\t6\t.\tG\tA\t29\tPASS\tNS=3\tGT\t./.\t./.\t./.\t./.\t./.\t0/1\t0/1\t0/1\t
21\t4\t.\tG\tA\t29\tPASS\tNS=3\tGT\t1/1\t0/0\t1/1\t0/0\t0/0\t1/1\t0/0\t1/1\t
'''
        vcf = StringIO(VCF_HEADER + vcf)
        snps = list(VCFReader(vcf).parse_snvs())

        recomb = _calc_recomb_rate(snps[0].record.samples,
                                   snps[1].record.samples,
                                   'ril_self')
        self.assertAlmostEqual(recomb, 0.0, 3)
        recomb = _calc_recomb_rate(snps[0].record.samples,
                                   snps[2].record.samples,
                                   'ril_self')
        self.assertAlmostEqual(recomb, 0.375, 3)
        recomb = _calc_recomb_rate(snps[0].record.samples,
                                   snps[2].record.samples,
                                   'test_cross')
        self.assertAlmostEqual(recomb, 0.5, 3)
        recomb = _calc_recomb_rate(snps[0].record.samples,
                                   snps[3].record.samples,
                                   'test_cross')
        assert recomb is None

        vcf = '''#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t1_14_1_gbs\t1_17_1_gbs\t1_18_4_gbs\t1_19_4_gbs\t1_26_1_gbs\t1_27_1_gbs1_2_2_gbs\t1_35_13_gbs\t1_3_2_gbs\t1_50_1_gbs\t1_59_1_gbs\t1_63_4_gbs\t1_6_2_gbs\t1_70_1_gbs\t1_74_1_gbs\t1_79_1_gbs\t1_7_2_gbs\t1_81_10_gbs\t1_86_1_gbs\t1_8_2_gbs\t1_91_2_gbs\t1_94_4_gbs\t2_107_1_gbs\t2_10_2_gbs\t2_116_1_gbs\t2_11_1_gbs\t2_125_2_gbs\t2_13_1_gbs\t2_16_3_gbs\t2_21_1_gbs\t2_22A_1_gbs\t2_24_2_gbs\t2_28_2_gbs\t2_31_2_gbs\t2_33_1_gbs\t2_39_3_gbs\t2_43_1_gbs2_5_1_gbs\t2_64_7_gbs\t2_67_2_gbs\t2_6_4_gbs\t2_84_2_gbs\t2_8_3_gbs\t2_95_2_gbs\t4_100B_4_gbs\t4_108_10_gbs\t4_110_11_gbs\t4_111_6_gbs\t4_115B_2_gbs\t4_11B_3_gbs\t4_123B_2_gbs\t4_127_6_gbs\t4_131_1_gbs\t4_136B_3_gbs\t4_136_10_T1_gbs\t4_138B_2_gbs\t4_26_11_gbs\t4_28_4_gbs\t4_33_2_gbs\t4_35_1_gbs\t4_38_2_gbs\t4_39_2_gbs\t4_41B_2_gbs\t4_42_11_gbs\t4_45_2_gbs\t4_53_2_gbs\t4_5_5_gbs\t4_62_4_gbs\t4_64B_1_gbs\t4_65_5_gbs\t4_66_2_gbs\t4_71_2_gbs\t4_72_1_gbs\t4_77_1_gbs\t4_7B_1_gbs\t4_7_2_gbs\t4_81B_2_gbs\t4_82B_4_gbs\t4_85_1_gbs\t4_95_1_gbs\t4_9_1_gbs\t5_14B_1_gbs\t5_15B_1_gbs\t5_18_1_gbs\t5_22_2_gbs\t5_24_2_gbs\t5_25_2_gbs\t5_32_3_gbs\t5_33B_4_gbs\t5_34B_2_gbs\t5_3_1_gbs\t5_40B_2_gbs\t5_49B_2_T1_gbs\t5_57_1_gbs\t5_58_1_gbs\t5_66_1_gbs\t5_80B_2_gbs\tMU_16_5_gbs\tV_196_2_gbs\t1\t2
s7\t4039693\tS7_4039693\tT\tG\t.\tPASS\tIV0=F\tGT\t0/0\t0/0\t0/0\t1/1\t0/0\t1/1\t1/1\t1/1\t1/1\t0/0\t0/0\t0/0\t1/1\t0/0\t0/0\t1/1\t0/0\t1/1\t0/0\t0/0\t0/0\t1/1\t0/0\t0/0\t0/0\t0/0\t1/1\t1/1\t0/0\t0/0\t0/0\t0/0\t1/1\t0/0\t1/1\t0/0\t0/0\t1/1\t1/1\t0/0\t1/1\t1/1\t1/1\t0/0\t1/1\t1/1\t1/1\t0/0\t1/1\t1/1\t0/0\t0/0\t0/0\t0/0\t0/0\t1/1\t0/0\t0/0\t./.\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t1/1\t0/0\t0/0\t1/1\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t1/1\t0/0\t1/1\t0/0\t0/0\t0/0\t1/1\t1/1\t1/1\t1/1\t1/1\t0/0\t1/1\t0/0\t0/0\t0/0\t0/0\t1/1
s7\t4028261\tS7_4028261\tC\tT\t.\tPASS\tIV0=F\tGT\t1/1\t1/1\t./.\t0/0\t1/1\t0/0\t./.\t0/0\t0/0\t1/1\t1/1\t1/1\t0/0\t1/1\t1/1\t0/0\t1/1\t0/0\t1/1\t1/1\t1/1\t0/0\t1/1\t1/1\t1/1\t1/1\t0/0\t0/0\t1/1\t1/1\t1/1\t1/1\t0/0\t1/1\t0/0\t1/1\t1/1\t0/0\t0/0\t1/1\t0/0\t0/0\t0/0\t1/1\t0/0\t0/0\t0/0\t0/0\t0/0\t./.\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t0/1\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t1/1\t0/0
'''
        vcf = StringIO(VCF_HEADER + vcf)
        reader = VCFReader(vcf)
        snps = list(reader.parse_snvs())

        recomb = _calc_recomb_rate(snps[0].record.samples,
                                   snps[1].record.samples,
                                   'ril_self')
        self.assertAlmostEqual(recomb, 0.8187, 3)

    def test_recomb_rate_along_chrom(self):
        vcf = FREEBAYES_VCF_PATH

        res = calc_recomb_rates_along_chroms(vcf, pop_type='test_cross')
        assert not list(res)

        samples = ['sample05_gbs', 'sample06_gbs', 'sample07_gbs',
                   'sample08_gbs', 'sample09_gbs']
        res = calc_recomb_rates_along_chroms(vcf, pop_type='test_cross',
                                             samples=samples)
        assert not list(res)

class _LDStatsCacheTest(unittest.TestCase):
    def test_ld_stats(self):
        stats = _LDStatsCache()
        stats.set_stat(2, 1, 'result')
        assert stats.get_stat(1, 2) == 'result'
        stats.del_lower_than(10)
        try:
            stats.get_stat(1, 2)
            self.fail('KeyError expected')
        except KeyError:
            pass


if __name__ == "__main__":
    # import sys; sys.argv = ['', 'RecombRateTest.test_recomb_rate']
    unittest.main()
