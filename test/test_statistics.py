import unittest
from os.path import join
from tempfile import NamedTemporaryFile
from subprocess import check_output

from vcf import Reader
from vcf_crumbs.statistics import (get_snpcaller_name, VARSCAN, GATK,
                                   calculate_maf_dp, FREEBAYES, VcfStats,
                                   HOM_REF, VCFcomparisons, _AlleleCounts2D,
                                   HOM_ALT, HET, HOM)

from vcf_crumbs.utils import TEST_DATA_DIR, BIN_DIR

VARSCAN_VCF_PATH = join(TEST_DATA_DIR, 'sample.vcf.gz')
REF_PATH = join(TEST_DATA_DIR, 'sample_ref.fasta')
GATK_VCF_PATH = join(TEST_DATA_DIR, 'gatk_sample.vcf.gz')
FREEBAYES_VCF_PATH = join(TEST_DATA_DIR, 'freebayes_sample.vcf.gz')
FREEBAYES_MULTI_VCF_PATH = join(TEST_DATA_DIR, 'freebayes_multisample.vcf.gz')


class TestVcfStats(unittest.TestCase):
    def  test_vcf_stats(self):
        vcf_stats = VcfStats(VARSCAN_VCF_PATH,
                             min_samples_for_heterozigosity=2)
        assert vcf_stats.mafs_dp().count == 153
        assert vcf_stats.mafs_dp()[63] == 5
        assert vcf_stats.mafs_dp('pepo').count == 29
        assert vcf_stats.gt_quals(HET)[21] == 2
        assert vcf_stats.gt_quals(HOM)[3] == 25
        assert vcf_stats.gt_quals(HET).count == 53
        assert (0.28 - vcf_stats.heterozigosity_for_sample('pepo')) < 0.01
        assert vcf_stats.het_by_snp[0] == 46


class AlleleCount2DTest(unittest.TestCase):
    def test_allele_count2d(self):
        allelecount = _AlleleCounts2D()
        allelecount.add(2, 3, (0, 0), 25)
        allelecount.add(2, 3, (1, 0), 25)
        allelecount.add(2, 3, (1, 1), 25)
        allelecount.add(2, 3, (0, 0), 25)
        allelecount.add(2, 3, (0, 1), 50)
        allelecount.add(2, 3, (2, 2), 25)
        allelecount.add(2, 3, (1, 0), 75)
        allelecount.add(2, 4, (1, 0), 75)

        assert allelecount.get_gt_count(2, 3, HOM_ALT) == 2
        assert allelecount.get_gt_count(2, 3, HOM_REF) == 2
        assert allelecount.get_gt_count(2, 3, HET) == 3
        assert allelecount.get_avg_gt_qual(2, 3, HOM_ALT) == 25
        assert allelecount.get_avg_gt_qual(2, 3, HOM_REF) == 25
        assert allelecount.get_avg_gt_qual(2, 3, HET) == 50

        allelecount.get_gt_depths_for_coverage(5)


class SnvStatTests(unittest.TestCase):

    def test_get_snpcaller(self):
        assert get_snpcaller_name(Reader(filename=VARSCAN_VCF_PATH)) == \
                                    VARSCAN
        assert get_snpcaller_name(Reader(filename=GATK_VCF_PATH)) == GATK

        assert get_snpcaller_name(Reader(filename=FREEBAYES_VCF_PATH)) == \
                                                                FREEBAYES

    def test_calc_maf(self):
        #varscan
        reader = Reader(filename=VARSCAN_VCF_PATH)
        snp = reader.next()
        maf = calculate_maf_dp(snp, vcf_variant=VARSCAN)
        assert 0.52 < maf['all'] < 0.53
        assert maf['upv196'] == 1

        #gatk
        reader = Reader(filename=GATK_VCF_PATH)
        snp = reader.next()
        maf = calculate_maf_dp(snp, vcf_variant=GATK)
        assert 0.7 < maf['all'] < 0.72
        assert 0.7 < maf['hib_amarillo'] < 0.72

        #freebayes
        reader = Reader(filename=FREEBAYES_VCF_PATH)
        snp = reader.next()
        maf = calculate_maf_dp(snp, vcf_variant=FREEBAYES)
        assert maf == {'all': 1.0, 'pep': 1.0}


class VCFcomparisonsTest(unittest.TestCase):
    def test_calculate_statistics(self):
        #with freebayes
        reader = Reader(filename=FREEBAYES_VCF_PATH)
        vcf_to_compare = VCFcomparisons(FREEBAYES_VCF_PATH)
        stats = vcf_to_compare.calculate_statistics(reader)
        assert stats['common'] == 944
        assert stats['uncalled'] == 0
        assert stats['different'] == 0
        assert stats['common_snps_prc'] == 100

        #with varscan
        reader = Reader(filename=VARSCAN_VCF_PATH)
        vcf_to_compare = VCFcomparisons(VARSCAN_VCF_PATH, samples=['mu16'])
        stats = vcf_to_compare.calculate_statistics(reader, samples=['mu16'])
        assert stats['common'] == 107
        assert stats['uncalled'] == 69
        assert stats['different'] == 0
        assert stats['common_snps_prc'] == 100

    def test_compare_vcfs_samples(self):
        binary = join(BIN_DIR, 'compare_vcfs_samples')
        assert 'usage' in check_output([binary, '-h'])
        samples_fhand = NamedTemporaryFile()
        samples_fhand.write('mu16\n')
        samples_fhand.flush()

        cmd = [binary, VARSCAN_VCF_PATH, '-s', samples_fhand.name,
               '-r', samples_fhand.name, '-v', VARSCAN_VCF_PATH]
        stats = check_output(cmd)
        result = 'common_snps_prc : 100.0\ndifferent : 0\ncommon : 107\n'
        result += 'uncalled : 69\n'
        assert stats == result


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'AlleleDepthsTests']
    unittest.main()
