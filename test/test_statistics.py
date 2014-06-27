import unittest
from os.path import join
import sys
from tempfile import NamedTemporaryFile

from vcf import Reader
from vcf_crumbs.statistics import (calc_density_per_chrom,
                                   get_snpcaller_name, VARSCAN, GATK,
                                   calculate_maf, FREEBAYES, VcfStats,
                                   calc_n_bases_in_chrom_with_snp, HOM_REF,
                                   VCFcomparisons, _AlleleCounts2D, HOM_ALT,
                                   HET, HOM, VcfStats_old)

from vcf_crumbs.utils import TEST_DATA_DIR, BIN_DIR
from subprocess import check_call, CalledProcessError, check_output

from crumbs.utils.file_utils import TemporaryDir
from crumbs.statistics import IntCounter
from crumbs.plot import draw_scatter


VARSCAN_VCF_PATH = join(TEST_DATA_DIR, 'sample.vcf.gz')
REF_PATH = join(TEST_DATA_DIR, 'sample_ref.fasta')
GATK_VCF_PATH = join(TEST_DATA_DIR, 'gatk_sample.vcf.gz')
FREEBAYES_VCF_PATH = join(TEST_DATA_DIR, 'freebayes_sample.vcf.gz')
FREEBAYES_MULTI_VCF_PATH = join(TEST_DATA_DIR, 'freebayes_multisample.vcf.gz')


class TestVcfStats(unittest.TestCase):
    def test_vcf_stats_old(self):
        vcf_stats = VcfStats_old(VARSCAN_VCF_PATH, gq_threshold=0)
        assert vcf_stats.het_by_sample == \
                        {'upv196': IntCounter({'num_gt': 90, 'num_het': 19}),
                         'pepo': IntCounter({'num_gt': 29, 'num_het': 8}),
                         'mu16': IntCounter({'num_gt': 107, 'num_het': 26})}
        assert vcf_stats.snv_quals.count == 0
        vcf_stats = VcfStats_old(VARSCAN_VCF_PATH, gq_threshold=25)
        assert vcf_stats.het_by_sample == \
                        {'upv196': IntCounter({'num_gt': 90, 'num_het': 7}),
                         'pepo': IntCounter({'num_gt': 29, 'num_het': 4}),
                         'mu16': IntCounter({'num_gt': 107, 'num_het': 18})}
        assert len(vcf_stats.samples) == 3
        assert len(vcf_stats.mafs) == 4
        assert vcf_stats.mafs['all'][63] == 11
        assert vcf_stats.mafs['pepo'].count == 29
        assert vcf_stats.snps_per_chromosome['CUUC00423_TC01'] == 26
        assert len(vcf_stats.call_data[HOM_REF]['x']) == 95
        assert vcf_stats.genotype_qualities.count == 226
        assert vcf_stats.genotype_qualities[3] == 25
        assert vcf_stats.variable_gt_per_snp[100] == 50

    def  test_vcf_stats(self):
        vcf_stats = VcfStats(VARSCAN_VCF_PATH,
                             min_samples_for_heterozigosity=2)
        assert vcf_stats.mafs().count == 153
        assert vcf_stats.mafs()[63] == 11
        assert vcf_stats.mafs('pepo').count == 29
        assert vcf_stats.gt_quals(HET)[21] == 2
        assert vcf_stats.gt_quals(HOM)[3] == 25
        assert vcf_stats.gt_quals(HET).count == 53
        assert (0.28 - vcf_stats.heterozigotes_by_sample('pepo')) < 0.01
        assert vcf_stats.het_by_snp[0] == 46


class AlleleCount2DTest(unittest.TestCase):
    def test_allele_count2d(self):
        allelecount = _AlleleCounts2D()
        allelecount.add(2, 3, '0/0', 25)
        allelecount.add(2, 3, '1/0', 25)
        allelecount.add(2, 3, '1/1', 25)
        allelecount.add(2, 3, '0/0', 25)
        allelecount.add(2, 3, '0/1', 50)
        allelecount.add(2, 3, '2/2', 25)
        allelecount.add(2, 3, '1/0', 75)
        allelecount.add(2, 4, '1/0', 75)

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

    def test_calc_densities(self):
        vcf_stats = VcfStats(VARSCAN_VCF_PATH)
        densities = calc_density_per_chrom(vcf_stats.snps_per_chromosome,
                                           open(REF_PATH))
        assert densities['CUUC00355_TC01'] == 3.74

    def test_calc_n_bases_in_chrom_with_snp(self):
        vcf_stats = VcfStats(VARSCAN_VCF_PATH)
        counts = vcf_stats.snps_per_chromosome
        n_bases = calc_n_bases_in_chrom_with_snp(counts, open(REF_PATH))
        assert n_bases == 16393

    def test_calc_maf(self):
        #varscan
        reader = Reader(filename=VARSCAN_VCF_PATH)
        snp = reader.next()
        maf = calculate_maf(snp, vcf_variant=VARSCAN)
        assert 0.52 < maf['all'] < 0.53
        assert maf['upv196'] == 1

        #gatk
        reader = Reader(filename=GATK_VCF_PATH)
        snp = reader.next()
        maf = calculate_maf(snp, vcf_variant=GATK)
        assert 0.7 < maf['all'] < 0.72
        assert 0.7 < maf['hib_amarillo'] < 0.72

        #freebayes
        reader = Reader(filename=FREEBAYES_VCF_PATH)
        snp = reader.next()
        maf = calculate_maf(snp, vcf_variant=FREEBAYES)
        assert maf == {'all': 1.0, 'pep': 1.0}

    def test_scatter_calldata(self):
        vcf_stats = VcfStats(FREEBAYES_VCF_PATH)
        fhand = NamedTemporaryFile(suffix='.png')
        draw_scatter(vcf_stats.call_data.values(), fhand, xlim=(0, 100),
                     ylim=(0, 100))
        #raw_input(fhand.name)

    def test_counts_distribution_in_genotype(self):
        vcf_stats = VcfStats(VARSCAN_VCF_PATH)
        results_dp11 = {'0/0': IntCounter({7: 8, 8: 3, 11: 3, 6: 2, 10: 1}),
                        '0/1': IntCounter({4: 4, 5: 4, 3: 3}),
                        '1/1': IntCounter({2: 2})}
        assert vcf_stats.counts_distribution_in_gt[11] == results_dp11

        vcf_stats = VcfStats(FREEBAYES_VCF_PATH)
        results_dp11 = {
        '1/1/1/1/1/1/1/1/1/1/1/1/1/1/1/1/1/1/1/1': IntCounter({0: 7, 1: 2,
                                                               2: 2}),
         '0/0/0/0/0/0/0/0/0/0/0/0/0/0/0/1/1/1/1/1': IntCounter({5: 3}),
         '0/0/0/0/0/0/0/0/0/0/0/0/0/0/0/0/0/0/1/1': IntCounter({8: 10}),
         '0/0/0/0/0/0/0/0/0/0/0/0/0/0/0/0/0/1/1/1': IntCounter({7: 5}),
         '0/0/0/0/0/0/0/0/0/0/0/0/0/0/0/0/1/1/1/1': IntCounter({6: 2}),
         '0/0/0/0/0/0/0/0/0/0/0/0/0/0/1/1/1/1/1/1': IntCounter({4: 4}),
         '0/0/0/0/0/0/0/0/0/0/1/1/1/1/1/1/1/1/1/1': IntCounter({3: 2, 0: 1})}
#                         '0/1': IntCounter({8: 10, 7: 5, 4: 4, 5: 3, 3: 2, 6:
#                                                                    2, 0: 1}),
#                         '1/1': IntCounter({0: 7, 1: 2, 2: 2})}
        assert vcf_stats.counts_distribution_in_gt[11] == results_dp11
        #raw_input(fhand.name)


class StatBinTests(unittest.TestCase):

    def test_draw_snv_stats_bin(self):
        binary = join(BIN_DIR, 'draw_snv_stats')
        tempdir = TemporaryDir()
        cmd = [binary, '-r', REF_PATH, '-o', tempdir.name, VARSCAN_VCF_PATH,
               '-d', '10', '-d', '20']
        stderr = NamedTemporaryFile()
        stdout = NamedTemporaryFile()
        try:
            check_call(cmd, stderr=stderr, stdout=stdout)
            #raw_input()
        except CalledProcessError:
            sys.stderr.write(open(stderr.name).read())
            sys.stdout.write(open(stdout.name).read())
        finally:
            tempdir.close()
        #FREEBAYES
        tempdir = TemporaryDir()
        cmd = [binary, '-r', REF_PATH, '-o', tempdir.name, FREEBAYES_VCF_PATH,
               '-d', '10']
        stderr = NamedTemporaryFile()
        stdout = NamedTemporaryFile()
        try:
            check_call(cmd, stderr=stderr, stdout=stdout)
        except CalledProcessError:
            sys.stderr.write(open(stderr.name).read())
            sys.stdout.write(open(stdout.name).read())
        finally:
            tempdir.close()


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
    import sys;sys.argv = ['', 'AlleleCount2DTest', 'TestVcfStats']
    unittest.main()
