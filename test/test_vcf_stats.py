import unittest
from os.path import join
import sys

from vcf import Reader
from vcf_crumbs.vcf_stats import (calc_density_per_chrom, get_data_from_vcf,
                                  get_snpcaller_name, VARSCAN, GATK,
                                  calculate_maf, FREEBAYES, VcfStats,
                                  calc_n_bases_in_chrom_with_snp)

from vcf_crumbs.utils import TEST_DATA_DIR, BIN_DIR
from subprocess import check_call, CalledProcessError
from tempfile import NamedTemporaryFile
from crumbs.utils.file_utils import TemporaryDir
from crumbs.statistics import IntCounter
from cyvcf.parser import HOM_REF
from bam_crumbs.plot import draw_scatter

VARSCAN_VCF_PATH = join(TEST_DATA_DIR, 'sample.vcf.gz')
REF_PATH = join(TEST_DATA_DIR, 'sample_ref.fasta')
GATK_VCF_PATH = join(TEST_DATA_DIR, 'gatk_sample.vcf.gz')
FREEBAYES_VCF_PATH = join(TEST_DATA_DIR, 'freebayes_sample.vcf.gz')
FREEBAYES_MULTI_VCF_PATH = join(TEST_DATA_DIR, 'freebayes_multisample.vcf.gz')


class TestVcfStats(unittest.TestCase):
    def test_vcf_stats(self):
        vcf_stats = VcfStats(VARSCAN_VCF_PATH, gq_threshold=0)
        assert vcf_stats.het_by_sample == \
                        {'upv196': IntCounter({'num_gt': 90, 'num_het': 19}),
                         'pepo': IntCounter({'num_gt': 29, 'num_het': 8}),
                         'mu16': IntCounter({'num_gt': 107, 'num_het': 26})}
        vcf_stats = VcfStats(VARSCAN_VCF_PATH, gq_threshold=25)
        assert vcf_stats.het_by_sample == \
                        {'upv196': IntCounter({'num_gt': 90, 'num_het': 7}),
                         'pepo': IntCounter({'num_gt': 29, 'num_het': 4}),
                         'mu16': IntCounter({'num_gt': 107, 'num_het': 18})}
        assert len(vcf_stats.samples) == 3
        assert len(vcf_stats.mafs) == 4
        assert vcf_stats.mafs['all'][63] == 11
        assert vcf_stats.snps_per_chromosome['CUUC00423_TC01'] == 26
        assert len(vcf_stats.call_data[HOM_REF]['x']) == 95
        assert vcf_stats.genotype_qualities.count == 226
        assert vcf_stats.genotype_qualities[3] == 25
        assert vcf_stats.variable_gt_per_snp[100] == 50


class SnvStatTests(unittest.TestCase):

    def test_get_snpcaller(self):
        assert get_snpcaller_name(Reader(filename=VARSCAN_VCF_PATH)) == \
                                    VARSCAN
        assert get_snpcaller_name(Reader(filename=GATK_VCF_PATH)) == GATK

        assert get_snpcaller_name(Reader(filename=FREEBAYES_VCF_PATH)) == \
                                                                FREEBAYES

    def test_get_data(self):
        data = get_data_from_vcf(VARSCAN_VCF_PATH, gq_threshold=0)
        assert data['samples'] == set(['upv196', 'pepo', 'mu16'])
        assert data['het_by_sample'] == \
                        {'upv196': IntCounter({'num_gt': 90, 'num_het': 19}),
                         'pepo': IntCounter({'num_gt': 29, 'num_het': 8}),
                         'mu16': IntCounter({'num_gt': 107, 'num_het': 26})}

        data = get_data_from_vcf(VARSCAN_VCF_PATH, gq_threshold=25)
        assert data['het_by_sample'] == \
                        {'upv196': IntCounter({'num_gt': 90, 'num_het': 7}),
                         'pepo': IntCounter({'num_gt': 29, 'num_het': 4}),
                         'mu16': IntCounter({'num_gt': 107, 'num_het': 18})}

        assert len(data['variable_gt_per_snp']) == 81

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


class StatBinTests(unittest.TestCase):

    def test_draw_snv_stats_bin(self):
        binary = join(BIN_DIR, 'draw_snv_stats')
        tempdir = TemporaryDir()
        cmd = [binary, '-r', REF_PATH, '-o', tempdir.name, VARSCAN_VCF_PATH]
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
        cmd = [binary, '-r', REF_PATH, '-o', tempdir.name, FREEBAYES_VCF_PATH]
        stderr = NamedTemporaryFile()
        stdout = NamedTemporaryFile()
        try:
            check_call(cmd, stderr=stderr, stdout=stdout)
        except CalledProcessError:
            sys.stderr.write(open(stderr.name).read())
            sys.stdout.write(open(stdout.name).read())
        finally:
            tempdir.close()


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'SnvStatTests.test_calc_n_bases_in_chrom_with_snp']
    unittest.main()
