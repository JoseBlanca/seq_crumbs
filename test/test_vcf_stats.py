import unittest
from os.path import join
from tempfile import NamedTemporaryFile

from vcf import Reader
from vcf_crumbs.vcf_stats import (calc_density_per_chrom, draw_histogram,
                                  get_data_from_vcf, get_snpcaller_name,
                                  VARSCAN, GATK, calculate_maf, draw_scatter)
from vcf_crumbs.utils.path_utils import TEST_DATA_DIR

VARSCAN_VCF_PATH = join(TEST_DATA_DIR, 'sample.vcf.gz')
REF_PATH = join(TEST_DATA_DIR, 'sample_ref.fasta')
GATK_VCF_PATH = join(TEST_DATA_DIR, 'gatk_sample.vcf.gz')


class SnvStatTests(unittest.TestCase):

    def test_get_snpcaller(self):
        assert get_snpcaller_name(Reader(filename=VARSCAN_VCF_PATH)) == \
                                    VARSCAN
        assert get_snpcaller_name(Reader(filename=GATK_VCF_PATH)) == GATK

    def test_get_data(self):
        data = get_data_from_vcf(VARSCAN_VCF_PATH)
        #print data

    def test_calc_densities(self):
        data = get_data_from_vcf(VARSCAN_VCF_PATH)
        densities = calc_density_per_chrom(data['snps_per_chromo'],
                                           open(REF_PATH))
        assert densities['CUUC00355_TC01'] == 3.74

    def test_calc_maf(self):
        #varscan
        reader = Reader(filename=VARSCAN_VCF_PATH)
        snp = reader.next()
        maf = calculate_maf(snp, snpcaller=VARSCAN)
        assert 0.52 < maf < 0.53
        #gatk
        reader = Reader(filename=GATK_VCF_PATH)
        snp = reader.next()
        maf = calculate_maf(snp, snpcaller=GATK)
        assert 0.7 < maf < 0.72

    def test_do_histogram(self):
        values = [1, 2, 3, 1, 2, 3, 2, 3, 2, 3, 2, 1, 4]
        fhand = NamedTemporaryFile(suffix='.png')
        draw_histogram(values, fhand, bins=100, title='hola')
        #raw_input(fhand.name)

    def test_draw_scatter(self):
        data = get_data_from_vcf(VARSCAN_VCF_PATH)
        fhand = NamedTemporaryFile(suffix='.png')
        draw_scatter(data['call_data'].values(), fhand, xlim=0, ylim=0)
        #raw_input(fhand.name)

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
