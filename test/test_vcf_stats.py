import unittest
from os.path import join


from vcf import Reader
from vcf_crumbs.vcf_stats import (calc_density_per_chrom, get_data_from_vcf,
                                  get_snpcaller_name, VARSCAN, GATK,
                                  calculate_maf)

from vcf_crumbs.utils import TEST_DATA_DIR

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

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
