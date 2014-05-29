from os.path import join, dirname
import unittest
from tempfile import NamedTemporaryFile
from subprocess import check_output

from vcf import Reader

from vcf_crumbs.filters import (GenotypesInSamplesFilter, AlleleNumberFilter,
                                MissingGenotypesFilter, remove_low_quality_gt)
from vcf_crumbs.utils import TEST_DATA_DIR
from vcf_crumbs.vcf_stats import (VARSCAN,  get_call_data, GQ,
                                  get_snpcaller_name, GT)


VCF_PATH = join(TEST_DATA_DIR, 'sample.vcf.gz')
VCF_INDEL_PATH = join(TEST_DATA_DIR, 'sample_indel.vcf.gz')
VARI_VCF_PATH = join(TEST_DATA_DIR, 'vari_filter.vcf')


class FiltersTest(unittest.TestCase):
    def xtest_vcfparser(self):
        reader = Reader(open(VCF_INDEL_PATH))
        for snp in reader:
            print "***************"
            print snp.REF
            print type(snp.ALT[0])
            print('snv type: ' + snp.var_subtype)
            print('alleles: ' + ','.join([unicode(al) for al in snp.alleles]))
            print 'info', snp.INFO
            for call in snp.samples:
                if call.gt_bases is None:
                    continue
                print('call alleles: ' + ','.join(call.gt_alleles))
                print call.data
                print call.site
            print 'pos', snp.POS
            print('start: ' + str(snp.start))
            print 'end', snp.end
            print snp.is_snp

    def test_gt_in_samples(self):
        records = Reader(filename=VCF_PATH)
        rec1 = records.next()
        filter_ = GenotypesInSamplesFilter([0], ['mu16'])
        filter_.vcf_variant = VARSCAN
        assert filter_(rec1)

        filter_ = GenotypesInSamplesFilter([0], ['pepo', 'mu16'], 2)
        filter_.vcf_variant = VARSCAN
        assert not filter_(rec1)

    def test_allele_number(self):
        records = Reader(filename=VCF_PATH)
        rec1 = records.next()
        filter_ = AlleleNumberFilter(2)
        filter_.vcf_variant = VARSCAN
        assert filter_(rec1)

        filter_ = AlleleNumberFilter(1)
        filter_.vcf_variant = VARSCAN
        assert not filter_(rec1)

    def test_missing_genotypes(self):
        records = Reader(filename=VCF_PATH)
        rec1 = records.next()
        filter_ = MissingGenotypesFilter(2)
        filter_.vcf_variant = VARSCAN
        assert filter_(rec1)

        filter_ = MissingGenotypesFilter(0)
        filter_.vcf_variant = VARSCAN
        assert not filter_(rec1)

    def test_genotype_quality(self):
        records = Reader(filename=VCF_PATH)
        vcf_variant = get_snpcaller_name(records)
        rec1 = records.next()
        for call in remove_low_quality_gt(rec1, 20, vcf_variant):
            if call.sample != 'upv196':
                assert call.gt_bases is None


class BinaryTest(unittest.TestCase):

    def test_bin_record_filters(self):
        binary = join(dirname(__file__), '..', 'bin', 'run_vcf_record_filters')
        assert 'usage' in check_output([binary, '-h'])

        out_fhand = NamedTemporaryFile()
        in_fpath = VCF_PATH
        samples_fhand = NamedTemporaryFile()
        samples_fhand.write('mu16\n')
        samples_fhand.flush()
        cmd = [binary, in_fpath, '-o', out_fhand.name, '-g', '0', '-s',
               samples_fhand.name, '-t', '20']
        check_output(cmd)
        reader = Reader(filename=out_fhand.name)
        vcf_variant = get_snpcaller_name(reader)
        for snp in reader:
            gt = snp.genotype('mu16').gt_type
            assert gt == 0 or gt == None
            for call in snp:
                data = get_call_data(call, vcf_variant)
                if data[GQ] < 20:
                    assert data[GT] is None

if __name__ == "__main__":
    unittest.main()
