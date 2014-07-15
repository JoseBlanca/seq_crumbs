from os.path import join, dirname
import unittest
from tempfile import NamedTemporaryFile
from subprocess import check_output

from vcf_crumbs.snv import VCFReader

from vcf_crumbs.filters import (PASSED, FILTERED_OUT, group_in_filter_packets,
                                CallRateFilter, BiallelicFilter, IsSNPFilter,
                                GenotypeQualFilter, ObsHetFilter)
from vcf_crumbs.utils import TEST_DATA_DIR


VCF_PATH = join(TEST_DATA_DIR, 'sample.vcf.gz')
VCF_INDEL_PATH = join(TEST_DATA_DIR, 'sample_indel.vcf.gz')
VARI_VCF_PATH = join(TEST_DATA_DIR, 'vari_filter.vcf')


class FiltersTest(unittest.TestCase):
    def test_group_in_filter_packets(self):
        items = list(range(10))
        packets = list(group_in_filter_packets(items, 4))
        res = [{FILTERED_OUT: [], PASSED: (0, 1, 2, 3)},
               {FILTERED_OUT: [], PASSED: (4, 5, 6, 7)},
               {FILTERED_OUT: [], PASSED: (8, 9)}]
        assert res == packets

    @staticmethod
    def eval_prop_in_packet(packet, prop):
        eval_prop = lambda snps: [getattr(snp, prop)for snp in snps]
        packet = {PASSED: eval_prop(packet[PASSED]),
                  FILTERED_OUT: eval_prop(packet[FILTERED_OUT])}
        return packet

    @staticmethod
    def filter_vcf(vcf_fhand, filter_):
        snps = VCFReader(vcf_fhand, min_calls_for_pop_stats=1).parse_snvs()
        packet = list(group_in_filter_packets(snps, 10))[0]
        filtered_packet = filter_(packet)
        return filtered_packet

    def test_missing_genotypes(self):
        packet = self.filter_vcf(open(VCF_PATH),
                                 filter_=CallRateFilter(min_calls=2))
        res = self.eval_prop_in_packet(packet, 'num_called')
        expected = {FILTERED_OUT: [0, 1, 1, 1, 1], PASSED: [2, 2, 2, 2, 2]}
        assert res == expected

        packet = self.filter_vcf(open(VCF_PATH),
                                 filter_=CallRateFilter(min_calls=1))
        res = self.eval_prop_in_packet(packet, 'num_called')
        expected = {FILTERED_OUT: [0], PASSED: [2, 2, 2, 2, 1, 2, 1, 1, 1]}
        assert res == expected

        packet = self.filter_vcf(open(VCF_PATH),
                                 filter_=CallRateFilter(min_calls=1,
                                                        reverse=True))
        res = self.eval_prop_in_packet(packet, 'num_called')
        expected = {PASSED: [0], FILTERED_OUT: [2, 2, 2, 2, 1, 2, 1, 1, 1]}
        assert res == expected

        packet = self.filter_vcf(open(VCF_PATH),
                                 filter_=CallRateFilter(min_calls=3))
        res = self.eval_prop_in_packet(packet, 'num_called')
        expected = {PASSED: [], FILTERED_OUT: [2, 2, 2, 2, 0, 1, 2, 1, 1, 1]}
        assert res == expected

    def test_biallelic(self):
        packet = self.filter_vcf(open(VCF_PATH),
                                 filter_=BiallelicFilter())
        res = self.eval_prop_in_packet(packet, 'alleles')
        assert not(res[FILTERED_OUT])
        assert len(res[PASSED]) == 10

        packet = self.filter_vcf(open(VCF_INDEL_PATH),
                                 filter_=BiallelicFilter())
        res = self.eval_prop_in_packet(packet, 'alleles')
        assert len(res[FILTERED_OUT]) == 1
        assert len(res[PASSED]) == 6

    def test_is_snp(self):
        packet = self.filter_vcf(open(VCF_PATH),
                                 filter_=IsSNPFilter())
        res = self.eval_prop_in_packet(packet, 'is_snp')
        assert not(res[FILTERED_OUT])
        assert res[PASSED] == [True] * 10

        packet = self.filter_vcf(open(VCF_INDEL_PATH),
                                 filter_=IsSNPFilter())
        res = self.eval_prop_in_packet(packet, 'is_snp')
        assert res[FILTERED_OUT] == [False] * 7
        assert not res[PASSED]

    def test_gt_qual(self):
        packet = self.filter_vcf(open(VCF_PATH),
                                 filter_=GenotypeQualFilter(20))
        res = self.eval_prop_in_packet(packet, 'qual')
        assert res[FILTERED_OUT] == [None] * 10
        assert not res[PASSED]

        fpath = join(TEST_DATA_DIR, 'freebayes_al_depth.vcf')
        packet = self.filter_vcf(open(fpath),
                                 filter_=GenotypeQualFilter(1))
        res = self.eval_prop_in_packet(packet, 'qual')
        assert len(res[FILTERED_OUT]) == 1
        assert len(res[PASSED]) == 4

    def test_obs_het(self):
        packet = self.filter_vcf(open(VCF_PATH),
                                 filter_=ObsHetFilter(min_het=0.5))
        res = self.eval_prop_in_packet(packet, 'obs_het')
        assert res[FILTERED_OUT] == [0.0, 0.0, 0.0, None, 0.0, 0.0, 0.0]
        assert res[PASSED] == [0.5, 1.0, 0.5]

        packet = self.filter_vcf(open(VCF_PATH),
                                 filter_=ObsHetFilter(max_het=0.5))
        res = self.eval_prop_in_packet(packet, 'obs_het')
        assert res[FILTERED_OUT] == [None, 1.0]
        assert res[PASSED] == [0.0, 0.0, 0.0, 0.5, 0.5, 0.0, 0.0, 0.0]

        packet = self.filter_vcf(open(VCF_PATH),
                                 filter_=ObsHetFilter(min_het=0.1,
                                                      max_het=0.9))
        res = self.eval_prop_in_packet(packet, 'obs_het')
        assert res[FILTERED_OUT] == [0.0, 0.0, 0.0, None, 1.0, 0.0, 0.0, 0.0]
        assert res[PASSED] == [0.5, 0.5]


# TODO fix binary
class BinaryTest(unittest.TestCase):

    def test_bin_record_filters(self):
        binary = join(dirname(__file__), '..', 'bin', 'filter_snvs')
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
                gq = data.get(GQ, None)
                if gq is None or gq < 20:
                    assert not call.called

if __name__ == "__main__":
    unittest.main()
