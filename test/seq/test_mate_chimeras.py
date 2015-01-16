'''
Created on 2014 uzt 10

@author: peio
'''
import unittest
import os.path
from subprocess import check_output
from tempfile import NamedTemporaryFile


from crumbs.seq.mate_chimeras import (classify_mapped_reads, classify_chimeras,
                                      calculate_distance_distribution)
from crumbs.utils.bin_utils import BIN_DIR
from crumbs.utils.test_utils import TEST_DATA_DIR
from crumbs.utils.tags import NON_CHIMERIC, CHIMERA, UNKNOWN
from crumbs.seq.seq import get_name
from crumbs.mapping import map_with_bwamem, map_process_to_sortedbam


class FilterByMappingType(unittest.TestCase):
    def test_classify_paired_reads(self):
        index_fpath = os.path.join(TEST_DATA_DIR, 'ref_example.fasta')
        #Non chimeric
        query1 = '>seq1 1:N:0:GATCAG\nGGGATCGCAGACCCATCTCGTCAGCATGTACCCTTGCTACATTGAACTT\n'
        query2 = '>seq1 2:N:0:GATCAG\nAGGAGGGATCGGGCACCCACGGCGCGGTAGACTGAGGCCTTCTCGAACT\n'
        #Chimeric
        query3 = '>seq2 1:N:0:GATCAG\nAAGTTCAATGTAGCAAGGGTACATGCTGACGAGATGGGTCTGCGATCCC\n'
        query4 = '>seq2 2:N:0:GATCAG\nACGTGGATGCGGCGACGGCCCTACGGCACATACTGTTATTAGGGTCACT\n'
        #unknown
        query5 = '>seq3 1:N:0:GATCAG\nAGTGACCCTAATAACAGTATGTGCCGTAGGGCCGTCGCCGCATCCACGT\n'
        query6 = '>seq3 2:N:0:GATCAG\nGTCGTGCGCAGCCATTGAGACCTTCCTAGGGTTTTCCCCATGGAATCGG\n'

        query = query1 + query2 + query5 + query6 + query3 + query4
        in_fhand = NamedTemporaryFile()
        in_fhand.write(query)
        in_fhand.flush()

        bam_fhand = NamedTemporaryFile(suffix='.bam')
        extra_params = ['-a', '-M']
        bwa = map_with_bwamem(index_fpath, interleave_fpath=in_fhand.name,
                              extra_params=extra_params)
        map_process_to_sortedbam(bwa, bam_fhand.name, key='queryname')
        result = classify_mapped_reads(bam_fhand, mate_distance=2000)
        for pair, kind in result:
            if kind == NON_CHIMERIC:
                assert 'seq1' in get_name(pair[0])
            elif kind == UNKNOWN:
                assert 'seq3' in get_name(pair[0])
            elif kind == CHIMERA:
                assert 'seq2' in get_name(pair[0])
            else:
                self.fail()

    def test_filter_chimeras(self):
        index_fpath = os.path.join(TEST_DATA_DIR, 'ref_example.fasta')
        #Non chimeric
        query1 = '>seq1 1:N:0:GATCAG\nGGGATCGCAGACCCATCTCGTCAGCATGTACCCTTGCTACATTGAACTT\n'
        query2 = '>seq1 2:N:0:GATCAG\nAGGAGGGATCGGGCACCCACGGCGCGGTAGACTGAGGCCTTCTCGAACT\n'
        #Chimeric
        query3 = '>seq2 1:N:0:GATCAG\nAAGTTCAATGTAGCAAGGGTACATGCTGACGAGATGGGTCTGCGATCCC\n'
        query4 = '>seq2 2:N:0:GATCAG\nACGTGGATGCGGCGACGGCCCTACGGCACATACTGTTATTAGGGTCACT\n'
        #unknown
        query5 = '>seq3 1:N:0:GATCAG\nAGTGACCCTAATAACAGTATGTGCCGTAGGGCCGTCGCCGCATCCACGT\n'
        query6 = '>seq3 2:N:0:GATCAG\nGTCGTGCGCAGCCATTGAGACCTTCCTAGGGTTTTCCCCATGGAATCGG\n'

        query = query1 + query2 + query5 + query6 + query3 + query4
        in_fhand = NamedTemporaryFile()
        in_fhand.write(query)
        in_fhand.flush()

        #classify_chimeras function
        out_fhand = NamedTemporaryFile()
        chimeras_fhand = NamedTemporaryFile()
        unknown_fhand = NamedTemporaryFile()
        classify_chimeras(in_fhand, index_fpath, mate_distance=2000,
                        out_fhand=out_fhand, chimeras_fhand=chimeras_fhand,
                        unknown_fhand=unknown_fhand)
        out_fhand.flush()
        chimeras_fhand.flush()
        unknown_fhand.flush()
        assert 'seq1' in open(out_fhand.name).next()
        assert 'seq2' in open(chimeras_fhand.name).next()
        assert 'seq3' in open(unknown_fhand.name).next()

    def test_filter_chimeras_bin(self):
        index_fpath = os.path.join(TEST_DATA_DIR, 'ref_example.fasta')
        #Non chimeric
        query1 = '>seq1 1:N:0:GATCAG\nGGGATCGCAGACCCATCTCGTCAGCATGTACCCTTGCTACATTGAACTT\n'
        query2 = '>seq1 2:N:0:GATCAG\nAGGAGGGATCGGGCACCCACGGCGCGGTAGACTGAGGCCTTCTCGAACT\n'
        #Chimeric
        query3 = '>seq2 1:N:0:GATCAG\nAAGTTCAATGTAGCAAGGGTACATGCTGACGAGATGGGTCTGCGATCCC\n'
        query4 = '>seq2 2:N:0:GATCAG\nACGTGGATGCGGCGACGGCCCTACGGCACATACTGTTATTAGGGTCACT\n'
        #unknown
        query5 = '>seq3 1:N:0:GATCAG\nAGTGACCCTAATAACAGTATGTGCCGTAGGGCCGTCGCCGCATCCACGT\n'
        query6 = '>seq3 2:N:0:GATCAG\nGTCGTGCGCAGCCATTGAGACCTTCCTAGGGTTTTCCCCATGGAATCGG\n'

        query = query1 + query2 + query5 + query6 + query3 + query4
        in_fhand = NamedTemporaryFile()
        in_fhand.write(query)
        in_fhand.flush()

        filter_chimeras_bin = os.path.join(BIN_DIR, 'classify_chimeras')
        assert 'usage' in check_output([filter_chimeras_bin, '-h'])
        chimeras_fhand = NamedTemporaryFile()
        unknown_fhand = NamedTemporaryFile()
        out_fhand = NamedTemporaryFile()
        cmd = [filter_chimeras_bin, in_fhand.name, '-r', index_fpath]
        cmd.extend(['-c', chimeras_fhand.name, '-u', unknown_fhand.name,
                    '-s', '2000', '-o', out_fhand.name])
        check_output(cmd, stdin=in_fhand)
        assert 'seq1' in open(out_fhand.name).next()
        assert 'seq2' in open(chimeras_fhand.name).next()
        assert 'seq3' in open(unknown_fhand.name).next()


class DrawDistanceDistribution(unittest.TestCase):
    def test_calculate_mp_distance_distribution(self):
        index_fpath = os.path.join(TEST_DATA_DIR, 'ref_example.fasta')
        query1 = '>seq1 1:N:0:GATCAG\n'
        query1 += 'GGGATCGCAGACCCATCTCGTCAGCATGTACCCTTGCTACATTGAACTT\n'
        query2 = '>seq1 2:N:0:GATCAG\n'
        query2 += 'AGGAGGGATCGGGCACCCACGGCGCGGTAGACTGAGGCCTTCTCGAACT\n'
        # Chimeric
        query3 = '>seq2 1:N:0:GATCAG\n'
        query3 += 'AAGTTCAATGTAGCAAGGGTACATGCTGACGAGATGGGTCTGCGATCCC\n'
        query4 = '>seq2 2:N:0:GATCAG\n'
        query4 += 'ACGTGGATGCGGCGACGGCCCTACGGCACATACTGTTATTAGGGTCACT\n'
        # unknown
        query5 = '>seq3 1:N:0:GATCAG\n'
        query5 += 'AGTGACCCTAATAACAGTATGTGCCGTAGGGCCGTCGCCGCATCCACGT\n'
        query6 = '>seq3 2:N:0:GATCAG\n'
        query6 += 'GTCGTGCGCAGCCATTGAGACCTTCCTAGGGTTTTCCCCATGGAATCGG\n'

        query = query1 + query2 + query5 + query6 + query3 + query4
        in_fhand = NamedTemporaryFile()
        in_fhand.write(query)
        in_fhand.flush()
        stats = calculate_distance_distribution(in_fhand, index_fpath,
                                                   max_clipping=0.05)
        assert stats['outies'][1776] == 1
        assert stats['innies'][82] == 1
        assert stats['others'][1417] == 1

    def test_draw_distance_distribution_bin(self):
        index_fpath = os.path.join(TEST_DATA_DIR, 'ref_example.fasta')
        # Non chimeric
        query1 = '>seq1 1:N:0:GATCAG\n'
        query1 += 'GGGATCGCAGACCCATCTCGTCAGCATGTACCCTTGCTACATTGAACTT\n'
        query2 = '>seq1 2:N:0:GATCAG\n'
        query2 += 'AGGAGGGATCGGGCACCCACGGCGCGGTAGACTGAGGCCTTCTCGAACT\n'
        # Chimeric
        query3 = '>seq2 1:N:0:GATCAG\n'
        query3 += 'AAGTTCAATGTAGCAAGGGTACATGCTGACGAGATGGGTCTGCGATCCC\n'
        query4 = '>seq2 2:N:0:GATCAG\n'
        query4 += 'ACGTGGATGCGGCGACGGCCCTACGGCACATACTGTTATTAGGGTCACT\n'
        # unknown
        query5 = '>seq3 1:N:0:GATCAG\n'
        query5 += 'AGTGACCCTAATAACAGTATGTGCCGTAGGGCCGTCGCCGCATCCACGT\n'
        query6 = '>seq3 2:N:0:GATCAG\n'
        query6 += 'GTCGTGCGCAGCCATTGAGACCTTCCTAGGGTTTTCCCCATGGAATCGG\n'

        query = query1 + query2 + query5 + query6 + query3 + query4
        in_fhand = NamedTemporaryFile()
        in_fhand.write(query)
        in_fhand.flush()

        distribution_fhand = NamedTemporaryFile()
        draw_distances_distribution_bin = os.path.join(BIN_DIR,
                                               'draw_pair_distance_distribution')
        assert 'usage' in check_output([draw_distances_distribution_bin, '-h'])
        cmd = [draw_distances_distribution_bin, '-r', index_fpath, '-o',
               distribution_fhand.name, in_fhand.name]
        check_output(cmd)
        # raw_input(distribution_fhand.name)

if __name__ == "__main__":
    #import sys; sys.argv = ['', 'DrawDistanceDistribution']
    unittest.main()
