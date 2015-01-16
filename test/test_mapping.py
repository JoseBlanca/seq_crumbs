# Copyright 2012 Jose Blanca, Peio Ziarsolo, COMAV-Univ. Politecnica Valencia
# This file is part of seq_crumbs.
# seq_crumbs is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.

# seq_crumbs is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR  PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with seq_crumbs. If not, see <http://www.gnu.org/licenses/>.

import unittest
import subprocess
import os.path
from tempfile import NamedTemporaryFile

from crumbs.utils.test_utils import TEST_DATA_DIR
from crumbs.mapping import (get_or_create_bowtie2_index, _bowtie2_index_exists,
                            map_with_bowtie2, get_or_create_bwa_index,
                            _bwa_index_exists, map_with_bwamem,
                            map_process_to_bam, sort_fastx_files,
                            map_with_tophat)
from crumbs.utils.file_utils import TemporaryDir
from crumbs.utils.bin_utils import get_binary_path
from crumbs.seq.seq import get_name
import pysam


# pylint: disable=R0201
# pylint: disable=R0904
# pylint: disable=W0402
# pylint: disable=C0111


class Bowtie2Test(unittest.TestCase):
    def test_get_or_create_index(self):
        db_name = 'arabidopsis_genes'
        seq_fpath = os.path.join(TEST_DATA_DIR, db_name)
        assert _bowtie2_index_exists(seq_fpath)

        directory = TemporaryDir()
        index_fpath = get_or_create_bowtie2_index(seq_fpath, directory.name)
        expected_index = os.path.join(directory.name,
                                      os.path.basename(db_name))
        assert index_fpath == expected_index
        assert _bowtie2_index_exists(index_fpath)

        # already exists
        index_fpath = get_or_create_bowtie2_index(seq_fpath, directory.name)
        assert index_fpath == expected_index
        assert _bowtie2_index_exists(index_fpath)
        directory.close()

    def test_map_with_bowtie2(self):
        reference_fpath = os.path.join(TEST_DATA_DIR, 'arabidopsis_genes')
        reads_fpath = os.path.join(TEST_DATA_DIR, 'arabidopsis_reads.fastq')
        directory = TemporaryDir()
        index_fpath = get_or_create_bowtie2_index(reference_fpath,
                                                  directory.name)
        bam_fhand = NamedTemporaryFile(suffix='.bam')
        bowtie2 = map_with_bowtie2(index_fpath, unpaired_fpath=reads_fpath)
        map_process_to_bam(bowtie2, bam_fhand.name)
        directory.close()

        #With paired_fpahts option
        reference_fpath = os.path.join(TEST_DATA_DIR, 'arabidopsis_genes')
        forward_fpath = os.path.join(TEST_DATA_DIR, 'arabidopsis_reads.fastq')
        reverse_fpath = NamedTemporaryFile().name
        paired_fpaths = (forward_fpath, reverse_fpath)
        directory = TemporaryDir()
        index_fpath = get_or_create_bowtie2_index(reference_fpath,
                                                  directory.name)
        bam_fhand = NamedTemporaryFile(suffix='.bam')
        bowtie2 = map_with_bowtie2(index_fpath, paired_fpaths=paired_fpaths)
        map_process_to_bam(bowtie2, bam_fhand.name)
        directory.close()

    def test_rev_compl_fragmented_reads(self):
        index_fpath = os.path.join(TEST_DATA_DIR, 'ref_example.fasta')

        #with unpaired_reads
        query_f = '>seq1\nAAGTTCAATGTAGCAAGGGTACATGCTGACGAGATGGGTCTGCGATCCCTG'
        query_f += 'AGGACACCCAGTCTCCCGGGAGTCTTTTCCAAGGTGTGCTCCTGATCGCCGTGTTA\n'

        query_r = '>seq2\nTAACACGGCGATCAGGAGCACACCTTGGAAAAGACTCCCGGGAGACTGGGTG'
        query_r += 'TCCTCAGGGATCGCAGACCCATCTCGTCAGCATGTACCCTTGCTACATTGAACTT\n'

        query = query_f + query_r
        in_fhand = NamedTemporaryFile()
        in_fhand.write(query)
        in_fhand.flush()

        bam_fhand = NamedTemporaryFile(suffix='.bam')
        bowtie2 = map_with_bowtie2(index_fpath, extra_params=['-a', '-f'],
                                   unpaired_fpath=in_fhand.name)
        map_process_to_bam(bowtie2, bam_fhand.name)
        samfile = pysam.Samfile(bam_fhand.name)
        #for aligned_read in samfile:
        #    print aligned_read

        #with paired_reads.
        #f is reversed r is direct
        query1 = '>seq10 f\nGGGATCGCAGACCCATCTCGTCAGCATGTACCCTTGCTACATTGAACTT'
        query1 += '\n'
        query2 = '>seq10 r\nATGTAATACGGGCTAGCCGGGGATGCCGACGATTAAACACGCTGTCATA'
        query2 += 'GTAGCGTCTTCCTAGGGTTTTCCCCATGGAATCGGTTATCGTGATACGTTAAATTT\n'
        #f is direct, r is reversed
        query3 = '>seq11 f\nAAGTTCAATGTAGCAAGGGTACATGCTGACGAGATGGGTCTGCGATCCC'
        query3 += '\n'
        query4 = '>seq11 r\nAAATTTAACGTATCACGATAACCGATTCCATGGGGAAAACCCTAGGAAG'
        query4 += 'ACGCTACTATGACAGCGTGTTTAATCGTCGGCATCCCCGGCTAGCCCGTATTACAT\n'

        query_f = query1 + query3
        query_r = query2 + query4

        f_fhand = NamedTemporaryFile()
        f_fhand.write(query_f)
        f_fhand.flush()
        r_fhand = NamedTemporaryFile()
        r_fhand.write(query_r)
        r_fhand.flush()
        paired_fpaths = (f_fhand.name, r_fhand.name)

        bam_fhand = NamedTemporaryFile(suffix='.bam')
        bowtie2 = map_with_bowtie2(index_fpath, extra_params=['-a', '-f'],
                                   paired_fpaths=paired_fpaths)
        map_process_to_bam(bowtie2, bam_fhand.name)
        samfile = pysam.Samfile(bam_fhand.name)
        #for aligned_read in samfile:
        #    print aligned_read

    def test_tophat(self):
        reference_fpath = os.path.join(TEST_DATA_DIR, 'arabidopsis_genes')
        reads_fpath = os.path.join(TEST_DATA_DIR, 'arabidopsis_reads.fastq')
        directory = TemporaryDir()
        index_fpath = get_or_create_bowtie2_index(reference_fpath,
                                                  directory.name)
        map_with_tophat(index_fpath, directory.name,
                        unpaired_fpath=reads_fpath)
        os.path.exists(os.path.join(directory.name, 'accepted_hits.bam'))
        directory.close()

    def test_tophat_paired(self):
        reference_fpath = os.path.join(TEST_DATA_DIR, 'arabidopsis_genes')
        reads_1_fpath = os.path.join(TEST_DATA_DIR, 'reads_1.fastq')
        reads_2_fpath = os.path.join(TEST_DATA_DIR, 'reads_2.fastq')
        try:
            directory = TemporaryDir()
            index_fpath = get_or_create_bowtie2_index(reference_fpath,
                                                      directory.name)
            map_with_tophat(index_fpath, directory.name,
                            paired_fpaths=[reads_1_fpath, reads_2_fpath])
            os.path.exists(os.path.join(directory.name, 'accepted_hits.bam'))
            self.fail('runtimeError expected')
        except RuntimeError:
            pass
        finally:
            directory.close()

        try:
            directory = TemporaryDir()
            index_fpath = get_or_create_bowtie2_index(reference_fpath,
                                                      directory.name)
            map_with_tophat(index_fpath, directory.name,
                            paired_fpaths=[reads_1_fpath, reads_2_fpath],
                            mate_inner_dist=350, mate_std_dev=50)
            os.path.exists(os.path.join(directory.name, 'accepted_hits.bam'))
        finally:
            directory.close()


class Bwa2Test(unittest.TestCase):
    def test_get_or_create_index(self):
        db_name = 'arabidopsis_genes'
        seq_fpath = os.path.join(TEST_DATA_DIR, db_name)
        assert not _bwa_index_exists(seq_fpath)

        directory = TemporaryDir()
        index_fpath = get_or_create_bwa_index(seq_fpath, directory.name)
        expected_index = os.path.join(directory.name,
                                      os.path.basename(db_name))
        assert index_fpath == expected_index
        assert _bwa_index_exists(index_fpath)

        # already exists
        index_fpath = get_or_create_bwa_index(seq_fpath, directory.name)
        assert index_fpath == expected_index
        assert _bwa_index_exists(index_fpath)
        directory.close()

    def test_map_with_bwa(self):
        reference_fpath = os.path.join(TEST_DATA_DIR, 'arabidopsis_genes')
        reads_fpath = os.path.join(TEST_DATA_DIR, 'arabidopsis_reads.fastq')
        directory = TemporaryDir()
        index_fpath = get_or_create_bwa_index(reference_fpath, directory.name)
        bam_fhand = NamedTemporaryFile(suffix='.bam')
        bwa = map_with_bwamem(index_fpath, unpaired_fpath=reads_fpath)
        map_process_to_bam(bwa, bam_fhand.name)
        out = subprocess.check_output([get_binary_path('samtools'), 'view',
                                       bam_fhand.name])
        assert  'TTCTGATTCAATCTACTTCAAAGTTGGCTTTATCAATAAG' in out

        directory.close()

    def test_rev_compl_fragmented_reads(self):
        index_fpath = os.path.join(TEST_DATA_DIR, 'ref_example.fasta')

        #with paired_reads.
        #f is reversed r is direct
        query1 = '>seq10 f\nGGGATCGCAGACCCATCTCGTCAGCATGTACCCTTGCTACATTGAACTT'
        query1 += '\n'
        query2 = '>seq10 r\nATGTAATACGGGCTAGCCGGGGATGCCGACGATTAAACACGCTGTCATA'
        query2 += 'GTAGCGTCTTCCTAGGGTTTTCCCCATGGAATCGGTTATCGTGATACGTTAAATTT\n'
        #f is direct, r is reversed
        query3 = '>seq11 f\nAAGTTCAATGTAGCAAGGGTACATGCTGACGAGATGGGTCTGCGATCCC'
        query3 += '\n'
        query4 = '>seq11 r\nAAATTTAACGTATCACGATAACCGATTCCATGGGGAAAACCCTAGGAAG'
        query4 += 'ACGCTACTATGACAGCGTGTTTAATCGTCGGCATCCCCGGCTAGCCCGTATTACAT\n'

        #f is fragmented in two reference sequences. r mapps completely
        query7 = '>seq4 f\nCAAATCATCACCAGACCATGTCCGATCCCGGGAGTCTTTTCCAAGGTGTGC'
        query7 += 'TCTTTATCCGGCCCTTGCTCAAGGGTATGTTAAAACGGCAAGAGCTGCCTGAGCGCG\n'
        query8 = '>seq4 r\nTGTTCTGCAATCGATACAACGATCGAATTTAATCTGAGTAACTGCCAATTC'
        query8 += 'TGAGTAATATTATAGAAAGT\n'

        query_f = query1 + query3 + query7
        query_r = query2 + query4 + query8

        f_fhand = NamedTemporaryFile()
        f_fhand.write(query_f)
        f_fhand.flush()
        r_fhand = NamedTemporaryFile()
        r_fhand.write(query_r)
        r_fhand.flush()
        paired_fpaths = (f_fhand.name, r_fhand.name)

        bam_fhand = NamedTemporaryFile(suffix='.bam')
        bwa = map_with_bwamem(index_fpath, paired_fpaths=paired_fpaths)
        map_process_to_bam(bwa, bam_fhand.name)
        samfile = pysam.Samfile(bam_fhand.name)
        #for aligned_read in samfile:
        #    print aligned_read


class SortSeqsFileTest(unittest.TestCase):
    def test_sort_by_position_in_ref(self):
        index_fpath = os.path.join(TEST_DATA_DIR, 'ref_example.fasta')

        #with fasta format
        query1 = '>seq1\nGAGAATTAAGCCTATCTGGAGAGCGGTACCAACAGGGAAACACCGACTCA\n'
        query2 = '>seq2\nTAACAGTATGTGCCGTAGGGCCGTCGCCGCATCCACGTTATCGGAAGGGC\n'
        query3 = '>seq3\nTACGGCCGTCCCCCTGCTGCTTATCATCAGGCGACGATAGTCAGCTCCGC\n'
        query4 = '>seq4\nTGCAGAGACCGACATGCGAAAGGAGTGACTATCACCGTCAATGGCGTGCC\n'
        query5 = '>seq5\nAATAAATAATCTGGGTATGTACTCGGAGTCTACGTAAGCGCGCTTAAATT\n'
        query6 = '>seq6\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n'
        query = query6 + query1 + query2 + query3 + query4 + query5
        in_fhand = NamedTemporaryFile()
        in_fhand.write(query)
        in_fhand.flush()

        sorted_names = []
        for seq in sort_fastx_files(in_fhand, 'coordinate', index_fpath):
            sorted_names.append(get_name(seq))
        expected_names = ['seq2', 'seq3', 'seq1', 'seq5', 'seq4', 'seq6']
        assert sorted_names == expected_names
        #it fails because bwa somehow gives a position to an unmapped seq

        #with fastq format
        query1 += '+\n??????????????????????????????????????????????????\n'
        query2 += '+\n??????????????????????????????????????????????????\n'
        query3 += '+\n??????????????????????????????????????????????????\n'
        query4 += '+\n??????????????????????????????????????????????????\n'
        query5 += '+\n??????????????????????????????????????????????????\n'
        query6 += '+\n??????????????????????????????????????????????????\n'
        query = query6 + query1 + query2 + query3 + query4 + query5
        in_fhand = NamedTemporaryFile()
        in_fhand.write(query)
        in_fhand.flush()

        sorted_names = []
        for seq in sort_fastx_files(in_fhand, 'coordinate', index_fpath):
            sorted_names.append(get_name(seq))
        expected_names = ['seq2', 'seq3', 'seq1', 'seq5', 'seq4', 'seq6']
        assert sorted_names == expected_names

        #sort by sequence
        sorted_names = []
        for seq in sort_fastx_files([in_fhand], key='seq', directory=None,
                     max_items_in_memory=None, tempdir=None):
            sorted_names.append(get_name(seq))
        expected_names = ['seq6', 'seq5', 'seq1', 'seq2', 'seq3', 'seq4']
        assert sorted_names == expected_names

    def test_add_rg_to_bam(self):
        reference_fpath = os.path.join(TEST_DATA_DIR, 'arabidopsis_genes')
        reads_fpath = os.path.join(TEST_DATA_DIR, 'arabidopsis_reads.fastq')
        directory = TemporaryDir()
        index_fpath = get_or_create_bwa_index(reference_fpath, directory.name)
        bam_fhand = NamedTemporaryFile(suffix='.bam')
        lib_name = 'aa'
        log_fhand = NamedTemporaryFile()
        readgroup = {'ID': lib_name, 'PL': 'illumina', 'LB': lib_name,
                     'SM': '{0}_illumina_pe'.format(lib_name), 'PU': '0'}
        bwa = map_with_bwamem(index_fpath, unpaired_fpath=reads_fpath,
                              readgroup=readgroup, log_fpath=log_fhand.name)
        map_process_to_bam(bwa, bam_fhand.name)
        out = subprocess.check_output([get_binary_path('samtools'), 'view',
                                       '-h', bam_fhand.name], stderr=log_fhand)
        assert '@RG\tID:aa' in out
        assert 'TTCTGATTCAATCTACTTCAAAGTTGGCTTTATCAATAAG' in out

        directory.close()
if __name__ == '__main__':
    # import sys;sys.argv = ['', 'SortSeqsFileTest.test_add_rg_to_bam']
    unittest.main()
