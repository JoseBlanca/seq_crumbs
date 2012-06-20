'''
Created on 18/06/2012

@author: jose
'''

import unittest
import os.path
from tempfile import NamedTemporaryFile

from crumbs.sff_extract import SffExtractor
from crumbs.tests.utils import TEST_DATA_DIR

# pylint: disable=R0201
# pylint: disable=R0904


class SffExtractTest(unittest.TestCase):
    'It tests the SFF extraction functionality.'

    def test_seqs_in_sff(self):
        'It extracts an sff file'

        sff_fpath = os.path.join(TEST_DATA_DIR, '10_454_reads.sff')

        # No trim
        seqs = SffExtractor([open(sff_fpath, 'rb')]).seqs
        seqs = list(seqs)
        assert len(seqs) == 10
        assert str(seqs[0].seq).startswith('tcagGGTCTACATGTTGGTTAACCCGTACTGAT')

        # trimmed
        seqs = SffExtractor([open(sff_fpath, 'rb')], trim=True).seqs
        seqs = list(seqs)
        assert len(seqs) == 10
        assert str(seqs[0].seq).startswith('GGTCTACATGTTGGTTAACCCGTACTGATTTGA')

        # min left_clip and no trim
        seqs = SffExtractor([open(sff_fpath, 'rb')], min_left_clip=5).seqs
        seqs = list(seqs)
        assert len(seqs) == 10
        assert str(seqs[0].seq).startswith('tcaggGTCTACATGTTGGTTAACCCGTACTGAT')

        # min left_clip and trim
        seqs = SffExtractor([open(sff_fpath, 'rb')], min_left_clip=5,
                            trim=True).seqs
        seqs = list(seqs)
        assert len(seqs) == 10
        assert str(seqs[0].seq).startswith('GTCTACATGTTGGTTAACCCGTACTGAT')

        # empty file
        empty_fhand = NamedTemporaryFile()
        try:
            seqs = SffExtractor([open(empty_fhand.name, 'rb')]).seqs
            list(seqs)
            self.fail('ValueError expected.')
        except ValueError:
            pass

        # Wrong file type
        fasta_fhand = NamedTemporaryFile()
        fasta_fhand.write('>a_seq\n' + 'ACTG' * 30 + '\n')
        fasta_fhand.flush()
        try:
            seqs = SffExtractor([open(fasta_fhand.name, 'rb')]).seqs
            list(seqs)
        except ValueError:
            pass

    def test_check_nucl_counts(self):
        'It checks that the nucleotide freqs are all below the given threshold'
        sff_fpath = os.path.join(TEST_DATA_DIR, '10_454_reads.sff')

        extractor = SffExtractor([open(sff_fpath, 'rb')])
        seqs = extractor.seqs
        seqs = list(seqs)
        assert len(seqs) == 10
        assert extractor.clip_advice[sff_fpath][0] == 5
        assert extractor.clip_advice[sff_fpath][1] == 'A'

        extractor = SffExtractor([open(sff_fpath, 'rb')], min_left_clip=4,
                            trim=False)
        seqs = extractor.seqs
        seqs = list(seqs)
        assert len(seqs) == 10
        assert extractor.clip_advice[sff_fpath][0] == 5

        extractor = SffExtractor([open(sff_fpath, 'rb')], min_left_clip=4,
                            trim=True)
        seqs = extractor.seqs
        seqs = list(seqs)
        assert len(seqs) == 10
        assert extractor.clip_advice[sff_fpath][0] == 5

        extractor = SffExtractor([open(sff_fpath, 'rb')], min_left_clip=5,
                            trim=True)
        seqs = extractor.seqs
        seqs = list(seqs)
        assert len(seqs) == 10
        assert not extractor.clip_advice[sff_fpath]

if __name__ == '__main__':
    #import sys;sys.argv = ['', 'SffExtractTest.test_items_in_gff']
    unittest.main()
