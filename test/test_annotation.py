'''
Created on 2012 urr 26

@author: peio
'''
import unittest
import os

from crumbs.annotation import EstscanOrfAnnotator
from crumbs.utils.test_utils import TEST_DATA_DIR
from crumbs.seqio import read_seqrecords


class AnnotationTest(unittest.TestCase):
    'Test annotator classes'

    def test_orf_annotator(self):
        'It tests orf annotator'
        fpath = os.path.join(TEST_DATA_DIR, 'orf_test.fasta')
        estscan_matrix = os.path.join(TEST_DATA_DIR,
                                      'Arabidopsis_thaliana.smat')
        seq_records = list(read_seqrecords([open(fpath)]))
        orf_annotator = EstscanOrfAnnotator(estscan_matrix)
        seq_records = orf_annotator(seq_records)
        orf1 = seq_records[0].features[0]
        orf2 = seq_records[1].features[0]
        assert orf1.strand == 1
        assert orf1.location.start.position == 0
        assert orf1.location.end.position == 541
        assert orf2.strand == -1
        assert orf2.location.start.position == 0
        assert orf2.location.end.position == 541
        assert not seq_records[2].features

if __name__ == "__main__":
    # import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
