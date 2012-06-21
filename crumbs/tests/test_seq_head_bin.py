'''
Created on 2012 eka 21

@author: peio
'''
import unittest
import os
from tempfile import NamedTemporaryFile
from subprocess import check_output

from crumbs.tests.utils import BIN_DIR

# pylint: disable=R0201
# pylint: disable=R0904


class SeqHead(unittest.TestCase):
    'It tests the seq head binary'

    def test_seq_head(self):
        'It tests the seq head'
        head_bin = os.path.join(BIN_DIR, 'seq_head')
        assert check_output([head_bin]).startswith('usage')

        #get one seq
        fasta_fhand = NamedTemporaryFile()
        fasta_fhand.write('>seq\nACTA\n>seq2\nACTA\n>seq3\nACTA\n')
        fasta_fhand.flush()

        result = check_output([head_bin, '-n', '1', fasta_fhand.name])
        assert result == '>seq\nACTA\n'

if __name__ == '__main__':
    #import sys;sys.argv = ['', 'SffExtractTest.test_items_in_gff']
    unittest.main()
