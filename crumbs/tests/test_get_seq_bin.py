'''
Created on 2012 eka 22

@author: peio
'''
import unittest
import os
from tempfile import NamedTemporaryFile
from subprocess import check_output, CalledProcessError, check_call

from crumbs.tests.utils import BIN_DIR


# pylint: disable=R0201
# pylint: disable=R0904


class GetHead(unittest.TestCase):
    'It tests the seq head binary'

    def test_seq_head(self):
        'It tests the seq head'
        get_seq_bin = os.path.join(BIN_DIR, 'get_seqs')
        assert check_output([get_seq_bin]).startswith('usage')

        #get one seq
        fasta_fhand = NamedTemporaryFile()
        fasta_fhand.write('>seq\nACTA\n>seq2\nACTA\n>seq3\nACTA\n')
        fasta_fhand.flush()

        #head
        result = check_output([get_seq_bin, '-t', 'head', '-n', '1',
                               fasta_fhand.name])
        assert result == '>seq\nACTA\n'

        #random
        result = check_output([get_seq_bin, '-t', 'random', '-n', '1',
                               fasta_fhand.name])
        assert 'ACTA' in result

        #bad_type
        stderr = NamedTemporaryFile()
        try:
            check_call([get_seq_bin, '-t', 'am', '-n', '1',
                          fasta_fhand.name], stderr=stderr)
            self.fail('Error expectd')
        except CalledProcessError:
            assert 'usage' in open(stderr.name).read()


if __name__ == '__main__':
    #import sys;sys.argv = ['', 'SffExtractTest.test_items_in_gff']
    unittest.main()
