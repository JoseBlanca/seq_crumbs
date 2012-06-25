'''
Created on 21/06/2012

@author: jose
'''

import unittest
from StringIO import StringIO

from crumbs.seqio import guess_format
from crumbs.utils.exceptions import UnknownFormatError

# pylint: disable=R0201
# pylint: disable=R0904


class GuessFormatTest(unittest.TestCase):
    'It tests the function that guess the sequence format'

    def test_fasta(self):
        'It guess fasta formats'
        fhand = StringIO('>seq\nACTC\n')
        assert guess_format(fhand) == 'fasta'

        # qual
        fhand = StringIO('>seq\n10 20\n')
        assert guess_format(fhand) == 'qual'

    def test_unkown(self):
        'It tests unkown formats'
        fhand = StringIO('xseq\nACTC\n')
        try:
            guess_format(fhand)
            self.fail('UnknownFormatError expected')
        except UnknownFormatError:
            pass

    def test_empty_file(self):
        'It guesses the format of an empty file'
        fhand = StringIO()
        try:
            guess_format(fhand)
            self.fail('UnknownFormatError expected')
        except UnknownFormatError:
            pass

    def test_fastq(self):
        'It guesses the format for the solexa and illumina fastq'

        txt = '@HWI-EAS209_0006_FC706VJ:5:58:5894:21141#ATCACG/1\n'
        txt += 'TTAATTGGTAAATAAATCTCCTAATAGCTTAGATNTTACCTTNNNNNNNNNNTAGTTTCT\n'
        txt += '+HWI-EAS209_0006_FC706VJ:5:58:5894:21141#ATCACG/1\n'
        txt += 'efcfffffcfeefffcffffffddf`feed]`]_Ba_^__[YBBBBBBBBBBRTT\]][]\n'
        fhand = StringIO(txt)
        assert guess_format(fhand) == 'fastq-illumina'

        fhand = StringIO('@HWI-EAS209\n@')
        try:
            assert guess_format(fhand) == 'fasta'
            self.fail('UnknownFormatError expected')
        except UnknownFormatError:
            pass

if __name__ == '__main__':
    #import sys;sys.argv = ['', 'SffExtractTest.test_items_in_gff']
    unittest.main()
