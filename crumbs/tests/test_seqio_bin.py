'''
Created on 2012 eka 25

@author: peio
'''
import unittest
import os
from subprocess import check_output, CalledProcessError
from tempfile import NamedTemporaryFile


from crumbs.tests.utils import BIN_DIR

FASTA = ">seq1\natctagtc\n>seq2\natctagtc\n>seq3\natctagtc\n"
QUAL = ">seq1\n30 30 30 30 30 30 30 30\n>seq2\n30 30 30 30 30 30 30 30\n"
QUAL += ">seq3\n30 30 30 30 30 30 30 30\n"
FASTQ = '@seq1\natcgt\n+\n?????\n@seq2\natcgt\n+\n?????\n@seq3\natcgt\n+\n?????\n'

class SeqioBinTest(unittest.TestCase):
    'It test the seqio binary'

    @staticmethod
    def _make_fhand(content=None):
        'It makes temporary fhands'
        if content is None:
            content = ''
        fhand = NamedTemporaryFile()
        fhand.write(content)
        fhand.flush()
        return fhand

    def test_seqio_bin(self):
        'It test the seqio binary'
        seqio_bin = os.path.join(BIN_DIR, 'seqio')
        assert check_output([seqio_bin]).startswith('usage')

        #get one se
        fasta_fhand = self._make_fhand(FASTA)
        qual_fhand = self._make_fhand(QUAL)
        fastq_fhand = self._make_fhand(FASTQ)

        # fasta-qual to fastq
        out_fhand = NamedTemporaryFile()
        check_output([seqio_bin, '-o', out_fhand.name, '-f', 'fastq',
                      fasta_fhand.name, qual_fhand.name])
        assert "@seq1\natctagtc\n+" in  open(out_fhand.name).read()

        #fastq to fast-qual
        fasta_out_fhand = NamedTemporaryFile()
        qual_out_fhand = NamedTemporaryFile()
        check_output([seqio_bin, '-o', fasta_out_fhand.name,
                      qual_out_fhand.name, '-f', 'fasta', fastq_fhand.name])
        assert ">seq1\natcgt" in  open(fasta_out_fhand.name).read()
        assert ">seq1\n30 30 30 30 30" in  open(qual_out_fhand.name).read()

        #bad_format_fasta
        bad_fasta_fhand = self._make_fhand(FASTA + 'asdsa')
        out_fhand = NamedTemporaryFile()
        stderr = NamedTemporaryFile()
        try:
            check_output([seqio_bin, '-o', out_fhand.name, '-f', 'fastq',
                         bad_fasta_fhand.name, qual_fhand.name], stderr=stderr)
            self.fail('error expected')
        except CalledProcessError:
            assert 'Sequence length and number' in open(stderr.name).read()

        #bad_format_fastq
        bad_fastq_fhand = self._make_fhand(FASTQ + 'aklsjhdas')
        fasta_out_fhand = NamedTemporaryFile()
        qual_out_fhand = NamedTemporaryFile()
        stderr = NamedTemporaryFile()
        try:
            check_output([seqio_bin, '-o', fasta_out_fhand.name,
                          qual_out_fhand.name, '-f', 'fasta',
                           bad_fastq_fhand.name], stderr=stderr)
            self.fail('error expected')
        except CalledProcessError:
            assert 'Lengths of sequence and qualit' in open(stderr.name).read()

        #malformed fastq to fastq
        bad_fastq_fhand = self._make_fhand(FASTQ + 'aklsjhdas')
        fastq_out_fhand = NamedTemporaryFile()
        stderr = NamedTemporaryFile()
        try:
            check_output([seqio_bin, '-o', fastq_out_fhand.name, '-f',
                          'fastq-illumina', bad_fastq_fhand.name],
                          stderr=stderr)
            self.fail('error expected')
        except CalledProcessError:
            assert 'Lengths of sequence and qualit' in open(stderr.name).read()





if __name__ == '__main__':
    #import sys;sys.argv = ['', 'SffExtractTest.test_items_in_gff']
    unittest.main()
