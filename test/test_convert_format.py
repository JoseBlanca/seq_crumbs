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
import os.path
from tempfile import NamedTemporaryFile
from subprocess import check_output, CalledProcessError

from crumbs.utils.bin_utils import BIN_DIR

# pylint: disable=R0201
# pylint: disable=R0904

FASTA = ">seq1\natctagtc\n>seq2\natctagtc\n>seq3\natctagtc\n"
QUAL = ">seq1\n30 30 30 30 30 30 30 30\n>seq2\n30 30 30 30 30 30 30 30\n"
QUAL += ">seq3\n30 30 30 30 30 30 30 30\n"
FASTQ = '@seq1\natcgt\n+\n?????\n@seq2\natcgt\n+\n?????\n@seq3\natcgt\n+\n'
FASTQ += '?????\n'


def _make_fhand(content=None):
    'It makes temporary fhands'
    if content is None:
        content = ''
    fhand = NamedTemporaryFile()
    fhand.write(content)
    fhand.flush()
    return fhand


class FastaQualToFastqTest(unittest.TestCase):
    'It test the convert_format binary'

    def test_fastaqual_to_fasta(self):
        'It test fastaqual to fastq binary'
        fasta_fhand = _make_fhand(FASTA)
        qual_fhand = _make_fhand(QUAL)

        seqio_bin = os.path.join(BIN_DIR, 'fastaqual_to_fastq')
        assert 'usage' in check_output([seqio_bin, '-h'])

        out_fhand = NamedTemporaryFile()
        check_output([seqio_bin, '-o', out_fhand.name,
                      fasta_fhand.name, qual_fhand.name])
        assert "@seq1\natctagtc\n+" in  open(out_fhand.name).read()


class SeqioBinTest(unittest.TestCase):
    'It test the convert_format binary'

    def test_seqio_bin(self):
        'It test the seqio binary'
        seqio_bin = os.path.join(BIN_DIR, 'convert_format')
        assert 'usage' in check_output([seqio_bin, '-h'])

        # get one seq
        fasta_fhand = _make_fhand(FASTA)
        qual_fhand = _make_fhand(QUAL)
        fastq_fhand = _make_fhand(FASTQ)

        # fasta to fastq should fail
        out_fhand = NamedTemporaryFile()
        stderr = NamedTemporaryFile()
        try:
            print check_output([seqio_bin, '-o', out_fhand.name,
                          fasta_fhand.name, '-f', 'fastq' ], stderr=stderr)
            self.fail('Error expected')
        except CalledProcessError:
            assert 'No qualities available' in open(stderr.name).read()

        # bad_format_fastq
        bad_fastq_fhand = _make_fhand(FASTQ + 'aklsjhdas')
        fasta_out_fhand = NamedTemporaryFile()
        qual_out_fhand = NamedTemporaryFile()
        stderr = NamedTemporaryFile()
        try:
            print check_output([seqio_bin, '-o', fasta_out_fhand.name,
                                '-f', 'fasta', bad_fastq_fhand.name],
                               stderr=stderr)
            self.fail('error expected')
        except CalledProcessError:
            assert 'Lengths of sequence and qualit' in open(stderr.name).read()

        # fastq to fastq
        fastq_fhand = _make_fhand(FASTQ)
        fastq_fhand2 = _make_fhand(FASTQ)
        fastq_out_fhand = NamedTemporaryFile()
        stderr = NamedTemporaryFile()
        check_output([seqio_bin, '-o', fastq_out_fhand.name, '-f',
                      'fastq-illumina', fastq_fhand.name, fastq_fhand2.name],
                     stderr=stderr)
        out_fastq = open(fastq_out_fhand.name).read()
        assert '+\n^^^^^\n@seq1\natcgt' in  out_fastq

        # test stdin
        fasta_out_fhand = NamedTemporaryFile()
        check_output([seqio_bin, '-o', fasta_out_fhand.name, '-f', 'fasta'],
                      stdin=open(fastq_fhand.name))
        assert ">seq1\natcgt" in  open(fasta_out_fhand.name).read()

    def test_version(self):
        'It can return its version number'
        guess_bin = os.path.join(BIN_DIR, 'convert_format')
        stderr = NamedTemporaryFile()
        check_output([guess_bin, '--version'], stderr=stderr)
        assert 'from seq_crumbs version:' in open(stderr.name).read()


if __name__ == '__main__':
    # import sys;sys.argv = ['', 'SeqioBinTest.test_seqio_bin']
    unittest.main()
