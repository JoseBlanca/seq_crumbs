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
import os
from subprocess import Popen, PIPE
from tempfile import NamedTemporaryFile
from cStringIO import StringIO

from crumbs.utils import (TemporaryDir, rel_symlink, wrap_in_buffered_reader,
                          check_process_finishes, popen, guess_format,
                          fhand_is_seekable)
from crumbs.exceptions import (ExternalBinaryError, MissingBinaryError,
                               UnknownFormatError)


# pylint: disable=R0201
# pylint: disable=R0904

class UtilsTest(unittest.TestCase):
    'It tests the miscelaneous utilities.'

    def test_rel_symlink(self):
        'It tests various cases of rel symlinks'
        tempdir = TemporaryDir()
        try:
            hola = os.path.join(tempdir.name, 'hola')
            os.mkdir(hola)
            caracola = os.path.join(tempdir.name, 'caracola')
            rel_symlink(hola, caracola)
            assert os.path.exists(caracola)

            fname = os.path.join(hola, 'fname')
            open(fname, 'w')
            caracola2 = os.path.join(tempdir.name, 'caracola2')
            rel_symlink(fname, caracola2)
            assert os.path.exists(caracola2)

            path2 = os.path.join(tempdir.name, 'dir1', 'dir2')
            os.makedirs(path2)
            caracola3 = os.path.join(path2, 'caracola3')
            rel_symlink(hola, caracola3)
            assert os.path.exists(caracola3)
        finally:
            tempdir.close()

    def test_check_process(self):
        'It checks that a process finishes OK'

        binary = 'ls'
        process = Popen([binary, '/directorio_que_no_existe'], stderr=PIPE)
        try:
            check_process_finishes(process, binary)
            self.fail('ExternalBinaryErrorExpected')
        except ExternalBinaryError, error:
            assert binary in str(error)

        stderr = NamedTemporaryFile()
        process = Popen([binary, '/directorio_que_no_existe'], stderr=stderr)
        try:
            check_process_finishes(process, binary, stderr=stderr)
        except ExternalBinaryError, error:
            assert binary in str(error)

        cmd = [binary]
        process = Popen(cmd, stdout=PIPE)
        check_process_finishes(process, binary)

    def test_popen(self):
        'It checks that we can create a process'

        try:
            popen(['bad_binary'])
            self.fail()
        except MissingBinaryError:
            pass

        popen(['ls'], stdout=PIPE)

    def test_io_open(self):
        'It checks the file wrapping in ReaderBuffers'

        # a standard file
        fhand = NamedTemporaryFile()
        fhand.write('hola')
        fhand.flush()
        fhand2 = open(fhand.name)
        fhand3 = wrap_in_buffered_reader(fhand2)
        assert fhand3.peek(10) == 'hola'


class GuessFormatTest(unittest.TestCase):
    'It tests the function that guess the sequence format'

    def test_fasta(self):
        'It guess fasta formats'
        fhand = StringIO('>seq\nACTC\n')
        assert guess_format(fhand) == 'fasta'

        # qual
        fhand = StringIO('>seq\n10 20\n')
        assert guess_format(fhand) == 'qual'

        # qual
        qual = ">seq1\n30 30 30 30 30 30 30 30\n>seq2\n30 30 30 30 30 30 30"
        qual += " 30\n>seq3\n30 30 30 30 30 30 30 30\n"

        fhand = StringIO(qual)
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

    def test_is_seekable(self):
        'It tests wether the fhands are seekable or not'

        # StringIO
        fhand = StringIO('hola')
        assert fhand_is_seekable(fhand)

        # standard file
        fhand = NamedTemporaryFile()
        fhand.seek(0)
        assert fhand_is_seekable(fhand)

        # a wrapped BufferedReader
        fhand2 = wrap_in_buffered_reader(fhand)
        assert fhand_is_seekable(fhand2)

        class NonSeekable(object):
            'Just for testing'
            pass

        assert not fhand_is_seekable(NonSeekable())

        class NonSeekable2(object):
            def seek(self):
                pass

            def seekable(self):
                return False

        assert not fhand_is_seekable(NonSeekable2())

if __name__ == '__main__':
    #import sys;sys.argv = ['', 'SffExtractTest.test_items_in_gff']
    unittest.main()
