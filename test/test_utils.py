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
from subprocess import Popen, PIPE, check_output, CalledProcessError
from tempfile import NamedTemporaryFile

from crumbs.utils.file_utils import (TemporaryDir, rel_symlink,
                                     wrap_in_buffered_reader)
from crumbs.utils.bin_utils import check_process_finishes, popen, BIN_DIR

from crumbs.exceptions import ExternalBinaryError, MissingBinaryError
from crumbs.utils.test_utils import TEST_DATA_DIR
from crumbs.utils.tags import ERROR_ENVIRON_VARIABLE
from crumbs.seqio import guess_seq_type


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
        fhand3 = wrap_in_buffered_reader(fhand2, force_wrap=True)
        assert fhand3.peek(10) == 'hola'

    @staticmethod
    def test_guess_seq_type():
        'It checks that we can guess the type of the seqs in a file'
        fpath = os.path.join(TEST_DATA_DIR, 'arabidopsis_genes')
        guess_seq_type(open(fpath))


class ErrorHandlingTest(unittest.TestCase):
    'It tests the handling of the unexpected errors'
    def test_error_handling(self):
        'We handle the errors'
        cat_bin = os.path.join(BIN_DIR, 'cat_seqs')

        os.environ[ERROR_ENVIRON_VARIABLE] = 'Fail'
        stderr = NamedTemporaryFile(suffix='.error')
        try:
            check_output([cat_bin, '-h'], stderr=stderr)
            self.fail('Error expected')
        except CalledProcessError:
            pass
        finally:
            del os.environ[ERROR_ENVIRON_VARIABLE]
            stderr.close()

if __name__ == '__main__':
    #import sys;sys.argv = ['', 'ErrorHandlingTest.test_error_handling']
    unittest.main()
