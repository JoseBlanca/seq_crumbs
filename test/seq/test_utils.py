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
import sys
import hashlib
from subprocess import Popen, PIPE, check_output, CalledProcessError
from tempfile import NamedTemporaryFile
from cStringIO import StringIO

from crumbs.utils.file_utils import (TemporaryDir, rel_symlink,
                                     wrap_in_buffered_reader)
from crumbs.utils.bin_utils import check_process_finishes, popen, BIN_DIR

from crumbs.exceptions import ExternalBinaryError, MissingBinaryError
from crumbs.utils.test_utils import TEST_DATA_DIR
from crumbs.utils.tags import ERROR_ENVIRON_VARIABLE
from crumbs.seq.seqio import guess_seq_type
from crumbs.settings import get_setting
from crumbs.seq.utils.file_formats import (get_format, FILEFORMAT_INVENTORY,
                                           set_format)
from crumbs.utils.sqlite_utils import SqliteCache



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

    def test_get_format_fhand(self):
        "It checks the get/set format functions"
        #file fhand
        fhand = NamedTemporaryFile()
        fhand.write('>seq\natgctacgacta\n')
        fhand.flush()
        name = fhand.name
        id_ = id(fhand)

        file_format = get_format(fhand)
        assert FILEFORMAT_INVENTORY[(id_, name)] == file_format
        num_keys = len(FILEFORMAT_INVENTORY)

        file_format = get_format(fhand)
        assert FILEFORMAT_INVENTORY[(id_, name)] == file_format
        assert len(FILEFORMAT_INVENTORY) == num_keys

        fhand = NamedTemporaryFile()
        set_format(fhand, 'fasta')

        assert 'fasta' == get_format(fhand)

    def test_get_format_stringio(self):
        "It checks the get/set format functions"
        #stiongIO
        stringIO_fhand = StringIO('>seq\natgctacgacta\n')

        striongIOhash = hashlib.sha224(stringIO_fhand.getvalue()[:100]).hexdigest()
        id_ = id(stringIO_fhand)

        file_format = get_format(stringIO_fhand)
        assert FILEFORMAT_INVENTORY[(id_, striongIOhash)] == file_format


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


class SettingsTest(unittest.TestCase):
    'It tests the get_settings function'
    def test_get_settings(self):
        'We get the settings'
        kmer_size = get_setting('DEFAULT_KMER_SIZE')
        assert kmer_size

    def test_environment_settings(self):
        'It sets settings with the environment variables'
        code = 'from crumbs.settings import get_setting\n'
        code += 'print get_setting("DEFAULT_KMER_SIZE")\n'

        test_script = NamedTemporaryFile()
        test_script.write(code)
        test_script.flush()

        environ = {key: val for key, val in os.environ.items()}
        environ['SEQ_CRUMBS_DEFAULT_KMER_SIZE'] = '1'
        cmd = [sys.executable, test_script.name]
        out = check_output(cmd, env=environ)
        assert out == '1\n'


class SqliteTest(unittest.TestCase):
    def test_init(self):
        fhand = NamedTemporaryFile()
        with SqliteCache(fhand.name) as sqlitecache:
            assert sqlitecache['seq1'] is None
            sqlitecache['seq1'] = 23
            assert sqlitecache['seq1'] == 23
            assert 'seq1' in sqlitecache
            assert 'seq2' not in sqlitecache

            # reopen the same cache file
            sqlitecache2 = SqliteCache(fhand.name)
            assert sqlitecache['seq1'] == 23
            sqlitecache2.close()

        fhand = NamedTemporaryFile()
        with SqliteCache(fhand.name) as sqlitecache:
            assert sqlitecache['seq1'] is None
            sqlitecache['seq1'] = {'1': 23, '3': 2}
            assert sqlitecache['seq1'] == {'1': 23, '3': 2}
            sqlitecache['seq1'] = {'1': 22, '3': 2}
            result = str(sqlitecache)
            assert result == "key\t1\t3\nseq1\t22\t2\n" in result
            sqlitecache['seq2'] = {'1': 22, '3': 2}
            dump = list(sqlitecache.dump())
            assert dump == [(u'seq1', {'1': 22, '3': 2}),
                            (u'seq2', {'1': 22, '3': 2})]


if __name__ == '__main__':
    #import sys;sys.argv = ['', 'UtilsTest.test_get_format_stringio']
    unittest.main()
