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

import tempfile
import shutil
import io
from gzip import GzipFile
import os.path

from Bio.bgzf import BgzfWriter

from crumbs.utils.tags import BGZF, GZIP, BZIP2
try:
    from crumbs.utils import BZ2File
except ImportError:
    pass
from crumbs.exceptions import OptionalRequirementError


def wrap_in_buffered_reader(fhand, force_wrap=False):
    '''It wraps the given file in a peekable BufferedReader.

    If the file is seekable it doesn't do anything.
    '''
    if not force_wrap and fhand_is_seekable(fhand):
        return fhand
    else:
        fhand = io.open(fhand.fileno(), mode='rb')  # with text there's no peek

    return fhand


def fhand_is_seekable(fhand):
    'It returns True if the fhand is seekable'
    try:
        try:
            # The stdin stream in some instances has seek, has no seekable
            fhand.tell()
        except IOError:
            return False
        try:
            if fhand.seekable():
                return True
            else:
                return False
        except AttributeError:
            return True
    except AttributeError:
        return False


def peek_chunk_from_file(fhand, chunk_size):
    'It returns the beginning of a file without moving the pointer'
    if fhand_is_seekable(fhand):
        fhand.seek(0)
        chunk = fhand.read(chunk_size)
        fhand.seek(0)
    else:
        chunk = fhand.peek(chunk_size)
    return chunk


BZIP_ERROR = 'bz2file (Python < 3.3) or bz2 module is required to work with '
BZIP_ERROR += 'bzip2 files'


def uncompress_if_required(fhand):
    'It returns a uncompressed handle if required'
    magic = peek_chunk_from_file(fhand, 2)
    if magic == '\037\213':
        fhand = GzipFile(fileobj=fhand)
    elif magic == 'BZ':
        try:
            fhand = BZ2File(fhand)
        except NameError:
            raise OptionalRequirementError(BZIP_ERROR)
    return fhand


def compress_fhand(fhand, compression_kind=None):
    'Compresses the file if required'
    if compression_kind == BGZF:
        if fhand_is_seekable(fhand):
            fhand = BgzfWriter(fileobj=fhand)
        else:
            raise RuntimeError('bgzf is only available for seekable files')
    elif compression_kind == GZIP:
        fhand = GzipFile(fileobj=fhand)
    elif compression_kind == BZIP2:
        mode = 'w' if 'w' in fhand.mode else 'r'
        try:
            fhand = BZ2File(fhand, mode=mode)
        except NameError:
            raise OptionalRequirementError(BZIP_ERROR)
    return fhand


class TemporaryDir(object):
    '''This class creates temporary directories '''
    def __init__(self, suffix='', prefix='', directory=None):
        '''It initiates the class.'''
        self._name = tempfile.mkdtemp(suffix=suffix, prefix=prefix,
                                      dir=directory)

    @property
    def name(self):
        'Returns path to the dict'
        return self._name

    def close(self):
        '''It removes the temp dir'''
        if os.path.exists(self._name):
            shutil.rmtree(self._name)


def _common_base(path1, path2):
    'it return the common and uncommon part of both strings'
    common = []
    path1 = path1.split(os.sep)
    path2 = path2.split(os.sep)

    for index, dir_ in enumerate(path1):
        try:
            if dir_ == path2[index]:
                common.append(dir_)
            else:
                break
        except IndexError:
            break

    uncommon1 = path1[len(common):]
    uncommon2 = path2[len(common):]
    return uncommon1, uncommon2


def _rel_path(path1, path2):
    'it return relative paths from path2 to path1'
    uncommon1, uncommon2 = _common_base(path1, path2)

    path = ['..'] * (len(uncommon2) - 1)
    path.extend(uncommon1)
    return os.sep.join(path)


def rel_symlink(path1, path2):
    'It makes the relative to path1 from within dirname(path2)'
    # The working directory can not be used because it is not thread safe
    # working directory is a global environment variable shared by all threads
    path1 = os.path.abspath(path1)
    path2 = os.path.abspath(path2)
    rel_path1 = _rel_path(path1, path2)

    fname2 = os.path.split(path2)[-1]

    temp_dir = tempfile.mkdtemp()
    try:
        temp_link_path = os.path.join(temp_dir, fname2)

        os.symlink(rel_path1, temp_link_path)
        shutil.move(temp_link_path, path2)
    finally:
        os.rmdir(temp_dir)
