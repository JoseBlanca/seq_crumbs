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

'''
Created on 21/06/2012

@author: jose
'''

import sys
import os
import tempfile
import shutil
from subprocess import check_call, Popen
import platform

from crumbs.third_party import cgitb
from crumbs.exceptions import (UnknownFormatError, FileNotFoundError,
                                     WrongFormatError, TooManyFiles,
                                     MalformedFile, SampleSizeError,
                                     ExternalBinaryError, MissingBinaryError)

STDIN = 'stdin'
STDOUT = 'stdout'
INFILES = 'infiles'
OUTFILE = 'output'


def main(funct):
    'The main function of a script'
    if len(sys.argv) == 1:
        sys.argv = sys.argv + ['-h']

    argv = sys.argv
    if '--error_log' in argv:
        error_fpath_index = argv.index('--error_log') + 1
        error_fpath = argv[error_fpath_index]
    else:
        binary = sys.argv[0]
        error_fpath = binary + '.error'

    stderr = sys.stderr
    try:
        return(funct())
    except FileNotFoundError, error:
        stderr.write(str(error) + '\n')
        return 2
    except UnknownFormatError, error:
        stderr.write(str(error) + '\n')
        return 3
    except WrongFormatError, error:
        stderr.write(str(error) + '\n')
        return 4
    except TooManyFiles, error:
        stderr.write(str(error) + '\n')
        return 5
    except MalformedFile, error:
        stderr.write(str(error) + '\n')
        return 6
    except SampleSizeError, error:
        stderr.write(str(error) + '\n')
        return 7
    except ExternalBinaryError, error:
        stderr.write(str(error) + '\n')
        return 8
    except MissingBinaryError, error:
        stderr.write(str(error) + '\n')
        return 9
    except Exception as error:
        msg = 'An unexpected error happened.\n'
        msg += 'The seq crumbs developers would appreciate your feedback\n'
        msg += 'Please send them the error log: '
        msg += error_fpath + '\n\n'
        msg += str(error)
        stderr.write(msg)
        hook = cgitb.Hook(display=0, format='text', logfpath=error_fpath)
        hook.handle()
        fhand = open(error_fpath, 'a')
        fhand.write('\nThe command was:\n' + ' '.join(sys.argv) + '\n')
        fhand.close()
        raise


def get_inputs_from_args(parsed_args):
    'It returns the input fhand'
    in_fpaths = getattr(parsed_args, INFILES)
    if in_fpaths == STDIN:
        in_fhands = [sys.stdin]
    else:
        in_fhands = []
        for in_fpath in in_fpaths:
            if os.path.exists(in_fpath):
                in_fhand = open(in_fpath, 'rt')
            else:
                raise FileNotFoundError('A file was not found: ' + in_fpath)
            in_fhands.append(in_fhand)
    return in_fhands


def get_output_from_args(parsed_args):
    'It returns the out_fhand'
    out_fpath = getattr(parsed_args, OUTFILE)
    if out_fpath == STDOUT:
        out_fhand = sys.stdout
    else:
        out_fhand = open(out_fpath, 'w')
    return out_fhand


class TemporaryDir(object):
    '''This class creates temporary directories '''
    def __init__(self, suffix='', prefix='', directory=None):
        '''It initiates the class.'''
        self._name = tempfile.mkdtemp(suffix=suffix, prefix=prefix,
                                      dir=directory)

    def get_name(self):
        'Returns path to the dict'
        return self._name
    name = property(get_name)

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
    'It makes the relative symlink'
    path1 = os.path.abspath(path1)
    path2 = os.path.abspath(path2)
    fname2 = os.path.split(path2)[-1]
    rel_path1 = _rel_path(path1, path2)

    #we need a temp dir to be thread safe when creating the link
    #we cannot use chdir, that is not threadsafe
    temp_dir = tempfile.mkdtemp()
    try:
        temp_link_path = os.path.join(temp_dir, fname2)

        cmd = ['ln', '-s', '-T', rel_path1, temp_link_path]
        check_call(cmd)
        cmd2 = ['mv', temp_link_path, path2]
        check_call(cmd2)
    finally:
        os.rmdir(temp_dir)


def check_process_finishes(process, binary, stdout=None, stderr=None):
    'It checks that the given process finishes OK, otherwise raises an Error'

    stdout_msg, stderr_msg = process.communicate()

    returncode = process.returncode
    if returncode == 0:
        return
    elif returncode is None:
        msg = 'The process for {:s} is still running with PID'
        msg = msg.format(binary, process.PID)
        raise RuntimeError(msg)

    if stdout and not stdout_msg:
        stdout.flush()
        stdout_msg = open(stdout.name).read()
    if stderr and not stderr_msg:
        stderr.flush()
        stderr_msg = open(stderr.name).read()
    msg = '{:s} had a problem running\n'.format(binary)
    if stdout_msg:
        msg += 'stdout:\n{:s}\n'.format(stdout_msg)
    if stderr_msg:
        msg += 'stderr:\n{:s}\n'.format(stderr_msg)
    raise ExternalBinaryError(msg)


def popen(*args, **kwargs):
    'It spawns a process with subprocess.Popen'
    try:
        return Popen(*args, **kwargs)
    except OSError, error:
        msg = 'The binary "%s" is not in the path.' % args[0][0]
        if 'No such file' in str(error):
            raise MissingBinaryError(msg)


def get_binary_path(binary_name):
    '''It return the path to the proper binary. It looks on platform and
    architecture to decide it.

    Fails if there is not binary for that architecture
    '''
    system = platform.system().lower()
    if system == 'windows':
        binary_name += '.exe'
    arch = platform.architecture()[0]

    join = os.path.join

    module_path = os.path.split(__file__)[0]
    third_party_path = join(module_path, 'third_party', 'bin')
    if not os.path.exists(third_party_path):
        msg = 'Third party bin directory not found, please fixme.'
        raise MissingBinaryError(msg)

    binary_path = os.path.abspath(join(third_party_path, system, arch,
                                       binary_name))

    if os.path.exists(binary_path):
        return binary_path
    elif arch == '64bit':
        arch = '32bit'
        binary_path = os.path.abspath(join(third_party_path, system, arch,
                                           binary_name))
        if os.path.exists(binary_path):
            return binary_path

    # At this point there is not available binary for the working platform
    msg = '{} not available for this platform: {}'.format(binary_name, system)
    raise MissingBinaryError(msg)
