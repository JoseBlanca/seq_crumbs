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

import sys
import os
from subprocess import Popen
import platform
import argparse
from gzip import GzipFile

from Bio.bgzf import BgzfWriter

from crumbs.third_party import cgitb
from crumbs.exceptions import (UnknownFormatError, FileNotFoundError,
                               WrongFormatError, TooManyFiles,
                               MalformedFile, SampleSizeError,
                               ExternalBinaryError, MissingBinaryError,
                               IncompatibleFormatError,
                               UndecidedFastqVersionError, MaxNumReadsInMem,
                               PairDirectionError)
from crumbs.utils.file_utils import (wrap_in_buffered_reader,
                                     uncompress_if_required, fhand_is_seekable)
from crumbs.utils.seq_utils import guess_format
from crumbs.settings import (SUPPORTED_OUTPUT_FORMATS, USE_EXTERNAL_BIN_PREFIX,
                              EXTERNAL_BIN_PREFIX, ADD_PATH_TO_EXT_BIN)
from crumbs.utils.tags import OUTFILE, GUESS_FORMAT


BIN_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..',
                                       '..', 'bin'))


def main(funct):
    'The main function of a script'
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
        return 3
    except UnknownFormatError, error:
        stderr.write(str(error) + '\n')
        return 4
    except WrongFormatError, error:
        stderr.write(str(error) + '\n')
        return 5
    except TooManyFiles, error:
        stderr.write(str(error) + '\n')
        return 6
    except MalformedFile, error:
        stderr.write(str(error) + '\n')
        return 7
    except SampleSizeError, error:
        stderr.write(str(error) + '\n')
        return 8
    except ExternalBinaryError, error:
        stderr.write(str(error) + '\n')
        return 9
    except MissingBinaryError, error:
        stderr.write(str(error) + '\n')
        return 10
    except IncompatibleFormatError, error:
        stderr.write(str(error) + '\n')
        return 11
    except UndecidedFastqVersionError, error:
        stderr.write(str(error) + '\n')
        return 12
    except MaxNumReadsInMem, error:
        stderr.write(str(error) + '\n')
        return 13
    except PairDirectionError, error:
        stderr.write(str(error) + '\n')
        return 14
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
    if USE_EXTERNAL_BIN_PREFIX:
        binary_name = EXTERNAL_BIN_PREFIX + binary_name

    if not ADD_PATH_TO_EXT_BIN:
        # I have to check if the bynary is on my current directory.
        # If it is there use it, else assumes that it is on the path
        if os.path.exists(os.path.join(os.getcwd(), binary_name)):
            return os.path.join(os.getcwd(), binary_name)
        return binary_name

    system = platform.system().lower()
    if system == 'windows':
        binary_name += '.exe'
    arch = platform.architecture()[0]

    join = os.path.join

    module_path = os.path.split(__file__)[0]
    third_party_path = join(module_path, '..', 'third_party', 'bin')
    third_party_path = os.path.abspath(third_party_path)

    if not os.path.exists(third_party_path):
        msg = 'Third party bin directory not found, please fix me.'
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


def create_io_argparse(**kwargs):
    'It returns a parser with several inputs and one output'
    parser = argparse.ArgumentParser(**kwargs)

    parser.add_argument('input', help='Sequence input files to process',
                        default=sys.stdin, nargs='*',
                        type=argparse.FileType('rt'))

    parser.add_argument('-t', '--in_format', help='Format of the input files',
                        default=GUESS_FORMAT)

    parser.add_argument('-o', '--outfile', default=sys.stdout, dest=OUTFILE,
                        help='Sequence output file to process',
                        type=argparse.FileType('wt'))
    group = parser.add_mutually_exclusive_group()
    group.add_argument('-z ', '--gzip', action='store_true',
                       help='Compress the output in gzip format')
    group.add_argument('-Z ', '--bgzf', action='store_true',
                       help='Compress the output in bgzf format')
    return parser


def create_basic_argparse(**kwargs):
    'It returns a cmd parser with inputs, output and format'
    parser = create_io_argparse(**kwargs)
    parser = argparse.ArgumentParser(parents=[parser], add_help=False)
    parser.add_argument('-f', '--out_format', dest='out_format',
                        help='output file format',
                        choices=SUPPORTED_OUTPUT_FORMATS)
    return parser


def create_basic_process_argparse(**kwargs):
    'It returns a cmd parser with inputs, output, format, num_processes'
    parser = create_basic_argparse(**kwargs)
    parser = argparse.ArgumentParser(parents=[parser], add_help=False)
    parser.add_argument('-p', '--processes', dest='processes', type=int,
                        help='Num. of processes to use', default=1)
    return parser


def parse_basic_args(parser):
    'It parses the command line and it returns a dict with the arguments.'
    parsed_args = parser.parse_args()
    # we have to wrap the file in a BufferedReader to allow peeking into stdin
    wrapped_fhands = []
    # if input is stdin it will be a fhand not a list of fhands.
    # we have to convert to a list
    in_fhands = parsed_args.input
    if not isinstance(in_fhands, list):
        in_fhands = [in_fhands]
    for fhand in in_fhands:
        fhand = wrap_in_buffered_reader(fhand)
        fhand = uncompress_if_required(fhand)
        wrapped_fhands.append(fhand)

    in_format = parsed_args.in_format

    out_fhand = getattr(parsed_args, OUTFILE)

    if parsed_args.bgzf:
        if fhand_is_seekable(out_fhand):
            out_fhand = BgzfWriter(fileobj=out_fhand)
        else:
            parser.error('bgzf is incompatible with STDOUT.')
    elif parsed_args.gzip:
        out_fhand = GzipFile(fileobj=out_fhand)

    out_format = parsed_args.out_format
    # The default format is the same as the first file
    if not out_format:
        if in_format == GUESS_FORMAT:
            out_format = guess_format(wrapped_fhands[0])
        else:
            out_format = in_format
    # The original fhands should be stored, because otherwise they would be
    # closed
    args = {'out_fhand': out_fhand, 'in_fhands': wrapped_fhands,
            'out_format': out_format, 'original_in_fhands': in_fhands,
            'in_format': in_format}
    return args, parsed_args


def parse_basic_process_args(parser):
    'It parses the command line and it returns a dict with the arguments.'
    args, parsed_args = parse_basic_args(parser)
    args['processes'] = parsed_args.processes
    return args, parsed_args
