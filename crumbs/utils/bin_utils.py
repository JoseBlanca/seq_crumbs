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


from crumbs.third_party import cgitb
from crumbs.exceptions import (UnknownFormatError, FileNotFoundError,
                               WrongFormatError, TooManyFiles,
                               MalformedFile, SampleSizeError,
                               ExternalBinaryError, MissingBinaryError,
                               IncompatibleFormatError,
                               UndecidedFastqVersionError, MaxNumReadsInMem,
                               PairDirectionError, InterleaveError)
from crumbs.utils.file_utils import (wrap_in_buffered_reader,
                                     uncompress_if_required, compress_fhand)
from crumbs.utils.seq_utils import guess_format
from crumbs.settings import get_setting
from crumbs.utils.tags import (OUTFILE, GUESS_FORMAT, BGZF, GZIP,
                               ERROR_ENVIRON_VARIABLE)
from crumbs import __version__ as version


BIN_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..',
                                       '..', 'bin'))


def main(funct):
    'The main function of a script'
    argv = sys.argv
    if '--error_log' in argv:
        error_fpath_index = argv.index('--error_log') + 1
        error_fpath = argv[error_fpath_index]
    else:
        binary = os.path.split(sys.argv[0])[-1]
        error_fpath = binary + '.error'

    stderr = sys.stderr
    try:
        # This code is required to test the error handling
        fail = os.environ.get(ERROR_ENVIRON_VARIABLE, None)
        if fail:
            raise RuntimeError('Generating a test error')
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
    except InterleaveError, error:
        stderr.write(str(error) + '\n')
        return 15
    except KeyboardInterrupt, error:
        stderr.write('Program stopped by user request\n')
        return 16
    except Exception as error:
        msg = 'An unexpected error happened.\n'
        msg += 'The seq_crumbs developers would appreciate your feedback.\n'
        try:
            fail = os.environ.get(ERROR_ENVIRON_VARIABLE, None)
            if fail:
                # error handling debugging
                fhand = sys.stderr
                error_fpath = None
            else:
                fhand = open(error_fpath, 'a')
        except IOError:
            # the fpath for the error is not writable
            fhand = sys.stderr
            error_fpath = None

        if error_fpath:
            msg += 'Please send them the error log'
            msg += ': ' + error_fpath + '\n\n'
        msg += str(error)
        stderr.write(msg)

        if error_fpath:
            hook = cgitb.Hook(display=0, format='text', logfpath=error_fpath)
            hook.handle()

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
    if get_setting('USE_EXTERNAL_BIN_PREFIX'):
        binary_name = get_setting('EXTERNAL_BIN_PREFIX') + binary_name

    if not get_setting('ADD_PATH_TO_EXT_BIN'):
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

    parser.add_argument('input', default=sys.stdin, nargs='*',
                        help='Sequence input files to process (default STDIN)',
                        type=argparse.FileType('rt'))

    parser.add_argument('-t', '--in_format', default=GUESS_FORMAT,
                       help='Format of the input files (default: %(default)s)')

    parser.add_argument('-o', '--outfile', default=sys.stdout, dest=OUTFILE,
                        help='Sequence output file (default: STDOUT)',
                        type=argparse.FileType('wt'))

    parser.add_argument('--version', action='version',
                        version=build_version_msg())

    group = parser.add_mutually_exclusive_group()
    group.add_argument('-z ', '--gzip', action='store_true',
                       help='Compress the output in gzip format')
    group.add_argument('-Z ', '--bgzf', action='store_true',
                       help='Compress the output in bgzf format')
    return parser


def build_version_msg():
    'It creates a message with the version.'
    bin_name = os.path.split(sys.argv[0])[-1]
    version_msg = bin_name + ' from seq_crumbs version: ' + version
    return version_msg


def create_basic_argparse(**kwargs):
    'It returns a cmd parser with inputs, output and format'
    parser = create_io_argparse(**kwargs)
    parser = argparse.ArgumentParser(parents=[parser], add_help=False)
    parser.add_argument('-f', '--out_format', dest='out_format',
                        help='output file format (default: same as input)',
                        choices=get_setting('SUPPORTED_OUTPUT_FORMATS'))
    return parser


def create_basic_parallel_argparse(**kwargs):
    'It returns a cmd parser with inputs, output, format, num_processes'
    parser = create_basic_argparse(**kwargs)
    parser = argparse.ArgumentParser(parents=[parser], add_help=False)
    parser.add_argument('-p', '--processes', dest='processes', type=int,
                        help='Num. of processes to use (default: %(default)s)',
                        default=1)
    return parser


def create_filter_argparse(add_reverse=True, **kwargs):
    'It returns a cmd parser for the filter executables'
    parser = create_basic_parallel_argparse(**kwargs)
    parser = argparse.ArgumentParser(parents=[parser], add_help=False)
    if add_reverse:
        parser.add_argument('-r', '--reverse', action='store_true',
                            help='Reverses the filtering')
    parser.add_argument('-e', '--filtered_file',
                        help='Filtered out sequences output file',
                        type=argparse.FileType('wt'))
    return parser


def create_trimmer_argparse(**kwargs):
    'It returns a cmd parser for the filter executables'
    parser = create_basic_parallel_argparse(**kwargs)
    parser = argparse.ArgumentParser(parents=[parser], add_help=False)
    parser.add_argument('-m', '--mask', dest='mask', action='store_true',
                        help='Do not trim, only mask by lowering the case')
    return parser


def get_requested_compression(parsed_args):
    'It looks in the selected options and return the selected compression kind'
    comp_kind = None
    bgzf = getattr(parsed_args, 'bgzf', False)
    gzip = getattr(parsed_args, 'gzip', False)
    if bgzf:
        comp_kind = BGZF
    elif gzip:
        comp_kind = GZIP
    return comp_kind


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

    # We have to add the one_line to the fastq files in order to get the
    # speed improvements of the seqitems
    in_format = parsed_args.in_format
    if 'fastq' in in_format:
        guessed_in_format = guess_format(wrapped_fhands[0])
        if '-one_line' in guessed_in_format:
            in_format += '-one_line'
    else:
        guessed_in_format = None

    out_fhand = getattr(parsed_args, OUTFILE)

    comp_kind = get_requested_compression(parsed_args)
    if isinstance(out_fhand, list):
        new_out_fhands = []
        for out_f in out_fhand:
            try:
                out_f = compress_fhand(out_f, compression_kind=comp_kind)
            except RuntimeError, error:
                parser.error(error)

            new_out_fhands.append(out_f)
        out_fhand = new_out_fhands
    else:
        try:
            out_fhand = compress_fhand(out_fhand, compression_kind=comp_kind)
        except RuntimeError, error:
            parser.error(error)

    out_format = parsed_args.out_format
    # The default format is the same as the first file
    if not out_format:
        if in_format == GUESS_FORMAT:
            if not guessed_in_format:
                guessed_in_format = guess_format(wrapped_fhands[0])
            out_format = guessed_in_format
        else:
            out_format = in_format
    # The original fhands should be stored, because otherwise they would be
    # closed
    args = {'out_fhand': out_fhand, 'in_fhands': wrapped_fhands,
            'out_format': out_format, 'original_in_fhands': in_fhands,
            'in_format': in_format}
    return args, parsed_args


def parse_basic_parallel_args(parser):
    'It parses the command line and it returns a dict with the arguments.'
    args, parsed_args = parse_basic_args(parser)
    args['processes'] = parsed_args.processes
    return args, parsed_args


def parse_filter_args(parser, add_reverse=True):
    'It parses the command line and it returns a dict with the arguments.'
    args, parsed_args = parse_basic_parallel_args(parser)
    if add_reverse:
        args['reverse'] = parsed_args.reverse
    args['filtered_fhand'] = parsed_args.filtered_file
    return args, parsed_args


def parse_trimmer_args(parser):
    'It parses the command line and it returns a dict with the arguments.'
    args, parsed_args = parse_basic_parallel_args(parser)
    args['mask'] = parsed_args.mask
    return args, parsed_args
