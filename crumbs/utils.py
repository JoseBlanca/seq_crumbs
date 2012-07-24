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
import argparse
import io
import cStringIO
from array import array
import itertools
from multiprocessing import Pool
from gzip import GzipFile

from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio.bgzf import BgzfWriter

from crumbs.third_party import cgitb
from crumbs.exceptions import (UnknownFormatError, FileNotFoundError,
                               WrongFormatError, TooManyFiles,
                               MalformedFile, SampleSizeError,
                               ExternalBinaryError, MissingBinaryError,
                               IncompatibleFormatError,
                               UndecidedFastqVersionError)
from crumbs.settings import (OUTFILE, SUPPORTED_OUTPUT_FORMATS, GUESS_FORMAT,
                             CHUNK_TO_GUESS_FASTQ_VERSION,
                             SEQS_TO_GUESS_FASTQ_VERSION,
                             LONGEST_EXPECTED_ILLUMINA_READ)


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
                        help='format of the output file',
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


def uncompress_if_required(fhand):
    'It returns a uncompressed handle if required'
    magic = _peek_chunk_from_file(fhand, 2)
    if magic == '\037\213':
        fhand = GzipFile(fileobj=fhand)
    return fhand


def parse_basic_process_args(parser):
    'It parses the command line and it returns a dict with the arguments.'
    args, parsed_args = parse_basic_args(parser)
    args['processes'] = parsed_args.processes
    return args, parsed_args


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


def _peek_chunk_from_file(fhand, chunk_size):
    'It returns the beginning of a file without moving the pointer'
    if fhand_is_seekable(fhand):
        fhand.seek(0)
        chunk = fhand.read(chunk_size)
        fhand.seek(0)
    else:
        chunk = fhand.peek(chunk_size)
    return chunk


def _get_some_qual_and_lengths(fhand, force_file_as_non_seek):
    'It returns the quality characters and the lengths'
    seqs_to_peek = SEQS_TO_GUESS_FASTQ_VERSION
    chunk_size = CHUNK_TO_GUESS_FASTQ_VERSION

    lengths = array('I')
    seqs_analyzed = 0
    if fhand_is_seekable(fhand) and not force_file_as_non_seek:
        fmt_fhand = fhand
    else:
        chunk = _peek_chunk_from_file(fhand, chunk_size)
        fmt_fhand = cStringIO.StringIO(chunk)

    try:
        for seq in FastqGeneralIterator(fmt_fhand):
            qual = [ord(char) for char in seq[2]]
            sanger_chars = [q for q in qual if q < 64]
            if sanger_chars:
                fhand.seek(0)
                return None, True     # no quals, no lengths, is_sanger
            lengths.append(len(qual))
            seqs_analyzed += 1
            if seqs_analyzed > seqs_to_peek:
                break
    except ValueError:
        raise UnknownFormatError('Malformed fastq')
    finally:
        fhand.seek(0)
    return lengths, None     # quals, lengths, don't know if it's sanger


def _guess_fastq_version(fhand, force_file_as_non_seek):
    '''It guesses the format of fastq files.

    It ignores the solexa fastq version.
    '''
    lengths, is_sanger = _get_some_qual_and_lengths(fhand,
                                                    force_file_as_non_seek)
    if is_sanger:
        return 'fastq'
    elif is_sanger is False:
        return 'fastq-illumina'
    n_long_seqs = [l for l in lengths if l > LONGEST_EXPECTED_ILLUMINA_READ]
    if n_long_seqs:
        msg = 'It was not possible to guess the format of '
        if hasattr(fhand, 'name'):
            msg += 'the file ' + fhand.name
        else:
            msg += 'a file '
        msg = '\n. The quality values could be Illumina, but there are '
        msg += 'sequences longer than %i bp.'
        msg %= LONGEST_EXPECTED_ILLUMINA_READ
        raise UndecidedFastqVersionError(msg)
    else:
        return 'fastq-illumina'


def guess_format(fhand):
    '''It guesses the format of the sequence file.

    It does ignore the solexa fastq version.
    '''
    return _guess_format(fhand, force_file_as_non_seek=False)


def _guess_format(fhand, force_file_as_non_seek):
    '''It guesses the format of the sequence file.

    This function is just for testing forcing the fhand as non-seekable.
    It does ignore the solexa fastq version.
    '''
    chunk_size = 1024
    chunk = _peek_chunk_from_file(fhand, chunk_size)
    if not chunk:
        raise UnknownFormatError('The file is empty')
    lines = chunk.splitlines()
    if chunk.startswith('>'):
        if lines[1].startswith('>'):
            raise UnknownFormatError('Malformed fasta')
        else:
            first_item = lines[1].strip().split()[0]
            if first_item.isdigit():
                return 'qual'
            else:
                return 'fasta'
    elif chunk.startswith('@'):
        return _guess_fastq_version(fhand, force_file_as_non_seek)
    elif chunk.startswith('LOCUS'):
        return 'genbank'
    elif chunk.startswith('ID'):
        return 'embl'
    raise UnknownFormatError('Sequence file of unknown format.')


def process_sequences(seq_packets, map_functions, processes=1,
                      keep_order=False):
    'It processes de sequences. It can do it in parallel'
    if processes > 1:
        pool = Pool(processes=processes)
        mapper = pool.imap if keep_order else pool.imap_unordered
    else:
        mapper = itertools.imap

    for map_function in map_functions:
        seq_packets = mapper(map_function, seq_packets)

    return seq_packets
