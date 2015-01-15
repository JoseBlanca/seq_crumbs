# Copyright 2013 Jose Blanca, Peio Ziarsolo, COMAV-Univ. Politecnica Valencia
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

import cStringIO
from array import array
import itertools
import hashlib

from crumbs.utils.optional_modules import FastqGeneralIterator
from crumbs.settings import get_setting
from crumbs.utils.file_utils import fhand_is_seekable, peek_chunk_from_file
from crumbs.exceptions import (UnknownFormatError, UndecidedFastqVersionError,
                               FileIsEmptyError)

FILETYPE = 'file'
STRINGIOTYPE = 'stringio'
OTHERTYPE = 'othertype'


def _get_some_qual_and_lengths(fhand, force_file_as_non_seek):
    'It returns the quality characters and the lengths'
    seqs_to_peek = get_setting('SEQS_TO_GUESS_FASTQ_VERSION')
    chunk_size = get_setting('CHUNK_TO_GUESS_FASTQ_VERSION')

    lengths = array('I')
    seqs_analyzed = 0
    if fhand_is_seekable(fhand) and not force_file_as_non_seek:
        fmt_fhand = fhand
        chunk = fmt_fhand.read(chunk_size)
        fhand.seek(0)
    else:
        chunk = peek_chunk_from_file(fhand, chunk_size)
        fmt_fhand = cStringIO.StringIO(chunk)

    try:
        for seq in FastqGeneralIterator(fmt_fhand):
            qual = [ord(char) for char in seq[2]]
            sanger_chars = [q for q in qual if q < 64]
            if sanger_chars:
                fhand.seek(0)
                return None, True, chunk  # no quals, no lengths, is_sanger
            lengths.append(len(qual))
            seqs_analyzed += 1
            if seqs_analyzed > seqs_to_peek:
                break
    except ValueError:
        msg = 'The file is Fastq, but the version is difficult to guess'
        raise UndecidedFastqVersionError(msg)
    finally:
        fhand.seek(0)
    return lengths, None, chunk  # don't know if it's sanger


def _guess_fastq_version(fhand, force_file_as_non_seek):
    '''It guesses the format of fastq files.

    It ignores the solexa fastq version.
    '''
    lengths, is_sanger, chunk = _get_some_qual_and_lengths(fhand,
                                                        force_file_as_non_seek)
    if is_sanger:
        fmt = 'fastq'
    elif is_sanger is False:
        fmt = 'fastq-illumina'
    else:
        fmt = None
    if fmt:
        return fmt

    longest_expected_illumina = get_setting('LONGEST_EXPECTED_ILLUMINA_READ')
    n_long_seqs = [l for l in lengths if l > longest_expected_illumina]
    if n_long_seqs:
        msg = 'It was not possible to guess the format of '
        if hasattr(fhand, 'name'):
            msg += 'the file ' + fhand.name
        else:
            msg += 'a file '
        msg = '\n. The quality values could be Illumina, but there are '
        msg += 'sequences longer than %i bp.'
        msg %= longest_expected_illumina
        raise UndecidedFastqVersionError(msg)
    else:
        return 'fastq-illumina'

FILEFORMAT_INVENTORY = {}


def _get_instance_type(fhand):
    try:
        name = fhand.name
    except AttributeError:
        name = None
    if name:
        return FILETYPE

    if 'getvalue' in dir(fhand):
        return STRINGIOTYPE
    else:
        return OTHERTYPE


def _get_fhand_id(fhand):
    'get unique id of the fhand'
    fhand_type = _get_instance_type(fhand)
    fhand_id = id(fhand)
    if fhand_type == FILETYPE:
        key = (fhand_id, fhand.name)
    elif fhand_type == STRINGIOTYPE:
        hash_ = hashlib.sha224(fhand.getvalue()[:100]).hexdigest()
        key = (fhand_id, hash_)
    else:
        key = fhand_id
    return key


def get_format(fhand):
    'It gets the format or it looks in the inventory'
    id_ = _get_fhand_id(fhand)
    try:
        file_format = FILEFORMAT_INVENTORY[id_]
    except KeyError:
        file_format = None

    if file_format is None:
        file_format = _guess_format(fhand, force_file_as_non_seek=False)
        FILEFORMAT_INVENTORY[id_] = file_format

    return file_format


def set_format(fhand, file_format):
    'It sets the file format in the global inventory variable'
    id_ = _get_fhand_id(fhand)
    if id_ in FILEFORMAT_INVENTORY:
        msg = 'The given instance already setted its file format'
        raise RuntimeError(msg)
    FILEFORMAT_INVENTORY[id_] = file_format


def _guess_format(fhand, force_file_as_non_seek):
    '''It guesses the format of the sequence file.

    This function is just for testing forcing the fhand as non-seekable.
    It does ignore the solexa fastq version.
    '''
    chunk_size = 2048
    chunk = peek_chunk_from_file(fhand, chunk_size)
    if not chunk:
        raise FileIsEmptyError('The file is empty')
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
