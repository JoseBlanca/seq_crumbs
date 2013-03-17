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

from Bio.SeqIO.QualityIO import FastqGeneralIterator

from crumbs.settings import get_setting
from crumbs.utils.file_utils import fhand_is_seekable, peek_chunk_from_file
from crumbs.exceptions import UnknownFormatError, UndecidedFastqVersionError


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
        raise UnknownFormatError('Malformed fastq')
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

    # onle line fastq? All seq in just one line?
    lines = [l for l in itertools.islice(chunk.splitlines(), 5)]
    if len(lines) != 5:
        one_line = ''
    else:
        if (not lines[0].startswith('@') or
            not lines[2].startswith('+') or
            not lines[4].startswith('@') or
            len(lines[1]) != len(lines[3])):
            one_line = '-multiline'
        else:
            one_line = ''

    if fmt:
        return fmt + one_line

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
        return 'fastq-illumina' + one_line


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
    chunk = peek_chunk_from_file(fhand, chunk_size)
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
