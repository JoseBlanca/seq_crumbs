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

import re
import cStringIO
from array import array
import itertools
from multiprocessing import Pool

from Bio.SeqIO.QualityIO import FastqGeneralIterator

from crumbs.exceptions import UnknownFormatError, UndecidedFastqVersionError
from crumbs.settings import (CHUNK_TO_GUESS_FASTQ_VERSION,
                             SEQS_TO_GUESS_FASTQ_VERSION,
                             LONGEST_EXPECTED_ILLUMINA_READ)
from crumbs.utils.file_utils import fhand_is_seekable, peek_chunk_from_file


def uppercase_length(string):
    'It returns the number of uppercase characters found in the string'
    return len(re.findall("[A-Z]", string))


def _get_some_qual_and_lengths(fhand, force_file_as_non_seek):
    'It returns the quality characters and the lengths'
    seqs_to_peek = SEQS_TO_GUESS_FASTQ_VERSION
    chunk_size = CHUNK_TO_GUESS_FASTQ_VERSION

    lengths = array('I')
    seqs_analyzed = 0
    if fhand_is_seekable(fhand) and not force_file_as_non_seek:
        fmt_fhand = fhand
    else:
        chunk = peek_chunk_from_file(fhand, chunk_size)
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


def process_seq_packets(seq_packets, map_functions, processes=1,
                        keep_order=False):
    'It processes the SeqRecord packets'
    if processes > 1:
        pool = Pool(processes=processes)
        mapper = pool.imap if keep_order else pool.imap_unordered
    else:
        mapper = itertools.imap

    for map_function in map_functions:
        seq_packets = mapper(map_function, seq_packets)

    return seq_packets
