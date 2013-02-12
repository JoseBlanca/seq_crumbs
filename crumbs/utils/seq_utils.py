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
from copy import deepcopy

from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from crumbs.exceptions import UnknownFormatError, UndecidedFastqVersionError
from crumbs.settings import get_setting
from crumbs.utils.file_utils import fhand_is_seekable, peek_chunk_from_file
from crumbs.utils.tags import (UPPERCASE, LOWERCASE, SWAPCASE, SEQITEM,
                               SEQRECORD)

# pylint: disable=R0903
# pylint: disable=C0111


def replace_seq_same_length(seqrecord, seq_str):
    'It replaces the str with another of equal length keeping the annots.'
    annots = seqrecord.letter_annotations
    seqrecord.letter_annotations = {}
    alphabet = seqrecord.seq.alphabet
    seqrecord.seq = Seq(seq_str, alphabet)
    seqrecord.letter_annotations = annots
    return seqrecord


def copy_seqrecord(seqrec, seq=None, name=None, id_=None):
    '''Given a seqrecord it returns a new seqrecord with seq or qual changed.

    This is necessary because our SeqWithQuality is inmutable
    '''
    if seq is None:
        seq = seqrec.seq
    if id_ is  None:
        id_ = seqrec.id
    if name is None:
        name = seqrec.name

    # the letter annotations
    let_annot = {annot: v for annot, v in seqrec.letter_annotations.items()}

    # the rest of parameters
    description = seqrec.description
    dbxrefs = seqrec.dbxrefs[:]
    features = seqrec.features[:]  # the features are not copied
    annotations = deepcopy(seqrec.annotations)

    # the new sequence
    new_seq = SeqRecord(seq=seq, id=id_, name=name, description=description,
                        dbxrefs=dbxrefs, features=features,
                        annotations=annotations, letter_annotations=let_annot)

    return new_seq


def uppercase_length(string):
    'It returns the number of uppercase characters found in the string'
    return len(re.findall("[A-Z]", string))


def get_uppercase_segments(string):
    '''It detects the unmasked regions of a sequence

    It returns a list of (start, end) tuples'''
    start = 0
    for is_upper, group in itertools.groupby(string, lambda x: x.isupper()):
        group = list(group)
        end = start + len(group) - 1
        if is_upper:
            yield start, end
        start = end + 1


class ChangeCase(object):
    'It changes the sequence case.'

    def __init__(self, action):
        'The initiator'
        if action not in (UPPERCASE, LOWERCASE, SWAPCASE):
            msg = 'Action should be: uppercase, lowercase or invertcase'
            raise ValueError(msg)
        self.action = action

    def __call__(self, seqrecords):
        'It changes the case of the seqrecords.'
        action = self.action
        processed_seqs = []
        for seqrecord in seqrecords:
            str_seq = str(seqrecord.seq)
            if action == UPPERCASE:
                str_seq = str_seq.upper()
            elif action == LOWERCASE:
                str_seq = str_seq.lower()
            elif action == SWAPCASE:
                str_seq = str_seq.swapcase()
            else:
                raise NotImplementedError()
            seqrecord = replace_seq_same_length(seqrecord, str_seq)
            processed_seqs.append(seqrecord)
        return processed_seqs


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
            one_line = ''
        else:
            one_line = '-one_line'

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


def append_to_description(seqrecord, text):
    'it appends the text to the seqrecord description'
    desc = seqrecord.description
    if desc == '<unknown description>' or desc == seqrecord.id:
        desc = ''
    desc += text
    seqrecord.description = desc


class _FunctionRunner(object):
    'a class to join all the mapper functions in a single function'
    def __init__(self, map_functions):
        'Class initiator'
        self.map_functions = map_functions

    def __call__(self, seq_packet):
        'It runs all the map_functions for each seq_packet '
        processed_packet = seq_packet
        for map_function in self.map_functions:
            processed_packet = map_function(processed_packet)
        return processed_packet


def process_seq_packets(seq_packets, map_functions, processes=1,
                        keep_order=True):
    'It processes the SeqRecord packets'
    if processes > 1:
        workers = Pool(processes=processes)
        mapper = workers.imap if keep_order else workers.imap_unordered
    else:
        workers = None
        mapper = itertools.imap
    run_functions = _FunctionRunner(map_functions)

    seq_packets = mapper(run_functions, seq_packets)

    return seq_packets, workers


def get_title(seq):
    'Given a seq it returns the title'
    # TODO remove this check when everything is adapted to the new system
    if 'SeqRecord' in seq.__class__.__name__:
        seq_class = SEQRECORD
    else:
        seq_class = seq.kind
        seq = seq.object

    if seq_class == SEQITEM:
        title = seq.lines[0][1:]
    elif seq_class == SEQRECORD:
        title = seq.id + ' ' + seq.description
    else:
        msg = 'Do not know how to guess title form this seq class'
        raise NotImplementedError(msg)
    return title


def get_name(seq):
    if 'SeqRecord' in seq.__class__.__name__:
        seq_class = SEQRECORD
    else:
        seq_class = seq.kind
        seq = seq.object
    if seq_class == SEQITEM:
        name = seq.name
    elif seq_class == SEQRECORD:
        name = seq.id
    return name


def get_file_format(seq):
    seq_class = seq.kind

    if seq_class == SEQITEM:
        fmt = seq.file_format
    elif seq_class == SEQRECORD:
        fmt = None
    return fmt
