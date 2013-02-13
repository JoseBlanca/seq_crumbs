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


from itertools import chain, tee
from shutil import copyfileobj
from tempfile import NamedTemporaryFile
import cStringIO
from collections import namedtuple

from Bio import SeqIO
from Bio.SeqIO import FastaIO
from Bio.SeqIO import QualityIO
from Bio.Alphabet import IUPAC

from crumbs.exceptions import (MalformedFile, error_quality_disagree,
                               UnknownFormatError, IncompatibleFormatError)
from crumbs.iterutils import length, group_in_packets
from crumbs.utils.file_utils import rel_symlink
from crumbs.utils.seq_utils import guess_format, peek_chunk_from_file
from crumbs.utils.tags import (GUESS_FORMAT, SEQS_PASSED, SEQS_FILTERED_OUT,
                               SEQITEM, SEQRECORD)
from crumbs.settings import get_setting


def clean_seqrecord_stream(seqs):
    'It removes the empty seqs and fixes the descriptions.'
    for seq in seqs:
        if seq and str(seq.seq):
            if seq.description == '<unknown description>':
                seq.description = ''
            yield seq


def write_seqrecords(seqs, fhand=None, file_format='fastq'):
    'It writes a stream of sequences to a file'
    file_format = _remove_one_line(file_format)

    if fhand is None:
        fhand = NamedTemporaryFile(suffix='.' + file_format.replace('-', '_'))
    seqs = clean_seqrecord_stream(seqs)
    SeqIO.write(seqs, fhand, file_format)
    return fhand


def write_seqrecord_packets(fhand, seq_packets, file_format='fastq',
                            workers=None):
    'It writes to file a stream of SeqRecord lists'
    try:
        write_seqrecords(chain.from_iterable(seq_packets), fhand,
                         file_format=file_format)
    except BaseException:
        if workers is not None:
            workers.terminate()
        raise


def write_seq_packets(fhand, seq_packets, file_format='fastq', workers=None):
    'It writes to file a stream of seq lists'
    file_format = _remove_one_line(file_format)
    try:
        write_seqs(chain.from_iterable(seq_packets), fhand,
                   file_format=file_format)
    except BaseException:
        if workers is not None:
            workers.terminate()
        raise


def write_filter_packets(passed_fhand, filtered_fhand, filter_packets,
                         file_format='fastq', workers=None):
    'It writes the filter stream into passed and filtered out sequence files'
    file_format = _remove_one_line(file_format)

    if filtered_fhand is None:
        seq_packets = (p[SEQS_PASSED] for p in filter_packets)
        return write_seq_packets(fhand=passed_fhand, seq_packets=seq_packets,
                                 file_format=file_format, workers=workers)
    for packet in filter_packets:
        try:
            write_seqs(packet[SEQS_PASSED], fhand=passed_fhand,
                       file_format=file_format)
            write_seqs(packet[SEQS_FILTERED_OUT], fhand=filtered_fhand,
                       file_format=file_format)
        except BaseException:
            if workers is not None:
                workers.terminate()
            raise


def title2ids(title):
    '''It returns the id, name and description as a tuple.

    It takes the title of the FASTA file (without the beginning >)
    '''
    items = title.strip().split()
    name = items[0]
    id_ = name
    if len(items) > 1:
        desc = ' '.join(items[1:])
    else:
        desc = ''
    return id_, name, desc


def read_seq_packets(fhands, size=get_setting('PACKET_SIZE'), out_format=None,
                     file_format=GUESS_FORMAT, prefered_seq_classes=None):
    '''It yields Seqitems in packets of the given size.'''
    seqs = read_seqs(fhands, file_format, out_format=out_format,
                     prefered_seq_classes=prefered_seq_classes)
    return group_in_packets(seqs, size)


def read_seqrecord_packets(fhands, size=get_setting('PACKET_SIZE'),
                           file_format=GUESS_FORMAT):
    '''It yields SeqRecords in packets of the given size.'''
    seqs = read_seqrecords(fhands, file_format=file_format)
    return group_in_packets(seqs, size)


def read_seqrecords(fhands, file_format=GUESS_FORMAT):
    'it returns an iterator of seqrecords'
    seq_iters = []
    for fhand in fhands:
        if file_format == GUESS_FORMAT or file_format is None:
            fmt = guess_format(fhand)
        else:
            fmt = file_format

        fmt = _remove_one_line(fmt)

        if fmt in ('fasta', 'qual') or 'fastq' in fmt:
            title = title2ids
        if fmt == 'fasta':
            seq_iter = FastaIO.FastaIterator(fhand, title2ids=title)
        elif fmt == 'qual':
            seq_iter = QualityIO.QualPhredIterator(fhand, title2ids=title)
        elif fmt == 'fastq' or fmt == 'fastq-sanger':
            seq_iter = QualityIO.FastqPhredIterator(fhand, title2ids=title)
        elif fmt == 'fastq-solexa':
            seq_iter = QualityIO.FastqSolexaIterator(fhand, title2ids=title)
        elif fmt == 'fastq-illumina':
            seq_iter = QualityIO.FastqIlluminaIterator(fhand, title2ids=title)
        else:
            seq_iter = SeqIO.parse(fhand, fmt)
        seq_iters.append(seq_iter)
    return chain.from_iterable(seq_iters)


def _remove_one_line(file_format):
    'It removes the one-line from the format'
    if file_format and file_format.endswith('-one_line'):
        file_format = file_format[:-9]
    return file_format


def seqio(in_fhands, out_fhand, out_format, copy_if_same_format=True):
    'It converts sequence files between formats'
    if out_format not in get_setting('SUPPORTED_OUTPUT_FORMATS'):
        raise IncompatibleFormatError("This output format is not supported")

    in_formats = [_remove_one_line(guess_format(fhand)) for fhand in in_fhands]

    if len(in_fhands) == 1 and in_formats[0] == out_format:
        if copy_if_same_format:
            copyfileobj(in_fhands[0], out_fhand)
        else:
            rel_symlink(in_fhands[0].name, out_fhand.name)
    else:
        seqs = read_seqrecords(in_fhands)
        try:
            SeqIO.write(seqs, out_fhand, out_format)
        except ValueError, error:
            if error_quality_disagree(error):
                raise MalformedFile(str(error))
            if 'No suitable quality scores' in str(error):
                msg = 'No qualities available to write output file'
                raise IncompatibleFormatError(msg)
            raise
    out_fhand.flush()


def fastaqual_to_fasta(seq_fhand, qual_fhand, out_fhand):
    'It converts a fasta and a qual file into a fastq format file'
    seqrecords = SeqIO.QualityIO.PairedFastaQualIterator(seq_fhand, qual_fhand)
    try:
        SeqIO.write(seqrecords, out_fhand.name, 'fastq')
    except ValueError, error:
        if error_quality_disagree(error):
            raise MalformedFile(str(error))
        raise
    out_fhand.flush()


def _count_seqs_in_fasta(fhand):
    'It counts the seqs in a fasta file'
    count = 0
    for line in fhand:
        if line[0] == '>':
            count += 1
    return count


def count_seqs_in_files(fhands, file_format=GUESS_FORMAT):
    'It counts the seqs in the given files'
    count = 0
    for fhand in fhands:
        if file_format == GUESS_FORMAT or file_format is None:
            file_format = guess_format(fhand)
        else:
            file_format = file_format

        if file_format == 'fasta':
            count += _count_seqs_in_fasta(fhand)
        elif 'fastq' in file_format:
            count += length(QualityIO.FastqGeneralIterator(fhand))
        else:
            count += length(read_seqrecords([fhand]))
    return count


def guess_seq_type(fhand):
    '''it guesses the file's seq type'''

    rna = set(IUPAC.ambiguous_rna.letters)
    dna = set(IUPAC.ambiguous_dna.letters)
    rna_dna = rna.union(dna)

    protein = set(IUPAC.extended_protein.letters)
    only_prot = list(protein.difference(rna_dna))

    chunk_size = 1024
    chunk = peek_chunk_from_file(fhand, chunk_size)
    if not chunk:
        raise UnknownFormatError('The file is empty')
    fhand_ = cStringIO.StringIO(chunk)
    total_letters = 0
    nucleotides = 0
    for seqrec in read_seqrecords([fhand_]):
        for letter in str(seqrec.seq):
            total_letters += 1
            if letter in ('gcatnuGCATNU'):
                nucleotides += 1
            if letter in only_prot:
                return 'prot'
    nucl_freq = nucleotides / total_letters
    if nucl_freq > 0.8:
        return 'nucl'

    raise RuntimeError('unable to guess the seq type')


def _get_name_from_lines(lines):
    'It returns the name and the chunk from a list of names'
    name = lines[0].split()[0][1:]
    return name


SeqWrapper = namedtuple('Seq', ['kind', 'object', 'file_format'])
SeqItem = namedtuple('SeqItem', ['name', 'lines', 'length'])


def _itemize_fasta(fhand):
    'It returns the fhand divided in chunks, one per seq'

    lines = []
    for line in fhand:
        if not line or line.isspace():
            continue
        if line.startswith('>'):
            if lines:
                yield SeqItem(_get_name_from_lines(lines), lines,
                              len(lines[1]) - 1)
                lines = []
        lines.append(line)
        if len(lines) == 1 and not lines[0].startswith('>'):
            raise RuntimeError('Not a valid fasta file')
    else:
        if lines:
            yield SeqItem(_get_name_from_lines(lines), lines,
                          len(lines[1]) - 1)


def _get_name_from_chunk_fastq(lines):
    'It returns the name and the chunk from a list of names'
    if len(lines) != 4:
        raise RuntimeError('Malformed fastq')
    if not lines[0].startswith('@'):
        raise RuntimeError('Not a valid fastq file: not start with @')
    if not lines[1].startswith('+'):
        raise RuntimeError('Too complex fastq for this function')
    if len(lines[1]) != len(lines[3]):
        raise RuntimeError('Qual has different length than seq')
    name = lines[0].split()[0][1:]
    return name


def _itemize_fastq(fhand):
    'It returns the fhand divided in chunks, one per seq'
    blobs = group_in_packets(fhand, 4)
    return (SeqItem(_get_name_from_lines(lines), lines, len(lines[1]) - 1) for lines in blobs)


def assing_kind_to_seqs(kind, seqs, file_format):
    'It puts each seq into a NamedTuple named Seq'
    return (SeqWrapper(kind, seq, file_format) for seq in seqs)


def _read_seqitems(fhands, file_format):
    'it returns an iterator of seq items (tuples of name and chunk)'
    seq_iters = []
    for fhand in fhands:
        if file_format == GUESS_FORMAT or file_format is None:
            file_format = guess_format(fhand)
        else:
            file_format = file_format

        if file_format == 'fasta':
            seq_iter = _itemize_fasta(fhand)
        elif file_format in ('fastq-one_line', 'fastq-illumina-one_line'):
            seq_iter = _itemize_fastq(fhand)
        else:
            msg = 'Format not supported by the itemizers: ' + file_format
            raise NotImplementedError(msg)
        seq_iter = assing_kind_to_seqs(SEQITEM, seq_iter, file_format)
        seq_iters.append(seq_iter)
    return chain.from_iterable(seq_iters)


def _write_seqitems(items, fhand, file_format):
    'It writes one seq item (tuple of name and string)'
    for seq in items:
        seqitems_fmt = _remove_one_line(seq.file_format)
        if file_format and 'fastq' in seqitems_fmt and 'fasta' in file_format:
            seq_lines = seq.object.lines
            fhand.write('>' + seq_lines[0][1:] + seq_lines[1])
        elif file_format and seqitems_fmt != file_format:
            msg = 'Input and output file formats do not match, you should not '
            msg += 'use SeqItems: ' + str(seq.file_format) + ' '
            msg += str(file_format)
            raise RuntimeError(msg)
        else:
            fhand.write(''.join(seq.object.lines))


def write_seqs(seqs, fhand=None, file_format=None):
    'It writes the given sequences'
    if fhand is None:
        fhand = NamedTemporaryFile(suffix='.' + file_format.replace('-', '_'))

    file_format = _remove_one_line(file_format)
    seqs, seqs2 = tee(seqs)
    try:
        seq = seqs2.next()
    except StopIteration:
        # No sequences to write, so we're done
        return
    seq_class = seq.kind
    if seq_class == SEQITEM:
        _write_seqitems(seqs, fhand, file_format)
    elif seq_class == SEQRECORD:
        seqs = (seq.object for seq in seqs)
        write_seqrecords(seqs, fhand, file_format)
    else:
        raise ValueError('Unknown class for seq: ' + seq_class)
    return fhand


def read_seqs(fhands, file_format, out_format=None, prefered_seq_classes=None):
    'It returns a stream of seqs in different codings: seqrecords, seqitems...'

    if not prefered_seq_classes:
        prefered_seq_classes = [SEQITEM, SEQRECORD]

    if file_format == GUESS_FORMAT:
        in_format = guess_format(fhands[0])
    else:
        in_format = file_format

    if out_format not in (None, GUESS_FORMAT):

        if in_format != out_format:
            if SEQITEM in prefered_seq_classes:
                # seqitems is incompatible with different input and output
                # formats
                prefered_seq_classes.pop(prefered_seq_classes.index(SEQITEM))

    if not prefered_seq_classes:
        msg = 'No valid seq class left or prefered'
        raise ValueError(msg)

    for seq_class in prefered_seq_classes:
        if seq_class == SEQITEM:
            try:
                return _read_seqitems(fhands, in_format)
            except NotImplementedError:
                continue
        elif seq_class == SEQRECORD:
            try:
                seqs = read_seqrecords(fhands, in_format)
                return assing_kind_to_seqs(SEQRECORD, seqs, None)
            except NotImplementedError:
                continue
        else:
            raise ValueError('Unknown class for seq: ' + seq_class)
    raise RuntimeError('We should not be here, fixme')
