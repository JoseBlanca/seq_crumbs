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


from itertools import chain
from shutil import copyfileobj
from tempfile import NamedTemporaryFile
import cStringIO

from Bio import SeqIO
from Bio.SeqIO import FastaIO
from Bio.SeqIO import QualityIO
from Bio.Alphabet import IUPAC

from crumbs.exceptions import (MalformedFile, error_quality_disagree,
                               UnknownFormatError)
from crumbs.iterutils import length, group_in_packets
from crumbs.utils.file_utils import rel_symlink
from crumbs.utils.seq_utils import guess_format, peek_chunk_from_file
from crumbs.utils.tags import (GUESS_FORMAT, SEQS_PASSED, SEQS_FILTERED_OUT)
from crumbs.settings import PACKET_SIZE


def clean_seq_stream(seqs):
    'It removes the empty seqs and fixes the descriptions.'
    for seq in seqs:
        if seq and str(seq.seq):
            if seq.description == '<unknown description>':
                seq.description = ''
            yield seq


def write_seqrecords(seqs, fhand=None, file_format='fastq'):
    'It writes a stream of sequences to a file'
    if fhand is None:
        fhand = NamedTemporaryFile(suffix='.' + file_format.replace('-', '_'))
    seqs = clean_seq_stream(seqs)
    SeqIO.write(seqs, fhand, file_format)
    return fhand


def write_seq_packets(fhand, seq_packets, file_format='fastq', workers=None):
    'It writes to file a stream of SeqRecord lists'
    try:
        write_seqrecords(chain.from_iterable(seq_packets), fhand,
                         file_format=file_format)
    except BaseException:
        if workers is not None:
            workers.terminate()
        raise


def write_filter_packets(passed_fhand, filtered_fhand, filter_packets,
                         file_format='fastq', workers=None):
    'It writes the filter stream into passed and filtered out sequence files'
    if filtered_fhand is None:
        seq_packets = (p[SEQS_PASSED] for p in filter_packets)
        return write_seq_packets(fhand=passed_fhand, seq_packets=seq_packets,
                                 file_format=file_format, workers=workers)
    for packet in filter_packets:
        try:
            write_seqrecords(packet[SEQS_PASSED], fhand=passed_fhand,
                             file_format=file_format)
            write_seqrecords(packet[SEQS_FILTERED_OUT], fhand=filtered_fhand,
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


def read_seq_packets(fhands, size=PACKET_SIZE, file_format=GUESS_FORMAT):
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


def seqio(in_fhands, out_fhands, out_format, copy_if_same_format=True):
    'It converts sequence files between formats'

    in_formats = [guess_format(fhand) for fhand in in_fhands]

    if (len(in_formats) == 1 and in_formats[0] == out_format and
        hasattr(in_fhands[0], 'name')):
        if copy_if_same_format:
            copyfileobj(in_fhands[0], out_fhands[0])
        else:
            rel_symlink(in_fhands[0].name, out_fhands[0].name)

    elif len(in_fhands) == 1 and len(out_fhands) == 1:
        try:
            SeqIO.convert(in_fhands[0], in_formats[0], out_fhands[0],
                          out_format)
        except ValueError as error:
            if error_quality_disagree(error):
                raise MalformedFile(str(error))
            raise
    elif (len(in_fhands) == 1 and len(out_fhands) == 2 and
          out_format == 'fasta'):
        try:
            for seq in read_seqrecords([in_fhands[0]]):
                SeqIO.write([seq], out_fhands[0], out_format)
                SeqIO.write([seq], out_fhands[1], 'qual')
        except ValueError, error:
            if error_quality_disagree(error):
                raise MalformedFile(str(error))
            raise
    elif (len(in_fhands) == 2 and len(out_fhands) == 1 and
          in_formats == ['fasta', 'qual']):
        seq_records = SeqIO.QualityIO.PairedFastaQualIterator(in_fhands[0],
                                                              in_fhands[1])
        try:
            SeqIO.write(seq_records, out_fhands[0].name, out_format)
        except ValueError, error:
            if error_quality_disagree(error):
                raise MalformedFile(str(error))
            raise
    elif (len(in_fhands) == 2 and len(out_fhands) == 2 and
          in_formats == ['fasta', 'qual'] and out_format == 'fasta'):
        if copy_if_same_format:
            copyfileobj(in_fhands[0], out_fhands[0])
            copyfileobj(in_fhands[1], out_fhands[1])
        else:
            rel_symlink(in_fhands[0].name, out_fhands[0].name)
            rel_symlink(in_fhands[1].name, out_fhands[1].name)
    else:
        raise RuntimeError('Please fixme, we should not be here')

    for out_fhand in out_fhands:
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
