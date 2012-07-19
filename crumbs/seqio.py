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
from itertools import chain
from shutil import copyfileobj


from Bio import SeqIO
from Bio.SeqIO import FastaIO
from Bio.SeqIO import QualityIO

from crumbs.exceptions import MalformedFile
from crumbs.iterutils import length, group_in_packets
from crumbs.utils import rel_symlink, guess_format
from crumbs.settings import GUESS_FORMAT

PACKET_SIZE = 1000


def clean_seq_stream(seqs):
    'It removes the empty seqs and fixes the descriptions.'
    for seq in seqs:
        if seq and str(seq.seq):
            if seq.description == '<unknown description>':
                seq.description = ''
            yield seq


def write_seqrecords(fhand, seqs, file_format='fastq'):
    'It writes a stream of sequences to a file'
    seqs = clean_seq_stream(seqs)
    SeqIO.write(seqs, fhand, file_format)
    fhand.flush()


def write_seqrecord_packets(fhand, seq_packets, file_format='fastq'):
    'It writes to file a stream of seqcuence lists'
    write_seqrecords(fhand, chain.from_iterable(seq_packets),
                     file_format=file_format)


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


def read_seqrecords_in_packets(fhands, size=PACKET_SIZE,
                               file_format=GUESS_FORMAT):
    '''it returns the sequences in packets of the given size. It returns a
    iterator of seq iterators'''
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
            msg = 'Lengths of sequence and quality values differs'
            if msg in str(error):
                raise MalformedFile(str(error))
            raise
    elif (len(in_fhands) == 1 and len(out_fhands) == 2 and
          out_format == 'fasta'):
        try:
            for seq in read_seqrecords([in_fhands[0]]):
                SeqIO.write([seq], out_fhands[0], out_format)
                SeqIO.write([seq], out_fhands[1], 'qual')
        except ValueError, error:
            msg = 'Sequence length and number of quality scores disagree'
            if msg in str(error):
                raise MalformedFile(str(error))
            raise
#        try:
#            SeqIO.convert(in_fhands[0], in_formats[0], out_fhands[0], 'fasta')
#            SeqIO.convert(in_fhands[0], in_formats[0], out_fhands[1], 'qual')
#        except ValueError as error:
#            msg = 'Lengths of sequence and quality values differs'
#            if msg in str(error):
#                raise MalformedFile(str(error))
#            raise

    elif (len(in_fhands) == 2 and len(out_fhands) == 1 and
          in_formats == ['fasta', 'qual']):
        seq_records = SeqIO.QualityIO.PairedFastaQualIterator(in_fhands[0],
                                                              in_fhands[1])
        try:
            SeqIO.write(seq_records, out_fhands[0].name, out_format)
        except ValueError, error:
            msg = 'Sequence length and number of quality scores disagree'
            if msg in str(error):
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
