'''
Created on 21/06/2012

@author: jose
'''
from itertools import chain

from shutil import copyfileobj

from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio import SeqIO

from crumbs.utils.exceptions import UnknownFormatError, MalformedFile


def write_seqrecords(fhand, seqs, file_format='fastq'):
    'It writes a stream of sequences to a file'
    SeqIO.write(seqs, fhand, file_format)


def read_seqrecords(fhands):
    'it returns an iterator of seqrecords'
    return chain(*[SeqIO.parse(fhand, guess_format(fhand)) for fhand in fhands])


def _guess_fastq_format(fhand):
    '''It guesses the format of fastq files.

    It ignores the solexa fastq version.
    '''
    try:
        for seq in FastqGeneralIterator(fhand):
            qual = seq[2]
            for letter in qual:
                num = ord(letter)
                if num < 64:
                    fhand.seek(0)
                    return 'fastq'
                elif num > 73:
                    fhand.seek(0)
                    return 'fastq-illumina'
    except ValueError:
        raise UnknownFormatError('Malformed fastq')
    msg = 'It is fastq file, but we do not know if sanger or illumina'
    raise UnknownFormatError(msg)


def guess_format(fhand):
    '''It guesses the format of the sequence file.

    It does ignore the solexa fastq version.
    '''
    fhand.seek(0)
    chunk = fhand.read(1024)
    fhand.seek(0)
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
        return _guess_fastq_format(fhand)
    elif chunk.startswith('LOCUS'):
        return 'genbank'
    elif chunk.startswith('ID'):
        return 'embl'
    raise UnknownFormatError('Sequence file of unknown format.')


def seqio(in_fhands, out_fhands, out_format):
    'It converts sequence files between formats'

    in_formats = [guess_format(fhand) for fhand in in_fhands]

    if len(in_fhands) == 1 and len(out_fhands) == 1:
        try:
            SeqIO.convert(in_fhands[0].name, in_formats[0],
                          out_fhands[0].name, out_format)
        except ValueError as error:
            msg = 'Lengths of sequence and quality values differs'
            if msg in str(error):
                raise MalformedFile(str(error))
            raise
    elif (len(in_fhands) == 1 and len(out_fhands) == 2 and
          out_format == 'fasta'):
        try:
            SeqIO.convert(in_fhands[0].name, in_formats[0],
                          out_fhands[0].name, 'fasta')
            SeqIO.convert(in_fhands[0].name, in_formats[0],
                          out_fhands[1].name, 'qual')
        except ValueError as error:
            msg = 'Lengths of sequence and quality values differs'
            if msg in str(error):
                raise MalformedFile(str(error))
            raise

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
        #out_fhands[0].write(in_fhands[0].read())
        copyfileobj(in_fhands[0], out_fhands[0])
        copyfileobj(in_fhands[1], out_fhands[1])

    else:
        raise RuntimeError('Please fixme, we should not be here')

    for out_fhand in out_fhands:
        out_fhand.flush()
