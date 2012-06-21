'''
Created on 21/06/2012

@author: jose
'''
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio import SeqIO

from crumbs.utils.exceptions import UnknownFormatError


def write_seqrecords(fhand, seqs):
    'It writes a stream of sequences to a file'
    SeqIO.write(seqs, fhand, 'fastq')


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
            return 'fasta'
    elif chunk.startswith('@'):
        return _guess_fastq_format(fhand)
    elif chunk.startswith('LOCUS'):
        return 'genbank'
    elif chunk.startswith('ID'):
        return 'embl'
    raise UnknownFormatError('Sequence file of unknown format.')
