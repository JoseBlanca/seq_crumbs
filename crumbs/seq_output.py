'''
Created on 18/06/2012

@author: jose
'''

from Bio import SeqIO


def write_seqrecords(fhand, seqs):
    'It writes a stream of sequences to a file'
    SeqIO.write(seqs, fhand, 'fastq')
