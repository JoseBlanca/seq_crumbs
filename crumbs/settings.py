'''
Created on 19/07/2012

@author: jose
'''

STDIN = 'stdin'
STDOUT = 'stdout'
INFILES = 'infiles'
OUTFILE = 'output'
SUPPORTED_OUTPUT_FORMATS = ['fasta', 'fastq', 'fastq-illumina']

# number of sequences to analyze in the fastq version guessing of a seekable
# file
SEQS_TO_GUESS_FASTQ_VERSION = 1000
# number of bytes to analyze in the fastq version guessing of a non-seekable
# file
CHUNK_TO_GUESS_FASTQ_VERSION = 50000
# maximum length expected for an Illumina read
LONGEST_EXPECTED_ILLUMINA_READ = 250
# Input format tag when we want to guess
GUESS_FORMAT = 'guess'
