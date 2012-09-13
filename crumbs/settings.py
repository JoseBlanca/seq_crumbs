'''
Created on 19/07/2012

@author: jose
'''

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

SUPPORTED_OUTPUT_FORMATS = ['fasta', 'fastq', 'fastq-illumina']

# number of sequences to process in a chunk. Lenght of the sequence list to
# hold in memory
PACKET_SIZE = 1000

# number of sequences to analyze in the fastq version guessing of a seekable
# file
SEQS_TO_GUESS_FASTQ_VERSION = 1000

# number of bytes to analyze in the fastq version guessing of a non-seekable
# file
CHUNK_TO_GUESS_FASTQ_VERSION = 50000

# maximum length expected for an Illumina read
LONGEST_EXPECTED_ILLUMINA_READ = 250


# 454 FLX mate pair linker
FLX_LINKER = 'GTTGGAACCGAAAGGGTTTGAATTCAAACCCTTTCGGTTCCAAC'
# Titanium mate pair linker. It could be found forward or reverse
TITANIUM_LINKER = 'TCGTATAACTTCGTATAATGTATGCTATACGAAGTTATTACG'
TITANIUM_LINKER_REV = 'CGTAATAACTTCGTATAGCATACATTATACGAAGTTATACGA'
FWD_454_LINKERS = [FLX_LINKER, TITANIUM_LINKER]
LINKERS = [SeqRecord(Seq(FLX_LINKER), id='flx_linker'),
           SeqRecord(Seq(TITANIUM_LINKER), id='titanium_linker')]


## Use this to modify how get_binary path works
# if need to modify the binary's name
USE_EXTERNAL_BIN_PREFIX = False
# prefix to add to the binary name
EXTERNAL_BIN_PREFIX = 'crumbs_'
# mark True if need the path or assumes that is on the path
ADD_PATH_TO_EXT_BIN = True

# how many reads can be hold in memory by default
DEFAULT_SEQS_IN_MEM_LIMIT = 500000
