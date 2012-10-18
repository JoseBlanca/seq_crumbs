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
USE_EXTERNAL_BIN_PREFIX = True
# prefix to add to the binary name
EXTERNAL_BIN_PREFIX = 'crumbs_'
# mark True if need the path or assumes that is on the path
ADD_PATH_TO_EXT_BIN = False

# how many reads can be hold in memory by default
DEFAULT_SEQS_IN_MEM_LIMIT = 500000

# max width of a line of an ASCII plot
MAX_WIDTH_ASCII_PLOT = 100

# default minimum number of bins in an histogram
MIN_BINS = 20
# default maximum number of bins in an histogram
MAX_BINS = 500
MEAN_VALUES_IN_BIN = 10000

# default number of location to plot in a nucleotide frequency plot
DEF_PLOT_FREQS_UP_TO_BASE = 40

#when 2 match parts are in this distance they are merges as just one matchpart
DEFAULT_IGNORE_ELONGATION_SHORTER = 3

# default kmer size to do the kmer stats
DEFAULT_KMER_SIZE = 20
