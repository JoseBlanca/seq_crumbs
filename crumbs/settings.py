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

import tempfile
import os

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

_SUPPORTED_OUTPUT_FORMATS = ['fasta', 'fastq', 'fastq-illumina']

# number of sequences to process in a chunk. Lenght of the sequence list to
# hold in memory
_PACKET_SIZE = 1000

# number of sequences to analyze in the fastq version guessing of a seekable
# file
_SEQS_TO_GUESS_FASTQ_VERSION = 1000

# number of bytes to analyze in the fastq version guessing of a non-seekable
# file
_CHUNK_TO_GUESS_FASTQ_VERSION = 50000

# maximum length expected for an Illumina read
_LONGEST_EXPECTED_ILLUMINA_READ = 250


# 454 FLX mate pair linker
_FLX_LINKER = 'GTTGGAACCGAAAGGGTTTGAATTCAAACCCTTTCGGTTCCAAC'
# Titanium mate pair linker. It could be found forward or reverse
_TITANIUM_LINKER = 'TCGTATAACTTCGTATAATGTATGCTATACGAAGTTATTACG'
_TITANIUM_LINKER_REV = 'CGTAATAACTTCGTATAGCATACATTATACGAAGTTATACGA'
_FWD_454_LINKERS = [_FLX_LINKER, _TITANIUM_LINKER]
_LINKERS = [SeqRecord(Seq(_FLX_LINKER), id='flx_linker'),
            SeqRecord(Seq(_TITANIUM_LINKER), id='titanium_linker')]


# # Use this to modify how get_binary path works
# if need to modify the binary's name
_USE_EXTERNAL_BIN_PREFIX = False
# prefix to add to the binary name
_EXTERNAL_BIN_PREFIX = 'crumbs_'
# mark True if need the path or assumes that is on the path
_ADD_PATH_TO_EXT_BIN = True

_THIRD_PART_JAVA_DIR = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                                   'third_party', 'java')
_PICARD_TOOLS_DIR = os.path.join(_THIRD_PART_JAVA_DIR, 'picard-tools')

# how many reads can be hold in memory by default
_DEFAULT_SEQS_IN_MEM_LIMIT = 500000

# max width of a line of an ASCII plot
_MAX_WIDTH_ASCII_PLOT = 100

# default minimum number of bins in an histogram
_MIN_BINS = 20
# default maximum number of bins in an histogram
_MAX_BINS = 500
_MEAN_VALUES_IN_BIN = 10000

# default number of location to plot in a nucleotide frequency plot
_DEF_PLOT_FREQS_UP_TO_BASE = 40

# when 2 match parts are in this distance they are merges as just one matchpart
_DEFAULT_IGNORE_ELONGATION_SHORTER = 3

# default kmer size to do the kmer stats
_DEFAULT_KMER_SIZE = 20

# trimest polyannotator
_POLYA_ANNOTATOR_MIN_LEN = 5
_POLYA_ANNOTATOR_MISMATCHES = 1

# quality trim
_DEFAULT_QUALITY_TRIM_TRESHOLD = 25
_DEFAULT_QUALITY_TRIM_WINDOW = 5

# dust score parameters
_DUST_WINDOWSIZE = 64
_DUST_WINDOWSTEP = 32
_DEFATULT_DUST_THRESHOLD = 7

_TEMP_DIR = None

# min_mapq to use as a filter for maped reads
_DEFAULT_MIN_MAPQ = 0


class _Settings(dict):
    '''A class that stores the seq_crumbs settings.'''
    def __init__(self):
        'It inits the class'
        super(_Settings, self).__init__()
        self.load_settings()
        tempfile.tempdir = self.__getitem__('TEMP_DIR')

    def load_settings(self):
        'It loads the settings defined in this module'
        for key, val in globals().viewitems():
            if not key.isupper():
                continue
            key = key[1:]  # strip the underscore
            super(_Settings, self).__setitem__(key, val)

        # Are there any environment variable to update the settings?
        for key, value in os.environ.items():
            if key.startswith('SEQ_CRUMBS_'):
                key = key[11:]
                if key in self.viewkeys():
                    value = type(key)(value)
                    super(_Settings, self).__setitem__(key, value)

_settings = _Settings()


def get_settings():
    'It returns the settings'
    # load the settings defined in this module
    return _settings


def get_setting(key):
    'It returns the value for one setting'
    return _settings[key]
