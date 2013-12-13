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

# pylint: disable=W0232
# pylint: disable=R0903
# pylint: disable=C0111


def error_quality_disagree(error_msg):
    'It checks if a Biopython ValueError is due to a malformed quality file'
    msg = 'Sequence length and number of quality scores disagree'
    msg2 = 'Lengths of sequence and quality values differs'
    if msg in str(error_msg) or msg2 in str(error_msg):
        return True
    else:
        return False


class UnknownFormatError(Exception):
    'Raised when the format of a sequence file cannot be guessed'
    pass


class FileIsEmptyError(Exception):
    pass


class WrongFormatError(Exception):
    'Raised when the file have an unexpected format'
    pass


class UndecidedFastqVersionError(Exception):
    'The file is Fastq, but the version is difficult to guess'
    pass


class FileNotFoundError(Exception):
    'The file does not exists'
    pass


class TooManyFiles(Exception):
    'Too many files given'
    pass


class MalformedFile(Exception):
    'The input file is malformed'
    pass


class SampleSizeError(Exception):
    'The asked sample size is bigger than the population'
    pass


class ExternalBinaryError(Exception):
    'A call to an external binary failed'
    pass


class MissingBinaryError(Exception):
    'A call to an external binary failed becase the binary is not in the path'
    pass


class IncompatibleFormatError(Exception):
    'Asked output format is incompatible with input file formats'
    pass


class MaxNumReadsInMem(Exception):
    'Number of reads in memory reached the limit'
    pass


class PairDirectionError(Exception):
    'Error in parser trying to get direction from description'
    pass


class InterleaveError(Exception):
    'Error in the alternation of forward and reverse reads'
    pass


class OptionalRequirementError(Exception):
    'An optional module is not present'
    pass


class ItemsNotSortedError(Exception):
    'A list or generator of items is not sorted as expected'
    pass
