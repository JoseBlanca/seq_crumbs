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


class UnknownFormatError(Exception):
    'Raised when the format of a sequence file cannot be guessed'
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
