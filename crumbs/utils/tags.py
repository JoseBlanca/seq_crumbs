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

STDIN = 'stdin'
STDOUT = 'stdout'
INFILES = 'infiles'
OUTFILE = 'output'
# Input format tag when we want to guess
GUESS_FORMAT = 'guess'
PROCESSED_SEQS = 'processed_seqs'
PROCESSED_PACKETS = 'processed_packets'
YIELDED_SEQS = 'yielded_seqs'
UPPERCASE = 'upper'
LOWERCASE = 'lower'
SWAPCASE = 'swap'
FWD = 'fwd'
REV = 'rev'
TRIMMING_RECOMMENDATIONS = 'trimming_recommendations'
VECTOR = 'vector'
QUALITY = 'quality'
OTHER = 'other'
TRIMMING_KINDS = [VECTOR, QUALITY, OTHER]
ELONGATED = 'elongated'
SUBJECT = 'subject'
QUERY = 'query'
START = 0
END = 1
BGZF = 'bgzf'
GZIP = 'gzip'

NUCL = 'nucl'
PROT = 'prot'

ERROR_ENVIRON_VARIABLE = 'seq_crumbs_error_debugging'

FIVE_PRIME = "5'"
THREE_PRIME = "3'"

REVERSE = -1
FORWARD = 1
