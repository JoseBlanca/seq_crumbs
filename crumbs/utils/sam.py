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

from operator import or_

# pylint: disable=C0111


IS_PAIRED = 0x0001           # the read is paired in sequencing
IS_IN_PROPER_PAIR = 0x0002   # the read is mapped in a proper pair
IS_UNMAPPED = 0x0004         # the query sequence itself is unmapped
MATE_IS_UNMAPPED = 0x0008    # the mate is unmapped
STRAND = 0x0010              # strand of the query (1 for reverse)
MATE_STRAND = 0x0020         # strand of the mate
IS_FIRST_IN_PAIR = 0x0040    # the read is the first read in a pair
IS_SECOND_IN_PAIR = 0x0080   # the read is the second read in a pair
IS_NOT_PRIMARY = 0x0100      # the alignment is not primary
FAILED_QUALITY = 0x0200      # the read fails platform/vendor quality checks
IS_DUPLICATE = 0x0400        # the read is either a PCR or an optical duplicate


SAM_FLAG_BITS = [IS_PAIRED, IS_IN_PROPER_PAIR, IS_UNMAPPED, MATE_IS_UNMAPPED,
                 STRAND, MATE_STRAND, IS_FIRST_IN_PAIR, IS_SECOND_IN_PAIR,
                 IS_NOT_PRIMARY, FAILED_QUALITY, IS_DUPLICATE]


def bit_tags_to_int_flag(bit_tags):
    'It returns the integer corresponding to the given list of tags'
    return reduce(or_, bit_tags)


def int_flag_to_bit_tags(flag):
    'It returns a list with the indexes of the bits set to 1 in the given flag'
    return [num for num in SAM_FLAG_BITS if num & flag]


def bit_tag_is_in_int_flag(bit_flag, int_flag):
    return bit_flag in int_flag_to_bit_tags(int_flag)
