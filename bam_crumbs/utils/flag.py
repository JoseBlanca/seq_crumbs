
from operator import or_

# pylint: disable=C0111

SAM_FLAG_BITS = {
    'is_paired': 0x0001,  # the read is paired in sequencing
    'is_in_proper_pair': 0x0002,  # the read is mapped in a proper pair
    'is_unmapped': 0x0004,  # the query sequence itself is unmapped
    'mate_is_unmapped': 0x0008,  # the mate is unmapped
    'strand': 0x0010,  # strand of the query (1 for reverse)
    'mate_strand': 0x0020,  # strand of the mate
    'is_first_in_pair': 0x0040,  # the read is the first read in a pair
    'is_second_in_pair': 0x0080,  # the read is the second read in a pair
    'is_not_primary': 0x0100,  # the alignment is not primary
    'failed_quality': 0x0200,  # the read fails platform/vendor quality checks
    'is_duplicate': 0x0400,  # the read is either a PCR or an optical duplicate
}

SAM_FLAGS = {number: flag for flag, number in SAM_FLAG_BITS.viewitems()}

SAM_FLAG_TAGS = list(sorted(SAM_FLAG_BITS.viewkeys()))
SAM_FLAG_BINARIES = list(sorted(SAM_FLAG_BITS.viewvalues()))


def create_flag(bit_tags):
    'It returns the integer corresponding to the bitwise or of the tags'
    return reduce(or_, [SAM_FLAG_BITS[t] for t in bit_tags])
